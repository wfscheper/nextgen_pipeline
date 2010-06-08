import re
import os
import logging

from glob import glob
from ruffus import files, follows, inputs, mkdir, regex, transform

from utils import call, pmsg, unpaired_re, unpaired_strings


GENOME_ANALYSIS_JAR = '/opt/GenomeAnalysisTK/GenomeAnalysisTK.jar'

CMD_DICT = {
    # commands
    'gatk': 'java -Xmx4g -Djava.io.tmpdir=/tmp/$USER -jar ' + GENOME_ANALYSIS_JAR,
    'picard': 'java -Xmx4g -jar /opt/picard-tools/',
    'filterIndels': '/usr/local/bin/filterSingleSampleCalls.pl',
    # data files
    'exome': '/usr/local/share/nextgen_resources/hg19_capture.interval_list',
    'genome': '/usr/local/share/nextgen_resources/human_g1k_v37.fasta',
    'dbsnp': '/usr/local/share/nextgen_resources/dbsnp_130_b37.rod',
    # filters
    'standard_filter': '\"QUAL < 30.0 || AB > 0.75 && DP > 40 || QD < 5.0 || HRun > 5 || ' + \
        'SB > -0.10\"',
    'hard_validate_filter': '\"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\"',
}

# Calculate Covariates for Quality Score Recalibration
def count_covariates_generator():
    cwd = os.getcwd()
    for file in glob('%s/sorted_bam/*.bam' % cwd):
        filename = '%s/recal_data/%s' % (cwd, unpaired_strings['recal_data'] % \
                unpaired_re.search(file).groupdict())
        yield [file, filename]

@follows(mkdir('recal_data'))
@files(count_covariates_generator)
def count_covariates(input_file, output_file):
    '''Run CounCovariates on files in sorted_bam/'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Count Covariates', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s -T CountCovariates -R %(genome)s -D %(dbsnp)s -standard ' % cmd_dict + \
            '-I %(infile)s -recalFile %(outfile)s' % cmd_dict
    call(gatk_cmd)

# Apply Quality Score Recalibration

# Generate QS Recalibration Dataplots

# Find candidate intervals for realignment

# Realign around possible indels

# Call Indels

# Filter Indels

# Call SNPs

# Create InDel mask

# Filter SNPs

stages_dict = {
    'count_covariates': count_covariates,
    'default': count_covariates,
}

