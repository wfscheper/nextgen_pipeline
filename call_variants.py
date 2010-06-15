import re
import os

from glob import glob
from ruffus import files, follows, inputs, mkdir, regex, transform

from utils import call, pmsg, unpaired_re, unpaired_strings


GENOME_ANALYSIS_JAR = '/opt/GenomeAnalysisTK/GenomeAnalysisTK.jar'
MEM = 'Xmx4g'

CMD_DICT = {
    # commands
    'gatk': 'java -%s -Djava.io.tmpdir=/tmp/$USER -jar %s' % (MEM, GENOME_ANALYSIS_JAR),
    'picard': 'java -%s -jar /opt/picard-tools/' % (MEM),
    'filterIndels': '/usr/local/bin/filterSingleSampleCalls.pl',
    'samtools': '/usr/local/bin/samtools',
    # data files
    'exome': '/usr/local/share/nextgen_resources/hg19_capture.interval_list',
    'genome': '/usr/local/share/nextgen_resources/human_g1k_v37.fasta',
    'dbsnp': '/usr/local/share/nextgen_resources/dbsnp_130_b37.rod',
    # filters
    'standard_filter': '\"QUAL < 30.0 || AB > 0.75 && DP > 40 || QD < 5.0 || HRun > 5 || ' + \
        'SB > -0.10\"',
    'hard_validate_filter': '\"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\"',
}

# Find candidate intervals for realignment
@follows(recalibrate_quality_scores, mkdir('indel_intervals'))
@transform(recalibrate_quality_scores, regex(r'^(.*)/recal_bam/(.*).recalibrated.bam$'),
        r'\1/indel_intervals/\2.intervals')
def create_intervals(input_file, output_file):
    '''Determine indel candidate intervals'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Creating Intervals', input_file, output_file)
    gatk_cmd = '%(gatk)s ' + \
            '-T RealignerTargetCreator ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(infile)s ' + \
            '-o %(outfile)s'
    call(gatk_cmd % cmd_dict)

# Realign around possible indels
@follows(create_intervals, mkdir('cleaned_bam'))
@transform(create_intervals, regex(r'^(.*)/indel_intervals/(.*).intervals$'),
        inputs([r'\1/recal_bam/\2.recalibrated.bam', r'\1/indel_intervals/\2.intervals']),
        r'\1/cleaned_bam/\2.cleaned.bam')
def local_realignment(input_files, output_file):
    '''Realign reads around candidate indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['recal_bam'] = input_files[0]
    cmd_dict['indel_intervals'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('Local Realignment', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s ' + \
            '-T IndelRealigner ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(recal_bam)s ' + \
            '-targetIntervals %(indel_intervals)s ' + \
            '-mrl 20000 ' + \
            '--output %(outfile)s'
    call(gatk_cmd % cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd % cmd_dict)

# Call Indels
@follows(local_realignment, mkdir('indels'))
@transform(local_realignment, regex(r'^(.*)/cleaned_bam/(.*).cleaned.bam$'),
        r'\1/indels/\2_indels.raw.bed',
        r'\1/indels/\2.detailed.bed')
def indel_genotyping(input_file, output_file, details_file):
    '''Call Indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['details'] = details_file
    pmsg('Indel Genotyping', input_file, output_file)
    gatk_cmd = '%(gatk)s -T IndelGenotyperV2 -R %(genome)s -I %(infile)s ' + \
            '-O %(outfile)s -o %(details)s --verbose'
    call(gatk_cmd % cmd_dict)

# Filter Indels

# Call SNPs

# Create InDel mask

# Filter SNPs

stages_dict = {
    'create_intervals': create_intervals,
    'local_realignment': local_realignment,
    'default': local_realignment
}

