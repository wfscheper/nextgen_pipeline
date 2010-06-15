import os

from glob import glob
from ruffus import files, follows, inputs, mkdir, regex, transform

from utils import call, pmsg, unpaired_re, unpaired_strings, CMD_DICT


# Calculate Covariates for Quality Score Recalibration
def count_covariates_generator():
    cwd = os.getcwd()
    for file in glob('%s/prepped_bam/*.bam' % cwd):
        filename = '%s/recal_data/%s' % (cwd, unpaired_strings['recal_data'] % \
                unpaired_re.search(file).groupdict())
        yield [file, filename]

def count_covariates(input_file, output_file):
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Count Covariates', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s -T CountCovariates -R %(genome)s -D %(dbsnp)s -standard ' + \
            '-I %(infile)s -recalFile %(outfile)s -dP illumina'
    call(gatk_cmd % cmd_dict)

@follows(mkdir('recal_data'))
@files(count_covariates_generator)
def sorted_count_covariates(input_file, output_file):
    '''Run CounCovariates on files in sorted_bam/'''
    count_covariates(input_file, output_file)

# Apply Quality Score Recalibration
@follows(sorted_count_covariates, mkdir('recal_bam'))
@transform(sorted_count_covariates, regex(r'^(.*)/recal_data/(.*).prepped.csv$'),
        inputs([r'\1/recal_data/\2.prepped.csv', r'\1/prepped_bam/\2.prepped.bam']),
        r'\1/recal_bam/\2.recalibrated.bam')
def recalibrate_quality_scores(input_files, output_file):
    '''Apply Recalibrated QSs to BAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['recal_data'] = input_files[0]
    cmd_dict['bam'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('Table Recalibration', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s -T TableRecalibration -R %(genome)s -I %(bam)s ' + \
            '-recalFile %(recal_data)s -outputBam %(outfile)s -dP illumina'
    call(gatk_cmd % cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd % cmd_dict)

# Generate QS Recalibration Covariates
@follows(recalibrate_quality_scores)
@transform(recalibrate_quality_scores, regex(r'^(.*)/recal_bam/(.*).recalibrated.bam$'),
        r'\1/recal_data/\2.recalibrated.csv')
def recal_count_covariates(input_file, output_file):
    '''Run CountCovariates on files in recal_bam/'''
    count_covariates(input_file, output_file)

stages_dict = {
    'count_covariates': sorted_count_covariates,
    'recalibrate_quality': recalibrate_quality_scores,
    'recount_covariates': recal_count_covariates,
    'default': recal_count_covariates
}

