'''
Recalibrates quality scores and generates data needed for checking what effect
recalibration had on quality scores. 

Depends on a properly aligned, sorted and indexed bam file.

Required for depth of coverage and variant calling.
'''
import os

from glob import iglob as glob
from ruffus import files, follows, inputs, mkdir, regex, transform

from utils import call, pmsg, unpaired_re, unpaired_strings, CMD_DICT


# Calculate Covariates for Quality Score Recalibration
def count_covariates_generator():
    cwd = os.getcwd()
    for file in glob('%s/prepped/*.bam' % cwd):
        filename = '%s/recal_data/%s' % (cwd, unpaired_strings['recal_data'] % \
                unpaired_re.search(file).groupdict())
        yield [file, filename]

def count_covariates(input_file, output_file):
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    print 'Cleaning up intermeidate files'
    call('rm -rf bam/ clipped/ sam/ sorted/', {})
    pmsg('Count Covariates', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s ' + \
            '-T CountCovariates ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(infile)s ' + \
            '-recalFile %(outfile)s ' + \
            '-standard '
    call(gatk_cmd, cmd_dict)

@follows(mkdir('recal_data'))
@files(count_covariates_generator)
def sorted_count_covariates(input_file, output_file):
    '''Run CounCovariates on files in sorted/'''
    count_covariates(input_file, output_file)

# Apply Quality Score Recalibration
@follows(sorted_count_covariates, mkdir('recalibrated'))
@transform(sorted_count_covariates, regex(r'^(.*)/recal_data/(.*).prepped.csv$'),
        inputs([r'\1/recal_data/\2.prepped.csv', r'\1/prepped/\2.prepped.bam']),
        r'\1/recalibrated/\2.recalibrated.bam')
def recalibrate_quality_scores(input_files, output_file):
    '''Apply Recalibrated QSs to BAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['recal_data'] = input_files[0]
    cmd_dict['bam'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('QS Recalibration', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s ' + \
            '-T TableRecalibration ' + \
            '-R %(genome)s ' + \
            '-I %(bam)s ' + \
            '-recalFile %(recal_data)s ' + \
            '-outputBam %(outfile)s'
    call(gatk_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict)

# Generate QS Recalibration Covariates
@follows(recalibrate_quality_scores)
@transform(recalibrate_quality_scores, regex(r'^(.*)/recalibrated/(.*).recalibrated.bam$'),
        r'\1/recal_data/\2.recalibrated.csv')
def recal_count_covariates(input_file, output_file):
    '''Run CountCovariates on files in recalibrated/'''
    count_covariates(input_file, output_file)

# Find candidate intervals for realignment
@follows(mkdir('indel_intervals'))
@transform(recalibrate_quality_scores, regex(r'^(.*)/recalibrated/(.+)\.recalibrated\.bam$'),
        r'\1/indel_intervals/\2.intervals')
def create_intervals(input_file, output_file):
    '''Determine indel candidate intervals'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Interval Creation', input_file, output_file)
    gatk_cmd = '%(gatk)s ' + \
            '-T RealignerTargetCreator ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(infile)s ' + \
            '-o %(outfile)s'
    call(gatk_cmd, cmd_dict)

# Realign around possible indels
@follows(create_intervals, mkdir('realigned'))
@transform(create_intervals, regex(r'^(.*)/indel_intervals/(.+)\.intervals$'),
        inputs([r'\1/recalibrated/\2.recalibrated.bam', r'\1/indel_intervals/\2.intervals']),
        r'\1/realigned/\2.realigned.bam')
def local_realignment(input_files, output_file):
    '''Realign reads around candidate indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['recalibrated'] = input_files[0]
    cmd_dict['indel_intervals'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('Local Realignment', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s ' + \
            '-T IndelRealigner ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(recalibrated)s ' + \
            '-targetIntervals %(indel_intervals)s ' + \
            '--output %(outfile)s'
    call(gatk_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict)

# Fix mate info post realignment
@follows(local_realignment, mkdir('realigned_fixmate'))
@transform(local_realignment, regex(r'^(.*)/realigned/(.+)\.realigned\.bam$'),
        r'\1/realigned_fixmate/\2.realigned_fixmate.bam')
def fix_mate_realigned(input_file, output_file):
    '''Fix mate info post-realignment'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Fix Mate Info', cmd_dict['infile'], cmd_dict['outfile'])
    picard_cmd = '%(picard)s FixMateInformation ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'SO=coordinate ' + \
            'MAX_RECORDS_IN_RAM=7500000 ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict)

stages_dict = {
    'count_covariates': sorted_count_covariates,
    'recalibrate_quality': recalibrate_quality_scores,
    'recount_covariates': recal_count_covariates,
    'create_intervals': create_intervals,
    'local_realignment': local_realignment,
    'fix_mate_realigned': fix_mate_realigned,
    'default': fix_mate_realigned,
}

