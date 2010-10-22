'''
Recalibrates quality scores and generates data needed for checking what effect
recalibration had on quality scores.

Depends on a properly aligned, sorted and indexed bam file.

Required for depth of coverage and variant calling.
'''
import os

from glob import iglob as glob
from ruffus import check_if_uptodate, files, follows, inputs, jobs_limit, mkdir, regex, transform

from utils import CMD_DICT, call, check_if_clean, pmsg, filename_re


def call_count_covariates(input_file, output_file):
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Count Covariates', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s --analysis_type CountCovariates ' + \
            '--reference_sequence %(reference)s ' + \
            '--DBSNP %(dbsnp)s ' + \
            '--input_file %(infile)s ' + \
            '--recal_file %(outfile)s ' + \
            '--standard_covs ' + \
            '--num_threads 2 '
    call(gatk_cmd, cmd_dict)

def create_intervals_generator():
    for infile in glob('deduped/*.bam'):
        outfile = '%(line)s_s_%(lane)s.intervals' % filename_re.search(infile).groupdict()
        yield [infile, 'intervals/%s' % outfile]

@jobs_limit(1)
@files(['bam/', 'clipped/', 'sam/', 'sorted/'], None)
@check_if_uptodate(check_if_clean)
def clean_up(input_files, output_file):
    '''Clean up intermediate files'''
    print('Cleaning up intermeidate files: %s' % ', '.join(input_files))
    call('rm -rf %s' % ' '.join(input_files), {}, is_logged=False)

# Find candidate intervals for realignment
@follows(clean_up, mkdir('intervals', 'logs'))
@files(create_intervals_generator)
def create_intervals(input_file, output_file):
    '''Determine indel candidate intervals'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Interval Creation', input_file, output_file)
    gatk_cmd = '%(gatk)s --analysis_type RealignerTargetCreator ' + \
            '--reference_sequence %(reference)s ' + \
            '--DBSNP %(dbsnp)s ' + \
            '--input_file %(infile)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)

# Realign around possible indels
@follows(mkdir('realigned'))
@transform(create_intervals,
           regex(r'^intervals/(.+)\.intervals$'),
           inputs([r'deduped/\1.deduped.bam', r'intervals/\1.intervals']),
           r'realigned/\1.realigned.bam')
def local_realignment(input_files, output_file):
    '''Realign reads around candidate indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['bam_file'] = input_files[0]
    cmd_dict['indel_intervals'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('Local Realignment', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s --analysis_type IndelRealigner ' + \
            '--reference_sequence %(reference)s ' + \
            '--DBSNP %(dbsnp)s ' + \
            '--input_file %(bam_file)s ' + \
            '--targetIntervals %(indel_intervals)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)

# Fix mate info post realignment
@follows(mkdir('fixmate'))
@transform(local_realignment,
           regex(r'^realigned/(.+)\.realigned\.bam$'),
           r'fixmate/\1.fixmate.bam')
def fix_mate_realigned(input_file, output_file):
    '''Fix mate info post-realignment'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Fix Mate Info', cmd_dict['infile'], cmd_dict['outfile'])
    picard_cmd = '%(picard)s FixMateInformation ' + \
            'INPUT=%(infile)s ' + \
            'OUTPUT=%(outfile)s ' + \
            'SORT_ORDER=coordinate ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict, is_logged=False)

# Calculate Covariates for Quality Score Recalibration
@follows(mkdir('covariates'))
@transform(fix_mate_realigned,
           regex(r'^fixmate/(.+)\.fixmate.bam'),
           r'covariates/\1.precalibration.csv')
def count_covariates(input_file, output_file):
    '''Run CounCovariates on files in sorted/'''
    call_count_covariates(input_file, output_file)

# Apply Quality Score Recalibration
@follows(mkdir('recalibrated'))
@transform(count_covariates,
           regex(r'^covariates/(.+).precalibration.csv$'),
           inputs([r'covariates/\1.precalibration.csv', r'fixmate/\1.fixmate.bam']),
           r'recalibrated/\1.recalibrated.bam')
def recalibrate_quality_scores(input_files, output_file):
    '''Apply Recalibrated QSs to BAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['recal_data'] = input_files[0]
    cmd_dict['bam'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('QS Recalibration', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s --analysis_type TableRecalibration ' + \
            '--reference_sequence %(reference)s ' + \
            '--input_file %(bam)s ' + \
            '--recal_file %(recal_data)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict, is_logged=False)

# Generate QS Recalibration Covariates
@transform(recalibrate_quality_scores,
           regex(r'^recalibrated/(.+).recalibrated.bam$'),
           r'covariates/\1.postcalibration.csv')
def recount_covariates(input_file, output_file):
    '''Run CountCovariates on files in recalibrated/'''
    call_count_covariates(input_file, output_file)

stages_dict = {
    'clean_up': clean_up,
    'intervals': create_intervals,
    'realignment': local_realignment,
    'fixmate': fix_mate_realigned,
    'count': count_covariates,
    'recalibrate': recalibrate_quality_scores,
    'recount': recount_covariates,
    'default': recount_covariates,
}

