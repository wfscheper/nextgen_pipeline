'''
Recalibrates quality scores and generates data needed for checking what effect
recalibration had on quality scores. 

Depends on a properly aligned, sorted and indexed bam file.

Required for depth of coverage and variant calling.
'''
import os

from glob import iglob as glob
from ruffus import check_if_uptodate, files, follows, inputs, jobs_limit, mkdir, regex, transform

from utils import call, check_if_clean, pmsg, unpaired_re, unpaired_strings, CMD_DICT


def count_covariates_generator():
    for file in glob('prepped/*.bam'):
        filename = 'covariates/%s' % (unpaired_strings['recal_data'] % \
                unpaired_re.search(file).groupdict())
        yield [file, filename]

def call_count_covariates(input_file, output_file):
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Count Covariates', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s ' + \
            '-T CountCovariates ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(infile)s ' + \
            '-recalFile %(outfile)s ' + \
            '-standard '
    call(gatk_cmd, cmd_dict)

@jobs_limit(1)
@files(["bam/", "clipped/", "deduped/", "header/", "sam/", "sorted/"], None)
@check_if_uptodate(check_if_clean)
def clean_up(input_files, output_file):
    '''Clean up intermediate files'''
    print('Cleaning up intermeidate files: %s' % ", ".join(input_files))
    call('rm -rf %s' % " ".join(input_files), {}, is_logged=False)

# Calculate Covariates for Quality Score Recalibration
@follows(clean_up, mkdir('covariates'))
@files(count_covariates_generator)
def count_covariates(input_file, output_file):
    '''Run CounCovariates on files in sorted/'''
    call_count_covariates(input_file, output_file)

# Apply Quality Score Recalibration
@follows(mkdir('recalibrated'))
@transform(count_covariates,
           regex(r'^covariates/(.+).prepped.csv$'),
           inputs([r'covariates/\1.prepped.csv', r'prepped/\1.prepped.bam']),
           r'recalibrated/\1.recalibrated.bam')
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
            '-o %(outfile)s'
    call(gatk_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict)

# Generate QS Recalibration Covariates
@transform(recalibrate_quality_scores,
           regex(r'^recalibrated/(.+).recalibrated.bam$'),
           r'covariates/\1.recalibrated.csv')
def recount_covariates(input_file, output_file):
    '''Run CountCovariates on files in recalibrated/'''
    call_count_covariates(input_file, output_file)

# Find candidate intervals for realignment
@follows(recount_covariates, mkdir('intervals'))
@transform(recalibrate_quality_scores,
           regex(r'^recalibrated/(.+)\.recalibrated\.bam$'),
           r'intervals/\1.intervals')
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
@follows(mkdir('realigned'))
@transform(create_intervals,
           regex(r'^intervals/(.+)\.intervals$'),
           inputs([r'recalibrated/\1.recalibrated.bam', r'intervals/\1.intervals']),
            r'realigned/\1.realigned.bam')
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
            '-o %(outfile)s'
    call(gatk_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict)

# Fix mate info post realignment
@follows(local_realignment, mkdir('fixmate'))
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
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'SO=coordinate ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict)

stages_dict = {
    'clean_up': clean_up,
    'count': count_covariates,
    'recalibrate': recalibrate_quality_scores,
    'recount': recount_covariates,
    'intervals': create_intervals,
    'realignment': local_realignment,
    'fixmate': fix_mate_realigned,
    'default': fix_mate_realigned,
}

