import logging
import os
import re
import subprocess


DEBUG_COMMAND = 'touch %(outfile)s'
DEBUG = False

CMD_DICT = {
    'threads': '4',
    # commands
    'bwa': '`which bwa`',
    'clip_reads': '`which clip_se_reads.py`',
    'filterIndels': '`which filterSingleSampleCalls.pl`',
    'gatk': '`which gatk.sh`',
    'makeIndelMask': '/usr/bin/env python `which makeIndelMask.py`',
    'picard': '`which picard.sh`',
    'sampl': '/usr/bin/perl `which samtools.pl`',
    'samtools': '`which samtools`',
    # data files
    'dbsnp': '../resources/dbsnp_130_b37.rod',
    'exome': '../resources/hg19_capture.interval_list',
    'genome': '../resources/human_g1k_v37.fasta',
    'header_template': '../resources/header.template',
    'primers': '../resources/ceph_raindance_1_primers_hg19.bed',
    'se_target': '../resources/ceph_target_regions_hg19.interval_list',
    # filters
    'standard_filter': '\"QUAL < 30.0 || AB > 0.75 && DP > 40 || QD < 5.0 || HRun > 5 || ' + \
        'SB > -0.10\"',
    'hard_to_validate_filter': '\"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\"',
}

paired_re = re.compile(r'(?P<line>\w+)_s_(?P<lane>\d+)_(?P<pair>[12])(?P<rest>.*)')
unpaired_re = re.compile(r'(?P<line>\w+)_s_(?P<lane>\d+)(?P<ext>.*)')
read_group_re = re.compile(
    r'^.*/(?P<read_group>(?P<sample>[a-zA-Z0-9]+)_(?P<run_barcode>[a-zA-Z0-9]+))' + \
    '(_\d+)?_s_(?P<lane>\d+).*$'
)
paired_strings = {
    'sequence': '%(line)s_s_%(lane)s_%(pair)s_sequence.txt',
    'fastq': '%(line)s_s_%(lane)s_%(pair)s_sequence.fastq',
    'sai': '%(line)s_s_%(lane)s_%(pair)s.sai',
}
unpaired_strings = {
    'sai': '%(line)s_s_%(lane)s.sai',
    'sam': '%(line)s_s_%(lane)s.sam',
    'bam': '%(line)s_s_%(lane)s.bam',
    'recal_data': '%(line)s_s_%(lane)s.prepped.csv',
    'intervals': '%(line)s_s_%(lane)s.intervals',
}

logger = logging.getLogger('main')

def call(command, command_dict, is_logged=True):
    command_line = command % command_dict
    logger.debug('Calling command line: %s' % command_line)
    if DEBUG:
        command_line = DEBUG_COMMAND % command_dict
    if is_logged:
        logfile = os.path.split(command_dict['outfile'])[-1]
        logfile = os.path.join('logs', logfile)
        subprocess.call(command_line, stdout=open(logfile, 'w'), stderr=subprocess.STDOUT,
                        shell=True)
    else:
        subprocess.call(command_line, shell=True)

def check_if_clean(input_files, output_files):
    '''Check if intermediate files from previous pipeline are still around'''
    for file in input_files:
        if os.path.exists(file):
            return True, "Need to clean out intermediate file: %s" % file
    return False, "Intermediate files cleaned out"

def log_info(func):
    def new_func(*args, **kwargs):
        logger.debug(func.__name__)
        logger.debug('inputs: %s' % str(args[0]))
        logger.debug('outputs: %s' % str(args[1]))
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func

def pmsg(msg, input, output):
    msgs = {
        'msg': msg,
        'in': input,
        'out': output,
    }
    print '%(msg)s with input %(in)s and output %(out)s' % msgs

def saicmp(x, y):
    '''Compare function for moving sai files to front of list'''
    if x.endswith('sai') and not y.endswith('sai'):
        return - 1
    elif y.endswith('sai') and not x.endswith('sai'):
        return 1
    else:
        return cmp(x, y)

