import logging
import re
import subprocess

from run import DEBUG


DEBUG_COMMAND = 'touch %(outfile)s'

CMD_DICT = {
    'threads': '4',
    # commands
    'bwa': '/usr/local/bin/bwa',
    'clip_reads': '/usr/local/bin/clip_se_reads.py',
    'filterIndels': '/usr/local/bin/filterSingleSampleCalls.pl',
    'gatk': '/usr/local/bin/gatk.sh',
    'makeIndelMask': '/usr/local/bin/makeIndelMask.py',
    'picard': '/usr/local/bin/picard.sh',
    'sampl': '/usr/bin/perl /usr/local/bin/samtools.pl',
    'samtools': '/usr/local/bin/samtools',
    # data files
    'dbsnp': '../resources/dbsnp_130_b37.rod',
    'exome': '../resources/hg19_capture.interval_list',
    'genome': '../resources/human_g1k_v37.fasta',
    'header_template': '../resources/header.template',
    'header_tmp': '/tmp/header-%(read_group)s',
    'primers': '../resources/ceph_raindance_1_primers_hg19.bed',
    # filters
    'standard_filter': '\"QUAL < 30.0 || AB > 0.75 && DP > 40 || QD < 5.0 || HRun > 5 || ' + \
        'SB > -0.10\"',
    'hard_to_validate_filter': '\"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\"',
}

paired_re = re.compile(r'(?P<line>\w+)_s_(?P<lane>\d+)_(?P<pair>[12])(?P<rest>.*)')
unpaired_re = re.compile(r'(?P<line>\w+)_s_(?P<lane>\d+)(?P<ext>.*)')
read_group_re = re.compile(
    r'^.*(?P<read_group>(?P<sample>[a-zA-Z0-9]+)_(?P<run_barcode>[a-zA-Z0-9]+))' + \
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
    print 'Running %(msg)s with input %(in)s and output %(out)s' % msgs

def call(command, command_dict):
    command_line = command % command_dict
    logger.debug('Calling command line: %s' % command_line)
    if DEBUG: command_line = DEBUG_COMMAND % command_dict
    subprocess.call(command_line, shell=True)

def saicmp(x,y):
    '''Compare function for moving sai files to front of list'''
    if x.endswith('sai') and not y.endswith('sai'):
        return -1
    elif y.endswith('sai') and not x.endswith('sai'):
        return 1
    else:
        return cmp(x,y)

