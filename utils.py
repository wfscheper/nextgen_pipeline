import logging
import re
import subprocess

fastq_re = re.compile('fastq')
paired_re = re.compile('(?P<line>\w+)_s_(?P<lane>\d+)_(?P<pair>[12])(?P<rest>.*)')
unpaired_re = re.compile('(?P<line>\w+)_s_(?P<lane>\d+)(?P<ext>.*)')

paired_strings = {
    'sequence': '%(line)s_s_%(lane)s_%(pair)s_sequence.txt',
    'fastq': '%(line)s_s_%(lane)s_%(pair)s_sequence.fastq',
    'sai': '%(line)s_s_%(lane)s_%(pair)s.sai',
}
unpaired_strings = {
    'sam': '%(line)s_s_%(lane)s.sam',
    'bam': '%(line)s_s_%(lane)s.bam',
    'recal_data': '%(line)s_s_%(lane)s.recal_data.csv',
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

def call(cmdline):
    logger.debug('Calling command line: %s' % cmdline)
    subprocess.call(cmdline, shell=True)
