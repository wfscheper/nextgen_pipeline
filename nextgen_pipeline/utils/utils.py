import logging
import os
import re
import subprocess


DEBUG_COMMAND = 'touch %(outfile)s'
CMD_DICT = {
    # defaults
    'bwa': '`which bwa`',
    'filter_indels': '`which filterSingleSampleCalls.pl`',
    'gatk': '`which gatk.sh`',
    'make_indel_mask': '`makeIndelMask.py`',
    'samtools': '`which samtools`',
    'threads': '1',
    'window_size': '200'
}

filename_re = re.compile(r'(?P<line>\w+)_s_(?P<lane>\d+)(?P<ext>.*)')
read_group_re = re.compile(
    r'^.*/(?P<read_group>(?P<sample>[a-zA-Z0-9]+)_(?P<run_barcode>[a-zA-Z0-9]+))' + \
    '(_\d+)?_s_(?P<lane>\d+).*$'
)

logger = logging.getLogger('main')

def call(command, command_dict, is_logged=True, is_debug=False):
    command_line = command % command_dict
    logger.debug('Calling command line: %s' % command_line)
    if is_debug:
        command_line = DEBUG_COMMAND % command_dict
    if is_logged:
        logfile = os.path.split(command_dict['outfile'])[-1] + '.log'
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

