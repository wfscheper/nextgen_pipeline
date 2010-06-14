import logging
import re
import subprocess

fastq_re = re.compile('fastq')
paired_re = re.compile(r'(?P<line>\w+)_s_(?P<lane>\d+)_(?P<pair>[12])(?P<rest>.*)')
unpaired_re = re.compile(r'(?P<line>\w+)_s_(?P<lane>\d+)(?P<ext>.*)')
read_group_re = re.compile(r'^(?P<sample>[a-zA-Z0-9]+)_(?P<run_barcode>[a-zA-Z0-9]+)(_\d+)?_s_(?P<lane>\d+)$')

paired_strings = {
    'sequence': '%(line)s_s_%(lane)s_%(pair)s_sequence.txt',
    'fastq': '%(line)s_s_%(lane)s_%(pair)s_sequence.fastq',
    'sai': '%(line)s_s_%(lane)s_%(pair)s.sai',
}
unpaired_strings = {
    'sam': '%(line)s_s_%(lane)s.sam',
    'bam': '%(line)s_s_%(lane)s.bam',
    'recal_data': '%(line)s_s_%(lane)s.original.csv',
}

GENOME_ANALYSIS_JAR = '/opt/GenomeAnalysisTK/GenomeAnalysisTK.jar'
MEM = 'Xmx4g'

CMD_DICT = {
    'threads': '4',
    # commands
    'bwa': '/usr/local/bin/bwa',
    'gatk': 'java -%s -Djava.io.tmpdir=/tmp/$USER -jar %s' % (MEM, GENOME_ANALYSIS_JAR),
    'filterIndels': '/usr/local/bin/filterSingleSampleCalls.pl',
    'picard': '/usr/local/bin/picard.sh',
    'sampl': '/usr/bin/perl /usr/local/bin/samtools.pl',
    'samtools': '/usr/local/bin/samtools',
    # data files
    'exome': '../resources/hg19_capture.interval_list',
    'genome': '../resources/human_g1k_v37.fasta',
    'dbsnp': '../resources/dbsnp_130_b37.rod',
    'header_template': '../resources/header.template',
    'header_tmp': '/tmp/header-%(read_group)s',
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
