import time
import subprocess
import re
import os
import sys
import logging
from Bio import SeqIO
from glob import glob
from shutil import move, copy

from ruffus import *

from zipper import zip, unzip
from utils import log_info, call

cdict = {
    'bwa': '/usr/local/bin/bwa',
    'samtools': '/usr/local/bin/samtools',
    'ref': '../resources/human_g1k_v37.fasta',
    'threads': '4',
    'sampl': '/usr/bin/perl /usr/local/bin/samtools.pl'
}

fastq_regex = re.compile('fastq')
paired_re = re.compile('(?P<line>[\d_]+)_s_(?P<lane>\d+)_(?P<pair>[12])(?P<rest>.*)')
unpaired_re = re.compile('(?P<line>[\d_]+)_s_(?P<lane>\d+)(?P<ext>.*)')

paired_strings = {
    'sequence': '%(line)s_s_%(lane)s_%(pair)s_sequence.txt',
    'fastq': '%(line)s_s_%(lane)s_%(pair)s_sequence.fastq',
    'sai': '%(line)s_s_%(lane)s_%(pair)s.sai',
}
unpaired_strings = {
    'sam': '%(line)s_s_%(lane)s.sam',
    'bam': '%(line)s_s_%(lane)s.bam',
}

logger = logging.getLogger('main')

def pmsg(msg, input, output):
    msgs = {
        'msg': msg,
        'in': input,
        'out': output,
    }
    print 'Running %(msg)s with input %(in)s and output %(out)s' % msgs

# Copy sequence from staging area
def sequence_copy_generator():
    cwd = os.getcwd()
    for file in glob('/media/thumper1/nextgen/staging_area/*_sequence.txt'):
        filename = paired_strings['sequence'] % paired_re.search(file).groupdict()
        yield [file, '%s/fastq/%s.fastq.gz' % (cwd, filename.strip('.txt'))]

@follows(mkdir('fastq'))
@files(sequence_copy_generator)
def copy_sequences(input_file, output_file):
    """Copy sequence files from staging area on thumper1"""
    cmd_dict = cdict.copy()
    cmd_dict['in_file'] = input_file
    cmd_dict['out_file'] = output_file
    pmsg('Sequence Copy', cmd_dict['in_file'], cmd_dict['out_file'])
    SeqIO.convert(input_file, 'fastq-illumina', output_file.strip('.gz'), 'fastq-sanger')
    zip(output_file)

# Convert fastq to sai
def fastq_to_sai_generator():
    cwd = os.getcwd()
    for file in glob('%s/fastq/*' % cwd):
        filename = '%s/sai/%s' % (cwd, paired_strings['sai'] % paired_re.search(file).groupdict())
        yield [file, filename]

@follows(copy_sequences, mkdir('sai'))
@files(fastq_to_sai_generator)
def fastq_to_sai(input_file, output_file):
    '''Convert FASTQ files to SAI files.'''
    cmd_dict = cdict.copy()
    cmd_dict['in_file'] = input_file
    cmd_dict['out_file'] = output_file
    pmsg('FASTQ to SAI', cmd_dict['in_file'], cmd_dict['out_file'])
    bwacmd = '%(bwa)s aln -t %(threads)s %(ref)s %(in_file)s > %(sai)s' % cmd_dict
    call(bwacmd)

# Merge paired ends to SAM
@follows(fastq_to_sai, mkdir('sam'))
@transform(fastq_to_sai,
            regex(r'^(.+)/sai/([\d_]+)_s_(\d)_1.sai$'),
            inputs([r'\1/sai/\2_s_\3_2.sai',
                   r'\1/sai/\2_s_\3_1.sai',
                   r'\1/fastq/\2_s_\3_1_sequence.fastq.gz',
                   r'\1/fastq/\2_s_\3_2_sequence.fastq.gz']),
            r'\1/sam/\2_s_\3.sam')
def paired_ends_to_sam(input_files, output_file):
    '''Convert SAI files and FASTQ files to SAM files.'''

    def saicmp(x,y):
        '''Compare function for moving sai files to front of list'''
        if x.endswith('sai') and not y.endswith('sai'):
            return -1
        elif y.endswith('sai') and not x.endswith('sai'):
            return 1
        else:
            return cmp(x,y)

    cmd_dict = cdict.copy()
    assert type(input_files) is type([])
    pmsg('SAM generation', ', '.join(input_files), output_file.strip('.gz'))
    # sort input files
    input_files.sort(cmp=saicmp)
    # Run bwa to merge paired ends into SAM file
    cmd_dict['in_files'] = ' '.join(input_files)
    cmd_dict['out_file'] = output_file.strip('.gz')
    cmdstr = '%(bwa)s sampe %(ref)s %(in_files)s > %(out_file)s' % cmd_dict
    call(cmdstr)

## Convert filtered SAM files to BAM files
@follows(paired_ends_to_sam, mkdir('bam'))
@transform(paired_ends_to_sam, regex(r'^(.*)/sam/(.*).sam$'), r'\1/bam/\2.bam')
def sam_to_bam(input_file, output_file):
    '''Convert SAM files to BAM files.'''
    cmd_dict = cdict.copy()
    finfo = unpaired_re.search(input_file).groupdict()
    cmd_dict['infile'] = input_file
    cmd_dict['ofile'] = output_file
    pmsg('SAM to BAM', cmd_dict['infile'], cmd_dict['ofile'])
    samcmd = '%(samtools)s import %(ref)s.fai %(infile)s %(ofile)s' % cmd_dict
    call(samcmd)
    
# Sort BAM file by name
@follows(sam_to_bam, mkdir('namesorted_bam'))
@transform(sam_to_bam, regex(r'^(.*)/bam/(.*).bam$'), r'\1/namesorted_bam/\2.namesorted.bam')
def namesort_bam(input_file, output_file):
    '''Sort BAM files by name.'''
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['outprefix'] = os.path.splitext(output_file)[0]
    pmsg('BAM Name Sort', cmd_dict['infile'], cmd_dict['outfile'])
    samcmd = '%(samtools)s sort -n %(infile)s %(outprefix)s' % cmd_dict
    call(samcmd)

# Run samtools fixmate on namesorted BAM file
@follows(namesort_bam, mkdir('fixmate_bam'))
@transform(namesort_bam, regex(r'^(.*)/namesorted_bam/(.*).namesorted.bam$'), r'\1/fixmate_bam/\2.fixmate.bam')
def fixmate_bam(input_file, output_file):
    '''Fix mate info in BAM file.'''
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Fix Mate Info', cmd_dict['infile'], cmd_dict['outfile'])
    samcmd = '%(samtools)s fixmate %(infile)s %(outfile)s' % cmd_dict
    call(samcmd)

# Sort BAM file
@follows(fixmate_bam, mkdir('sorted_bam'))
@transform(fixmate_bam, regex(r'^(.*)/fixmate_bam/(.*).fixmate.bam'), r'\1/sorted_bam/\2.sorted.bam')
def sort_bam(input_file, output_file):
    '''Sort BAM files by coordinate.'''
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = os.path.splitext(output_file)[0]
    cmd_dict['outprefix'] = cmd_dict['outfile']
    pmsg('BAM Coord Sort', cmd_dict['infile'], cmd_dict['outfile'])
    samcmd = '%(samtools)s sort %(infile)s %(outprefix)s' % cmd_dict
    call(samcmd)

# Create index from BAM - creates BAI files
@follows(sort_bam)
@transform(sort_bam, regex(r'^(.*).sorted.bam'), r'\1.sorted.bai')
def bam_index(input_file, output_file):
    '''Index BAM file and create a BAI file.'''
    pmsg('Create BAM Index', input_file.strip('.gz'), output_file.strip('.gz'))
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file.strip('.gz')
    idxcmd = '%(samtools)s index %(infile)s' % cmd_dict
    call(idxcmd)

stages_dict = {
    'copy_sequences': copy_sequences,
    'align_sequences': fastq_to_sai,
    'make_sam': paired_ends_to_sam,
    'make_bam': sam_to_bam,
    'namesort_bam': namesort_bam,
    'fixmate_bam': fixmate_bam,
    'sort_bam': sort_bam,
    'index_bam': bam_index,
    'default': bam_index,
}

