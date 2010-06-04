import time
import subprocess
import re
import os
import sys
import logging
from glob import glob
from shutil import move, copy

from ruffus import *

from zipper import zip, unzip
from utils import log_info, call

cdict = {
    'bwa': '/usr/local/bin/bwa',
    'samtools': '/usr/local/bin/samtools',
    'ref': '/media/thumper0/genome/human_g1k_v37.fasta',
    'threads': '4',
    'sampl': '/usr/bin/perl /usr/local/bin/samtools.pl'
}

fastq_regex = re.compile('fastq')
paired_re = re.compile('(?P<line>[\d_]+)_s_(?P<lane>\d+)_(?P<pair>[12])(?P<rest>.*)')
unpaired_re = re.compile('(?P<line>[\d_]+)_s_(?P<lane>\d+)(?P<ext>.*)')

paired_strings = {
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

# Convert fastq to sai
def fastq_to_sai_generator():
    cwd = os.getcwd()
    for item in [[file, '%s/sai/%s' % (cwd, paired_strings['sai'] % paired_re.search(file).groupdict())] for file in glob('%s/fastq/*' % cwd)]:
        yield item

@follows(mkdir('sai'))
@files(fastq_to_sai_generator)
def fastq_to_sai(input_file, output_file):
    '''Convert FASTQ files to SAI files.'''
    cmd_dict = cdict.copy()
    cmd_dict['fastq'] = input_file.strip('.gz')
    cmd_dict['sai'] = output_file.strip('.gz')
    cmd_dict['in_file'] = input_file
    cmd_dict['out_file'] = output_file
    pmsg('FASTQ to SAI', cmd_dict['in_file'], cmd_dict['out_file'])
    #unzip(input_file)
    bwacmd = '%(bwa)s aln -t %(threads)s %(ref)s %(in_file)s > %(sai)s' % cmd_dict
    call(bwacmd)
    #zip(cmd_dict['fastq'])

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
    #unzipped = [file.strip('.gz') for file in input_files]
    pmsg('SAM generation', ', '.join(input_files), output_file.strip('.gz'))
    # Unzip input files
    #unzip(*input_files)
    # Get unzipped filenames
    #unzipped.sort(cmp=saicmp)
    input_files.sort(cmp=saicmp)
    # Run bwa to merge paired ends into SAM file
    #cmd_dict['in_files'] = ' '.join(unzipped)
    cmd_dict['in_files'] = ' '.join(input_files)
    cmd_dict['out_file'] = output_file.strip('.gz')
    cmdstr = '%(bwa)s sampe %(ref)s %(in_files)s > %(out_file)s' % cmd_dict
    call(cmdstr)
    #unzipped.append(cmd_dict['out_file'])
    # Zip up unzipped files.
    #zip(*unzipped)
    #zip(cmd_dict['out_file'])

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
    #unzip(input_file)
    samcmd = '%(samtools)s import %(ref)s.fai %(infile)s %(ofile)s' % cmd_dict
    call(samcmd)
    #zip(cmd_dict['infile'])
    
# Sort BAM file
@follows(sam_to_bam, mkdir('sorted_bam'))
@transform(sam_to_bam, regex(r'^(.*)/bam/(.*).bam'), r'\1/sorted_bam/\2.sorted.bam')
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
    #zip(output_file.strip('.gz'))

"""
# Sort BAM file by name
@follows(sam_to_bam, mkdir('bam_namesorted'))
@transform(sam_to_bam, regex(r'^(.+)/bam/(.+)$'), r'\1/bam_namesorted/\2')
def bam_namesort(input_file, output_file):
    '''Sort BAM files by name.'''
    pmsg('BAM name sort', input_file.strip('.gz'), output_file.strip('.gz'))
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file.strip('.gz')
    cmd_dict['outfile'] = output_file.strip('.gz')
    cmd_dict['outprefix'] = os.path.splitext(cmd_dict['outfile'])[0]
    unzip(input_file)
    samcmd = '%(samtools)s sort -n %(infile)s %(outprefix)s' % cmd_dict
    call(samcmd)
    zip(cmd_dict['infile'], cmd_dict['outfile'])

# Add fields to BAM
@follows(bam_namesort, mkdir('fixmate'))
@transform(bam_namesort, regex(r'^(.+)/bam_namesorted/(.+)$'), r'\1/fixmate/\2')
def bam_fixmate(input_file, output_file):
    '''Add fields to BAM (samtools fixmate).'''
    pmsg('BAM fixmate', input_file.strip('.gz'), output_file.strip('.gz'))
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file.strip('.gz')
    cmd_dict['outfile'] = output_file.strip('.gz')
    unzip(input_file)
    samcmd = '%(samtools)s fixmate %(infile)s %(outfile)s' % cmd_dict
    call(samcmd)
    zip(cmd_dict['infile'], cmd_dict['outfile'])

# Add MD tags
@follows(sort_bam, mkdir('bamfinal'))
@transform(sort_bam, regex(r'^(.+)/sorted_bam/(.+).sorted.bam.gz$'), r'\1/bamfinal/\2.final.bam.gz')
def add_mdtags(input_file, output_file):
    '''Add MD tags to BAM'''
    pmsg('Add MD tags', input_file.strip('.gz'), output_file.strip('.gz'))
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file.strip('.gz')
    cmd_dict['outfile'] = output_file.strip('.gz')
    unzip(input_file)
    samcmd = '%(samtools)s calmd -b %(infile)s %(ref)s > %(outfile)s' % cmd_dict
    call(samcmd)
    zip(cmd_dict['infile'], cmd_dict['outfile'])

# Merge sorted BAMs from different lanes
@follows(add_mdtags, mkdir('merged_bam'))
@collate(add_mdtags, regex(r'^(.*)/bamfinal/(\d+)_s_(.*).final.bam.gz$'), r'\1/merged_bam/\2.merged.bam.gz')
def merge_bam(input_files, output_file):
    '''Merge all BAM files for this line into a single BAM file.'''
    cmd_dict = cdict.copy()
    cmd_dict['ofile'] = output_file.strip('.gz')
    cmd_dict['infiles'] = ' '.join([file.strip('.gz') for file in input_files])
    pmsg('Merge BAMs', cmd_dict['infiles'], cmd_dict['ofile'])
    unzip(*input_files)
    bamcmd = '%(samtools)s merge %(ofile)s %(infiles)s' % cmd_dict
    call(bamcmd)
    to_zip = [file.strip('.gz') for file in input_files]
    to_zip.append(output_file.strip('.gz'))
    zip(*to_zip)

# Sort merged BAM file
@transform(merge_bam, regex(r'^(.*)/merged_bam/(.*).bam.gz'), r'\1/merged_bam/\2.sorted.bam.gz')
def sort_merged_bam(input_file, output_file):
    '''Sort merged BAM files by coordinate.'''
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = os.path.splitext(input_file)[0]
    cmd_dict['outfile'] = os.path.splitext(output_file)[0]
    cmd_dict['outprefix'] = os.path.splitext(cmd_dict['outfile'])[0]
    pmsg('BAM Merged Sort', cmd_dict['infile'], cmd_dict['outfile'])
    unzip(input_file)
    samcmd = '%(samtools)s sort %(infile)s %(outprefix)s' % cmd_dict
    call(samcmd)
    zip(cmd_dict['infile'], cmd_dict['outfile'])    

# Create Pileup from BAM
@follows(bam_index, mkdir('pileup'))
@transform(sort_merged_bam, regex(r'^(.*)/merged_bam/(.*).merged.sorted.bam.gz$'), r'\1/pileup/\2.pileup.gz')
def build_pileup(input_file, output_file):
    '''Build a pileup from merged BAM file.'''
    pmsg('Build pileup', input_file.strip('.gz'), output_file.strip('.gz'))
    unzip(input_file)
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file.strip('.gz')
    cmd_dict['ofile'] = output_file.strip('.gz')
    samcmd = '%(samtools)s pileup -c -f %(ref)s %(infile)s > %(ofile)s' % cmd_dict
    call(samcmd)
    zip(cmd_dict['infile'], cmd_dict['ofile'])

# Clean Pileup
@transform(build_pileup, regex(r'^(.*)/pileup/(\d+).pileup.gz$'), r'\1/pileup/\2.pileup.clean.txt.gz')
def clean_pileup(input_file, output_file):
    '''Clean up the pileup.'''
    pmsg('Clean pileup', input_file.strip('.gz'), output_file.strip('.gz'))
    unzip(input_file)
    cmd_dict = cdict.copy()
    cmd_dict['infile'] = input_file.strip('.gz')
    cmd_dict['ofile'] = output_file.strip('.gz')
#
#Original    samcmd = '%(sampl)s varFilter -d 5 -D 70 -Q 25 -q 40 %(infile)s > %(ofile)s' % cmd_dict
    samcmd = '%(sampl)s varFilter -d 5 -D 70 -Q 25 -q 40 %(infile)s > %(ofile)s' % cmd_dict
    call(samcmd)
    zip(cmd_dict['infile'], cmd_dict['ofile'])

# Analyze Pileup
@transform(clean_pileup, regex(r'(.*)/pileup/(.*).pileup.gz'), r'\1/analysis/\2.analysis.txt')
def analyze_pileup(input_files, output_files):
    '''Analyze Pileup file.'''
    time.sleep(1)
"""
