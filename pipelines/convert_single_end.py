'''
Prepares bam files from staged sequence.txt files

Required for any pipelines that use the Genome Analysis Toolkit.

Currnetly, this means quality score recalibration and variant calling.
'''
import gzip
import os
from Bio import SeqIO
from glob import iglob as glob

from ruffus import follows, files, inputs, jobs_limit, mkdir, regex, transform

from zipper import zip
from utils import call, unpaired_re, pmsg, read_group_re, saicmp, CMD_DICT


# Copy sequence from staging area
def copy_sequence_generator():
    cwd = os.getcwd()
    for infile in glob('./staging_area/*'):
        outfile = '%(line)s_s_%(lane)s_sequence.fastq.gz' % \
                unpaired_re.search(infile).groupdict()
        yield [infile, '%s/fastq/%s' % (cwd, outfile)]

@follows(mkdir('fastq'))
@files(copy_sequence_generator)
def copy_sequence(input_file, output_file):
    """Copy sequence files from staging area on thumper1"""
    cmd_dict = CMD_DICT.copy()
    cmd_dict['outfile'] = output_file
    cmd_dict['outfile_prefix'] = output_file.strip('.gz')

    if input_file.endswith('gz'):
        infile = gzip.open(input_file)
    else:
        infile = open(input_file)

    pmsg('Copying sequence file', input_file, cmd_dict['outfile_prefix'])
    SeqIO.convert(infile, 'fastq-illumina', cmd_dict['outfile_prefix'], 'fastq-sanger')
    pmsg('Compressing sequence file', cmd_dict['outfile_prefix'], cmd_dict['outfile'])
    zip(cmd_dict['outfile_prefix'])

@jobs_limit(2)
@follows(mkdir('sai'))
@transform(copy_sequence, regex(r'^(.*)/fastq/(.*)\.fastq\.gz$'), r'\1/sai/\2.sai')
def align_sequence(input_file, output_file):
    '''Align sequence files to reference genome'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Aligning sequence', cmd_dict['infile'], cmd_dict['outfile'])
    bwa_cmd = '%(bwa)s aln -t %(threads)s %(genome)s %(infile)s > %(outfile)s'
    call(bwa_cmd, cmd_dict)

@follows(mkdir('sam'))
@transform(align_sequence, regex(r'^(.*)/sai/(.*)\.sai$'),
        inputs([r'\1/fastq/\2.fastq.gz',
            r'\1/sai/\2.sai']),
        r'\1/sam/\2.sam')
def create_sam(input_files, output_file):
    '''Convert fastq+sai files to sam file'''
    cmd_dict = CMD_DICT.copy()
    assert type(input_files) is type([])
    pmsg('SAM generation', ', '.join(input_files), output_file.strip('.gz'))
    # sort input files
    input_files.sort(cmp=saicmp)
    # Run bwa to merge paired ends into SAM file
    cmd_dict['infiles'] = ' '.join(input_files)
    cmd_dict['outfile'] = output_file
    bwa_cmd = '%(bwa)s samse %(genome)s %(infiles)s > %(outfile)s'
    call(bwa_cmd, cmd_dict)

## Convert filtered SAM files to BAM files
@follows(mkdir('bam'))
@transform(create_sam, regex(r'^(.*)/sam/(.+)\.sam$'), r'\1/bam/\2.bam')
def sam_to_bam(input_file, output_file):
    '''Convert SAM files to BAM files.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('SAM to BAM', cmd_dict['infile'], cmd_dict['outfile'])
    sam_cmd = '%(samtools)s import %(genome)s.fai %(infile)s %(outfile)s'
    call(sam_cmd, cmd_dict)

stages_dict = {
    'copy_sequence': copy_sequence,
    'align_sequence': align_sequence,
    'create_sam': create_sam,
    'sam_to_bam': sam_to_bam,
    'default': sam_to_bam,
}

