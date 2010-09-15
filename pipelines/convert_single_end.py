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


def copy_sequence_generator():
    for infile in glob('staging_area/*'):
        outfile = '%(line)s_s_%(lane)s.fastq.gz' % \
                unpaired_re.search(infile).groupdict()
        yield [infile, 'fastq/%s' % (outfile)]

# Copy sequence from staging area
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
@transform(copy_sequence, regex(r'^fastq/(.+)\.fastq\.gz$'), r'sai/\1.sai')
def align_sequence(input_file, output_file):
    '''Align sequence files to reference genome'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Aligning sequence', cmd_dict['infile'], cmd_dict['outfile'])
    bwa_cmd = '%(bwa)s aln -t %(threads)s %(genome)s %(infile)s > %(outfile)s'
    call(bwa_cmd, cmd_dict)

@follows(mkdir('sam'))
@transform(align_sequence,
           regex(r'^sai/(.+)\.sai$'),
           inputs([r'fastq/\1.fastq.gz', r'sai/\1.sai']),
           r'sam/\1.sam')
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

@follows(mkdir('sorted'))
@transform(create_sam, regex(r'^sam/(.+)\.sam$'), r'sorted/\1.sorted.sam')
def coordinate_sort_sam(input_file, output_file):
    '''Sort SAM file by coordinates'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('SAM coordinate sort', cmd_dict['infile'], cmd_dict['outfile'])
    picard_cmd = '%(picard)s SortSam ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'SO=coordinate ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

@follows(mkdir('clipped'))
@transform(coordinate_sort_sam,
           regex(r'^sorted/(.+)\.sorted\.sam$'),
           r'clipped/\1.clipped.sam')
def clip_reads(input_file, output_file):
    '''Clip reads around primers'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Read clipping', cmd_dict['infile'], cmd_dict['outfile'])
    clip_cmd = '%(clip_reads)s --out %(outfile)s %(primers)s %(infile)s'
    call(clip_cmd, cmd_dict)

@follows(mkdir('bam'))
@transform(clip_reads, regex(r'^clipped/(.+)\.clipped\.sam$'), r'bam/\1.clipped.bam')
def sam_to_bam(input_file, output_file):
    '''Convert SAM file to BAM'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('SAM to BAM', cmd_dict['infile'], cmd_dict['outfile'])
    picard_cmd = '%(picard)s SamFormatConverter ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

# Update header with missing data
@follows(mkdir('prepped'))
@transform(sam_to_bam, regex(r'^bam/(.+)\.clipped\.bam$'), r'prepped/\1.prepped.bam')
def fix_header(input_file, output_file):
    '''Fix header info'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict.update(read_group_re.match(input_file).groupdict())
    cmd_dict['header_tmp'] = '/tmp/header_%(read_group)s_%(lane)s' % cmd_dict
    open(cmd_dict['header_tmp'], 'w').write(
        open(cmd_dict['header_template'], 'r').read() % cmd_dict
    )
    picard_cmd = '%(picard)s ReplaceSamHeader ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'HEADER=%(header_tmp)s ' + \
            'MAX_RECORDS_IN_RAM=500000 ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)
    os.remove(cmd_dict['header_tmp'])
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict)

@follows(mkdir('coverage'))
@transform(bam_index, regex(r'^prepped/(.+)\.prepped\.bam\.bai$'), r'coverage/\1.coverage')
def calculate_coverage(input_file, output_file):
    '''Calculate coverage statistics'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = os.path.splitext(input_file)[0]
    cmd_dict['outfile'] = output_file
    pmsg('Coverage calculations', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s -T DepthOfCoverage ' + \
            '-I %(infile)s ' + \
            '-L %(se_target)s ' + \
            '-R %(genome)s ' + \
            '-o %(outfile)s ' + \
            '--minMappingQuality 10 ' + \
            '--minBaseQuality 10 ' + \
            ''#'-omitBaseOutput'
    call(gatk_cmd, cmd_dict)

stages_dict = {
    'copy_sequence': copy_sequence,
    'align_sequence': align_sequence,
    'create_sam': create_sam,
    'sort_sam': coordinate_sort_sam,
    'clip_reads': clip_reads,
    'sam_to_bam': sam_to_bam,
    'fix_header': fix_header,
    'calculate_coverage': calculate_coverage,
    'default': calculate_coverage,
}

