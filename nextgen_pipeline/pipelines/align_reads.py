'''
Prepares bam files from staged sequence.txt files

Required for any pipelines that use the Genome Analysis Toolkit.

Currnetly, this means quality score recalibration and variant calling.
'''
import gzip
import os
from Bio import SeqIO
from bz2 import BZ2File
from glob import iglob, glob

from ruffus import follows, files, inputs, merge, mkdir, regex, transform

from ..utils import CMD_DICT, call, pmsg, read_group_re

def copy_sequence_generator():
    for in_file in iglob('staging_area/*'):
        out_file = os.path.split(in_file)[-1]
        out_file = out_file.split(os.path.extsep)[0] + '.fastq.gz'
        out_file = os.path.join('fastq', out_file)
        yield [in_file, out_file]

# Copy sequence from staging area
@follows(mkdir('fastq', 'logs'))
@files(copy_sequence_generator)
def copy_sequence(input_file, output_file):
    '''Copy sequence files from staging area'''
    GZIP_HEADER = '\x1f\x8b'
    BZIP_HEADER = 'BZ'

    pmsg('Copying sequence files', input_file, output_file)
    # check if this is actually a gzipped file
    header = open(input_file).read(2)
    if header == GZIP_HEADER:
        input_file_handle = gzip.open(input_file, 'rb')
    elif header == BZIP_HEADER:
        input_file_handle = BZ2File(input_file, 'r')
    else:
        input_file_handle = open(input_file, 'rb')
    output_file_handle = gzip.open(output_file, 'wb')

    # check whether this is a illumina or sanger fastq file
    try:
        SeqIO.convert(input_file_handle, 'fastq-illumina', output_file_handle, 'fastq-sanger')
    except ValueError as e:
        # check if this is a quality score problem
        if e.args != ('Invalid character in quality string',):
            raise e
        input_file_handle.seek(0)
        output_file_handle.seek(0)
        output_file_handle.writelines(input_file_handle.readlines())
    finally:
        input_file_handle.close()
        output_file_handle.close()

@follows(mkdir('sai'), mkdir('logs'))
@transform(copy_sequence, regex(r'^fastq/(.+)_sequence\.fastq\.gz$'), r'sai/\1.sai')
def fastq_to_sai(input_file, output_file):
    '''Convert FASTQ files to SAI files.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Aligning sequences', cmd_dict['infile'], cmd_dict['outfile'])
    bwacmd = '%(bwa)s aln -t %(threads)s -f %(outfile)s %(reference)s %(infile)s'
    call(bwacmd, cmd_dict)

# Merge paired ends to SAM
@follows(mkdir('sam'))
@transform(fastq_to_sai, regex(r'^sai/(\w+)_s_(\d)(_1)?\.sai$'),
           inputs([r'sai/\1_s_\2*.sai', r'fastq/\1_s_\2*.fastq.gz']),
           r'sam/\1_s_\2.sam')
def make_sam(input_files, output_file):
    '''Convert SAI files and FASTQ files to SAM files.'''

    def saicmp(x, y):
        '''Compare function for moving sai files to front of list'''
        if x.endswith('sai') and not y.endswith('sai'):
            return - 1
        elif y.endswith('sai') and not x.endswith('sai'):
            return 1
        else:
            return cmp(x, y)

    cmd_dict = CMD_DICT.copy()
    assert type(input_files) is type([])
    pmsg('Generating SAM file', ', '.join(input_files), output_file)
    # sort input files
    input_files.sort(cmp=saicmp)
    # Run bwa to merge into SAM file
    cmd_dict['infiles'] = ' '.join(input_files)
    cmd_dict['outfile'] = output_file
    cmd_dict.update(read_group_re.match(input_files[-1]).groupdict())
    bwa_cmd = '%(bwa)s %(sam_type)s -i %(read_group)s -m %(sample)s -l %(read_group)s ' + \
            '-p illumina -f %(outfile)s %(reference)s %(infiles)s'
    call(bwa_cmd, cmd_dict)

    # if this is single ended, then we need to sort and clip the sam file
    if cmd_dict['sam_type'] == 'samse':
        # sort sam file
        cmd_dict['infile'] = output_file
        cmd_dict['outfile'] = os.path.splitext(output_file)[0] + '.sorted.sam'
        pmsg('Sorting SAM file', cmd_dict['infile'], cmd_dict['outfile'])
        sam_cmd = '%(samtools)s view -hSu %(infile)s | ' + \
                '%(samtools)s sort -o - %(outfile)s | ' + \
                '%(samtools)s view -h -o %(outfile)s -'
        call(sam_cmd, cmd_dict)
        call('rm %(infile)s', cmd_dict, is_logged=False)
        # clip reads
        cmd_dict['infile'] = cmd_dict['outfile']
        cmd_dict['outfile'] = output_file
        pmsg('Clipping reads', cmd_dict['infile'], cmd_dict['outfile'])
        clip_cmd = '%(clip_se_reads)s -o %(outfile)s %(primers)s %(infile)s'
        call(clip_cmd, cmd_dict)
        call('rm %(infile)s', cmd_dict, is_logged=False)

# Sort BAM file
@follows(mkdir('sorted'))
@transform(make_sam, regex(r'^sam/(.+)\.sam'), r'sorted/\1.sorted.bam')
def sort_bam(input_file, output_file):
    '''Sort BAM files by coordinate.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['outprefix'] = os.path.splitext(cmd_dict['outfile'])[0]
    pmsg('Cooridinate sorting BAM file', cmd_dict['infile'], cmd_dict['outfile'])
    sam_cmd = '%(samtools)s view -hSu %(infile)s | %(samtools)s sort - %(outprefix)s'
    call(sam_cmd, cmd_dict)

# Remove duplicates
@follows(mkdir('deduped'))
@transform(sort_bam, regex(r'^sorted/(.+)\.sorted\.bam$'), r'deduped/\1.deduped.bam')
def remove_duplicates(input_file, output_file):
    '''Remove duplicates from BAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    if cmd_dict['sam_type'] == 'sampe':
        cmd_dict['metrics'] = output_file.rstrip('bam') + 'metrics'
        pmsg('Removing duplicates', input_file, output_file)
        sam_cmd = '%(samtools)s rmdup %(infile)s %(outfile)s'
        call(sam_cmd, cmd_dict)
    else:
        pmsg('Linking files', input_file, output_file)
        call('ln %(infile)s %(outfile)s', cmd_dict, is_logged=False)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd, cmd_dict, is_logged=False)

@follows(remove_duplicates, mkdir('coverage'))
@merge('deduped/*.bam', 'coverage/merged.coverage')
def calculate_coverage(input_files, output_file):
    '''Calculate coverage statistics'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infiles'] = ' '.join([ '--input_file %s' % f for f in input_files ])
    cmd_dict['outfile'] = output_file
    pmsg('Coverage calculations', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s --analysis_type DepthOfCoverage ' + \
            '--reference_sequence %(reference)s ' + \
            '--intervals %(target_interval)s ' + \
            '--minMappingQuality 10 ' + \
            '--minBaseQuality 10 ' + \
            '--omitDepthOutputAtEachBase ' + \
            '%(infiles)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)

stages_dict = {
    'copy': copy_sequence,
    'align': fastq_to_sai,
    'sam': make_sam,
    'sort': sort_bam,
    'dedupe': remove_duplicates,
    'coverage': calculate_coverage,
    'default': calculate_coverage,
}

