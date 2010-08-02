'''
Prepares bam files from staged sequence.txt files

Required for any pipelines that use the Genome Analysis Toolkit.

Currnetly, this means quality score recalibration and variant calling.
'''
import os
from Bio import SeqIO
from glob import iglob as glob

from ruffus import follows, files, inputs, jobs_limit, mkdir, regex, transform

from zipper import zip
from utils import call, paired_re, paired_strings, pmsg, read_group_re, saicmp, CMD_DICT


# Copy sequence from staging area
def copy_sequence_generator():
    cwd = os.getcwd()
    for file in glob('../staging_area/*_sequence.txt'):
        filename = unpaired_strings['sequence'] % unpaired_re.match(file).groupdict()
        yield [file, '%s/fastq/%s.fastq.gz' % (cwd, filename.rstrip('.txt'))]

@follows(mkdir('fastq'))
@files(copy_sequence_generator)
@jobs_limit(2)
def copy_sequences(input_file, output_file):
    """Copy sequence files from staging area on thumper1"""
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['outfile_prefix'] = output_file.rstrip('.gz')
    pmsg('Sequence Copy', input_file, cmd_dict['outfile_prefix'])
    SeqIO.convert(cmd_dict['infile'], 'fastq-illumina', cmd_dict['outfile_prefix'], 'fastq-sanger')
    pmsg('Compressing file', cmd_dict['outfile_prefix'], cmd_dict['outfile'])
    zip(cmd_dict['outfile_prefix'])
    os.remove(cmd_dict['outfile_prefix'])

# Convert fastq to sai
def fastq_to_sai_generator():
    cwd = os.getcwd()
    for file in glob('%s/fastq/*' % cwd):
        filename = '%s/sai/%s' % (cwd, unpaired_strings['sai'] % unpaired_re.match(file).groupdict())
        yield [file, filename]

@jobs_limit(2)
@follows(copy_sequences, mkdir('sai'))
#@files(fastq_to_sai_generator)
@transform(copy_sequences, regex(r'^(.*)/fastq/(.*)\.fastq\.gz$'), r'\1/sai/\2.sai')
def fastq_to_sai(input_file, output_file):
    '''Convert FASTQ files to SAI files.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('FASTQ to SAI', cmd_dict['infile'], cmd_dict['outfile'])
    bwacmd = '%(bwa)s aln -t %(threads)s %(genome)s %(infile)s > %(outfile)s'
    call(bwacmd, cmd_dict)

# Merge paired ends to SAM
@follows(fastq_to_sai, mkdir('sam'))
@transform(fastq_to_sai,
            regex(r'^(.*)/sai/(\w+)_s_(\d+)_sequence\.sai$'),
                   r'\1/sai/\2_s_\3_sequence.sai',
                   r'\1/fastq/\2_s_\3_sequence.fastq.gz',
            r'\1/sam/\2_s_\3.sam')
def paired_ends_to_sam(input_files, output_file):
    '''Convert SAI files and FASTQ files to SAM files.'''
    cmd_dict = CMD_DICT.copy()
    assert type(input_files) is type([])
    pmsg('SAM generation', ', '.join(input_files), output_file.strip('.gz'))
    # sort input files
    input_files.sort(cmp=saicmp)
    # Run bwa to merge paired ends into SAM file
    cmd_dict['infiles'] = ' '.join(input_files)
    cmd_dict['outfile'] = output_file.strip('.gz')
    bwa_cmd = '%(bwa)s samse %(genome)s %(infiles)s > %(outfile)s'
    call(bwa_cmd, cmd_dict)

## Convert filtered SAM files to BAM files
@follows(paired_ends_to_sam, mkdir('bam'))
@transform(paired_ends_to_sam, regex(r'^(.*)/sam/(.*)\.sam$'), r'\1/bam/\2.bam')
def sam_to_bam(input_file, output_file):
    '''Convert SAM files to BAM files.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('SAM to BAM', cmd_dict['infile'], cmd_dict['outfile'])
    sam_cmd = '%(samtools)s import %(genome)s.fai %(infile)s %(outfile)s'
    call(sam_cmd, cmd_dict)

# Sort BAM file
@follows(fixmate_bam, mkdir('sorted_bam'))
@transform(fixmate_bam, regex(r'^(.*)/fixmate_bam/(.*)\.fixmate\.bam'), \
        r'\1/sorted_bam/\2.sorted.bam')
def sort_bam(input_file, output_file):
    '''Sort BAM files by coordinate.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['outprefix'] = os.path.splitext(cmd_dict['outfile'])[0]
    pmsg('BAM Coord Sort', cmd_dict['infile'], cmd_dict['outfile'])
    picard_cmd = '%(picard)s SortSam ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'SO=coordinate ' + \
            'MAX_RECORDS_IN_RAM=50000000 ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

@follows(mkdir('clipped_bam'))
@transform(sort_bam, regex(r'^(.*)/sorted_bam/(.*)\.sorted\.bam$'),
        r'\1/clipped_bam/\2.sorted.sam')
def uncompress_bam(input_file, output_file):
    '''Uncompress sorted BAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmesg('Uncompressing BAM', cmd_dict['infile'], cmd_dict['outfile'])
    sam_cmd = '%(samtools)s view -h %(infile)s > %(outfile)s'
    call(sam_cmd, cmd_dict)

@transform(uncompress_bam, regex(r'^(.*)/clipped_bam/(.*).sorted\.sam$'),
        r'\1/clipped_bam/\2.clipped.sam')
def clip_sam(input_file, output_file):
    '''Clip reads in SAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmesg('Clpping reads', cmd_dict['infile'], cmd_dict['outfile'])
    clip_cmd = '%(cliptool)s --out %(outfile)s %(primer_file)s %(infile)s'
    call(clip_cmd, cmd_dict)

# Update header with missing data
@follows(mkdir('prepped_bam'))
@transform(clip_sam, regex(r'^(.*)/clipped_bam/(.*)\.clipped\.sam$'),
        r'\1/prepped_bam/\2.prepped.bam')
def fix_header(input_file, output_file):
    '''Fix header info'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['read_group'] = os.path.split(input_file)[1].rstrip('.clipped.sam')
    cmd_dict.update(read_group_re.match(cmd_dict['read_group']).groupdict())
    cmd_dict['header_tmp'] = cmd_dict['header_tmp'] % cmd_dict
    open(cmd_dict['header_tmp'], 'w').write(
        open(cmd_dict['header_template'], 'r').read() % cmd_dict
    )
    picard_cmd = '%(picard)s ReplaceSamHeader ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s' + \
            'HEADER=%(header_tmp)s ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)
    os.remove(cmd_dict['header_tmp'])

# Create index from BAM - creates BAI files
@follows(fix_header)
@transform(fix_header, regex(r'^(.*)/prepped_bam/(.*)\.bam'), r'\1/prepped_bam/\2.bam.bai')
def bam_index(input_file, output_file):
    '''Index BAM file and create a BAI file.'''
    pmsg('Create BAM Index', input_file, output_file)
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    sam_cmd = '%(samtools)s index %(infile)s'
    call(sam_cmd, cmd_dict)

stages_dict = {
    'copy_sequences': copy_sequences,
    'align_sequences': fastq_to_sai,
    'make_sam': paired_ends_to_sam,
    'make_bam': sam_to_bam,
    'sort_bam': sort_bam,
    'clip_reads': clip_reads,
    'fix_header': fix_header,
    'index_bam': bam_index,
    'default': bam_index,
}

