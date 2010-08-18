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
    for file in glob('./staging_area/*'):
        filename = '%(line)s_s_%(lane)s_%(pair)s_sequence.fastq.gz' % \
                paired_re.search(file).groupdict()
        yield [file, '%s/fastq/%s' % (cwd, filename)]

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

@jobs_limit(2)
@follows(mkdir('sai'))
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
@follows(mkdir('sam'))
@transform(fastq_to_sai, regex(r'^(.*)/sai/(\w+)_s_(\d+)_1_sequence\.sai$'),
           inputs([r'\1/sai/\2_s_\3_2_sequence.sai',
                   r'\1/sai/\2_s_\3_1_sequence.sai',
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

    cmd_dict = CMD_DICT.copy()
    assert type(input_files) is type([])
    pmsg('SAM generation', ', '.join(input_files), output_file.strip('.gz'))
    # sort input files
    input_files.sort(cmp=saicmp)
    # Run bwa to merge paired ends into SAM file
    cmd_dict['infiles'] = ' '.join(input_files)
    cmd_dict['outfile'] = output_file.strip('.gz')
    bwa_cmd = '%(bwa)s sampe %(genome)s %(infiles)s > %(outfile)s'
    call(bwa_cmd, cmd_dict)

## Convert filtered SAM files to BAM files
@follows(mkdir('bam'))
@transform(paired_ends_to_sam, regex(r'^(.*)/sam/(.+)\.sam$'), r'\1/bam/\2.bam')
def sam_to_bam(input_file, output_file):
    '''Convert SAM files to BAM files.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('SAM to BAM', cmd_dict['infile'], cmd_dict['outfile'])
    sam_cmd = '%(samtools)s import %(genome)s.fai %(infile)s %(outfile)s'
    call(sam_cmd, cmd_dict)

# Sort BAM file
@follows(mkdir('sorted'))
@transform(sam_to_bam, regex(r'^(.*)/bam/(.+)\.bam$'), r'\1/sorted/\2.sorted.bam')
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
            'MAX_RECORDS_IN_RAM=7500000 ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

# Remove duplicates
@follows(mkdir('deduped'))
@transform(sort_bam, regex(r'^(.*)/sorted/(.*)\.sorted\.bam$'),
        r'\1/deduped/\2.deduped.bam')
def remove_duplicates(input_file, output_file):
    '''Remove duplicates from BAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['metrics'] = output_file.rstrip('bam') + 'metrics'
    pmsg('Remove duplicates', input_file, output_file)
    picard_cmd = '%(picard)s MarkDuplicates ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'M=%(metrics)s ' + \
            'REMOVE_DUPLICATES=true' + \
            'MAX_RECORDS_IN_RAM=7500000 ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

# Update header with missing data
@follows(mkdir('prepped'))
@transform(remove_duplicates, regex(r'^(.*)/deduped/(.*)\.deduped\.bam$'), r'\1/prepped/\2.prepped.bam')
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
            'MAX_RECORDS_IN_RAM=7500000 ' + \
            'VALIDATION_STRINGENCY=SILENT'

    call(picard_cmd, cmd_dict)
    os.remove(cmd_dict['header_tmp'])

# Create index from BAM - creates BAI files
@transform(fix_header, regex(r'^(.*)/prepped/(.+)\.bam'), r'\1/prepped/\2.bam.bai')
def bam_index(input_file, output_file):
    '''Index BAM file and create a BAI file.'''
    pmsg('Create BAM Index', input_file, output_file)
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    picard_cmd = '%(picard)s BuildBamIndex ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'MAX_RECORDS_IN_RAM=7500000 ' + \
            'OVERWRITE=true ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

@follows(mkdir('coverage'))
@transform(bam_index, regex(r'^(.*)/prepped/(.+)\.prepped\.bam\.bai$'), r'\1/coverage/\2.coverage')
def calculate_coverage(input_file, output_file):
    '''Calculate coverage statistics'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = os.path.splitext(input_file)[0]
    cmd_dict['outfile'] = output_file
    pmsg('Coverage calculations', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s -T DepthOfCoverage ' + \
            '-I %(infile)s ' + \
            '-R %(genome)s ' + \
            '-L %(exome)s ' + \
            '-o %(outfile)s ' + \
            '--minMappingQuality 10 ' + \
            '--minBaseQuality 10 ' + \
            ''#'-omitBaseOutput'
    call(gatk_cmd, cmd_dict)

stages_dict = {
    'copy_sequences': copy_sequences,
    'align_sequences': fastq_to_sai,
    'make_sam': paired_ends_to_sam,
    'make_bam': sam_to_bam,
    'sort_bam': sort_bam,
    'fix_header': fix_header,
    'index_bam': bam_index,
    'calculate_coverage': calculate_coverage,
    'default': calculate_coverage,
}

