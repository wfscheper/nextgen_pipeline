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


def copy_sequence_generator():
    for file in glob('staging_area/*'):
        filename = '%(line)s_s_%(lane)s_%(pair)s.fastq.gz' % \
                paired_re.search(file).groupdict()
        yield [file, 'fastq/%s' % (filename)]

# Copy sequence from staging area
@follows(mkdir('fastq'))
@files(copy_sequence_generator)
@jobs_limit(2)
def copy_sequence(input_file, output_file):
    '''Copy sequence files from staging area'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['outfile_prefix'] = output_file.rstrip('.gz')
    pmsg('Copying sequence files', input_file, cmd_dict['outfile_prefix'])
    try:
        SeqIO.convert(cmd_dict['infile'], 'fastq-illumina', cmd_dict['outfile_prefix'], 'fastq-sanger')
    except ValueError:
        call('cp %(infile)s %(outfile_prefix)s', cmd_dict, is_logged=False)
    pmsg('Compressing file', cmd_dict['outfile_prefix'], cmd_dict['outfile'])
    zip(cmd_dict['outfile_prefix'])

@follows(mkdir('sai'), mkdir('logs'))
@transform(copy_sequence, regex(r'^fastq/(.+)\.fastq\.gz$'), r'sai/\1.sai')
def fastq_to_sai(input_file, output_file):
    '''Convert FASTQ files to SAI files.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    bwacmd = '%(bwa)s aln -t %(threads)s %(genome)s %(infile)s > %(outfile)s'
    pmsg('Aligning sequences', cmd_dict['infile'], cmd_dict['outfile'])
    call(bwacmd, cmd_dict)

# Merge paired ends to SAM
@follows(mkdir('sam'))
@transform(fastq_to_sai, regex(r'^sai/(\w+)_s_(\d+)_1\.sai$'),
           inputs([r'sai/\1_s_\2_2.sai',
                   r'sai/\1_s_\2_1.sai',
                   r'fastq/\1_s_\2_1.fastq.gz',
                   r'fastq/\1_s_\2_2.fastq.gz']),
            r'sam/\1_s_\2.sam')
def paired_ends_to_sam(input_files, output_file):
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
    # Run bwa to merge paired ends into SAM file
    cmd_dict['infiles'] = ' '.join(input_files)
    cmd_dict['outfile'] = output_file
    bwa_cmd = '%(bwa)s sampe -f %(outfile)s %(genome)s %(infiles)s'
    call(bwa_cmd, cmd_dict)

## Convert filtered SAM files to BAM files
@follows(mkdir('bam'))
@transform(paired_ends_to_sam, regex(r'^sam/(.+)\.sam$'), r'bam/\1.bam')
def sam_to_bam(input_file, output_file):
    '''Convert SAM files to BAM files.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Converting SAM file to BAM', cmd_dict['infile'], cmd_dict['outfile'])
    sam_cmd = '%(samtools)s import %(genome)s.fai %(infile)s %(outfile)s'
    call(sam_cmd, cmd_dict)

# Sort BAM file
@follows(mkdir('sorted'))
@transform(sam_to_bam, regex(r'^bam/(.+)\.bam$'), r'sorted/\1.sorted.bam')
def sort_bam(input_file, output_file):
    '''Sort BAM files by coordinate.'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['outprefix'] = os.path.splitext(cmd_dict['outfile'])[0]
    pmsg('Cooridinate sorting BAM file', cmd_dict['infile'], cmd_dict['outfile'])
    picard_cmd = '%(picard)s SortSam ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'SO=coordinate ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

# Remove duplicates
@follows(mkdir('deduped'))
@transform(sort_bam, regex(r'^sorted/(.+)\.sorted\.bam$'), r'deduped/\1.deduped.bam')
def remove_duplicates(input_file, output_file):
    '''Remove duplicates from BAM file'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['metrics'] = output_file.rstrip('bam') + 'metrics'
    pmsg('Removing duplicates', input_file, output_file)
    picard_cmd = '%(picard)s MarkDuplicates ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'M=%(metrics)s ' + \
            'REMOVE_DUPLICATES=true ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)

# Update header with missing data
@follows(mkdir('header'))
@transform(remove_duplicates, regex(r'^deduped/(.+)\.deduped\.bam$'), r'header/\1.header.bam')
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
    pmsg('Fixing header', input_file, output_file)
    picard_cmd = '%(picard)s ReplaceSamHeader ' + \
            'I=%(infile)s ' + \
            'O=%(outfile)s ' + \
            'HEADER=%(header_tmp)s ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(picard_cmd, cmd_dict)
    os.remove(cmd_dict['header_tmp'])

@follows(mkdir('prepped'))
@transform(fix_header, regex('^header/(.+)\.header\.bam$'), r'prepped/\1.prepped.bam')
def tag_reads(input_file, output_file):
    '''Tag reads with read group'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['samfile'] = os.path.splitext(output_file)[0] + '.sam'
    cmd_dict.update(read_group_re.match(input_file).groupdict())
    pmsg('Tagging reads', input_file, output_file)
    command = '%(samtools)s view -h %(infile)s | ' + \
            'sed "/^[^@]/s|$|\tPG:Z:bwa\tRG:Z:%(read_group)s|" > %(samfile)s'
    call(command, cmd_dict, is_logged=False)
    pmsg('Compressing SAM file', cmd_dict['samfile'], output_file)
    command = '%(picard)s SamFormatConverter ' + \
            'I=%(samfile)s ' + \
            'O=%(outfile)s ' + \
            'CREATE_INDEX=true ' + \
            'VALIDATION_STRINGENCY=SILENT'
    call(command, cmd_dict)
    os.remove(cmd_dict['samfile'])

@follows(mkdir('coverage'))
@transform(tag_reads, regex(r'^prepped/(.+)\.prepped\.bam$'), r'coverage/\1.coverage')
def calculate_coverage(input_file, output_file):
    '''Calculate coverage statistics'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
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
    'copy': copy_sequence,
    'align': fastq_to_sai,
    'sam': paired_ends_to_sam,
    'bam': sam_to_bam,
    'sort': sort_bam,
    'dedupe': remove_duplicates,
    'fix_header': fix_header,
    'tag': tag_reads,
    'coverage': calculate_coverage,
    'default': calculate_coverage,
}

