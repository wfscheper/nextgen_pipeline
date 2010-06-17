import re
import os

from glob import glob
from ruffus import files, follows, inputs, mkdir, regex, transform

from utils import CMD_DICT, call, pmsg, unpaired_re, unpaired_strings


# Find candidate intervals for realignment
def create_intervals_generator():
    cwd = os.getcwd()
    for file in glob('%s/recal_bam/*.bam' % cwd):
        filename = '%s/indel_intervals/%s' % (cwd, unpaired_strings['intervals'] % \
                unpaired_re.search(file).groupdict())
        yield [file, filename]

@follows(mkdir('indel_intervals'))
@files(create_intervals_generator)
def create_intervals(input_file, output_file):
    '''Determine indel candidate intervals'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Interval Creation', input_file, output_file)
    gatk_cmd = '%(gatk)s ' + \
            '-T RealignerTargetCreator ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(infile)s ' + \
            '-o %(outfile)s'
    call(gatk_cmd % cmd_dict)

# Realign around possible indels
@follows(create_intervals, mkdir('cleaned_bam'))
@transform(create_intervals, regex(r'^(.*)/indel_intervals/(.+)\.intervals$'),
        inputs([r'\1/recal_bam/\2.recalibrated.bam', r'\1/indel_intervals/\2.intervals']),
        r'\1/cleaned_bam/\2.cleaned.bam')
def local_realignment(input_files, output_file):
    '''Realign reads around candidate indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['recal_bam'] = input_files[0]
    cmd_dict['indel_intervals'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('Local Realignment', ', '.join(input_files), output_file)
    gatk_cmd = '%(gatk)s ' + \
            '-T IndelRealigner ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-I %(recal_bam)s ' + \
            '-targetIntervals %(indel_intervals)s ' + \
            '-mrl 20000 ' + \
            '--output %(outfile)s'
    call(gatk_cmd % cmd_dict)
    samtools_cmd = '%(samtools)s index %(outfile)s'
    call(samtools_cmd % cmd_dict)

# Call Indels
@follows(local_realignment, mkdir('indels'))
@transform(local_realignment, regex(r'^(.*)/cleaned_bam/(.*)\.cleaned\.bam$'),
        r'\1/indels/\2_indels.detailed.bed',
        r'\1/indels/\2_indels.raw.bed')
def indel_genotyping(input_file, details_file, raw_file):
    '''Call Indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['detailsfile'] = details_file
    cmd_dict['rawfile'] = raw_file
    pmsg('Indel Genotyping', input_file, ', '.join([details_file, raw_file]))
    gatk_cmd = '%(gatk)s ' + \
            '-T IndelGenotyperV2 ' + \
            '-R %(genome)s ' + \
            '-I %(infile)s ' + \
            '-O %(rawfile)s ' + \
            '-o %(detailsfile)s ' + \
            '--verbose'
    call(gatk_cmd % cmd_dict)

# Filter Indels
@follows(indel_genotyping, mkdir('filtered_calls'))
@transform(indel_genotyping, regex(r'^(.*)/indels/(.*)\.detailed\.bed$'),
        r'\1/filtered_calls/\2.filtered.bed')
def filter_indels(input_file, output_file):
    '''Filter indel calls'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Filter Indels', input_file, output_file)
    filter_cmd = '%(filterIndels)s ' + \
            '--calls %(infile)s ' + \
            '--max_cons_av_mm 3.0 ' + \
            '--max_cons_nqs_av_mm 0.5 ' + \
            '--mode ANNOTATE ' + \
            '> %(outfile)s'
    call(filter_cmd % cmd_dict)

# Call SNPs
@follows(local_realignment, mkdir('snps'))
@transform(local_realignment, regex(r'^(.*)/cleaned_bam/(.*)\.cleaned\.bam$'),
        r'\1/snps/\2_snps.raw.vcf')
def snp_genotyping(input_file, output_file):
    '''Call SNP variants'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('SNP Genotyping', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s ' + \
            '-T UnifiedGenotyper ' + \
            '-R %(genome)s ' + \
            '-D %(dbsnp)s ' + \
            '-L %(exome)s ' + \
            '-I %(infile)s ' + \
            '-varout %(outfile)s ' + \
            '-stand_call_conf 30.0'
    call(gatk_cmd % cmd_dict)

# Create InDel mask
@follows(indel_genotyping)
@transform(indel_genotyping, regex(r'^(.*)/indels/(.*)\.detailed\.bed$'),
        inputs([r'\1/indels/\2.raw.bed']),
        r'\1/indels/\2.mask.bed')
def create_indel_mask(input_file, output_file):
    '''Create Indel Mask'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file[0]
    cmd_dict['outfile'] = output_file
    pmsg('Indel Masking', cmd_dict['infile'], cmd_dict['outfile'])
    python_cmd = '%(makeIndelMask)s ' + \
            '%(infile)s ' + \
            '10 ' + \
            '%(outfile)s'
    call(python_cmd % cmd_dict)

# Filter SNPs
@follows(snp_genotyping, create_indel_mask, mkdir('filtered_calls'))
@transform(snp_genotyping, regex(r'^(.*)/snps/(.*)_snps.raw.vcf$'),
        inputs([r'\1/snps/\2_snps.raw.vcf',
                r'\1/indels/\2_indels.mask.bed']),
        r'\1/filtered_calls/\2_snps.filtered.vcf')
def filter_snps(input_files, output_file):
    '''Filter SNP Calls'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['snpfile'] = input_files[0]
    cmd_dict['maskfile'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('Filter SNPs', ', '.join(input_files), cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s ' + \
            '-T VariantFiltration ' + \
            '-R %(genome)s ' + \
            '-o %(outfile)s ' + \
            '-B variant,VCF,%(snpfile)s ' + \
            '-B mask,Bed,%(maskfile)s ' + \
            '--maskName InDel ' + \
            '--clusterWindowSize 10 ' + \
            '--filterExpression %(standard_filter)s ' + \
            '--filterName \"StandardFilter\" ' + \
            '--filterExpression %(hard_to_validate_filter)s ' + \
            '--filterName \"HardToValidateFilter\" '
    call(gatk_cmd % cmd_dict)

stages_dict = {
    'create_intervals': create_intervals,
    'local_realignment': local_realignment,
    'indel_genotyping': indel_genotyping,
    'filter_indels': filter_indels,
    'snp_genotyping': snp_genotyping,
    'create_indel_mask': create_indel_mask,
    'filter_snps': filter_snps,
    'default': filter_snps
}

