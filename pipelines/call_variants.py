'''
Pipeline for calling variants.

Depends on recalibration of quality scores being successfully run, and
generating an indexed, recalibrated bam file.
'''
import re
import os

from glob import iglob as glob
from ruffus import files, follows, inputs, mkdir, regex, transform

from utils import CMD_DICT, call, pmsg, unpaired_re, unpaired_strings


def indel_genoytping_generator():
    cwd = os.getcwd()
    for infile in glob('%s/realigned_fixmate/*.bam' % cwd):
        details_file = '%(line)s_s_%(lane)s_indels.detailed.bed' % \
                unpaired_re.search(infile).groupdict()
        raw_file = '%(line)s_s_%(lane)s_indels.raw.bed' % \
                unpaired_re.search(infile).groupdict()
        yield [infile, '%s/indels/%s' % (cwd, details_file), '%s/indels/%s' % (cwd, raw_file)]

# Call Indels
@follows(mkdir('indels'))
@files(indel_genoytping_generator)
def indel_genotyping(input_file, details_file, raw_file):
    '''Call Indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['detailsfile'] = details_file
    cmd_dict['outfile'] = raw_file
    pmsg('Indel Genotyping', input_file, ', '.join([details_file, raw_file]))
    gatk_cmd = '%(gatk)s ' + \
            '-T IndelGenotyperV2 ' + \
            '-R %(genome)s ' + \
            '-I %(infile)s ' + \
            '-O %(outfile)s ' + \
            '-o %(detailsfile)s ' + \
            '--verbose'
    call(gatk_cmd, cmd_dict)

# Call SNPs
def snp_genoytping_generator():
    cwd = os.getcwd()
    for infile in glob('%s/realigned_fixmate/*.bam' % cwd):
        outfile = '%(line)s_s_%(lane)s_snps.raw.vcf' % \
                unpaired_re.search(infile).groupdict()
        yield [infile, '%s/snps/%s' % (cwd, outfile)]

@follows(mkdir('snps'))
@files(snp_genoytping_generator)
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
    call(gatk_cmd, cmd_dict)

# Filter Indels
@transform(indel_genotyping, regex(r'^(.*)/indels/(.+)\.detailed\.bed$'),
        inputs([r'\1/indels/\2.detailed.bed']),
        r'\1/indels/\2.filtered.bed')
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
    call(filter_cmd, cmd_dict)

# Create InDel mask
@follows(filter_indels)
@transform(indel_genotyping, regex(r'^(.*)/indels/(.*)\.detailed\.bed$'),
        inputs([r'\1/indels/\2.raw.bed']), r'\1/indels/\2.mask.bed')
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
    call(python_cmd, cmd_dict)

# Filter SNPs
@follows(snp_genotyping, create_indel_mask)
@transform(snp_genotyping, regex(r'^(.*)/snps/(.*)_snps.raw.vcf$'),
        inputs([r'\1/snps/\2_snps.raw.vcf',
                r'\1/indels/\2_indels.mask.bed']),
        r'\1/snps/\2_snps.filtered.vcf')
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
    call(gatk_cmd, cmd_dict)

stages_dict = {
    'indel_genotyping': indel_genotyping,
    'filter_indels': filter_indels,
    'snp_genotyping': snp_genotyping,
    'create_indel_mask': create_indel_mask,
    'filter_snps': filter_snps,
    'default': filter_snps
}

