'''
Pipeline for calling variants.

Depends on recalibration of quality scores being successfully run, and
generating an indexed, recalibrated bam file.
'''
import re
import os

from glob import iglob as glob
from ruffus import check_if_uptodate, files, follows, inputs, jobs_limit, mkdir, regex, transform

from utils import CMD_DICT, call, check_if_clean, pmsg, unpaired_re, unpaired_strings


def indel_genoytping_generator():
    for infile in glob('recalibrated/*.bam'):
        raw_file = '%(line)s_s_%(lane)s.indels_raw.vcf' % \
                unpaired_re.search(infile).groupdict()
        yield [infile, 'indels/%s' % (raw_file)]

def snp_genoytping_generator():
    for infile in glob('recalibrated/*.bam'):
        outfile = '%(line)s_s_%(lane)s.snps_raw.vcf' % \
                unpaired_re.search(infile).groupdict()
        yield [infile, 'snps/%s' % (outfile)]

@jobs_limit(1)
@files(["fixmate", "intervals/", "prepped/", "realigned/"], None)
@check_if_uptodate(check_if_clean)
def clean_up(input_files, output_file):
    print('Cleaning up intermeidate files: %s' % ", ".join(input_files))
    call('rm -rf %s' % " ".join(input_files), {}, is_logged=False)

# Call Indels
@follows(clean_up, mkdir('indels'))
@files(indel_genoytping_generator)
def indel_genotyping(input_file, output_file):
    '''Call Indels'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    cmd_dict['outfile_bed'] = os.path.splitext(output_file)[0] + ".bed"
    pmsg('Indel Genotyping', input_file, output_file)
    gatk_cmd = '%(gatk)s ' + \
            '--analysis_type IndelGenotyperV2 ' + \
            '--reference_sequence %(reference)s ' + \
            '--intervals %(target_interval)s ' + \
            '--input_file %(infile)s ' + \
            '--window_size %(window_size)s ' + \
            '--out %(outfile)s ' + \
            '--bedOutput %(outfile_bed)s'
    call(gatk_cmd, cmd_dict)

# Call SNPs
@jobs_limit(1)
@follows(clean_up, mkdir('snps'))
@files(snp_genoytping_generator)
def snp_genotyping(input_file, output_file):
    '''Call SNP variants'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('SNP Genotyping', cmd_dict['infile'], cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s ' + \
            '--analysis_type UnifiedGenotyper ' + \
            '--reference_sequence %(reference)s ' + \
            '--DBSNP %(dbsnp)s ' + \
            '--intervals %(target_interval)s ' + \
            '--input_file %(infile)s ' + \
            '--out %(outfile)s ' + \
            '--standard_min_confidence_threshold_for_calling 30.0 ' + \
            '--annotation AlleleBalance ' + \
            '--annotation DepthOfCoverage ' + \
            '--annotation HomopolymerRun ' + \
            '--annotation MappingQualityZero ' + \
            '--annotation QualByDepth ' + \
            '--annotation RMSMappingQuality ' + \
            '--annotation HaplotypeScore ' + \
            '--num_threads 8'
    call(gatk_cmd, cmd_dict)

# Filter Indels
@transform(indel_genotyping,
           regex(r'^indels/(.+)\.indels_raw.vcf$'),
           inputs(r'indels/\1.indels_raw.bed'),
           r'indels/\1.indels_filtered.bed')
def filter_indels(input_file, output_file):
    '''Filter indel calls'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Filter Indels', input_file, output_file)
    filter_cmd = '%(filter_indels)s ' + \
            '--calls %(infile)s ' + \
            '--max_cons_av_mm 3.0 ' + \
            '--max_cons_nqs_av_mm 0.5 ' + \
            '--mode ANNOTATE ' + \
            '> %(outfile)s'
    call(filter_cmd, cmd_dict, is_logged=False)

# Create InDel mask
@follows(filter_indels)
@transform(indel_genotyping,
           regex(r'^indels/(.+)\.indels_raw\.vcf$'),
           inputs(r'indels/\1.indels_raw.bed'),
           r'indels/\1.indels_mask.bed')
def create_indel_mask(input_file, output_file):
    '''Create Indel Mask'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infile'] = input_file
    cmd_dict['outfile'] = output_file
    pmsg('Indel Masking', cmd_dict['infile'], cmd_dict['outfile'])
    python_cmd = '%(make_indel_mask)s ' + \
            '%(infile)s ' + \
            '10 ' + \
            '%(outfile)s'
    call(python_cmd, cmd_dict, is_logged=False)

# Filter SNPs
@follows(snp_genotyping, create_indel_mask)
@transform(snp_genotyping, regex(r'^snps/(.+).snps_raw.vcf$'),
        inputs([r'snps/\1.snps_raw.vcf', r'indels/\1.indels_mask.bed']),
        r'snps/\1.snps_filtered.vcf')
def filter_snps(input_files, output_file):
    '''Filter SNP Calls'''
    cmd_dict = CMD_DICT.copy()
    cmd_dict['snpfile'] = input_files[0]
    cmd_dict['maskfile'] = input_files[1]
    cmd_dict['outfile'] = output_file
    pmsg('Filter SNPs', ', '.join(input_files), cmd_dict['outfile'])
    gatk_cmd = '%(gatk)s ' + \
            '--analysis_type VariantFiltration ' + \
            '--reference_sequence %(reference)s ' + \
            '--out %(outfile)s ' + \
            '--rodBind:variant,VCF %(snpfile)s ' + \
            '--rodBind:mask,Bed %(maskfile)s ' + \
            '--maskName InDel ' + \
            '--clusterWindowSize 10 ' + \
            '--filterExpression %(standard_filter)s ' + \
            '--filterName \"StandardFilter\" ' + \
            '--filterExpression %(hard_to_validate_filter)s ' + \
            '--filterName \"HardToValidateFilter\"'
    call(gatk_cmd, cmd_dict)

stages_dict = {
    'clean': clean_up,
    'call_indels': indel_genotyping,
    'filter_indels': filter_indels,
    'call_snps': snp_genotyping,
    'create_mask': create_indel_mask,
    'filter_snps': filter_snps,
    'default': filter_snps
}

