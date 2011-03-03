'''
Pipeline for calling variants.

Depends on recalibration of quality scores being successfully run, and
generating an indexed, recalibrated bam file.
'''
import os

from glob import glob, iglob
from ruffus import check_if_uptodate, files, follows, inputs, jobs_limit, mkdir, regex, transform

from ..utils import CMD_DICT, call, check_if_clean, pmsg, filename_re, read_group_re


def indel_genoytping_generator():
    '''
    Generates tuples of input files and their resulting output file for indel
    genotyper

    Groups files by sample name and returns a list of those files, and the
    output file they should generate: <sample name>.indels_raw.vcf
    '''
    infiles_by_sample = {}
    for infile in iglob('recalibrated/*.bam'):
        sample = read_group_re.match(infile).groupdict()['sample']
        infiles_by_sample[sample] = infiles_by_sample.get(sample, []) + [infile]
    for sample, infiles in infiles_by_sample.items():
        yield [ infiles, 'indels/%s.indels_raw.vcf' % sample ]

@jobs_limit(1)
@files(['fixmate/', 'intervals/', 'deduped/', 'realigned/'], None)
@check_if_uptodate(check_if_clean)
def clean_up(input_files, output_file):
    '''Clean up intermediate files from recalibration stage'''
    print('Cleaning up intermeidate files: %s' % ', '.join(input_files))
    call('rm -rf %s' % ' '.join(input_files), {}, is_logged=False)

# Call SNPs
@jobs_limit(1)
@follows(clean_up, mkdir('snps', 'logs'))
@files(glob('recalibrated/*.bam'), 'snps/merged.snps_raw.vcf')
def snp_genotyping(input_files, output_file):
    '''Call SNP variants'''
    pmsg('SNP Genotyping', ', '.join(input_files), output_file)
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infiles'] = ' '.join([ '--input_file %s' % f for f in input_files ])
    cmd_dict['outfile'] = output_file
    gatk_cmd = '%(gatk)s ' + \
            '--analysis_type UnifiedGenotyper ' + \
            '--reference_sequence %(reference)s ' + \
            '--DBSNP %(dbsnp)s ' + \
            '--intervals %(target_interval)s ' + \
            '--standard_min_confidence_threshold_for_calling 50 ' + \
            '--standard_min_confidence_threshold_for_emitting 30 ' + \
            '--annotation AlleleBalance ' + \
            '--annotation DepthOfCoverage ' + \
            '--annotation HomopolymerRun ' + \
            '--annotation MappingQualityZero ' + \
            '--annotation QualByDepth ' + \
            '--annotation RMSMappingQuality ' + \
            '--annotation HaplotypeScore ' + \
            '--num_threads 7 ' + \
            '%(infiles)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)

# Call Indels
@follows(snp_genotyping, mkdir('indels'))
@files(indel_genoytping_generator)
def indel_genotyping(input_files, output_file):
    '''Call Indels'''
    pmsg('Indel Genotyping', ', '.join(input_files), output_file)
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infiles'] = ' '.join([ '--input_file %s' % f for f in input_files ])
    cmd_dict['outfile'] = output_file
    gatk_cmd = '%(gatk)s ' + \
            '--analysis_type IndelGenotyperV2 ' + \
            '--reference_sequence %(reference)s ' + \
            '--intervals %(target_interval)s ' + \
            '--window_size %(window_size)s ' + \
            '%(infiles)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)

# Create InDel mask
@follows(indel_genotyping)
@files(glob('indels/*.indels_raw.vcf'), 'indels/indels_mask.bed')
def create_indel_mask(input_files, output_file):
    '''Create Indel Mask'''
    pmsg('Create Indel Mask', ', '.join(input_files), output_file)
    cmd_dict = CMD_DICT.copy()
    cmd_dict['infiles'] = ' '.join([ '--rodBind:%s,vcf %s' % (os.path.basename(f.split('.')[0]), f)
                                    for f in input_files ])
    cmd_dict['outfile'] = output_file
    gatk_cmd = '%(gatk)s ' + \
            '--analysis_type CreateTriggerTrack ' + \
            '--reference_sequence %(reference)s ' + \
            '%(infiles)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)

# Apply basic filters to snp calls
@follows(create_indel_mask)
@files(['snps/merged.snps_raw.vcf', 'indels/indels_mask.bed'], 'snps/merged.snps_filtered.vcf')
def filter_snps(input_files, output_file):
    '''Apply basic filters to SNPs'''
    pmsg('Apply basic filters to SNPs', ', '.join(input_files), output_file)
    cmd_dict = CMD_DICT.copy()
    cmd_dict['snps'], cmd_dict['mask'] = input_files
    cmd_dict['outfile'] = output_file
    gatk_cmd = '%(gatk)s ' + \
            '--analysis_type VariantFiltration ' + \
            '--reference_sequence %(reference)s ' + \
            '--clusterWindowSize 10 ' + \
            '--filterExpression %(hard_to_validate_filter)s ' + \
            '--filterName \"HARD_TO_VALIDATE\" ' + \
            '--rodBind:mask,bed %(mask)s ' + \
            '--maskName InDel ' + \
            '--rodBind:variant,vcf %(snps)s ' + \
            '--out %(outfile)s'
    call(gatk_cmd, cmd_dict)

# Define stages
stages_dict = {
    'clean': clean_up,
    'call_snps': snp_genotyping,
    'call_indels': indel_genotyping,
    'create_mask': create_indel_mask,
    'filter_snps': filter_snps,
    'default': filter_snps,
}
