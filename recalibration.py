import re
import os
import logging

from glob import glob
from ruffus import *

from zipper import zip, unzip
from utils import log_info, call

STANDARD_FILTER = '\"QUAL < 30.0 || AB > 0.75 && DP > 40 || QD < 5.0 || HRun > 5 || SB > -0.10\"'
HARD_TO_VALIDATE_FILTER = '\"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\"'

cdict = {
    # commands
    'gatk': 'java -Xmx4g -Djava.io.tmpdir=/tmp/$USER -jar /opt/GenomeAnalysisTK/GenomeAnalysisTK.jar',
    'picard': 'java -Xmx4g -jar /opt/picard-tools/',
    'filterIndels': '/usr/local/bin/filterSingleSampleCalls.pl'
    # data files
    'exome': '/usr/local/share/nextgen_resources/hg19_capture.interval_list',
    'genome': '/usr/local/share/nextgen_resources/human_g1k_v37.fasta',
    'dbsnp': '/usr/local/share/nextgen_resources/dbsnp_130_b37.rod',
    # gatk tools
    'CountCovariates': '-T CountCovariates -standard',
    'TableRecalibration': '-T TableRecalibration',
    'RealignerTargetCreator': '-T RealignerTargetCreator',
    'IndelRealigner': '-T IndelRealigner',
    'IndelGenotyper': '-T IndelGenotyperV2 --verbose -mnr 1000000',
    'UnifiedGenotyper': '-T UnifiedGenotyper -stand',
    'VariantFiltration': \
            '-T VariantFiltration ' + \
            '--maskName InDel ' + \
            '-window 10 ' + \
            '-filter ' + STANDARD_FILTER + ' -filterName \"StandardFiler\"' + \
            '-filter ' + HARD_TO_VALIDATE_FILTER  + ' -filterName \"HARD_TO_VALIDATE\"',
}

# Calculate Covariates for Quality Score Recalibration

# Apply Quality Score Recalibration

# Generate QS Recalibration Dataplots

# Find candidate intervals for realignment

# Realign around possible indels

# Call Indels

# Filter Indels

# Call SNPs

# Create InDel mask

# Filter SNPs

stages_dict = {
    'count_covariates': count_covariates,
    'default': count_covariates,
}

