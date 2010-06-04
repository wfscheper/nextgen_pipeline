#!/usr/bin/env python

import datetime
import sys
import multiprocessing
from optparse import OptionParser

from ruffus import *

from steps import *
from log import *

pipeline_stages = {
    'sai': fastq_to_sai,
    'sam': paired_ends_to_sam,
    'bam': sam_to_bam,
    'sort_bam': sort_bam,
    'index_bam': bam_index,
}

def show_pipeline_stage_help():
    print "The pipeline stage you selected does not exist."
    print "Please choose one of the following options:"
    for stage, fn in pipeline_stages.items():
        print '\t%s:\t%s' % (stage, fn.__doc__)
    sys.exit(0)

parser = OptionParser()
parser.add_option('-s', '--stage', dest='stage', help='stage of pipeline to run')
parser.add_option('-t', '--threads', action='store', type='int', dest='threads',
                  default=0, help='number of threads to use')
parser.add_option('-f', '--force', action='store_true', dest='force_run', 
                  default=False, help='Force pipeline to run stage')
parser.add_option('-l', '--log', dest='log', default=None, help='path to log file')

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    
    ncpus = multiprocessing.cpu_count()

    if options.stage:
        if options.stage not in pipeline_stages.keys():
            show_pipeline_stage_help()
        else:
            start_stage = pipeline_stages[options.stage]
    else:
        start_stage = bam_index
    
    if options.threads:
        NUM_JOBS = options.threads if options.threads <= ncpus else ncpus
    else:
        NUM_JOBS = ncpus/2

    if options.log:
        (head, tail) = os.path.split(options.log)
        if os.path.exists(head):
            log_fn = options.log
        else:
            print "Unable to write to that log file."
            sys.exit(0)
    else:
        cwd = os.getcwd()
        now = datetime.datetime.now()
        ts = now.strftime('%Y-%m-%d_%H%M%S')
        log_fn = '%s/%s.%s.log' % (cwd, os.path.split(sys.argv[0])[1], ts)

    logger = quick_start_log(log_fn=log_fn)

    logger.debug('pipeline_run: %d jobs' % NUM_JOBS)

    if options.force_run:
        pipeline_run([start_stage], forcedtorun_tasks=[start_stage], multiprocess=NUM_JOBS, logger=logger)
    else:
        pipeline_run([start_stage], multiprocess=NUM_JOBS, logger=logger)
