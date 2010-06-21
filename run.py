#!/usr/bin/env python

import datetime
import os
import sys
import multiprocessing

from glob import iglob as glob
from optparse import OptionParser

from log import quick_start_log
from ruffus import pipeline_run


DEBUG = True
def show_pipeline_stage_help():
    print "The pipeline stage you selected does not exist."
    print "Please choose one of the following options (* default stage):"
    for pipeline, stages_dict in pipeline_stages.items():
        print '\t%s' % (pipeline)
        for stage, fn in stages_dict.items():
            if stage is not 'default':
                if fn == stages_dict['default']:
                    print '\t\t*%s:\t%s' % (stage, fn.__doc__)
                else:
                    print '\t\t%s:\t%s' % (stage, fn.__doc__)
    sys.exit(0)

parser = OptionParser()
parser.add_option('-p', '--pipeline', dest='pipeline', default=None, help='pipeline to run')
parser.add_option('-s', '--stage', dest='stage', help='stage of pipeline to run')
parser.add_option('-t', '--threads', action='store', type='int', dest='threads',
                  default=0, help='number of threads to use')
parser.add_option('-f', '--force', action='store_true', dest='force_run', 
                  default=False, help='Force pipeline to run stage')
parser.add_option('-l', '--log', dest='log', default=None, help='path to log file')

if __name__ == '__main__':
    pipeline_stages = {}
    (options, args) = parser.parse_args()
    
    ncpus = multiprocessing.cpu_count()

    if options.pipeline is None:
        parser.error('Must specify a pipeline')
    else:
        try:
            pipeline = __import__(options.pipeline, globals(), locals(), ['*'])
        except ImportError as e:
            parser.error('Unable to load pipeline named %s' % options.pipeline)
        pipeline_stages.update({options.pipeline: pipeline.stages_dict})

    if options.stage:
        if options.stage not in pipeline_stages[options.pipeline].keys():
            show_pipeline_stage_help()
        else:
            start_stage = pipeline_stages[options.pipeline][options.stage]
    else:
        start_stage = pipeline_stages[options.pipeline]['default']
    
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
        pipeline_run([start_stage], forcedtorun_tasks=[start_stage], multiprocess=NUM_JOBS,
                logger=logger)
    else:
        pipeline_run([start_stage], multiprocess=NUM_JOBS, logger=logger)

