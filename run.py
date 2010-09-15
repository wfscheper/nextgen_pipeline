#!/usr/bin/env python
import datetime
import os
import sys
import multiprocessing

from glob import iglob as glob
from optparse import OptionParser
from ruffus import pipeline_run, pipeline_printout_graph

from log import quick_start_log


PIPELINE_PATH = 'pipelines.'

def show_pipeline_help():
    print "Please choose one of the following pipelines:"
    for pipeline in glob('%s/pipelines/*.py' % sys.path[0]):
        if '__init__' in pipeline: continue
        pipeline = os.path.splitext(pipeline.rpartition(os.path.sep)[-1])[0]
        print '\t%s' % (pipeline)
    sys.exit(0)

def show_pipeline_stage_help():
    print "The pipeline stage you selected does not exist."
    print "Please choose one of the following stages (* default stage):"
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
parser.add_option('--graph', dest='print_graph', default=False, action='store_true',
        help='Print a graph of the pipeline rather than running it')

if __name__ == '__main__':
    pipeline_stages = {}
    (options, args) = parser.parse_args()

    # get number of cpus on machine
    ncpus = multiprocessing.cpu_count()

    # load the pipeline(s) request by the user
    if options.pipeline:
        pipelines = [ pipeline.strip() for pipeline in options.pipeline.split(',') ]
    else:
        from pipelines import pipelines

    for pipeline in pipelines:
        try:
            pipeline = __import__(PIPELINE_PATH + pipeline, globals(), locals(), ['*'])
        except (ImportError, TypeError) as e:
            # either no pipeline was requested or a missing/non-existant
            # pipeline was chosen
            show_pipeline_help()
        pipeline_stages.update({options.pipeline: pipeline.stages_dict})

        # did the user specify a stage
        if options.stage:
            if options.stage not in pipeline_stages[options.pipeline].keys():
                # missing or non-existant stage chosen
                show_pipeline_stage_help()
            else:
                start_stage = pipeline_stages[options.pipeline][options.stage]
        else:
            # user did not specify a stage, use default
            start_stage = pipeline_stages[options.pipeline]['default']

        # user specified thread count, capped at the number of cpus
        if options.threads:
            NUM_JOBS = options.threads if options.threads <= ncpus else ncpus
        else:
            NUM_JOBS = ncpus/2

        # user speicified log file
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

        # user said to force running of stage
        if options.print_graph:
            pipeline_printout_graph('pipeline.jpg', 'jpg', [start_stage],
                    forcedtorun_tasks=([start_stage] if options.force_run else []))
        else:
            pipeline_run([start_stage], forcedtorun_tasks=([start_stage] if options.force_run else
                []), multiprocess=NUM_JOBS, logger=logger)

