import datetime
import multiprocessing
import os
import sys

from ConfigParser import SafeConfigParser
from optparse import OptionParser

from ruffus import pipeline_run, pipeline_printout_graph

from nextgen_pipeline.pipelines import pipelines
from nextgen_pipeline.utils import CMD_DICT
from nextgen_pipeline.utils.log import quick_start_log


PIPELINE_PATH = 'nextgen_pipeline.pipelines.'

def show_pipeline_help(pipelines):
    print("Please choose one of the following pipelines:")
    for pipeline in pipelines:
        print('\t%s' % (pipeline))
    sys.exit(1)

def show_pipeline_stage_help(pipeline_stages):
    print("The pipeline stage you selected does not exist.")
    print("Please choose one of the following stages (* default stage):")
    for pipeline, stages_dict in pipeline_stages.items():
        print '\t%s' % (pipeline)
        for stage, fn in stages_dict.items():
            if stage is not 'default':
                if fn == stages_dict['default']:
                    print('\t\t*%s:\t%s' % (stage, fn.__doc__))
                else:
                    print('\t\t%s:\t%s' % (stage, fn.__doc__))
    sys.exit(1)

def run_pipeline():
    #==============================================================================================
    # Config command line parser
    #==============================================================================================
    parser = OptionParser('nextgen_pipeline [option]... pipeline...')
    parser.add_option('-s', '--stage', dest='stage', help='stage of pipeline to run')
    parser.add_option('-t', '--threads', action='store', type='int', dest='threads',
                      default=None, help='number of threads to use')
    parser.add_option('-j', '--jobs', action='store', type='int', dest='jobs',
                      default=None, help='number of jobs to use')
    parser.add_option('-f', '--force', action='store_true', dest='force_run',
                      default=False, help='Force pipeline to run stage')
    parser.add_option('-c', '--config', dest='config', default='pipeline.cfg',
                      help='config file, defaults to pipeline.cfg')
    parser.add_option('-l', '--log', dest='log', default=None, help='path to log file')
    parser.add_option('--graph', dest='print_graph', default=False, action='store_true',
                      help='Print a graph of the pipeline rather than running it')

    (options, args) = parser.parse_args()
    if len(args) < 1:
        args = pipelines

    # Configuration parsing
    if not os.path.isfile(options.config):
        parser.error('Could not find config file: %s' % options.config)
    config = SafeConfigParser()
    config.read(options.config)

    for section in config.sections():
        for option in config.options(section):
            CMD_DICT[option] = config.get(section, option)

    # get number of cpus on machine
    ncpus = multiprocessing.cpu_count()

    # load the pipeline(s) request by the user
    pipeline_stages = {}
    for arg in args:
        try:
            pipeline = __import__(PIPELINE_PATH + arg, globals(), locals(), ['*'])
        except ImportError, TypeError:
            # either no pipeline was requested or a missing/non-existant
            # pipeline was chosen
            print(parser.usage + '\n')
            show_pipeline_help(pipelines)
        pipeline_stages.update({arg: pipeline.stages_dict})

        # did the user specify a stage
        if options.stage:
            if options.stage not in pipeline_stages[arg].keys():
                # missing or non-existant stage chosen
                show_pipeline_stage_help()
            else:
                start_stage = pipeline_stages[arg][options.stage]
        else:
            # user did not specify a stage, use default
            start_stage = pipeline_stages[arg]['default']

        # user specified job count, capped at the number of cpus
        if options.jobs:
            NUM_JOBS = options.jobs if options.jobs <= ncpus else ncpus
        else:
            NUM_JOBS = config.getint('jobs') \
                    if config.has_option('Threading', 'jobs') else ncpus / 2

        # user specified log file
        if options.log:
            (head, _) = os.path.split(options.log)
            if os.path.exists(head):
                log_fn = options.log
            else:
                print "Unable to write to that log file."
                sys.exit(1)
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
