import logging, subprocess

logger = logging.getLogger('main')

def log_info(func):
    def new_func(*args, **kwargs):
        logger.debug(func.__name__)
        logger.debug('inputs: %s' % str(args[0]))
        logger.debug('outputs: %s' % str(args[1]))
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func

def call(cmdline):
    logger.debug('Calling command line: %s' % cmdline)
    subprocess.call(cmdline, shell=True)
