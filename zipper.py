import logging
import subprocess

logger = logging.getLogger('main')

pigz = '/usr/local/bin/pigz'

def zip(*files, **kwargs):
    default = {
        'threads': 4
    }
    default.update(kwargs)
    for file in files:
        logger.debug('Zipping file: %s ...' % file)
        subprocess.call('%s -f -p %s %s' % (pigz, str(default['threads']), file), shell=True)

def unzip(*files, **kwargs):
    default = {
        'threads': 4
    }
    default.update(kwargs)
    for file in files:
        logger.debug('Unzipping file: %s ...' % file)
        subprocess.call('%s -d -f -p %s %s' % (pigz, str(default['threads']), file), shell=True)
