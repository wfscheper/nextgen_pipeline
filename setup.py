'''
Created on Mar 3, 2011

@author: wfscheper
'''
from setuptools import setup, find_packages

setup(name='nextgen_pipeline',
      version='1.0',
      description='Nextgen Pipeline build with ruffus',
      author='Walter Scheper',
      author_email='scheper@unc.edu',
      url='http://github.com/kmlong/nextgen_pipeline',
      packages=find_packages(exclude=['*.tests']),
      package_data={
          '': ['*.diff']
      },
      install_requires=['ruffus >= 2.2', 'biopython >= 1.55'],
      provides=['nextgen_pipeline (1.0)'],
      entry_points={'console_scripts': ['nextgen_pipeline = nextgen_pipeline.run:run_pipeline']}
      )
