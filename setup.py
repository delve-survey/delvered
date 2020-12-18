#!/usr/bin/env python

from distutils.core import setup

setup(name='delvered',
      version='1.0',
      description='David Nidever DELVE reduction code',
      author='David Nidever',
      author_email='dnidever@montana.edu',
      url='https://github.com/dnidever/delvered',
      package_dir = {'': 'python'},
      packages=['delvered'],
      scripts=['bin/make_delvered_summary_table','bin/query_delvered_summary_table','bin/delve_archive_search','bin/parse_archive_search'],
      requires=['numpy','astropy','scipy','dlnpyutils']
)
