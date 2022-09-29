.. image:: https://github.com/timcera/hspf_utils/actions/workflows/python-package.yml/badge.svg
    :alt: Tests
    :target: https://github.com/timcera/hspf_utils/actions/workflows/python-package.yml
    :height: 20

.. image:: https://img.shields.io/coveralls/github/timcera/hspf_utils
    :alt: Test Coverage
    :target: https://coveralls.io/r/timcera/hspf_utils?branch=master
    :height: 20

.. image:: https://img.shields.io/pypi/v/hspf_utils.svg
    :alt: Latest release
    :target: https://pypi.python.org/pypi/hspf_utils/
    :height: 20

.. image:: https://img.shields.io/pypi/l/hspf_utils.svg
    :alt: BSD-3 clause license
    :target: https://pypi.python.org/pypi/hspf_utils/
    :height: 20

.. image:: https://img.shields.io/pypi/dd/hspf_utils.svg
    :alt: hspf_utils downloads
    :target: https://pypi.python.org/pypi/hspf_utils/
    :height: 20

.. image:: https://img.shields.io/pypi/pyversions/hspf_utils
    :alt: PyPI - Python Version
    :target: https://pypi.org/project/hspf_utils/
    :height: 20

hspf_utils - Quick Guide
========================
The hspf_utils is a command line script and Python library of utilities to work
with the Hydrological Simulation Program - FORTRAN (HSPF).  Uses pandas
(http://pandas.pydata.org/) or numpy (http://numpy.scipy.org) for any heavy
lifting.

Requirements
------------
* tstoolbox - Time-series toolbox; collected and installed by 'pip' or
  'easy_install' command.
* hspfbintoolbox - Utility to extract time-series from HSFP binary output
  files; collected and installed by 'pip' or 'easy_install' command.

Installation
------------
Should be as easy as running ``pip install hspf_utils`` or
``easy_install hspf_utils`` at any command line.

Usage - Command Line
--------------------
Just run 'hspf_utils --help' to get a list of subcommands::

  usage: hspf_utils [-h] {about,detailed,summary,mapping,parameters} ...

  positional arguments:
    {about,detailed,summary,mapping,parameters}
      about               Display version number and system information.
      detailed            Develops a detailed water balance.
      summary             Develops a summary water balance.
      mapping             Develops a csv file appropriate for joining to a GIS
                          layer.
      parameters          Develops a table of parameter values.

  optional arguments:
    -h, --help            show this help message and exit

For the subcommands that output data it is printed to the screen and you can
then redirect to a file.

Usage - API
-----------
You can use all of the command line subcommands as functions.  The function
signature is identical to the command line subcommands.  The return is always
a PANDAS DataFrame.  Input can be a CSV or TAB separated file, or a PANDAS
DataFrame and is supplied to the function via the 'input_ts' keyword.

Simply import hspf_utils::

  from hspf_utils import hspf_utils

  # Then you could call the functions

  # Once you have a PANDAS DataFrame you can use that as input to other
  # hspf_utils functions.
  ntsd = hspf_utils.aggregate(statistic='mean', agg_interval='daily', input_ts=ntsd)
