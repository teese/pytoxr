#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pytoxr contains tools for the analysis of data from ToxR experiments

Copyright (C) 2016  Mark George Teese

This software is licensed under the permissive MIT License.
"""

from setuptools import setup, find_packages

classifiers = """\
Development Status :: Experimental
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Classifier: Operating System :: OS Independent
"""

setup(name='pytoxr',
      author='Mark Teese',
      author_email='mark.teese /-at-/ tum /-dot-/ de',
      license='MIT',
      packages=find_packages(),
      classifiers=classifiers.splitlines(),
      platforms=['ALL'],
      keywords=["ToxR", "transmembrane", "TOXCAT", "TMD", "homodimer",
                "GALLEX", "AraTM", "BacTH"],
      version='0.0.2')