#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pytoxr contains tools for the analysis of data from ToxR experiments

Copyright (C) 2016  Mark George Teese

This software is licensed under the permissive MIT License.
"""

from setuptools import setup, find_packages
from os import path
from codecs import open

# grab the long_description from the readme file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, "readme.rst")) as f:
    long_description = f.read()


classifiers = """\
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
"""

setup(name='pytoxr',
	author="Mark Teese",
    url="https://github.com/teese/pytoxr",
    download_url = 'https://github.com/teese/pytoxr/archive/0.0.7.tar.gz',
    author_email="mark.teese@checkmytumhomepage.de",
    description = "Tools for the analysis of data from ToxR experiments.",
    long_description=long_description,
	long_description_content_type='text/x-rst',
    license='MIT',
    packages=find_packages(),
    classifiers=classifiers.splitlines(),
    keywords="ToxR transmembrane TOXCAT TMDhomodimer GALLEX AraTM BacTH",
	project_urls={'LangoschLab':'http://cbp.wzw.tum.de/index.php?id=9', "TU_Munich":"https://www.tum.de"},
    install_requires=["pandas", "numpy", "matplotlib", "scipy", "seaborn", "eccpy"],
    version='0.0.7')