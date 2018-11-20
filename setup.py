#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

__doc__ = "Detecting poly-adenylation tails and primer sequences in next-generation sequencing reads"
__version__ = '0.1'

# import detectpolya

setup(name='detectpolya',
    version = __version__,
    description='Detecting poly-adenylation tails and primer sequences in next-generation sequencing reads',
    author='Olivier M. F. Martin',
    author_email='oliviermfmartin@gmail.com',
    url='https://github.com/oliviermfmartin/detectpolya',
    packages=['detectpolya'],
 	install_requires = ['progressbar','collections', 'itertools', 'numpy', 'HTSeq', 'fuzzysearch', 'Bio.Seq', 'Bio', 'warnings']
 	)
