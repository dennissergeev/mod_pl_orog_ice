# -*- coding: utf-8 -*-
"""
Paths to data

Depends on path.py package
"""
import os
from path import Path

# Top-level directory containing code and data
topdir = Path('.').abspath().parent

# Text directory
textdir = topdir/'text'

# PL tracking results
trackdir = topdir/'data'/'tracks'

# Output directories
plotdir = topdir/'figures'

# UM
extdir = Path('/media')/os.getenv('USER')/'Elements'/'phd'
umdatadir = extdir/'modelling'/'UM'/'exp_results'

# Wild cards and path templates
FNAME_MASK = 'umnsa_*'
PATH_MASK = umdatadir/'{fcst_init}'/'{um_res}'/'{idir}'/'{subdir}'
PATH_MASK /= FNAME_MASK