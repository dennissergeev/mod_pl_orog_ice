# -*- coding: utf-8 -*-
"""
Common objects for MPLOSI paper
"""
import arke


UM_TIME_FMT = '%Y%m%dT%H%MZ'

um_res = arke.units.GrdStep('km2p2')

hours = dict(STARS72=[15, 27, 39],
             STARS77=[24, 36, 48])
# hours = dict(STARS72=[39],
#              STARS77=[48])
# hours = dict(STARS72=[19, 27],
#              STARS77=[])
pl_cases = dict(
    STARS72='20070404T1200Z',
    STARS77='20080129T1200Z'
)

fcst_strs = [*pl_cases.values()]

runs = ('ctrl', 'nosva', 'sva200', 'ice76n', 'ice82n')

rad = 150