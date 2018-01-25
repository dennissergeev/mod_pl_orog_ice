# -*- coding: utf-8 -*-
import math


def _make_dir(d, force=False):
    if d.exists():
        if force:
            d.rmtree_p()
            d.makedirs()
        # else:
        #     warnings.warn('Directory {d} exists!'.format(d=d))
    else:
        d.makedirs()


def unit_format(value, unit='1'):
    if value == 1:
        if unit == '1':
            string = ''
        else:
            string = f'${str(unit).replace(" ", "$ $")}$'
            if '%' in string:
                string = string.replace("%", "\%")
    else:
        exp = math.floor(math.log10(value))  # use math instead of numpy
        base = value / 10**exp
        if exp == 0 or exp == 1:
            string = r'${0}$'.format(value)
        elif exp == -1:
            string = r'${0:0.1f}$'.format(value)
        else:
            if int(base) == 1:
                string = r'$10^{{{0:d}}}$'.format(int(exp))
            else:
                string = r'${0:d}\times10^{{{1:d}}}$'.format(int(base),
                                                             int(exp))
        if not unit == '1':
            string += r' ${}$'.format(unit)
    return string
