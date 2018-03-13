#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot MetUM output from 1 or more experiments on the same figure
  - horizontal cross-sections (maps)
"""
# Standard packages
import argparse
# import cartopy.crs as ccrs
from datetime import timedelta, datetime
import iris
import daiquiri, logging  # NOQA
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
from path import Path
import pandas as pd
import re
# from scipy.ndimage.filters import gaussian_filter
import sys
# My packages
import arke
from arke.io import extract_as_single_cube
from arke.plot import prepare_map
# Local files
# from cart_func import lcc_map_grid  # , best_ticks
import misc_utils as misc
import mypaths
from um_cube_params import CUBE_PARAMS, vert_coord_by_units
from um_cube_diags import compute

# Settings for saving figures
fmt = 'png'
svfigkw = dict(format=fmt, dpi=100, bbox_inches='tight')
imgname_mask = '{vrbls}_{subdir}_{time}{track}.{fmt}'
# Dates and times
# DEFAULT_FCST = '20070120T1200Z'  # '20130324_1200'
# DATETIME_START = '2007-01-20-12:00'  # '2013-03-24-12:00'
# DATETIME_END = '2007-01-22-12:00'  # '2013-03-26-12:00'
DEFAULT_START = '0h'
DEFAULT_NTIMES = '48h'
DEFAULT_FREQ = '1H'
DEFAULT_GRD_DELTA = 'km2p2'
# File path masks
PLOTDIR = Path(mypaths.plotdir)
FNAME_MASK = 'umnsa_*'  # 'um*rholev_*-proc.nc'
TOPDIR = Path(mypaths.umdatadir)
PATH_MASK = TOPDIR / '{fcst_init}' / '{um_res}' / '{idir}' / '{subdir}'
PATH_MASK /= FNAME_MASK
LOGPATH = Path(__file__).dirname() / 'logs'
SCRIPT = Path(__file__).basename().splitext()[0]


def parse_args(args=None):
    ap = argparse.ArgumentParser(SCRIPT,
                                 description=__doc__,
                                 formatter_class=argparse.
                                 ArgumentDefaultsHelpFormatter,
                                 epilog='Example of use: coming soon')

    ap.add_argument('-v', '--verbose', action='count',
                    default=0, help='Verbosity (-v, -vv, etc)')

    ap_exp = ap.add_argument_group(title='Experiments and variables')
    ap_exp.add_argument('-n', '--names', type=str,
                        required=True,
                        help=('Variables @ levels separated by commas'
                              '\n(e.g. vort@950hPa,temp@850hPa,slp)'))
    ap_exp.add_argument('-r', '--runs', type=str,
                        required=True,
                        help='Experiment names, comma-sep (e.g. ctrl,nosva)')
    ap_exp.add_argument('--umres', type=str, nargs='+',
                        default=[DEFAULT_GRD_DELTA],
                        help='MetUM grid spacing')
    ap_exp.add_argument('--fcst', required=True, nargs='+',
                        type=misc.valid_date,
                        help='Forecast reference time (YYYYmmddTHHMMZ)')

    ap_time = ap.add_argument_group('Time coordinate selection')
    ap_time.add_argument('--start', type=str,
                         default=DEFAULT_START,
                         help='Start period')
    ap_time.add_argument('--ntimes', type=str,
                         default=DEFAULT_NTIMES,
                         help='Number of time slices')
    ap_time.add_argument('--freq', type=str,
                         default=DEFAULT_FREQ,
                         help='Time frequency')

    ag_sub = ap.add_argument_group(title='Subset')
    ag_tr = ag_sub.add_argument_group(title='Subset along a trajectory')
    ag_tr.add_argument('--track_fname', type=str,
                       help='Name of text file with coordinates')

    ag_sub = ag_sub.add_mutually_exclusive_group()
    ag_sub.add_argument('-ra', '--radius', type=float,
                        help='Subset radius (km)')
    ag_sub.add_argument('-wh', '--width_height', type=str,
                        help='Subset width and height in grid points (w,h)')
    ag_sub.add_argument('--rim', type=str,
                        help='Cut rim zone around domain in grid points (w)')
    ag_sub.add_argument('-xy', '--xy', type=str,
                        help='Subset output data by indices (x0,x1,y0,y1)')
    ag_sub.add_argument('-ll', '--lonlat', type=str,
                        help=('Subset output data by longitude and latitude'
                              '(lon0,lon1,lat0,lat1)'))

    ap_plt = ap.add_argument_group(title='Plotting settings')
    ap_plt.add_argument('-g', '--geoax', action='store_true', default=False,
                        help='Plot results on the map')
    ap_plt.add_argument('-f', '--force', action='store_true',
                        default=False,
                        help='Overwrite output directory')
    ap_plt.add_argument('--gif', action='store_true', default=False,
                        help='Merge the saved figures into a .gif animation')

    return ap.parse_args(args)


# @profile
def prepare_data(cubelist_in, varnames,
                 idt, level_dict, h_subset):
    # if extract_before_diag:
    #    cl = arke.subset.extract_levels(cubelist, level_dict)
    if len(h_subset) > 0:
        cubelist = arke.subset.subset_cubelist(cubelist_in, h_subset)
    else:
        cubelist = cubelist_in
    # Preserve useful attributes
    add_attr = {k: cubelist[0].attributes[k]
                for k in ('run', 'um_res', 'subdir')}

    vrbls = []
    for name in varnames:
        ld = level_dict[name]
        if (CUBE_PARAMS[name].get('extract_level_before_diag', False)
           and 'value' in ld):
            cl = arke.subset.slice_cubelist(cubelist, ld['name'], ld['value'])
        else:
            cl = cubelist
        # print(f'\nName: {name}\n LD: {ld}\n\n')
        subdir_cnstr = iris.AttributeConstraint(subdir=ld['subdir'])
        if 'stash' in CUBE_PARAMS[name]:
            cnstr = iris.AttributeConstraint(STASH=CUBE_PARAMS[name]['stash'])
        else:
            cnstr = CUBE_PARAMS[name]['_cf_name']
        diag_method = CUBE_PARAMS[name].get('diag_method')
        if diag_method is None:
            # print(name)
            # print('AAA', cl)
            # print(cl.extract(cnstr))
            # print(idt, cl[0].attributes['subdir'], cl[0].attributes['run'])
            # _tmp = cl.extract(cnstr)
            # print('\nBy name:', _tmp)
            # print(f'{_tmp[0].attributes}\n')
            # print(f'By {subdir_cnstr}:', cl.extract(subdir_cnstr))
            # cube = cl.extract(cnstr)[0]
            # print('CCC', cl)
            # logger.debug(cl)
            try:
                cube = cl.extract_strict(cnstr)
            except:
                cube = extract_as_single_cube(cl, cnstr & subdir_cnstr)
            # print('DDD', cube)
            cube = cube.extract(
                   iris.Constraint(time=arke.coords.nearest_tval(cube, idt)))
        else:
            cube = compute(cl, diag_method, date_time=idt,
                           cnstr=subdir_cnstr,
                           )

        if isinstance(cube, iris.cube.Cube):
            if cube.ndim > 2:
                try:
                    zcoord = cube.coord(ld['name'])
                    idx = zcoord.nearest_neighbour_index(ld['value'])
                    constr = {zcoord.name(): zcoord.points[idx]}
                    cube = cube.extract(iris.Constraint(**constr))
                except iris.exceptions.CoordinateNotFoundError:
                    msg = 'Cube "{cube}" does not have z coordinate;'\
                          'Cannot extract by {name}={value}'
                    raise Exception(msg.format(cube=cube.name(), **ld))
            if cube.name() == 'unknown':
                cube.rename(CUBE_PARAMS[name].get('_cf_name'))
                cube.units = CUBE_PARAMS[name].get('_cf_units')
            cube.attributes.update(CUBE_PARAMS[name])
            cube.attributes.update(add_attr)  # bc it is empty after compute()
        else:
            for _cube in cube:
                if _cube.ndim > 2:
                    try:
                        zcoord = _cube.coord(ld['name'])
                        idx = zcoord.nearest_neighbour_index(ld['value'])
                        constr = {zcoord.name(): zcoord.points[idx]}
                        _cube = _cube.extract(iris.Constraint(**constr))
                    except iris.exceptions.CoordinateNotFoundError:
                        msg = 'Cube "{cube}" does not have z coordinate;'\
                              'Cannot extract by {name}={value}'
                        raise Exception(msg.format(cube=_cube.name(), **ld))

                if _cube.name() == 'unknown':
                    _cube.rename(CUBE_PARAMS[name].get('_cf_name'))
                    _cube.units = CUBE_PARAMS[name].get('_cf_units')
                _cube.attributes.update(CUBE_PARAMS[name])
                _cube.attributes.update(add_attr)
        vrbls.append(cube)

    return vrbls


def plot_on_same_fig(fig, axgr, mapkey, vrbls_lists, geoax=False):
    cax = axgr.cbar_axes[0]
    cbar_ready = False
    _dt = None
    for count, (ax, vrbls) in enumerate(zip(axgr.axes_all, vrbls_lists)):
        ttl = []
        for icube in vrbls:
            if isinstance(icube, (list, tuple)):
                if len(icube) == 2:
                    # Assume we are dealing with u, v winds
                    u, v = icube
                    scl = u.attributes.get('scl', 1.)
                    method = getattr(ax, u.attributes['plt_method'])
                    if geoax:
                        x, y = arke.grid.unrotate_lonlat_grids(u)
                    else:
                        # x, y = iris.analysis.cartography.get_xy_grids(u)
                        dx = u.attributes['um_res'].to_flt('km')
                        x = np.arange(u.shape[1]) * dx
                        y = np.arange(u.shape[0]) * dx
                        ax.set_xlabel('Distance along x-axis, km', fontsize=10)
                        ax.set_ylabel('Distance along y-axis, km', fontsize=10)
                        ax.set_xlim(x.min(), x.max())
                        ax.set_ylim(y.min(), y.max())
                    kw = u.attributes.get('plt_kw', {})
                    if u.attributes['plt_method'] == 'streamplot':
                        ws = (u.data**2 + v.data**2)**0.5
                        kw.update(linewidth=ws/ws.max())
                    if u.attributes['plt_method'] == 'quiver':
                        xydim = len(x.shape)
                        skip = [slice(None, None, 10)]
                        p = method(x[skip*xydim], y[skip*xydim],
                                   u.data[skip*2]*scl,
                                   v.data[skip*2]*scl,
                                   **mapkey, **kw)
                        ax.quiverkey(p, 0.7, 0.9, 5, r'$5 \frac{m}{s}$',
                                     labelpos='E', coordinates='figure')
                    else:
                        p = method(x, y, u.data*scl, v.data*scl,
                                   **mapkey, **kw)
                    if _dt is None:
                        _dt = arke.coords.get_cube_datetimes(u)[0]
                    ttl.append(u.attributes.get('_cf_name'))
                # ttl.append('\nand '.join([i.name() for i in icube]) +
                #            ' ({scl})'.format(scl=misc.scl2latex(scl)))
                else:
                    raise NotImplementedError()
            else:
                scl = icube.attributes.get('scl', 1.)
                scl_str = misc.unit_format(scl**(-1), icube.units)
                scl_str = scl_str.replace('hertz', 's^{-1}')
                method = getattr(ax, icube.attributes['plt_method'],
                                 'contourf')
                ma_out = icube.attributes.get('ma_out')
                if ma_out:
                    icube.data = np.ma.masked_outside(icube.data, *ma_out)
                if geoax:
                    x, y = arke.grid.unrotate_xy_grids(icube)
                else:
                    # x, y = iris.analysis.cartography.get_xy_grids(icube)
                    dx = icube.attributes['um_res'].to_flt('km')
                    x = np.arange(icube.shape[1]) * dx
                    y = np.arange(icube.shape[0]) * dx
                    ax.set_xlabel('Distance along x-axis, km', fontsize=10)
                    ax.set_ylabel('Distance along y-axis, km', fontsize=10)
                    ax.set_xlim(x.min(), x.max())
                    ax.set_ylim(y.min(), y.max())
                kw = icube.attributes.get('plt_kw', {})
                p = method(x, y, icube.data*scl, **mapkey, **kw)

                cbar = icube.attributes.get('colorbar')
                if isinstance(cbar, dict) and not cbar_ready:
                    cb = fig.colorbar(p, cax=cax, **cbar)
                    cb.ax.set_title(scl_str, fontsize=10)
                    cbar_ready = True
                    ttl.append(icube.attributes["_cf_name"])
                else:
                    _ttl = icube.attributes["_cf_name"]
                    if scl_str:
                        _ttl += f' ({scl_str})'
                    ttl.append(_ttl)

                clab = icube.attributes.get('clabels')
                if isinstance(clab, dict):
                    ax.clabel(p, **clab)
                if _dt is None:
                    _dt = arke.coords.get_cube_datetimes(icube)[0]
        if count == 0:
            txt = '\n'.join(ttl)
            at = AnchoredText(txt, prop=dict(size=12), frameon=False,
                              loc=3, bbox_to_anchor=(0., 1.),
                              bbox_transform=ax.transAxes)
            ax.add_artist(at)
        # ax.set_title('\n'.join(ttl))
        txt = f'{_dt:%b %d, %H%M}UTC'
        at = AnchoredText(txt, prop=dict(size=10), frameon=True, loc=2)
        ax.add_artist(at)
        exp_label = icube.attributes.get('run', 'test')
        txt = f'{exp_label}'
        at = AnchoredText(txt, prop=dict(size=10), frameon=True, loc=1)
        ax.add_artist(at)
        # plt.axis('off')
        # cb.remove()


def main(args=None):
    """ Main entry point of the script """
    args = parse_args(args)
    # Logging set up
    misc._make_dir(LOGPATH)
    LOGFILE = LOGPATH / '{}_{:%Y%m%d%H%M}.log'.format(SCRIPT, datetime.now())
    daiquiri.setup(
                   level=logging.DEBUG,
                   outputs=(daiquiri.output.File(LOGFILE),
                            daiquiri.output.Stream(sys.stdout))
                   )
    logger = daiquiri.getLogger(SCRIPT)
    logger.setLevel((5 - min(args.verbose, 4))*10)

    level_dict = dict()
    level_regex = re.compile('([0-9]+)([a-zA-Z]+)')
    lev_str = ''
    var_str = ''
    lev_flag = 0
    subdirs = set()
    varnames = []
    for name in args.names.split(','):
        if '@' in name:
            vrbl, level = name.split('@')
            varnames.append(vrbl)
            level_value, level_unit = level_regex.match(level).groups()
            level_dict[vrbl] = dict(value=float(level_value), unit=level_unit)
            lev_str += '{vrbl}{value:.0f}{unit}_'.format(vrbl=vrbl,
                                                         **level_dict[vrbl])
            var_str += '{vrbl}_'.format(vrbl=vrbl)
            lev_flag = 1
            level_dict[vrbl].update(vert_coord_by_units(level_unit))
            subdirs.add(level_dict[vrbl]['subdir'])
        else:
            varnames.append(name)
            level_dict[name] = vert_coord_by_units('')
            subdirs.add(level_dict[name]['subdir'])
            lev_str += '{name}_'.format(name=name)
            var_str += '{name}_'.format(name=name)

    logger.info(level_dict)

    um_res_list = [arke.units.GrdStep(i) for i in args.umres]
    fcst_strs = [i.strftime('%Y%m%dT%H%MZ') for i in args.fcst]
    runs = args.runs.split(',')
    all_exp = fcst_strs + [str(i) for i in um_res_list] + runs
    lab_mask = ''
    if len(fcst_strs) > 1:
        lab_mask += '{fcst} '
    if len(um_res_list) > 1:
        lab_mask += '{um_res} '
    if len(runs) > 1:
        lab_mask += '{run}'

    input_paths = []
    list_of_cubelists = []
    exp_labels = []
    for fcst in fcst_strs:
        for um_res in um_res_list:
            for run in runs:
                _joint_cubelist = iris.cube.CubeList()
                for subdir in subdirs:
                    input_path = PATH_MASK.format(um_res=um_res,
                                                  fcst_init=fcst,
                                                  idir=run,
                                                  subdir=subdir)
                    cubelist = iris.load(input_path,
                                         callback=arke.io.clean_call)
                    for cube in cubelist:
                        add_attr = dict(subdir=subdir,
                                        run=run,
                                        um_res=um_res)
                        cube.attributes.update(add_attr)
                    # print(cubelist[0].attributes)
                    logger.debug('{} loaded'.format(input_path))
                    exp_labels.append(lab_mask.format(fcst=fcst,
                                                      um_res=um_res,
                                                      run=run))
                    input_paths.append(input_path)
                    _joint_cubelist += cubelist
                list_of_cubelists.append(_joint_cubelist)
    logger.info(exp_labels)
    # for a in list_of_cubelists:
    #     print(a)

    if args.track_fname:
        tr_fname = args.track_fname
        try:
            tr_df = pd.read_csv(tr_fname, parse_dates=[0])
        except FileNotFoundError:
            tr_fname = ((Path(input_path).parent.parent/'track')
                        .glob(f'pl_loc*{tr_fname}*')[0])
            tr_df = pd.read_csv(tr_fname, parse_dates=[0])
        track_time = tr_df.time.values.astype('<M8[ms]').astype(datetime)

    subset = {}
    if args.track_fname and (args.xy or args.lonlat):
        raise argparse.ArgumentTypeError(('Use --radius/--width_height '
                                          'when subsetting along a track'))
    if args.radius:
        subset['method'] = 'radius'
        subset['r'] = args.radius
    if args.rim:
        subset['method'] = 'rim'
        subset['width'] = int(args.rim)
    if args.xy:
        subset['method'] = 'xy'
        subset['corners'] = [int(i) for i in args.xy.split(',')]
    if args.lonlat:
        subset['method'] = 'lonlat'
        subset['corners'] = [float(i) for i in args.lonlat.split(',')]
    if args.width_height:
        subset['method'] = 'wh'
        subset['w'], subset['h'] = [float(i)
                                    for i in args.width_height.split(',')]

    if args.track_fname:
        time_range = track_time
    else:
        if len(args.fcst) == 1:
            if args.start[-1] == 'h':
                tdelta_start = dict(hours=int(args.start[:-1]))
            if args.ntimes[-1] == 'h':
                tdelta_ntimes = dict(hours=int(args.ntimes[:-1]))
            t_start = args.fcst[0] + timedelta(**tdelta_start)
            t_end = args.fcst[0] + timedelta(**tdelta_ntimes)
            time_range = pd.date_range(start=t_start,
                                       end=t_end,
                                       freq=args.freq).to_pydatetime()
        else:
            # TODO: make more flexible instead of
            # using the first cube in each cubelist for now
            input_datetimes = [set(arke.coords.get_cube_datetimes(i[0]))
                               for i in list_of_cubelists]
            time_range = sorted(set.intersection(*input_datetimes))
            time_range = pd.date_range(start=time_range[0],
                                       end=time_range[-1],
                                       freq=args.freq).to_pydatetime()

    subdir_out = 'cf_{}'.format('_'.join(all_exp))
    subsubdir = var_str[:-1]
    if args.track_fname:
        subsubdir += '_track'
        track_str = '_track'
    else:
        track_str = ''
    output_dir = PLOTDIR / subdir_out / subsubdir
    if lev_flag == 1:
        output_dir /= lev_str[:-1]
    misc._make_dir(output_dir, args.force)

    for idt in time_range:
        logger.debug(f'Processing {idt:%Y-%m-%d %H:%M}')
        vrbls_all_exp = []
        # Prepare cubes
        for cubelist in list_of_cubelists:
            if args.track_fname:
                subset['ilon'], subset['ilat'] = (tr_df[tr_df.time == idt]
                                                  [['lon', 'lat']].values[0])
            # vrbls_all_exp.append(prepare_data(cubelist, varnames,
            #                                   idt, level_dict, {}, # subset,
            #                                   ))
            cl = prepare_data(cubelist, varnames,
                              idt, level_dict, {})
            if len(subset) > 0:
                vrbls_all_exp.append(arke.subset.subset_cubelist(cl, subset))
            else:
                vrbls_all_exp.append(cl)

        # logger.debug(vrbls_all_exp)

        fig, axgr, mapkey = prepare_map(vrbls_all_exp, geoax=args.geoax)
        plot_on_same_fig(fig, axgr, mapkey, vrbls_all_exp, geoax=args.geoax)

        imgname = imgname_mask.format(vrbls=lev_str,
                                      track=track_str,
                                      subdir=subdir_out,
                                      time=idt.strftime('%Y%m%d%H%M'),
                                      fmt=fmt)
        fig.savefig(output_dir / imgname, **svfigkw)
        logger.debug('Saved to {}'.format(output_dir / imgname))
        plt.close()
    if args.gif:
        fnames = output_dir / f'*.{fmt}'
        gifname = imgname_mask.format(vrbls=lev_str,
                                      track=track_str,
                                      subdir=subdir_out,
                                      time='',
                                      fmt='gif')
        misc.merge_to_gif(fnames, gifname, delay=25, resize=100,
                          output_dir=output_dir)


if __name__ == '__main__':
    sys.exit(main())
