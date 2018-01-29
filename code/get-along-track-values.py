# -*- coding: utf-8 -*-
"""
Load polar low positions and times and calculate fields within certain radius
"""
from datetime import datetime
import iris
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
#
from arke import coords, grid, numerics, units
#
import mypaths


iris.FUTURE.netcdf_promote = True
iris.FUTURE.cell_datetime_objects = True

um_res = units.GrdStep('km2p2')

FNAME_MASK = 'umnsa_*'
TOPDIR = mypaths.um_dir_ext
PATH_MASK = TOPDIR / '{fcst_str}' / '{um_res}' / '{run}' / '{subdir}'
PATH_MASK /= FNAME_MASK
TRACKDIR = mypaths.trackdir
OUTDIR = mypaths.trackdir.parent

RAD = 150
plevels = [850, 875, 900, 925, 950]

cf_names = dict(
    plev=['x_wind', 'y_wind', 'air_temperature'],
    sfc=['surface_upward_sensible_heat_flux']
)

varnames = dict(
    plev=['rel_vort'],
    sfc=['shf']
)  # , 'theta', 'thetae', 'brunt_vaisala_squared', 'relh', 'spechum', 'mixr']

unit_dict = dict(rel_vort='s-1', shf='W m-2')

pl_cases = dict(
    STARS72='20070404T1200Z',
    STARS77='20080129T1200Z'
               )
fcst_strs = [*pl_cases.values()]
runs = ['ctrl', 'nosva', 'sva200', 'ice76n', 'ice82n']


tracks = {}
for fcst_str in fcst_strs:
    tracks[fcst_str] = {}
    for run in runs:
        input_dir = TRACKDIR / fcst_str / run
        input_file = f'pl_loc.{fcst_str}.{um_res}.{run}.vort_mslp.txt'
        tracks[fcst_str][run] = pd.read_csv(input_dir / input_file,
                                            parse_dates=[0])


data3d = dict()
for fcst_str in tqdm(fcst_strs):
    data3d[fcst_str] = dict()
    for run in runs:
        data3d[fcst_str][run] = iris.cube.CubeList()
        for subdir, names in cf_names.items():
            input_files = PATH_MASK.format(fcst_str=fcst_str, um_res=um_res,
                                           run=run, subdir=subdir)
            data3d[fcst_str][run] += iris.load(input_files, constraints=names)


plevels = np.array(plevels, dtype=int)
plevel_str = f'{plevels.max()}-{plevels.min()}'
pconstr = iris.Constraint(pressure=lambda x: int(x.point) in plevels)
prep_kw = dict(lonlat2cart_kw=dict(dx=um_res.to_flt()))
prep_kw.update(prep_zcoord_kw=None, rm_surf_alt=None,
               rm_sigma=None, rm_aux_factories=None)


vrbls_along_track = dict()
for fcst_str in tqdm(fcst_strs):
    vrbls_along_track[fcst_str] = dict()
    for run in tqdm(runs, leave=False):
        tr_df = tracks[fcst_str][run]
        tr_datetimes = tr_df.time.values.astype('<M8[ms]').astype(datetime)
        npoints = len(tr_datetimes)

        along_track = {v: np.empty((npoints, ), np.float64)
                       for k in varnames for v in varnames[k]}
        for i, (idt, ilon, ilat) in tqdm(enumerate(zip(tr_datetimes,
                                                       tr_df.lon.values,
                                                       tr_df.lat.values)),
                                         leave=False):
            subdir = 'plev'
            sublist = (data3d[fcst_str][run]
                       .extract(cf_names[subdir]).extract(pconstr))
            cl = iris.cube.CubeList()
            for cube in sublist:
                _tc = iris.Constraint(time=coords.nearest_tval(cube, idt))
                cl.append(cube.extract(_tc))

            temp = cl.extract('air_temperature', strict=True)
            temp.data = np.ma.masked_outside(temp.data, 100, 350)
            u = cl.extract('x_wind', strict=True)
            v = cl.extract('y_wind', strict=True)
            u, v = grid.unrotate_uv(u, v)
            u.data = np.ma.masked_where(temp.data.mask, u.data)
            v.data = np.ma.masked_where(temp.data.mask, v.data)

            cubes = dict(u=u, v=v)
            for key in cubes.keys():
                cubes[key] = numerics.prepare_cube_on_model_levels(cubes[key],
                                                                   **prep_kw)

            AF = numerics.AtmosFlow(cartesian=True, **cubes)
            for var in varnames[subdir]:
                cube = numerics.replace_dimcoord(getattr(AF, var), temp)
                cube = grid.mask_cube_outside_circle_xy(cube, RAD, ilon, ilat,
                                                        dx=um_res.to_flt('km'),
                                                        return_mask=False)
                if cube.data.mask.all():
                    raise ValueError('all masked!!!')
                along_track[var][i] = np.nanmean(cube.data)

            subdir = 'sfc'
            var = 'shf'
            sublist = data3d[fcst_str][run].extract(cf_names[subdir])
            cube = sublist.extract_strict('surface_upward_sensible_heat_flux')
            tc = iris.Constraint(time=coords.nearest_tval(cube, idt))
            cube = cube.extract(tc)
            cube = grid.mask_cube_outside_circle_xy(cube, RAD, ilon, ilat,
                                                    dx=um_res.to_flt('km'),
                                                    return_mask=False)
            if cube.data.mask.all():
                raise ValueError('all masked!!!')
            along_track[var][i] = np.nanmean(cube.data)
        vrbls_along_track[fcst_str][run] = along_track


#
# Write results to HDF5 file
#
global_dict = dict(radius=RAD, pressure_levels=plevels, method='mean')
fname = (f'along_track.mean.{"_".join(fcst_strs)}.{"_".join(runs)}'
         '.{"_".join([v for k in varnames for v in varnames[k]])}'
         '.r{RAD}.p{plevel_str}.h5')
with h5py.File(OUTDIR/fname, 'w') as hf:
    for fcst_str in pl_cases.values():
        g = hf.create_group(fcst_str)
        for run in runs:
            g2 = g.create_group(run)
            for key, val in vrbls_along_track[fcst_str][run].items():
                da = g2.create_dataset(key, data=val)
                da.attrs.update(units=unit_dict[key])
    hf.attrs.update(global_dict)
