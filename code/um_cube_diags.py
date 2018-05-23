# -*- coding: utf-8 -*-
"""
Common functions to compute diags from the MetUM output
"""
# Standard packages
import datetime
import iris
import numpy as np
# My packages
from arke.coords import nearest_tval
from arke.io import extract_as_single_cube
from arke.grid import unrotate_uv
from arke.numerics import (AtmosFlow,
                           replace_dimcoord,
                           prepare_cube_on_model_levels)
import arke.met_const as mconst
# Local files
from um_cube_params import CUBE_PARAMS
from um_cube_params import wind_aliases_p as wind_aliases


def _get_required_cubes(cubelist, req_vars,
                        date_time=None, cnstr=iris.Constraint()):
    cubes = dict()
    for k, v in req_vars.items():
        cubes[k] = extract_as_single_cube(cubelist, v & cnstr)
        if isinstance(date_time, datetime.datetime):
            t_cnstr = iris.Constraint(time=nearest_tval(cubes[k],
                                                        date_time))
            cubes[k] = cubes[k].extract(t_cnstr)
    return cubes


def compute(cubelist, diag, date_time=None, cnstr=iris.Constraint(),
            unrotatewinds=True, **kw):
    """
    Compute field using the available cubes

    E.g. compute vorticity using wind components
    """
    kwargs = dict(cnstr=cnstr, date_time=date_time)
    if diag in ('vort_z', 'rvstretch', 'rvtilt', 'rvhadv', 'rvvadv',
                'div_h', 'wspd', 'wspd_h', 'egr',
                'stream', 'quiver',
                'frontgen',
                'dfmstretch', 'dfmshear'):
        try:
            cubes = _get_required_cubes(cubelist, wind_aliases, **kwargs)
        except:
            only_u_and_v = {k: v for k, v in wind_aliases.items() if k != 'w'}
            cubes = _get_required_cubes(cubelist, only_u_and_v, **kwargs)
        _u = cubes['u']

        cubes['v'] = cubes['v'].regrid(cubes['u'], iris.analysis.Linear())
        if unrotatewinds:
            cubes['u'], cubes['v'] = unrotate_uv(cubes['u'], cubes['v'])
        if cubes['u'].ndim < 3:
            kw.update(prep_zcoord_kw=None, rm_surf_alt=None,
                      rm_sigma=None, rm_aux_factories=None)
        for k in cubes.keys():
            dx = cubes[k].attributes['um_res'].to_flt()
            kw.update(lonlat2cart_kw=dict(dx=dx))
            cubes[k] = prepare_cube_on_model_levels(cubes[k], **kw)

        if 'w' in cubes:
            cubes['w'].rename('z_wind')
        else:
            cubes['w'] = None

    if diag in ('stream', 'quiver'):
        cube = (cubes['u'], cubes['v'])

    elif diag == 'vort_z':
        # Compute relative vorticity using arke.numerics package
        # AF = AtmosFlow(*[cubes[k] for k in wind_aliases.keys()])
        # cube = cubes['w']
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.rel_vort, _u)

    elif diag == 'rvstretch':
        # Compute relative vorticity budget stretching term
        AF = AtmosFlow(**cubes)
        cube = -1 * replace_dimcoord(AF.rel_vort_stretch, _u)

    elif diag == 'rvtilt':
        # Compute relative vorticity budget tilting term
        AF = AtmosFlow(**cubes)
        cube = -1 * replace_dimcoord(AF.rel_vort_tilt, _u)

    elif diag == 'rvhadv':
        # Compute relative vorticity budget horizontal advection term
        AF = AtmosFlow(**cubes)
        cube = -1 * replace_dimcoord(AF.rel_vort_hadv, _u)

    elif diag == 'rvvadv':
        # Compute relative vorticity budget vertical advection term
        AF = AtmosFlow(**cubes)
        cube = -1 * replace_dimcoord(AF.rel_vort_vadv, _u)

    elif diag == 'div_h':
        # Compute horizontal divergence using arke.numerics package
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.div_h, _u)

    elif diag == 'dfmstretch':
        # Compute stretching deformation using arke.numerics package
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.dfm_stretch, _u)

    elif diag == 'dfmshear':
        # Compute shearing deformation using arke.numerics package
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.dfm_shear, _u)

    elif diag == 'wspd':
        # Compute wind speed using arke.numerics package
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.wspd, _u)

    elif diag == 'wspd_h':
        # Compute horizontal wind speed using arke.numerics package
        cubes.pop('w', None)
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.wspd, _u)

    elif diag == 'relh':
        # Relative humdity
        req_vars = dict(temp='air_temperature', spechum='specific_humidity')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        AF = AtmosFlow(**cubes, cartesian=False)
        cube = AF.relh
        cube.convert_units('%')

    elif diag == 'theta':
        req_vars = dict(temp='air_temperature')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        AF = AtmosFlow(**cubes, cartesian=False)
        cube = AF.theta

    elif diag == 'slp_anom':
        slp_name = CUBE_PARAMS['slp']['_cf_name']
        cube = extract_as_single_cube(cubelist, slp_name & cnstr)
        cube = cube - cube.collapsed('time', iris.analysis.MEAN)
        if isinstance(date_time, datetime.datetime):
            t_cnstr = iris.Constraint(time=nearest_tval(cube,
                                                        date_time))
            cube = cube.extract(t_cnstr)
        cube.rename(CUBE_PARAMS[diag]['_cf_name'])
    elif diag == 'precip':
        rain = extract_as_single_cube(cubelist,
                                      'stratiform_rainfall_rate')
        snow = extract_as_single_cube(cubelist,
                                      'stratiform_snowfall_rate')
        cube = rain + snow
        if isinstance(date_time, datetime.datetime):
            t_cnstr = iris.Constraint(time=nearest_tval(cube,
                                                        date_time))
            cube = cube.extract(t_cnstr)
        cube /= iris.coords.AuxCoord(mconst.rho_l.data,
                                     units=mconst.rho_l.units)
        cube.convert_units('mm h-1')
        cube.rename(CUBE_PARAMS[diag]['_cf_name'])

    elif diag == 'thetae':
        req_vars = dict(temp='air_temperature', spechum='specific_humidity')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        AF = AtmosFlow(**cubes, cartesian=False)
        cube = AF.thetae

    elif diag == 'SSTmT':
        t = extract_as_single_cube(cubelist,
                                   'air_temperature')
        sst = extract_as_single_cube(cubelist,
                                     'surface_temperature')
        sst = sst.regrid(t, iris.analysis.Linear())
        sst.data = np.ma.masked_less(sst.data, 271.35)
        cubes = []
        for cube in (t, sst):
            if isinstance(date_time, datetime.datetime):
                t_cnstr = iris.Constraint(time=nearest_tval(cube,
                                                            date_time))
                cubes.append(cube.extract(t_cnstr))
        cube = cubes[1] - cubes[0]
        cube.rename(CUBE_PARAMS[diag]['_cf_name'])

    elif diag == 'SSTmTheta':
        req_vars = dict(temp='air_temperature')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        AF = AtmosFlow(**cubes, cartesian=False)
        t = AF.theta
        sst = extract_as_single_cube(cubelist,
                                     'surface_temperature')
        sst = sst.regrid(t, iris.analysis.Linear())
        sst.data = np.ma.masked_less(sst.data, 271.35)
        cubes = []
        for cube in (t, sst):
            if isinstance(date_time, datetime.datetime):
                t_cnstr = iris.Constraint(time=nearest_tval(cube,
                                                            date_time))
                cubes.append(cube.extract(t_cnstr))
        cube = cubes[1] - cubes[0]
        cube.rename(CUBE_PARAMS[diag]['_cf_name'])

    elif diag == 'bvfreq':
        req_vars = dict(temp='air_temperature')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        AF = AtmosFlow(**cubes, cartesian=False)
        cube = AF.brunt_vaisala_squared**0.5
        cube.rename('brunt_vaisala_frequency_in_air')
        cube.convert_units('s-1')

    elif diag == 'bvfreqsq':
        req_vars = dict(temp='air_temperature')  # , pres='air_pressure')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        AF = AtmosFlow(**cubes, cartesian=False)
        cube = AF.brunt_vaisala_squared

    elif diag == 'meanbvfreqsq':
        req_vars = dict(temp='air_temperature')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        AF = AtmosFlow(**cubes, cartesian=False)
        cube = AF.brunt_vaisala_squared
        cube = cube.collapsed('pressure', iris.analysis.MEAN)

    elif diag == 'egr':
        req_vars = dict(temp='air_temperature')
        cubes.pop('w', None)
        cubes.update(_get_required_cubes(cubelist, req_vars, **kwargs))
        cubes['temp'] = prepare_cube_on_model_levels(cubes['temp'], **kw)
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.eady_growth_rate, _u)
        cube.convert_units('day-1')

    elif diag == 'gradT':
        req_vars = dict(temp='air_temperature')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        _t = [v for v in cubes.values()][0].copy()
        for k in cubes.keys():
            dx = cubes[k].attributes['um_res'].to_flt()
            if cubes[k].ndim < 3:
                kw.update(prep_zcoord_kw=None, rm_surf_alt=None,
                          rm_sigma=None, rm_aux_factories=None)
            kw.update(lonlat2cart_kw=dict(dx=dx))
            cubes[k] = prepare_cube_on_model_levels(cubes[k], **kw)
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.hgradmag('air_temperature'), _t)
        cube.rename('gradient_magnitude_of_air_temperature')
        cube.convert_units('K m-1')

    elif diag == 'gradThetae':
        req_vars = dict(temp='air_temperature')
        cubes = _get_required_cubes(cubelist, req_vars, **kwargs)
        _t = [v for v in cubes.values()][0].copy()
        for k in cubes.keys():
            dx = cubes[k].attributes['um_res'].to_flt()
            if cubes[k].ndim < 3:
                kw.update(prep_zcoord_kw=None, rm_surf_alt=None,
                          rm_sigma=None, rm_aux_factories=None)
            kw.update(lonlat2cart_kw=dict(dx=dx))
            cubes[k] = prepare_cube_on_model_levels(cubes[k], **kw)
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.hgradmag(alias='thetae'), _t)
        cube.rename('gradient_magnitude_of_equivalent_potential_temperature')
        cube.convert_units('K m-1')

    elif diag == 'frontgen':
        req_vars = dict(temp='air_temperature')
        cubes.pop('w', None)
        cubes.update(_get_required_cubes(cubelist, req_vars, **kwargs))
        cubes['temp'] = prepare_cube_on_model_levels(cubes['temp'], **kw)
        AF = AtmosFlow(**cubes)
        cube = replace_dimcoord(AF.kinematic_frontogenesis, _u)
        cube.convert_units('K (100km)-1 (3h)-1')

    if isinstance(cube, iris.cube.Cube):
        try:
            cube.coord('pressure').convert_units('hPa')
        except iris.exceptions.CoordinateNotFoundError:
            pass

    return cube
