# -*- coding: utf-8 -*-
"""
Common settings for UM output plotting
"""
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from plot_utils import div_cmap


vert_coord_params = [dict(name='pressure', units='hPa', subdir='plev'),
                     dict(name='level_height', units='m', subdir='mlev'),
                     dict(name='', units='', subdir='sfc')]


def vert_coord_by_units(units):
    return [{k: v for k, v in d.items() if k != 'units'}
            for d in vert_coord_params if d['units'] == units][0]


wind_aliases = OrderedDict((('u', 'x_wind'),
                            ('v', 'y_wind'),
                            ('w', 'upward_air_velocity')))
wind_aliases_p = OrderedDict((('u', 'x_wind'),
                              ('v', 'y_wind'),
                              ('w', 'lagrangian_tendency_of_air_pressure')))

clev101 = list(np.linspace(-1, 1, 9))
clev101.remove(0)
clev101 = np.array(clev101)

w_cmap = div_cmap(mincol='chartreuse', maxcol='fuchsia',
                  midcol=(1, 1, 1, 0.75),
                  under='darkgreen', over='deeppink')

omega_cmap = w_cmap

_cmap = plt.cm.PuOr_r(np.linspace(0.5, 1, 128))
vort_pos_cmap = mcolors.LinearSegmentedColormap.from_list('Or', _cmap)
vort_pos_cmap.set_over('#330000')
vort_pos_cmap.set_under('w', alpha=0)

vort_cmap = plt.cm.PuOr_r
vort_cmap.set_over('#330000')
vort_cmap.set_under('#000033')

d_cmap = plt.cm.PuOr
d_cmap.set_over('#000033')
d_cmap.set_under('#330000')

pv_cmap = plt.cm.BrBG_r
pv_cmap.set_over('#330000')
pv_cmap.set_under('#000033')

precip_cmap = plt.cm.cool
precip_cmap.set_over('#330055')

temp_cmap = plt.cm.YlGnBu_r
temp_cmap.set_over('#FFFFEE')
temp_cmap.set_under('#222233')

cmap_max = plt.cm.plasma_r
cmap_max.set_under('w', alpha=0)

cmap_max2 = plt.cm.viridis
cmap_max2.set_under('#000000')

cmap_front = plt.cm.plasma_r
cmap_front.set_under('#DDDDDD')
cmap_front.set_over('#FFFFFF')

CUBE_PARAMS = {
'stream': {
         '_cf_name': 'streamline_plot',
         'diag_method': 'stream',
         'extract_level_before_diag': True,
         'plt_method': 'streamplot',
         'plt_kw': dict(density=1, color='#222222'),
         },
'quiver': {
         '_cf_name': 'quiver_plot',
         'diag_method': 'quiver',
         'extract_level_before_diag': True,
         'plt_method': 'quiver',
         'plt_kw': dict(),
         },
'vort': {
         '_cf_name': 'atmospheric_relative_vorticity',
         'diag_method': 'vort_z',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=[1, 5, 10, 20], cmap=vort_pos_cmap, extend='both'),
         'exp_diff_kw': dict(levels=[10, 25, 50, 100]),
         'colorbar': {},
         'scl': 1e4,
         'mark': 'max',
        },
# 'vort': {
#          '_cf_name': 'atmospheric_relative_vorticity',
#          'diag_method': 'vort_z',
#          'diag_kw': {},
#          'extract_level_before_diag': True,
#          'plt_method': 'contourf',
#          'plt_kw': dict(levels=clev101*20, cmap=vort_cmap, extend='both'),
#          'exp_diff_kw': dict(levels=[10, 25, 50, 100]),
#          'colorbar': {},
#          'scl': 1e4,
#          'mark': 'max',
#         },
'rv': {
         '_cf_name': 'atmospheric_relative_vorticity',
         'diag_method': 'vort_z',
         'extract_level_before_diag': True,
         'plt_method': 'contour',
         'plt_kw': dict(levels=[5], colors='C1', linewidths=1.5),
         'scl': 1e4,
        },
'div': {
         '_cf_name': 'divergence_of_wind',
         'diag_method': 'div_h',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=clev101*5, cmap=d_cmap, extend='both'),
         'exp_diff_kw': dict(levels=[-100, -50, -25, -10, -5], linestyles='-'),
         'colorbar': {},
         'scl': 1e4
        },
'rvstretch': {
         '_cf_name': 'stretching_term_of_atmosphere_relative_vorticity_budget',
         'diag_method': 'rvstretch',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=clev101*10, cmap=vort_cmap, extend='both'),
         'colorbar': {},
         'scl': 1e7
        },
'rvtilt': {
         '_cf_name': 'tilting_term_of_atmosphere_relative_vorticity_budget',
         'diag_method': 'rvtilt',
         'diag_kw': {},
         'extract_level_before_diag': False,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=clev101*10, cmap=vort_cmap, extend='both'),
         'colorbar': {},
         'scl': 1e7
        },
'rvhadv': {
         '_cf_name': ('horizontal_advection_term'
                      '_of_atmosphere_relative_vorticity_budget'),
         'diag_method': 'rvhadv',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=clev101*10, cmap=vort_cmap, extend='both'),
         'colorbar': {},
         'scl': 1e7
        },
'rvvadv': {
         '_cf_name': ('vertical_advection_term'
                      '_of_atmosphere_relative_vorticity_budget'),
         'diag_method': 'rvvadv',
         'diag_kw': {},
         'extract_level_before_diag': False,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=clev101*10, cmap=vort_cmap, extend='both'),
         'colorbar': {},
         'scl': 1e7
        },
'dfmstretch': {
         'diag_method': 'dfmstretch',
         '_cf_name': 'stretching_deformation_2d',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=clev101*20, cmap=vort_cmap, extend='both'),
         'colorbar': {},
         'scl': 1e4
        },
'dfmshear': {
         '_cf_name': 'shearing_deformation_2d',
         'diag_method': 'dfmshear',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=clev101*20, cmap=vort_cmap, extend='both'),
         'colorbar': {},
         'scl': 1e4
        },
'wspd': {
         '_cf_name': 'wind_speed',
         'diag_method': 'wspd',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=np.arange(3, 43, 3), cmap='magma_r', extend='max'),
         'exp_diff_kw': dict(levels=np.arange(18, 33, 3)),
         'colorbar': {},
        },
'wspd_h': {
         '_cf_name': 'horizontal_wind_speed',
         'diag_method': 'wspd_h',
         'diag_kw': {},
         'extract_level_before_diag': True,
         'plt_method': 'contourf',
         'plt_kw': dict(levels=np.arange(3, 43, 3), cmap='magma_r', extend='max'),
         'exp_diff_kw': dict(levels=np.arange(18, 33, 3)),
         'colorbar': {},
        },
 'pv': {
          '_cf_name': 'potential_vorticity_of_atmosphere_layer',
          '_cf_units': 'm2 s-1 K kg-1',
          'stash': 'm01s15i229',
          'extract_level_before_diag': True,
          'plt_method': 'contour',
          'plt_kw': dict(levels=np.arange(2, 10, 2), colors='#000FFF', linewidths=1.5),
          'clabels': dict(inline=1, fmt='%1.0f', fontsize=12),
          'scl': 1e6,
        },
# 'pv': {
#          '_cf_name': 'potential_vorticity_of_atmosphere_layer',
#          '_cf_units': 'm2 s-1 K kg-1',
#          'stash': 'm01s15i229',
#          'extract_level_before_diag': True,
#          'plt_method': 'contourf',
#          'plt_kw': dict(levels=np.arange(-5, 6, 1), cmap=pv_cmap, extend='both'),
#          'exp_diff_kw': dict(levels=[1, 5, 10]),
#          'colorbar': {},
#          'scl': 1e6,
#          'mark': 'max',
#        },
'w': {
      '_cf_name': 'upward_air_velocity',
      'plt_method': 'contourf',
      'extract_level_before_diag': True,
      'plt_kw': dict(cmap=w_cmap, levels=clev101, extend='both'),
      'colorbar': {},
     },
'omega': {
      '_cf_name': 'lagrangian_tendency_of_air_pressure',
      'plt_method': 'contourf',
      'extract_level_before_diag': True,
      'plt_kw': dict(cmap=omega_cmap, levels=clev101*2, extend='both'),
      'colorbar': {},
     },
'slp': {
        '_cf_name': 'air_pressure_at_sea_level',
        'plt_method': 'contour',
        'plt_kw': dict(levels=np.arange(952, 1053, 2.), colors='#FF0000', linewidths=0.5),
        'clabels': dict(inline=1, fmt='%1.0f', fontsize=10, colors='#FF0000'),
        'exp_diff_kw': dict(levels=[980, 990, 1000, 1010, 1020]),
        'scl': 1e-2,
        'mark': 'min',
        },
'mslp': {
        '_cf_name': 'air_pressure_at_sea_level',
        'plt_method': 'contour',
        'plt_kw': dict(levels=np.arange(952, 1053, 1.), colors='#000FFF', linewidths=0.5),
        'clabels': dict(inline=1, fmt='%1.0f', fontsize=10, colors='#000FFF'),
        'exp_diff_kw': dict(levels=[980, 990, 1000, 1010, 1020]),
        'scl': 1e-2,
        'mark': 'min',
        'filter': False,
       },
'shf': {
        '_cf_name': 'surface_upward_sensible_heat_flux',
        'plt_method': 'contourf',
        'plt_kw': dict(levels=np.arange(50, 400, 50), cmap='Reds', extend='max'),
        'exp_diff_kw': dict(levels=[200, 300, 400]),
        'colorbar': {}
       },
'lhf': {
        '_cf_name': 'surface_upward_latent_heat_flux',
        'plt_method': 'contourf',
        'plt_kw': dict(levels=np.arange(50, 400, 50), cmap='Reds', extend='max'),
        # 'plt_kw': dict(levels=np.arange(50, 1000, 50), colors='C0', linewidths=0.5), #, cmap='Reds'),
        'exp_diff_kw': dict(levels=[200, 300, 400]),
        'colorbar': {},
        # 'clabels': dict(inline=1, fmt='%1.0f', fontsize=10, colors='C0'),
       },
'lwtoa': {
          '_cf_name': 'toa_outgoing_longwave_flux',
          'plt_method': 'pcolormesh',
          'plt_kw': dict(cmap=plt.cm.gray_r, vmin=140, vmax=230),
          # 'colorbar': {extend='neither', ticks=[140, 170, 200, 230]},
         },
'rainrate': {
           '_cf_name': 'stratiform_rainfall_rate',
           'plt_method': 'contourf',
           'plt_kw': dict(cmap=precip_cmap, extend='max', levels=[0.1, 1.0, 1.5, 1.75, 2.0]),
           'exp_diff_kw': dict(levels=[1.0, 2.0], extend='max'),
           'colorbar': {},
           'scl': 1e3,
          },
'snowrate': {
           '_cf_name': 'stratiform_snowfall_rate',
           'plt_method': 'contourf',
           'plt_kw': dict(cmap=precip_cmap, extend='max', levels=[0.1, 1.0, 1.5, 1.75, 2.0]),
           'exp_diff_kw': dict(levels=[1.0, 2.0], extend='max'),
           'colorbar': {},
           'scl': 1e3,
          },
'precip': {
           '_cf_name': 'stratiform_precipitation_rate',
           'diag_method': 'precip',
           'plt_method': 'contourf',
           'plt_kw': dict(cmap=precip_cmap, extend='max', levels=[1, 2, 5, 10]),
           'exp_diff_kw': dict(levels=[1.0, 2.0], extend='max'),
           'colorbar': {},
          },
'theta': {
          '_cf_name': 'air_potential_temperature',
          'diag_method': 'theta',
          'extract_level_before_diag': True,
          'plt_method': 'contour',
          'plt_kw': dict(levels=np.arange(200, 400, 5), colors='C1'), #cmap=temp_cmap, extend='both'),
          'clabels': dict(inline=1, fmt='%1.0f', fontsize=12, colors='C1'),
          # 'plt_method': 'contourf',
          # 'plt_kw': dict(levels=np.arange(200, 400, 2), colors='#222222'), #cmap=temp_cmap, extend='both'),
          # 'clabels': dict(inline=1, fmt='%1.0f', fontsize=12, colors='#222222'),
          #'plt_kw': dict(levels=np.arange(260, 275, 2), cmap=temp_cmap, extend='both'),
          #'colorbar': {},
          'exp_diff_kw': dict(levels=[264, 267]),
         },
'thetae': {
          '_cf_name': 'equivalent_potential_temperature',
          'diag_method': 'thetae',
          'extract_level_before_diag': True,
          'plt_method': 'contourf',
          # 'plt_kw': dict(levels=np.arange(200, 400, 5), colors='C9'), #cmap=temp_cmap, extend='both'),
          # 'clabels': dict(inline=1, fmt='%1.0f', fontsize=12, colors='C9'),
          'plt_kw': dict(levels=np.arange(250, 281, 5), cmap=temp_cmap, extend='both'),
          'colorbar': {},
          'exp_diff_kw': dict(levels=[264, 267]),
         },
'temp': {
          '_cf_name': 'air_temperature',
          'extract_level_before_diag': True,
          'plt_method': 'contour',
          'plt_kw': dict(levels=np.arange(200, 400, 5), colors='C1'),
          'clabels': dict(inline=1, fmt='%1.0f', fontsize=12, colors='C1'),
          # 'plt_method': 'contourf',
          # 'plt_kw': dict(levels=np.arange(225, 246, 3), cmap=temp_cmap, extend='both'),
          # 'colorbar': {},
          'exp_diff_kw': dict(levels=[250, 255, 260, 265]),
          'ma_out': (200, 350),
         },
'gradT': {
          '_cf_name': 'gradient_magnitude_of_air_temperature',
          'diag_method': 'gradT',
          'scl': 1e3,
          'plt_method': 'contourf',
          'extract_level_before_diag': True,
          'plt_kw': dict(levels=np.arange(0.1, 0.51, 0.1), cmap=cmap_max, extend='both'),
          'colorbar': {},
         },
'gradThetae': {
          '_cf_name': 'gradient_magnitude_of_equivalent_potential_temperature',
          'diag_method': 'gradThetae',
          'scl': 1e3,
          'plt_method': 'contourf',
          'extract_level_before_diag': True,
          'plt_kw': dict(levels=np.arange(0.1, 0.51, 0.1), cmap=cmap_max, extend='both'),
          'colorbar': {},
         },
'frontgen': {
          '_cf_name': 'kinematic_frontogenesis',
          'diag_method': 'frontgen',
          'scl': 1,
          'plt_method': 'contourf',
          'extract_level_before_diag': True,
          'plt_kw': dict(levels=[10, 100, 200, 500, 1000], cmap=cmap_front, extend='both'),
          'colorbar': {},
         },
'egr': {
          '_cf_name': 'eady_growth_rate',
          'diag_method': 'egr',
          'plt_method': 'contourf',
          'plt_kw': dict(levels=np.arange(-5, 6, 1), cmap=w_cmap, extend='both'),
          'colorbar': {},
         },
'bvfreq': {
          '_cf_name': 'brunt_vaisala_frequency_in_air',
          'diag_method': 'bvfreq',
          'plt_method': 'contourf',
          'plt_kw': dict(levels=np.arange(0, 0.0251, 0.0025), cmap=cmap_max2, extend='both'),
          'colorbar': {},
         },
'bvfreqsq': {
          '_cf_name': 'square_of_brunt_vaisala_frequency_in_air',
          'diag_method': 'bvfreqsq',
          'plt_method': 'contourf',
          'plt_kw': dict(levels=clev101, cmap=w_cmap, extend='both'),
          'colorbar': {},
          'scl': 1e4,
         },
'meanbvfreqsq': {
          '_cf_name': 'square_of_brunt_vaisala_frequency_in_air',
          'diag_method': 'meanbvfreqsq',
          'plt_method': 'contourf',
          'plt_kw': dict(levels=clev101*5, cmap=w_cmap, extend='both'),
          'colorbar': {},
          'scl': 1e4,
         },
'dtempdp': {
          '_cf_name': 'derivative_of_air_temperature_wrt_pressure',
          'plt_method': 'contourf',
          'plt_kw': dict(cmap='RdBu_r'),
          'colorbar': {},
          'exp_diff_kw': dict(levels=np.arange(-0.2, 0.21, 0.5)),
         },
'spechum': {
          '_cf_name': 'specific_humidity',
          'plt_method': 'contour',
          'extract_level_before_diag': True,
          'plt_kw': dict(levels=np.arange(1, 10, 1), colors='#000FFF', linewidths=0.5), #, cmap='Reds'),
          'scl': 1e3,
          'clabels': dict(inline=1, fmt='%1.0f', fontsize=12, colors='#000FFF'),
         },
'relh': {
          '_cf_name': 'relative_humidity',
          'diag_method': 'relh',
          'plt_method': 'contourf',
          'extract_level_before_diag': True,
          'plt_kw': dict(levels=np.linspace(0, 100, 5), cmap='viridis_r'),
          'colorbar': {},
         },
'gh': {
       '_cf_name': 'geopotential_height',
       'extract_level_before_diag': True,
       'plt_method': 'contour',
       'plt_kw': dict(colors='magenta', levels=np.arange(0, 10000, 20)),
       'clabels': dict(inline=1, fmt='%1.0f', fontsize=12, colors='magenta'),
#       'scl': 1e-1,
       'exp_diff_kw': dict(levels=np.arange(0, 10000, 100)),
         },
'seaice': {
           '_cf_name': 'sea_ice_area_fraction',
           'plt_method': 'contour',
           'plt_kw': dict(levels=[0.15], colors='#333333', linewidths=3, alpha=0.75),
           'exp_diff_kw': dict(levels=[0.15], linestyles='--'),
          },
'orog': {
           '_cf_name': 'surface_altitude',
           'plt_method': 'contourf',
           'plt_kw': dict(cmap='cubehelix', alpha=0.5, levels=np.arange(0, 1500, 200)),
          },
'blh': {
        '_cf_name': 'atmosphere_boundary_layer_thickness',
        'plt_method': 'contourf',
        'plt_kw': dict(levels=np.arange(250, 3100, 250), cmap='viridis', extend='both'),
        'exp_diff_kw': dict(levels=[2000, 3000]),
        'colorbar': {},
       },
'SSTmT': {
           '_cf_name': 'sst_minus_t',
           'extract_level_before_diag': True,
           'diag_method': 'SSTmT',
           'plt_method': 'contour',
           'plt_kw': dict(levels=np.arange(40, 51, 1), cmap='plasma', linestyles=['--'], linewidths=2),
           'clabels': dict(inline=1, fmt='%1.0f', fontsize=12),
           'exp_diff_kw': dict(levels=[43], linestyles='--'),
          },
'SSTmTheta': {
           '_cf_name': 'sst_minus_theta',
           'extract_level_before_diag': True,
           'diag_method': 'SSTmTheta',
           'plt_method': 'contourf',
           'plt_kw': dict(levels=clev101*10, cmap='RdBu_r', extend='both'),
          },
'slp_anom' : {
              '_cf_name': 'air_pressure_anomaly',
              'diag_method': 'slp_anom',
              'plt.method': 'contourf',
              'plt_kw': dict(cmap='RdBu_r', levels=np.arange(-20, 22, 2)),
              'exp_diff_kw': dict(linestyles='-', levels=np.arange(-20, 0, 2)),
              'colorbar': {},
              'scl': 1e2,
             },
}
