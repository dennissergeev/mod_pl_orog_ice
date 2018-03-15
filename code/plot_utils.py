# -*- coding: utf-8 -*-
"""
Various plotting functions
"""
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
import string

import arke

import misc_utils as misc


cc = plt.rcParams['axes.prop_cycle']
# Common plotting settings
CBARKW = dict(orientation='vertical')
AXGR_KW = dict(axes_pad=0.2,
               cbar_location='right',
               cbar_mode='single',
               cbar_pad=0.2,
               cbar_size='3%')
mstart_kw = dict(marker='x', ms=8)
mfin_kw = dict(marker='o', mfc='w', ms=8)
# Common cartopy settings
map_kw = dict(transform=ccrs.PlateCarree())
COAST = dict(scale='50m', alpha=0.5,
             edgecolor='#333333', facecolor='#AAAAAA')
clon = 10
clat = 76
extent = [-14, 33, 68, 83]
LCC_KW = dict(clon=clon, clat=clat, coast=COAST, extent=extent,
              ticks=None)



def plotter(fig, axgr, mapkey, vrbls2plot, geoax=False):
    ax_labels = iter(string.ascii_lowercase)
    cax = axgr.cbar_axes[0]
    cbar_ready = False
    for count, (axrow, (fcst_str, times_list)) in enumerate(zip(axgr.axes_row, vrbls2plot.items())):
        for nn, (ax, vrbls) in enumerate(zip(axrow, times_list)):
            ttl = []
            _dt = None
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
                            dx = u.attributes['um_res'].to_flt('km')
                            x = np.arange(u.shape[1]) * dx
                            y = np.arange(u.shape[0]) * dx
                            ax.set_xlabel('Distance along x-axis, km') # , fontsize=10)
                            ax.set_ylabel('Distance along y-axis, km') # , fontsize=10)
                            ax.set_xlim(x.min(), x.max())
                            ax.set_ylim(y.min(), y.max())
                        kw = u.attributes.get('plt_kw', {})
                        if u.attributes['plt_method'] == 'streamplot':
                            ws = (u.data**2 + v.data**2)**0.5
                            kw.update(linewidth=ws/ws.max())
                        if u.attributes['plt_method'] == 'quiver':
                            xydim = len(x.shape)
                            skip = [slice(None, None, 20)]
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
                        ax.set_xlabel('Distance along x-axis, km') # , fontsize=10)
                        ax.set_ylabel('Distance along y-axis, km') # , fontsize=10)
                        ax.set_xlim(x.min(), x.max())
                        ax.set_ylim(y.min(), y.max())
                    kw = icube.attributes.get('plt_kw', {})
                    p = method(x, y, icube.data*scl, **mapkey, **kw)

                    cbar = icube.attributes.get('colorbar')
                    if isinstance(cbar, dict) and not cbar_ready:
                        CBARKW.update(cbar)
                        cb = fig.colorbar(p, cax=cax, **CBARKW)
                        cb.ax.set_title(scl_str, fontsize='xx-large')
                        cb.ax.tick_params(labelsize='x-large')
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
                        fcst_per = icube.coord('forecast_period')[0]
                        fcst_per.convert_units('hours')
            ax_lab = next(ax_labels)
            txt = (f'({ax_lab}) {_dt:%b %d %Y}\n{_dt:%H%M}UTC'
                   f' (+{fcst_per.points[0]:.0f}h)')
            at = AnchoredText(txt, frameon=True, loc=2, prop=dict(size='large')) 
            ax.add_artist(at)


def div_cmap(numcolors=11, name='bwr_div_cmap',
             mincol='blue', midcol='white', maxcol='red',
             under=None, midcol_alpha=0, over=None):
    """
    Create a custom diverging colormap with three colors

    Default is blue to transparent white to red with 11 colors.
    Colors can be specified in any way understandable
    by matplotlib.colors.ColorConverter.to_rgb()
    """
    c_max = mpl.colors.colorConverter.to_rgba(maxcol)
    c_min = mpl.colors.colorConverter.to_rgba(mincol)
    c_mid = mpl.colors.colorConverter.to_rgba(midcol, alpha=midcol_alpha)
    cmap = mpl.colors.LinearSegmentedColormap.from_list(name=name,
                                                        colors=[c_min,
                                                                c_mid, c_max],
                                                        N=numcolors)
    if under is not None:
        cmap.set_under(under)
    if over is not None:
        cmap.set_over(over)
    return cmap