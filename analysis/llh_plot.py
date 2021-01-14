# usual
import os
import sys
import copy
import numpy as np
import scipy.optimize as opt
from scipy.interpolate import interp1d, interp2d, griddata

# plotting
import matplotlib as mpl
mpl.use('agg')
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plotcorner(llh_arr5d, params,
               interp=True, interp_n=25,
               figsize=(8,6),
               cmap=plt.cm.viridis,
               show=True, 
               save=False, outname='test.pdf'):
    '''
    Corner plot creator for 5-dimensional LLH array from
    POCAM PPC simulations for low-energy. 
    '''
    
    ### get parameters
    grid = np.round(params['grid'], 3)

    # data / scan
    scan_dict  = params['scan_dict']
    scan_N     = scan_dict['n']
    scan_Nph   = scan_dict['nph']
    scan_abs   = scan_dict['abs']
    scan_sca   = scan_dict['sca']
    scan_dome  = scan_dict['domeff']
    scan_p0    = scan_dict['p0']
    scan_p1    = scan_dict['p1']
    
    # truth
    truth_dict = params['truth_dict']
    t_N        = truth_dict['n']
    t_Nph      = truth_dict['nph']
    t_abs      = truth_dict['abs']
    t_sca      = truth_dict['sca']
    t_dome     = truth_dict['domeff']
    t_p0       = truth_dict['p0']
    t_p1       = truth_dict['p1']

    # get llh min 5d-index
    min_idx = np.hstack(np.where(np.amin(llh_arr5d) == llh_arr5d))

    # prepare scan dimensions in proper order
    plot_ts  = [t_abs, t_sca, t_dome, t_p0, t_p1]
    plot_xs  = [scan_abs, scan_sca, scan_dome, scan_p0, scan_p1]
    plot_xl  = ['Absorption coefficient', 'Scattering coefficient', 
                'DOM efficiency', 'P0', 'P1']
    plot_lim = [[np.min(i), np.max(i)] for i in plot_xs]
    
    # get input shape from parameters with dimension > 1
    dims  = [i for i, x in enumerate(plot_xs) if len(x) > 1]
    ndim  = len(dims)
    
    # create corner plot
    # dimensionality = ndim, plots = 2*ndim - 1
    # example: 2d --> 1x 2d-grid, 2x 1d-curve = 3 plots
    fig = plt.figure()

    # whether to interpolate the LLH profile with n points
    n = interp_n

    # iterate through dimensions for y-axis
    cntr = 1
    axs  = []
    for i, dim1 in enumerate(dims):
            
        # iterate through dimensions for x-axis
        for j, dim2 in enumerate(dims):
                    
            # get diagonal entries
            if i == j:
                            
                # create subplot
                ax = plt.subplot(ndim, ndim, cntr)
        
                # get scan x interval
                x = plot_xs[dim2] # or dim2 as dim1=dim2 for diagonals
                                  
                # get scan y interval
                # need adjustable 5d-index, use slice method directly
                llh_x = np.zeros((len(x),))
                for xii, xi in enumerate(x):
                    # get part of llh array that is variable given point xi
                    ix  = (grid[:,dim2] == xi) # or dim1, dim1 = dim2 for diagonals
                    zi  = llh_arr5d[ix]
                    # get minimum llh for this point
                    ix2  = np.where(zi == np.amin(zi))
                    llhi = zi[ix2][0]
                    # and put to array
                    llh_x[xii] = llhi
                
                # put to plot array
                y = llh_x
                
                # plot discrete points
                ax.scatter(x, y, color='cornflowerblue', s=8, alpha=1)
                
                # whether to interpolate
                if interp:
                    xn = np.linspace(np.min(x), np.max(x), n)
                    f  = interp1d(x, y, kind='quadratic')
                    x  = xn
                    y  = f(xn)
                    
                # plot
                ax.plot(x, y, color='cornflowerblue', alpha=0.5)
                
                # add truth vertical lines
                ax.axvline(plot_ts[dim1], color='red', ls=':', alpha=0.5) # or dim2, dim1 = dim2 for diagonals
                
                # add llh min vertical lines
                ix = np.where(y == np.amin(y))
                ax.axvline(x[ix][0], color='cyan', ls=':', alpha=0.5) # or dim2, dim1 = dim2 for diagonals            
                
                # add labels
                if i == 0:
                    ax.set_ylabel('LLH')
                if i == len(dims) - 1:
                    ax.set_xlabel(plot_xl[dim1]) # or dim2, dim1 = dim2 for diagonals
                
                # store axis for sharing of columns
                axs.append(ax)
                
                # make aspect fit llh profile colormeshes
                divider = make_axes_locatable(ax)
                cax     = divider.append_axes('right', size='5%', pad=0.05)
                cax.axis('off')
                
            # get off diagonal triangle entries
            if i > j:
                
                # create subplot
                ax = plt.subplot(ndim, ndim, cntr, sharex=axs[j])
                
                # get x and y intervals
                x = plot_xs[dim2]
                y = plot_xs[dim1]
                # get z profile
                # for every pixel (xi, yj) choose remaining parameters so that llh is minimal    
                llh_xy = np.zeros((len(y), len(x)))
                for xii, xi in enumerate(x):
                    for yii, yi in enumerate(y):
                        # get part of llh array that is variable given pixel
                        ix  = (grid[:,dim2] == xi) & (grid[:,dim1] == yi)
                        zi  = llh_arr5d[ix]
                        # get minimum llh for this pixel
                        if len(zi) > 0:
                            ix2  = np.where(zi == np.amin(zi))
                            llhi = zi[ix2][0]
                        else:
                            llhi = -1
                            
                        # put to grid
                        llh_xy[yii, xii] = llhi
                
                # put to plot array
                z = llh_xy
                
                # mask invalid 
                mz = (z >= 0)
                
                # create pcolormesh from xyz data
                dx = np.abs((np.unique(x)[1:] - np.unique(x)[:-1])/2)[0]
                dy = np.abs((np.unique(y)[1:] - np.unique(y)[:-1])/2)[0]
                xl = np.unique(np.append(x - dx, sorted(x + dx)[-1]))
                yl = np.unique(np.append(y - dy, sorted(y + dy)[-1]))
                
                # find range of color mapping
                vmin = np.min(z[mz])
                vmax = np.max(z)
                
                # set invalid values to white
                cmap.set_under('white')
                
                # plot mc truth
                tx = plot_ts[dim2]
                ty = plot_ts[dim1]
                ax.scatter(tx, ty, color='r', s=100, label='MC truth', marker='*', zorder=100)
                
                # interpolate
                if interp:
                    # get rid of invalid pixels
                    zn     = copy.deepcopy(z)
                    zn[zn<0] = np.nan
                    xx, yy = np.meshgrid(x, y)
                    arr    = np.ma.masked_invalid(zn)
                    x1     = xx[~arr.mask]
                    y1     = yy[~arr.mask]
                    narr   = arr[~arr.mask]
                    znew   = griddata((x1, y1), narr.ravel(), (xx, yy), method='cubic')
                    
                    # define new interpolation range
                    xn   = np.linspace(x.min(), x.max(), n)
                    yn   = np.linspace(y.min(), y.max(), n)
                    dxn  = np.abs((np.unique(xn)[1:] - np.unique(xn)[:-1])/2)[0]
                    dyn  = np.abs((np.unique(yn)[1:] - np.unique(yn)[:-1])/2)[0]
                    xxn  = np.unique(np.append(xn - dxn, sorted(xn + dxn)[-1]))
                    yyn  = np.unique(np.append(yn - dyn, sorted(yn + dyn)[-1]))
                    xxx, yyy = np.meshgrid(xxn, yyn)
                    xxxm, yyym = np.meshgrid(xn, yn)
                    
                    # scipy interpolate 2d
                    f    = interp2d(np.unique(x), np.unique(y), znew, kind='cubic')
                    zzz  = f(xn, yn)
                    
                    # find range of color mapping
                    vmin = np.min(zzz[zzz >= 0])
                    vmax = np.max(zzz)
                    
                    #else:
                    # plot interpolated colormesh 
                    im = ax.pcolormesh(xxx, yyy, zzz, vmin=vmin, vmax=vmax, cmap=cmap)

                    # plot llh min
                    idx = np.where(zzz == np.min(zzz[zzz >= 0]))
                    lx  = xxx[idx]
                    ly  = yyy[idx]
                    plt.scatter(lx, ly, color='cyan', s=150, label='LLH min', marker='*')

                    # contours
                    zmin = np.min(zzz[zzz >= 0])    
                    sigmas = np.array([1,3])
                    cs = ax.contour(xxxm, yyym, zzz, levels=zmin + 2.7 * sigmas, 
                                    colors=['white'], linestyles=['--'], label='Confidence level')
                    # contour labels
                    fmt  = {}
                    strs = [r'$%i\sigma$' %sig for sig in sigmas]
                    for l, s in zip(cs.levels, strs):
                        fmt[l] = s
                    ax.clabel(cs, cs.levels, inline=True, fmt=fmt)
                
                else:
                    # plot colormesh
                    im = ax.pcolormesh(xl, yl, z, vmin=vmin, vmax=vmax, cmap=cmap)
                    
                    # plot llh min
                    if len(z[mz]) > 0 and np.max(z[mz]) > 0:
                        index_1, index_2 = np.where(z == np.min(z[mz]))
                        lx = x[int(index_2[0])]
                        ly = y[int(index_1[0])]
                        ax.scatter(lx, ly, color='cyan', s=150, label='LLH min', marker='*')
                    
                # add colorbar
                norm    = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=False)
                divider = make_axes_locatable(ax)
                cax     = divider.append_axes('right', size='5%', pad=0.05)
                cb      = mpl.colorbar.ColorbarBase(cax, format='%i', norm=norm, cmap=cmap)
                cb.set_ticks(np.linspace(vmin, vmax, 5))
                
                # add labels
                # to y-axis on leftmost column
                if j == 0:
                    ax.set_ylabel(plot_xl[dim1])
                # to x-axis on lowest row
                if i == len(dims)-1:
                    ax.set_xlabel(plot_xl[dim2])
                    
            # set xlims
            ax.set_xlim(plot_lim[dim2])
                    
            # advance plot counter
            cntr += 1      

    # add legend
    legend_elements = [Line2D([0], [0], color='cornflowerblue', lw=2, label='1D LLH'),
                    Line2D([0], [0], color='r', lw=2, ls=':', label='1D MC truth'),
                    Line2D([0], [0], color='cyan', lw=2, ls=':', label='1D LLH min'),
                    Patch(facecolor='#31678e', label='2D LLH'),
                    Line2D([0], [0], marker='*', color='r', ls='', label='2D MC truth'),
                    Line2D([0], [0], marker='*', color='cyan', ls='', label='2D LLH min')]
    ax = plt.subplot(ndim, ndim, ndim)
    ax.legend(handles=legend_elements, loc='center', frameon=False)
    plt.axis('off')

    plt.tight_layout()
    
    if save:
        fig.savefig(outname)
        
    if show:
        plt.show()
