import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def plot_phase_diagram(fig,ax,phlist,steprefine=10):
    """Bulid a colormap plot of the phase diagram constructed with the
    phases given in the list phlist.

    By default, the maximum possible temperature and pressure range
    possible is used.

    steprefine = use a step in temperature and pressure that is
    steprefine times lower than the smallest step in any of the
    phases. A higher steprefine means a slower plot but smoother
    contours.
    """

    ## determine plot limits and step
    Tmin = pmin = np.inf
    Tmax = pmax = -np.inf
    Tstep = pstep = np.inf
    for ph in phlist:
        Tmin = min(Tmin,ph.Tmin)
        pmin = min(pmin,ph.pmin)
        Tmax = max(Tmax,ph.Tmax)
        pmax = max(pmax,ph.pmax)
        Tstep = min(Tstep,ph.Tstep)
        pstep = min(pstep,ph.pstep)

    ## build the grid
    p = np.arange(pmin,pmax,pstep / steprefine)
    T = np.arange(Tmin,Tmax,Tstep / steprefine)
    pm,Tm = np.meshgrid(p,T)

    ## build the minimum G map and the minimum id phase map
    Gmlist = []
    Gmmin = np.full_like(pm,np.inf)
    Gmmin2 = np.full_like(pm,np.inf)
    Gmid = -np.ones_like(pm)
    for i,ph in enumerate(phlist):
        Gmlist.append(ph.G(pm,Tm))
        mask = (Gmlist[-1] < Gmmin) & (~np.isnan(Gmlist[-1]))
        Gmmin2[mask] = Gmmin[mask]

        Gmmin[mask] = Gmlist[-1][mask]
        Gmid[mask] = i

        mask = (~mask) & (Gmlist[-1] < Gmmin2) & (~np.isnan(Gmlist[-1]))
        Gmmin2[mask] = Gmlist[-1][mask]

    ## put back the nans
    for G in Gmlist:
        Gmid[np.isnan(G)] = -1

    ## list of colors
    clist = plt.rcParams['axes.prop_cycle'].by_key()['color']

    norm = mcolors.Normalize(vmin=0, vmax=10, clip=True)
    ## build the contour plot, one per phase
    for i,ph in enumerate(phlist):
        Gtemp = np.ones_like(Gmid)
        Gtemp[Gmid != i] = 0

        # contour line
        levels = np.array([0.5])
        ax.contour(pm,Tm,Gtemp,levels,colors='k',linestyles='solid')

        # colormap
        Gtemp = Gmmin2 - Gmmin
        Gtemp[Gmid != i] = np.nan
        cmap = mcolors.LinearSegmentedColormap.from_list('custom',[(0,'#ffffff'),(1,clist[i])],N=256)
        cnt = ax.contourf(pm, Tm, Gtemp, levels=256, cmap=cmap, norm=norm)

    ## colorbar
    cmap = mcolors.LinearSegmentedColormap.from_list('custom',[(0,'#ffffff'),(1,'#000000')],N=256)
    cnt = cm.ScalarMappable(norm=norm,cmap=cmap)
    fig.colorbar(cnt,orientation='vertical',ax=ax,pad=0.01)

    ## hatching for nans
    Gtemp = np.ones_like(Gmid)
    Gtemp[Gmid != -1] = 0
    ax.contourf(pm, Tm, Gtemp, hatches=[None,'X'], levels=[-0.5,0.5,1.5], colors='none', cmap=None)

    ## axes
    ax.set_xlabel(r"Pressure (GPa)")
    ax.set_ylabel("Temperature (K)")

    return fig, ax
