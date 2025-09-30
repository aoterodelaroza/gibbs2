import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def plot_phase_diagram(fig,ax,phlist,colorbar=True,steprefine=10):
    """Bulid a colormap plot of the phase diagram constructed with the
    phases given in the Phase list phlist and add it to the provided
    figure and axes.

    By default, the widest temperature and pressure range possible is
    used.

    colorbar = add a vertical colorbar to the plot.

    steprefine = use a step in temperature and pressure that is
    steprefine times smaller than the smallest step in any of the
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

    ## flag the nans
    for G in Gmlist:
        Gmid[np.isnan(G)] = -1

    ## list of colors
    clist = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ## the colorbar goes up to 10 kJ/mol
    norm = mcolors.Normalize(vmin=0, vmax=10, clip=True)

    ## build the contour plot, one per phase
    for i,ph in enumerate(phlist):
        Gtemp = np.ones_like(Gmid)
        Gtemp[Gmid != i] = 0

        # contour lines
        levels = np.array([0.5])
        ax.contour(pm,Tm,Gtemp,levels,colors='k',linestyles='solid')

        # colormap
        Gtemp = Gmmin2 - Gmmin
        Gtemp[Gmid != i] = np.nan
        cmap = mcolors.LinearSegmentedColormap.from_list('custom',[(0,'#ffffff'),(1,clist[i % len(clist)])],N=256)
        cnt = ax.contourf(pm, Tm, Gtemp, levels=256, cmap=cmap, norm=norm)

    ## colorbar
    if colorbar:
        cmap = mcolors.LinearSegmentedColormap.from_list('custom',[(0,'#ffffff'),(1,'#000000')],N=256)
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),orientation='vertical',ax=ax,pad=0.01)

    ## hatching for nans
    Gtemp = np.ones_like(Gmid)
    Gtemp[Gmid != -1] = 0
    ax.contourf(pm, Tm, Gtemp, hatches=[None,'X'], levels=[-0.5,0.5,1.5], colors='none', cmap=None)

    ## axes
    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("Temperature (K)")

    return fig, ax

def plot_enthalpy_pressure(fig,ax,phlist,steprefine=10):
    """Bulid an enthalpy-pressure plot constructed with the phases
    given in the StaticPhase list phlist.

    By default, the maximum pressure range possible is used.

    steprefine = use a step in pressure that is steprefine times lower
    than the smallest step in any of the phases.
    """

    ## determine plot limits and step
    ## pick as reference the phase with the largest pressure range
    idref = -1
    prange = 0
    pstep = np.inf
    for i,ph in enumerate(phlist):
        pstep = min(pstep,ph.pstep)
        dp = ph.pmax - ph.pmin
        if dp > prange:
            idref = i
            prange = dp
            pmin = ph.pmin
            pmax = ph.pmax

    ## build the grid
    p = np.arange(pmin,pmax,pstep / steprefine)

    ## plot the enthalpies
    for i,ph in enumerate(phlist):
        y = ph.H(p)
        mask = ~np.isnan(y)
        x = p[mask]
        y = y[mask] - phlist[idref].H(x)
        ax.plot(x,y,'-',label=ph._name)

    ## axes
    ax.set_xlim(pmin,pmax)
    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("Enthalpy (kJ/mol)")
    ax.grid()
    ax.legend()

    return fig, ax

def plot_barh_stablephase(fig,ax,phlist,y,thinlines=True,steprefine=10):
    """Bulid a horizontal bar (barh) plot indicating the most table
    phase as a function of pressure in the StaticPhase list phlist.
    Place the plot at height y.

    By default, the maximum pressure range possible is used.

    thinlines = indicate the domain of each phase with thin lines
    under the bar plot.

    steprefine = use a step in pressure that is steprefine times lower
    than the smallest step in any of the phases.
    """

    ## determine plot limits and step
    prange = 0
    pmin = np.inf
    pmax = -np.inf
    pstep = np.inf
    for i,ph in enumerate(phlist):
        pmin = min(ph.pmin,pmin)
        pmax = max(ph.pmax,pmax)
        pstep = min(pstep,ph.pstep)

    ## list of colors
    clist = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ## build the grid
    plist = np.arange(pmin,pmax,pstep / steprefine)
    Hlist = np.zeros((len(plist),len(phlist)))
    for i,ph in enumerate(phlist):
        Hlist[:,i] = np.array(ph.H(plist))

    ## height of the boxes
    ylo, yhi = ax.get_ylim()
    yh = (yhi - ylo) / 30

    width = []
    left = []
    color = []
    iminlast = np.nan
    for i,p in enumerate(plist):
        imin = np.nanargmin(Hlist[i,:])
        if np.isnan(iminlast):
            iminlast = imin
            plast = p
        elif imin != iminlast:
            left.append(plast)
            width.append(p - plast)
            color.append(clist[iminlast % len(clist)])
            iminlast = imin
            plast = p

    ## the actual bar
    ax.set_axisbelow(True)
    ax.barh(y, width=width, height=yh, left=left, color=color, edgecolor='k')

    ## the thin lines under the bar
    if thinlines:
        ypos = y - 0.5 * yh
        for i,ph in enumerate(phlist):
            ypos = ypos - 0.25 * yh
            ax.hlines(ypos,ph.pmin,ph.pmax,colors=clist[i % len(clist)])

    return fig, ax
