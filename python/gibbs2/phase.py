import numpy as np
from scipy.interpolate import CloughTocher2DInterpolator
import copy

class Phase:
    """Thermodynamic properties of a phase calculated by gibbs2. This class
    is derived from the information written to a .eos file and is used
    to create plots and performing tasks not done in gibbs2 by default."""

    def __init__(self,name,p,T,V,G):
        """Initialize phase given name and 2D grid of thermodynamic
        properties."""

        ## assign fields
        self._name = name
        self._plist = np.copy(p)
        self._Tlist = np.copy(T)
        self._Vlist = np.copy(V)
        self._Glist = np.copy(G)

        ## set up the interpolants
        self._Ginterp = CloughTocher2DInterpolator(list(zip(self._plist, self._Tlist)), self._Glist)
        self._Vinterp = CloughTocher2DInterpolator(list(zip(self._plist, self._Tlist)), self._Vlist)

    @classmethod
    def fromlines(cls,name,strl):
        """Initialize phase from a list of strings with the
        thermodynamic properties obtained from reading the eos file
        (see examples)."""

        x = np.genfromtxt(strl)
        p = x[:,0]
        T = x[:,1]
        V = x[:,2]
        G = x[:,4]
        return cls(name,p,T,V,G)

    @classmethod
    def fromfile(cls,filename):
        """Initialize phase from an eos file. If several phases are
        present, use the first one."""

        return read_phases_from_eos_file(filename)[0]

    ## interpolated thermodynamic properties
    def G(self,p,T):
        """Returns the pT-interpolated Gibbs energy."""
        return self._Ginterp((p,T))
    def V(self,p,T):
        """Returns the pT-interpolated volume."""
        return self._Vinterp((p,T))

    ## some properties: min, max, step temperature and pressure, name
    @property
    def Tmin(self):
        return np.min(self._Tlist)
    @property
    def Tmax(self):
        return np.max(self._Tlist)
    @property
    def pmin(self):
        return np.min(self._plist)
    @property
    def pmax(self):
        return np.max(self._plist)
    @property
    def Tstep(self):
        return np.min(np.diff(np.unique(np.sort(self._Tlist))))
    @property
    def pstep(self):
        return np.min(np.diff(np.unique(np.sort(self._plist))))
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self,value):
        self._name = value

    ## operations
    def __mul__(self,other):
        """Multiply by scalar: all extensive properties are
        multiplied, all intensive properties are left alone. Update
        the name."""
        ph = copy.deepcopy(self)
        ph._Vlist *= other
        ph._Glist *= other
        ph._Ginterp = CloughTocher2DInterpolator(list(zip(ph._plist, ph._Tlist)), ph._Glist)
        ph._Vinterp = CloughTocher2DInterpolator(list(zip(ph._plist, ph._Tlist)), ph._Vlist)

        if len(self.name.split()) > 1:
            sname = f"({self.name})"
        else:
            sname = self.name
        ph._name = f"{other} * {sname}"
        return ph
    __rmul__ = __mul__

    def __add__(self,other):
        """Add two phases. The pressure range is the intersection of
        the pressure ranges and the extensive properties are added."""
        if other == 0:
            return self
        else:
            pmin = np.min([np.min(self._plist),np.min(other._plist)])
            pmax = np.max([np.max(self._plist),np.max(other._plist)])
            pstep = np.min([self.pstep,other.pstep])
            Tmin = np.min([np.min(self._Tlist),np.min(other._Tlist)])
            Tmax = np.max([np.max(self._Tlist),np.max(other._Tlist)])
            Tstep = np.min([self.Tstep,other.Tstep])
            p = np.arange(pmin,pmax+pstep,pstep)
            T = np.arange(Tmin,Tmax+Tstep,Tstep)
            pm,Tm = np.meshgrid(p,T)
            pm = pm.flatten()
            Tm = Tm.flatten()

            ## fill the Gibbs energy and volume
            G = self._Ginterp((pm,Tm)) + other._Ginterp((pm,Tm))
            V = self._Vinterp((pm,Tm)) + other._Vinterp((pm,Tm))

            ## remove the nans
            mask = ~(np.isnan(G) | np.isnan(V))
            pm = pm[mask]
            Tm = Tm[mask]
            V = V[mask]
            G = G[mask]

            if len(self.name.split()) > 1:
                sname = f"({self.name})"
            else:
                sname = self.name
            if len(other.name.split()) > 1:
                oname = f"({other.name})"
            else:
                oname = other.name

            return Phase(f"{sname} + {oname}",pm,Tm,V,G)
    __radd__ = __add__

    def __or__(self,other):
        """The operation A | B builds the thermodynamically stable
        phase in the union (p,T) range of A and B."""

        pmin = np.min([np.min(self._plist),np.min(other._plist)])
        pmax = np.max([np.max(self._plist),np.max(other._plist)])
        pstep = np.min([self.pstep,other.pstep])
        Tmin = np.min([np.min(self._Tlist),np.min(other._Tlist)])
        Tmax = np.max([np.max(self._Tlist),np.max(other._Tlist)])
        Tstep = np.min([self.Tstep,other.Tstep])
        p = np.arange(pmin,pmax+pstep,pstep)
        T = np.arange(Tmin,Tmax+Tstep,Tstep)
        pm,Tm = np.meshgrid(p,T)
        pm = pm.flatten()
        Tm = Tm.flatten()

        ## fill the Gibbs energy
        Gs = self._Ginterp((pm,Tm))
        Go = other._Ginterp((pm,Tm))
        G = np.fmin(Gs,Go)

        ## fill the others with a mask
        V = np.empty_like(G)
        mask = (G == Gs)
        V[mask] = self.V(pm[mask],Tm[mask])
        mask = ~mask
        V[mask] = other.V(pm[mask],Tm[mask])

        ## remove the nans
        mask = ~(np.isnan(G) | np.isnan(V))
        pm = pm[mask]
        Tm = Tm[mask]
        V = V[mask]
        G = G[mask]
        if len(self.name.split()) > 1:
            sname = f"({self.name})"
        else:
            sname = self.name
        if len(other.name.split()) > 1:
            oname = f"({other.name})"
        else:
            oname = other.name

        return Phase(f"{sname} | {oname}",pm,Tm,V,G)
    __ror__ = __or__

    ## representation functions
    def __str__(self):
        return f"Phase {self._name} with {len(self._plist)} pressures and {len(self._Tlist)} temperatures."

def read_phases_from_eos_file(eosfile):
    """Read a .eos file generated by gibbs2 and return a list of phases containing
    all calculated thermodynamic properties."""

    phlist = []
    name = None
    with open(eosfile,"r") as f:
        for line in f:
            if line.startswith("# Phase"):
                if name is not None:
                    phlist.append(Phase.fromlines(name,strl))
                name = line.split()[-1]
                strl = []
            if name is not None and not line.startswith("#"):
                strl.append(line)
    if name is not None:
        phlist.append(Phase.fromlines(name,strl))

    return phlist
