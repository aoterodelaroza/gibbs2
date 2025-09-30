import numpy as np
from scipy.interpolate import interp1d
import copy

class StaticPhase:
    """Thermodynamic properties of a static phase calculated by
    gibbs2. This class is derived from the information written to a
    .eos_static file and is used to create plots and performing tasks
    not done in gibbs2 by default."""

    ## constructors
    def __init__(self,name,p,E,H,V):
        """Initialize staticphase given name, pressure, energy,
        enthalpy, and volume."""

        ## assign fields
        self._name = name
        self._plist = np.copy(p)
        self._Elist = np.copy(E)
        self._Hlist = np.copy(H)
        self._Vlist = np.copy(V)

        ## unique temperatures and pressures to within 2 decimal places
        self._pkeys = np.unique(np.sort(np.round(self._plist,decimals=2)))

        ## set up the interpolant
        self._Einterp = interp1d(self._plist,self._Elist,'cubic',bounds_error=False,fill_value=np.nan)
        self._Hinterp = interp1d(self._plist,self._Hlist,'cubic',bounds_error=False,fill_value=np.nan)
        self._Vinterp = interp1d(self._plist,self._Vlist,'cubic',bounds_error=False,fill_value=np.nan)

    @classmethod
    def fromlines(cls,name,strl):
        """Initialize staticphase from a list of strings with the
        thermodynamic properties obtained from reading the eos_static
        file (see examples)."""

        x = np.genfromtxt(strl)
        p = x[:,0]
        E = x[:,1]
        H = x[:,2]
        V = x[:,3]
        return cls(name,p,E,H,V)

    @classmethod
    def fromfile(cls,filename):
        """Initialize staticphase from an eos_static file. If several
        phases are present, use the first one."""

        return read_staticphases_from_eos_file(filename)[0]

    ## interpolated thermodynamic properties
    def H(self,p):
        """Returns the pressure-interpolated static enthalpy."""
        return self._Hinterp(p)
    def E(self,p):
        """Returns the pressure-interpolated static energy."""
        return self._Einterp(p)
    def V(self,p):
        """Returns the pressure-interpolated volume."""
        return self._Vinterp(p)

    ## min, max, step temperature and pressure
    @property
    def pmin(self):
        return np.min(self._plist)
    @property
    def pmax(self):
        return np.max(self._plist)
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
        ph._Elist *= other
        ph._Hlist *= other
        ph._Einterp = interp1d(ph._plist,ph._Elist,'cubic',bounds_error=False,fill_value=np.nan)
        ph._Hinterp = interp1d(ph._plist,ph._Hlist,'cubic',bounds_error=False,fill_value=np.nan)
        ph._Vinterp = interp1d(ph._plist,ph._Vlist,'cubic',bounds_error=False,fill_value=np.nan)

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
            pmin = max(self.pmin,other.pmin)
            pmax = min(self.pmax,other.pmax)
            pstep = min(self.pstep,other.pstep)
            p = np.arange(other.pmin,other.pmax,other.pstep)
            E = self.E(p) + other.E(p)
            H = self.H(p) + other.H(p)
            V = self.V(p) + other.V(p)

            if len(self.name.split()) > 1:
                sname = f"({self.name})"
            else:
                sname = self.name
            if len(other.name.split()) > 1:
                oname = f"({other.name})"
            else:
                oname = other.name

            return StaticPhase(f"{sname} + {oname}",p,E,H,V)
    __radd__ = __add__

    def __and__(self,other):
        """The operation A & B builds the thermodynamically stable
        phase in the intersection pressure range of A and B."""
        pmin = max(self.pmin,other.pmin)
        pmax = min(self.pmax,other.pmax)
        pstep = min(self.pstep,other.pstep)
        p = np.arange(other.pmin,other.pmax,other.pstep)

        ## fill the enthalpy
        Hs = self.H(p)
        Ho = other.H(p)
        H = np.minimum(Hs,Ho)

        ## fill the others with a mask
        V = np.empty_like(H)
        E = np.empty_like(H)
        mask = (H == Hs)
        V[mask] = self.V(p[mask])
        E[mask] = self.E(p[mask])
        mask = ~mask
        V[mask] = other.V(p[mask])
        E[mask] = other.E(p[mask])

        if len(self.name.split()) > 1:
            sname = f"({self.name})"
        else:
            sname = self.name
        if len(other.name.split()) > 1:
            oname = f"({other.name})"
        else:
            oname = other.name

        return StaticPhase(f"{sname} & {oname}",p,E,H,V)

    __rand__ = __and__


    ## representation functions
    def __str__(self):
        return f"Phase {self._name} with {len(self._pkeys)} pressures ({np.min(self._plist)} GPa -> {np.max(self._plist)} GPa)."

def read_staticphases_from_eos_file(eosfile):
    """Read a .eos_static file generated by gibbs2 and return a list
    of phases containing all static properties."""

    phlist = []
    name = None
    reading = False
    strl = []
    with open(eosfile,"r") as f:
        for line in f:
            if line.startswith("#") and name is not None and reading:
                phlist.append(StaticPhase.fromlines(name,strl))
                reading = False
                name = None
                strl = []
            if line.startswith("# Phase"):
                name = line.split()[-1]
            if name is not None and not line.startswith("#") and len(line) > 0:
                reading = True
                strl.append(line)
    if name is not None:
        phlist.append(StaticPhase.fromlines(name,strl))

    return phlist
