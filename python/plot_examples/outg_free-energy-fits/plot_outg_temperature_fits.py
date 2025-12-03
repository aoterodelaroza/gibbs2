import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

## Make a plot of the Fvib(V;T) for all phases at the given
## temperature from the info in the gibbs2 output file.

####
gibbs2out="pirssonite.outg" ## output file containing the data
tref=300 ## temperature to plot
iplot=3 ## plot column (1=F,2=Fvib,3=Fvib-F0,4=S,5=Cv)
####

pnames = []
data = []

inblock = False
name = None
with open(gibbs2out,"r") as fin:
    for line in fin:
        if (line.startswith("* Thermodynamic properties")):
            inblock = True
        if inblock and line.startswith("+ Phase"):
            name = line.split()[-1]
        if inblock and name and line.startswith("# Temperature ="):
            if (np.abs(float(line.split()[-1]) - tref) < 1e-2):
                strlist = []
                fin.readline()
                fin.readline()
                for line2 in fin:
                    if line2.startswith("#"):
                        break
                    strlist.append(line2)

                stream = StringIO(''.join(strlist))
                data.append(np.genfromtxt(stream))
                pnames.append(name)
                name = None
        if len(line) == 0:
            inblock = False

fig, ax = plt.subplots()
labels=["V(A^3)","F(Ha)","Fvib(Ha)","Fvib-F0(Ha)","S(Ha/K)","CV(Ha/K)"]

for idx, ff in enumerate(pnames):
    xx = data[idx]
    v = xx[:,0] * 0.52917720859**3
    y = xx[:,2]
    ax.plot(v,y,'-o',label=ff)
ax.set_xlabel(labels[0])
ax.set_ylabel(labels[iplot])
ax.grid()
ax.legend()
fig.tight_layout()
fig.savefig("plot_outg.pdf")

