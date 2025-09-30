## append the gibbs2 python module to the path
import sys
sys.path.append('/home/alberto/git/gibbs2/python')

##
import gibbs2 as g2
import matplotlib.pyplot as plt

phlist = g2.read_staticphases_from_eos_file("till.eos_static")

fig,ax = plt.subplots()
fig,ax = g2.plot_enthalpy_pressure(fig,ax,phlist,steprefine=10)
fig,ax = g2.plot_barh_stablephase(fig,ax,phlist,-29,steprefine=10)

ax.set_xlim([0,30])

fig.tight_layout()
fig.savefig("enthalpy_pressure.pdf")

