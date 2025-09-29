## append the gibbs2 python module to the path
import sys
sys.path.append('/home/alberto/git/gibbs2/python')

##
import gibbs2 as g2
import matplotlib.pyplot as plt

phlist = g2.read_phases_from_eos_file("silica.eos")
fig,ax = plt.subplots()

fig,ax = g2.plot_phase_diagram(fig,ax,phlist,steprefine=10)

ax.set_ylim([0,2500])

fig.tight_layout()
fig.savefig("phase_diagram.png")
