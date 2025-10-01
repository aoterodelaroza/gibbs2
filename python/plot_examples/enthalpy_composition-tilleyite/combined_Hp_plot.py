## append the gibbs2 python module to the path
import sys
sys.path.append('/home/alberto/git/gibbs2/python')

##
import gibbs2 as g2
import matplotlib.pyplot as plt

co2I   = g2.StaticPhase.fromfile("co2-I.eos_static")
co2II  = g2.StaticPhase.fromfile("co2-II.eos_static")
co2III = g2.StaticPhase.fromfile("co2-III.eos_static")
co2V   = g2.StaticPhase.fromfile("co2-V.eos_static")
co2 = (co2I & co2II & co2III & co2V)
co2.name = r"CO$_2$"

cao = g2.StaticPhase.fromfile("cao-B1.eos_static")
cao.name = r"CaO"

coesite    = g2.StaticPhase.fromfile("coesite.eos_static")
quartz     = g2.StaticPhase.fromfile("quartz.eos_static")
stishovite = g2.StaticPhase.fromfile("stishovite.eos_static")
sio2 = (coesite & quartz & stishovite)
sio2.name = r"SiO$_2$"

calcite       = g2.StaticPhase.fromfile("calcite.eos_static")
aragonite     = g2.StaticPhase.fromfile("aragonite.eos_static")
caco3VI       = g2.StaticPhase.fromfile("caco3-VI.eos_static")
caco3VII      = g2.StaticPhase.fromfile("caco3-VII.eos_static")
postaragonite = g2.StaticPhase.fromfile("post-aragonite.eos_static")
caco3 = (calcite & aragonite & caco3VI & caco3VII & postaragonite)
caco3.name = r"CaCO$_3$"

ca2sico3 = g2.StaticPhase.fromfile("ca2sico3.eos_static")
ca2sico3.name = r"Ca$_2$Si(CO$_3$)$_4$"

ptillhi = g2.StaticPhase.fromfile("post-till-hi.eos_static")
ptilllo = g2.StaticPhase.fromfile("post-till-lo.eos_static")
tillhi  = g2.StaticPhase.fromfile("till-hi.eos_static")
tilllo  = g2.StaticPhase.fromfile("till-lo.eos_static")
till = (ptillhi | ptilllo | tillhi | tilllo)
till.name = "Tilleyite"

phlist = [till + 6 * co2,
          2 * ca2sico3 + cao,
          5 * cao + 2 * sio2 + 8 * co2,
          5 * caco3 + 2 * sio2 + 3 * co2,
          ]
fig,ax = plt.subplots()
fig,ax = g2.plot_enthalpy_pressure(fig,ax,phlist,steprefine=10,idref=0)
ax.set_xlim([0,30])
ax.set_title(r"Relative enthalpy per Ca$_5$Si$_2$C$_8$O$_{25}$ formula unit")
fig.tight_layout()
fig.savefig("composition.pdf")

