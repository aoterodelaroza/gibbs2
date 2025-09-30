## append the gibbs2 python module to the path
import sys
sys.path.append('/home/alberto/git/gibbs2/python')

##
import gibbs2 as g2
import matplotlib.pyplot as plt

calcite   = g2.StaticPhase.fromfile("calcite.eos_static")
aragonite = g2.StaticPhase.fromfile("aragonite.eos_static")
cao       = g2.StaticPhase.fromfile("cao-B1.eos_static")
co2I      = g2.StaticPhase.fromfile("co2-I.eos_static")
co2II     = g2.StaticPhase.fromfile("co2-II.eos_static")
co2III    = g2.StaticPhase.fromfile("co2-III.eos_static")
co2V      = g2.StaticPhase.fromfile("co2-V.eos_static")

caco3 = calcite & aragonite
caco3.name = "CaCO3"

cao_co2 = cao + (co2I & co2II & co2III & co2V)
cao_co2.name = "CaO + CO2"

fig,ax = plt.subplots()
fig,ax = g2.plot_enthalpy_pressure(fig,ax,[caco3,cao_co2],steprefine=10)
fig,ax = g2.plot_barh_stablephase(fig,ax,[caco3,cao_co2],-29,steprefine=10)

fig.tight_layout()
fig.savefig("enthalpy_phase-composition.pdf")
