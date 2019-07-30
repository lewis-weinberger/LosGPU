"""
Plot the time taken for Arepo projection, comparing GPU and CPU implementations.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter
mpl.rcParams["axes.linewidth"] = 1.25
mpl.rcParams["text.usetex"]    = True
mpl.rcParams["font.family"]    = "DejaVu Sans"
mpl.rcParams["font.serif"]     = "cm"

# OPTIONS
dark_theme = False

if dark_theme:
    # DARK THEME
    mpl.rcParams["savefig.facecolor"] = "#000000"
    mpl.rcParams["savefig.edgecolor"] = "#000000"
    mpl.rcParams["text.color" ]      = "#FFFFFF"
    mpl.rcParams["axes.labelcolor" ] = "#FFFFFF"
    mpl.rcParams["axes.facecolor"]   = "#000000"
    mpl.rcParams["axes.edgecolor"]   = "#FFFFFF"
    mpl.rcParams["xtick.color" ]     = "#FFFFFF"
    mpl.rcParams["ytick.color" ]     = "#FFFFFF"

f = "times.txt"
data = np.loadtxt(f, unpack=True)
CPU_mask = data[2] == 0
GPU_mask = data[2] == 1

res_CPU = data[0,CPU_mask]
t_CPU   = data[1,CPU_mask]/3600
res_GPU = data[0,GPU_mask]
t_GPU   = data[1,GPU_mask]/3600

fig = plt.figure(figsize=(8, 8))
ax  = fig.add_axes([0, 0, 0.4, 0.4])

ax.plot(res_CPU, t_CPU, c="coral", ls="-", marker="o", lw=2.5, label='LosGPU: 1 CPU',
        markersize=8, markeredgecolor="none")
ax.plot(res_GPU, t_GPU, c="lightseagreen", ls="-", marker="o", lw=2.5, 
        label='LosGPU: 1 CPU, 1 GPU', markersize=8, markeredgecolor="none")

ax.set_xlabel("Number of sightlines", labelpad=10)
ax.set_ylabel("Extraction Time [hr]", labelpad=10)

ax.legend(loc='upper left', frameon=False, numpoints=1, fontsize=12)

#ax.set_yscale('log')
ax.set_xscale('log')
ax.xaxis.set_major_formatter(ScalarFormatter())
#ax.yaxis.set_major_formatter(ScalarFormatter())

bufx = 1.15
bufy = 1.15
ax.set_xlim(res_CPU.min()/bufx, res_CPU.max()*bufx)
ax.set_ylim(-0.1, t_CPU.max()*bufy)

plt.savefig("time_comparison.pdf", bbox_inches='tight')
