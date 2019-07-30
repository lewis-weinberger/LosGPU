"""
Generates a parameter file for random sightlines.
"""
import numpy as np

BOXLEN    = 10 # cMpc/h
NUMLOS    = 100
NBINS     = 2048
LENLOS    = 10  # cMpc/h
direction = (1, 0, 1)
output    = "output_{}.los".format(NUMLOS)

x = np.random.rand(NUMLOS)*BOXLEN
y = np.random.rand(NUMLOS)*BOXLEN
z = np.random.rand(NUMLOS)*BOXLEN

fname = "new_example.params"
with open(fname, "w") as f:
    f.write("# NUMLOS NBINS LENLOS/[cMpc/h] x_dir y_dir z_dir output_path\n")
    f.write("{}    {}    {:.3f}    {:.3f}    {:.3f}    {:.3f}   {}\n".format(NUMLOS,
           NBINS, LENLOS, direction[0], direction[1], direction[2], output))
    f.write("# x_mid/[cMpc/h] y_mid/[cMpc/h] z_mid/[cMpc/h]\n")
    for i in range(NUMLOS):
        f.write("{:.3f}    {:.3f}    {:.3f}\n".format(x[i], y[i], z[i]))
