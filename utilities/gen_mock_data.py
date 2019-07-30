"""
A script to generate mock data. Creates the appropriate directory structure
and HDF5 files to mimic GADGET output, essentially creating random density,
velocity and temperature fields. 

Note: care should be taken when using this script!
The generated data can easily become VERY large on disk!
e.g. with N = 512, total size of snapshot on disk is ~6 GB
     with N = 2048, total size on disk is ~380 GB

Note: if generating a large mock data set, spread the data over more files
to minimise memory footprint during generation. 
e.g. with N = 512, use NUMFILES = 4 (each file ~1.5 GB)
     with N = 2048, use NUMFILES = 256 (each file ~1.5 GB)
"""

import numpy as np
import h5py as h5
import os


# PARAMETERS TO TWEAK
N        = 512                # Cube root number of particles
BOXLEN   = 10000.0            # Box side length (in ckpc/h) (NB: float not int)
NUMFILES = 4                  # Split the simulation over multiple files
SNAPNUM  = 15                 # Mock snapshot number
PATH     = "../"              # Directory to save mock data into


################################################################################
# DO NOT EDIT BELOW THIS LINE ##################################################
################################################################################


# GENERATE DATA
N3       = N**3               # Total number of particles
NperFile = int(N3/NUMFILES)   # Number of particles per file
REDSHIFT = 5.756              # Cosmological parameter
OMEGA0   = 0.308              # Cosmological parameter
OMEGAL   = 0.692              # Cosmological parameter
SMALLH   = 0.678              # Cosmological parameter
ATIME    = 1/(1+REDSHIFT)     # Cosmological parameter

# FUNCTION TO SAVE DATA TO HDF5 FILE
def gen_sub_file(parent_dir, snapnum, nfile):
    fname = parent_dir + "/snap_{:03d}.{:d}.hdf5".format(snapnum, nfile)
    print("Saving to: ", fname)
    
    # GENERATE RANDOM VALUES
    coords = BOXLEN*np.random.rand(NperFile,3)
    masses = 9.967319e-6*np.ones(shape=NperFile)
    vels   = 150*np.random.rand(NperFile,3)
    U      = 250*np.random.rand(NperFile)
    h      = 30.0*np.random.rand(NperFile)+50.0
    Ne     = 1.08*np.ones(shape=NperFile)
    NHI    = np.random.rand(NperFile)

    # MIX SIGNS FOR VELOCITIES
    signs  = np.random.choice([-1,1], size=NperFile*3).reshape(NperFile,3)
    vels   = np.multiply(signs,vels)

    # CONVERT TO CORRECT DATATYPE
    coords = coords.astype(np.float32)
    masses = masses.astype(np.float32)
    vels   = vels.astype(np.float32)
    U      = U.astype(np.float32)
    h      = h.astype(np.float32)
    Ne     = Ne.astype(np.float32)
    NHI    = NHI.astype(np.float32)

    # CREATE HDF5 FILE
    with h5.File(fname, 'w') as f:
        h5header = f.create_group("Header")
        h5part0  = f.create_group("PartType0")

        h5header.attrs['NumFilesPerSnapshot'] = NUMFILES
        h5header.attrs['BoxSize']             = BOXLEN 
        h5header.attrs['Redshift']            = REDSHIFT 
        h5header.attrs['Omega0']              = OMEGA0 
        h5header.attrs['OmegaLambda']         = OMEGAL 
        h5header.attrs['HubbleParam']         = SMALLH
        h5header.attrs['Time']                = ATIME
        h5header.attrs['NumPart_ThisFile']    = NperFile

        h5part0.create_dataset("Coordinates", data=coords)
        h5part0.create_dataset("Velocities", data=vels)
        h5part0.create_dataset("Masses", data=masses)
        h5part0.create_dataset("InternalEnergy", data=U)
        h5part0.create_dataset("SmoothingLength", data=h)
        h5part0.create_dataset("ElectronAbundance", data=Ne)
        h5part0.create_dataset("NeutralHydrogenAbundance", data=NHI)


# CREATE DIRECTORIES AND SAVE DATA TO FILES
print("Generating mock data...")
print("Particles per file: ", NperFile)
print("Estimated size on disk: ~{:.3f} GB".format(4*11*N3*1e-9))

directory = PATH + "mock_data"
if not os.path.exists(directory):
    os.makedirs(directory)

subdirectory = directory + "/snapdir_{:03d}".format(SNAPNUM)
if not os.path.exists(subdirectory):
    os.makedirs(subdirectory)

for i in range(NUMFILES):
    gen_sub_file(subdirectory, SNAPNUM, i)

print("Finished!")
print(os.listdir(PATH))
print(os.listdir(directory))
print(os.listdir(subdirectory))
