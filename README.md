### LosGPU
**LosGPU** is a GPU accelerated sightline extraction code designed to be run on the [Sherwood](https://www.nottingham.ac.uk/astronomy/sherwood/) simulation suite. 
Since Sherwood is a [GADGET](https://wwwmpa.mpa-garching.mpg.de/gadget/)-based simulation, LosGPU can be used on other GADGET simulations (possibly requiring some tweaking to use appropriate particle fields).

#### Requirements

**LosGPU** is written in C++ using OpenACC directives, so it requires a compiler such as PGI's `pgc++`. It also requires an [MPI](https://www.mpi-forum.org/) library, and the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) library.

#### Installation:
To compile the code simply `cd` into the `src` directory and run `make`. This will create an executable **LosGPU** in the main directory. The Makefile uses the HDF5 compiler wrapper (`h5c++`) which in turn uses the MPI wrapper (`mpicxx`) which finally should use the PGI compiler (`pgc++`) in order to take advantage of the OpenACC directives for acceleration. If you do not have the PGI compiler, the Makefile will need to be tweaked. Note: although `gcc` supports some OpenACC directives, their coverage is still not extensive, so using PGI's compiler is highly recommended.

The provided Makefile has three target configuration options:

1. `gpu`: targeting GPU acceleration,
2. `cpu`: targeting multicore acceleration\*,
3. `test`: for testing (default option).

If you run `make` or `make test` then the default setup is to compile for the "test" configuration. 
Otherwise to compile for the accelerated configurations run either `make gpu` or `make cpu`. 
To get rid of the object files and the exectuble, run `make clean`.

\* Although `multicore` does compile, I haven't been able to run it successfully. I get an error about `MP_BLIST`, which is an environment variable that the PGI runtime uses to determine the binding between OpenMP threads and physical cores.

#### Usage:
The executable takes 3 command line arguments:

1. Simulation directory path (one level up from snapdir\_XXX),
2. Snapshot number,
3. Parameter file.

The parameter file determines what sightlines will be extracted from the simulation volume. Some example versions are provided in **example\_params**. The format of the parameter file is as follows:

- Line 0: comment line (ignored)
- Line 1: NUMLOS NBINS LENLOS/[cMpc/h] x\_dir y\_dir z\_dir output\_path
- Line 2: comment line (ignore)
- Remaining lines: x\_mid/[cMpc/h] y\_mid/[cMpc/h] z\_mid/[cMpc/h]

For example, to extract sightlines from a simulation at `/path/to/sim/`, snapshot number 15, with the parameter file **example.params**, run the following:

`mpirun -np 1 ./LosGPU /path/to/sim 15 example.params`

Information about the simulation and the sightline parameters will be printed during execution, so you can check it's running properly.

#### Output
Once extracted, the sightlines are outputted to a binary formatted file. This can be read with the provided Python utility script **read\_los.py**.

### Utility Scripts
There are some Python 3 utility scripts for automating some processes like generating parameter files or mock data. These include:

1. `gen_mock_data.py`: a script for generating mock data to test the code on,
2. `gen_test.py`: a script to generate a parameter file with random coordinates,
3. `read_los.py`: a script to read the binary output file.

I haven't tested if these run properly with Python 2.

---

### A NOTE ON PARALLELIZATION
The code is currently parallelized in two layers. Firstly, for each sightline the code needs to loop over all the particles to add their (possible) contributions. However since this is only additive, it doesn't matter what order it is done, so we can divide the particles over different processes and then reduce the final result at the end. This level of parallelism is achieved with MPI, where each MPI task handles some subset of all the particles, and we do a reduction at the end to get the full sightlines. This level of parallelization affects the time taken to calculate each individual sightline. 

The second level of parallelization is across the different sightlines. Note these are all independent of each other, and so can be run embarassingly parallel. This is accelerated using the GPU, where the many "cores" of a GPU can handle individual sightlines. For example, modern Nvidia GPUs like the Tesla P100 have thousands of "CUDA cores", and so they can calculate many many sightlines at once. When combined, the two parallelisms give a massive speed-up compared to running in serial.

Disk IO is performed in serial by the root MPI task, and data is then broadcasted to other tasks. This IO scheme is used to prevent file-system bottlenecks; if the application is run with many MPI tasks trying to read files on disk simulataneously, this could cause issues. For a small number of tasks the safer scheme may be very slightly slower depending on the hardware/filesystem used, I haven't checked.
