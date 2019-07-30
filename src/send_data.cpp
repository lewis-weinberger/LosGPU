/* Copyright 2018 Lewis Weinberger.                                           */
/* GPU accelerate sightline extraction for Sherwood (GADGET) snapshots.       */

#include "mpi.h"
#include "./global_vars.h"

void send_header(MPI_Comm comm) {
  /* Send data from header file */
  MPI_Bcast(&BOXLEN, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&REDSHIFT, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&OMEGA0, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&OMEGALAMBDA, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&SMALLH, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&ATIME, 1, MPI_DOUBLE, 0, comm);
}

void send_params(MPI_Comm comm, int rank) {
  /* Send data from parameter file */
  MPI_Bcast(&NUMLOS, 1, MPI_INT, 0, comm);
  MPI_Bcast(&NBINS, 1, MPI_INT, 0, comm);
  MPI_Bcast(&NUMSEG_TOT, 1, MPI_INT, 0, comm);

  /* Allocate memory on non-root processes */
  if (rank != 0 ) {
    NUMSEG      = new int[NUMLOS];
    NUMSEG_PREV = new int[NUMLOS];
    xseg        = new double[NUMSEG_TOT];
    yseg        = new double[NUMSEG_TOT];
    zseg        = new double[NUMSEG_TOT];
    xlos_mid    = new double[NUMLOS];
    ylos_mid    = new double[NUMLOS];
    zlos_mid    = new double[NUMLOS];
    xlos        = new double[NUMLOS*NBINS];
    ylos        = new double[NUMLOS*NBINS];
    zlos        = new double[NUMLOS*NBINS];
    rhoker      = new double[NUMLOS*NBINS];
    tempker     = new double[NUMLOS*NBINS];
    velker      = new double[NUMLOS*NBINS];
  }

  MPI_Bcast(NUMSEG, NUMLOS, MPI_INT, 0, comm);
  MPI_Bcast(NUMSEG_PREV, NUMLOS, MPI_INT, 0, comm);
  MPI_Bcast(xseg, NUMSEG_TOT, MPI_DOUBLE, 0, comm);
  MPI_Bcast(yseg, NUMSEG_TOT, MPI_DOUBLE, 0, comm);
  MPI_Bcast(zseg, NUMSEG_TOT, MPI_DOUBLE, 0, comm);
  
  MPI_Bcast(&LENLOS, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(&drbin, 1, MPI_DOUBLE, 0, comm);

  MPI_Bcast(dir, 3, MPI_DOUBLE, 0, comm);
  MPI_Bcast(xlos_mid, NUMLOS, MPI_DOUBLE, 0, comm);
  MPI_Bcast(ylos_mid, NUMLOS, MPI_DOUBLE, 0, comm);
  MPI_Bcast(zlos_mid, NUMLOS, MPI_DOUBLE, 0, comm);
  
  MPI_Bcast(xlos, NUMLOS*NBINS, MPI_DOUBLE, 0, comm);
  MPI_Bcast(ylos, NUMLOS*NBINS, MPI_DOUBLE, 0, comm);
  MPI_Bcast(zlos, NUMLOS*NBINS, MPI_DOUBLE, 0, comm);
  
  MPI_Bcast(rhoker, NUMLOS*NBINS, MPI_DOUBLE, 0, comm);
  MPI_Bcast(velker, NUMLOS*NBINS, MPI_DOUBLE, 0, comm);
  MPI_Bcast(tempker, NUMLOS*NBINS, MPI_DOUBLE, 0, comm);
}

void send_data(MPI_Comm comm, int rank) {
  /* Send particle data */
  int typesize;
  int count = 7;
  int blocklengths[7] = {3, 3, 1, 1, 1, 1, 1};
  MPI_Aint displacements[7] = {0, 12, 24, 28, 32, 36, 40};
  MPI_Datatype types[7] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
                   MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
  MPI_Datatype P_data_type;
  MPI_Type_create_struct(count, blocklengths, displacements, types, 
                         &P_data_type);
  MPI_Type_commit(&P_data_type);
  MPI_Type_size(P_data_type, &typesize);

  MPI_Bcast(&NUMPARTICLES, 1, MPI_INT, 0, comm);
  
  /* Allocate memory on non-root processes */
  if (rank != 0) 
    P = new particle_data[NUMPARTICLES];

  MPI_Bcast(P, NUMPARTICLES, P_data_type, 0, comm);
}
