/* Copyright 2018 Lewis Weinberger.                                           */
/* GPU accelerated line-of-sight extraction for Sherwood (GADGET) snapshots.  */

#if defined(_OPENACC)
#include "openacc.h"
#endif
#include "mpi.h"
#include <math.h>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <sstream>

#include "./prototypes.h"
#include "./global_vars.h"

using std::cerr;
using std::cout;
using std::endl;
using std::invalid_argument;
using std::string;

double *frho = nullptr, *fvel = nullptr, *ftemp = nullptr;

int main(int argc, char **argv) {

  /* MPI variables */
  MPI_Comm comm;
  int ntasks, rank;
  int *npart, *rnpart;

  /* Initialise MPI */
  MPI_Init(NULL, NULL);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &ntasks);
  MPI_Comm_rank(comm, &rank);
  npart  = new int[ntasks]; // number of particles each task processes
  rnpart = new int[ntasks]; // cumulative number across processes

  /* Process runtime arguments */
  if (argc != 4) {
    if (rank == 0) {
      cerr << "Incorrect argument(s)." << endl;
      cerr << "Usage: "
           << argv[0]
           << " [snapshot path] [snapshot number] [parameter file]"
           << endl;
    }
    MPI_Finalize();
    return 1;
  }

  /* Start time */
  double t_start = MPI_Wtime();

  /* Print fancy figlet title */
  if (rank == 0)
    print_title();

  string snap_dir(argv[1]);
  int snapshot = std::stoi(argv[2]);
  string param_path(argv[3]);

  /* Process snapshot file path */
  std::stringstream path;
  path << snap_dir << "snapdir_" << std::setfill('0') << std::setw(3)
       << snapshot << "/snap_"   << std::setfill('0') << std::setw(3)
       << snapshot;
  string snap_path = path.str();

  /* Load in simulation snapshot header info on root process */
  int rerr = 0;
  if (rank == 0)
    cout << endl << "Reading snapshot header file: " << snap_dir << endl;
    try {
      read_snapshot_header(snap_path, rank);
    } catch (invalid_argument err) {
        cout << err.what();
        rerr = 1;
      }
  MPI_Bcast(&rerr, 1, MPI_INT, 0, comm);
  if (rerr == 1) {
    MPI_Finalize();
    return 1;
  }
  /* Send data to other processes */
  send_header(comm);

  /* Read sightline parameter file on root process */
  if (rank == 0)
    cout << endl << "Reading parameter file: " << param_path << endl;
    try {
      read_params(param_path, rank);
    } catch (invalid_argument err) {
        cout << err.what();
        rerr = 1;
      }
  MPI_Bcast(&rerr, 1, MPI_INT, 0, comm);
  if (rerr == 1) {
    MPI_Finalize();
    return 1;
  }
  /* Send data to other processes */
  send_params(comm, rank);

  #if defined(_OPENACC)
  /* Associate GPU devices with MPI tasks */
  int ngpus = acc_get_num_devices(acc_device_nvidia);
  int gpunum;
  if (ngpus == ntasks) {
    gpunum = rank % ntasks;
    acc_set_device_num(gpunum, acc_device_nvidia);
  } else {
      /* If there aren't as many GPUs as MPI tasks, then don't use GPUs */
      acc_set_device_type(acc_device_host);
    }
  #endif

  /* Timing statistics */
  double t_init = MPI_Wtime() - t_start;
  double t_file = 0.0, t_extract = 0.0;

  /* GPU initialization */
  #if defined(_OPENACC)
  #pragma acc enter data copyin(NBINS, drbin, NUMLOS,    \
                                BOXLEN, SMALLH, ATIME,   \
                                dir[0:3],                \
                                xlos[0:NUMLOS*NBINS],    \
                                ylos[0:NUMLOS*NBINS],    \
                                zlos[0:NUMLOS*NBINS],    \
                                NUMSEG[0:NUMLOS],        \
                                NUMSEG_PREV[0:NUMLOS],   \
                                xseg[0:NUMSEG_TOT],      \
                                yseg[0:NUMSEG_TOT],      \
                                zseg[0:NUMSEG_TOT],      \
                                rhoker[0:NUMLOS*NBINS],  \
                                velker[0:NUMLOS*NBINS],  \
                                tempker[0:NUMLOS*NBINS])
  #endif
    
  /* Loop over snapshot files */
  for (int i = 0; i < NUMFILES; i++) {

    /* Load snapshot data on root process*/
    double t_f0 = MPI_Wtime();
    if (rank == 0)
      read_snapshot(snap_path, i, rank);

    /* Send data to other processes */
    send_data(comm, rank);
    t_file += MPI_Wtime() - t_f0;

    /* Decomposition of particles across MPI tasks (and therefore GPUs) */
    int tmp  = 0;
    rnpart[0] = 0;
    for (int n = 0; n < ntasks; n++) {
      if (n < ntasks - 1)
        npart[n] = NUMPARTICLES/ntasks;
      else 
        npart[n] = NUMPARTICLES - (ntasks - 1)*(NUMPARTICLES/ntasks);

      if (n > 0) {
        tmp += npart[n - 1];
        rnpart[n] = tmp;
      }
    }

    /* Print decomposition info */
    if (rank == 0) {
      cout << "Using " << ntasks << " MPI task(s)." << endl;
      #if defined(_OPENACC)
      cout << "Using " << ngpus << " GPU device(s)." << endl;
      #endif
      cout << "Decomposing " << NUMPARTICLES << " particles:"  << endl;
      cout << "npart = {";
      for (int n = 0; n < ntasks; n++) {
        if (n == ntasks - 1)
          cout << npart[n];
        else  
          cout << npart[n] << ", ";
      }
      cout << "}" << endl;
      cout << "rnpart = {";
      for (int n = 0; n < ntasks; n++) {
        if (n == ntasks - 1)
          cout << rnpart[n];
        else  
          cout << rnpart[n] << ", ";
      }
      cout << "}" << endl;
      cout << "Extracting sightlines from file..." << endl;
    }
    
    double t_e0 = MPI_Wtime();

    /* Parallel kernels region */
    #if defined(_OPENACC)
    #pragma acc kernels present(NBINS, drbin, NUMLOS,    \
                                BOXLEN, SMALLH, ATIME,   \
                                dir[0:3],                \
                                xlos[0:NUMLOS*NBINS],    \
                                ylos[0:NUMLOS*NBINS],    \
                                zlos[0:NUMLOS*NBINS],    \
                                NUMSEG[0:NUMLOS],        \
                                NUMSEG_PREV[0:NUMLOS],   \
                                xseg[0:NUMSEG_TOT],      \
                                yseg[0:NUMSEG_TOT],      \
                                zseg[0:NUMSEG_TOT],      \
                                rhoker[0:NUMLOS*NBINS],  \
                                velker[0:NUMLOS*NBINS],  \
                                tempker[0:NUMLOS*NBINS]) \
                         copyin(P[0:NUMPARTICLES],       \
                                NUMPARTICLES,            \
                                npart[0:ntasks],         \
                                rnpart[0:ntasks])
    #pragma acc loop independent
    #endif
    /* Loop over sightlines */
    for (int losnum = 0; losnum < NUMLOS; losnum++) {

      /* Loop over particles (each task takes a subset of particles) */
      for (int partnum = rnpart[rank]; 
               partnum < rnpart[rank] + npart[rank]; 
               partnum++) {
        
        /* Particle coordinates */
        double xx, yy, zz, hh;
        xx   = P[partnum].Pos[0]; /* ckpc/h */
        yy   = P[partnum].Pos[1];
        zz   = P[partnum].Pos[2];
        hh   = P[partnum].h * 0.5; /* factor of 0.5 from kernel definition */
        
        /* Check if particle contributes to sightline  */
        double sx, sy, sz, dist, distx, disty, distz, vecp;
        dist = 2.0 * hh + 10.0; /* value larger than 2*hh */

        /* Loop over all segments of sightline */
        for (int k = 0; k < NUMSEG[losnum]; k++) {

          /* a coordinate position along this line segment */
          sx = xseg[NUMSEG_PREV[losnum]+k];
          sy = yseg[NUMSEG_PREV[losnum]+k];
          sz = zseg[NUMSEG_PREV[losnum]+k];

          /* determine the nearest distance from the particle to this line segment */
          vecp  = (sx - xx) * dir[0] + (sy - yy) * dir[1] + (sz - zz) * dir[2];
          distx = fabs(sx - xx - vecp*dir[0]);
          disty = fabs(sy - yy - vecp*dir[1]);
          distz = fabs(sz - zz - vecp*dir[2]);
            
          /* periodicity */
          if (distx > BOXLEN/2) distx = BOXLEN - distx;
          if (disty > BOXLEN/2) disty = BOXLEN - disty;
          if (distz > BOXLEN/2) distz = BOXLEN - distz;

          if (distx*distx + disty*disty + distz*distz < dist*dist)
            dist = sqrt(distx*distx + disty*disty + distz*distz);
        }

        /* dist now minimum distance between the particle and the sightline */
        if (dist <= 2.0*hh) {
          /* Particle properties */
          double hinv, hinv3, vel, mass, temp;
          vel  = P[partnum].Vel[0]*dir[0] + P[partnum].Vel[1]*dir[1]
                 + P[partnum].Vel[2]*dir[2];
          mass = P[partnum].Mass;
          temp = P[partnum].U;

          hinv = 1.0/hh;
          hinv3 = hinv*hinv*hinv;

          /* Conversion factors from internal units */
          double rscale, vscale, mscale, escale, dscale, mu;
          rscale = (KPC*ATIME)/SMALLH;            /* length to cm */
          vscale = sqrt(ATIME);                   /* velocity to kms^-1 */
          mscale = (1.0e10*SOLAR_MASS)/SMALLH;    /* mass to g */
          escale = 1.0e10;                        /* (km/s)^2 to (cm/s)^2 */
          dscale = mscale/(rscale*rscale*rscale); /* density to g cm^-3 */

          /* Convert to appropriate units */
          mu    = 1.0/(XH*(0.75 + P[partnum].Ne) + 0.25); /* mean molecular weight */
          temp *= ((GAMMA-1.0)*mu*HMASS*AMU)/BOLTZMANN;   /* K/escale */

          /* Loop over cells within sightline */
          #if defined(_OPENACC)
          #pragma acc loop independent
          #endif
          for (int cell = 0; cell < NBINS; cell++) {
            double kernel, dx, dy, dz, dr2, dr;
            int cindex = losnum*NBINS + cell;

            /* distance from particle to cell */
            dx = fabs(xx - xlos[cindex]);
            dy = fabs(yy - ylos[cindex]);
            dz = fabs(zz - zlos[cindex]);

            /* periodicity */
            if (dx > BOXLEN/2) dx = BOXLEN - dx;
            if (dy > BOXLEN/2) dy = BOXLEN - dy;
            if (dz > BOXLEN/2) dz = BOXLEN - dz;

            dr2 = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr2);

            /* particle contributes to cell if within smoothing length */
            if (dr <= 2.0*hh) {
              double q;
              q = dr*hinv; /* q = r/h */

              kernel = (q <= 1.0) ? 0.318309886 - 0.238732414 * q * q * (2.0 - q)
                                  : 0.079577471 * (2.0 - q) * (2.0 - q) * (2.0 - q);
              kernel *= hinv3;
              kernel  = dscale * kernel * mass;
              rhoker[cindex]  += kernel * XH;                 /* g cm^-3      */
              velker[cindex]  += vscale * kernel * vel  * XH; /* g cm^-3 km/s */
              tempker[cindex] += escale * kernel * temp * XH; /* g cm^-3 K    */
            } /* particle contribution to given cell */
          } /* loop over cells */
        } /* contributing particles */
      } /* loop over particles */
    } /* loop over sightlines */
   
    t_extract += MPI_Wtime() - t_e0;

    /* deallocate particle data for this file on host */
    delete[] P;

    cout << "[Task " << rank << "] finished particles " 
         << rnpart[rank] << " to " << rnpart[rank] + npart[rank] 
         << "!" << endl;

  } /* Loop over files */
  
  /* Remove data from GPU device(s) */
  #if defined(_OPENACC)
  #pragma acc exit data delete(NBINS, drbin, NUMLOS,    \
                               BOXLEN, SMALLH, ATIME,   \
                               dir[0:3],                \
                               xlos[0:NUMLOS*NBINS],    \
                               ylos[0:NUMLOS*NBINS],    \
                               zlos[0:NUMLOS*NBINS],    \
                               NUMSEG[0:NUMLOS],        \
                               NUMSEG_PREV[0:NUMLOS],   \
                               xseg[0:NUMSEG_TOT],      \
                               yseg[0:NUMSEG_TOT],      \
                               zseg[0:NUMSEG_TOT])      \
                       copyout(rhoker[0:NUMLOS*NBINS],  \
                               velker[0:NUMLOS*NBINS],  \
                               tempker[0:NUMLOS*NBINS]) 
  #endif

  /* Reduce onto rank 0 and save */
  if (rank == 0) {
    frho = new double[NUMLOS*NBINS];
    fvel = new double[NUMLOS*NBINS];
    ftemp = new double[NUMLOS*NBINS];
  }
  MPI_Reduce(rhoker, frho, NUMLOS*NBINS, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(velker, fvel, NUMLOS*NBINS, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(tempker, ftemp, NUMLOS*NBINS, MPI_DOUBLE, MPI_SUM, 0, comm);

  /* Write to file */
  int werr = 0;
  if (rank == 0) {
    try {
      write_to_file();
    } catch (invalid_argument err) {
        cout << err.what();
        werr = 1;
      }
  } 
  MPI_Bcast(&werr, 1, MPI_INT, 0, comm);
  if (werr == 1) {
    MPI_Finalize();
    return 1;
  }

  deallocate();

  /* Print timing statistics */
  double t_end = MPI_Wtime();
  if (rank == 0) {
    cout << "Wallclock statistics: [seconds]" << endl;
    cout << "Total time = " << t_end - t_start << endl;
    cout << "Snapshot loading = " << t_file << endl;
    cout << "Extraction = " << t_extract << endl;
    cout << "Initialization = " << t_init << "   :-)" << endl << endl;
  }

  /* Finish */
  MPI_Finalize();
  return 0;
}
