#ifndef _GLOBAL_VARS_H_
#define _GLOBAL_VARS_H_

/* Copyright 2018 Lewis Weinberger.                                           */
/* GPU accelerated sightline extraction for Sherwood (GADGET) snapshots.      */

#include <vector>
#include <string>

/* Sightline parameters */
extern int NUMLOS, NBINS;
extern double LENLOS, drbin;
extern std::string output_path;
extern double *xlos_mid, *ylos_mid, *zlos_mid;
extern double dir[3];
extern double *xlos, *ylos, *zlos;
extern double *rhoker, *velker, *tempker;
extern int *NUMSEG, *NUMSEG_PREV;
extern int NUMSEG_TOT;
extern double *xseg, *yseg, *zseg;

/* Particle data */
extern int NUMFILES, NUMPARTICLES;
extern double BOXLEN, REDSHIFT, OMEGA0, OMEGALAMBDA, SMALLH, ATIME;
struct particle_data {
  /* Position and velocity vectors */
  float Pos[3], Vel[3];
  /* Mass, internal energy, smoothing length, something, something */
  float Mass, U, h, NH0, Ne;
};
extern particle_data *P;

/* output */
extern double *frho, *fvel, *ftemp;

/* Numbers */
#define  PI          3.14159265358979323846
#define  GAMMA       (5.0/3.0)
#define  GRAVITY     6.67384e-8
#define  BOLTZMANN   1.3806488e-16
#define  C           2.99792458e10
#define  AMU         1.66053886e-24 /* 1 a.m.u */
#define  MPC         3.08568025e24
#define  KPC         3.08568025e21
#define  SOLAR_MASS  1.989e33
#define  HMASS       1.00794  /* Hydrogen mass in a.m.u. */
#define  XH          0.76  /* hydrogen fraction by mass */
#define  OMEGAB      0.0482

#endif // _GLOBAL_VARS_H_
