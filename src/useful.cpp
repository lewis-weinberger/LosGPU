/* Copyright 2018 Lewis Weinberger.                                           */
/* GPU accelerated sightline extraction for Sherwood (GADGET) snapshots.      */

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include "./global_vars.h"
#include "./prototypes.h"

using std::cout;
using std::endl;
using std::string;
using std::stoi;
using std::stod;
using std::vector;
using std::getline;

/* Global variables */
int NUMLOS, NBINS;
double LENLOS, drbin;
string output_path;
double *xlos_mid, *ylos_mid, *zlos_mid;
double dir[3];
double *xlos, *ylos, *zlos;
double *rhoker, *velker, *tempker;
int *NUMSEG, *NUMSEG_PREV;
int NUMSEG_TOT;
double *xseg, *yseg, *zseg;

void read_params(const string& param_path, int rank) {
  /* Read in parameter file */
  std::ifstream pfile(param_path);
  string line;

  if (pfile.is_open()) {
    /* First line contains a comment */
    getline(pfile, line);

    /* Second line contains general sightline info */
    getline(pfile, line);
    vector<string> sep = split(line);
    NUMLOS      = stoi(sep[0]);
    NBINS       = stoi(sep[1]);
    LENLOS      = stod(sep[2]);
    dir[0]      = stod(sep[3]);
    dir[1]      = stod(sep[4]);
    dir[2]      = stod(sep[5]);
    output_path = sep[6];

    /* Third line contains a comment */
    getline(pfile, line);

    /* Allocate memory for sightline arrays */
    xlos_mid = new double[NUMLOS];
    ylos_mid = new double[NUMLOS];
    zlos_mid = new double[NUMLOS];
    NUMSEG   = new int[NUMLOS];

    int linenum = 0;
    while (getline(pfile, line)) {
      /* Check to see if incorrect number of sightlines given */
      if (linenum == NUMLOS) {
        if (rank == 0) {
          cout << "NUMLOS = " << NUMLOS << endl;
          cout << "Given  = " << linenum << endl;
        }
        throw std::invalid_argument("Number of sightlines incorrect!\n");
      }

      /* The remaining lines contain sightline midpoint coordinates */
      vector<string> sep = split(line);
      xlos_mid[linenum] = stod(sep[0]);
      ylos_mid[linenum] = stod(sep[1]);
      zlos_mid[linenum] = stod(sep[2]);
      linenum++;
    }
    pfile.close();
  } else {
      throw std::invalid_argument("Unable to open parameter file!\n");
    }

  /* Normalize direction vector */
  double mag_dir = 0.0;
  for (int i = 0; i < 3; i++) {
    mag_dir += dir[i]*dir[i];
  }
  mag_dir = std::sqrt(mag_dir);

  if (mag_dir <= 1e-16)
    throw std::invalid_argument("Please give a real (non-zero) direction!\n");

  for (int i = 0; i < 3; i++) {
    dir[i] /= mag_dir;
  }

//  double max_len2 = 3*BOXLEN*BOXLEN*1e-6;
//  if (LENLOS*LENLOS > max_len2)
//    throw std::invalid_argument("Please give a non-overlapping length!\n");

  /* Bin length along direction */
  drbin = LENLOS/NBINS*1e3; /* units of ckpc/h */

  /* Sanity prints */
  if (rank == 0) {
    cout << "Sightline info:" << endl;
    cout << "NUMLOS    = " << NUMLOS << endl;
    cout << "NBINS     = " << NBINS << endl;
    cout << "LENLOS    = " << LENLOS << " cMpc/h" << endl;
    cout << "drbin     = " << drbin << " ckpc/h" << endl;
    cout << "Direction = (" << dir[0] << ", " << dir[1] << ", " << dir[2]
         << ")" << endl;
    cout << "Output file path: " << output_path << endl << endl;
    cout << "Example midpoint coordinates: " << endl;
    cout << "xlos_mid[0] = " << xlos_mid[0] << " cMpc/h" << endl;
    cout << "ylos_mid[0] = " << ylos_mid[0] << " cMpc/h" << endl;
    cout << "zlos_mid[0] = " << zlos_mid[0] << " cMpc/h" << endl << endl;
  }

  /* Generate sightline position data */
  int nb2 = NBINS/2;
  double xstart, ystart, zstart;

  /* Allocate memory for sightline data */
  xlos    = new double[NUMLOS*NBINS];
  ylos    = new double[NUMLOS*NBINS];
  zlos    = new double[NUMLOS*NBINS];
  rhoker  = new double[NUMLOS*NBINS];
  velker  = new double[NUMLOS*NBINS];
  tempker = new double[NUMLOS*NBINS];

  for (int i = 0; i < NUMLOS; i++) {
    /* Starting coordinates in ckpc/h - note check periodicity */
    xstart = wrap<double>(xlos_mid[i]*1e3 - nb2 * drbin * dir[0], BOXLEN);
    ystart = wrap<double>(ylos_mid[i]*1e3 - nb2 * drbin * dir[1], BOXLEN);
    zstart = wrap<double>(zlos_mid[i]*1e3 - nb2 * drbin * dir[2], BOXLEN);

    /* Segment data */
    NUMSEG[i] = 1;

    /* Loop over all bins in the given direction */
    for (int j = 0; j < NBINS; j++) {
      /* wrap around periodic boundaries */
      xlos[i*NBINS + j] = wrap<double>(xstart + j * drbin * dir[0], BOXLEN);
      ylos[i*NBINS + j] = wrap<double>(ystart + j * drbin * dir[1], BOXLEN);
      zlos[i*NBINS + j] = wrap<double>(zstart + j * drbin * dir[2], BOXLEN);

      /* if we cross a boundary, we must have another line segment */
      if (j > 0) {
        if (abs(xlos[i*NBINS + j] - xlos[i*NBINS + j - 1]) > drbin*dir[0]
         || abs(ylos[i*NBINS + j] - ylos[i*NBINS + j - 1]) > drbin*dir[1]
         || abs(zlos[i*NBINS + j] - zlos[i*NBINS + j - 1]) > drbin*dir[2])
          NUMSEG[i] += 1;
      }

      /* Zero-initialise */
      rhoker[i*NBINS + j] = 0.0;
      velker[i*NBINS + j] = 0.0;
      tempker[i*NBINS + j] = 0.0;
    }
  }

  /* Determine segment points */
  NUMSEG_PREV = new int[NUMLOS];
  NUMSEG_PREV[0] = 0;
  NUMSEG_TOT = NUMSEG[0];
  for (int i = 1; i < NUMLOS; i++) {
    /* Cumulative number of line segments */
    NUMSEG_PREV[i] = NUMSEG[i-1] + NUMSEG_PREV[i-1];
    /* Total number of points stored */
    NUMSEG_TOT += NUMSEG[i];
  }
  xseg = new double[NUMSEG_TOT];
  yseg = new double[NUMSEG_TOT];
  zseg = new double[NUMSEG_TOT];

  int n;
  double diffx, diffy, diffz;
  for (int i = 0; i < NUMLOS; i++) {
    /* add  point for first segment */
    xseg[NUMSEG_PREV[i]] = xlos[i*NBINS];
    yseg[NUMSEG_PREV[i]] = ylos[i*NBINS];
    zseg[NUMSEG_PREV[i]] = zlos[i*NBINS];
    
    n = 1;
    for (int j = 1; j < NBINS; j++) {
      diffx = abs(xlos[i*NBINS + j] - xlos[i*NBINS + j - 1]);
      diffy = abs(ylos[i*NBINS + j] - ylos[i*NBINS + j - 1]);
      diffz = abs(zlos[i*NBINS + j] - zlos[i*NBINS + j - 1]);

      /* If we cross a boundary the difference will be a box length,
       * so we can save the new point as belonging to the next segment */
      if (diffx > drbin*dir[0] || diffy > drbin*dir[1] || diffz > drbin*dir[2]) {
        xseg[NUMSEG_PREV[i]+n] = xlos[i*NBINS + j];
        yseg[NUMSEG_PREV[i]+n] = ylos[i*NBINS + j];
        zseg[NUMSEG_PREV[i]+n] = zlos[i*NBINS + j];
        n++;
      }
    }

  }

  /* Sanity prints */
  if (rank == 0) {
    cout << "Example generated coordinates:" << endl;
    cout << "xlos[0][0]         = " << xlos[0] << endl;
    cout << "xlos[0][NBINS-1]   = " << xlos[NBINS-1] << endl;
    cout << "ylos[0][0]         = " << ylos[0] << endl;
    cout << "ylos[0][NBINS-1]   = " << ylos[NBINS-1] << endl;
    cout << "zlos[0][0]         = " << zlos[0] << endl;
    cout << "zlos[0][NBINS-1]   = " << zlos[NBINS-1] << endl << endl;
  }
}

template<class T> T wrap(T coord, T len) {
  /* Wraps coordinates or grid indices (type-dependent) based on box
   * periodicity */

  if (coord < 0.0) {
    coord += len;
  } else if (coord >= len) {
      coord -= len;
    }

  return coord;
}

vector<string> split(const string& s) {
  /* Split string into vector of words. Based on Koenig & Moo (2000) page 88 */
  vector<string> ret;
  typedef string::size_type string_size;
  string_size i = 0;

  while (i != s.size()) {
    while (i != s.size() && std::isspace(s[i]))
      ++i;

    string_size j = i;
    while (j != s.size() && !std::isspace(s[j]))
      ++j;

    if (i != j) {
      ret.push_back(s.substr(i, j - i));
      i = j;
    }
  }
  return ret;
}
