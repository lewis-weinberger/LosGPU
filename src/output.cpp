/* Copyright 2018 Lewis Weinberger.                                           */
/* GPU accelerated sightline extraction for Sherwood (GADGET) snapshots.      */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ios>
#include <stdexcept>
#include "./global_vars.h"
#include "./prototypes.h"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::vector;

void write_to_file() {
  /* For density normalization */
  double H0, rhoc, critH;
  H0    = 1.0e7/MPC; /* 100 km s^-1 Mpc^-1 in cgs */
  rhoc  = 3.0*(H0*SMALLH)*(H0*SMALLH)/(8.0*PI*GRAVITY); /* g cm^-3 */
  critH = rhoc*OMEGAB*XH/(ATIME*ATIME*ATIME); /* g cm^-3*/

  /* Normalise velocity and temperature */
  for (int index = 0; index < NBINS*NUMLOS; index++) {
    fvel[index]  /= frho[index]; /* km/s */
    ftemp[index] /= frho[index]; /* K    */
    frho[index]  /= critH;
  }

  cout << "Saving to: " << output_path << endl << endl;

  /* Open file in binary mode */
  ofstream outfile(output_path, std::ios::binary);
  if (outfile.is_open()) {

    /* Write data to file */
    outfile.write((char*)&NBINS, sizeof(int));
    outfile.write((char*)&NUMLOS, sizeof(int));
    outfile.write((char*)&LENLOS, sizeof(double));

    outfile.write((char*)NUMSEG, sizeof(int)*NUMLOS);

    outfile.write((char*)xlos_mid, sizeof(double)*NUMLOS);
    outfile.write((char*)ylos_mid, sizeof(double)*NUMLOS);
    outfile.write((char*)zlos_mid, sizeof(double)*NUMLOS);

    outfile.write((char*)xlos, sizeof(double)*NUMLOS*NBINS);
    outfile.write((char*)ylos, sizeof(double)*NUMLOS*NBINS);
    outfile.write((char*)zlos, sizeof(double)*NUMLOS*NBINS);

    outfile.write((char*)frho, sizeof(double)*NUMLOS*NBINS);
    outfile.write((char*)fvel, sizeof(double)*NUMLOS*NBINS);
    outfile.write((char*)ftemp, sizeof(double)*NUMLOS*NBINS);

    /* Close file */
    outfile.close();
  } else {
      throw(std::invalid_argument("Unable to open output file!\n"));
    }

  /* Deallocate memory (only allocated on root task) */
  delete frho;
  delete fvel;
  delete ftemp;
}

void deallocate() {
  /* Deallocates remaining dynamic arrays */
  delete rhoker;
  delete velker;
  delete tempker;
  delete xlos;
  delete ylos;
  delete zlos;
  delete xlos_mid;
  delete ylos_mid;
  delete zlos_mid;
  delete NUMSEG;
  delete NUMSEG_PREV;
  delete xseg;
  delete yseg;
  delete zseg;
}

void print_title() {
  /* Prints nice FIGLET title */
  vector<string> title;
  title.push_back(" ___      _______  _______  _______  _______  __   __ ");
  title.push_back("|   |    |       ||       ||       ||       ||  | |  |");
  title.push_back("|   |    |   _   ||  _____||    ___||    _  ||  | |  |");
  title.push_back("|   |    |  | |  || |_____ |   | __ |   |_| ||  |_|  |");
  title.push_back("|   |___ |  |_|  ||_____  ||   ||  ||    ___||       |");
  title.push_back("|       ||       | _____| ||   |_| ||   |    |       |");
  title.push_back("|_______||_______||_______||_______||___|    |_______|");

  cout << endl;
  for (vector<string>::iterator it = title.begin(); it != title.end(); it++) {
    cout << *it <<endl;
  }
  cout << endl;
}
