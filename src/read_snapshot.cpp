/* Copyright 2018 Lewis Weinberger.                                           */
/* GPU accelerated sightline extraction for Sherwood (GADGET) snapshots.      */

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <H5Cpp.h>
#include "./global_vars.h"
#include "./prototypes.h"

using std::cout;
using std::endl;
using std::string;
using H5::H5File;
using H5::Group;
using H5::Attribute;
using H5::DataType;
using H5::DataSet;
using H5::DataSpace;

/* Global variables */
double BOXLEN, REDSHIFT, OMEGA0, OMEGALAMBDA, SMALLH, ATIME;
int NUMFILES, NUMPARTICLES;
particle_data *P;

void read_snapshot_header(const string& snap_path, int rank) {
  /* Opens the header of the first snapshot file */

  string filepath = snap_path + ".0.hdf5";
  if (rank == 0)
    cout << "Opening header from file: " << filepath << endl;

  /* Open file */
  H5File *file = new H5File(filepath.c_str(), H5F_ACC_RDONLY);

  /* Open header group */
  Group *header = new Group(file->openGroup("/Header"));

  /* "NumFilesPerSnapshot" */
  Attribute *attr = new Attribute(header->openAttribute("NumFilesPerSnapshot"));
  DataType *dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &NUMFILES);
  delete attr;
  delete dtype;

  /* "BoxSize" */
  attr  = new Attribute(header->openAttribute("BoxSize"));
  dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &BOXLEN);
  delete attr;
  delete dtype;

  /* "Redshift" */
  attr  = new Attribute(header->openAttribute("Redshift"));
  dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &REDSHIFT);
  delete attr;
  delete dtype;

  /* "Omega0" */
  attr  = new Attribute(header->openAttribute("Omega0"));
  dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &OMEGA0);
  delete attr;
  delete dtype;

  /* "OmegaLambda" */
  attr  = new Attribute(header->openAttribute("OmegaLambda"));
  dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &OMEGALAMBDA);
  delete attr;
  delete dtype;

  /* "HubbleParam" */
  attr  = new Attribute(header->openAttribute("HubbleParam"));
  dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &SMALLH);
  delete attr;
  delete dtype;

  /* "HubbleParam" */
  attr  = new Attribute(header->openAttribute("Time"));
  dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &ATIME);
  delete attr;
  delete dtype;

  /* Close header group */
  delete header;

  /* Close file */
  delete file;

  /* Sanity prints */
  if (rank == 0) {
    cout << "Header info:" << endl;
    cout << "NUMFILES    = " << NUMFILES << endl;
    cout << "BOXLEN      = " << BOXLEN << " ckpc/h" << endl;
    cout << "REDSHIFT    = " << REDSHIFT << endl;
    cout << "OMEGA0      = " << OMEGA0 << endl;
    cout << "OMEGALAMBDA = " << OMEGALAMBDA << endl;
    cout << "SMALLH      = " << SMALLH << endl;
    cout << "ATIME       = " << ATIME << endl;
  }
}

void read_snapshot(const string& snap_path, int i, int rank) {
  string id = std::to_string(i);
  string filepath = snap_path + "." + id + ".hdf5";
  if (rank == 0)
    cout << endl << "Opening snapshot file: " << filepath << endl;

  /* Open file */
  H5File *file = new H5File(filepath.c_str(), H5F_ACC_RDONLY);

  /* Open header group */
  Group *header = new Group(file->openGroup("/Header"));

  /* "NumPart_ThisFile" */
  int NumPart[6];
  Attribute *attr = new Attribute(header->openAttribute("NumPart_ThisFile"));
  DataType *dtype = new DataType(attr->getDataType());
  attr->read(*dtype, &NumPart);
  NUMPARTICLES = NumPart[0];
  if (rank == 0)
    cout << "NUMPARTICLES = " << NUMPARTICLES << endl;
  delete attr;
  delete dtype;
  delete header;

  /* Allocate data structure for this file */
  P = new particle_data[NUMPARTICLES];

  /* Now open dataset and read particle data */
  Group *gas = new Group(file->openGroup("/PartType0"));

  /* Buffer to read into for coordinates and velocities*/
  float *buffer;
  buffer = new float[3*NUMPARTICLES];

  /* Array size for HDF5 memory space */
  hsize_t dim[1];

  /* "Coordinates" */
  DataSet *data = new DataSet(gas->openDataSet("Coordinates"));
  DataSpace *filespace = new DataSpace(data->getSpace());
  dim[0] = 3*NUMPARTICLES;
  DataSpace *memspace = new DataSpace(1, dim);
  dtype = new DataType(data->getDataType());
  data->read(buffer, *dtype, *memspace, *filespace);
  delete data;
  delete dtype;
  delete filespace;
  delete memspace;
  for (int n = 0; n < NUMPARTICLES; n++) {
    P[n].Pos[0] = buffer[3*n];
    P[n].Pos[1] = buffer[3*n + 1];
    P[n].Pos[2] = buffer[3*n + 2];
  }

  /* "Velocities" */
  data = new DataSet(gas->openDataSet("Velocities"));
  filespace = new DataSpace(data->getSpace());
  dim[0] = 3*NUMPARTICLES;
  memspace = new DataSpace(1, dim);
  dtype = new DataType(data->getDataType());
  data->read(buffer, *dtype, *memspace, *filespace);
  delete data;
  delete dtype;
  delete filespace;
  delete memspace;
  for (int n = 0; n < NUMPARTICLES; n++) {
    P[n].Vel[0] = buffer[3*n];
    P[n].Vel[1] = buffer[3*n + 1];
    P[n].Vel[2] = buffer[3*n + 2];
  }

  /* Adjust buffer for other variables */
  delete buffer;
  buffer = new float[NUMPARTICLES];

  /* "Masses" */
  data = new DataSet(gas->openDataSet("Masses"));
  filespace = new DataSpace(data->getSpace());
  dim[0] = NUMPARTICLES;
  memspace = new DataSpace(1, dim);
  dtype = new DataType(data->getDataType());
  data->read(buffer, *dtype, *memspace, *filespace);
  delete data;
  delete dtype;
  delete filespace;
  delete memspace;
  for (int n = 0; n < NUMPARTICLES; n++) {
    P[n].Mass = buffer[n];
  }

  /* "InternalEnergy" */
  data = new DataSet(gas->openDataSet("InternalEnergy"));
  filespace = new DataSpace(data->getSpace());
  dim[0] = NUMPARTICLES;
  memspace = new DataSpace(1, dim);
  dtype = new DataType(data->getDataType());
  data->read(buffer, *dtype, *memspace, *filespace);
  delete data;
  delete dtype;
  delete filespace;
  delete memspace;
  for (int n = 0; n < NUMPARTICLES; n++) {
    P[n].U = buffer[n];
  }

  /* "SmoothingLength" */
  data = new DataSet(gas->openDataSet("SmoothingLength"));
  filespace = new DataSpace(data->getSpace());
  dim[0] = NUMPARTICLES;
  memspace = new DataSpace(1, dim);
  dtype = new DataType(data->getDataType());
  data->read(buffer, *dtype, *memspace, *filespace);
  delete data;
  delete dtype;
  delete filespace;
  delete memspace;
  for (int n = 0; n < NUMPARTICLES; n++) {
    P[n].h = buffer[n];
  }

  /* "ElectronAbundance" */
  data = new DataSet(gas->openDataSet("ElectronAbundance"));
  filespace = new DataSpace(data->getSpace());
  dim[0] = NUMPARTICLES;
  memspace = new DataSpace(1, dim);
  dtype = new DataType(data->getDataType());
  data->read(buffer, *dtype, *memspace, *filespace);
  delete data;
  delete dtype;
  delete filespace;
  delete memspace;
  for (int n = 0; n < NUMPARTICLES; n++) {
    P[n].Ne = buffer[n];
  }

  /* "NeutralHydrogenAbundance" */
  data = new DataSet(gas->openDataSet("NeutralHydrogenAbundance"));
  filespace = new DataSpace(data->getSpace());
  dim[0] = NUMPARTICLES;
  memspace = new DataSpace(1, dim);
  dtype = new DataType(data->getDataType());
  data->read(buffer, *dtype, *memspace, *filespace);
  delete data;
  delete dtype;
  delete filespace;
  delete memspace;
  for (int n = 0; n < NUMPARTICLES; n++) {
    P[n].NH0 = buffer[n];
  }

  /* deallocate buffer */
  delete buffer;

  /* Close group */
  delete gas;

  /* Close file */
  delete file;

  /* Sanity prints */
  if (rank == 0) {
    cout << "Example particle data:" << endl;
    cout << "P[0].Pos[0] = " << P[0].Pos[0] << " ckpc/h" << endl;
    cout << "P[0].Pos[1] = " << P[0].Pos[1] << " ckpc/h" << endl;
    cout << "P[0].Pos[2] = " << P[0].Pos[2] << " ckpc/h" << endl;
    cout << "P[0].Vel[1] = " << P[0].Vel[0] << " km/s" << endl;
    cout << "P[0].Vel[2] = " << P[0].Vel[1] << " km/s" << endl;
    cout << "P[0].Vel[3] = " << P[0].Vel[2] << " km/s" << endl;
    cout << "P[0].Mass   = " << P[0].Mass << " 1e10 Msun/h" << endl;
    cout << "P[0].U      = " << P[0].U << " (km/s)^2" << endl;
    cout << "P[0].h      = " << P[0].h << " ckpc/h" << endl;
    cout << "P[0].NH0    = " << P[0].NH0 << endl;
    cout << "P[0].Ne     = " << P[0].Ne << endl << endl;
  }
}
