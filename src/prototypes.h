#ifndef _PROTOTYPES_H_
#define _PROTOTYPES_H_

/* Copyright 2018 Lewis Weinberger.                                           */
/* GPU accelerated sightline extraction for Sherwood (GADGET) snapshots.      */

#include "mpi.h"
#include <string>
#include <vector>

/* useful.cpp */
void read_params(const std::string&, int);
std::vector<std::string> split(const std::string&);
template<class T> T wrap(T, T);

/* read_snapshot.cpp */
void read_snapshot_header(const std::string&, int);
void read_snapshot(const std::string&, int, int);

/* send_data.cpp */
void send_header(MPI_Comm);
void send_params(MPI_Comm, int rank);
void send_data(MPI_Comm, int rank);

/* output.cpp */
void write_to_file(void);
void deallocate();
void print_title();

#endif // _PROTOTYPES_H_
