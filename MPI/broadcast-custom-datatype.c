#include <stdio.h>
#include <stddef.h>
#include "mpi.h"

typedef struct {
  int max_iter;
  double t0;
  double tf;
  double xmin;
} Pars;

int main(int argc,char **argv) {

  MPI_Init(&argc,&argv);
  int rank;
  int root = 0;  // define the root process
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); // identify the rank

  Pars pars;
  if ( rank == root ) {
    pars.max_iter = 10;
    pars.t0 = 0.0;
    pars.tf = 1.0;
    pars.xmin = -5.0;
  } 

  int nitems = 4;
  MPI_Datatype types[nitems];
  MPI_Datatype mpi_par;  // give my new type a name
  MPI_Aint offsets[nitems];  // an array for storing the element offsets
  int blocklengths[nitems];

  types[0] = MPI_INT; offsets[0] = offsetof(Pars,max_iter);blocklengths[0] = 1;
  types[1] = MPI_DOUBLE; offsets[1] = offsetof(Pars,t0);blocklengths[1] = 1;
  types[2] = MPI_DOUBLE; offsets[2] = offsetof(Pars,tf);blocklengths[2] = 1;
  types[3] = MPI_DOUBLE; offsets[3] = offsetof(Pars,xmin);blocklengths[3] = 1;

  MPI_Type_create_struct(nitems,blocklengths,offsets,types,&mpi_par);
  MPI_Type_commit(&mpi_par);

  MPI_Bcast(&pars,1,mpi_par,root,MPI_COMM_WORLD);

  printf("Hello from rank %d; my max_iter value is %d\n",rank,pars.max_iter);

  MPI_Finalize();

  return 0;
}


