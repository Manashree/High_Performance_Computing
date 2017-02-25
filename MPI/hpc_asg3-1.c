#include "mpi.h"
#include <stdio.h>
#include <stddef.h>

typedef struct {
  int max_iter;
  double t0;
  double tf;
  double xmax[12];
  double xmin;
} Pars;

int main(int argc, char *argv[])
{
    int myid, numprocs, left, right;
    Pars buffer, buffer2;
    MPI_Request request;
    MPI_Status status;
 
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    right = (myid + 1) % numprocs;
    left = myid - 1;
    if (left < 0)
        left = numprocs - 1;

    //initialize the send buffer
    buffer.max_iter = myid;
    buffer.t0 = 3.14*myid;
    buffer.tf = 1.67*myid;
    buffer.xmin = 2.55*myid;
    for(int i=0;i<12;i++) {
      buffer.xmax[i] = 2.7*myid;
    }

    int nitems = 5;
    MPI_Datatype types[nitems];
    MPI_Datatype mpi_par;  	
    MPI_Aint offsets[nitems];
    int blocklengths[nitems];

    types[0] = MPI_INT;    
	offsets[0] = offsetof(Pars,max_iter); 
	blocklengths[0] = 1;
    types[1] = MPI_DOUBLE; 
	offsets[1] = offsetof(Pars,t0);       
	blocklengths[1] = 1;
    types[2] = MPI_DOUBLE; 
	offsets[2] = offsetof(Pars,tf);       
	blocklengths[2] = 1;
    types[3] = MPI_DOUBLE; 
	offsets[3] = offsetof(Pars,xmax);     
	blocklengths[3] = 12;
    types[4] = MPI_DOUBLE; 
	offsets[4] = offsetof(Pars,xmin);     
	blocklengths[4] = 1;

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_par);
    MPI_Type_commit(&mpi_par);

    MPI_Sendrecv(&buffer, 1, mpi_par, left, 123,&buffer2, 1, mpi_par, right, 123, MPI_COMM_WORLD, &status);

    printf(" Process %d received %d\n",myid,buffer2.max_iter);

    MPI_Finalize();
    return 0;
}
/*
 *  Process 0 received 1
 Process 9 received 10
 Process 10 received 11
 Process 8 received 9
 Process 15 received 0
 Process 14 received 15
 Process 1 received 2
 Process 5 received 6
 Process 12 received 13
 Process 7 received 8
 Process 13 received 14
 Process 11 received 12
 Process 4 received 5
 Process 6 received 7
 Process 3 received 4
 Process 2 received 3

 * */

