#include "mpi.h"
#include <stdio.h>
#include <stddef.h>

int proc_map(int i, int no_procs){
    return (i % (no_procs - 1)) + 1;
}

int main(int argc, char** argv){
    int nproc, rank; 
    const int nrows = 4;
    const int ncols = 4;
    MPI_Status Stat;
 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
    if (rank == 0) {
        int mat[nrows][ncols], vec[ncols], res[ncols];
	for (int i = 0; i < nrows; i++) {
	    for (int j = 0; j < ncols; j++) {
		mat[i][j] = i*j;
	    }	
	    vec[i] = 2;
	}

        for (int i = 0; i < nrows; i++) {
	    int proc_idx = proc_map(i, nproc);
            MPI_Send(vec, ncols, MPI_INT, proc_idx, 99, MPI_COMM_WORLD);
            MPI_Send(mat[i], ncols, MPI_INT, proc_idx, 100*(i+1), MPI_COMM_WORLD);
        }
        for (int i = 0; i < nrows; i++) {
            int proc_idx = proc_map(i, nproc);
            MPI_Recv(&res[i], 1, MPI_INT, proc_idx, i, MPI_COMM_WORLD, &Stat);
            printf("P%d : res[%d]t= %dn", rank, i, res[i]);
        }
    } else {
        int vec_recv[ncols], row_recv[nrows];
        for (int i = 0; i < nrows; i++) {
            int proc_idx = proc_map(i, nproc);
            if (rank == proc_idx) {
		MPI_Recv(vec_recv, ncols, MPI_INTEGER, 0, 99, MPI_COMM_WORLD, &Stat);
                MPI_Recv(row_recv, ncols, MPI_INT, 0, 100*(i+1), MPI_COMM_WORLD, &Stat);
                int sum = 0;
                for (int j = 0; j < ncols; j++){
                    sum += (row_recv[j] * vec_recv[j]);
                }
                MPI_Send(&sum, 1, MPI_INT, 0, i, MPI_COMM_WORLD);
            }
        }
    }
 
    MPI_Finalize();
    return 0;
}