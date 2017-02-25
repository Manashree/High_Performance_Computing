/* Write an embarrassingly parallel code using MPI which computes the first 10 digits of pi
* (i.e. 3.14159….). Run it on bigred2 and submit the strong scaling performance plot.*/

#include "mpi.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 1E2
#define d 1E-2

double generateRandom(){
    return (double)rand() / (double)(RAND_MAX+1.0);
}

int main(int argc, char** argv){

    int rank, numprocs, parameter, sum, pointsInCircle,j;
    int root = 0;
    double pi, x_coord, y_coord;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // identify the rank
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    srand((int)time(0));
	
    for(j = rank; j < N; j += numprocs) {
		x_coord = generateRandom();
		y_coord = generateRandom();

		/* if dart lands in circle, increment score */
		if((x_coord*x_coord) + (y_coord*y_coord) < 1.0) {
			sum++;
		}
    }

	MPI_Reduce(&sum, &pointsInCircle, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(rank == root) {
		// calculate pi
		pi = 4.0 * d * pointsInCircle;
		printf("Pi is: %0.2f : ", pi);
	}
	
	MPI_Finalize();
	return 0;
}

/**
 *  time aprun -n 64 ./a.out 
Pi is: 3.1424358400
Application 6989551 resources: utime ~25s, stime ~10s, Rss ~9776, inblocks ~18078, outblocks ~31288

real	0m3.467s
user	0m0.292s
sys	0m0.032s

 * > time aprun -n 32 ./a.out 
Pi is: 3.1412992000
Application 6989552 resources: utime ~5s, stime ~2s, Rss ~3944, inblocks ~8695, outblocks ~15644

real	0m2.431s
user	0m0.296s
sys	0m0.024s

 *  time aprun -n 16 ./a.out 
Pi is: 3.1407033600
Application 6989554 resources: utime ~5s, stime ~1s, Rss ~3944, inblocks ~7446, outblocks ~15643

real	0m2.416s
user	0m0.288s
sys	0m0.036s
 *  time aprun -n 8 ./a.out 

Pi is: 3.1411152000
Application 6989555 resources: utime ~5s, stime ~1s, Rss ~3944, inblocks ~6821, outblocks ~15642

real	0m2.763s
user	0m0.308s
sys	0m0.012s

 *  time aprun -n 4 ./a.out 
Pi is: 3.1411227200
Application 6989560 resources: utime ~4s, stime ~0s, Rss ~3944, inblocks ~6509, outblocks ~15642

real	0m3.226s
user	0m0.288s
sys	0m0.052s
manarao@aprun8:~/codes-3> time aprun -n 2 ./a.out 
Pi is: 3.1414593600
Application 6989574 resources: utime ~4s, stime ~0s, Rss ~3944, inblocks ~6352, outblocks ~15642

real	0m4.197s
user	0m0.304s
sys	0m0.016s
manarao@aprun8:~/codes-3> time aprun -n 1 ./a.out 
Pi is: 3.1415392000
Application 6989586 resources: utime ~4s, stime ~0s, Rss ~3944, inblocks ~6274, outblocks ~15642

real	0m5.841s
user	0m0.288s
sys	0m0.052s

 * */
