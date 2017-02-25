/*
Write an MPI code to integrate cos(x) ∗ sin (x/2) on the interval between 0 and p to an
accuracy of at least 1e-6. Plot the strong scaling of the code from 1 to 64 processes.
*/
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846
#define N 1E6
/*Midpoint Rule
 * delta_x = (b-a)/n
 * x_i is the midpoint
 * h is height
 * */
float f_x(float x);
float integrate(float l, float h);
int main(int argc,char** argv){
	float a,b,delta_x,score=0.0,sbuffer=0.0,rbuffer=0.0;
	int rank,size,i;

	MPI_Status status;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	a=0.0;
	b=PI*0.5;
	delta_x =(b-a)/(size-1);
	if(rank==0){
		for(i=1;i<size;i++){
			MPI_Recv(&rbuffer, 1, MPI_FLOAT, i, 999, MPI_COMM_WORLD, &status);
			score += rbuffer;
		}
		printf("Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: %0.10f \n",score);
	}else{
		sbuffer = integrate((delta_x * (rank-1)), delta_x/N);
		MPI_Send(&sbuffer, 1, MPI_FLOAT, 0, 999, MPI_COMM_WORLD);
	}
	MPI_Finalize();
}
float f_x(float x){
	//return(cos(x));
	return(cos(x)*sin(x/2.0));
}
float integrate(float l, float h) {
	float x=0.0;
	for(int i=0;i<N;i++){
		x += f_x(l+(i+0.5)*h)*h;
	}
	return(x);
}
/*
 *   time aprun -n 64 ./a.out 
Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: 0.2761529386 
Application 6990416 resources: utime ~16s, stime ~10s, Rss ~9804, inblocks ~19430, outblocks ~34786

real	0m3.151s
user	0m0.324s
sys	0m0.052s

 *  time aprun -n 32 ./a.out 
Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: 0.2761054933 
Application 6990419 resources: utime ~2s, stime ~2s, Rss ~4888, inblocks ~9381, outblocks ~17393

real	0m2.556s
user	0m0.312s
sys	0m0.064s

 *  time aprun -n 16 ./a.out 
Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: 0.2761433125 
Application 6990420 resources: utime ~1s, stime ~1s, Rss ~4592, inblocks ~8132, outblocks ~17392

real	0m2.385s
user	0m0.336s
sys	0m0.048s

 *  time aprun -n 8 ./a.out 
Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: 0.2761515379 
Application 6990421 resources: utime ~0s, stime ~1s, Rss ~3892, inblocks ~7507, outblocks ~17391

real	0m2.228s
user	0m0.328s
sys	0m0.028s

 *  time aprun -n 4 ./a.out 
Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: 0.2760576010 
Application 6990423 resources: utime ~0s, stime ~0s, Rss ~3892, inblocks ~7195, outblocks ~17391

real	0m2.275s
user	0m0.308s
sys	0m0.064s

 *  time aprun -n 2 ./a.out 
Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: 0.2761739492 
Application 6990424 resources: utime ~0s, stime ~0s, Rss ~3892, inblocks ~7038, outblocks ~17391

real	0m2.285s
user	0m0.320s
sys	0m0.052s

 *  time aprun -n 1 ./a.out 
Definite integral of cos(x)*sin(x/2) from 0 to pi/2 is: 0.0000000000 
Application 6990425 resources: utime ~0s, stime ~0s, Rss ~3892, inblocks ~6960, outblocks ~17391

real	0m2.177s
user	0m0.308s
sys	0m0.068s
 * 
 * */
