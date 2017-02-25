/*
Parallel sorting by regular sampling(distributed quicksort algorithm)
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#define N 27

/*
p: no of proc
n: no of elements in array

1. divide data into 'p' parts and send to 'p' processors
2. each processor sorts its own data usin local seq. quicksort
3. processor selects data at index 0, n/(p-1)^2, .... n(p-1)/(p^2) as samples
4. root will gather and sort these samples using seq quicksort.
5. root selects (p-1) pivot elements from the sample list in 4, and broadcast these.
5. All other processes will split their own samples into 'p' parts.
6. Each process i keeps its ith partition and sends the jth partition to process j, for all j ne i
7. Each process merges its P partitions into a single list.
*/

/* Comparison function. Receives two generic (void) pointers to the items under comparison. */
int compare_ints(const void* p, const void* q)
{
    int x = *(const int*)p;
    int y = *(const int*)q;

    /* Avoid return x - y, which can cause undefined behaviour
       because of signed integer overflow. */
    if(x < y)
	return -1; // Return -1 if you want ascending, 1 if you want descending order.
    else if(x > y)
	return 1; // Return 1 if you want ascending, -1 if you want descending order.

    return 0;
}

void print_array(int* arr){
	int l = sizeof(arr)/sizeof(arr[0]);
	for(int i = 0; i<l; i++){
		printf("%d ",arr[i]);
	}
}

int main(int argc, char* argv[])
{
    int rank, size,length=0;
    int myData[N] = { 15, 46, 48, 93, 39, 6, 72, 91, 14, 36, 69, 40, 89, 61, 97, 12, 21, 54, 53, 97, 84, 58, 32, 27, 33, 72, 20 };

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int myDataLengths[size];
    int myDataStarts[size];

    // communication buffer used for determination of pivot values
    int pivotbuffer[size * size],tempbuffer[size-1],segLenSent[size-1],segLenRecv[size-1];
    int pivotbufferSize;

	int recvbuffer[N]; 
	int recvLengths[size];
	int recvStarts[size];
	
	int classStart[size];
	int classLength[size];
	int dataindex=0;
	
    for(int i = 0; i < size; i++) {
		myDataLengths[i] = N / size;
		myDataStarts[i] = i * N / size;
    }
    myDataLengths[size - 1] += (N % size);

    // if root scatter else receive data
    if(rank == 0) {
		MPI_Scatter(&myData, myDataLengths[0], MPI_INT, MPI_IN_PLACE, myDataLengths[rank], MPI_INT,0, MPI_COMM_WORLD);
    } else {
		MPI_Scatter(&myData, myDataLengths[rank], MPI_INT, myData, myDataLengths[rank], MPI_INT,0, MPI_COMM_WORLD);
    }

    // All processors sort their piece of the data using cstdlib::quicksort
    qsort(myData, myDataLengths[rank], sizeof(int), compare_ints);

	// All processors collect regular samples from sorted list; an offset to the myData[] index
	for(int index = 0; index<size; index++){
		pivotbuffer[index]= myData[index*myDataLengths[rank]/size];
    }
	
	// Step 4a - root gets all data & append arrays
	if(rank == 0) {
		MPI_Gather(MPI_IN_PLACE, size, MPI_INT, &pivotbuffer, size, MPI_INT, 0, MPI_COMM_WORLD);
		length = sizeof(pivotbuffer)/sizeof(pivotbuffer[0]);
		printf("Before sorting samples: ");
		for(int i = 0; i<length; i++){
			printf("%d ",pivotbuffer[i]);
		}
		printf("\n");
		
		qsort(pivotbuffer, length, sizeof(int), compare_ints);
		printf("After sorting samples: ");
		for(int i = 0; i<length; i++){
			printf("%d ",pivotbuffer[i]);
		}
		printf("\n");
		// All processors collect regular samples from sorted list; an offset to the myData[] index
		for(int index = 0; index<size; index++){
			tempbuffer[index]= pivotbuffer[(index+1)*size];
		}
	}else{
		MPI_Gather(&pivotbuffer, size, MPI_INT, &pivotbuffer, size, MPI_INT, 0, MPI_COMM_WORLD);
	}
	MPI_Bcast(&tempbuffer, size-1, MPI_INT, 0,MPI_COMM_WORLD);
	length = sizeof(tempbuffer)/sizeof(tempbuffer[0]);
	printf("After broadcast: ");
	for(int i = 0; i<length; i++){
		printf("%d ",tempbuffer[i]);
	}
	printf("\n");
	//start partitioning list

	
	for(int classindex=0; classindex<size-1; classindex++){
		classStart[classindex] = dataindex;
		classLength[classindex]=0;

		// as long as dataindex refers to data in the current class
		while((dataindex< myDataLengths[rank]) && (myData[dataindex]<=tempbuffer[classindex])){
			classLength[classindex]++;
			dataindex++;
		}		
	}
	
	// set Start and Length for last class
	classStart[size-1] = dataindex;
	classLength[size-1] = myDataLengths[rank] - dataindex;
	
	// PHASE V:  All ith classes are gathered by processor i 
//http://stackoverflow.com/questions/15049190/difference-between-mpi-allgather-and-mpi-alltoall-functions
	
	for(int iprocessor=0; iprocessor<size; iprocessor++){
		MPI_Gather(&classLength[iprocessor], 1, MPI_INT, recvLengths, size, MPI_INT, 0, MPI_COMM_WORLD);
		printf("After broadcast of LENGTH for processor %d- length gathered is: \n",rank);
		length = sizeof(recvLengths)/sizeof(recvLengths[0]);
		for(int i = 0; i<length; i++){
			printf("%d ",recvLengths[i]);
		}
		printf("\n");
		//calculate displacements
		if (rank == iprocessor){
			recvStarts[0]=0;
			for(int i=1;i<size; i++){
				recvStarts[i] = recvStarts[i-1]+recvLengths[i-1];
			}
		}		
		MPI_Alltoallv(&myData,classLength,classStart,MPI_INT,&recvbuffer,recvLengths,recvStarts,MPI_INT,MPI_COMM_WORLD);
	}
	
	printf("After broadcast of partitions for processor %d- ",rank);
	length = sizeof(recvbuffer)/sizeof(recvbuffer[0]);
	for(int i = 0; i<length; i++){
		printf("%d ",recvbuffer[i]);
	}
	
    MPI_Finalize();
    return 0;
}
