/*
Heat distribution - This code is based on a simplified two-dimensional heat equation. The initial temperature is computed to be high in the middle & zero at the boundaries. 
The approach to solving this problem is to divide the area into a fine mesh of points, h[i][j].	Temperatures at an inside point are taken to be the average of temperatures of four neighboring	points.
The program should terminate the computation after a certain precision of computation is reached.
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define COLS 10000
#define ROWS 10000
#define TEMP 50.
#define DEBUG 1
#define EPS 1e-3
#define I_FIX 2
#define J_FIX 2
#define BEGIN 123
#define HEAT 234
#define NONE 0
#define LGHOST 789
#define RGHOST 987
#define WDONE 999
#define FDONE 997
#define MDONE 111 

/* PSEUDOCODE

For MASTER:
1. Compute number of rows for each processor
2. For each worker send rows computed as above, add extra rows to first few processors. The maximum number of rows owned by any process will be at most 1 more than the number of rows owned by a process with the least number of rows.
3. Tell each worker who its neighbours(left & right) are.
4. Send offset, neighbours and row count computed above to each worker.
5. Collect max values from each worker till max value <ETA.
6. Once max value <ETA, stop and print matrix.

For workers:
1. receive offset, rows, left, right
2. allocate memory for the matrices
3. Compute if the heat source-offset exists in the number of rows sent to the worker if it fits in the boundaries of rows for this worker, set the heatFlag which will help other functions to pinpoint the heat source
4. Communicate border rows with neighbors
5. call compute_new_values to update values of grid points
6. Compute maximum value and send to MASTER
7. Wait for proceedFlag from MASTER
8. send this worker's portion of final results back to master if proceedFlag is 0
9. continue with the computation, create copy of matrix if proceedFlag is 1

*/

/*
finds max value in matrix
ignore row0,row=ROWS as they are ghost rows
ignore col=0 and col=COLS as they be border columns hence 0.
*/
double max_abs(int nrows, double** m1, double** m2){
    double max_val = DBL_MIN;
    double temp =0.;
    for (int i = 1; i < nrows-1; i++)
        for (int j = 1; j < COLS-1; j++){
        	temp = fabs(m1[i][j] - m2[i][j]);
            if (temp > max_val) {
                max_val = temp;
            }
        }
    return max_val;
}

void print_matrix(double** matrix,int nrows){
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < COLS; j++)
            printf("%f ",matrix[i][j]);
        printf("\n");
    }
}

/*
copy values from source to destination in matrix
ignore row0,row=ROWS as they are ghost rows
ignore col=0 and col=COLS as they be border columns hence 0.
fix heat source co-ordinates according to heat flag(+1 is done to accomodate for padding for ghost rows)
*/
void copy_matrix(int nrows, double** dest, double** source, int heat, int heat_source_I_FIX) {
    for (int i = 1; i < nrows-1; i++)
        for (int j = 1; j < COLS-1; j++)
            dest[i][j] = source[i][j];
    if(heat == 1)
    	dest[heat_source_I_FIX+1][J_FIX] = TEMP;
}

// Pass number of rows if each processor is calling individually according to its allocation
double** alloc_matrix(int nrows){
    double** matrix;
    matrix = (double**) malloc(nrows * sizeof(double *));
    matrix[0] = (double*) malloc(nrows * COLS * sizeof(double));
    for (int i = 0; i < nrows; i++)
        matrix[i] = matrix[0] + i*COLS;
    return matrix;
}

/*
Compute new values from a_old to a_new in matrix
ignore row=0,row=ROWS as they are ghost rows
ignore row=1 for first processor as it is the border,row=ROWS at rank=size-1 as it is the last border row
ignore col=0 and col=COLS as they be border columns hence 0.
fix heat source co-ordinates according to heat flag(+1 is done to accomodate for padding for ghost rows)
manage left and right for cells adjacent to border rows
*/
void compute_new_values(int nrows, double** old_matrix, double** new_matrix, int heat, int heat_source_I_FIX,int psuedo_rank,int size){
	double left=0., right=0., top=0.,bottom=0.;
    for (int i = 1; i < nrows-1; i++){
        if(psuedo_rank==1 && i==1){
            //printf("Don't touch row 1 for processor %d\n",psuedo_rank);
            continue;
        }
        if(psuedo_rank==size-1 && i==nrows-2){
            //printf("Don't touch row %d for processor %d\n",i, psuedo_rank);
            continue;
        }
        for (int j = 1; j < COLS-1; j++){
        	//printf("i=%d,j=%d for rank %d\n",i,j,psuedo_rank);
        	top = old_matrix[i-1][j]; 
        	bottom = old_matrix[i+1][j];
            left =  (j == 1) ? 0. :old_matrix[i][j-1];
            right = (j == (COLS-2)) ? 0. : old_matrix[i][j+1];
            new_matrix[i][j] = 0.25 * (left + right + top + bottom);
        }
    }
    if(heat == 1)
    	new_matrix[heat_source_I_FIX+1][J_FIX] = TEMP;
}

//initialize in matrix
void init_matrix(int nrows, double** matrix, int heat, int heat_source_I_FIX){
    for (int i = 0; i < nrows; i++)
        for (int j = 0; j < COLS; j++) {
            matrix[i][j] = 0.;
        }
    if(heat == 1){
    	matrix[heat_source_I_FIX+1][J_FIX] = TEMP;
    }
}

double max_in_array(double a[], int num_elements){
   int i;
   double max = DBL_MIN;
   for (i=1; i<num_elements; i++){
	 if (a[i]>max){
	    max=a[i];
	 }
   }
   return(max);
}


int main(int argc, char *argv[]) {
	int rank,size,i=0, avgRows=0, remainder=0;
	int left=0,right=0,rows=0, MASTER=0, heat_source=0;
	int r_cnt=0, offset =0, temp=0;
	double **final_matrix = alloc_matrix(ROWS); 
	init_matrix(ROWS, final_matrix,0,0);
	double cuu_max = 0.;
	MPI_Status status,status1,status2,status3,status4,status5;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Request send_req, rec_req, send_req_L, rec_req_L, send_req_R, rec_req_R;
	double max_values[size];
	int offset_arr[size];


	for (i = 0;i<size;i++){
    	max_values[i]=0.;  
    	offset_arr[i]=0;
    }

	if(rank==MASTER){
		//compute number of rows for each processor
		avgRows = ROWS/(size-1);

		//calculate remaining rows
		remainder = ROWS%(size-1);
		offset = 0;

		for(i = 1; i < size; i++){
			//for each worker send rows, add extra rows to first few processors
			/*
				The maximum number of rows owned by any process will
				be at most 1 more than the number of rows owned by a
				process with the least number of rows
			*/
			rows = (i <= remainder) ? avgRows+1 : avgRows;
			
			//save offset for each worker computed in earlier loop iteration
			offset_arr[i] = offset;

			//Tell each worker who its neighbors(left & right) are
			//special case: worker 1 and last worker
			//row 0 can't exchange with row 0-1
			if (i == 1) 
                left = NONE;
	        else
	            left = i - 1;
	        //last row can't exchange with last+1.
	        if (i == size-1)
	            right = NONE;
	        else
	            right = i + 1;

	        //send all re info to each worker
			MPI_Send(&rows, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&left, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&right, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&offset, 1, MPI_INT, i, BEGIN, MPI_COMM_WORLD);
			//printf("Sending to process: %d rows= %d with neighbours as left= %d right= %d\n",i,rows,left,right);
			//compute offset for next worker
			offset = offset + rows;

		}

		// continue working till heat is < eta
        while(1){
        	//printf("Waiting for receive from all processors\n");
	         //receive max heat value from workers in array
	        for(i = 1; i < size; i++){
	         	MPI_Recv(&max_values[i], 1, MPI_DOUBLE, i, WDONE, MPI_COMM_WORLD, &status);
	         	//printf("Received %f from rank: %d\n", i, max_values[i]);
	        }
	        cuu_max = max_in_array(max_values, size);
	        //printf("Current max is %f\n",cuu_max);
	        if( (cuu_max < EPS) && (cuu_max !=0) ){
	        	temp=0;
                int k,vals;
	        	printf("Current max is below the stopping eta\n");
	        	
	        	/*
				If stopping condition is met, send stopping flag to all workers.
				Receive row counts from respective workers 
				According to those row counts, add rows in the final matrix
	        	*/
				for(i = 1; i < size; i++){
	         		// get rows from all
	         		/* Wait for results from all worker tasks */
	         		MPI_Send(&temp, 1, MPI_INT, i, MDONE, MPI_COMM_WORLD);
					MPI_Recv(&r_cnt, 1, MPI_INT, i, FDONE, MPI_COMM_WORLD, &status);
					MPI_Status fStat;
	         		MPI_Recv(&final_matrix[offset_arr[i]][0], r_cnt*COLS, MPI_DOUBLE, i, FDONE, MPI_COMM_WORLD, &fStat);			         		
					//printf("Added entries in final_matrix for processor %d\n", i);
	        	}
	        	printf("Final Matrix is: \n");
	        	print_matrix(final_matrix,ROWS);
	        	MPI_Finalize();	
	        	break;
	        } else {
	        	/*
				If stopping condition is not met, ask workers to continue working. 
	        	*/
	        	temp=1;
	        	int i = 1;
	        	for(i = 1; i < size; i++){
					MPI_Send(&temp, 1, MPI_INT, i , MDONE, MPI_COMM_WORLD);
				}						
	        }
         }
	} else { 
			//if not MASTER
		    int heatFlag = 0;	
		    int l_index = 0, r_index = 0, proceedFlag = 0,k=0;	
		    int heat_source_proc=0, heat_source_I_FIX=0,start=0, end=0;
		 	// receive offset, rows, left, right
		    MPI_Status rw_status, l_status, r_status, o_status, pf_status;

			MPI_Recv(&rows, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &rw_status);
			MPI_Recv(&left, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &l_status);
			MPI_Recv(&right, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &r_status);
			MPI_Recv(&offset, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &o_status);
			//printf("At process:%d rows=%d with neighbours as left=%d right=%d\n",rank,rows,left,right);
			//printf("Processor:%d offset:%d",rank, offset);
			
			//allocate memory for the matrices
			double **a_old = alloc_matrix(rows+2); 
		    double **a_new = alloc_matrix(rows+2);
		    double **a_final = alloc_matrix(rows);
		    double max_diff = 1.0;

		    /*
		    Compute if the heat source-offset exists in the number of rows sent to the worker
			if it fits in the boundaries of rows for this worker, set the heatFlag which will
			help other functions to pinpoint the heat source
		    */
			heat_source_I_FIX = I_FIX - offset;

			if ((0 <= heat_source_I_FIX) && (heat_source_I_FIX < rows))
				heatFlag = 1;
			else 
				heatFlag = 0;
			//printf("At processor %d with heatFlag %d at co-ordinates[%d][%d]\n",rank,heatFlag,heat_source_I_FIX,J_FIX);
			//printf("Initializing matrices at processor %d\n", rank);
		    init_matrix(rows+2, a_old, heatFlag, heat_source_I_FIX); //initialize the matrices
		    init_matrix(rows+2, a_new, heatFlag, heat_source_I_FIX); 
		    init_matrix(rows, a_old, heatFlag, heat_source_I_FIX);
			
			start = 1; //index to send to left neighbor
			end = rows; //index to send to right neighbor

			r_index = rows+1;	//index to receive from right neighbor
			l_index = 0;		//index to receive from left neighbor
			
			 // Communicate border rows with neighbors
	        while(1){
		        //printf("left: %d right: %d",left,right);
		        //printf("Processor %d is sending max value: %f\n",rank,max_diff);
		        MPI_Isend(&max_diff,1,MPI_DOUBLE, MASTER, WDONE, MPI_COMM_WORLD, &send_req);
		        //if left neighbor exists i.e., is not 0 communicate border rows
		        if (left != NONE){
		        	//printf("Processor %d is sending row at %d to processor %d\n",rank,start,left);
		        	MPI_Isend(a_old[start], COLS, MPI_FLOAT, left, RGHOST, MPI_COMM_WORLD, &send_req_L);
		        	//printf("Processor %d is receiving row at %d from processor %d\n",rank,l_index,left);
		            MPI_Irecv(a_old[l_index], COLS, MPI_FLOAT, left, LGHOST, MPI_COMM_WORLD, &rec_req_L);

		            MPI_Wait(&send_req_L, &status2);
		            MPI_Wait(&rec_req_L, &status3);
	         	}

	         	//if right neighbor exists i.e., is not 0 communicate border rows
		        if (right != NONE){
		        	//printf("Processor %d is sending row at %d to processor %d\n",rank,end,right);
		        	MPI_Isend(a_old[end], COLS, MPI_FLOAT, right, LGHOST, MPI_COMM_WORLD, &send_req_R);
		        	//printf("Processor %d is receiving row at %d from processor %d\n",rank,r_index,right);
		            MPI_Irecv(a_old[r_index], COLS, MPI_FLOAT, right, RGHOST, MPI_COMM_WORLD, &rec_req_R);

		            MPI_Wait(&send_req_R, &status4);
		            MPI_Wait(&rec_req_R, &status5);
	         	}

	         	//call compute_new_values to update values of grid points
	         	compute_new_values(rows+2, a_old, a_new, heatFlag, heat_source_I_FIX,rank,size);
	         	
	         	//calculate the maximum absolute differences among pairwise
				// differences of old and new matrix elements
	         	max_diff = max_abs(rows+2, a_new, a_old);
	         	
	         	//Compute maximum value and send to MASTER
	         	MPI_Wait(&send_req, &status1);

	         	//Wait for proceedFlag from MASTER
	         	MPI_Recv(&proceedFlag, 1, MPI_INT, MASTER, MDONE, MPI_COMM_WORLD, &pf_status);
	         	//printf("Processor %d received proceedFlag: \n",rank, proceedFlag);
	         	
	         	//send this worker's portion of final results back to master if proceedFlag is 0
		        if(proceedFlag !=1){
					for (int i = 1; i < rows+1; i++)
				        for (int j = 0; j < COLS; j++)
				            a_final[i-1][j] = a_old[i][j];
				    //printf("-------------------Sending to MASTER------------------------\n");
		        	MPI_Send(&rows, 1, MPI_INT, MASTER, FDONE, MPI_COMM_WORLD);
		        	MPI_Send(a_final[0], rows*COLS, MPI_DOUBLE, MASTER, FDONE, MPI_COMM_WORLD);
	    		    MPI_Finalize();
	    		    break;
		        }else{
		        	//continue with the computation, create copy of matrix
		        	copy_matrix(rows+2, a_old, a_new, heatFlag, heat_source_I_FIX);
		        }
		    }
		}
	}
