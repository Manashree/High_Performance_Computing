#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
        const int N = 20;
        int nthreads,threadid,i;
        double a[N], b[N], result[N];
        /* Some initializations */
        int chunk = 5;
        for(i = 0; i < N; i++) {
                a[i] = i * 1.0;
                b[i] = i * 2.0;
        }
        #pragma omp parallel private(threadid){
            // fork
            threadid = omp_get_thread_num();
            #pragma omp for schedule(dynamic, chunk)
            for(i = 0; i < N; i++) {
		        result[i] = a[i] + b[i];
		        printf("Thread id: %d working on index %d \n", threadid, i);
            }
        } // join
        printf("TEST result last index result[19]= %g\n", result[N-1]);
        return 0;
}
