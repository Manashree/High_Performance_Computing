#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]){
        const int N = 100;
        int nthreads, threadid, i;
        double a[N], b[N], c[N];

        // Initialize
        for (i=0; i < N; i++) {
                a[i] = 1.0*i;
                b[i] = 2.0*i;
        }
        //fork-join
        #pragma omp parallel shared(a,b,c,nthreads) private(i,threadid){
                threadid = omp_get_thread_num();
                if (threadid == 0) {
                        nthreads = omp_get_num_threads();
                        printf("Number of threads = %d\n", nthreads);
                }
                printf(" My threadid %d\n",threadid);
                #pragma omp for 
                for (i=0; i<N; i++){
                        c[i] = a[i] + b[i];
                }
        } // join  
        printf(" TEST c[99] = %g\n",c[N-1]);
        return 0;
}
