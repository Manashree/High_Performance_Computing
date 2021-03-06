#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

int main(){
  const int N = 100;
  int x[N], i, max_x, min_x, sum,sum2;
  float mean,mean2,var;
  max_x = 0;
  min_x = 100;
  sum = 0;
  sum2 = 0;

  /* initialize x */
  srand(time(NULL)); // Initialize random variable seed
  #pragma omp parallel for 
  for(i = 0; i < N; i++) { 
    x[i] = rand();
  }


#pragma omp parallel private(i) shared(x){
#pragma omp sections{
  /* fork 3 different threads */
      { 
        for(i = 0; i < N; i++) { /* find min & max of x */
          if (x[i] > max_x) max_x = x[i];
          if (x[i] < min_x) min_x = x[i];
        }
        printf("The max of x = %d\n", max_x);
        printf("The min of x = %d\n", min_x);
      }
      
      #pragma omp section{ /* calculate the mean of x */
        for(i = 0; i < N; i++)
          sum = sum + x[i];
        mean = sum/N;
        printf("Mean of x = %f\n", mean);
      }
#pragma omp section {
        for(i = 0; i < N; i++)
          sum2 = sum2 + x[i]*x[i];
        mean2 = sum2/N;
      }
  }
}
  var = mean2 - mean*mean;
  printf("Variance of x = %f\n",var);
  return 0;
}
