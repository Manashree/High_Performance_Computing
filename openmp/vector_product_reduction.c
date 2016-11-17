#include <stdio.h>
#include <omp.h> 

/* Calculate inner product of two vectors */

int main(){
  int const N=100;
  int i, k;
  double a[N], b[N]; /* declare vectors a and b */
  double dot_prod = 0.0; 

  for(i = 0; i < N; i++) {
    a[i] = 3.14;
    b[i] = 6.67;
  }

#pragma omp parallel{ 
/*
Reduction: 
 * Specifies that one or more variables that are private to each thread are the subject of 
 * a reduction operation at the end of the parallel region. Syntax : reduction(operation:var)
 * dot_prod is the reduction variable and hence cannot be declared shared because the 
 * threads would overwrite the value of dot_prod and it also 
 * cannot be declared private as private variables don't persist outside of parallel region. 
 * It has to be specified as a reduction operation as it is performed on individual values 
 * from each thread 
*/
  #pragma omp for reduction(+:dot_prod)
    for(i = 0; i<N; i++) /* compute dot  product */
      dot_prod = dot_prod + a[i] * b[i];
  }
  printf("Inner product of a[] and b[] is %g\n", dot_prod);
  return 0;
}
