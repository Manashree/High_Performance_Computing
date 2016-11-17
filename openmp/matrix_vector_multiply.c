#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char* argv[])
{
    const int n = 3;
    int i, j, sum;
    double b[n] = { 1, 2, 3 };
    double a[n][n] = { { 1, 2, 3 }, { 1, 2, 3 }, { 1, 2, 3 } };
    double c[n] = { 0, 0, 0 };

    printf("Input matrix A\n");
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
                printf("%f", a[i][j]);
        }
        printf("\n");
    }
    printf("Input vector B\n");
    for (i = 0; i < 3; i++) {
        printf("%f", b[i]);
        printf("\n");
    }
   #pragma omp parallel{
        double sum2 = 0;
        int i2, j2;
        for (i2 = 0; i2 < n; i2++) {
        #pragma omp for
            for (j2 = 0; j2 < n; j2++) {
                sum2 = sum2 + (a[i2][j2] * b[j2]);
            }
        #pragma omp critical{
                c[i2] = c[i2] + sum2;
                sum2 = 0;
            }
        }
    }
    printf("Result is:\n");
    for (i = 0; i < n; i++) {
        printf("%f", c[i]);
        printf("\n");
    }
}
