#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "time_benchmark.h"
#include "linearAlg.h"
//inv(X'*X)*X'*Y

int main()
{
    start_tracking_time();

    Matrix* a = la_initRandMatrix(5, 10);
    Matrix* b = la_initRandMatrix(10, 5);

    Matrix* c = la_multiplyMM(a, b);
    la_freeMatrix(a);
    la_freeMatrix(b);
    la_inverse(c);
    
    Vector* v = la_initRandVector(5);
    Vector* v1 = la_multiplyMV(c, v);
    la_freeVector(v);
    la_printVector(v1);

    end_tracking_time();
    printf("\ntime: %.5f ms\n", get_time_taken());

    printf("fin");
    return 0;
}


