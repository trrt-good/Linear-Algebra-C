#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include "..\\LinearAlg\\linearAlg.h"
//inv(X'*X)*X'*Y

int main()
{
    LARGE_INTEGER tps;
    LARGE_INTEGER t1, t2;
    float timeDiff;
    QueryPerformanceFrequency(&tps);
    QueryPerformanceCounter(&t1);

    srand(time(NULL));

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

    QueryPerformanceCounter(&t2);
    timeDiff = (t2.QuadPart - t1.QuadPart) * 1000.0 / tps.QuadPart;
    printf("\ntime: %.5f ms\n", timeDiff);

    printf("fin");
    return 0;
}


