#ifndef LINEAR_ALG
#define LINEAR_ALG

#define LA_MAX_INVERSE_SIZE 10

typedef struct matrix
{
    int rows;
    int columns;
    float** entries;
} Matrix;

typedef struct vector
{
    int rows;
    float* entries;
} Vector;

//la for linear algebra

//=============== vectors ===============

void la_printVector(Vector *vector);
void laa_writeVectorBin(float *vector, int rows, FILE* filePointer);
void laa_readVectorBin(float* destVector, FILE* filePointer);
Vector *la_cloneVector(Vector *vector);
Vector *la_initVector(int rows, float intialValue);
Vector *la_initVectorArray(float *entries, int rows);
Vector *la_initZerosVector(int rows);
Vector *la_initRandVector(int rows);
void la_freeVector(Vector *vector);

// vector math:

float la_dot(Vector *a, Vector *b);
void la_multiplyMVTo(Matrix *a, Vector *b, Vector *destination);
Vector *la_multiplyMV(Matrix *a, Vector *b);
void la_addVectorsTo(Vector *a, Vector *b, Vector *destination);
Vector *la_addVectors(Vector *a, Vector *b);
void la_subtractVectorsTo(Vector *a, Vector *b, Vector *destination);
Vector *la_subtractVectors(Vector *a, Vector *b);

//=============== matrices ===============

void la_printMatrix(Matrix *matrix);
void laa_writeMatrixBin(float** matrix, int rows, int columns, FILE* filePointer);
void laa_readMatrixBin(float** destMatrix, FILE* filePointer);
Matrix *la_cloneMatrix(Matrix *matrix);
Matrix *la_initMatrix(int rows, int columns, float initialValue);
Matrix *la_initMatrixArray(float **entries, int rows, int columns);
Matrix *la_initZerosMatrix(int rows, int columns);
Matrix *la_initRandMatrix(int rows, int columns);
void la_freeMatrix(Matrix *matrix);

// matrix math:

void la_addToMatrix(Matrix *a, Matrix *b);
void la_subtractFromMatrix(Matrix *a, Matrix *b);
Matrix *la_addMatrix(Matrix *a, Matrix *b);
void la_addMatrixTo(Matrix *a, Matrix *b, Matrix *destination);
Matrix *la_subtractMatrix(Matrix *a, Matrix *b);
void la_subtractMatrixTo(Matrix *a, Matrix *b, Matrix *destination);
Matrix *la_transposed(Matrix *matrix);
Matrix *la_multiplyMM(Matrix *a, Matrix *b);
void la_multiplyMMTo(Matrix *a, Matrix *b, Matrix *destination);
void la_inverse(Matrix *matrix);
Matrix *la_getInverse(Matrix *matrix);

//=====================================================================================================================================
//=================================================== Linear algebra with arrays ======================================================
//=====================================================================================================================================

//=============== vectors ===============

void laa_printVector(float *vector, int columns);
float *laa_allocVector(int rows, float initialValue);
float *laa_allocVectorRaw(int rows);
float *laa_allocRandVector(int rows);
float *laa_cloneVector(float *cloneVector, int rows);
void laa_freeVector(float *vector);
void laa_setVector(float *vector, int rows, float value);
void laa_copyVectorValues(float *copyFrom, float *pasteTo, int rows);

// vector math:

float laa_dot(float *a, float *b, int rows);
void laa_addVectorsTo(float *a, float *b, float *destination, int rows);
void laa_subtractVectorsTo(float *a, float *b, float *destination, int rows);
float *laa_addVectors(float *a, float *b, int rows);
float *laa_subtractVectors(float *a, float *b, int rows);
void laa_multiplyMVTo(float **matrix, int rows, int columns, float *vector, float *destination);
float *laa_multiplyMV(float **matrix, int rows, int columns, float *vector);

//=============== matrices ===============

void laa_printMatrix(float **matrix, int rows, int columns);
float **laa_allocMatrix(int rows, int columns, float initialValue);
float **laa_allocMatrixRaw(int rows, int columns);
float **laa_allocRandMatrix(int rows, int columns);
float **laa_cloneMatrix(float **values, int rows, int columns);
float** cloneMatrix2dArray(int rows, int columns, float values[rows][columns]);
void laa_freeMatrix(float **matrix, int rows);
void laa_setMatrix(float **matrix, int rows, int columns, float value);
void laa_copyMatrixValues(float **copyFrom, float **pasteTo, int rows, int columns);

// matrix math:

void laa_addMatrixTo(float **a_matrixVals, float **b_matrixVals, float **destination, int rows, int columns);
float **laa_addMatrix(float **a_matrixVals, float **b_matrixVals, int rows, int columns);
void laa_subtractMatrixTo(float **a_matrixVals, float **b_matrixVals, float **destination, int rows, int columns);
float **laa_subtractMatrix(float **a_matrixVals, float **b_matrixVals, int rows, int columns);
float **laa_matrixDifference(float **a_matrixVals, float **b_matrixVals, int rows, int columns);
void laa_transposeTo(float **matrix, float **destination, int rows, int columns);
float **laa_transpose(float **matrix, int rows, int columns);
void laa_multiplyMMTo(float **a_matrixVals, int a_rows, int a_columns, float **b_matrixVals, int b_rows, int b_columns, float **destination);
void laa_multiplyMMReplace(float **a_matrixVals, int a_rows, int a_columns, float **b_matrixVals, int b_rows, int b_columns, float **destination);
float **laa_multiplyMM(float **a_matrixVals, int a_rows, int a_columns, float **b_matrixVals, int b_rows, int b_columns);
void laa_inverse(float **matrixEntries, int rows, int columns);
float **laa_getInverse(float **matrixEntries, int rows, int columns);

#endif