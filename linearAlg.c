#include "linearAlg.h"

//=============== vectors ===============
void la_printVector(Vector *vector)
{
    laa_printVector(vector->entries, vector->rows);
}

Vector *la_cloneVector(Vector *vector)
{
    return la_initVectorArray(laa_cloneVector(vector->entries, vector->rows), vector->rows);
}

Vector *la_initVector(int rows, float intialValue)
{
    Vector *newVector = (Vector *)malloc(sizeof(Vector));
    newVector->rows = rows;
    newVector->entries = laa_allocVector(rows, intialValue);
    return newVector;
}

Vector *la_initVectorArray(float *entries, int rows)
{
    Vector *newVector = (Vector *)malloc(sizeof(Vector));
    newVector->rows = rows;
    newVector->entries = entries;
    return newVector;
}

Vector *la_initZerosVector(int rows)
{
    return la_initVector(rows, 0);
}

Vector *la_initRandVector(int rows, float from, float to)
{
    return la_initVectorArray(laa_allocRandVector(rows, from, to), rows);
}

void la_freeVector(Vector *vector)
{
    free(vector->entries);
    free(vector);
}

// vector math:

float la_dot(Vector *a, Vector *b)
{
    if (a->rows != b->rows)
        printf("cannot preform dot product. Vectors must be the same length");

    return laa_dot(a->entries, b->entries, a->rows);
}

void la_multiplyMVTo(Matrix *a, Vector *b, Vector *destination)
{
    laa_multiplyMVTo(a->entries, a->rows, a->columns, b->entries, destination->entries);
}

Vector *la_multiplyMV(Matrix *a, Vector *b)
{
    return la_initVectorArray(laa_multiplyMV(a->entries, a->rows, a->columns, b->entries), a->columns);
}

void la_addVectorsTo(Vector *a, Vector *b, Vector *destination)
{
    laa_addVectorTo(a->entries, b->entries, destination->entries, a->rows);
}

Vector *la_addVectors(Vector *a, Vector *b)
{
    return la_initVectorArray(laa_addVectors(a->entries, b->entries, a->rows), a->rows);
}

void la_subtractVectorsTo(Vector *a, Vector *b, Vector *destination)
{
    laa_subtractVectorsTo(a->entries, b->entries, destination->entries, a->rows);
}

Vector *la_subtractVectors(Vector *a, Vector *b)
{
    return la_initVectorArray(laa_subtractVectors(a->entries, b->entries, a->rows), a->rows);
}

//=============== matrices ===============

void la_printMatrix(Matrix *matrix)
{
    laa_printMatrix(matrix->entries, matrix->rows, matrix->columns);
}

Matrix *la_cloneMatrix(Matrix *matrix)
{
    Matrix *newMatrix = (Matrix *)malloc(sizeof(Matrix));

    newMatrix->rows = matrix->rows;
    newMatrix->columns = matrix->columns;

    newMatrix->entries = laa_cloneMatrix(matrix->entries, matrix->rows, matrix->columns);
    return newMatrix;
}

Matrix *la_initMatrix(int rows, int columns, float initialValue)
{
    Matrix *newMatrix = (Matrix *)malloc(sizeof(Matrix));
    newMatrix->rows = rows;
    newMatrix->columns = columns;

    newMatrix->entries = laa_allocMatrix(rows, columns, initialValue);
    return newMatrix;
}

Matrix *la_initMatrixArray(float **entries, int rows, int columns)
{
    Matrix *newMatrix = (Matrix *)malloc(sizeof(Matrix));
    newMatrix->rows = rows;
    newMatrix->columns = columns;
    newMatrix->entries = entries;
    return newMatrix;
}

Matrix *la_initZerosMatrix(int rows, int columns)
{
    int i, j;
    Matrix *newMatrix = (Matrix *)malloc(sizeof(Matrix));
    newMatrix->rows = rows;
    newMatrix->columns = columns;

    newMatrix->entries = laa_allocMatrix(rows, columns, 0);
    return newMatrix;
}

Matrix *la_initRandMatrix(int rows, int columns, float from, float to)
{
    int i, j;
    Matrix *newMatrix = (Matrix *)malloc(sizeof(Matrix));
    newMatrix->rows = rows;
    newMatrix->columns = columns;

    newMatrix->entries = laa_allocRandMatrix(rows, columns, from, to);
    return newMatrix;
}

void la_freeMatrix(Matrix *matrix)
{
    int i;
    for (i = 0; i < matrix->rows; i++)
    {
        free(matrix->entries[i]);
    }
    free(matrix->entries);
    free(matrix);
    return;
}

// matrix math:

/**
 * @brief adds matrix b to matrix a. Matrix b is unaffected, but matrix a
 * will be equal to a+b. This does not create any new matrices
 * @param a a pointer to the matrix to recieve the operation (end result is a+b)
 * @param b a pointer to matrix to be added to matrix a
 */
void la_addToMatrix(Matrix *a, Matrix *b)
{
    if (a->columns != b->columns | a->rows != b->rows)
    {
        printf("Matrices have incompatable dimensions! cannot add");
        exit(0);
    }
    laa_addMatrixTo(a->entries, b->entries, a->entries, a->rows, a->columns);
}

/**
 * @brief subtracts matrix b from matrix a. Matrix b is unaffected, but matrix a
 * will be equal to a-b. This does not create any new matrices
 * @param a pointer to the matrix to recieve the operation (end result is a-b)
 * @param b pointer to the matrix to be subtracted from matrix a
 */
void la_subtractFromMatrix(Matrix *a, Matrix *b)
{
    if (a->columns != b->columns | a->rows != b->rows)
    {
        printf("Matrices have incompatable dimensions! cannot add");
        exit(0);
    }
    laa_subtractMatrixTo(a->entries, b->entries, a->entries, a->rows, a->columns);
}

Matrix *la_addMatrix(Matrix *a, Matrix *b)
{
    if (a->columns != b->columns | a->rows != b->rows)
    {
        printf("Matrices have incompatable dimensions! cannot add");
        exit(0);
    }
    la_initMatrixArray(laa_addMatrix(a->entries, b->entries, a->rows, a->columns), a->rows, a->columns);
}

void la_addMatrixTo(Matrix *a, Matrix *b, Matrix *destination)
{
    if (a->columns != b->columns | a->rows != b->rows)
    {
        printf("Matrices have incompatable dimensions! cannot add");
        exit(0);
    }
    laa_addMatrixTo(a->entries, b->entries, destination->entries, a->rows, b->rows);
}

Matrix *la_subtractMatrix(Matrix *a, Matrix *b)
{
    if (a->columns != b->columns | a->rows != b->rows)
    {
        printf("Matrices have incompatable dimensions! cannot add");
        exit(0);
    }
    la_initMatrixArray(laa_subtractMatrix(a->entries, b->entries, a->rows, a->columns), a->rows, a->columns);
}

void la_subtractMatrixTo(Matrix *a, Matrix *b, Matrix *destination)
{
    if (a->columns != b->columns | a->rows != b->rows)
    {
        printf("Matrices have incompatable dimensions! cannot add");
        exit(0);
    }
    laa_subtractMatrixTo(a->entries, b->entries, destination->entries, a->rows, a->columns);
}

/**
 * @brief creates a new matrix that is equal to the transpose of matrix
 * @param matrix pointer to the matrix that will be used to create
 * a transposed version of the matrix
 * @return Matrix* a pointer to the new matrix
 */
Matrix *la_transposed(Matrix *matrix)
{
    return la_initMatrixArray(laa_transpose(matrix->entries, matrix->rows, matrix->columns), matrix->columns, matrix->rows);
}

/**
 * @brief creates a new matrix equal to a*b
 * @param a pointer to a matrix
 * @param b pointer to a matrix
 * @return Matrix* a pointer to the new matrix equal to a*b
 */
Matrix *la_multiplyMM(Matrix *a, Matrix *b)
{
    return la_initMatrixArray(laa_multiplyMM(a->entries, a->rows, a->columns, b->entries, b->rows, b->columns), a->rows, b->columns);
}

void la_multiplyMMTo(Matrix *a, Matrix *b, Matrix *destination)
{
    if (a->rows != destination->rows || b->columns != destination->columns)
        laa_multiplyMMTo(a->entries, a->rows, a->columns, b->entries, b->rows, b->columns, destination->entries);
}

/**
 * @brief inverses the given matrix
 * @param matrix the matrix to be inversed
 * @return Matrix* the same pointer to the given matrix
 */
void la_inverse(Matrix *matrix)
{
    laa_inverse(matrix->entries, matrix->rows, matrix->columns);
}

Matrix *la_getInverse(Matrix *matrix)
{
    return la_initMatrixArray(laa_getInverse(matrix->entries, matrix->rows, matrix->columns), matrix->rows, matrix->columns);
}

//=====================================================================================================================================
//=================================================== Linear algebra with arrays ======================================================
//=====================================================================================================================================

#pragma region laa - linear algebra with arrays

//=============== vectors ===============

void laa_printVector(float *vector, int columns)
{
    int i;
    printf("[ ");
    for (i = 0; i < columns - 1; i++)
    {
        printf("%.13f, ", vector[i]);
    }
    printf("%.13f ]\n", vector[columns - 1]);
}

void laa_writeVectorBin(float *vector, int rows, FILE *filePointer)
{
    int elements_written;
    elements_written = fwrite(&rows, sizeof(int), 1, filePointer);
    elements_written += fwrite(vector, sizeof(float), rows, filePointer);
    if (elements_written != rows + 1)
    {
        printf("error while writing vector");
        exit(1);
    }
}

void laa_readVectorBin(float *destVector, FILE *filePointer)
{
    int rows;
    fread(&rows, sizeof(int), 1, filePointer);
    fread(destVector, sizeof(float), rows, filePointer);
}

/**
 * @brief allocates a vector
 *
 * @param rows
 * @param initialValue
 * @return float*
 */
float *laa_allocVector(int rows, float initialValue)
{
    int i;
    float *vector = (float *)malloc(rows * sizeof(float));
    if (!vector)
        exit(1);
    if (initialValue == 0.0)
        memset(vector, 0, rows*sizeof(float));
    else
        for (i = 0; i < rows; i++)
            vector[i] = initialValue;
    return vector;
}

/**
 * @brief allocates a vector with length rows
 * without initializing any of the values
 *
 * @param rows size of the vector
 * @return float* pointer to the vector
 */
float *laa_allocVectorRaw(int rows)
{
    return (float *)malloc(rows * sizeof(float));
}

/**
 * @brief allocates a vector of the specified length
 *
 * @param rows length of the vector
 * @return float* pointer to the vector in memory
 */
float *laa_allocRandVector(int rows, float from, float to)
{
    int i;
    float *vector = (float *)malloc(rows * sizeof(float));
    if (!vector)
        exit(1);
    for (i = 0; i < rows; i++)
    {
        vector[i] = ((float)rand() / (float)(RAND_MAX))*(to-from)+from;
    }
    return vector;
}

float *laa_cloneVector(float *cloneVector, int rows)
{
    int i;
    float *vector = (float *)malloc(rows * sizeof(float));
    if (!vector)
        exit(1);
    for (i = 0; i < rows; i++)
    {
        vector[i] = cloneVector[i];
    }
    return vector;
}

/**
 * @brief frees a vector from memory
 *
 * @param vector pointer to the vector in memory
 */
void laa_freeVector(float *vector)
{
    free(vector);
    return;
}

void laa_setVectorTo(float *vector, int rows, float value)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        vector[i] = value;
    }
}

void laa_setVectorToRand(float *vector, int rows)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        vector[i] = ((float)rand()) / ((float)RAND_MAX) - 0.5;
    }
}

void laa_copyVectorValues(float *copyFrom, float *pasteTo, int rows)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        pasteTo[i] = copyFrom[i];
    }
}

int laa_maxIndexValue(float* vector, int rows)
{
    int max = 0;
    for (int i = 0; i < rows; i ++)
    {
        if (vector[i] > vector[max])
            max = i;
    }
    return max;
}

// vector math:

/**
 * @brief dot product between arrays a and b
 *
 * @param a float array
 * @param b float array
 * @param rows the length of the arrays a and b
 * @return float dot product value
 */
float laa_dot(float *a, float *b, int rows)
{
    int i;
    float sum = 0;
    for (i = rows; i--;)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

void laa_addVectorTo(float *a, float *b, float *destination, int rows)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        destination[i] = a[i] + b[i];
    }
}

void laa_subtractVectorsTo(float *a, float *b, float *destination, int rows)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        destination[i] = a[i] - b[i];
    }
}

float *laa_addVectors(float *a, float *b, int rows)
{
    int i;
    float *newVector = laa_allocVectorRaw(rows);
    for (i = 0; i < rows; i++)
    {
        newVector[i] = a[i] + b[i];
    }
    return newVector;
}

float *laa_subtractVectors(float *a, float *b, int rows)
{
    int i;
    float *newVector = laa_allocVectorRaw(rows);
    for (i = 0; i < rows; i++)
    {
        newVector[i] = a[i] - b[i];
    }
    return newVector;
}

/**
 * @brief MV: Matrix-Vector multiplciation. for simplicity, it is
 * assumed that the number of rows of the vector is equal to the
 * number of columns of the matrix (however this means that there
 * can be memory issues if the vector is too small)
 * @param vector
 * @param matrix
 * @param rows rows of the matrix
 * @param columns columns of the matrix, equal to rows of the vector
 */
void laa_multiplyMVTo(float **matrix, int rows, int columns, float *vector, float *destination)
{
    if (destination == vector)
    {
        printf("destination of the matrix-vector multiplication cannot be the vector itself\n");
        exit(1);
    }

    int i;
    for (i = 0; i < rows; i++)
    {
        destination[i] = laa_dot(vector, matrix[i], columns);
    }
}

float *laa_multiplyMV(float **matrix, int rows, int columns, float *vector)
{
    int i;
    float *newVector = laa_allocVectorRaw(rows);
    for (i = 0; i < rows; i++)
    {
        newVector[i] = laa_dot(vector, matrix[i], columns);
    }
    return newVector;
}

void laa_multiplyVSTo(float* vector, int rows, float multiplier, float* destination)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        destination[i] = vector[i]*multiplier;
    }
}

float* laa_multiplyVS(float* vector, int rows, float multiplier)
{
    int i;
    float* new_vector = malloc(sizeof(float)*rows);
    for (i = 0; i < rows; i ++)
    {
        new_vector[i] = vector[i]*multiplier;
    }
    return new_vector;
}

//=============== matrices ===============

/**
 * @brief prints a 2d array
 *
 * @param matrix pointer to the array in memory
 * @param rows rows of the matrix
 * @param columns columns of the matrix
 */
void laa_printMatrix(float **matrix, int rows, int columns)
{
    int i, j;
    if (rows <= 0 | columns <= 0)
    {
        printf("invalid matrix size! cannot print");
    }

    printf("\n      ");
    for (j = 0; j < columns; j++)
    {
        printf("%20d ", j);
    }
    for (i = 0; i < rows; i++)
    {
        printf("\n%5d |", i);
        for (j = 0; j < columns; j++)
        {
            printf("%20.13f", matrix[i][j]);
        }
        printf("    |");
    }
    printf("\n\n");
}

void laa_writeMatrixBin(float **matrix, int rows, int columns, FILE *filePointer)
{
    int i = 0;
    fwrite(&rows, sizeof(int), 1, filePointer);
    fwrite(&columns, sizeof(int), 1, filePointer);
    for (i = 0; i < rows; i++)
    {
        fwrite(matrix[i], sizeof(float), columns, filePointer);
    }
}

void laa_readMatrixBin(float **destMatrix, FILE *filePointer)
{
    int i = 0, rows, columns;
    fread(&rows, sizeof(int), 1, filePointer);
    fread(&columns, sizeof(int), 1, filePointer);
    for (i = 0; i < rows; i++)
    {
        fread(destMatrix[i], sizeof(float), columns, filePointer);
    }
}

/**
 * @brief allocates an array of arrays which has
 * dimensions rows x columns. The length of the pointer
 * of pointer array is the rows and the length of
 * each single pointer array is the columns
 *
 * @param rows number of rows of the matrix
 * @param columns number of columns of the matrix
 * @param initialValue the initial value for all entries in the matrix
 * @return float** pointer to the array of arrays
 */
float **laa_allocMatrix(int rows, int columns, float initialValue)
{
    int i, j;

    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < columns; j++)
        {
            matrixEntries[i][j] = initialValue;
        }
    }
    return matrixEntries;
}

/**
 * @brief allocates an array of arrays wihtout
 * intitalizing values, which has dimensions rows
 * x columns. The length of the pointer
 * of pointer array is the rows and the length of
 * each single pointer array is the columns.
 *
 * O(rows)
 *
 * @param rows number of rows of the matrix
 * @param columns number of columns of the matrix
 */
float **laa_allocMatrixRaw(int rows, int columns)
{
    int i, j;

    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
    {
        printf("allocation error! 1");
        exit(1);
    }
        

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            printf("allocation error! 2");
            exit(1);
        }
    }
    return matrixEntries;
}

/**
 * @brief allocates a matrix of size rows x columns in memory
 * with random float values between 0.0 and 1.0
 *
 * @param rows rows of the matrix
 * @param columns columns of the matrix
 * @return float** pointer to the matrix in memory
 */
float **laa_allocRandMatrix(int rows, int columns, float from, float to)
{
    int i, j;

    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < columns; j++)
        {
            matrixEntries[i][j] = ((float)rand() / (float)(RAND_MAX))*(to-from)+from;
        }
    }
    return matrixEntries;
}

/**
 * @brief clones a matrix given the rows, columns
 * and the matrix entries.
 *
 * @param rows rows of the matrix that will be cloned
 * @param columns columns of the matrix to be cloned
 * @param values values of the matrix to be cloned
 * @return float** the cloned matrix
 */
float **laa_cloneMatrix(float **values, int rows, int columns)
{
    int i, j;
    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < columns; j++)
        {
            matrixEntries[i][j] = values[i][j];
        }
    }
    return matrixEntries;
}

float **cloneMatrix2dArray(int rows, int columns, float values[rows][columns])
{
    int i, j;

    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < columns; j++)
        {
            matrixEntries[i][j] = values[i][j];
        }
    }
    return matrixEntries;
}

/**
 * @brief frees an array of arrays from memory.
 *
 * @param matrix the pointer to the matrix in memory
 * @param rows number of rows of the matrix, or length of the pointer to pointer array
 */
void laa_freeMatrix(float **matrix, int rows)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
    return;
}

void laa_setMatrixTo(float **matrix, int rows, int columns, float value)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            matrix[i][j] = value;
        }
    }
}

void laa_setMatrixToRand(float **matrix, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            matrix[i][j] = ((float)rand() / (float)(RAND_MAX)) - 0.5;
        }
    }
}

/**
 * @brief copies values of matrix to
 *
 * @param rows
 * @param columns
 * @param copyFrom
 * @param pasteTo
 */
void laa_copyMatrixValues(float **copyFrom, float **pasteTo, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            pasteTo[i][j] = copyFrom[i][j];
        }
    }
}

// math section:

/**
 * @brief adds a_matrix to b_matrix
 *
 * @param rows number of rows in both a and b
 * @param columns number of columns in both a and b
 * @param a_matrixVals values of matrix a
 * @param b_matrixVals values of matrix b
 * @param destination WARNING: MUST ALREADY BE ALLOCATED.
 * location for the output values, which can be the same
 * pointer as either matrix a or b
 */
void laa_addMatrixTo(float **a_matrixVals, float **b_matrixVals, float **destination, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            destination[i][j] = a_matrixVals[i][j] + b_matrixVals[i][j];
        }
    }
}

/**
 * @brief returns the pointer to a new matrix that is the
 * sum of matrix a and matrix b
 *
 * @param a_matrixVals the values of matrix a
 * @param b_matrixVals the values of matrix b
 * @param destination pointer to the destination matrix
 * @param rows
 * @param columns
 * @return pointer to the matrix a+b
 */
float **laa_addMatrix(float **a_matrixVals, float **b_matrixVals, int rows, int columns)
{
    int i, j;

    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < columns; j++)
        {
            matrixEntries[i][j] = a_matrixVals[i][j] + b_matrixVals[i][j];
        }
    }
    return matrixEntries;
}

/**
 * @brief subtracts b_matrix from a_matrix (a-b)
 *
 * @param rows number of rows in both a and b
 * @param columns number of columns in both a and b
 * @param a_matrixVals values of matrix a
 * @param b_matrixVals values of matrix b
 * @param destination WARNING: MUST ALREADY BE ALLOCATED.
 * location for the output values, which can be the same
 * pointer as either matrix a or b
 */
void laa_subtractMatrixTo(float **a_matrixVals, float **b_matrixVals, float **destination, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            destination[i][j] = a_matrixVals[i][j] - b_matrixVals[i][j];
        }
    }
}

float **laa_subtractMatrix(float **a_matrixVals, float **b_matrixVals, int rows, int columns)
{
    int i, j;

    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < columns; j++)
        {
            matrixEntries[i][j] = a_matrixVals[i][j] - b_matrixVals[i][j];
        }
    }
    return matrixEntries;
}

/**
 * @brief returns the pointer to a new matrix that is the
 * difference of matrix a and matrix b (a-b)
 *
 * @param a_matrixVals the values of matrix a
 * @param b_matrixVals the values of matrix b
 * @param destination pointer to the destination matrix
 * @param rows
 * @param columns
 * @return pointer to the matrix a-b
 */
float **laa_matrixDifference(float **a_matrixVals, float **b_matrixVals, int rows, int columns)
{
    int i, j;

    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < columns; j++)
        {
            matrixEntries[i][j] = a_matrixVals[i][j] - b_matrixVals[i][j];
        }
    }
    return matrixEntries;
}

/**
 * @brief transposes the given matrix into the specified destination,
 * which must be different than the pointer to the matrix to be transposed.
 * Make sure that the result location is big enough to hold the transposed
 * matrix (otherwise fear a memory issue)
 *
 * @param rows number of rows of the matrix to be transposed
 * @param columns number of columns of the matrix to be transposed
 * @param matrix a pointer to the locaion of the matrixto be transposed.
 * @param resultLocation the location for the transposed matrix
 * to go (must already be allocated)
 */
void laa_transposeTo(float **matrix, float **destination, int rows, int columns)
{
    int i, j;
    for (i = 0; i < columns; i++)
    {
        for (j = 0; j < rows; j++)
        {
            destination[i][j] = matrix[j][i];
        }
    }
}

float **laa_transpose(float **matrix, int rows, int columns)
{
    int i, j;

    float **matrixEntries = (float **)malloc(columns * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < columns; i++)
    {
        matrixEntries[i] = (float *)malloc(rows * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < rows; j++)
        {
            matrixEntries[i][j] = matrix[j][i];
        }
    }
    return matrixEntries;
}

void laa_multiplyMSTo(float** matrixVals, int rows, int cols, float multiplier, float** destination)
{
    int i, j;
    for (i = 0; i < rows; i ++)
    {
        for (j = 0; j < cols; j++)
        {
            destination[i][j] = matrixVals[i][j]*multiplier;
        }
    }
}

float** laa_multiplyMS(float** matrixVals, int rows, int cols, float multiplier)
{
    int i, j;
    float **matrixEntries = (float **)malloc(rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < rows; i++)
    {
        matrixEntries[i] = (float *)malloc(cols * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < cols; j++)
        {
            matrixEntries[i][j] = matrixVals[i][j]*multiplier;
        }
    }
    return matrixEntries;
}

/**
 * @brief
 *
 * @param a_matrixVals
 * @param a_rows
 * @param a_columns
 * @param b_matrixVals
 * @param b_rows
 * @param b_columns
 * @param destination
 */
void laa_multiplyMMTo(float **a_matrixVals, int a_rows, int a_columns, float **b_matrixVals, int b_rows, int b_columns, float **destination)
{
    int i, j, k;
    float sum = 0;
    if (a_columns != b_rows)
    {
        printf("matrices must have compatable dimensions!\n");
        exit(0);
    }
    else if (destination == a_matrixVals | destination == b_matrixVals)
    {
        printf("destination cannot be equal to one of the multipliers! please use a different matrix to store the product result or use the laa_multiplyMMReplace method");
        exit(0);
    }

    for (i = 0; i < a_rows; i++)
    {
        for (j = 0; j < b_columns; j++)
        {
            // dot product:
            for (k = 0; k < a_columns; k++)
            {
                sum += a_matrixVals[i][k] * b_matrixVals[k][j];
            }
            destination[i][j] = sum;
            sum = 0;
        }
    }
}

/**
 * @brief preforms matrix multplication but supports the result being able to
 * be stored in one of the multpliers.
 *
 * @param a_rows
 * @param a_columns
 * @param a_matrixVals
 * @param b_rows
 * @param b_columns
 * @param b_matrixVals
 * @param destination
 */
void laa_multiplyMMReplace(float **a_matrixVals, int a_rows, int a_columns, float **b_matrixVals, int b_rows, int b_columns, float **destination)
{
    int i, j, k;
    float sum = 0;
    float buffer[a_columns];
    if (a_columns != b_rows)
    {
        printf("matrices must have compatable dimensions!\n");
        exit(0);
    }

    for (i = 0; i < a_rows; i++)
    {
        sum = 0;
        for (k = 0; k < a_columns; k++)
        {
            buffer[k] = a_matrixVals[i][k];
            sum += a_matrixVals[i][k] * b_matrixVals[k][0];
        }
        destination[i][0] = sum;
        for (j = 1; j < b_columns; j++)
        {
            sum = 0;
            for (k = 0; k < a_columns; k++)
            {
                sum += buffer[k] * b_matrixVals[k][j];
            }
            destination[i][j] = sum;
        }
    }
}

float **laa_multiplyMM(float **a_matrixVals, int a_rows, int a_columns, float **b_matrixVals, int b_rows, int b_columns)
{
    int i, j, k;
    float sum = 0;
    if (a_columns != b_rows)
    {
        printf("matrices must have compatable dimensions!\n");
        exit(1);
    }

    float **matrixEntries = (float **)malloc(a_rows * sizeof(float *));
    if (!matrixEntries)
        exit(1);

    for (i = 0; i < a_rows; i++)
    {
        matrixEntries[i] = (float *)malloc(b_columns * sizeof(float));
        if (!matrixEntries[i])
        {
            for (; i > 0; i--)
                free(matrixEntries[i - 1]);
            free(matrixEntries);
            exit(1);
        }
        for (j = 0; j < b_columns; j++)
        {
            // dot product:
            for (k = 0; k < a_columns; k++)
            {
                sum += a_matrixVals[i][k] * b_matrixVals[k][j];
            }
            matrixEntries[i][j] = sum;
            sum = 0;
        }
    }
    return matrixEntries;
}

/**
 * @brief inverses the static matrix
 */
void laa_inverse(float **matrixEntries, int rows, int columns)
{
    static float invMatrix[LA_MAX_INVERSE_SIZE][2 * LA_MAX_INVERSE_SIZE];
    int i, j, k;
    if (rows != columns)
    {
        printf("No inverse of non-square matrix!");
        exit(0);
    }
    if (rows > LA_MAX_INVERSE_SIZE)
    {
        printf("Matrix too large to calculate inverse. Increase the LA_MAX_INVERSE_SIZE");
        exit(0);
    }
    float ratio;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            invMatrix[i][j] = matrixEntries[i][j];
            if (i == j)
            {
                invMatrix[i][j + columns] = 1;
            }
            else
            {
                invMatrix[i][j + columns] = 0;
            }
        }
    }

    /* Applying Gauss Jordan Elimination */
    for (i = 0; i < rows; i++)
    {
        if (invMatrix[i][i] == 0.0)
        {
            printf("Mathematical error while preforming inverse!");
            exit(0);
        }

        for (j = 0; j < columns; j++)
        {
            if (i != j)
            {
                ratio = invMatrix[j][i] / invMatrix[i][i];
                for (k = 0; k < 2 * columns; k++)
                {
                    invMatrix[j][k] = invMatrix[j][k] - ratio * invMatrix[i][k];
                }
            }
        }
    }

    /* Row Operation to Make Principal Diagonal to 1 */
    for (i = 0; i < rows; i++)
    {
        for (j = columns; j < 2 * columns; j++)
        {
            matrixEntries[i][j - columns] = invMatrix[i][j] / invMatrix[i][i];
        }
    }
}

float **laa_getInverse(float **matrixEntries, int rows, int columns)
{
    static float invMatrix[LA_MAX_INVERSE_SIZE][2 * LA_MAX_INVERSE_SIZE];
    int i, j, k;
    if (rows != columns)
    {
        printf("No inverse of non-square matrix!");
        exit(0);
    }
    if (rows > LA_MAX_INVERSE_SIZE)
    {
        printf("Matrix too large to calculate inverse. Increase the LA_MAX_INVERSE_SIZE");
        exit(0);
    }
    float ratio;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            invMatrix[i][j] = matrixEntries[i][j];
            if (i == j)
            {
                invMatrix[i][j + columns] = 1;
            }
            else
            {
                invMatrix[i][j + columns] = 0;
            }
        }
    }

    /* Applying Gauss Jordan Elimination */
    for (i = 0; i < rows; i++)
    {
        if (invMatrix[i][i] == 0.0)
        {
            printf("Mathematical error while preforming inverse!");
            exit(0);
        }

        for (j = 0; j < columns; j++)
        {
            if (i != j)
            {
                ratio = invMatrix[j][i] / invMatrix[i][i];
                for (k = 0; k < 2 * columns; k++)
                {
                    invMatrix[j][k] = invMatrix[j][k] - ratio * invMatrix[i][k];
                }
            }
        }
    }

    float **newMatrix = (float **)malloc(rows * sizeof(float *));
    if (!newMatrix)
        exit(1);

    /* Row Operation to Make Principal Diagonal to 1 */
    for (i = 0; i < rows; i++)
    {
        newMatrix[i] = (float *)malloc(columns * sizeof(float));
        if (!newMatrix[i])
        {
            for (; i > 0; i--)
                free(newMatrix[i - 1]);
            free(newMatrix);
            exit(1);
        }
        for (j = columns; j < 2 * columns; j++)
        {
            newMatrix[i][j - columns] = invMatrix[i][j] / invMatrix[i][i];
        }
    }
    return newMatrix;
}

#pragma endregion