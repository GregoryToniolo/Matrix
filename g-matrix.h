#include <stdint.h>
#include <stdlib.h>

typedef struct {
    double* elements;
    uint32_t rows;
    uint32_t columns;
}Matrix;

void Matrix_init(Matrix* matrix)
{
    if (!matrix)
    {
        return;
    }
    else if (matrix->rows == 0 || matrix->columns == 0)
    {
        return;
    }

    matrix->elements = (double*) malloc(sizeof(double[matrix->rows * matrix->columns]));

    if (matrix->elements == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix->rows * matrix->columns; i++) 
    {
        matrix->elements[i] = 0;
    }
}

void Matrix_free(Matrix* matrix)
{
    free(matrix->elements);
}

void Matrix_determinant(Matrix* matrix, long double* result); //algorithms for determinant are stupid difficult and mostly run on O(n^3) so probably won't implement it 
void Matrix_inverse(Matrix* matrix, Matrix* result); //uses determinant to check whether the matrix can be inverted

//initialize all matrices with matrix_init before calling this function
void Matrix_add(Matrix* matrix1, Matrix* matrix2, Matrix* result)
{
    if (!matrix1 || !matrix2 || !result)
    {
        return;
    }
    else if ((matrix1->columns != matrix2->columns) || (matrix2->columns != result->columns))
    {
        return;
    }
    else if ((matrix1->rows != matrix2->rows) || (matrix2->rows != result->rows)) {
        return;
    }
    else if (matrix1->elements == NULL || matrix2->elements == NULL || result->elements == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->elements[j + i * result->columns] = matrix1->elements[j + i * result->columns] + matrix2->elements[j + i * result->columns];
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_subtract(Matrix* matrix1, Matrix* matrix2, Matrix* result)
{
    if (!matrix1 || !matrix2 || !result)
    {
        return;
    }
    else if ((matrix1->columns != matrix2->columns) || (matrix2->columns != result->columns))
    {
        return;
    }
    else if ((matrix1->rows != matrix2->rows) || (matrix2->rows != result->rows)) {
        return;
    }
    else if (matrix1->elements == NULL || matrix2->elements == NULL || result->elements == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->elements[j + i * result->columns] = matrix1->elements[j + i * result->columns] - matrix2->elements[j + i * result->columns];
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_hadamard_product(Matrix* matrix1, Matrix* matrix2, Matrix* result)
{
    if (!matrix1 || !matrix2 || !result)
    {
        return;
    }
    else if ((matrix1->columns != matrix2->columns) || (matrix2->columns != result->columns))
    {
        return;
    }
    else if ((matrix1->rows != matrix2->rows) || (matrix2->rows != result->rows)) {
        return;
    }
    else if (matrix1->elements == NULL || matrix2->elements == NULL || result->elements == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->elements[j + i * result->columns] = matrix1->elements[j + i * result->columns] * matrix2->elements[j + i * result->columns];
        }
    }
}

void Matrix_hadamard_division(Matrix* matrix1, Matrix* matrix2, Matrix* result)
{
    if (!matrix1 || !matrix2 || !result)
    {
        return;
    }
    else if ((matrix1->columns != matrix2->columns) || (matrix2->columns != result->columns))
    {
        return;
    }
    else if ((matrix1->rows != matrix2->rows) || (matrix2->rows != result->rows)) {
        return;
    }
    else if (matrix1->elements == NULL || matrix2->elements == NULL || result->elements == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->elements[j + i * result->columns] = matrix1->elements[j + i * result->columns] / matrix2->elements[j + i * result->columns];
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_scalar(Matrix* matrix, double scalar, Matrix* result)
{
    if (!matrix || !result)
    {
        return;
    }
    else if ((matrix->columns != result->columns) || (matrix->rows != result->columns))
    {
        return;
    }
    else if (matrix->elements == NULL || result->elements == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < result->columns; j++) 
        {
            result->elements[j + i * result->columns] = matrix->elements[j + i * matrix->columns] * scalar;
        }
    }

}

//initialize all matrices with matrix_init before calling this function
void Matrix_dot(Matrix* matrix1, Matrix* matrix2, Matrix* result)
{
    if (!matrix1 || !matrix2 || !result)
    {
        return;
    }
    else if (matrix1->elements == NULL || matrix2->elements == NULL || result->elements == NULL)
    {
        return;
    }
    else if ((matrix1->columns != matrix2->rows) || (matrix1->rows != result->rows) || (matrix2->columns != result->columns))
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++)
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            for (int k = 0; k < matrix2->rows; k++)
            {
                result->elements[j + i * matrix2->columns] += matrix1->elements[k + i * matrix2->rows] * matrix2->elements[j + k * matrix2->columns];
            }
        }
    }

}
