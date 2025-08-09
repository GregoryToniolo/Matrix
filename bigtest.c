#include "g-matrix.h"
#include <stdio.h>

void printMatrix(Matrix* a)
{
    for (int i = 0; i < a->rows; i++) 
    {
        for (int j = 0; j < a->columns; j++) 
        {
            printf("%f ", a->elements[j + i * a->columns]);
        }
        printf("\n");
    }
}

int main(void)
{
    Matrix a;
    a.rows = 4;
    a.columns = 1;
    Matrix b;
    b.rows = 1;
    b.columns = 4;
    Matrix c;
    c.rows = 1;
    c.columns = 1;


    Matrix_init(&a);
    Matrix_init(&b);
    Matrix_init(&c);

    for (int i = 0; i < 4; i++) 
    {
        a.elements[i] = i;
        b.elements[i] = i;
    }

    Matrix_dot(&b, &a, &c);


    printMatrix(&c);

    Matrix_free(&a);
    Matrix_free(&b);
    Matrix_free(&c);
}
