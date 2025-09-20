#include "g-matrix.h"
#include <iostream>
#define PRINTMATRIX(c)\
    for (int i = 0; i < c.rows; i++)\
    {\
        for (int j = 0; j < c.columns; j++)\
        {\
            std::cout << c.entries[i][j] << ' ';\
        }\
        std::cout << '\n';\
    }\

int main()
{
    Matrix a(4, 4);
    Matrix b(4, 4);

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            a.entries[i][j] = i;
        }
    }

    for (int i = 0; i < b.rows; i++)
    {
        for (int j = 0; j < b.columns; j++)
        {
            b.entries[i][j] = j;
        }
    }

    Matrix c = Matrix_dot(a, b);

    std::cout << "here\n";

    PRINTMATRIX(a)
    PRINTMATRIX(b)
    PRINTMATRIX(c)
    
    return 0;
}
