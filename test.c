#include <g-vectors.h>
#include <stdio.h>

int main()
{
    vector a = { (double[]) {1, 2, 3}, 3};
    vector b = { (double[]) {1, 2, 3}, 3};

    double product = dot_product(&a, &b);

    a = vector_subtract(&a, &b);

    for (int i = 0; i < a.size; i++)
    {
        printf("%f ", a.components[i]);
    }

    printf("%f\n", product);

    return 0;
}
