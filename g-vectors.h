typedef struct
{
    double *components;
    int size;
}vector;

//performs dot product (weighted sum) between 2 vectors
double dot_product(const vector* Vec1,const vector* Vec2)
{
    double scalar = 0;

    for (int i = 0; i < Vec1->size; i++)
    {
        scalar += Vec1->components[i] * Vec2->components[i];
    }

    return scalar;
}

//performs the product between a vector and a scalar
vector scalar_product(vector* vector1,const double* scalar)
{
    for (int i = 0; i < vector1->size; i++)
    {
        vector1->components[i] *= *scalar;
    }

    vector result = *vector1;

    return result;
}

//performs addition between 2 vectors
vector vector_add(const vector* Vec1, const vector* Vec2)
{
    vector result = *Vec1;

    for (int i = 0; i < Vec1->size; i++)
    {
        result.components[i] = Vec1->components[i] + Vec2->components[i];
    }

    return result;
}

//performs subtraction between 2 vectors
vector vector_subtract(const vector* Vec1, const vector* Vec2)
{
    vector result = *Vec1;

    for (int i = 0; i < Vec1->size; i++)
    {
        result.components[i] = Vec1->components[i] - Vec2->components[i];
    }

    return result;
}
