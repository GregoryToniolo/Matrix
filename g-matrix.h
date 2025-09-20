#pragma once

#ifndef __cplusplus
#include <stdlib.h>
#include <stdio.h>
#else
#include <stdexcept>
#endif

#define PI 3.1415265359
#define TWOPI 6.28318530718
#define PI_OVER_2 1.57079632679
#define PI_SQUARED 9.86960440109

float fsqrt(float n)
{
    float x = (n > 100) ? n/2 : n / 32;

    for (int i = 0; i < 18; i++)
    {
        x = 0.5f * (x + n/x);
    }

    return x;
}

float sine_deg(float theta)
{
    while (theta >= 360)
    {
        theta -= 360;
    }

    if (theta <= 180)
    {
        return (4 * theta * (180 - theta)) / (40500 - theta * (180 - theta));
    }

    theta -= 180;
    return -1 * (4 * theta * (180 - theta)) / (40500 - theta * (180 - theta));
}

float sine_rad(float theta)
{
    while (theta >= TWOPI)
    {
        theta -= TWOPI;
    }

    if (theta <= PI)
    {
        return (16 * theta * (PI - theta)) / (5 * PI_SQUARED - 4 * theta * (PI - theta));
    }

    theta -= PI_OVER_2;
    return -1 * (16 * theta * (PI - theta)) / (5 * PI_SQUARED - 4 * theta * (PI - theta));
}

float sine_pi_over_n(int theta)
{
    if (theta < 0)
    {
        return  -1 * (16 * (theta - 1)) / (5 * theta * theta - 4 * theta + 4);
    }

    return (16 * (theta - 1)) / (5 * theta * theta - 4 * theta + 4);
}

float cosine_deg(float theta)
{
    return sine_deg(theta + 90);
}

float cosine_rad(float theta)
{
    return sine_rad(theta + PI_OVER_2);
}

float tangent_deg(float theta)
{
    return sine_deg(theta) / cosine_deg(theta);
}

float tangent_rad(float theta)
{
    return sine_rad(theta) / cosine_rad(theta);
}

float secant_deg(float theta)
{
    return 1 / cosine_deg(theta);
}

float secant_rad(float theta)
{
    return 1 / cosine_rad(theta);
}

float cotangent_deg(float theta)
{
    return cosine_deg(theta) / sine_deg(theta);
}

float cotangent_rad(float theta)
{
    return cosine_rad(theta) / sine_rad(theta);
}

float cosecant_deg(float theta)
{
    return 1 / sine_deg(theta);
}

float cosecant_rad(float theta)
{
    return 1 / sine_rad(theta);
}

#ifndef MATRIX_EMBEDDED
#ifdef __cplusplus

struct Matrix
{
    const unsigned int rows, columns;
    float** entries;
    Matrix(unsigned int lines, unsigned int cols) : rows(lines), columns(cols)
    {
        if (lines == 0 || cols == 0)
        {
            throw std::invalid_argument("matrix rows and columns must be nonzero");
        }

        entries = new float*[rows];

        for (int i = 0; i < rows; i++)
        {
            entries[i] = new float[columns]();
        }
    }

    //copy constructor
    Matrix(const Matrix& other) : rows(other.rows), columns(other.columns)
    {
        entries = new float*[rows];

        for (int i = 0; i < rows; i++)
        {
            entries[i] = new float[columns]();
        }

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                entries[i][j] = other.entries[i][j];
            }
        }
    }

    float Get(const unsigned int row, const unsigned int column) const
    {
        if (row < rows && row >= 0 && column < columns && column >= 0)
        {
            return entries[row][column];
        }
        throw std::out_of_range("tried to access element not in matrix");
    }

    ~Matrix()
    {
        for (int i = 0; i < rows; i++)
        {
            delete[] entries[i];
            entries[i] = nullptr;
        }
        delete[] entries;
        entries = nullptr;
    }
};

Matrix Matrix_add(const Matrix& a, const Matrix& b)
{
    if (a.rows != b.rows || a.columns != b.columns)
    {
        throw std::invalid_argument("invalid number of rows or columns in input matrix");
    }
    Matrix result(a.rows, b.columns);

    for (int i = 0; i < result.rows; i++)
    {
        for (int j = 0; j < result.columns; j++)
        {
            result.entries[i][j] = a.entries[i][j] + b.entries[i][j];
        }
    }

    return result;
}

Matrix Matrix_subtract(const Matrix& a, const Matrix& b)
{
    if (a.rows != b.rows || a.columns != b.columns)
    {
        throw std::invalid_argument("invalid number of rows or columns in input matrix");
    }
    Matrix result(a.rows, b.columns);

    for (int i = 0; i < result.rows; i++)
    {
        for (int j = 0; j < result.columns; j++)
        {
            result.entries[i][j] = a.entries[i][j] - b.entries[i][j];
        }
    }

    return result;
}

Matrix Matrix_hadamard_product(const Matrix& a, const Matrix& b)
{
    if (a.rows != b.rows || a.columns != b.columns)
    {
        throw std::invalid_argument("invalid number of rows or columns in input matrix");
    }
    Matrix result(a.rows, b.columns);

    for (int i = 0; i < result.rows; i++)
    {
        for (int j = 0; j < result.columns; j++)
        {
            result.entries[i][j] = a.entries[i][j] * b.entries[i][j];
        }
    }

    return result;
}

Matrix Matrix_hadamard_division(const Matrix& a, const Matrix& b)
{
    if (a.rows != b.rows || a.columns != b.columns)
    {
        throw std::invalid_argument("invalid number of rows or columns in input matrix");
    }
    Matrix result(a.rows, b.columns);

    for (int i = 0; i < result.rows; i++)
    {
        for (int j = 0; j < result.columns; j++)
        {
            result.entries[i][j] = a.entries[i][j] / b.entries[i][j];
        }
    }

    return result;
}

Matrix Matrix_scalar(const Matrix& a, const float b)
{
    Matrix result(a.rows, a.columns);

    for (int i = 0; i < result.rows; i++)
    {
        for (int j = 0; j < result.columns; j++)
        {
            result.entries[i][j] = a.entries[i][j] * b;
        }
    }

    return result;
}

Matrix Matrix_dot(const Matrix& a, const Matrix& b)
{
    if (a.columns != b.rows)
    {
        throw std::invalid_argument("invalid number of rows or columns in input matrix");
    }

    Matrix result(a.rows, b.columns);

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < b.columns; j++)
        {
            for (int k = 0; k < b.rows; k++)
            {
                result.entries[i][j] += a.entries[i][k] * b.entries[k][j];
            }
        }
    }

    return result;
}

#else

void terminate(const char* error)
{
    fprintf(stderr, "Error [%s] has occurred now aborting.\n", error);
    exit(EXIT_FAILURE);
}

typedef struct {
    float** entries;
    unsigned int rows;
    unsigned int columns;
}Matrix;

void Matrix_init(Matrix* matrix)
{
    if (matrix == NULL)
    {
        terminate("cannot initialize null pointer");
    }
    else if (matrix->rows == 0 || matrix->columns == 0)
    {
        terminate("matrix rows and columns must be nonzero");
    }
    
    matrix->entries = (float**) malloc(sizeof(float*[matrix->rows]));

    for (int i = 0; i < matrix->rows; i++)
    {
        matrix->entries[i] = (float*) malloc(sizeof(float[matrix->columns]));
    }

    if (matrix->entries == NULL)
    {
        terminate("matrix initialization has failed");
        return;
    }

    for (int i = 0; i < matrix->rows; i++) 
    {
        for (int j = 0; j < matrix->rows; j++) 
        {
            matrix->entries[i][j] = 0;
        }
    }
}

void Matrix_free(Matrix* matrix)
{
    for (int i = 0; i < matrix->rows; i++)
    {
        free(matrix->entries[i]);
    }
    free(matrix->entries);
}

//initialize all matrices with matrix_init before calling this function
void Matrix_add(Matrix* a, Matrix* b, Matrix* result)
{
    if (!a || !b || !result)
    {
        terminate("cannot perform operations on null pointer");
    }
    else if ((a->columns != b->columns) || (b->columns != result->columns))
    {
        terminate("invalid number of columns in input or result matrix");
    }
    else if ((a->rows != b->rows) || (b->rows != result->rows)) 
    {
        terminate("invalid number of rows in input or result matrix");
    }
    else if (a->entries == NULL || b->entries == NULL || result->entries == NULL)
    {
        terminate("cannot perform operations on matrix with null entries");
    }

    for (int i = 0; i < a->rows; i++) 
    {
        for (int j = 0; j < b->columns; j++) 
        {
            result->entries[i][j] = a->entries[i][j] + b->entries[i][j];
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_subtract(Matrix* a, Matrix* b, Matrix* result)
{
    if (!a || !b || !result)
    {
        terminate("cannot perform operations on null pointer");
    }
    else if ((a->columns != b->columns) || (b->columns != result->columns))
    {
        terminate("invalid number of columns in input or result matrix");
    }
    else if ((a->rows != b->rows) || (b->rows != result->rows)) 
    {
        terminate("invalid number of rows in input or result matrix");
    }
    else if (a->entries == NULL || b->entries == NULL || result->entries == NULL)
    {
        terminate("cannot perform operations on matrix with null entries");
    }

    for (int i = 0; i < a->rows; i++) 
    {
        for (int j = 0; j < b->columns; j++) 
        {
            result->entries[i][j] = a->entries[i][j] - b->entries[i][j];
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_hadamard_product(Matrix* a, Matrix* b, Matrix* result)
{
    if (!a || !b || !result)
    {
        terminate("cannot perform operations on null pointer");
    }
    else if ((a->columns != b->columns) || (b->columns != result->columns))
    {
        terminate("invalid number of columns in input or result matrix");
    }
    else if ((a->rows != b->rows) || (b->rows != result->rows)) 
    {
        terminate("invalid number of rows in input or result matrix");
    }
    else if (a->entries == NULL || b->entries == NULL || result->entries == NULL)
    {
        terminate("cannot perform operations on matrix with null entries");
    }

    for (int i = 0; i < a->rows; i++) 
    {
        for (int j = 0; j < b->columns; j++) 
        {
            result->entries[i][j] = a->entries[i][j] * b->entries[i][j];
        }
    }
}

void Matrix_hadamard_division(Matrix* a, Matrix* b, Matrix* result)
{
    if (!a || !b || !result)
    {
        terminate("cannot perform operations on null pointer");
    }
    else if ((a->columns != b->columns) || (b->columns != result->columns))
    {
        terminate("invalid number of columns in input or result matrix");
    }
    else if ((a->rows != b->rows) || (b->rows != result->rows)) 
    {
        terminate("invalid number of rows in input or result matrix");
    }
    else if (a->entries == NULL || b->entries == NULL || result->entries == NULL)
    {
        terminate("cannot perform operations on matrix with null entries");
    }

    for (int i = 0; i < a->rows; i++) 
    {
        for (int j = 0; j < b->columns; j++) 
        {
            result->entries[i][j] = a->entries[i][j] / b->entries[i][j];
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_scalar(Matrix* matrix, float scalar, Matrix* result)
{
    if (!matrix || !result)
    {
        terminate("cannot perform operations on null pointer");
    }
    else if ((matrix->columns != result->columns) || (matrix->rows != result->rows))
    {
        terminate("invalid number of rows or columns in input or result matrix");
    }
    else if (matrix->entries == NULL || result->entries == NULL)
    {
        terminate("cannot perform operations on matrix with null entries");
    }

    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < result->columns; j++) 
        {
            result->entries[i][j] = matrix->entries[i][j] * scalar;
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_dot(Matrix* a, Matrix* b, Matrix* result)
{
    if (!a || !b || !result)
    {
        terminate("cannot perform operations on null pointer");
    }
    else if (a->entries == NULL || b->entries == NULL || result->entries == NULL)
    {
        terminate("cannot perform operations on matrix with null entries");
    }
    else if ((a->columns != b->rows) || (a->rows != result->rows) || (b->columns != result->columns))
    {
        terminate("invalid number of rows or columns in input or result matrix");
    }

    for (int i = 0; i < a->rows; i++)
    {
        for (int j = 0; j < b->columns; j++) 
        {
            for (int k = 0; k < b->rows; k++)
            {
                result->entries[i][j] += a->entries[i][k] * b->entries[k][j];
            }
        }
    }

}
#endif

#endif

typedef struct
{
    float x, y;
}Vec2;

typedef struct
{
    float x, y, z;
}Vec3;

typedef struct
{
    float x, y, z, w;
}Vec4;

typedef struct
{
    float entries[2][2];
}Mat2;

typedef struct
{
    float entries[3][3];
}Mat3;

typedef struct
{
    float entries[4][4];
}Mat4;

Vec2 Vec2_add(Vec2 a, Vec2 b)
{
    return (Vec2) {a.x + b.x, a.y + b.y};
}

Vec2 Vec2_subtract(Vec2 a, Vec2 b)
{
    return (Vec2) {a.x - b.x, a.y - b.y};
}

Vec2 Vec2_scalar(Vec2 a, float b)
{
    return (Vec2) {a.x * b, a.y * b};
}

float Vec2_dot(Vec2 a, Vec2 b)
{
    return a.x * b.x + a.y * b.y;
}

float Vec2_length(Vec2 a)
{
    return fsqrt(a.x * a.x + a.y * a.y);
}

float Vec2_distance(Vec2 a, Vec2 b)
{
    return fsqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

Vec2 Vec3_dehomogenize(Vec3 a)
{
    return (Vec2) {a.x / a.z, a.y / a.z};
}

Vec2 Vec3_project_orthogonal(Vec3 a)
{
    return (Vec2) {a.x, a.y};
}

Vec2 Vec3_project_perspective(Vec3 a, float distance)
{
    return (Vec2) {a.x /(distance - a.z) , a.y/(distance - a.z)};
}

Vec3 Vec2_homogenize(Vec2 a)
{
    return (Vec3) {a.x, a.y, 1};
}

Vec3 Vec3_add(Vec3 a, Vec3 b)
{
    return (Vec3) {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 Vec3_subtract(Vec3 a, Vec3 b)
{
    return (Vec3) {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 Vec3_scalar(Vec3 a, float b)
{
    return (Vec3) {a.x * b, a.y * b, a.z * b};
}

float Vec3_dot(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

float Vec3_length(Vec3 a)
{
    return fsqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

float Vec3_distance(Vec3 a, Vec3 b)
{
    return fsqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

Vec3 Vec4_dehomogenize(Vec4 a)
{
    return (Vec3) {a.x / a.w, a.y / a.w, a.z / a.w};
}

Vec3 Vec4_project_orthogonal(Vec4 a)
{
    return (Vec3) {a.x, a.y, a.z};
}

Vec3 Vec4_project_perspective(Vec4 a, float distance)
{
    return (Vec3) {a.x / (distance - a.w), a.y / (distance - a.w), a.z / (distance - a.w)};
}

Vec4 Vec3_homogenize(Vec3 a)
{
    return (Vec4) {a.x, a.y, a.z, 1};
}

Vec4 Vec4_add(Vec4 a, Vec4 b)
{
    return (Vec4) {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

Vec4 Vec4_subtract(Vec4 a, Vec4 b)
{
    return (Vec4) {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

Vec4 Vec4_scalar(Vec4 a, float b)
{
    return (Vec4) {a.x * b, a.y * b, a.z * b, a.w * b};
}

float Vec4_dot(Vec4 a, Vec4 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

float Vec4_length(Vec4 a)
{
    return fsqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
}

float Vec4_distance(Vec4 a, Vec4 b)
{
    return fsqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z) + (a.w - b.w) * (a.w - b.w));
}

Mat2 Mat2_add(Mat2 a, Mat2 b)
{
    return (Mat2) {{
        {a.entries[0][0] + b.entries[0][0], a.entries[0][1] + b.entries[0][1]},
        {a.entries[1][0] + b.entries[1][0], a.entries[1][1] + b.entries[1][1]}
    }};
}

Mat2 Mat2_subtract(Mat2 a, Mat2 b)
{
    return (Mat2) {{
        {a.entries[0][0] - b.entries[0][0], a.entries[0][1] - b.entries[0][1]},
        {a.entries[1][0] - b.entries[1][0], a.entries[1][1] - b.entries[1][1]}
    }};
}

Mat2 Mat2_hadamard(Mat2 a, Mat2 b)
{
    return (Mat2) {{
        {a.entries[0][0] * b.entries[0][0], a.entries[0][1] * b.entries[0][1]},
        {a.entries[1][0] * b.entries[1][0], a.entries[1][1] * b.entries[1][1]}
    }};
}

Mat2 Mat2_scalar(Mat2 a, float b)
{
    return (Mat2) {{
        {a.entries[0][0] * b, a.entries[0][1] * b},
        {a.entries[1][0] * b, a.entries[1][1] * b}
    }};
}

Mat2 Mat2_dot(Mat2 a, Mat2 b)
{
    return (Mat2) {{
        {a.entries[0][0] * b.entries[0][0] + a.entries[0][1] * b.entries[1][0],a.entries[0][0] * b.entries[0][1] + a.entries[0][1] * b.entries[1][1]},
        {a.entries[1][0] * b.entries[0][0] + a.entries[1][1] * b.entries[1][0],a.entries[1][0] * b.entries[0][1] + a.entries[1][1] * b.entries[1][1]},
    }};
}

Mat2 Mat2_power(Mat2 a, unsigned int exp)
{
    Mat2 b = {{
        {1, 0},
        {0, 1}
    }};

    while (exp)
    {
        if (exp % 2 == 1)
        {
            b = Mat2_dot(b, a);
        }
        a =  Mat2_dot(a, a);
        exp >>= 1;
    }

    return b;
}

Vec2 Vec2_transform(Mat2 mat, Vec2 vec)
{
    return (Vec2) {mat.entries[0][0] * vec.x + mat.entries[0][1] * vec.y, mat.entries[1][0] * vec.x + mat.entries[1][1] * vec.y};
}

Mat3 Mat3_add(Mat3 a, Mat3 b)
{
    return (Mat3) {{
        {a.entries[0][0] + b.entries[0][0], a.entries[0][1] + b.entries[0][1], a.entries[0][2] + b.entries[0][2]},
        {a.entries[1][0] + b.entries[1][0], a.entries[1][1] + b.entries[1][1], a.entries[1][2] + b.entries[1][2]},
        {a.entries[2][0] + b.entries[2][0], a.entries[2][1] + b.entries[2][1], a.entries[2][2] + b.entries[2][2]}
    }};
}

Mat3 Mat3_subtract(Mat3 a, Mat3 b)
{
    return (Mat3) {{
        {a.entries[0][0] - b.entries[0][0], a.entries[0][1] - b.entries[0][1], a.entries[0][2] - b.entries[0][2]},
        {a.entries[1][0] - b.entries[1][0], a.entries[1][1] - b.entries[1][1], a.entries[1][2] - b.entries[1][2]},
        {a.entries[2][0] - b.entries[2][0], a.entries[2][1] - b.entries[2][1], a.entries[2][2] - b.entries[2][2]}
    }};
}

Mat3 Mat3_hadamard(Mat3 a, Mat3 b)
{
    return (Mat3) {{
        {a.entries[0][0] * b.entries[0][0], a.entries[0][1] * b.entries[0][1], a.entries[0][2] * b.entries[0][2]},
        {a.entries[1][0] * b.entries[1][0], a.entries[1][1] * b.entries[1][1], a.entries[1][2] * b.entries[1][2]},
        {a.entries[2][0] * b.entries[2][0], a.entries[2][1] * b.entries[2][1], a.entries[2][2] * b.entries[2][2]}
    }};
}

Mat3 Mat3_scalar(Mat3 a, float b)
{
    return (Mat3) {{
        {a.entries[0][0] * b, a.entries[0][1] * b, a.entries[0][2] * b},
        {a.entries[1][0] * b, a.entries[1][1] * b, a.entries[1][2] * b},
        {a.entries[2][0] * b, a.entries[2][1] * b, a.entries[2][2] * b}
    }};
}

Mat3 Mat3_dot(Mat3 a, Mat3 b)
{
    Mat3 c = {{{0}}};
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                c.entries[i][j] += a.entries[i][k] * b.entries[k][j];
            }
        }
    }
    return c;
}

Mat3 Mat3_power(Mat3 a, unsigned int exp)
{
    Mat3 b = {{
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    }};

    while (exp)
    {
        if (exp % 2 == 1)
        {
            b = Mat3_dot(b, a);
        }
        a =  Mat3_dot(a, a);
        exp >>= 1;
    }

    return b;
}

Vec3 Vec3_transform(Mat3 mat, Vec3 vec)
{
    return (Vec3) {
        vec.x * mat.entries[0][0] + vec.y * mat.entries[0][1] + vec.z * mat.entries[0][2],
        vec.x * mat.entries[1][0] + vec.y * mat.entries[1][1] + vec.z * mat.entries[1][2],
        vec.x * mat.entries[2][0] + vec.y * mat.entries[2][1] + vec.z * mat.entries[2][2]
    };
}

Mat4 Mat4_add(Mat4 a, Mat4 b)
{
    return (Mat4) {{
        {a.entries[0][0] + b.entries[0][0], a.entries[0][1] + b.entries[0][1], a.entries[0][2] + b.entries[0][2], a.entries[0][3] + b.entries[0][3]},
        {a.entries[1][0] + b.entries[1][0], a.entries[1][1] + b.entries[1][1], a.entries[1][2] + b.entries[1][2], a.entries[1][3] + b.entries[1][3]},
        {a.entries[2][0] + b.entries[2][0], a.entries[2][1] + b.entries[2][1], a.entries[2][2] + b.entries[2][2], a.entries[2][3] + b.entries[2][3]},
        {a.entries[3][0] + b.entries[3][0], a.entries[3][1] + b.entries[3][1], a.entries[3][2] + b.entries[3][2], a.entries[3][3] + b.entries[3][3]},
    }};
}

Mat4 Mat4_subtract(Mat4 a, Mat4 b)
{
    return (Mat4) {{
        {a.entries[0][0] - b.entries[0][0], a.entries[0][1] - b.entries[0][1], a.entries[0][2] - b.entries[0][2], a.entries[0][3] - b.entries[0][3]},
        {a.entries[1][0] - b.entries[1][0], a.entries[1][1] - b.entries[1][1], a.entries[1][2] - b.entries[1][2], a.entries[1][3] - b.entries[1][3]},
        {a.entries[2][0] - b.entries[2][0], a.entries[2][1] - b.entries[2][1], a.entries[2][2] - b.entries[2][2], a.entries[2][3] - b.entries[2][3]},
        {a.entries[3][0] - b.entries[3][0], a.entries[3][1] - b.entries[3][1], a.entries[3][2] - b.entries[3][2], a.entries[3][3] - b.entries[3][3]},
    }};
}

Mat4 Mat4_hadamard(Mat4 a, Mat4 b)
{
    return (Mat4) {{
        {a.entries[0][0] * b.entries[0][0], a.entries[0][1] * b.entries[0][1], a.entries[0][2] * b.entries[0][2], a.entries[0][3] * b.entries[0][3]},
        {a.entries[1][0] * b.entries[1][0], a.entries[1][1] * b.entries[1][1], a.entries[1][2] * b.entries[1][2], a.entries[1][3] * b.entries[1][3]},
        {a.entries[2][0] * b.entries[2][0], a.entries[2][1] * b.entries[2][1], a.entries[2][2] * b.entries[2][2], a.entries[2][3] * b.entries[2][3]},
        {a.entries[3][0] * b.entries[3][0], a.entries[3][1] * b.entries[3][1], a.entries[3][2] * b.entries[3][2], a.entries[3][3] * b.entries[3][3]},
    }};
}

Mat4 Mat4_dot(Mat4 a, Mat4 b)
{
    Mat4 c = {{{0}}};
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                c.entries[i][j] += a.entries[i][k] * b.entries[k][j];
            }
        }
    }
    return c;
}

Mat4 Mat4_power(Mat4 a, unsigned int exp)
{
    Mat4 b = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};

    while (exp)
    {
        if (exp % 2 == 1)
        {
            b = Mat4_dot(b, a);
        }
        a =  Mat4_dot(a, a);
        exp >>= 1;
    }

    return b;
}

Vec4 Vec4_transform(Mat4 mat, Vec4 vec)
{
    return (Vec4) {
        vec.x * mat.entries[0][0] + vec.y * mat.entries[0][1] + vec.z * mat.entries[0][2] + vec.w * mat.entries[0][3],
        vec.x * mat.entries[1][0] + vec.y * mat.entries[1][1] + vec.z * mat.entries[1][2] + vec.w * mat.entries[1][3],
        vec.x * mat.entries[2][0] + vec.y * mat.entries[2][1] + vec.z * mat.entries[2][2] + vec.w * mat.entries[2][3],
        vec.x * mat.entries[3][0] + vec.y * mat.entries[3][1] + vec.z * mat.entries[3][2] + vec.w * mat.entries[3][3],
    };
}

Vec2 Vec2_rotate_deg(Vec2 vec, float theta)
{
    Mat2 rotation = {{
        {cosine_deg(theta), -sine_deg(theta)},
        {sine_deg(theta), cosine_deg(theta)}
    }};
    return Vec2_transform(rotation, vec);
}

Vec2 Vec2_rotate_rad(Vec2 vec, float theta)
{
    Mat2 rotation = {{
        {cosine_rad(theta), -sine_rad(theta)},
        {sine_rad(theta), cosine_rad(theta)}
    }};
    return Vec2_transform(rotation, vec);
}

Vec3 Vec3_rotateX_deg(Vec3 vec, float theta)
{
    Mat3 rotation = {{
        {1, 0, 0},
        {0, cosine_deg(theta), -sine_deg(theta)},
        {0, sine_deg(theta), cosine_deg(theta)}
    }};

    return Vec3_transform(rotation, vec);
}

Vec3 Vec3_rotateX_rad(Vec3 vec, float theta)
{
    Mat3 rotation = {{
        {1, 0, 0},
        {0, cosine_rad(theta), -sine_rad(theta)},
        {0, sine_rad(theta), cosine_rad(theta)}
    }};

    return Vec3_transform(rotation, vec);
}

Vec3 Vec3_rotateY_deg(Vec3 vec, float theta)
{
    Mat3 rotation = {{
        {cosine_deg(theta), 0, sine_deg(theta)},
        {0, 1, 0},
        {-sine_deg(theta), 0, cosine_deg(theta)},
    }};

    return Vec3_transform(rotation, vec);
}

Vec3 Vec3_rotateY_rad(Vec3 vec, float theta)
{
    Mat3 rotation = {{
        {cosine_rad(theta), 0, sine_rad(theta)},
        {0, 1, 0},
        {-sine_rad(theta), 0, cosine_rad(theta)},
    }};

    return Vec3_transform(rotation, vec);
}

Vec3 Vec3_rotateZ_deg(Vec3 vec, float theta)
{
    Mat3 rotation = {{
        {cosine_deg(theta), -sine_deg(theta), 0},
        {sine_deg(theta), cosine_deg(theta), 0},
        {0, 0, 1}
    }};

    return Vec3_transform(rotation, vec);
}

Vec3 Vec3_rotateZ_rad(Vec3 vec, float theta)
{
    Mat3 rotation = {{
        {cosine_rad(theta), -sine_rad(theta), 0},
        {sine_rad(theta), cosine_rad(theta), 0},
        {0, 0, 1}
    }};

    return Vec3_transform(rotation, vec);
}

Vec4 Vec4_rotateZW_deg(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {cosine_deg(theta), -sine_deg(theta), 0, 0},
        {sine_deg(theta), cosine_deg(theta), 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateZW_rad(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {cosine_rad(theta), -sine_rad(theta), 0, 0},
        {sine_rad(theta), cosine_rad(theta), 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateYW_deg(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {cosine_deg(theta), 0, -sine_deg(theta), 0},
        {0, 1, 0, 0},
        {sine_deg(theta), 0, cosine_deg(theta), 0},
        {0, 0, 0, 1}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateYW_rad(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {cosine_rad(theta), 0, -sine_rad(theta), 0},
        {0, 1, 0, 0},
        {sine_rad(theta), 0, cosine_rad(theta), 0},
        {0, 0, 0, 1}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateYZ_deg(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {cosine_deg(theta), 0, 0, -sine_deg(theta)},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {sine_deg(theta), 0, 0, cosine_deg(theta)},
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateYZ_rad(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {cosine_rad(theta), 0, 0, -sine_rad(theta)},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {sine_rad(theta), 0, 0, cosine_rad(theta)},
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateXW_deg(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {1, 0, 0, 0},
        {0, cosine_deg(theta), -sine_deg(theta), 0},
        {0, sine_deg(theta), cosine_deg(theta), 0},
        {0, 0, 0, 1}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateXW_rad(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {1, 0, 0, 0},
        {0, cosine_rad(theta), -sine_rad(theta), 0},
        {0, sine_rad(theta), cosine_rad(theta), 0},
        {0, 0, 0, 1}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateXZ_deg(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {1, 0, 0, 0},
        {0, cosine_deg(theta), 0, -sine_deg(theta)},
        {0, 0, 1, 0},
        {0, sine_deg(theta), 0, cosine_deg(theta)}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateXZ_rad(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {1, 0, 0, 0},
        {0, cosine_rad(theta), 0, -sine_rad(theta)},
        {0, 0, 1, 0},
        {0, sine_rad(theta), 0, cosine_rad(theta)}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateXY_deg(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, cosine_deg(theta), -sine_deg(theta)},
        {0, 0, sine_deg(theta), cosine_deg(theta)}
    }};
    return Vec4_transform(rotation, vec);
}

Vec4 Vec4_rotateXY_rad(Vec4 vec, float theta)
{
    Mat4 rotation = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, cosine_rad(theta), -sine_rad(theta)},
        {0, 0, sine_rad(theta), cosine_rad(theta)}
    }};
    return Vec4_transform(rotation, vec);
}
