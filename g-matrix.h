#ifndef __cplusplus
#include <stdlib.h>
#endif

#ifndef NULL
#define NULL nullptr
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
typedef struct {
    float* entries;
    unsigned int rows;
    unsigned int columns;
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
#ifndef __cplusplus
    matrix->entries = (float*) malloc(sizeof(float[matrix->rows * matrix->columns]));
#else
    matrix->entries = new float[matrix->rows * matrix->columns];
#endif

    if (matrix->entries == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix->rows * matrix->columns; i++) 
    {
        matrix->entries[i] = 0;
    }
}

void Matrix_free(Matrix* matrix)
{
#ifndef  __cplusplus
    free(matrix->entries);
#else
    delete[] matrix->entries;
#endif
}

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
    else if (matrix1->entries == NULL || matrix2->entries == NULL || result->entries == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->entries[j + i * result->columns] = matrix1->entries[j + i * result->columns] + matrix2->entries[j + i * result->columns];
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
    else if (matrix1->entries == NULL || matrix2->entries == NULL || result->entries == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->entries[j + i * result->columns] = matrix1->entries[j + i * result->columns] - matrix2->entries[j + i * result->columns];
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
    else if (matrix1->entries == NULL || matrix2->entries == NULL || result->entries == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->entries[j + i * result->columns] = matrix1->entries[j + i * result->columns] * matrix2->entries[j + i * result->columns];
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
    else if (matrix1->entries == NULL || matrix2->entries == NULL || result->entries == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix1->rows; i++) 
    {
        for (int j = 0; j < matrix2->columns; j++) 
        {
            result->entries[j + i * result->columns] = matrix1->entries[j + i * result->columns] / matrix2->entries[j + i * result->columns];
        }
    }
}

//initialize all matrices with matrix_init before calling this function
void Matrix_scalar(Matrix* matrix, float scalar, Matrix* result)
{
    if (!matrix || !result)
    {
        return;
    }
    else if ((matrix->columns != result->columns) || (matrix->rows != result->columns))
    {
        return;
    }
    else if (matrix->entries == NULL || result->entries == NULL)
    {
        return;
    }

    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < result->columns; j++) 
        {
            result->entries[j + i * result->columns] = matrix->entries[j + i * matrix->columns] * scalar;
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
    else if (matrix1->entries == NULL || matrix2->entries == NULL || result->entries == NULL)
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
                result->entries[j + i * matrix2->columns] += matrix1->entries[k + i * matrix2->rows] * matrix2->entries[j + k * matrix2->columns];
            }
        }
    }

}
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

Vec2 Vec3_dehomogenize(Vec3 a)
{
    return (Vec2) {a.x / a.z, a.y / a.z};
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

Vec3 Vec4_dehomogenize(Vec4 a)
{
    return (Vec3) {a.x / a.w, a.y / a.w, a.z / a.w};
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
    Mat3 c = {{0}};
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
    Mat4 c = {{0}};
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

Vec4 Vec4_transform(Mat4 mat, Vec4 vec)
{
    return (Vec4) {
        vec.x * mat.entries[0][0] + vec.y * mat.entries[0][1] + vec.z * mat.entries[0][2] + vec.w * mat.entries[0][3],
        vec.x * mat.entries[1][0] + vec.y * mat.entries[1][1] + vec.z * mat.entries[1][2] + vec.w * mat.entries[1][3],
        vec.x * mat.entries[2][0] + vec.y * mat.entries[2][1] + vec.z * mat.entries[2][2] + vec.w * mat.entries[2][3],
        vec.x * mat.entries[3][0] + vec.y * mat.entries[3][1] + vec.z * mat.entries[3][2] + vec.w * mat.entries[3][3],
    };
}
