#include "g-matrix.h"
#include <stdio.h>

int main(void)
{
    for (int theta = 0; theta <= 360; theta++)
    {
        Mat2 a = {{
        {cosine_deg(theta), -sine_deg(theta)},
        {sine_deg(theta), cosine_deg(theta)}
        }};

        Vec2 b = {1, 0};
        printf("%f %f\n", (Vec2_transform(a, b)).x, (Vec2_transform(a, b)).y);
    }
}
