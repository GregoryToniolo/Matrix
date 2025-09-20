# Matrix

A simple library written in C/C++ for basic Matrix and vector operations
## License

[MIT](https://choosealicense.com/licenses/mit/)


## Features
- memory safe functionality
- matrix constructor and destructor functions
- dot product, scalar product, addition, subtraction, hadamard product, hadamard division
- Mat2, Mat3, Mat4, Vec2, Vec3 and Vec4 types for more efficient operations
- M x N matrix type
- fast square root using heron's algorithm
- sine, cosine, tangent, secant, cosecant and cotangent functions with versions for use with degrees and radians
- vector homogenization and dehomogenization
- orthogonal and perspective projection
- 2d, 3d and 4d rotation

## Extra
- to use only vec2-4 and mat2-4 types define MATRIX_EMBEDDED before including the header file in your C/C++ program.

## Installation

Linux

```bash
  git clone https://github.com/GregoryToniolo/Matrix && cd Matrix && sudo cp g-matrix.h /usr/include
```
Usage:

Include the "g-matrix.h" file in your C/C++ program.

## FAQ:
What about the determinant, Matrix inversion, Matrix diagonalization ...etc...etc?

Not adding those because they are too hard to implement and understand.
