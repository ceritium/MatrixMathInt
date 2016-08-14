# MatrixMathInt

Fork of https://github.com/codebendercc/MatrixMath

It is the same but for work with integers.

## FEATURES
Briefly, the functions provided by MatrixMath:
```
void MatrixPrint(int* A, int m, int n, String label);
void MatrixCopy(int* A, int n, int m, int* B);
void MatrixMult(int* A, int* B, int m, int p, int n, int* C);
void MatrixAdd(int* A, int* B, int m, int n, int* C);
void MatrixSubtract(int* A, int* B, int m, int n, int* C);
void MatrixTranspose(int* A, int m, int n, int* C);
int MatrixInvert(int* A, int n);
```
Matrices should be stored in row-major arrays, which is fairly standard. The user must keep track of array dimensions and send them to the functions; mistakes on dimensions will not be caught by the library.

It's worth pointing out that the MatrixInvert() function uses Gauss-Jordan elimination with partial pivoting. Partial pivoting is a compromise between a numerically unstable algorithm and full pivoting, which involves more searching and swapping matrix elements.

Also, the inversion algorithm stores the result matrix on top of the the input matrix, meaning no extra memory is allocated during inversion but your original matrix is gone.

## HOW TO IMPORT/INSTALL
Grab the source code below, and put in a folder called MatrixMathInt. [Apparently recent Playground changes have killed the zip file. Please edit this if there is a preferred way of making the source available.

~~Put the MatrixMathInt folder in "libraries\".~~

*Edit: clone this repo at https://github.com/buckbaskin/MatrixMath.git into your `libraries\` folder*

In the Arduino IDE, create a new sketch (or open one) and

select from the menubar "Sketch->Import Library->MatrixMathInt".

Once the library is imported, a "#include MatrixMathInt.h" line will appear at the top of your Sketch.

The MatrixMathDemo in the Examples folder demonstrates multiplication and inversion using the MatrixPrint() function to show results.
