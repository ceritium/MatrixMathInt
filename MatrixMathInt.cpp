/*
 *  MatrixMath.cpp Library for Matrix Math
 *
 *  Created by Charlie Matlack on 12/18/10.
 *  Modified from code by RobH45345 on Arduino Forums, algorithm from 
 *  NUMERICAL RECIPES: The Art of Scientific Computing.
 *
 *  Modified by buckbaskin starting on 2/7/16.
 *  For original see: http://playground.arduino.cc/Code/MatrixMath
 */

#include "MatrixMathInt.h"

#define NR_END 1

MatrixMathInt Matrix;			// Pre-instantiate

// Matrix Printing Routine
// Uses tabs to separate numbers under assumption printed float width won't cause problems
void MatrixMathInt::Print(int* input, int m, int n, String label){
	// input = input matrix (m x n)
	int i,j;
	Serial.println();
	Serial.println(label);
	for (i=0; i<m; i++){
		for (j=0;j<n;j++){
			Serial.print(input[n*i+j]);
		}
		Serial.println();
	}
}

void MatrixMathInt::Copy(int* input, int n, int m, int* out)
{
	int i, j, k;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			out[n*i+j] = input[n*i+j];
		}
}

//Matrix Multiplication Routine
// out = left*right
void MatrixMathInt::Multiply(int* left, int* right, int m, int p, int n, int* out)
{
	// left = input matrix (m x p)
	// right = input matrix (p x n)
	// m = number of rows in left
	// p = number of columns in left = number of rows in right
	// n = number of columns in right
	// out = output matrix = left*right (m x n)
	int i, j, k;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			out[n*i+j]=0;
			for (k=0;k<p;k++)
				out[n*i+j]= out[n*i+j]+left[p*i+k]*right[n*k+j];
		}
}


//Matrix Addition Routine
void MatrixMathInt::Add(int* left, int* right, int m, int n, int* out)
{
	// left = input matrix (m x n)
	// right = input matrix (m x n)
	// m = number of rows in left = number of rows in right
	// n = number of columns in left = number of columns in right
	// out = output matrix = left+right (m x n)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			out[n*i+j]=left[n*i+j]+right[n*i+j];
}


//Matrix Subtraction Routine
void MatrixMathInt::Subtract(int* left, int* right, int m, int n, int* out)
{
	// left = input matrix (m x n)
	// right = input matrix (m x n)
	// m = number of rows in left = number of rows in right
	// n = number of columns in left = number of columns in right
	// out = output matrix = left-right (m x n)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			out[n*i+j]=left[n*i+j]-right[n*i+j];
}


//Matrix Transpose Routine
void MatrixMathInt::Transpose(int* input, int m, int n, int* out)
{
	// input = input matrix (m x n)
	// m = number of rows in input
	// n = number of columns in input
	// out = output matrix = the transpose of input (n x m)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			out[m*j+i]=input[n*i+j];
}

void MatrixMathInt::Scale(int* input, int m, int n, float k, int* out)
{
	// input = input matrix (m x n)
	// m = number of rows in input
	// n = number of columns in input
	// out = output matrix = elementwise multiplication of k*input[][]
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			out[n*i+j] = input[n*i+j]*k;
}


//Matrix Inversion Routine
// * This function inverts a matrix based on the Gauss Jordan method.
// * Specifically, it uses partial pivoting to improve numeric stability.
// * The algorithm is drawn from those presented in 
//	 NUMERICAL RECIPES: The Art of Scientific Computing.
// * The function returns 1 on success, 0 on failure.
// * NOTE: The argument is ALSO the result matrix, meaning the input matrix is REPLACED
int MatrixMathInt::Invert(int* square, int n)
{
	// square = input matrix AND result matrix
	// n = number of rows = number of columns in A (n x n)
	int pivrow;		// keeps track of current pivot row
	int k,i,j;		// k: overall index along diagonal; i: row index; j: col index
	int pivrows[n]; // keeps track of rows swaps to undo at end
	float tmp;		// used for finding max value and making column swaps

	for (k = 0; k < n; k++)
	{
		// find pivot row, the row with biggest entry in current column
		tmp = 0;
		for (i = k; i < n; i++)
		{
			if (abs(square[i*n+k]) >= tmp)	// 'Avoid using other functions inside abs()?'
			{
				tmp = abs(square[i*n+k]);
				pivrow = i;
			}
		}

		// check for singular matrix
		if (square[pivrow*n+k] == 0.0f)
		{
			Serial.println("Inversion failed due to singular matrix");
			return 0;
		}

		// Execute pivot (row swap) if needed
		if (pivrow != k)
		{
			// swap row k with pivrow
			for (j = 0; j < n; j++)
			{
				tmp = square[k*n+j];
				square[k*n+j] = square[pivrow*n+j];
				square[pivrow*n+j] = tmp;
			}
		}
		pivrows[k] = pivrow;	// record row swap (even if no swap happened)

		tmp = 1.0f/square[k*n+k];	// invert pivot element
		square[k*n+k] = 1.0f;		// This element of input matrix becomes result matrix

		// Perform row reduction (divide every element by pivot)
		for (j = 0; j < n; j++)
		{
			square[k*n+j] = square[k*n+j]*tmp;
		}

		// Now eliminate all other entries in this column
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				tmp = square[i*n+k];
				square[i*n+k] = 0.0f;  // The other place where in matrix becomes result mat
				for (j = 0; j < n; j++)
				{
					square[i*n+j] = square[i*n+j] - square[k*n+j]*tmp;
				}
			}
		}
	}

	// Done, now need to undo pivot row swaps by doing column swaps in reverse order
	for (k = n-1; k >= 0; k--)
	{
		if (pivrows[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				tmp = square[i*n+k];
				square[i*n+k] = square[i*n+pivrows[k]];
				square[i*n+pivrows[k]] = tmp;
			}
		}
	}
	return 1;
}
