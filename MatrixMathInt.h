/*
 *  MatrixMath.h Library for Matrix Math
 *
 *  Created by Charlie Matlack on 12/18/10.
 *  Modified from code by RobH45345 on Arduino Forums, algorithm from 
 *  NUMERICAL RECIPES: The Art of Scientific Computing.
 *
 *  Modified by buckbaskin starting on 2/7/16.
 *  For original, see http://playground.arduino.cc/Code/MatrixMath
 */

#ifndef MatrixMath_h
#define MatrixMath_h

#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif

class MatrixMathInt
{
public:
    //MatrixMath();
    void Print(int* A, int m, int n, String label);
    void Copy(int* A, int n, int m, int* out);
    void Multiply(int* A, int* B, int m, int p, int n, int* out);
    void Add(int* A, int* B, int m, int n, int* out);
    void Subtract(int* A, int* B, int m, int n, int* out);
    void Transpose(int* A, int m, int n, int* out);
    void Scale(int* A, int m, int n, float k, int* out);
    int Invert(int* A, int n);
};

extern MatrixMathInt Matrix;
#endif
