#ifndef FMATH_H
#define FMATH_H

#include "fmatrix.h"
#include <cmath>

using namespace std;

// Math with tmatrix - Use with CAUTION: only with very small arrays, otherwise it will run slowly!

// Vector functions:
template<typename T>
inline T norm(const tmatrix<T> & a) {
	T n = 0;
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) n += pow(a.val[i], 2); // Only works for 1D (vector) fmatrix!
	return sqrt(n);
}

template<typename T>
inline T dot(const tmatrix<T> & a, const tmatrix<T> & b) {
	T n = 0;
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) n += a.val[i] * b.val[i]; // Only works for 1D (vector) fmatrix!
	return n;
}

template<typename T>
inline tmatrix<T> cross (const tmatrix<T> & a, const tmatrix<T> & b) { // Only works with 3D vectors (n1 = 3)!! 
	
	tmatrix<T> c(3);
	c.val[0] = a.val[1] * b.val[2] - a.val[2] * b.val[1];
	c.val[1] = a.val[2] * b.val[0] - a.val[0] * b.val[2];
	c.val[2] = a.val[0] * b.val[1] - a.val[1] * b.val[0];

	return c;
}


// SLOW (binary) operators:

template<typename T>
inline tmatrix<T> operator + (const tmatrix<T> & a, const tmatrix<T> & b) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] + b.val[i];
	return c;
}

template<typename T>
inline tmatrix<T> operator - (const tmatrix<T> & a, const tmatrix<T> & b) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] - b.val[i];
	return c;
}

template<typename T>
inline tmatrix<T> operator / (const tmatrix<T> & a, const tmatrix<T> & b) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] / b.val[i];
	return c;
}

template<typename T>
inline tmatrix<T> operator * (const tmatrix<T> & a, const tmatrix<T> & b)  {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] / b.val[i];
	return c;
}

template<typename T>
inline tmatrix<T> operator * (const tmatrix<T> & a, const T & b) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] * b;
	return c;
}

template<typename T>
inline tmatrix<T> operator * (const T & b, const tmatrix<T> & a) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] * b;
	return c;
}

template<typename T>
inline tmatrix<T> operator / (const tmatrix<T> & a, const T & b) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] / b;
	return c;
}

template<typename T>
inline tmatrix<T> operator / ( const  T & b, const  tmatrix<T> & a) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] / b;
	return c;
}

template<typename T>
inline tmatrix<T> operator + (const tmatrix<T> & a, const  T & b) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] + b;
	return c;
}

template<typename T>
inline tmatrix<T> operator + (const T & b, const tmatrix<T> & a) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] + b;
	return c;
}

template<typename T>
inline tmatrix<T> operator - (const tmatrix<T> & a, const T & b) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = a.val[i] - b;
	return c;
}

template<typename T>
inline tmatrix<T> operator - (const T & b, const tmatrix<T> & a) {
	tmatrix<T> c(a.n1, a.n2, a.n3);
	for (size_t i = 0; i < a.n1 * a.n2 * a.n3; i++) c.val[i] = b - a.val[i];
	return c;
}

// Trigonometric functions:
// Do we need it?



#endif // !FMATH_H

