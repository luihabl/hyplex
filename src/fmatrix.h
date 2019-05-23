#ifndef FMATRIX_H
#define FMATRIX_H

#include <iostream>
#include <initializer_list>
#include <cmath>
using namespace std;

// tmatrix stants for template matrix

template <class T>
struct tmatrix
{
	size_t n1, n2, n3;
	T * val;

	tmatrix(const size_t _n1, const size_t _n2, const size_t _n3);
	tmatrix(const size_t _n1, const size_t _n2);
	tmatrix(const size_t _n1);
	tmatrix(const tmatrix & other);
	tmatrix(initializer_list<T> list);
	tmatrix();
	~tmatrix();

	tmatrix<T> & operator = (const tmatrix<T>& other);
    tmatrix<T> & operator = (initializer_list<T> list);
	/*T& operator [](const size_t i) { return val[i]; }
	T operator [](const size_t i)const { return val[i]; }*/
	
	T first();
	T last();
	void set_zero();
	void set_all(const T a);
	void set_sequence(T offset = 0);
	void multipy_all(const T a);
	void getbox(tmatrix<T> & box, int ilower, int jlower, int iupper, int jupper);
	void setbox(tmatrix<T> & box, int ilower, int jlower, int iupper, int jupper);
	void setbox_value(T value, int ilower, int jlower, int iupper, int jupper);
	tmatrix<T> getnewbox(int ilower, int jlower, int iupper, int jupper);

    static tmatrix<T> zeros(int n1, int n2 = 1, int n3 = 1);
    static tmatrix<T> linspace(double x_0, double x_1, int size = 50);
    static tmatrix<T> logspace(double e_0, double e_1, int size = 50);
	
};


// Constructors, operator overloadings and destructors

template <class T>
tmatrix<T>::tmatrix(const size_t _n1, const size_t _n2, const size_t _n3)
{

	n1 = _n1;
	n2 = _n2;
	n3 = _n3;

	val = new T[n1 * n2 * n3];

	// remove the zero value initialization to increase code speed?
	//this->set_zero();
}

template <class T>
tmatrix<T>::tmatrix(const size_t _n1, const size_t _n2)
{

	n1 = _n1;
	n2 = _n2;
	n3 = 1;

	val = new T[n1 * n2 * n3];

	//this->set_zero();
}


template <class T>
tmatrix<T>::tmatrix(const size_t _n1)
{

	n1 = _n1;
	n2 = 1;
	n3 = 1;

	val = new T[n1 * n2 * n3];

}

template <class T>
tmatrix<T>::tmatrix(initializer_list<T> list)
{
	typename initializer_list<T>::iterator it = list.begin();

	n1 = list.size();
	n2 = 1;
	n3 = 1;

	val = new T[n1 * n2 * n3];

	for (size_t i = 0; i < n1; i++)
		val[i] = *it++;
}

template <class T>
tmatrix<T>::tmatrix() { val = new T[0]; }


template <class T>
tmatrix<T>::tmatrix(const tmatrix & other)
{

	n1 = other.n1;
	n2 = other.n2;
	n3 = other.n3;

	val = new T[n1 * n2 * n3];

	for (size_t i = 0; i < n1 * n2 * n3; i++) val[i] = other.val[i];
}

template <class T>
tmatrix<T> & tmatrix<T>::operator = (const tmatrix<T>& other)
{
	if (this != &other) { // protect against invalid self-assignment

		n1 = other.n1;
		n2 = other.n2;
		n3 = other.n3;

		delete[](val);
		val = new T[n1 * n2 * n3];

		//*val = *(other.val);
		for (size_t i = 0; i < n1 * n2 * n3; i++) val[i] = other.val[i];

	}
	// by convention, always return *this
	return *this;
}

template <class T>
tmatrix<T> & tmatrix<T>::operator = (initializer_list<T> list)
{
    typename initializer_list<T>::iterator it = list.begin();
    
    for (size_t i = 0; i < n1 * n2 * n3; i++)
        val[i] = *it++;
    
    return *this;
}


template <class T>
tmatrix<T>::~tmatrix()
{
	delete[] val;
}

// Methods

template <class T>
void tmatrix<T>::set_zero() {
	for (size_t i = 0; i < n1 * n2 * n3; i++) val[i] = 0;
}

template <class T>
void tmatrix<T>::set_all(const T a)
{
	for (size_t i = 0; i < n1 * n2 * n3; i++) val[i] = a;
}

template <class T>
void tmatrix<T>::set_sequence(T offset)
{
	for (size_t i = 0; i < n1 * n2 * n3; i++) val[i] = offset + (int) i;
}


template <class T>
T tmatrix<T>::first()
{
	return val[0];
}

template <class T>
T tmatrix<T>::last()
{
	return val[n1 * n2 * n3 - 1];
}


template <class T>
tmatrix<T> tmatrix<T>::zeros(int n1, int n2, int n3)
{
    tmatrix<T> m(n1, n2, n3);
    m.set_zero();
    return m;
}

template <class T>
tmatrix<T> tmatrix<T>::linspace(double x_0, double x_1, int size)
{
    tmatrix<T> m(size);
    for (size_t i = 0; i < size; i++) {
        m.val[i] = i * (x_1 - x_0) / (size - 1) + x_0;
    }
    
    return m;
}

template <class T>
tmatrix<T> tmatrix<T>::logspace(double e_0, double e_1, int size)
{
    tmatrix<T> m(size);
    for (size_t i = 0; i < size; i++) {
        m.val[i] = pow(10, i * (e_1 - e_0) / (size - 1) + e_0);
    }
    
    return m;
}

template <class T> 
void tmatrix<T>::getbox(tmatrix<T> & box, int ilower, int jlower, int iupper, int jupper){
	for (int i = 0; i < iupper - ilower + 1; i++)
	{
		for (int j = 0; j < jupper - jlower + 1; j++)
		{
			box.val[i * box.n2 + j] = val[(i + ilower) * n2 + (j + jlower)];
		}
	}
}

template <class T> 
void tmatrix<T>::setbox(tmatrix<T> & box, int ilower, int jlower, int iupper, int jupper){
	for (int i = 0; i < iupper - ilower + 1; i++)
	{
		for (int j = 0; j < jupper - jlower + 1; j++)
		{
			val[(i + ilower) * n2 + (j + jlower)] = box.val[i * box.n2 + j];
		}
	}
}

template <class T> 
void tmatrix<T>::setbox_value(T value, int ilower, int jlower, int iupper, int jupper){
	for (int i = 0; i < iupper - ilower + 1; i++)
	{
		for (int j = 0; j < jupper - jlower + 1; j++)
		{
			val[(i + ilower) * n2 + (j + jlower)] = value;
		}
	}
}

template <class T>
tmatrix<T> tmatrix<T>::getnewbox(int ilower, int jlower, int iupper, int jupper)
{
    tmatrix<T> m(iupper - ilower + 1, jupper - jlower + 1);

    for (int i = 0; i < iupper - ilower + 1; i++)
	{
		for (int j = 0; j < jupper - jlower + 1; j++)
		{
			m.val[i * m.n2 + j] = val[(i + ilower) * n2 + (j + jlower)];
		}
	}
	
	return m;
}

// Type definitions

typedef tmatrix<double> fmatrix;  //fmatrix stants for floating-point matrix
typedef tmatrix<int> imatrix;     //imatrix stants for integer matrix

#endif


