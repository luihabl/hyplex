#ifndef NUM_TOOLS_H
#define NUM_TOOLS_H

#include "fmatrix.h"

double interp(const fmatrix & data, double x);
fmatrix interp(const fmatrix & data, fmatrix & x);

imatrix sample_from_sequence_shuffle(int sample_size, int range);
imatrix sample_from_sequence_naive(int sample_size, int range);

void average_field(fmatrix & av_field, const fmatrix & field, int step);
int average_field_over_period(fmatrix & av_field, const fmatrix & field, int period, int total_period, int counter);

inline
void swap(double& a, double& b) {
    double temp = a;
    a = b;
    b = temp;
}

template <class T>
T clamp_n(T low, T hi, T val){
    if (val < low) {return low;}
    else if (val > hi) {return hi;}
    else {return val;}
}

template <class T>
bool hasnan(tmatrix<T> & m, int max_n=-1){
    max_n = (int) (max_n < 0 ? m.n1 * m.n2 * m.n3 : max_n);
    for(int i = 0; i < max_n; i++){
        if(isnan(m.val[i]) || isinf(m.val[i])){
            return true;
        }
    }
    return false;
}



#endif