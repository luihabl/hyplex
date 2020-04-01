#ifndef CLOCK_H
#define CLOCK_H

#include <chrono>
#include "date.h"
#include <iostream>

using namespace std;
using namespace std::chrono;


template <class T>
double tdiff(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count() / 1e3;
}

template <class T>
double tdiff_us(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count();
}

inline high_resolution_clock::time_point now(){
    return high_resolution_clock::now();
}



#endif
