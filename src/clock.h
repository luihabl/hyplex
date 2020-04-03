#ifndef CLOCK_H
#define CLOCK_H

#include <chrono>
#include <iostream>
#include <string>
#include <sstream>

#include "date.h"

using namespace std;
using namespace std::chrono;
using namespace date;

template <class T>
double tdiff(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count() / 1e3;
}

template <class T>
double tdiff_h(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count() / (1e6 * 60 * 60);
}

template <class T>
double tdiff_us(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count();
}

inline high_resolution_clock::time_point now(){
    return high_resolution_clock::now();
}

inline system_clock::time_point sys_now(){
    return system_clock::now();
} 

template <class T>
string time_to_string(T t){
    stringstream ss;
    ss << t;
    return ss.str();
}


inline string get_utc_datetime_string()
{
    system_clock::time_point t = std::chrono::system_clock::now();
    stringstream ss;
    ss << t;
    return ss.str();
}




#endif
