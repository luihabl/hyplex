#ifndef CLOCK_H
#define CLOCK_H

#include <chrono>
#include <iostream>
#include <string>
#include <sstream>

#include "date.h"
#include "fmatrix.h"

using namespace std;
using namespace std::chrono;
using namespace date;

namespace clk {

template <class T>
double tdiff_s(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count() / 1e6;
}

template <class T>
double tdiff_min(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count() / (1e6 * 60);
}


template <class T>
double tdiff_h(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count() / (1e6 * 60 * 60);
}

template <class T>
double tdiff_ms(T t1, T t2){
    auto duration = duration_cast<microseconds>(t2 - t1);
    return (double) duration.count() / (1e3);
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

inline void print_td(fmatrix & td, int n_td){
    for(int j=0; j<n_td; j++){
        cout << j << "-" <<  j+1 << ": " << td.val[j] << "\t";
    }
    cout << endl;
}

}

#endif
