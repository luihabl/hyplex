#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "fmatrix.h"
#include "config.h"
#include "util.h"

using namespace std;

void verbose_log(string message);
void print_info(int i, int step_offset, fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i, double v_cap, int step_interval);
void print_initial_info(double p_null_e, double p_null_i);
void print_dsmc_info(int i, int n_active_n, int step_interval, int n_steps);

void save_state(fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i,  int i, fmatrix & misc, string suffix);
void load_state(fmatrix & p_e, int & n_active_e, fmatrix & p_i, int & n_active_i, int & step_offset, fmatrix & misc, string suffix);
fmatrix load_csv(string file_path, char delim = ';', int cols = 2);

template <class T>
void save_to_csv(tmatrix<T> & m, string name = "out.csv", int nrows = -1, int ncols = -1)
{
    name = OUTPUT_PATH + name;
    nrows = nrows < 0 ? (int) m.n1 : nrows;
    ncols = ncols < 0 ? (int) m.n2 : ncols;

    ofstream file(name);
    
    for (int i = 0; i < nrows; i++) {
        
        for (int j = 0; j < ncols; j++) {
            ostringstream stream_obj;
            stream_obj << setprecision(16);
            stream_obj << m.val[i * m.n2 + j];
            file << stream_obj.str();
            if (j != ncols - 1) file << ",";
        }
        file << endl;
    }
    
    file.close();

    verbose_log("Saved " + name);
}


template <class T>
void print_fmatrix(tmatrix<T> & phi){
    for(size_t i = 0; i < phi.n1; i++){
        for(size_t j = 0; j < phi.n2; j++){
            printf("%-7.2f", (double) phi.val[i * phi.n2 + j]);
        }
        cout << endl << endl;
    }
}

#endif





