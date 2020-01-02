// Printing, saving and loading

#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include "fmatrix.h"
#include "config.h"
#include "state_info.h"
#include "input-output.h"

using namespace std;

void verbose_log(string message);
void print_info(state_info state, int step_interval);
void print_initial_info(double p_null_e, double p_null_i);
void print_dsmc_info(int i, int n_active_n, int step_interval, int n_steps);

void save_state(fmatrix & p_e, fmatrix & p_i, state_info & state);
void load_state(fmatrix & p_e, fmatrix & p_i, state_info & state);
void save_fields_snapshot(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, state_info & state, string suffix);
void save_series(map<string, fmatrix> & series, int & n_points, state_info state, string suffix);

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





