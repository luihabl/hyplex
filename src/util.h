#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "fmatrix.h"
#include "config.h"
#include "util.h"

#define GET_VARIABLE_NAME(Variable) (#Variable)


using namespace std;

double interp(const fmatrix & data, double x);
fmatrix interp(const fmatrix & data, fmatrix & x);
fmatrix load_csv(string file_path, char delim = ';');
imatrix sample_from_sequence_shuffle(int sample_size, int range);
imatrix sample_from_sequence_naive(int sample_size, int range);
// imatrix sample_from_sequence_swap(int sample_size, int range);
void print_info(int i, fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i,  int step_interval = 1000);
void average_field(fmatrix & av_field, const fmatrix & field, int i);
int clamp(int low, int hi, int val);
void verbose_log(string message);

inline
void swap(double& a, double& b) {
    double temp = a;
    a = b;
    b = temp;
}

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





