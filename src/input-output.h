// Printing, saving and loading

#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include <chrono>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <unordered_map>

#include "filesystem.hpp"

#include "clock.h"
#include "date.h"
#include "fmatrix.h"
#include "fields.h"
#include "state-info.h"
#include "input-output.h"
#include "configuration.h"
#include "exdir.h"
#include "diagnostics.h"

using namespace std;
using namespace std::chrono;

void verbose_log(string message, bool print);
void print_initial_info(double p_null_e, double p_null_i, configuration & config);
void print_dsmc_info(int i, int n_active_n, int step_interval, int n_steps);

void load_state(fmatrix & p_e, fmatrix & p_i, state_info & state, configuration & config, string filename="state.exdir");

// template<class T>
// void save_to_file(ghc::filesystem::path file_path, T content, ios::openmode mode=ios::app)
// {
//     ofstream file;
//     file.open(file_path, mode);
//     file << content;
//     file.close();
// }

fmatrix load_csv(string file_path, char delim = ';', int cols = 2);
template <class T>
void save_to_csv(tmatrix<T> & m, string output_path, string name = "out.csv", int nrows = -1, int ncols = -1)
{
    name = output_path + name;
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

    verbose_log("Saved " + name, true);
}


template <class T>
void print_fmatrix(tmatrix<T> & m, int n1_max=-1, int n2_max=-1){
    n1_max = n1_max < 0 ? m.n1 : n1_max;
    n2_max = n2_max < 0 ? m.n2 : n2_max;
    
    for(size_t i = 0; i < n1_max; i++){
        for(size_t j = 0; j < n2_max; j++){
            printf("%.2e  ", (double) m.val[i * m.n2 + j]);
        }
        cout << endl << endl;
    }
}

template <class T>
void save_fmatrix(tmatrix<T> & m, string filename, string dataname){
    exdir file(filename);
    file.write_dataset(dataname, m);
    verbose_log("Saved " + filename + "/" + dataname, true);
}

template <class T>
void load_fmatrix(tmatrix<T> & m, string filename, string dataname){
    exdir file(filename);
    file.read_dataset(dataname, m);
    verbose_log("Loaded " + filename + "/" + dataname, true);
}

// ------------------- output manager class ---------------------------------

class output_manager 
{
    private:
        system_clock::time_point start_utc, end_utc;
        ghc::filesystem::path output_path;
        string output_name;
        string job_name; 
        
        state_info & state;
        configuration & config;
        mesh_set & mesh;
        
        int verbosity, step_print_info, step_save_state, step_save_fields, step_save_series, step_save_fseries, step_update_metadata,
            rf_period_i, n_steps, print_timing_step, step_save_vdist, mpi_rank, mpi_size, n_mesh_x, n_mesh_y;
        double start_progress, dt;

        fmatrix wmesh_e_av, wmesh_i_av, phi_av;
        fmatrix td;
    

        void check_output_folder();
        string build_output_name(string filename_preffix);
    
    public:
        exdir file;
        output_manager(system_clock::time_point _start_utc, state_info & state, configuration & config, mesh_set & mesh);
        output_manager(string prefix, state_info & state, configuration & config, mesh_set & mesh);
        void save_state(fmatrix & p_e, fmatrix & p_i, bool force = false);
        void save_fields_snapshot(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, diagnostics & diag, mesh_set & mesh, string suffix, bool force = false);
        void save_series(diagnostics & diag, bool force = false); 
        void save_field_series(fmatrix & field, double conversion_constant, bool force = false);
        void fields_rf_average(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, diagnostics & diag, mesh_set & mesh);
        template <class T> void save_fmatrix(tmatrix<T> & m, ghc::filesystem::path datapath);
        void save_initial_data();
        void update_metadata(string status = "running", bool force = false);
        void print_info();
        void print_loop_timing(tmatrix<system_clock::time_point> & tp);
        void save_distributions(diagnostics & diag, bool force = false);
};

template <class T>
void output_manager::save_fmatrix(tmatrix<T> & m, ghc::filesystem::path datapath){
    file.write_dataset(datapath, m);
    string message = output_name / datapath;
    verbose_log("Saved " +  message, true);
}




#endif





