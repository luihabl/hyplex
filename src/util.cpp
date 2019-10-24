
#include "util.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>

#include "fmatrix.h"
#include "fmath.h"

using namespace std;
using namespace std::chrono;

void verbose_log(string message){
    #ifdef VERBOSE
    cout << message << endl;
    #endif
}

void print_dsmc_info(int i, int n_active_n, int step_interval, int n_steps){
    static high_resolution_clock::time_point t0;
    if(i==0) t0 = high_resolution_clock::now();
    if ((i + 1) % step_interval == 0 || i == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (i + 1) / n_steps));
        printf("Step: %-8d ", i + 1);
        printf("Active neutrals: %-8d ", n_active_n);
        printf("Loop time: %.2f ms ", (double) duration_cast<microseconds>(high_resolution_clock::now() - t0).count() / (1e3 * step_interval));
        printf("\n");
        t0 = high_resolution_clock::now();
    }
}

void print_info(int i, int step_offset, fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i, double v_cap, int step_interval)
{
    static high_resolution_clock::time_point t0;
    if(i==0) t0 = high_resolution_clock::now();
    if ((i + 1) % step_interval == 0 || i == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (i + 1 - step_offset) / N_STEPS));
        printf("Step: %-8d ", i + 1);
        printf("Active electrons: %-8d ", n_active_e);
        printf("Active ions: %-8d ", n_active_i);
        printf("Cap. voltage: %.4f V   ", v_cap);
        printf("Loop time: %.2f ms ", (double) duration_cast<microseconds>(high_resolution_clock::now() - t0).count() / (1e3 * step_interval));
        printf("\n");
        t0 = high_resolution_clock::now();
    }
}

void print_initial_info(double p_null_e, double p_null_i)
{
    verbose_log("\n ---- Simulation parameters ----");
    printf("Grid size:\t\t (%d, %d)\n", N_MESH_X, N_MESH_Y);
    printf("Number of steps:\t %d\n", N_STEPS);
    printf("Gas:\t\t\t %s\n", GAS_NAME.c_str());
    printf("P Null (e):\t\t %.4e\n", p_null_e);
    printf("P Null (I):\t\t %.4e\n\n", p_null_i);
}

fmatrix load_csv(string file_path, char delim, int cols)
{
    string line;
    ifstream file(file_path);
    
    int line_count = 0;
    while (getline(file, line)) line_count++;
    file.clear();
    file.seekg(0, ios::beg);
    
    fmatrix data(line_count, cols); 
    int i = 0;
    while (getline(file, line))
    {
        istringstream s(line);
        string field;
        int j = 0;
        while (getline(s, field, delim))
        {
            if (j > cols - 1) {
                cout << "csv exceed j limit!" << endl;
                break;
            }
            data.val[i * cols + j] = atof(field.c_str());
            j++;
        }
        i++;
    }
    
    if (!file.eof()) {
        cout << "Could not read file " << endl;
    }
    return data;
}

void save_state(fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i,  int i, fmatrix & misc, string suffix){
    
    
    fmatrix p_e_corrected = DX * p_e;
    for(int i = 0; i < n_active_e; i++){
        p_e_corrected.val[i * 6 + 3] = p_e_corrected.val[i * 6 + 3] / DT;
        p_e_corrected.val[i * 6 + 4] = p_e_corrected.val[i * 6 + 4] / DT;
        p_e_corrected.val[i * 6 + 5] = p_e_corrected.val[i * 6 + 5] / DT;
    }
    save_to_csv(p_e_corrected, "p_e" + suffix + ".csv", n_active_e);

    fmatrix p_i_corrected = DX * p_i;
    for(int i = 0; i < n_active_i; i++){
        p_i_corrected.val[i * 6 + 3] = p_i_corrected.val[i * 6 + 3] / DT;
        p_i_corrected.val[i * 6 + 4] = p_i_corrected.val[i * 6 + 4] / DT;
        p_i_corrected.val[i * 6 + 5] = p_i_corrected.val[i * 6 + 5] / DT;
    }
    save_to_csv(p_i_corrected, "p_i" + suffix + ".csv", n_active_i);
    
    fmatrix sim_state(1 + misc.n1);
   
    sim_state.val[0] = (double) i;
    for (int j = 1; j < (int) sim_state.n1; j++) sim_state.val[i] = misc.val[i - 1];
   
    save_to_csv(sim_state, "sim" + suffix + ".csv");
    
    verbose_log("State saved");
}

void load_state(fmatrix & p_e, int & n_active_e, fmatrix & p_i, int & n_active_i, int & step_offset, fmatrix & misc, string suffix){

    fmatrix p_i_load = load_csv("output/p_i" + suffix + ".csv",',', 6);
    fmatrix p_e_load = load_csv("output/p_e" + suffix + ".csv",',', 6);

    n_active_i = (int) p_i_load.n1;
    n_active_e = (int) p_e_load.n1;

    p_e_load = p_e_load / DX;
    for(int i = 0; i < n_active_e; i++){
        p_e_load.val[i * 6 + 3] = p_e_load.val[i * 6 + 3] * DT;
        p_e_load.val[i * 6 + 4] = p_e_load.val[i * 6 + 4] * DT;
        p_e_load.val[i * 6 + 5] = p_e_load.val[i * 6 + 5] * DT;
    }

    p_i_load = p_i_load / DX;
    for(int i = 0; i < n_active_i; i++){
        p_i_load.val[i * 6 + 3] = p_i_load.val[i * 6 + 3] * DT;
        p_i_load.val[i * 6 + 4] = p_i_load.val[i * 6 + 4] * DT;
        p_i_load.val[i * 6 + 5] = p_i_load.val[i * 6 + 5] * DT;
    }

    for(int i=0; i < n_active_i * 6; i++) p_i.val[i] = p_i_load.val[i];
    for(int i=0; i < n_active_e * 6; i++) p_e.val[i] = p_e_load.val[i];

    fmatrix sim_state = load_csv("output/sim" + suffix + ".csv",',', 1);
    for(int i=1; i < (int) sim_state.n1; i++) misc.val[i-1] = sim_state.val[i];
    step_offset = sim_state.val[0];

    verbose_log("State loaded: i: " + to_string(step_offset) + " Active electrons: " + to_string(n_active_e) + " Active ions: " + to_string(n_active_i));
}

void save_fields(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, string suffix){
        
    fmatrix phi_corrected = phi * (M_EL * pow(DX, 2))/(Q * pow(DT, 2));
    save_to_csv(phi_corrected, "phi" + suffix + ".csv");

    fmatrix dens_e = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_e / vmesh;
    save_to_csv(dens_e, "dens_e" + suffix + ".csv");

    fmatrix dens_i = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_i / vmesh;
    save_to_csv(dens_i, "dens_i" + suffix + ".csv");
}


