
#include "input-output.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <map>
#include <sys/stat.h>

#include "fmatrix.h"
#include "fmath.h"

#include "state-info.h"
// #include "h5io.h"
// #include "H5Cpp.h"
#include "exdir.h"

using namespace std;
using namespace std::chrono;
// using namespace H5;

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

void print_info(state_info state, int step_interval)
{
    static high_resolution_clock::time_point t0;
    if(state.step==state.step_offset) t0 = high_resolution_clock::now();
    if ((state.step + 1) % step_interval == 0 || state.step == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (state.step + 1 - state.step_offset) / N_STEPS));
        printf("Step: %-8d ", state.step + 1);
        printf("Active electrons: %-8d ", state.n_active_e);
        printf("Active ions: %-8d ", state.n_active_i);
        printf("Cap. voltage: %.4f V   ", state.phi_zero / K_PHI);
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

void save_state(fmatrix & p_e, fmatrix & p_i, state_info & state, string suffix){
    
    
    fmatrix p_e_corrected = DX * p_e;
    for(int i = 0; i < state.n_active_e; i++){
        p_e_corrected.val[i * 6 + 3] = p_e_corrected.val[i * 6 + 3] / DT;
        p_e_corrected.val[i * 6 + 4] = p_e_corrected.val[i * 6 + 4] / DT;
        p_e_corrected.val[i * 6 + 5] = p_e_corrected.val[i * 6 + 5] / DT;
    }

    fmatrix p_i_corrected =  DX * p_i;
    for(int i = 0; i < state.n_active_i; i++){
        p_i_corrected.val[i * 6 + 3] = p_i_corrected.val[i * 6 + 3] / DT;
        p_i_corrected.val[i * 6 + 4] = p_i_corrected.val[i * 6 + 4] / DT;
        p_i_corrected.val[i * 6 + 5] = p_i_corrected.val[i * 6 + 5] / DT;
    }

    exdir file(OUTPUT_PATH + "state" + suffix + ".exdir");
   
    p_i_corrected.n1 = state.n_active_i;
    file.write_dataset("/p_i", p_i_corrected);

    p_e_corrected.n1 = state.n_active_e;
    file.write_dataset("/p_e", p_e_corrected);

    map<string, double> state_attrs_double = {
        {"Time [s]", (double) state.step * DT},
        {"Capacitor voltage [V]", state.phi_zero / K_PHI},
        {"Capacitor charge [norm. C]", state.q_cap},
        {"Surface charge density [norm. C/m^2]", state.sigma_1}
    };

    file.write_attribute("/", state_attrs_double);

    map<string, int> state_attrs_int = {
        {"Step", state.step},
        {"Active ions", state.n_active_i},
        {"Removed ions", state.n_out_ob_i},
        {"Active electrons", state.n_active_e},
        {"Removed electrons", state.n_out_ob_e}
    };

    file.write_attribute("/", state_attrs_int);

    verbose_log("Saved state");
}

void load_state(fmatrix & p_e, fmatrix & p_i, state_info & state, string filename){

    // H5File file(INPUT_PATH + "state.h5", H5F_ACC_RDONLY);
    exdir file(INPUT_PATH + filename);

    // DataSet p_i_dataset = file.openDataSet("p_i");
    // DataSet p_e_dataset = file.openDataSet("p_e");

    YAML::Node attrs;
    file.read_all_attributes("/", attrs);

    state.step_offset = attrs["Step"].as<int>();
    state.q_cap = attrs["Capacitor charge [norm. C]"].as<double>();
    state.sigma_1 = attrs["Surface charge density [norm. C/m^2]"].as<double>();

    state.n_active_i = attrs["Active ions"].as<int>();
    state.n_out_ob_i = attrs["Removed ions"].as<int>();

    state.n_active_e = attrs["Active electrons"].as<int>();
    state.n_out_ob_e = attrs["Removed electron"].as<int>();

    // read_attribute(file, "Step", state.step_offset);
    // read_attribute(file, "Capacitor charge [norm. C]", state.q_cap);
    // read_attribute(file, "Surface charge density [norm. C/m^2]", state.sigma_1);
    
    // read_attribute(p_i_dataset, "Active ions", state.n_active_i);
    // read_attribute(p_i_dataset, "Removed ions", state.n_out_ob_i);

    // read_attribute(p_e_dataset, "Active electrons", state.n_active_e);
    // read_attribute(p_e_dataset, "Removed electrons", state.n_out_ob_e);

    fmatrix p_i_load = fmatrix::zeros(state.n_active_i, 6);
    fmatrix p_e_load = fmatrix::zeros(state.n_active_e, 6);

    file.read_dataset("/p_e", p_e_load);
    file.read_dataset("/p_i", p_i_load);

    // p_i_dataset.read(p_i_load.val, PredType::NATIVE_DOUBLE);
    // p_e_dataset.read(p_e_load.val, PredType::NATIVE_DOUBLE);

    p_e_load = p_e_load / DX;
    for(int i = 0; i < state.n_active_e; i++){
        p_e_load.val[i * 6 + 3] = p_e_load.val[i * 6 + 3] * DT;
        p_e_load.val[i * 6 + 4] = p_e_load.val[i * 6 + 4] * DT;
        p_e_load.val[i * 6 + 5] = p_e_load.val[i * 6 + 5] * DT;
    }

    p_i_load = p_i_load / DX;
    for(int i = 0; i < state.n_active_i; i++){
        p_i_load.val[i * 6 + 3] = p_i_load.val[i * 6 + 3] * DT;
        p_i_load.val[i * 6 + 4] = p_i_load.val[i * 6 + 4] * DT;
        p_i_load.val[i * 6 + 5] = p_i_load.val[i * 6 + 5] * DT;
    }

    for(int i=0; i < state.n_active_i * 6; i++) p_i.val[i] = p_i_load.val[i];
    for(int i=0; i < state.n_active_e * 6; i++) p_e.val[i] = p_e_load.val[i];

    verbose_log("Loaded state: Step: " + to_string(state.step_offset) + " Active electrons: " + to_string(state.n_active_e) + " Active ions: " + to_string(state.n_active_i));
}

void save_fields_snapshot(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, state_info & state, string suffix){
    
    exdir file(OUTPUT_PATH + "fields" + suffix + ".exdir");

    file.write_attribute("/", "Time [s]", (double) state.step * DT);
    file.write_attribute("/", "Step", (double) state.step);


    fmatrix phi_corrected = phi * (M_EL * pow(DX, 2))/(Q * pow(DT, 2));
    file.write_dataset("/phi", phi_corrected);

    fmatrix dens_e = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_e / vmesh;
    file.write_dataset("/dens_e", dens_e);

    fmatrix dens_i = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_i / vmesh;
    file.write_dataset("/dens_i", dens_i);

    verbose_log("Saved fields snapshot");
}

void save_series(map<string, fmatrix> & series, int & n_points, state_info state, string suffix){
    
    exdir file(OUTPUT_PATH + "series" + suffix + ".exdir");
    
    file.write_attribute("/", "Time [s]", (double) state.step * DT);
    file.write_attribute("/", "Step", (double) state.step);

    map<string, fmatrix>::iterator cursor;

    for (auto & element : series)
    {
        string path = "/" + element.first;
        file.write_dataset(path, element.second);
    }

    verbose_log("Saved time series");
}

void save_field_series(fmatrix & field, state_info state, double conversion_constant, string suffix){

    exdir file(OUTPUT_PATH + "fseries" + suffix + ".exdir");
    
    fmatrix field_corrected = conversion_constant * field;

    file.write_dataset("/" + to_string(state.step), field_corrected);
    file.write_attribute("/" + to_string(state.step), "Time [s]", (double) state.step * DT);
}

void save_fmatrix(fmatrix & m, string filename, string dataname){
    exdir file(filename);
    file.write_dataset("/" + dataname, m);
    verbose_log("Saved " + filename + "/" + dataname);
}

void load_fmatrix(fmatrix & m, string filename, string dataname){
    exdir file(filename);
    file.read_dataset("/dataname", m);
    verbose_log("Loaded " + filename + "/" + dataname);
}
