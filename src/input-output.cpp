
#include "input-output.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <sys/stat.h>

#include "fmatrix.h"
#include "fmath.h"

#include "state-info.h"
#include "exdir.h"
#include "configuration.h"
#include "fields.h"

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

void print_info(state_info state, int step_interval, configuration & config)
{
    static high_resolution_clock::time_point t0;
    if(state.step==state.step_offset) t0 = high_resolution_clock::now();
    if ((state.step + 1) % step_interval == 0 || state.step == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (state.step + 1 - state.step_offset) / config.i("time/n_steps")));
        printf("Step: %-8d ", state.step + 1);
        printf("Active electrons: %-8d ", state.n_active_e);
        printf("Active ions: %-8d ", state.n_active_i);
        printf("Cap. voltage: %.4f V   ", state.phi_zero / config.f("p/k_phi"));
        printf("Loop time: %.2f ms ", (double) duration_cast<microseconds>(high_resolution_clock::now() - t0).count() / (1e3 * step_interval));
        printf("\n");
        t0 = high_resolution_clock::now();
    }
}

void print_initial_info(double p_null_e, double p_null_i, configuration & config)
{
    verbose_log("\n ---- Simulation parameters ----");
    printf("Grid size:\t\t (%d, %d)\n", config.i("geometry/n_mesh_x"), config.i("geometry/n_mesh_y"));
    printf("Number of steps:\t %d\n", config.i("time/n_steps"));
    printf("Gas:\t\t\t %s\n", config.s("ions/gas_name").c_str());
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

void save_state(fmatrix & p_e, fmatrix & p_i, state_info & state, configuration & config, string suffix){
    
    const double dx = config.f("geometry/dx");
    const double dt = config.f("time/dt");
    string output_path = config.s("project/output_path");
    
    fmatrix p_e_corrected(state.n_active_e, 6);
    for(int i = 0; i < state.n_active_e; i++){
        p_e_corrected.val[i * 6 + 0] = p_e.val[i * 6 + 0] * dx;
        p_e_corrected.val[i * 6 + 1] = p_e.val[i * 6 + 1] * dx;
        p_e_corrected.val[i * 6 + 2] = p_e.val[i * 6 + 2] * dx;
        p_e_corrected.val[i * 6 + 3] = p_e.val[i * 6 + 3] * dx / dt;
        p_e_corrected.val[i * 6 + 4] = p_e.val[i * 6 + 4] * dx / dt;
        p_e_corrected.val[i * 6 + 5] = p_e.val[i * 6 + 5] * dx / dt;
    }

    fmatrix p_i_corrected(state.n_active_i, 6);
    for(int i = 0; i < state.n_active_i; i++){
        p_i_corrected.val[i * 6 + 0] = p_i.val[i * 6 + 0] * dx;
        p_i_corrected.val[i * 6 + 1] = p_i.val[i * 6 + 1] * dx;
        p_i_corrected.val[i * 6 + 2] = p_i.val[i * 6 + 2] * dx;
        p_i_corrected.val[i * 6 + 3] = p_i.val[i * 6 + 3] * dx / dt;
        p_i_corrected.val[i * 6 + 4] = p_i.val[i * 6 + 4] * dx / dt;
        p_i_corrected.val[i * 6 + 5] = p_i.val[i * 6 + 5] * dx / dt;
    }

    exdir file(output_path + "state" + suffix + ".exdir");
   
    file.write_dataset("/p_i", p_i_corrected);
    file.write_dataset("/p_e", p_e_corrected);

    map<string, double> state_attrs_double = {
        {"Time [s]", (double) state.step * dt},
        {"Capacitor voltage [V]", state.phi_zero / config.f("p/k_phi")},
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

void load_state(fmatrix & p_e, fmatrix & p_i, state_info & state, configuration & config, string filename){

    string input_path = config.s("project/input_path");
    double dx = config.f("geometry/dx");
    double dt = config.f("time/dt");


    exdir file(input_path + filename);

    YAML::Node attrs;
    file.read_all_attributes("/", attrs);

    state.step_offset = attrs["Step"].as<int>();
    state.q_cap = attrs["Capacitor charge [norm. C]"].as<double>();
    state.sigma_1 = attrs["Surface charge density [norm. C/m^2]"].as<double>();

    state.n_active_i = attrs["Active ions"].as<int>();
    state.n_out_ob_i = attrs["Removed ions"].as<int>();

    state.n_active_e = attrs["Active electrons"].as<int>();
    state.n_out_ob_e = attrs["Removed electrons"].as<int>();

    fmatrix p_i_load = fmatrix::zeros(state.n_active_i, 6);
    fmatrix p_e_load = fmatrix::zeros(state.n_active_e, 6);

    file.read_dataset("/p_e", p_e_load);
    file.read_dataset("/p_i", p_i_load);

    p_e_load = p_e_load / dx;
    for(int i = 0; i < state.n_active_e; i++){
        p_e_load.val[i * 6 + 3] = p_e_load.val[i * 6 + 3] * dt;
        p_e_load.val[i * 6 + 4] = p_e_load.val[i * 6 + 4] * dt;
        p_e_load.val[i * 6 + 5] = p_e_load.val[i * 6 + 5] * dt;
    }

    p_i_load = p_i_load / dx;
    for(int i = 0; i < state.n_active_i; i++){
        p_i_load.val[i * 6 + 3] = p_i_load.val[i * 6 + 3] * dt;
        p_i_load.val[i * 6 + 4] = p_i_load.val[i * 6 + 4] * dt;
        p_i_load.val[i * 6 + 5] = p_i_load.val[i * 6 + 5] * dt;
    }

    for(int i=0; i < state.n_active_i * 6; i++) p_i.val[i] = p_i_load.val[i];
    for(int i=0; i < state.n_active_e * 6; i++) p_e.val[i] = p_e_load.val[i];

    verbose_log("Loaded state: Step: " + to_string(state.step_offset) + " Active electrons: " + to_string(state.n_active_e) + " Active ions: " + to_string(state.n_active_i));
}

void save_fields_snapshot(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, mesh_set & mesh, state_info & state, configuration & config, string suffix)
{
    
    exdir file(config.s("project/output_path") + "fields" + suffix + ".exdir");

    double dx = config.f("geometry/dx");
    double dt = config.f("time/dt");
    double n_factor = config.f("particles/n_factor");
    double k_phi = config.f("p/k_phi");

    file.write_attribute("/", "Time [s]", (double) state.step * dt);
    file.write_attribute("/", "Step", (double) state.step);


    fmatrix phi_corrected = phi / k_phi;
    file.write_dataset("/phi", phi_corrected);

    fmatrix dens_e = (4 / pow(dx, 2)) *  n_factor * wmesh_e / mesh.v;
    file.write_dataset("/dens_e", dens_e);

    fmatrix dens_i = (4 / pow(dx, 2)) *  n_factor * wmesh_i / mesh.v;
    file.write_dataset("/dens_i", dens_i);
    
    file.create_group("/mesh");
    file.write_dataset("/mesh/mesh_x", mesh.x);
    file.write_dataset("/mesh/mesh_y", mesh.y);
    file.write_attribute("/mesh", "a_x", mesh.a_x);
    file.write_attribute("/mesh", "a_y", mesh.a_y);
    
    verbose_log("Saved fields snapshot");
}

void save_series(unordered_map<string, fmatrix> & series, int & n_points, state_info state, configuration & config, string suffix)
{
    
    exdir file(config.s("project/output_path") + "series" + suffix + ".exdir");

    double dt = config.f("time/dt");
    
    file.write_attribute("/", "Time [s]", (double) state.step * dt);
    file.write_attribute("/", "Step", (double) state.step);

    for (auto & element : series)
    {
        string path = "/" + element.first;
        file.write_dataset(path, element.second);
    }

    verbose_log("Saved time series");
}

void save_field_series(fmatrix & field, state_info state, configuration & config, double conversion_constant, string suffix)
{

    exdir file(config.s("project/output_path") + "fseries" + suffix + ".exdir");

    double dt = config.f("time/dt");
    
    fmatrix field_corrected = conversion_constant * field;

    file.write_dataset("/" + to_string(state.step), field_corrected);
    file.write_attribute("/" + to_string(state.step), "Time [s]", (double) state.step * dt);
}

