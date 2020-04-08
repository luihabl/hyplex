
#include "input-output.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <filesystem>
#include <sys/stat.h>

#include "yaml-cpp/yaml.h"

#include "fmatrix.h"
#include "fmath.h"

#include "state-info.h"
#include "exdir.h"
#include "configuration.h"
#include "fields.h"

#include "date.h"


using namespace std;
using namespace std::chrono;

void verbose_log(string message, bool print){
    if(print) cout << message << endl;
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

void print_info(state_info & state, int step_interval, configuration & config)
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
    verbose_log("\n ---- Simulation parameters ----", config.i("simulation/verbosity") >= 1);
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

void output_manager::save_state(fmatrix & p_e, fmatrix & p_i){
    check_output_folder();
    const double dx = config.f("geometry/dx");
    const double dt = config.f("time/dt");
    
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


    file.create_group("state");
   
    file.write_dataset("state/p_i", p_i_corrected);
    file.write_dataset("state/p_e", p_e_corrected);

    map<string, double> state_attrs_double = {
        {"Time [s]", (double) state.step * dt},
        {"Capacitor voltage [V]", state.phi_zero / config.f("p/k_phi")},
        {"Capacitor charge [norm. C]", state.q_cap},
        {"Surface charge density [norm. C/m^2]", state.sigma_1}
    };
    

    file.write_attribute("state", state_attrs_double);

    map<string, int> state_attrs_int = {
        {"Step", state.step},
        {"Active ions", state.n_active_i},
        {"Removed ions", state.n_out_ob_i},
        {"Active electrons", state.n_active_e},
        {"Removed electrons", state.n_out_ob_e}
    };

    file.write_attribute("state", state_attrs_int);

    verbose_log("Saved state", verbosity >= 1);
}

void load_state(fmatrix & p_e, fmatrix & p_i, state_info & state, configuration & config, string filename){

    filesystem::path input_path = config.s("project/input_path");
    input_path /= filename;

    double dx = config.f("geometry/dx");
    double dt = config.f("time/dt");

    exdir file(input_path, false);

    YAML::Node attrs;
    file.read_all_attributes("", attrs);

    state.step_offset = attrs["Step"].as<int>();
    state.q_cap = attrs["Capacitor charge [norm. C]"].as<double>();
    state.sigma_1 = attrs["Surface charge density [norm. C/m^2]"].as<double>();

    state.n_active_i = attrs["Active ions"].as<int>();
    state.n_out_ob_i = attrs["Removed ions"].as<int>();

    state.n_active_e = attrs["Active electrons"].as<int>();
    state.n_out_ob_e = attrs["Removed electrons"].as<int>();

    fmatrix p_i_load = fmatrix::zeros(state.n_active_i, 6);
    fmatrix p_e_load = fmatrix::zeros(state.n_active_e, 6);

    
    file.read_dataset("p_e", p_e_load);
    file.read_dataset("p_i", p_i_load);

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

    verbose_log("Loaded state: Step: " + to_string(state.step_offset) + " Active electrons: " + to_string(state.n_active_e) + " Active ions: " + to_string(state.n_active_i), true);
}

void output_manager::save_fields_snapshot(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, mesh_set & mesh, string suffix)
{
    check_output_folder();
    file.create_group("fields" + suffix);

    double dx = config.f("geometry/dx");
    double dt = config.f("time/dt");
    double n_factor = config.f("particles/n_factor");
    double k_phi = config.f("p/k_phi");

    file.write_attribute("fields", "Time [s]", (double) state.step * dt);
    file.write_attribute("fields", "Step", (double) state.step);


    fmatrix phi_corrected = phi / k_phi;
    file.write_dataset("fields/phi", phi_corrected);

    fmatrix dens_e = (4 / pow(dx, 2)) *  n_factor * wmesh_e / mesh.v;
    file.write_dataset("fields/dens_e", dens_e);

    fmatrix dens_i = (4 / pow(dx, 2)) *  n_factor * wmesh_i / mesh.v;
    file.write_dataset("fields/dens_i", dens_i);
    
    verbose_log("Saved fields snapshot", verbosity >= 1);
}

void output_manager::save_series(unordered_map<string, fmatrix> & series, int & n_points)
{
    check_output_folder();

    file.create_group("series");

    double dt = config.f("time/dt");
    
    file.write_attribute("series", "Time [s]", (double) state.step * dt);
    file.write_attribute("series", "Step", (double) state.step);

    for (auto & element : series)
    {
        string path = "series/" + element.first;
        file.write_dataset(path, element.second);
    }

    verbose_log("Saved time series", verbosity >= 1);
}

void output_manager::save_field_series(fmatrix & field, double conversion_constant)
{
    check_output_folder();
    file.create_group("field_series");

    double dt = config.f("time/dt");
    
    fmatrix field_corrected = conversion_constant * field;

    file.write_dataset("field_series/" + to_string(state.step), field_corrected);
    file.write_attribute("field_series/" + to_string(state.step), "Time [s]", (double) state.step * dt);
}

output_manager::output_manager(system_clock::time_point _start_utc, state_info & _state, configuration & _config, mesh_set & _mesh): state(_state), config(_config), mesh(_mesh){
    
    start_utc = _start_utc;
    output_path = config.s("project/output_path");
    job_name = config.s("p/job_name");

    output_name = build_output_name();
    
    file = exdir(output_path / output_name, false);    
    verbosity = config.i("simulation/verbosity");
}

string output_manager::build_output_name(){
    string initial_datetime = format("%Y-%m-%d", start_utc);    
    string filename_preffix = "hy_" + initial_datetime;

    if(config.s("p/job_name") != "") filename_preffix += "_" + config.s("p/job_name");

    string tmp_filename = filename_preffix + EXDIR_EXT;

    int n_file = 1;
    while (filesystem::exists(output_path / tmp_filename))
    {
        tmp_filename = filename_preffix + "(" + to_string(n_file) + ")" + EXDIR_EXT;
        n_file += 1;
    }
    return tmp_filename;
}

void output_manager::save_initial_data(){

    double dx = config.f("geometry/dx");
    fmatrix mesh_x = mesh.x * dx;
    fmatrix mesh_y = mesh.y * dx;

    file.create_group("mesh");
    file.write_dataset("mesh/x", mesh_x);
    file.write_dataset("mesh/y", mesh_y);

    YAML::Node config_node = YAML::LoadFile(config.filename);

    YAML::Node attributes = file.get_attributes("");

    attributes["metadata"]["version"] = GIT_VERSION;
    attributes["metadata"]["start_utc"] = time_to_string(start_utc);
    attributes["metadata"]["stop_utc"] = time_to_string(start_utc);
    attributes["metadata"]["initial_step"] = state.step_offset;
    attributes["config"] = config_node;
    attributes["config"]["simulation"]["job_name"] = config.s("p/job_name");
    
    file.write_attribute("", attributes);

    update_metadata();
}

void output_manager::update_metadata(string status){
    YAML::Node attributes = file.get_attributes("");
    attributes["metadata"]["stop_utc"] = time_to_string(sys_now());
    attributes["metadata"]["elapsed_hours"] = tdiff_h(start_utc, sys_now());
    attributes["metadata"]["status"] = status;
    attributes["metadata"]["progress"] = float (state.step - state.step_offset) / config.f("time/n_steps");
    file.write_attribute("", attributes);
}

void output_manager::check_output_folder(){
    if(!file.file_exists()){  
        file = exdir(output_path / output_name, false);    
        save_initial_data();
        update_metadata();
    }
}


