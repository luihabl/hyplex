
#include "input-output.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <unordered_map>

#include "filesystem.hpp"

#include <sys/stat.h>

#include "yaml-cpp/yaml.h"

#include "fmatrix.h"
#include "fmath.h"

#include "state-info.h"
#include "exdir.h"
#include "configuration.h"
#include "fields.h"
#include "num-tools.h"
#include "date.h"
#include "mpi.h"
#include "diagnostics.h"

using namespace std;
using namespace std::chrono;

namespace io {

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

void output_manager::print_info()
{
    
    int n_active_e_total, n_active_i_total;
    
//    MPI_Reduce(&state->n_active_e, &n_active_e_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&state->n_active_i, &n_active_i_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    n_active_e_total = state->n_active_e * mpi_size; // Shows the approximate number for better performance
    n_active_i_total = state->n_active_i * mpi_size;
    
    if(mpi_rank !=0 ) return;
    
    static high_resolution_clock::time_point t0;
    if(state->step==state->step_offset) t0 = high_resolution_clock::now();
    if ((state->step + 1) % step_print_info == 0 || state->step == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (state->step + 1 - state->step_offset) / config->i("time/n_steps")));
        printf("Step: %-8d ", state->step + 1);
        printf("Active electrons: %-8d ", n_active_e_total);
        printf("Active ions: %-8d ", n_active_i_total);
        printf("Cap. voltage: %.4f V   ", state->phi_zero / config->f("p/k_phi"));
        printf("Loop time: %.2f ms ", (double) duration_cast<microseconds>(high_resolution_clock::now() - t0).count() / (1e3 * step_print_info));
        printf("\n");
        t0 = high_resolution_clock::now();
    }
}

void print_initial_info(double p_null_e, double p_null_i, configuration & config)
{
    
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(mpi_rank !=0 ) return;
    
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

void output_manager::save_state(fmatrix & p_e, fmatrix & p_i, bool force){

    if(!state_save) return;
    if(!(force || state->step % step_save_state == 0)) return;

    check_output_folder();
    const double dx = config->f("geometry/dx");
    const bool save_p_i = config->b("diagnostics/state/save_p_i");
    const bool save_p_e = config->b("diagnostics/state/save_p_e");

    int n_active_e_total=0, n_active_i_total=0,
        n_out_ob_e_total=0, n_out_ob_i_total=0;

    MPI_Reduce(&state->n_active_e, &n_active_e_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&state->n_active_i, &n_active_i_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&state->n_out_ob_e, &n_out_ob_e_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&state->n_out_ob_i, &n_out_ob_i_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    imatrix n_elements = imatrix::zeros(mpi_size);
    imatrix displacement = imatrix::zeros(mpi_size);

    fmatrix p_e_corrected(n_active_e_total, 6);
    fmatrix p_i_corrected(n_active_i_total, 6);

    MPI_Gather(&state->n_active_e, 1, MPI_INT, n_elements.val, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < mpi_size; i++ ) n_elements.val[i] = n_elements.val[i] * 6;
    displacement.val[0] = 0;
    for (int i = 1; i < mpi_size; i++ ) displacement.val[i] = displacement.val[i-1] + n_elements.val[i-1];
    MPI_Gatherv(p_e.val, state->n_active_e * 6, MPI_DOUBLE, p_e_corrected.val,
                n_elements.val, displacement.val, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    MPI_Gather(&state->n_active_i, 1, MPI_INT, n_elements.val, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < mpi_size; i++ ) n_elements.val[i] = n_elements.val[i] * 6;
    displacement.val[0] = 0;
    for (int i = 1; i < mpi_size; i++ ) displacement.val[i] = displacement.val[i-1] + n_elements.val[i-1];
    MPI_Gatherv(p_i.val, state->n_active_i * 6, MPI_DOUBLE, p_i_corrected.val,
                n_elements.val, displacement.val, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(mpi_rank != 0) return;


    for(int i = 0; i < n_active_e_total; i++){
        p_e_corrected.val[i * 6 + 0] = p_e_corrected.val[i * 6 + 0] * dx;
        p_e_corrected.val[i * 6 + 1] = p_e_corrected.val[i * 6 + 1] * dx;
        p_e_corrected.val[i * 6 + 2] = p_e_corrected.val[i * 6 + 2] * dx;
        p_e_corrected.val[i * 6 + 3] = p_e_corrected.val[i * 6 + 3] * dx / dt;
        p_e_corrected.val[i * 6 + 4] = p_e_corrected.val[i * 6 + 4] * dx / dt;
        p_e_corrected.val[i * 6 + 5] = p_e_corrected.val[i * 6 + 5] * dx / dt;
    }


    for(int i = 0; i < n_active_i_total; i++){
        p_i_corrected.val[i * 6 + 0] = p_i_corrected.val[i * 6 + 0] * dx;
        p_i_corrected.val[i * 6 + 1] = p_i_corrected.val[i * 6 + 1] * dx;
        p_i_corrected.val[i * 6 + 2] = p_i_corrected.val[i * 6 + 2] * dx;
        p_i_corrected.val[i * 6 + 3] = p_i_corrected.val[i * 6 + 3] * dx / dt;
        p_i_corrected.val[i * 6 + 4] = p_i_corrected.val[i * 6 + 4] * dx / dt;
        p_i_corrected.val[i * 6 + 5] = p_i_corrected.val[i * 6 + 5] * dx / dt;
    }

    file.create_group("state");

    if (save_p_i) file.write_dataset("state/p_i", p_i_corrected);
    if (save_p_e) file.write_dataset("state/p_e", p_e_corrected);
    
    map<string, double> state_attrs_double;
    state_attrs_double = {
        {"Time [s]", (double) state->step * dt},
        {"Capacitor voltage [V]", state->phi_zero / config->f("p/k_phi")},
        {"Capacitor charge [norm. C]", state->q_cap},
        {"Surface charge density [norm. C/m^2]", state->sigma_1}
    };

    file.write_attribute("state", state_attrs_double);

    map<string, int> state_attrs_int;
    state_attrs_int = {
        {"Step", state->step},
        {"Active ions", n_active_i_total},
        {"Removed ions", n_out_ob_i_total},
        {"Active electrons", n_active_e_total},
        {"Removed electrons", n_out_ob_e_total}
    };

    file.write_attribute("state", state_attrs_int);

    verbose_log("Saved state", verbosity >= 1);
}

void load_state(fmatrix & p_e, fmatrix & p_i, state_info & state, configuration & config, string filename){
    
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ghc::filesystem::path input_path = config.s("project/input_path");
    input_path /= filename;

    double dx = config.f("geometry/dx");
    double dt = config.f("time/dt");

    exdir file(input_path, false);

    YAML::Node attrs;
    file.read_all_attributes("", attrs);

    state.step_offset = attrs["Step"].as<int>();
    state.q_cap = attrs["Capacitor charge [norm. C]"].as<double>();
    state.sigma_1 = attrs["Surface charge density [norm. C/m^2]"].as<double>();

    if (mpi_rank == 0){
        state.n_out_ob_i = attrs["Removed ions"].as<int>();
        state.n_out_ob_e = attrs["Removed electrons"].as<int>();
    }
    
    int n_total_e = attrs["Active electrons"].as<int>();
    int part_e = n_total_e / mpi_size;
    int mod_e =  n_total_e % mpi_size;

    int n_total_i = attrs["Active ions"].as<int>();
    int part_i = n_total_i / mpi_size;
    int mod_i =  n_total_i % mpi_size;

    state.n_active_i = part_i;
    state.n_active_e = part_e;
    if (mpi_rank == mpi_size - 1) {
        state.n_active_i = part_i + mod_i;
        state.n_active_e = part_e + mod_e;
    }
    
    fmatrix p_i_load;
    fmatrix p_e_load;
    
    if(mpi_rank == 0){
        p_i_load = fmatrix::zeros(n_total_i, 6);
        p_e_load = fmatrix::zeros(n_total_e, 6);
        file.read_dataset("p_e", p_e_load);
        file.read_dataset("p_i", p_i_load);
    }
    
    imatrix n_elements   = imatrix::zeros(mpi_size);
    imatrix displacement = imatrix::zeros(mpi_size);
    
    n_elements.set_all(part_i * 6);
    n_elements.val[mpi_size - 1] = (part_i + mod_i) * 6;
    displacement.val[0] = 0;
    for (int i = 1; i < mpi_size; i++ ) displacement.val[i] = displacement.val[i-1] + n_elements.val[i-1];
    MPI_Scatterv(p_i_load.val, n_elements.val, displacement.val, MPI_DOUBLE,
                 p_i.val, state.n_active_i * 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    n_elements.set_all(part_e * 6);
    n_elements.val[mpi_size - 1] = (part_e + mod_e) * 6;
    displacement.val[0] = 0;
    for (int i = 1; i < mpi_size; i++ ) displacement.val[i] = displacement.val[i-1] + n_elements.val[i-1];
    MPI_Scatterv(p_e_load.val, n_elements.val, displacement.val, MPI_DOUBLE,
                 p_e.val, state.n_active_e * 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for(int i = 0; i < state.n_active_e; i++){
        p_e.val[i * 6 + 0] = p_e.val[i * 6 + 0] / dx;
        p_e.val[i * 6 + 1] = p_e.val[i * 6 + 1] / dx;
        p_e.val[i * 6 + 2] = p_e.val[i * 6 + 2] / dx;
        p_e.val[i * 6 + 3] = p_e.val[i * 6 + 3] * dt / dx;
        p_e.val[i * 6 + 4] = p_e.val[i * 6 + 4] * dt / dx;
        p_e.val[i * 6 + 5] = p_e.val[i * 6 + 5] * dt / dx;
    }

    for(int i = 0; i < state.n_active_i; i++){
        p_i.val[i * 6 + 0] = p_i.val[i * 6 + 0] / dx;
        p_i.val[i * 6 + 1] = p_i.val[i * 6 + 1] / dx;
        p_i.val[i * 6 + 2] = p_i.val[i * 6 + 2] / dx;
        p_i.val[i * 6 + 3] = p_i.val[i * 6 + 3] * dt / dx;
        p_i.val[i * 6 + 4] = p_i.val[i * 6 + 4] * dt / dx;
        p_i.val[i * 6 + 5] = p_i.val[i * 6 + 5] * dt / dx;
    }

    verbose_log("Loaded state: Step: " + to_string(state.step_offset) + " Active electrons: " + to_string(state.n_active_e) + " Active ions: " + to_string(state.n_active_i), true);


}

void output_manager::save_fields_snapshot(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, diagnostics & diag, mesh_set & mesh, string suffix, bool force)
{
    if(!(force || state->step % step_save_fields == 0)) return;

    diag.reduce_izfield();

    if(mpi_rank != 0) return;

    check_output_folder();

    ghc::filesystem::path obj_path = "fields" + suffix;

    file.create_group(obj_path);

    const double dx = config->f("geometry/dx");
    const double dt = config->f("time/dt");
    const double n_factor = config->f("particles/n_factor");
    const double k_phi = config->f("p/k_phi");
    const double m_el = config->f("electrons/m_el");
    const double m_i = config->f("ugas/m_i");
    const double q = config->f("physical/q");

    file.write_attribute(obj_path, "Time [s]", (double) state->step * dt);
    file.write_attribute(obj_path, "Step", (double) state->step);

    fmatrix phi_corrected = phi / k_phi;
    file.write_dataset(obj_path / "phi", phi_corrected);

    fmatrix k_dens = mesh.v / ((4 / pow(dx, 2)) *  n_factor);

    fmatrix dens_e = wmesh_e / k_dens;
    file.write_dataset(obj_path / "dens_e", dens_e);

    fmatrix dens_i = wmesh_i / k_dens;
    file.write_dataset(obj_path / "dens_i", dens_i);


    ghc::filesystem::path ufield_group = obj_path / "u";

    file.create_group(ufield_group);
    
    fmatrix ufield_e_x_corrected = (dx / dt) * diag.ufield_e_x_global; 
    file.write_dataset(ufield_group / "u_e_x", ufield_e_x_corrected);

    fmatrix ufield_e_y_corrected = (dx / dt) * diag.ufield_e_y_global ; 
    file.write_dataset(ufield_group / "u_e_y", ufield_e_y_corrected);

    fmatrix ufield_i_x_corrected = (dx / dt) * diag.ufield_i_x_global; 
    file.write_dataset(ufield_group / "u_i_x", ufield_i_x_corrected);

    fmatrix ufield_i_y_corrected = (dx / dt) * diag.ufield_i_y_global; 
    file.write_dataset(ufield_group / "u_i_y", ufield_i_y_corrected);


    ghc::filesystem::path pfield_group = obj_path / "E_k";

    file.create_group(pfield_group);
    
    fmatrix kfield_e_corrected = (0.5 * m_el * (dx * dx) / (dt * dt * q)) * diag.kfield_e_global; 
    file.write_dataset(pfield_group / "E_e", kfield_e_corrected);

    fmatrix kfield_i_corrected = (0.5 * m_i * (dx * dx) / (dt * dt * q)) * diag.kfield_i_global; 
    file.write_dataset(pfield_group / "E_i", kfield_i_corrected);


    ghc::filesystem::path izfield_group = obj_path / "coll";

    file.create_group(izfield_group);

    fmatrix izfield_corrected = diag.izfield_global / (k_dens * dt * (double) step_save_fields);
    file.write_dataset(izfield_group / "ionization", diag.izfield_global);

    verbose_log("Saved fields snapshot", verbosity >= 1);
}



void output_manager::save_series(diagnostics & diag, bool force)
{
    if(!series_save) return;
    if(!(force || state->step % step_save_series == 0)) return;

    if(mpi_rank == 0){
        
        check_output_folder();

        file.create_group("series");

        file.write_attribute("series", "Time [s]", (double) state->step * dt);
        file.write_attribute("series", "Step", (double) state->step);
    
        for (int i = 0; i < (int) diag.gseries_keys.n1; i++){
            string path = "series/" + diag.gseries_keys.val[i];
            file.write_dataset(path, diag.gseries[diag.gseries_keys.val[i]]);
        }
    }
    
    for (int i = 0; i < (int) diag.lseries_keys.n1; i++){
        string path = "series/" + diag.lseries_keys.val[i];
        MPI_Reduce(diag.lseries[diag.lseries_keys.val[i]].val, diag.tmp_array.val, diag.series_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(mpi_rank == 0) file.write_dataset(path, diag.tmp_array);
    }
    
    if(mpi_rank != 0) return;
    
    verbose_log("Saved time series", verbosity >= 1);
}


output_manager::output_manager(system_clock::time_point _start_utc, state_info * _state, configuration * _config, mesh_set * _mesh): state(_state), config(_config), mesh(_mesh){
    
    start_utc = _start_utc;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    output_path = config->s("project/output_path");
    job_name = config->s("p/job_name");

    rf_period_i = config->i("p/rf_period_i");
    n_steps = config->i("time/n_steps");
    start_progress = config->f("diagnostics/rf_av/start_progress");
    print_timing_step = config->i("diagnostics/print_info/print_timing_step"); 
    dt = config->f("time/dt");
    
    step_file_refresh = config->i("diagnostics/output_file/file_refresh_step");
    step_print_info  = config->i("diagnostics/print_info/print_step");
    step_save_state  = config->i("diagnostics/state/save_step");
    step_save_fields = config->i("diagnostics/fields_snapshot/save_step");
    step_save_series = config->i("diagnostics/series/save_step");
    step_update_metadata = config->i("diagnostics/metadata/update_step");
    n_mesh_x = config->i("geometry/n_mesh_x");
    n_mesh_y = config->i("geometry/n_mesh_y");

    output_overwrite = config->b("diagnostics/output_file/overwrite");
    state_save = config->b("diagnostics/state/save");
    series_save = config->b("diagnostics/series/save");
    
    step_save_vdist = config->i("diagnostics/vdist/save_step");

    int nx = config->i("geometry/n_mesh_x");
    int ny = config->i("geometry/n_mesh_y");

    wmesh_e_av = fmatrix::zeros(nx, ny);
    wmesh_i_av = fmatrix::zeros(nx, ny);
    phi_av = fmatrix::zeros(nx, ny);
    td = fmatrix::zeros(50);
    td.set_zero();
    verbosity = config->i("simulation/verbosity");
    
    start_new_file();
}


output_manager::output_manager(string prefix, state_info * _state, configuration * _config, mesh_set * _mesh, bool dsmc): state(_state), config(_config), mesh(_mesh){
    
    start_utc = clk::sys_now();

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    output_path = config->s("project/output_path");
    job_name = config->s("p/job_name");

    verbosity = config->i("simulation/verbosity");

    if (dsmc) n_steps = config->i("time/n_steps_dsmc");
    else n_steps = config->i("time/n_steps");
    
    start_new_file(prefix);
}

void output_manager::start_new_file(string prefix){

    if(mpi_rank == 0) output_name = build_output_name(prefix);
    
    int str_size = (int) output_name.size();
    MPI_Bcast(&str_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank != 0) output_name.resize(str_size);
    
    MPI_Bcast(const_cast<char*>(output_name.data()), str_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    file = exdir(output_path / output_name, false);
    
    save_initial_data();
}

void output_manager::start_new_file(){    
    start_new_file("hy_" + format("%Y-%m-%d", start_utc));
}

void output_manager::refresh_file(){
    if(!output_overwrite && (state->step % step_file_refresh == 0)) start_new_file();
}

string output_manager::build_output_name(string filename_preffix){

    if(config->s("p/job_name") != "") filename_preffix += "_" + config->s("p/job_name");

    string tmp_filename = filename_preffix + EXDIR_EXT;

    int n_file = 1;
    while (ghc::filesystem::exists(output_path / tmp_filename))
    {
        tmp_filename = filename_preffix + "-" + to_string(n_file) + EXDIR_EXT;
        n_file += 1;
    }
    return tmp_filename;
}

void output_manager::save_initial_data(){
    
    if(mpi_rank !=0 ) return;

    double dx = config->f("geometry/dx");
    fmatrix mesh_x = mesh->x * dx;
    fmatrix mesh_y = mesh->y * dx;

    file.create_group("mesh");
    file.write_dataset("mesh/x", mesh_x);
    file.write_dataset("mesh/y", mesh_y);

    YAML::Node config_node = YAML::LoadFile(config->filename);

    YAML::Node attributes = file.get_attributes("");

    attributes["metadata"]["version"] = GIT_VERSION;
    attributes["metadata"]["start_utc"] = clk::time_to_string(start_utc);
    attributes["metadata"]["stop_utc"] = clk::time_to_string(start_utc);
    attributes["metadata"]["initial_step"] = state->step_offset;
    attributes["config"] = config_node;
    attributes["config"]["simulation"]["job_name"] = config->s("p/job_name");
    
    file.write_attribute("", attributes);

    update_metadata();
}

void output_manager::update_metadata(string status, bool force){

    if(!(force || state->step % step_update_metadata == 0)) return;
    if(mpi_rank != 0) return;

    YAML::Node attributes = file.get_attributes("");
    attributes["metadata"]["stop_utc"] = clk::time_to_string(clk::sys_now());
    attributes["metadata"]["elapsed_hours"] = clk::tdiff_h(start_utc, clk::sys_now());
    attributes["metadata"]["status"] = status;
    attributes["metadata"]["progress"] = (float) (state->step - state->step_offset) / (float) n_steps;
    file.write_attribute("", attributes);

    verbose_log("Metadata updated", verbosity >= 1);
}

void output_manager::check_output_folder(){
    if(!file.file_exists()){
        file = exdir(output_path / output_name, false);
        save_initial_data();
    }
}


void output_manager::fields_rf_average(fmatrix & phi, fmatrix & wmesh_e, fmatrix & wmesh_i, diagnostics & diag, mesh_set & mesh){
    
    if(mpi_rank != 0) return;

    int i = state->step - state->step_offset;
    if(i > start_progress * (double) n_steps){
        
        int counter_offset = n_steps % rf_period_i;
        int counter_av = ((i - counter_offset) % rf_period_i) + 1;

        if(i >= counter_offset){
            average_field(phi_av, phi, counter_av);
            average_field(wmesh_e_av, wmesh_e, counter_av);
            average_field(wmesh_i_av, wmesh_i, counter_av);
	    }

        if(counter_av == rf_period_i){
            save_fields_snapshot(phi_av, wmesh_e_av, wmesh_i_av, diag, mesh, "-rfav", true);
        }
    }
}


 void output_manager::print_loop_timing(tmatrix<system_clock::time_point> & tp){
     
    if(mpi_rank != 0) return;
    
    if(verbosity >= 2)
    {
        for(size_t j=0; j < tp.n1 - 1; j++){
            td.val[j] += clk::tdiff_ms(tp.val[j], tp.val[j+1]) / (double) print_timing_step;
        }

        if((state->step+1) % print_timing_step == 0){
            clk::print_td(td, (int) tp.n1 - 1);
            td.set_zero();
        }
    }
}

void output_manager::save_distributions(diagnostics & diag, bool force){

    if(!(force || state->step % step_save_vdist == 0)) return;

    check_output_folder();

    ghc::filesystem::path obj_path = "dist", obj_path_e = obj_path / "electrons", obj_path_i = obj_path / "ions";

    diag.reduce_distributions();
    
    if(mpi_rank == 0){
        file.create_group(obj_path);
        file.create_group(obj_path_e);
        file.create_group(obj_path_i);

        file.write_dataset(obj_path_i / "vx", diag.vdist_i_global_x);
        file.write_dataset(obj_path_i / "vy", diag.vdist_i_global_y);
        file.write_dataset(obj_path_e / "vx", diag.vdist_e_global_x);
        file.write_dataset(obj_path_e / "vy", diag.vdist_e_global_y);

        file.write_dataset(obj_path_i / "i_top", diag.top_dist_i_global);
        file.write_dataset(obj_path_i / "i_rhs", diag.rhs_dist_i_global);
        file.write_dataset(obj_path_e / "i_top", diag.top_dist_e_global);
        file.write_dataset(obj_path_e / "i_rhs", diag.rhs_dist_e_global);
    }
    
    if(mpi_rank != 0) return;

    file.write_attribute(obj_path, "Time [s]", (double) state->step * dt);
    file.write_attribute(obj_path, "Step", (double) state->step);

    file.write_attribute(obj_path_e, "n_v", diag.n_v_e);
    file.write_attribute(obj_path_e, "vlim_x", {diag.vlim_e.val[0], diag.vlim_e.val[1]});
    file.write_attribute(obj_path_e, "vlim_y", {diag.vlim_e.val[2], diag.vlim_e.val[3]});
    
    file.write_attribute(obj_path_i, "n_v", diag.n_v_i);
    file.write_attribute(obj_path_i, "vlim_x", {diag.vlim_i.val[0], diag.vlim_i.val[1]});
    file.write_attribute(obj_path_i, "vlim_y", {diag.vlim_i.val[2], diag.vlim_i.val[3]});
    
    verbose_log("Saved distributions", verbosity >= 1);
}



}