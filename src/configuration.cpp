#include "configuration.h"

#include <iomanip>
#include <typeinfo>
#include <string>
#include <unordered_map>
#include <cmath>
#include <cxxabi.h>

#include "fmatrix.h"
#include "yaml-cpp/yaml.h"

#define ERROR_MSG(Map, Type) "ERROR: Key not found in Configuration." #Map "(" Type "): "



void configuration::calculate_parameters(){

    string calc_key = "p/";

    // ---- dx, dy -----
//    double dx = f("geometry/dx");
//    double dy = f("geometry/dy");
    
    double dx = f("geometry/l_x") / ((double) i("geometry/n_mesh_x") - 1);
    double dy = f("geometry/l_y") / ((double) i("geometry/n_mesh_y") - 1);

    set("geometry/dx", dx);
    set("geometry/dy", dy);

    // ---- voltages ----
    double dt = f("time/dt");
    double q =  f("physical/q");
    double m_el = f("electrons/m_el");
    double n_factor = f("particles/n_factor");

    set(calc_key + "volt_0_norm", f("boundaries/volt_0") * q * pow(dt, 2) / (m_el * pow(dx, 2)));
    set(calc_key + "volt_1_norm", f("boundaries/volt_1") * q * pow(dt, 2) / (m_el * pow(dx, 2)));

    // ---- injection ----
    set(calc_key + "n_inj_el", f("electrons/i_el") * dt / (q * n_factor));
    set(calc_key + "n_inj_i", f("ions/i_i") * dt / (q * n_factor));
    set(calc_key + "n_inj_n", (1 - f("thruster/eta_p")) * (4.477962e17 * f("thruster/mfc")) * dt / f("particles/n_factor_dsmc"));
    
    // ---- norm constants ----

    double pi = f("physical/pi");

    set(calc_key + "k_phi", q * dt * dt / (m_el * dx * dx));
    set(calc_key + "k_q", f(calc_key + "k_phi") * q * n_factor / f("boundaries/c_cap"));

    set(calc_key + "gamma", 2 * n_factor * pow(q, 2) * pow(dt, 2) / (m_el * f("physical/eps_0") * pow(dx, 2)));
    set(calc_key + "alpha", (double) i("time/k_sub") * m_el / f("ugas/m_i"));

    set(calc_key + "omega_i", 2 * pi * f("thruster/freq") * dt);
    set(calc_key + "rf_period_i", round(2 * pi / f(calc_key + "omega_i")));

    set(calc_key + "k_inj_el",  f(calc_key + "n_inj_i") * sqrt(f("ugas/m_i") / (2 * pi * m_el)));

    // ---- cross-sections ----
    smatrix exc_path = ss("ugas/exc_path");
    set(calc_key + "n_exc", (int) exc_path.n1);
}

void configuration::set_job_name(string arg_job_name){
    set("p/job_name", arg_job_name == "" ? s("simulation/job_name") : arg_job_name);
}


void configuration::select_gas(){
    select_gas_on_map(_m);
    select_gas_on_map(_ms);
}

configuration::configuration(string _filename)
{
    filename = _filename;
    YAML::Node config_node = YAML::LoadFile(filename);
    flatten_map(config_node);
    select_gas();
    calculate_parameters();
}


void configuration::flatten_map(YAML::Node node, string parent_key)
{

    for (auto const& n : node) {
        YAML::Node local_node = n.second;
        string new_key = parent_key == "" ? n.first.as<string>() : parent_key + "/" + n.first.as<string>();

        if (local_node.IsMap()) 
            flatten_map(local_node, new_key);
        
        if (n.second.IsScalar()) 
            set_value(local_node, new_key);
        
        if (n.second.IsSequence())
            set_sequence(local_node, new_key);
    }
}

int configuration::i(string key)
{
    try{
        return to_int(_m.at(key));
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_m,"int") + key << endl;
        abort();
    }
}

double configuration::f(string key)
{
    try{
        return to_double(_m.at(key));
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_m, "double") + key << endl;
        abort();
    }
}

string configuration::s(string key)
{
    try{
        return _m.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_m, "string") + key << endl;
        abort();
    }
}

bool configuration::b(string key)
{
    try{
        return to_bool(_m.at(key));
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_m, "bool") + key << endl;
        abort();
    }
}

tmatrix<int> configuration::is(string key)
{
    try{
        return convert_seq<int>(_ms[key], to_int);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_ms, "int") + key << endl;
        abort();
    }
}

tmatrix<double> configuration::fs(string key)
{
    try{
        return convert_seq<double>(_ms[key], to_double);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_ms, "double") + key << endl;
        abort();
    }
}

tmatrix<string> configuration::ss(string key)
{
    try{
        return convert_seq<string>(_ms[key], pass_string);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_ms, "string") + key << endl;
        abort();
    }
}

tmatrix<bool> configuration::bs(string key)
{
    try{
        return convert_seq<bool>(_ms[key], to_bool);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_ms, "bool") + key << endl;
        abort();
    }
}


void configuration::set_value(YAML::Node node, string key)
{
    _m[key] = node.as<string>();
}

void configuration::set_sequence(YAML::Node node, string key){

    int size = (int) node.size();
    
    smatrix seq(size);
    for (int i = 0; i < size; i++) 
        seq.val[i] = node[i].as<string>();
    _ms[key] = seq;
}


string lowercase(string s)
{
    for (auto& c : s)
        c = tolower(c);
    return s;
}

int to_int(string str_val){
    int int_val;
    string e = "e";
    if (str_val.find(e) != string::npos) {
        double value_double = stod(str_val);
        int_val = (int)value_double;
    }
    else {
        int_val = stoi(str_val);
    }
    
    return int_val;
}

double to_double(string str_val){
    return stod(str_val);
}



bool to_bool(string str_val){
    str_val = lowercase(str_val);

    if (str_val == "true" || str_val == "1") {
        return true;
    }
    else {
        return false;
    }
}

string pass_string(string str_val){
    return str_val;
}




void configuration::print_m(){
    for (auto const& x : _m) {

        cout << std::left;
        cout << setw(35) << x.first << x.second << endl;
    } 
}

void configuration::print_ms(){
    for (auto const& x : _ms) {
        cout << x.first << endl;
        for(size_t i = 0; i < x.second.n1; i++){
            cout << std::left;
            cout << setw(35) << "  " << x.second.val[i] << endl;
        }
    }
}

void configuration::print_all(){
    print_m();
    print_ms();
}
