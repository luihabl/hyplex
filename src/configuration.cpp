#include "configuration.h"

#include <iomanip>
#include <typeinfo>
#include <string>
#include <unordered_map>
#include <cmath>
#include <cxxabi.h>

#include "fmatrix.h"
#include "yaml-cpp/yaml.h"

#define TAG_STR "!string"
#define TAG_DOUBLE "!double"
#define TAG_INT "!int"
#define TAG_BOOL "!bool"
#define ERROR_MSG(Map) "ERROR: Key not found in Configuration." #Map ": "

template<typename T> using dict = unordered_map<string, T>;
template<typename T> using dict_seq = unordered_map<string, tmatrix<T>>;


void configuration::calculate_parameters(){

    string calc_key = "p/";

    // ---- dx, dy -----
    double dx = f("geometry/l_x") / ((double) i("geometry/n_mesh_x") - 1); 
    double dy = f("geometry/l_y") / ((double) i("geometry/n_mesh_y") - 1);

    _map_double[calc_key + "dx"] = dx;
    _map_double[calc_key + "dy"] = dy;

    // ---- voltages ----
    double dt = f("time/dt");
    double q =  f("physical/q");
    double m_el = f("electrons/m_el");
    double n_factor = f("particles/n_factor");

    _map_double[calc_key + "volt_0_norm"] = f("boundaries/volt_0") * q * pow(dt, 2) / (m_el * pow(dx, 2));
    _map_double[calc_key + "volt_1_norm"] = f("boundaries/volt_1") * q * pow(dt, 2) / (m_el * pow(dx, 2));


    // ---- injection ----
    _map_double[calc_key + "n_inj_el"] = f("electrons/i_el") * dt / (q * n_factor);
    _map_double[calc_key + "n_inj_i"] = f("ions/i_i") * dt / (q * n_factor);
    _map_double[calc_key + "n_inj_n"] = (1 - f("thruster/eta_p")) * (4.477962e17 * f("thruster/mfc")) * dt / f("particles/n_factor_dsmc");
    
    // ---- norm constants ----

    double pi = f("physical/pi");

    _map_double[calc_key + "k_phi"] = q * dt * dt / (m_el * dx * dx);
    _map_double[calc_key + "k_q"] = f(calc_key + "k_phi") * q * n_factor / f("boundaries/c_cap");

    _map_double[calc_key + "gamma"] = 2 * n_factor * pow(q, 2) * pow(dt, 2) / (m_el * f("physical/eps_0") * pow(dx, 2));
    _map_double[calc_key + "alpha"] = (double) i("time/k_sub") * m_el / f("ugas/m_i");

    _map_double[calc_key + "omega_i"] = 2 * pi * f("thruster/freq") * dt;
    _map_double[calc_key + "rf_period_i"] = round(2 * pi / f(calc_key + "omega_i"));

    _map_double[calc_key + "k_inj_el"] =  f(calc_key + "n_inj_i") * sqrt(f("ugas/m_i") / (2 * pi * m_el));

    // ---- cross-sections ----
    tmatrix<string> exc_path = ss("ugas/exc_path");
    _map_int[calc_key + "n_exc"] = exc_path.n1;
}




template<typename T>
void configuration::select_gas_on_map(unordered_map<string, T> & map){
    
    unordered_map<string, T> map_temp;


    string gas_name = s("ions/gas_name");

    for (auto const& n : map) {
        string key = n.first;
        T value = n.second;

        string root_key = key.substr(0, key.find("/"));

        if (root_key == gas_name) {
            size_t root_key_len = root_key.length();
            key.erase(key.begin() + 0, key.begin() + root_key_len + 1);
            map_temp["ugas/" + key] = value;
        }
    }

    for (auto const& n : map_temp) {;
        map[n.first] = n.second;
    }
}




void configuration::select_gas(){
    select_gas_on_map(_map_double);
    select_gas_on_map(_map_int);
    select_gas_on_map(_map_string);
    select_gas_on_map(_map_bool);
    select_gas_on_map(_map_seq_double);
    select_gas_on_map(_map_seq_int);
    select_gas_on_map(_map_seq_string);
    select_gas_on_map(_map_seq_bool);
}

configuration::configuration(string filename)
{
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
        return _map_int.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_int) + key << endl;
        abort();
    }
}

double configuration::f(string key)
{
    try{
        return _map_double.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_double) + key << endl;
        abort();
    }
}

string configuration::s(string key)
{
    try{
        return _map_string.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_string) + key << endl;
        abort();
    }
}

bool configuration::b(string key)
{
    try{
        return _map_bool.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_bool) + key << endl;
        abort();
    }
}

tmatrix<int> configuration::is(string key)
{
    try{
        return _map_seq_int.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_seq_int) + key << endl;
        abort();
    }
}

tmatrix<double> configuration::fs(string key)
{
    try{
        return _map_seq_double.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_seq_double) + key << endl;
        abort();
    }
}

tmatrix<string> configuration::ss(string key)
{
    try{
        return _map_seq_string.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_seq_string) + key << endl;
        abort();
    }
}

tmatrix<bool> configuration::bs(string key)
{
    try{
        return _map_seq_bool.at(key);
    }
    catch (const out_of_range & e){
        cerr << ERROR_MSG(_map_seq_bool) + key << endl;
        abort();
    }
}


void configuration::set_value(YAML::Node node, string key)
{

    string tag = node.Tag();

    if (tag == TAG_INT) {
        _map_int[key] = get_int(node);
    }

    else if (tag == TAG_DOUBLE) {
        _map_double[key] = get_double(node);
    }

    else if (tag == TAG_STR) {
        _map_string[key] = get_string(node);
    }

    else if (tag == TAG_BOOL) {
        _map_bool[key] = get_bool(node);
    }
}

void configuration::set_sequence(YAML::Node node, string key){
    
    string tag = node.Tag();
    int size = node.size();
       
    if (tag == TAG_INT) {
        imatrix int_seq = imatrix::zeros(size);
        for (int i = 0; i < size; i++) 
            int_seq.val[i] = get_int(node[i]);
        _map_seq_int[key] = int_seq;
    }
    
    else if (tag == TAG_DOUBLE) {
        fmatrix double_seq = fmatrix::zeros(size);
        for (int i = 0; i < size; i++) 
            double_seq.val[i] = get_double(node[i]);
        _map_seq_double[key] = double_seq;
    }
    
    else if (tag == TAG_STR) {
        tmatrix<string> string_seq(size);
        for (int i = 0; i < size; i++) 
            string_seq.val[i] = get_string(node[i]);
        _map_seq_string[key] = string_seq;
    }
    
    else if (tag == TAG_BOOL) {
        tmatrix<bool> bool_seq(size);
        for (int i = 0; i < size; i++) 
            bool_seq.val[i] = get_bool(node[i]);
        _map_seq_bool[key] = bool_seq;
    }
}


string lowercase(string s)
{
    for (auto& c : s)
        c = tolower(c);
    return s;
}

int get_int(YAML::Node node){
    int int_val;

    string str_val = node.as<string>();

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

double get_double(YAML::Node node){
    string str_val = node.as<string>();
    return stod(str_val);
}

string get_string(YAML::Node node){
    return node.as<string>();
}

bool get_bool(YAML::Node node){
    string str_val = node.as<string>();
    str_val = lowercase(str_val);

    if (str_val == "true" || str_val == "1") {
        return true;
    }
    else {
        return false;
    }
}

template<typename T>
void configuration::print_map(dict<T> & map){
    int s;
    string types = abi::__cxa_demangle(typeid(T).name(), 0, 0, &s);
    if(typeid(T) == typeid(string)) types = "string";

    for (auto const& x : map) {

        cout << std::left;
        cout << setw(35) << x.first << "(" + types + ") " << x.second << endl;
    }

    
}

template<typename T>
void configuration::print_seq(dict_seq<T> & map){
    int s;
    string types = abi::__cxa_demangle(typeid(T).name(), 0, 0, &s);
    if(typeid(T) == typeid(string)) types = "string";
    for (auto const& x : map) {
        cout << x.first << endl;
        for(size_t i = 0; i < x.second.n1; i++){
            cout << std::left;
            cout << setw(35) << "  " <<  "(" + types + ") " << x.second.val[i] << endl;
        }
    }
}

void configuration::print_all(){
    print_map(_map_double);
    print_map(_map_int);
    print_map(_map_string);
    print_map(_map_bool);
    print_seq(_map_seq_double);
    print_seq(_map_seq_int);
    print_seq(_map_seq_string);
    print_seq(_map_seq_bool);
}