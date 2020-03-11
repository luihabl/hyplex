#ifndef CONFIGC_H
#define CONFIGC_H

#include <string>
#include <unordered_map>

#include "yaml-cpp/yaml.h"
#include "fmatrix.h"

using namespace std;

template<typename T> using dict = unordered_map<string, T>;
template<typename T> using dict_seq = unordered_map<string, tmatrix<T>>;

class configuration {

    private:
        
        dict<double>        _map_double;
        dict<int>           _map_int;
        dict<string>        _map_string;
        dict<bool>          _map_bool;

        dict_seq<double>    _map_seq_double;
        dict_seq<int>       _map_seq_int;
        dict_seq<string>    _map_seq_string;
        dict_seq<bool>      _map_seq_bool;


        void set_value(YAML::Node node, string key);
        void set_sequence(YAML::Node node, string key);
        void flatten_map(YAML::Node node, string parent_key = "");
        void calculate_parameters();
        void select_gas();
        template<typename T> void select_gas_on_map(unordered_map<string, T> & map);

    public:
        configuration(string filename);
        void print_all();
        template<typename T> void print_map(dict<T> & map);
        template<typename T> void print_seq(dict_seq<T> & map);
        
        int     i(string key);
        double  f(string key);
        string  s(string key);
        bool    b(string key);

        tmatrix<int>     is(string key);
        tmatrix<double>  fs(string key);
        tmatrix<string>  ss(string key);
        tmatrix<bool>    bs(string key);

};


int get_int(YAML::Node node);
bool get_bool(YAML::Node node);
double get_double(YAML::Node node);
string get_string(YAML::Node node);
string lowercase(string s);

#endif