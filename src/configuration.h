#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <unordered_map>

#include "yaml-cpp/yaml.h"
#include "fmatrix.h"

using namespace std;

class configuration {

    private:
        
        unordered_map<string, string> _m;
        unordered_map<string, smatrix> _ms;

        void set_value(YAML::Node node, string key);
        void set_sequence(YAML::Node node, string key);
        void flatten_map(YAML::Node node, string parent_key = "");
        void calculate_parameters();
        void select_gas();
        template<typename T> void select_gas_on_map(unordered_map<string, T> & map);

    public:
        configuration(string filename);
        void print_all();
        void print_m();
        void print_ms();
        
        int     i(string key);
        double  f(string key);
        string  s(string key);
        bool    b(string key);

        tmatrix<int>     is(string key);
        tmatrix<double>  fs(string key);
        tmatrix<string>  ss(string key);
        tmatrix<bool>    bs(string key);

        template<typename T> void set(string key, T val);
        template<typename T> void set_seq(string key, tmatrix<T> seq);

};

string lowercase(string s);
template <typename T, typename U> tmatrix<T> convert_seq(smatrix & in, U func);
bool to_bool(string str_val);
double to_double(string str_val);
int to_int(string str_val);
string pass_string(string str_val);

#endif