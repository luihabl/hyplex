#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <sstream>
#include <iomanip>
#include <unordered_map>

#include "yaml-cpp/yaml.h"
#include "filesystem.hpp" 
#include "fmatrix.h"

#define ERROR_MSG(Map, Type) "ERROR: Key not found in Configuration." #Map "(" Type "): "

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
        ghc::filesystem::path filename;
        configuration(ghc::filesystem::path filename);
        void print_all();
        void print_m();
        void print_ms();
        
        int     i(string key);
        double  f(string key);
        string  s(string key);
        bool    b(string key);
        // template<typename T> T get(string key);

        tmatrix<int>     is(string key);
        tmatrix<double>  fs(string key);
        tmatrix<string>  ss(string key);
        tmatrix<bool>    bs(string key);

        template<typename T> void set(string key, T val);
        template<typename T> void set_seq(string key, tmatrix<T> seq);

        void set_job_name(string arg_job_name);
};

string lowercase(string s);
bool to_bool(string str_val);
double to_double(string str_val);
int to_int(string str_val);
string pass_string(string str_val);

// template <typename T> T convert_to (const string & str)
// {
//     if constexpr(is_integral<T>::value)
//     {
//         istringstream ss(str);
        
//         if constexpr(is_same<T, string>::value){
//             return str;
//         }
        
//         if(str.find("e") != string::npos){
//             double num;
//             ss >> num;
//             return static_cast<T>(num);
//         }

//         T num;
//         ss >> num;
//         return num;
        
//     }
//     else
//     {
//         istringstream ss(str);
//         T num;
//         ss >> num;
//         return num;
//     }
// }

// template<typename T>
// T configuration::get(string key)
// {
//     try{
//         return convert_to<T>(_m.at(key));
//     }
//     catch (const out_of_range & e){
//         cerr << ERROR_MSG(_m, "type") + key << endl;
//         abort();
//     }
// }


template <typename T, typename U>
tmatrix<T> convert_seq(smatrix & in, U func)
{
    size_t size = in.n1 * in.n2 * in.n3;
    tmatrix<T> out(size);

    for(size_t i=0; i < size; i++){
        out.val[i] = func(in.val[i]);
    }

    return out;
}

template<typename T>
void configuration::set(string key, T val)
{
    ostringstream strs;
    strs << setprecision(numeric_limits<T>::digits10);
    strs << val;
    _m[key] = strs.str();
}



template<typename T>
void configuration::set_seq(string key, tmatrix<T> seq)
{
    if(typeid(T) == typeid(string)){
        _ms[key] = seq;
    }
    else
    {
        _ms[key] = convert_seq<T>(_ms[key], to_string);
    }
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


#endif
