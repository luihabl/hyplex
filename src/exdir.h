#ifndef EXDIR_H
#define EXDIR_H

#include <string>
#include <cstdlib>
#include <fstream>
#include <map>

#include "yaml-cpp/yaml.h"
#include "npy.h"

#define VERSION "1"

using namespace std;
class exdir
{
    private:
        string root_file_path;
        bool overwrite_file = false;
        
        int mkdir(string path);
        int rm_file(string path);
        int create_folder(string path, string type);
        void write_metadata(string path, string type, string version);
        void ensure_attr_exists(string object_path);
        bool is_type(string path, string type);
    
    public:
        exdir(string _file_path, bool _overwrite_file = false);
        void create_group(string group_path);
        template <class T> void write_dataset(string dataset_path, tmatrix<T> & data);
        template <class T> void write_attribute(string object_path, string key, T value);
        template <class T> void write_attribute(string object_path, tmatrix<string> & key, tmatrix<T> & value);
        template <class T> void write_attribute(string object_path, map<string, T> & dict);
        template <class T> void read_dataset(string dataset_path, tmatrix<T> & data);
        template <class T> void read_attribute(string object_path, string key, T & value);
        void get_dataset_size(string dataset_path, tmatrix<unsigned long> & shape);
        void read_all_attributes(string object_path, YAML::Node & node);
        void clean_attr_file(string object_path);

};

exdir::exdir(string _file_path, bool _overwrite_file){
    root_file_path = _file_path;
    overwrite_file = _overwrite_file;
    
    if(overwrite_file && is_type(root_file_path, "file")) rm_file(root_file_path);
    else if(overwrite_file && !is_type(root_file_path, "file")) cout << "Invalid file, cannot overwrite." << endl;
    
    create_folder(root_file_path, "file");
}

int exdir::mkdir(string path){
    string cmd = "mkdir -p " + path;
    return system(cmd.c_str());
}

int exdir::rm_file(string path){
    string cmd = "rm -r -f " + path;
    return system(cmd.c_str());
}


int exdir::create_folder(string path, string type){
    int dir = mkdir(path);
    write_metadata(path, type, VERSION);
    return dir;
}   


bool exdir::is_type(string path, string type){
    try { 
        YAML::Node node = YAML::LoadFile(path + "/exdir.yaml");
        return node["exdir"]["type"].as<string>() == type;
    
    } catch (const std::exception& e) {
        return false;
    }
    return false;
}


void exdir::write_metadata(string path, string type, string version){
   
    YAML::Emitter out;
    
    out << YAML::BeginMap;
        out << YAML::Key << "exdir";
        out << YAML::Value;
        out << YAML::BeginMap;
                out << YAML::Key << "type" << YAML::Value << YAML::DoubleQuoted << type;
                out << YAML::Key << "version" << YAML::Value << version;
        out << YAML::EndMap;
	out << YAML::EndMap;

    
    ofstream fout;
    fout.open(path + "/exdir.yaml");
    fout << out.c_str();
    fout.close();
       
}

void exdir::create_group(string group_path){
    create_folder(root_file_path + group_path, "group");
}

template <class T>
void exdir::write_dataset(string dataset_path, tmatrix<T> & data){
    
    create_folder(root_file_path + dataset_path, "dataset");
    
    tmatrix<unsigned long> shape = {data.n1, data.n2, data.n3};
    if (shape.val[2] == 1)  
    {
        shape.n1 -= 1;
        if (shape.val[1] == 1)
            shape.n1 -= 1;
    }
    
    npy::SaveArrayAsNumpy(root_file_path + dataset_path + "/data.npy", false, shape, data);
}

void exdir::clean_attr_file(string object_path){
    
    string attr_path = root_file_path + object_path + "/attributes.yaml";
    
    ofstream fs;
    fs.open(attr_path, ofstream::out | ofstream::trunc);
    fs.close();
}

void exdir::ensure_attr_exists(string object_path){
    
    string attr_path = root_file_path + object_path + "/attributes.yaml";
    
    ofstream fs;
    fs.open(attr_path, ofstream::app);
    fs.close();
}

template <class T>
void exdir::write_attribute(string object_path, string key, T value) {
    
    ensure_attr_exists(object_path);
 
    string attr_path = root_file_path + object_path + "/attributes.yaml";
       
    YAML::Node node = YAML::LoadFile(attr_path);
    node[key] = value;
    
    ofstream fout(attr_path);
    fout << node;
    fout.close();
    
}


template <class T>
void exdir::write_attribute(string object_path, tmatrix<string> & key, tmatrix<T> & value) {
    
    ensure_attr_exists(object_path);
    
    string attr_path = root_file_path + object_path + "/attributes.yaml";

    YAML::Node node = YAML::LoadFile(attr_path);
    
    for(int k = 0; k < key.n1; k++){
        node[key.val[k]] = value.val[k];
    }
   
    ofstream fout(attr_path);
    fout << node;
    fout.close();
}

template <class T>
void exdir::write_attribute(string object_path, map<string, T> & dict) {
    
    ensure_attr_exists(object_path);
    
    string attr_path = root_file_path + object_path + "/attributes.yaml";

    YAML::Node node = YAML::LoadFile(attr_path);
    
    for (auto const& element : dict)
    {
        node[element.first] = element.second;

    }
   
    ofstream fout(attr_path);
    fout << node;
    fout.close();
}




void exdir::get_dataset_size(string dataset_path, tmatrix<unsigned long> & shape){
    string data_path = root_file_path + dataset_path + "/data.npy";
    bool fortran_order = false;
    npy::GetArrayShape(data_path, shape, fortran_order);
}


template <class T>
void exdir::read_dataset(string dataset_path, tmatrix<T> & data){
    
    string data_path = root_file_path + dataset_path + "/data.npy";
    
    tmatrix<unsigned long> shape(3);
    
    bool fortran_order = false;
    npy::LoadArrayFromNumpy(data_path, shape, fortran_order, data);
}

template <class T>
void exdir::read_attribute(string object_path, string key, T & value) {
    string attr_path = root_file_path + object_path + "/attributes.yaml";
    
    YAML::Node node = YAML::LoadFile(attr_path);
    value = node[key].as<T>();       
}


void exdir::read_all_attributes(string object_path, YAML::Node & node) {
    string attr_path = root_file_path + object_path + "/attributes.yaml";
    node = YAML::LoadFile(attr_path);
}



#endif