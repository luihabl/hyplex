#ifndef EXDIR_H
#define EXDIR_H

#include <string>
#include <filesystem>
#include <cstdlib>
#include <fstream>
#include <map>

#include "yaml-cpp/yaml.h"
#include "fmatrix.h"
#include "npy.h"

#define VERSION "1"
#define ATTR_FILE "attributes.yaml"
#define METADATA_FILE "exdir.yaml"
#define DATA_FILE "data.npy"
#define EXDIR_EXT ".exdir"

using namespace std;
class exdir
{
    private:
        void create_folder(filesystem::path path, string type);
        void write_metadata(filesystem::path path, string type, string version);
        void ensure_attr_exists(filesystem::path object_path);
        template<class T> void save_to_file(filesystem::path file_path, T content, ios::openmode mode=ios::app);
        string build_file_name(string file_preffix);
        bool is_type(filesystem::path path, string type);
    
    public:
        filesystem::path file_path;
        
        exdir(filesystem::path _file_path, bool _overwrite_file=false);
        exdir(const exdir & other);
        exdir() {};
        exdir & operator = (const exdir & other);
        void create_group(filesystem::path group_path);
        template <class T> void write_dataset(filesystem::path dataset_path, tmatrix<T> & data);
        YAML::Node get_attributes(filesystem::path object_path);
        void read_all_attributes(filesystem::path object_path, YAML::Node & node);
        template <class T> void write_attribute(filesystem::path object_path, string key, T value);
        template <class T> void write_attribute(filesystem::path object_path, tmatrix<string> & key, tmatrix<T> & value);
        template <class T> void write_attribute(filesystem::path object_path, map<string, T> & dict);
        void write_attribute(filesystem::path object_path, YAML::Node & other_node);
        template <class T> void read_dataset(filesystem::path dataset_path, tmatrix<T> & data);
        template <class T> void read_attribute(filesystem::path object_path, string key, T & value);
        void get_dataset_size(filesystem::path dataset_path, tmatrix<unsigned long> & shape);
        void clean_attr_file(filesystem::path object_path);
        bool file_exists();
};

inline exdir::exdir(filesystem::path _file_path, bool _overwrite_file){
    file_path = _file_path;
    if(_overwrite_file) 
        filesystem::remove_all(file_path);

    create_folder(file_path, "file");
}


inline exdir & exdir::operator=(const exdir & other)
{
    if (this != &other) { // protect against invalid self-assignment
        file_path = other.file_path;
    }
    return *this;
}

inline exdir::exdir (const exdir & other)
{
    file_path = other.file_path;
}

inline void exdir::create_folder(filesystem::path path, string type){
    filesystem::create_directories(path);
    write_metadata(path, type, VERSION);
}   

template<class T>
void exdir::save_to_file(filesystem::path file_path, T content, ios::openmode mode)
{
    ofstream file;
    file.open(file_path, mode);
    file << content;
    file.close();
}


inline bool exdir::is_type(filesystem::path path, string type){
    try { 
        YAML::Node node = YAML::LoadFile(path / METADATA_FILE);
        return node["exdir"]["type"].as<string>() == type;
    
    } catch (const std::exception& e) {
        return false;
    }
    return false;
}

inline bool exdir::file_exists(){
    return filesystem::exists(file_path);
}

inline void exdir::write_metadata(filesystem::path path, string type, string version){
   
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
    fout.open(path / METADATA_FILE);
    fout << out.c_str();
    fout.close();      
}

inline void exdir::create_group(filesystem::path group_path){
    create_folder(file_path / group_path, "group");
}

template <class T>
inline void exdir::write_dataset(filesystem::path dataset_path, tmatrix<T> & data){
    
    create_folder(file_path / dataset_path, "dataset");
    
    tmatrix<unsigned long> shape = {data.n1, data.n2, data.n3};
    if (shape.val[2] == 1)  
    {
        shape.n1 -= 1;
        if (shape.val[1] == 1)
            shape.n1 -= 1;
    }
    
    npy::SaveArrayAsNumpy(file_path / dataset_path / DATA_FILE, false, shape, data);
}

inline void exdir::clean_attr_file(filesystem::path object_path){
    
    filesystem::path attr_path = file_path / object_path / ATTR_FILE;
    
    ofstream fs;
    fs.open(attr_path, ofstream::out | ofstream::trunc);
    fs.close();
}

inline void exdir::ensure_attr_exists(filesystem::path object_path){
    
    string attr_path = file_path / object_path / ATTR_FILE;
    
    ofstream fs;
    fs.open(attr_path, ofstream::app);
    fs.close();
}

template <class T>
inline void exdir::write_attribute(filesystem::path object_path, string key, T value) {
    
    ensure_attr_exists(object_path);
 
    string attr_path = file_path / object_path / ATTR_FILE;
       
    YAML::Node node = YAML::LoadFile(attr_path);
    node[key] = value;
    
    ofstream fout(attr_path);
    fout << node;
    fout.close();
    
}


template <class T>
inline void exdir::write_attribute(filesystem::path object_path, tmatrix<string> & key, tmatrix<T> & value) {
    
    ensure_attr_exists(object_path);
    
    filesystem::path attr_path = file_path / object_path / ATTR_FILE;

    YAML::Node node = YAML::LoadFile(attr_path);
    
    for(int k = 0; k < key.n1; k++){
        node[key.val[k]] = value.val[k];
    }
   
    ofstream fout(attr_path);
    fout << node;
    fout.close();
}

template <class T>
inline void exdir::write_attribute(filesystem::path object_path, map<string, T> & dict) {
    
    ensure_attr_exists(object_path);
    
    filesystem::path attr_path = file_path / object_path / ATTR_FILE;

    YAML::Node node = YAML::LoadFile(attr_path);
    
    for (auto const& element : dict)
    {
        node[element.first] = element.second;

    }
   
    ofstream fout(attr_path);
    fout << node;
    fout.close();
}

inline void exdir::write_attribute(filesystem::path object_path, YAML::Node & other_node) {
    
    ensure_attr_exists(object_path);
    
    filesystem::path attr_path = file_path / object_path / ATTR_FILE;

    YAML::Node node = YAML::LoadFile(attr_path);

    for (YAML::const_iterator it=other_node.begin();it!=other_node.end();++it) {
        node[it->first.as<string>()] = other_node[it->first.as<string>()];
    }

    save_to_file(attr_path, node, ios::trunc);
}

inline YAML::Node exdir::get_attributes(filesystem::path object_path){
    ensure_attr_exists(object_path);
    filesystem::path attr_path = file_path / object_path / ATTR_FILE;
    return YAML::LoadFile(attr_path);
}

inline void exdir::get_dataset_size(filesystem::path dataset_path, tmatrix<unsigned long> & shape){
    filesystem::path data_path = file_path / dataset_path / DATA_FILE;
    bool fortran_order = false;
    npy::GetArrayShape(data_path, shape, fortran_order);
}


template <class T>
inline void exdir::read_dataset(filesystem::path dataset_path, tmatrix<T> & data){
    
    filesystem::path data_path = file_path / dataset_path / DATA_FILE;

    cout << data_path << endl;
    
    tmatrix<unsigned long> shape(3);
    
    bool fortran_order = false;
    npy::LoadArrayFromNumpy(data_path, shape, fortran_order, data);
}

template <class T>
inline void exdir::read_attribute(filesystem::path object_path, string key, T & value) {
    filesystem::path attr_path = file_path / object_path / ATTR_FILE;
    
    YAML::Node node = YAML::LoadFile(attr_path);
    value = node[key].as<T>();       
}


inline void exdir::read_all_attributes(filesystem::path object_path, YAML::Node & node) {
    ensure_attr_exists(object_path);
    filesystem::path attr_path = file_path / object_path / ATTR_FILE;
    node = YAML::LoadFile(attr_path);
}



#endif
