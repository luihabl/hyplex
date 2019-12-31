#ifndef H5IO_H
#define H5IO_H

#include "H5Cpp.h"
#include "h5io.h"

using namespace H5;

DataSet create_dataset(H5File & file, fmatrix & data, string name, int rank);
DataSet edit_dataset(H5File & file, fmatrix & data, string name);

template <class T>
void write_attribute(T & obj, string name, double val){
    
    hsize_t dims[1] = { 1 };
    DataSpace attr_dataspace = DataSpace (1,  dims);
    
    Attribute attribute = obj.createAttribute( name, PredType::NATIVE_DOUBLE, attr_dataspace);
    
    double attr_data[1] = { val };
    attribute.write( PredType::NATIVE_DOUBLE, attr_data);
}

template <class T>
void write_attribute(T & obj, string name, int val){
    
    hsize_t dims[1] = { 1 };
    DataSpace attr_dataspace = DataSpace (1,  dims);
    
    Attribute attribute = obj.createAttribute( name, PredType::NATIVE_INT, attr_dataspace);
    
    int attr_data[1] = { val };
    attribute.write( PredType::NATIVE_INT, attr_data);
}

template <class T, class G>
void read_attribute(T & obj, string name, G & val){
    Attribute attribute = obj.openAttribute(name);
    DataType type = attribute.getDataType();
    attribute.read(type, &val);
}

#endif

