#include <iostream>
#include <string>
#include "H5Cpp.h"
#include "fmatrix.h" 

using namespace H5;
using namespace std;


DataSet create_dataset(H5File & file, fmatrix & data, string name, int rank){
   
    hsize_t dims[3];
    dims[0] = data.n1;
	dims[1] = data.n2;
    dims[2] = data.n3;
	
    DataSpace dataspace(rank, dims);
    DataSet dataset = file.createDataSet(name, PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(data.val, PredType::NATIVE_DOUBLE);

    return dataset;
}

DataSet edit_dataset(H5File & file, fmatrix & data, string name){
    DataSet dataset = file.openDataSet(name);
    dataset.write(data.val, PredType::NATIVE_DOUBLE);

    return dataset;
}
