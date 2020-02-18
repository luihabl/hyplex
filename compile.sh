#!/bin/bash
# mpicxx -std=c++11 src/*.cpp -Ilib/hypre-2.16.0/src/hypre/include -Ilib -Llib/hypre-2.16.0/src/hypre/lib -lHYPRE -O3 -march=native -Wall -D VERBOSE -o run

machine=$1

if [ "$machine" == "zoidberg" ]
then
    echo "Compiling for zoidberg"
    
    module load openmpi/4.0.2
    module load hdf5-openmpi-4.0.2/1.10.5
    module load hypre/2.11.1

    make machine=zoidberg

elif [ "$machine" == "" ]
then
    echo "Compiling for other machines"

    make machine=others
fi

