#!/bin/bash
mpicxx -std=c++11 src/*.cpp -Ilib/hypre-2.16.0/src/hypre/include -Ilib -Llib/hypre-2.16.0/src/hypre/lib -lHYPRE -O3 -march=native -Wall -D VERBOSE -o run
