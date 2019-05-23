/*
  Random number generator.
*/

#include "random-numbers.h"
#include <random>

using namespace std;

random_device r;
seed_seq seed{ r(), r(), r(),  r() };
mt19937 gen{ seed };

uniform_real_distribution<double> uniform(0.0, 1.0);
normal_distribution<double> normal(0.0, 1.0);


double r_unif(){
    return uniform(gen);
}

double r_norm(){
    return normal(gen);
}

double r_norm(double average, double sigma){
    return sigma * normal(gen) + average;
}


