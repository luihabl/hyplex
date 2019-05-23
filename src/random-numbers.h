//
//  random-numbers.h
//  pic-benchmark
//
//  Created by Lui Txai Calvoso Habl on 29/01/2019.
//  Copyright © 2019 Lui Txai Calvoso Habl. All rights reserved.
//

#ifndef RANDOM_NUMBERS_H
#define RANDOM_NUMBERS_H

#include <random>

using namespace std;

extern random_device r;
extern seed_seq seed;
extern mt19937 gen;

extern uniform_real_distribution<double> uniform;
extern normal_distribution<double> normal;

double r_unif();
double r_norm();
double r_norm(double average, double sigma);


#endif /* RANDOM_NUMBERS_H */
