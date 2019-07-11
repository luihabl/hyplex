
#include "util.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <chrono>

#include "fmatrix.h"
#include "fmath.h"
#include "mcc.h"
#include "random-numbers.h"

using namespace std;
using namespace std::chrono;

double interp(const fmatrix & data, double x)
{
	if (x <= data.val[0 * 2 + 0])
		return data.val[0 * 2 + 1];

	else if (x >= data.val[(data.n1 - 1) * 2 + 0])
		return data.val[(data.n1 - 1) * 2 + 1];

	else
	{
		int left_index = 0;
		while (x > data.val[(left_index + 1) * 2 + 0]) left_index++;

		double x0 = data.val[left_index * 2 + 0];
		double y0 = data.val[left_index * 2 + 1];
		double x1 = data.val[(left_index + 1) * 2 + 0];
		double y1 = data.val[(left_index + 1) * 2 + 1];

		return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
	}
}

fmatrix interp(const fmatrix & data, fmatrix & x)
{
	fmatrix interp_results(x.n1);
	for (size_t i = 0; i < interp_results.n1; i++)
		interp_results.val[i] = interp(data, x.val[i]);

	return interp_results;
}


// Identify what is going wrong here?
// imatrix sample_from_sequence_swap(int sample_size, int range){

//     imatrix sequence(range);
// 	sequence.set_sequence();
//     imatrix samples(sample_size);

//     int n_storage = range;
//     int num = 0;

//     for (int i = 0; i < sample_size; i++)
//     {
//         num = r_unif() * n_storage;
//         samples.val[i] = sequence.val[num];
//         sequence.val[num] = sequence.val[range - i - 1];
//         n_storage -= 1;
//     }

//     return samples;
// }

imatrix sample_from_sequence_shuffle(int sample_size, int range)
{
	imatrix shuffled_sequence(range);
	shuffled_sequence.set_sequence();

    shuffle(&shuffled_sequence.val[0], &shuffled_sequence.val[range], gen);
	imatrix samples(sample_size);

	for (int i = 0; i < sample_size; i++)
		samples.val[i] = shuffled_sequence.val[i];

	return samples;
}

bool is_in_set(imatrix & seq, int last_index) {
	for (int j = 0; j < last_index; j++)
		if (seq.val[j] == seq.val[last_index]) return true;
	return false;
}

imatrix sample_from_sequence_naive(int sample_size, int range)
{
	imatrix samples(sample_size);
    
	bool repeated = false;
	for (int i = 0; i < sample_size; i++)
	{
		do {
			samples.val[i] = (int) floor(r_unif() * (double) range);
			repeated = is_in_set(samples, i);
		} while (repeated && i > 0);
	}
	return samples;
}

fmatrix load_csv(string file_path, char delim)
{
    string line;
    ifstream file(file_path);
    
    int line_count = 0;
    while (getline(file, line)) line_count++;
    file.clear();
    file.seekg(0, ios::beg);
    
    fmatrix data(line_count, 2); //line_count - 1 because the first row it's assumed to be a header.
    size_t i = 0;
    while (getline(file, line))
    {
        istringstream s(line);
        string field;
        size_t j = 0;
        while (getline(s, field, delim))
        {
            if (j > 1) {
                cout << "csv exceed j limit!" << endl;
                break;
            }
            
            data.val[i * 2 + j] = atof(field.c_str());
            j++;
        }
        i++;
    }
    
    if (!file.eof()) {
        cout << "Could not read file " << endl;
    }
    
    return data;
}

void print_dsmc_info(int i, int n_active_n, int step_interval, int n_steps){
    static high_resolution_clock::time_point t0;
    if(i==0) t0 = high_resolution_clock::now();
    if ((i + 1) % step_interval == 0 || i == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (i + 1) / n_steps));
        printf("Step: %-8d ", i + 1);
        printf("Active neutrals: %-8d ", n_active_n);
        printf("Loop time: %.2f ms ", (double) duration_cast<microseconds>(high_resolution_clock::now() - t0).count() / (1e3 * step_interval));
        printf("\n");
        t0 = high_resolution_clock::now();
    }
}

void print_info(int i, fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i,  int step_interval)
{
    static high_resolution_clock::time_point t0;
    if(i==0) t0 = high_resolution_clock::now();
    if ((i + 1) % step_interval == 0 || i == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (i + 1) / N_STEPS));
        printf("Step: %-8d ", i + 1);
        printf("Active electrons: %-8d ", n_active_e);
        printf("Active ions: %-8d ", n_active_i);
        printf("Loop time: %.2f ms ", (double) duration_cast<microseconds>(high_resolution_clock::now() - t0).count() / (1e3 * step_interval));
        printf("\n");
        t0 = high_resolution_clock::now();
    }
}

void print_initial_info(double p_null_e, double p_null_i)
{
    verbose_log("\n ---- Simulation parameters ----");
    printf("Grid size:\t\t (%d, %d)\n", N_MESH_X, N_MESH_Y);
    printf("Number of steps:\t %d\n", N_STEPS);
    printf("P Null (e):\t\t %.4e\n", p_null_e);
    printf("P Null (I):\t\t %.4e\n\n", p_null_i);
}



void average_field(fmatrix & av_field, const fmatrix & field, int i)
{
    if(i > (N_STEPS - N_AVERAGE)) {
        double step = (double) i - (N_STEPS - N_AVERAGE);
        for (size_t j = 0; j < field.n1 * field.n2 * field.n3; j++) {
            av_field.val[j] = ((step - 1) * av_field.val[j] / step) + (field.val[j] / step);
        }
    }
}


int clamp(int low, int hi, int val){
    if (val < low) {return low;}
    else if (val > hi) {return hi;}
    else {return val;}
}

void verbose_log(string message){
    #ifdef VERBOSE
    cout << message << endl;
    #endif
}
