
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

fmatrix load_csv(string file_path, char delim, int cols)
{
    string line;
    ifstream file(file_path);
    
    int line_count = 0;
    while (getline(file, line)) line_count++;
    file.clear();
    file.seekg(0, ios::beg);
    
    fmatrix data(line_count, cols); //line_count - 1 because the first row it's assumed to be a header.
    int i = 0;
    while (getline(file, line))
    {
        istringstream s(line);
        string field;
        int j = 0;
        while (getline(s, field, delim))
        {
            if (j > cols - 1) {
                cout << "csv exceed j limit!" << endl;
                break;
            }
            data.val[i * cols + j] = atof(field.c_str());
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

void print_info(int i, fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i, double v_cap, int step_interval)
{
    static high_resolution_clock::time_point t0;
    if(i==0) t0 = high_resolution_clock::now();
    if ((i + 1) % step_interval == 0 || i == 0)
    {
        printf("[%05.2f%%] ", (double) (100.0 * (i + 1) / N_STEPS));
        printf("Step: %-8d ", i + 1);
        printf("Active electrons: %-8d ", n_active_e);
        printf("Active ions: %-8d ", n_active_i);
        printf("Cap. voltage: %.4f V   ", v_cap);
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
    printf("Gas:\t\t\t %s\n", GAS_NAME.c_str());
    printf("P Null (e):\t\t %.4e\n", p_null_e);
    printf("P Null (I):\t\t %.4e\n\n", p_null_i);
}



void average_field(fmatrix & av_field, const fmatrix & field, int step)
{
    // step is actually step = i - (N_STEPS - N_AVERAGE);
    for (size_t j = 0; j < field.n1 * field.n2 * field.n3; j++) {
        av_field.val[j] = ((step - 1) * av_field.val[j] / step) + (field.val[j] / step);
    }
}

void save_state(fmatrix & p_e, int n_active_e, fmatrix & p_i, int n_active_i, fmatrix & phi,  fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, string suffix){
    
    
    // fmatrix p_e_corrected = DX * p_e;
    // for(int i = 0; i < n_active_e; i++){
    //     p_e_corrected.val[i * 6 + 3] = p_e_corrected.val[i * 6 + 3] / DT;
    //     p_e_corrected.val[i * 6 + 4] = p_e_corrected.val[i * 6 + 4] / DT;
    //     p_e_corrected.val[i * 6 + 5] = p_e_corrected.val[i * 6 + 5] / DT;
    // }
    // save_to_csv(p_e_corrected, "p_e" + suffix + ".csv", n_active_e);

    // fmatrix p_i_corrected = DX * p_i;
    // for(int i = 0; i < n_active_i; i++){
    //     p_i_corrected.val[i * 6 + 3] = p_i_corrected.val[i * 6 + 3] / DT;
    //     p_i_corrected.val[i * 6 + 4] = p_i_corrected.val[i * 6 + 4] / DT;
    //     p_i_corrected.val[i * 6 + 5] = p_i_corrected.val[i * 6 + 5] / DT;
    // }
    // save_to_csv(p_i_corrected, "p_i" + suffix + ".csv", n_active_i);
    
    
    fmatrix phi_corrected = phi * (M_EL * pow(DX, 2))/(Q * pow(DT, 2));
    save_to_csv(phi_corrected, "phi" + suffix + ".csv");

    fmatrix dens_e = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_e / vmesh;
    save_to_csv(dens_e, "dens_e" + suffix + ".csv");

    fmatrix dens_i = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_i / vmesh;
    save_to_csv(dens_i, "dens_i" + suffix + ".csv");
    
    verbose_log("State saved");
}

void verbose_log(string message){
    #ifdef VERBOSE
    cout << message << endl;
    #endif
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
