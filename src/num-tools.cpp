#include "num-tools.h"

#include "random-numbers.h"
#include <algorithm>
#include <cmath>


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


void average_field(fmatrix & av_field, const fmatrix & field, int step)
{
    // step is actually step = i - (N_STEPS - N_AVERAGE);
    for (size_t j = 0; j < field.n1 * field.n2 * field.n3; j++) {
        av_field.val[j] = ((step - 1) * av_field.val[j] / step) + (field.val[j] / step);
    }
}
