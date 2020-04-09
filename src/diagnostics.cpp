
#include "diagnostics.h"
#include "fmatrix.h"


diagnostics::diagnostics(){}

void diagnostics::velocity_distribution(fmatrix & p, int & n_active, int & vcol, int n_v, double v_0, double v_1,  fmatrix & dmesh){

	dmesh.set_zero();

	double v_p, v_p_norm, s;
	int left_index;

	double k = (v_1 - v_0) / (double) (n_v - 1);

	for(int i = 0; i < n_active; i++){

	    v_p = p.val[i * 6 + vcol];

		if(v_p < v_0 || v_p > v_1) continue;
	
		v_p_norm = (v_p - v_0) / k;
		left_index = floor(v_p_norm);
		s = v_p_norm - (double) left_index;
		
		dmesh.val[left_index] = 1 - s;
		dmesh.val[left_index + 1] = s;
	}
}

