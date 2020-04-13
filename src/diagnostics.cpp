
#include "diagnostics.h"
#include "fmatrix.h"

#include <unordered_map>

using namespace std;

diagnostics::diagnostics(configuration & _config, state_info & _state): config(_config), state(_state){
	k_v  = config.f("p/k_v");
	dt = config.f("time/dt");
	k_phi = config.f("p/k_phi");
	n_factor = config.f("particles/n_factor");
	q = config.f("physical/q");

	int n_steps = config.f("time/n_steps");
	series["time"]              = fmatrix::zeros(n_steps);
    series["v_cap"]             = fmatrix::zeros(n_steps);
    series["n_active_e"]        = fmatrix::zeros(n_steps);
    series["n_active_i"]        = fmatrix::zeros(n_steps);
    series["I_out_ob_e"]        = fmatrix::zeros(n_steps);
    series["I_out_ob_i"]        = fmatrix::zeros(n_steps);
    series["I_out_thr_e"]       = fmatrix::zeros(n_steps);
    series["I_out_thr_i"]       = fmatrix::zeros(n_steps);
    series["I_in_thr_e"]        = fmatrix::zeros(n_steps);
    series["I_in_thr_i"]        = fmatrix::zeros(n_steps);
	n_points_series = 0;

	series_measure_step = config.i("diagnostics/series/measure_step");
}

void diagnostics::velocity_distribution(fmatrix & p, int & n_active, int vcol, int n_v, double v_0, double v_1,  fmatrix & dmesh){

	dmesh.set_zero();

	double v_p, v_p_norm, s;
	int left_index;

	v_0 *= k_v;
	v_1 *= k_v;

	double dv = (v_1 - v_0) / (double) (n_v - 1);

	for(int i = 0; i < n_active; i++){

	    v_p = p.val[i * 6 + vcol];

		if(v_p < v_0 || v_p > v_1) continue;
	
		v_p_norm = (v_p - v_0) / dv;
		left_index = floor(v_p_norm);
		s = v_p_norm - (double) left_index;
		
		dmesh.val[left_index] += 1 - s;
		dmesh.val[left_index + 1] += s;
	}
}

void diagnostics::update_series(double n_inj_el, double n_inj_i) {
	if(state.step % series_measure_step == 0){
		series["time"].val[n_points_series] = state.step * dt;
		series["v_cap"].val[n_points_series] = state.phi_zero / k_phi;
		series["n_active_i"].val[n_points_series] = state.n_active_i;
		series["n_active_e"].val[n_points_series] = state.n_active_e;
	
		series["I_out_ob_e"].val[n_points_series] = state.n_out_ob_e * n_factor * q / dt;
		series["I_out_ob_i"].val[n_points_series] = state.n_out_ob_i * n_factor * q / dt;
		series["I_out_thr_e"].val[n_points_series] = state.n_out_thr_e * n_factor * q / dt;
		series["I_out_thr_i"].val[n_points_series] = state.n_out_thr_i * n_factor * q / dt;
		series["I_in_thr_e"].val[n_points_series] = n_inj_el * n_factor * q / dt;
		series["I_in_thr_i"].val[n_points_series] = n_inj_i * n_factor * q / dt;
		n_points_series += 1;
    }
}

