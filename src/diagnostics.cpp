
#include "diagnostics.h"
#include "fmatrix.h"
#include "particles-in-mesh.h"
#include "fields.h"
#include "mpi.h"

#include <unordered_map>

using namespace std;

diagnostics::diagnostics(configuration & _config, state_info & _state): config(_config), state(_state){
	k_v  = config.f("p/k_v");
	dt = config.f("time/dt");
	k_phi = config.f("p/k_phi");
	n_factor = config.f("particles/n_factor");
	q = config.f("physical/q");

	n_v_e = config.i("diagnostics/vdist/electrons/n_v");
    n_v_i = config.i("diagnostics/vdist/ions/n_v");

	n_mesh_x = config.i("geometry/n_mesh_x");
	n_mesh_y = config.i("geometry/n_mesh_y");

	step_save_vdist = config.i("diagnostics/vdist/save_step");
	step_save_fields = config.i("diagnostics/fields_snapshot/save_step");

	dist_e_x = fmatrix::zeros(n_v_e);
	dist_e_y = fmatrix::zeros(n_v_e);
    dist_e_global_x = fmatrix::zeros(n_v_e);
	dist_e_global_y = fmatrix::zeros(n_v_e);
    dist_i_x = fmatrix::zeros(n_v_i);
	dist_i_y = fmatrix::zeros(n_v_i);
    dist_i_global_x = fmatrix::zeros(n_v_i);
	dist_i_global_y = fmatrix::zeros(n_v_i);

	pfield_e = fmatrix::zeros(n_mesh_x, n_mesh_y);
	pfield_i = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_e_x = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_e_y = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_i_x = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_i_y = fmatrix::zeros(n_mesh_x, n_mesh_y);
	izfield = fmatrix::zeros(n_mesh_x, n_mesh_y);

	pfield_e_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	pfield_i_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_e_x_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_e_y_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_i_x_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ffield_i_y_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	izfield_global = fmatrix::zeros(n_mesh_x, n_mesh_y);

    vlim_e = fmatrix({config.fs("diagnostics/vdist/electrons/vlim_x").val[0], config.fs("diagnostics/vdist/electrons/vlim_x").val[1],
              config.fs("diagnostics/vdist/electrons/vlim_y").val[0], config.fs("diagnostics/vdist/electrons/vlim_y").val[1]});

    vlim_i = fmatrix({config.fs("diagnostics/vdist/ions/vlim_x").val[0], config.fs("diagnostics/vdist/ions/vlim_x").val[1],
              config.fs("diagnostics/vdist/ions/vlim_y").val[0], config.fs("diagnostics/vdist/ions/vlim_y").val[1]});
    
    gseries_keys = tmatrix<string>({"time", "v_cap"});
        
    lseries_keys = tmatrix<string>({"n_active_e", "n_active_i", "I_out_ob_e", "I_out_ob_i",
                                    "I_out_thr_e", "I_out_thr_i", "I_in_thr_i", "I_in_thr_e"});
    
	series_measure_step = config.i("diagnostics/series/measure_step");
    
    initialize_series();
}

void diagnostics::initialize_series(){
    
    double n_steps = config.f("time/n_steps");
    double measure_step = config.f("diagnostics/series/measure_step");
    series_size = ceil(n_steps / measure_step);
    
    for(int i = 0; i < (int) lseries_keys.n1; i++)
        lseries[lseries_keys.val[i]] = fmatrix::zeros(series_size);
    
    for(int i = 0; i < (int) gseries_keys.n1; i++)
        gseries[gseries_keys.val[i]] = fmatrix::zeros(series_size);
    
    tmp_array = fmatrix::zeros(series_size);
    
    n_points_series = 0;
}

void diagnostics::velocity_distribution(fmatrix & p, int & n_active, int vcol, double v_0, double v_1,  fmatrix & dmesh){

	dmesh.set_zero();

	double v_p, v_p_norm, s;
	int left_index;

	v_0 *= k_v;
	v_1 *= k_v;

	double dv = (v_1 - v_0) / (double) (dmesh.n1 - 1);

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

void diagnostics::update_distributions(fmatrix & p_e, fmatrix & p_i){
	if(state.step % step_save_vdist != 0) return;
	
	velocity_distribution(p_i, state.n_active_i, 3, vlim_i.val[0], vlim_i.val[1], dist_i_x);
    MPI_Reduce(dist_i_x.val, dist_i_global_x.val, n_v_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    velocity_distribution(p_i, state.n_active_i, 4, vlim_i.val[2], vlim_i.val[3], dist_i_y);
    MPI_Reduce(dist_i_y.val, dist_i_global_y.val, n_v_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    velocity_distribution(p_e, state.n_active_e, 3, vlim_e.val[0], vlim_e.val[1], dist_e_x);
    MPI_Reduce(dist_e_x.val, dist_e_global_x.val, n_v_e, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    velocity_distribution(p_e, state.n_active_e, 4, vlim_e.val[2], vlim_e.val[3], dist_e_y);
    MPI_Reduce(dist_e_y.val, dist_e_global_y.val, n_v_e, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}


void diagnostics::update_series(double n_inj_el, double n_inj_i) {
	if(state.step % series_measure_step == 0){
        
		gseries["time"].val[n_points_series] = state.step * dt;
		gseries["v_cap"].val[n_points_series] = state.phi_zero / k_phi;
        
        lseries["I_in_thr_e"].val[n_points_series] = state.n_in_thr_e * n_factor * q / dt;
        lseries["I_in_thr_i"].val[n_points_series] = state.n_in_thr_i * n_factor * q / dt;
		lseries["n_active_i"].val[n_points_series] = state.n_active_i;
		lseries["n_active_e"].val[n_points_series] = state.n_active_e;
		lseries["I_out_ob_e"].val[n_points_series] = state.n_out_ob_e * n_factor * q / dt;
		lseries["I_out_ob_i"].val[n_points_series] = state.n_out_ob_i * n_factor * q / dt;
		lseries["I_out_thr_e"].val[n_points_series] = state.n_out_thr_e * n_factor * q / dt;
		lseries["I_out_thr_i"].val[n_points_series] = state.n_out_thr_i * n_factor * q / dt;

		n_points_series += 1;
    }
}



void diagnostics::update_ffield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, fmatrix & wmesh_e_global, fmatrix & wmesh_i_global, imatrix & lpos_e, imatrix & lpos_i){

	if(state.step % step_save_fields != 0) return;

	pic_operations::flux_field(p_e, state.n_active_e, ffield_e_x, ffield_e_y, mesh, lpos_e); 
	MPI_Reduce(ffield_e_x.val, ffield_e_x_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(ffield_e_y.val, ffield_e_y_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	pic_operations::flux_field(p_i, state.n_active_i, ffield_i_x, ffield_i_y, mesh, lpos_i); 
	MPI_Reduce(ffield_i_x.val, ffield_i_x_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(ffield_i_y.val, ffield_i_y_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}


void diagnostics::update_pfield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, fmatrix & wmesh_e_global, fmatrix & wmesh_i_global, imatrix & lpos_e, imatrix & lpos_i){

	if(state.step % step_save_fields != 0) return;

	pic_operations::pres_field(p_e, state.n_active_e, pfield_e, mesh, lpos_e); 
	MPI_Reduce(pfield_e.val, pfield_e_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	pic_operations::pres_field(p_i, state.n_active_i, pfield_i, mesh, lpos_i); 
	MPI_Reduce(pfield_i.val, pfield_i_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}



void diagnostics::update_izfield(mesh_set & mesh, fmatrix & p_i, imatrix & lpos_i, int n_iz){
	int n_start =  state.n_active_i - n_iz;
	izfield_global.set_zero();
	pic_operations::weight_n(p_i, n_start, state.n_active_i, izfield, mesh, lpos_i);
	MPI_Reduce(izfield.val, izfield_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
