
#include "diagnostics.h"
#include "fmatrix.h"
#include "particles-in-mesh.h"
#include "particles.h"
#include "fields.h"
#include "mpi.h"

#include <unordered_map>

#define HIGH 1e99
#define LOW -1e99

using namespace std;

diagnostics::diagnostics(configuration & _config, state_info & _state): config(_config), state(_state){
	k_v  = config.f("p/k_v");
	dt = config.f("time/dt");
	k_phi = config.f("p/k_phi");
	n_factor = config.f("particles/n_factor");
	q = config.f("physical/q");

	n_v_e = config.i("diagnostics/vdist/electrons/n_v");
    n_v_i = config.i("diagnostics/vdist/ions/n_v");

	n_e_iz = 0;

	n_mesh_x = config.i("geometry/n_mesh_x");
	n_mesh_y = config.i("geometry/n_mesh_y");
	
	double dx = config.f("geometry/dx");
	double dy = config.f("geometry/dy");
	x_max = ((double) n_mesh_x - 1);
	y_max = ((double) n_mesh_y - 1) * (dy / dx);

	step_save_vdist = config.i("diagnostics/vdist/save_step");
	step_save_fields = config.i("diagnostics/fields_snapshot/save_step");

    n_steps = config.i("time/n_steps");
	izfield_start_progress = config.f("diagnostics/iz_field/start_progress");

	steps_since_last_izfield_save = 0;

	vdist_e_x = fmatrix::zeros(n_v_e);
	vdist_e_y = fmatrix::zeros(n_v_e);
    vdist_e_global_x = fmatrix::zeros(n_v_e);
	vdist_e_global_y = fmatrix::zeros(n_v_e);
    vdist_i_x = fmatrix::zeros(n_v_i);
	vdist_i_y = fmatrix::zeros(n_v_i);
    vdist_i_global_x = fmatrix::zeros(n_v_i);
	vdist_i_global_y = fmatrix::zeros(n_v_i);


	top_dist_e = fmatrix::zeros(n_v_e);
	top_dist_i = fmatrix::zeros(n_v_i);
	top_dist_e_global = fmatrix::zeros(n_v_e);
	top_dist_i_global = fmatrix::zeros(n_v_i);

	rhs_dist_e = fmatrix::zeros(n_v_e);
	rhs_dist_i = fmatrix::zeros(n_v_i);
	rhs_dist_e_global = fmatrix::zeros(n_v_e);
	rhs_dist_i_global = fmatrix::zeros(n_v_i);

	vmin = fmatrix::zeros(6);
	vmax = fmatrix::zeros(6);

	kfield_e = fmatrix::zeros(n_mesh_x, n_mesh_y);
	kfield_i = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_e_x = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_e_y = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_i_x = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_i_y = fmatrix::zeros(n_mesh_x, n_mesh_y);
	izfield = fmatrix::zeros(n_mesh_x, n_mesh_y);

	kfield_e_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	kfield_i_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_e_x_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_e_y_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_i_x_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	ufield_i_y_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	izfield_global = fmatrix::zeros(n_mesh_x, n_mesh_y);

	wmesh_e = fmatrix::zeros(n_mesh_x, n_mesh_y);
	wmesh_i = fmatrix::zeros(n_mesh_x, n_mesh_y);
	wmesh_e_global = fmatrix::zeros(n_mesh_x, n_mesh_y);
	wmesh_i_global = fmatrix::zeros(n_mesh_x, n_mesh_y);

	p_select = fmatrix::zeros(100000, 6);
	p_e_removed = fmatrix::zeros(100000, 6);
	p_i_removed = fmatrix::zeros(100000,6);
	n_removed_e = n_removed_i = 0;

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

void diagnostics::dist(fmatrix & p, int & n_active, int col, double v_0, double v_1,  fmatrix & dmesh){

	double v_p, v_p_norm, s;
	int left_index;

	double dv = (v_1 - v_0) / (double) (dmesh.n1 - 1);

	for(int i = 0; i < n_active; i++){

	    v_p = p.val[i * 6 + col];

		if(v_p < v_0 || v_p > v_1) continue;
	
		v_p_norm = (v_p - v_0) / dv;
		left_index = floor(v_p_norm);
		s = v_p_norm - (double) left_index;
		
		dmesh.val[left_index] += 1 - s;
		dmesh.val[left_index + 1] += s;
	}
}

double diagnostics::av_div(fmatrix & p, int & n_active){
	double cos_theta = 0, vx = 0, vy = 0;
	for(int i = 0; i < n_active; i++){
		vx = p.val[i * 6 + 3];
		vy = p.val[i * 6 + 4];

		cos_theta += vx / sqrt(vx*vx + vy*vy);
	}
	return cos_theta / (double) n_active;
}


void diagnostics::update_velocity_distributions(fmatrix & p_e, fmatrix & p_i){
	if(state.step % step_save_vdist != 0) return;
	
	vdist_i_x.set_zero();
	vdist_i_y.set_zero();
	vdist_e_x.set_zero();
	vdist_e_y.set_zero();

	dist(p_i, state.n_active_i, 3, k_v * vlim_i.val[0], k_v * vlim_i.val[1], vdist_i_x);
    dist(p_i, state.n_active_i, 4, k_v * vlim_i.val[2], k_v * vlim_i.val[3], vdist_i_y);
    dist(p_e, state.n_active_e, 3, k_v * vlim_e.val[0], k_v * vlim_e.val[1], vdist_e_x);
    dist(p_e, state.n_active_e, 4, k_v * vlim_e.val[2], k_v * vlim_e.val[3], vdist_e_y);
}

void diagnostics::update_boundary_current_distributions(){

	int n_active_select = 0;

	// select particles at the top boundary
	vmin = {x_min, y_max, LOW, LOW, LOW, LOW};
	vmax = {x_max, HIGH, HIGH, HIGH, HIGH, HIGH};

	particle_operations::select_particles(p_i_removed, n_removed_i, p_select, n_active_select, vmin, vmax);
	dist(p_select, n_active_select, 0, 0, x_max, top_dist_i);

	particle_operations::select_particles(p_e_removed, n_removed_e, p_select, n_active_select, vmin, vmax);
	dist(p_select, n_active_select, 0, 0, x_max, top_dist_e);

	// add axial current
	// dont forget to divide these distributions by n*dt to give current 

	// select particles at the rhs boundary
	vmin = {x_max, y_min, LOW, LOW, LOW, LOW};
	vmax = {HIGH, y_max, HIGH, HIGH, HIGH, HIGH};

	particle_operations::select_particles(p_i_removed, n_removed_i, p_select, n_active_select, vmin, vmax);
	dist(p_select, n_active_select, 1, 0, y_max, rhs_dist_i);

	particle_operations::select_particles(p_e_removed, n_removed_e, p_select, n_active_select, vmin, vmax);
	dist(p_select, n_active_select, 1, 0, y_max, rhs_dist_e);
}

void diagnostics::reduce_distributions(){

	//velocity distributions
	MPI_Reduce(vdist_i_x.val, vdist_i_global_x.val, n_v_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(vdist_i_y.val, vdist_i_global_y.val, n_v_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(vdist_e_x.val, vdist_e_global_x.val, n_v_e, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(vdist_e_y.val, vdist_e_global_y.val, n_v_e, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	//boundary current distributions
	MPI_Reduce(top_dist_i.val, top_dist_i_global.val, n_v_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(rhs_dist_i.val, rhs_dist_i_global.val, n_v_i, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(top_dist_e.val, top_dist_e_global.val, n_v_e, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rhs_dist_e.val, rhs_dist_e_global.val, n_v_e, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	boundary_dist_set_zero();
}



void diagnostics::update_series() {
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


void diagnostics::update_internal_wmesh(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, imatrix & lpos_e, imatrix & lpos_i, bool force){
		
		if(!(force||state.step % step_save_fields == 0)) return;
		
		wmesh_e_global.set_zero();
		wmesh_i_global.set_zero();

	    pic_operations::weight(p_i, state.n_active_i, wmesh_i, mesh, lpos_i);
        pic_operations::weight(p_e, state.n_active_e, wmesh_e, mesh, lpos_e);
        
        MPI_Reduce(wmesh_i.val, wmesh_i_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(wmesh_e.val, wmesh_e_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void diagnostics::set_internal_wmesh(fmatrix & _wmesh_e_global, fmatrix & _wmesh_i_global){
	for (int i = 0; i < n_mesh_x * n_mesh_y; i++){
		wmesh_e_global.val[i]= _wmesh_e_global.val[i];
		wmesh_i_global.val[i]= _wmesh_i_global.val[i];
	}
}


void diagnostics::update_ufield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, imatrix & lpos_e, imatrix & lpos_i, bool force){

	if(!(force||state.step % step_save_fields == 0)) return;

	ufield_e_x_global.set_zero();
	ufield_e_y_global.set_zero();
	ufield_i_x_global.set_zero();
	ufield_i_y_global.set_zero();

	pic_operations::ufield(p_e, state.n_active_e, ufield_e_x, ufield_e_y, mesh, lpos_e); 
	MPI_Reduce(ufield_e_x.val, ufield_e_x_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(ufield_e_y.val, ufield_e_y_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for (int i = 0; i < n_mesh_x * n_mesh_y; i++){
		if(wmesh_e_global.val[i] > 0){
			ufield_e_x_global.val[i] = ufield_e_x_global.val[i] / wmesh_e_global.val[i];
			ufield_e_y_global.val[i] = ufield_e_y_global.val[i] / wmesh_e_global.val[i];
		}
		else {
			ufield_e_x_global.val[i] = 0;
			ufield_e_y_global.val[i] = 0;
		}
	}


	pic_operations::ufield(p_i, state.n_active_i, ufield_i_x, ufield_i_y, mesh, lpos_i); 
	MPI_Reduce(ufield_i_x.val, ufield_i_x_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(ufield_i_y.val, ufield_i_y_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for (int i = 0; i < n_mesh_x * n_mesh_y; i++){
		if(wmesh_i_global.val[i] > 0){
			ufield_i_x_global.val[i] = ufield_i_x_global.val[i] / wmesh_i_global.val[i];
			ufield_i_y_global.val[i] = ufield_i_y_global.val[i] / wmesh_i_global.val[i];
		}
		else {
			ufield_i_x_global.val[i] = 0;
			ufield_i_y_global.val[i] = 0;
		}
	}
}


void diagnostics::update_kfield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, imatrix & lpos_e, imatrix & lpos_i, bool force){

	if(!(force||state.step % step_save_fields == 0)) return;

	kfield_e_global.set_zero();
	kfield_i_global.set_zero();

	pic_operations::kfield(p_e, state.n_active_e, kfield_e, mesh, lpos_e); 
	MPI_Reduce(kfield_e.val, kfield_e_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	pic_operations::kfield(p_i, state.n_active_i, kfield_i, mesh, lpos_i); 
	MPI_Reduce(kfield_i.val, kfield_i_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	for (int i = 0; i < n_mesh_x * n_mesh_y; i++){
		if(wmesh_e_global.val[i] > 0)
			kfield_e_global.val[i] = kfield_e_global.val[i] / wmesh_e_global.val[i];
		else 
			kfield_e_global.val[i] = 0;
	}

	for (int i = 0; i < n_mesh_x * n_mesh_y; i++){
		if(wmesh_i_global.val[i] > 0)
			kfield_i_global.val[i] = kfield_i_global.val[i] / wmesh_i_global.val[i];
		else 
			kfield_i_global.val[i] = 0;
	}
}

void diagnostics::update_izfield(mesh_set & mesh, fmatrix & p_i, imatrix & lpos_i, int n_iz){
	if(state.step < izfield_start_progress * (double) n_steps) return;
	int n_start =  state.n_active_i - n_iz;
	izfield_global.set_zero();
	pic_operations::weight_n(p_i, n_start, state.n_active_i, izfield, mesh, lpos_i);
	MPI_Reduce(izfield.val, izfield_global.val, n_mesh_x * n_mesh_y, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void diagnostics::izfield_set_zero(){
	izfield_global.set_zero();
	izfield.set_zero();
}

void diagnostics::boundary_dist_set_zero(){
	top_dist_e.set_zero();
	top_dist_i.set_zero();
	rhs_dist_e.set_zero();
	rhs_dist_i.set_zero();
}

void diagnostics::update_all(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, imatrix & lpos_e, imatrix & lpos_i){
	update_series();
	update_velocity_distributions(p_e, p_i);
	update_boundary_current_distributions();
	update_internal_wmesh(mesh, p_e, p_i, lpos_e, lpos_i);
	update_ufield(mesh, p_e, p_i, lpos_e, lpos_i);
	update_kfield(mesh, p_e, p_i, lpos_e, lpos_i);
	update_izfield(mesh, p_i, lpos_i, n_e_iz);
}