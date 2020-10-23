#ifndef PARTICLES_H
#define PARTICLES_H

#include "fmatrix.h"
#include "configuration.h"
#include "fields.h"
#include "particles-in-mesh.h"


struct particle_operations{

	double dx, dy, y_thruster, dt, q, k_inj, v_sb, v_rf, m_el, freq, duty_cycle, n_factor, k_sub, i_i, alpha, omega_i,c_cap, temp_e, a_x, a_y, k_phi, pi;
	int n_mesh_x, n_mesh_y, n_thruster;
	pic_operations & pic;

	double x_max;
	double y_max;

	imatrix tbremoved;

	particle_operations(configuration & config, pic_operations & _pic);

	int add_particle_copy(fmatrix & p, int & n_active, imatrix & lpos, const int & i);
	void add_maxwellian_particles(fmatrix & p, int & n_active, const double temperature, const double mass, const size_t n_add);
	void add_maxwellian_particle_at_position(fmatrix & p, int & n_active, imatrix & lpos, const double temperature, const double mass, double x_pos, double y_pos, int lpos_x, int lpos_y);
	int add_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub);
	void add_maxwellian_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub);
	
	double balanced_injection(double old_n_inj, double rate_constant, fmatrix & wmesh_i, fmatrix & wmesh_e, int ill, int jll, int iur, int jur);
	double pulsed_injection(int i);
	double square_injection(int i);

	void boundaries_ob_count(fmatrix & p, int & n_active, imatrix & lpos, int & n_removed_ob, int & n_removed_thr, fmatrix & p_removed, int & n_removed, bool copy_removed);
	void copy_removed_particles(fmatrix & p, fmatrix & p_removed, int & n_remove);
	static void select_particles(fmatrix & p, int & n_active, fmatrix & p_select, int & n_active_select, fmatrix & vmin, fmatrix & vmax);
	void boundaries_e(fmatrix & p, int & n_active, imatrix & lpos, int n_out_i);
	void boundaries_e_cap(fmatrix & p, int & n_active, imatrix & lpos, int & n_out_e, double v_cap, fmatrix & phi, mesh_set & mesh);
	void boundaries_n(fmatrix & p, int & n_active, imatrix & lpos);
	void boundaries_n_pump(fmatrix & p, int & n_active, imatrix & lpos, double pump_prob);
	double cap_voltage(double voltage, int n_out_e, int n_out_i);
	void remove_particle(fmatrix & p, int & n_active, int i, imatrix & lpos);
	void remove_particles(imatrix & tbremoved, int & n_remove, fmatrix & p, int & n_active, imatrix & lpos);
	void add_for_removal(imatrix & tbremoved, int & n_remove, int i);
	void reflect_particle(fmatrix & p, int & n_active, int i, double x, double y, double vx, double vy);
	void reflect_particle_specular(fmatrix & p, int i, double x, double y, int boundary_number);
	void reflect_particle_diffuse(fmatrix & p, int & n_active, int i, double x, double y, int boundary_number);
	double find_e_crit(int n_out_i, imatrix & out, int n_out, fmatrix & p, int n_active);
	void find_phi_at_particles(fmatrix & phi_at_patricles, fmatrix & phi, mesh_set & mesh, imatrix & out, int n_out, fmatrix & p, int n_active, imatrix & lpos);
	
	void move_e(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
	void move_i(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
	void move_n(fmatrix & p, int & n_active, double k_sub);

};

inline
void particle_operations::remove_particle(fmatrix & p, int & n_active, int i, imatrix & lpos)
{
	// for (size_t j = 0; j < 6; j++)
	// {
	//   swap(p.val[i * 6 + j], p.val[(n_active - 1) * 6 + j]);
	// }

	// swap(lpos.val[i * 2 + 0], lpos.val[(n_active - 1) * 2 + 0]);
	// swap(lpos.val[i * 2 + 1], lpos.val[(n_active - 1) * 2 + 1]);

	p.val[i * 6 + 0] = p.val[(n_active - 1) * 6 + 0];
	p.val[i * 6 + 1] = p.val[(n_active - 1) * 6 + 1];
	p.val[i * 6 + 2] = p.val[(n_active - 1) * 6 + 2];
	p.val[i * 6 + 3] = p.val[(n_active - 1) * 6 + 3];
	p.val[i * 6 + 4] = p.val[(n_active - 1) * 6 + 4];
	p.val[i * 6 + 5] = p.val[(n_active - 1) * 6 + 5];

	lpos.val[i * 2 + 0] = lpos.val[(n_active - 1) * 2 + 0];
	lpos.val[i * 2 + 1] = lpos.val[(n_active - 1) * 2 + 1];

	n_active -= 1;
}

inline
void particle_operations::add_for_removal(imatrix & tbremoved, int & n_remove, int i){
	tbremoved.val[n_remove] = i;
    n_remove += 1;
}

inline
void particle_operations::remove_particles(imatrix & tbremoved, int & n_remove, fmatrix & p, int & n_active, imatrix & lpos){
	for (int i = n_remove - 1; i >= 0; i--)
    {
        remove_particle(p, n_active, tbremoved.val[i], lpos);
    }
}

#endif // !PARTICLES_H
