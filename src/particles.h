#ifndef PARTICLES_H
#define PARTICLES_H

#include "fmatrix.h"
#include "configuration.h"

int add_particle_copy(fmatrix & p, int & n_active, imatrix & lpos, const int & i);
void add_maxwellian_particles(fmatrix & p, int & n_active, const double temperature, const double mass, const size_t n_add, const double & dx, const double & dy, const double & dt, const int & n_mesh_x,  const int & n_mesh_y, const double & q);
void add_maxwellian_particle_at_position(fmatrix & p, int & n_active, imatrix & lpos, const double temperature, const double mass, double x_pos, double y_pos, int lpos_x, int lpos_y, const double & dt, const double & dx, const double & q);
void add_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub, const double & dx, const double & dy, const double & dt, const int & n_thruster, const double & q);
void add_maxwellian_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub, const double & dx, const double & dy, const double & dt, const int & n_thruster, const double & q);
double balanced_injection(double old_n_inj, double rate_constant, fmatrix & wmesh_i, fmatrix & wmesh_e, int ill, int jll, int iur, int jur);
double pulsed_injection(double k_inj, double v_sb, double v_rf, double temp_e, double omega_i, int i);
double square_injection(double alpha, double freq, double dt, double duty_cycle, int i, double n_factor, double i_i, double q);

void boundaries_ob_count(fmatrix & p, int & n_active, imatrix & lpos, int & n_removed_ob, int & n_removed_thr, configuration & config);
void boundaries_e(fmatrix & p, int & n_active, imatrix & lpos, int n_out_i, configuration & config);
void boundaries_e_cap(fmatrix & p, int & n_active, imatrix & lpos, int & n_out_e, double v_cap, fmatrix & phi, fmatrix & mesh_x, fmatrix & mesh_y, configuration & config);
void boundaries_n(fmatrix & p, int & n_active, imatrix & lpos, configuration & config);
double cap_voltage(double voltage, int n_out_e, int n_out_i, double n_factor, double q, double c_cap);
void remove_particle(fmatrix & p, int & n_active, int i, imatrix & lpos);
void reflect_particle(fmatrix & p, int & n_active, int i, double x, double y, double vx, double vy, int n_mesh_x,
						 int n_mesh_y, double dx, double dy);
double find_e_crit(int n_out_i, imatrix & out, int n_out, fmatrix & p, int n_active);
void find_phi_at_particles(fmatrix & phi_at_patricles, fmatrix & phi, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & out, 
						   int n_out, fmatrix & p, int n_active, imatrix & lpos, const double & a_x, const double & a_y, const double & dx, const double & dy, const double & k_phi);
void move_e(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
void move_i(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y, const double & alpha, const double & k_sub);
void move_n(fmatrix & p, int & n_active, double k_sub);

#endif // !PARTICLES_H
