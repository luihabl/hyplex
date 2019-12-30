#ifndef PARTICLES_H
#define PARTICLES_H

#include "fmatrix.h"

int add_particle_copy(fmatrix & p, int & n_active, imatrix & lpos, const int & i);
void add_maxwellian_particles(fmatrix & p, int & n_active, const double temperature, const double mass, const size_t n_add);
void add_maxwellian_particle_at_position(fmatrix & p, int & n_active, imatrix & lpos, const double temperature, const double mass, double x_pos, double y_pos, int lpos_x, int lpos_y);
void add_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub=1);
void add_maxwellian_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub);
double balanced_injection(double old_n_inj, double rate_constant, fmatrix & wmesh_i, fmatrix & wmesh_e, int ill, int jll, int iur, int jur);
double pulsed_injection(double k_inj, double v_sb, double v_rf, double temp_e, double omega_i, int i);
double square_injection(double alpha, double freq, double dt, double duty_cycle, int i);

void boundaries_ob_count(fmatrix & p, int & n_active, imatrix & lpos, int & n_removed);
void boundaries_e(fmatrix & p, int & n_active, imatrix & lpos, int n_out_i);
void boundaries_e_cap(fmatrix & p, int & n_active, imatrix & lpos, int & n_out_e, double v_cap, fmatrix & phi, fmatrix & mesh_x, fmatrix & mesh_y);
void boundaries_n(fmatrix & p, int & n_active, imatrix & lpos);
double cap_voltage(double voltage, int n_out_e, int n_out_i);
void remove_particle(fmatrix & p, int & n_active, int i, imatrix & lpos);
void reflect_particle(fmatrix & p, int & n_active, int i, double x, double y, double vx, double vy);
double find_e_crit(int n_out_i, imatrix & out, int n_out, fmatrix & p, int n_active);
void find_phi_at_particles(fmatrix & phi_at_patricles, fmatrix & phi, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & out, int n_out, fmatrix & p, int n_active, imatrix & lpos);

void move_e(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
void move_i(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
void move_n(fmatrix & p, int & n_active, double k_sub);

#endif // !PARTICLES_H
