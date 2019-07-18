#ifndef PARTICLES_H
#define PARTICLES_H

#include "fmatrix.h"

int add_particle_copy(fmatrix & p, int & n_active, imatrix & lpos, const int & i);
void add_maxwellian_particles(fmatrix & p, int & n_active, const double temperature, const double mass, const size_t n_add);
void add_maxwellian_particle_at_position(fmatrix & p, int & n_active, imatrix & lpos, const double temperature, const double mass, double x_pos, double y_pos, int lpos_x, int lpos_y);
void add_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub=1);
double balanced_injection(double old_n_inj, double rate_constant, fmatrix & wmesh_i, fmatrix & wmesh_e, int ill, int jll, int iur, int jur);

void boundaries_i(fmatrix & p, int & n_active, imatrix & lpos, int & n_removed);
void boundaries_e(fmatrix & p, int & n_active, imatrix & lpos, int n_out_i);
void boundaries_n(fmatrix & p, int & n_active, imatrix & lpos);
void remove_particle(fmatrix & p, int & n_active, int i, imatrix & lpos);
void reflect_particle(fmatrix & p, int & n_active, int i, double x, double y, double vx, double vy);
double find_e_crit(int n_out_i, imatrix & out, int n_out, fmatrix & p, int n_active);

void move_e(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
void move_i(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
void move_n(fmatrix & p, int & n_active, double k_sub);

#endif // !PARTICLES_H
