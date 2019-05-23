#ifndef PARTICLES_H
#define PARTICLES_H

#include "fmatrix.h"

void remove_particle(fmatrix & p, int & n_active, int position, imatrix & lpos);
void add_maxwellian_particles(fmatrix & p, int & n_active, const double temperature, const double mass, const size_t n_add);
void add_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const size_t n_add);
void add_maxwellian_particle_at_position(fmatrix & p, int & n_active, imatrix & lpos, const double temperature, const double mass, double x_pos, double y_pos, int lpos_x, int lpos_y);
void boundaries(fmatrix & p, int & n_active, imatrix & lpos);
void move_e(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
void move_i(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y);
int add_particle_copy(fmatrix & p, int & n_active, imatrix & lpos, const int & i);

#endif // !PARTICLES_H
