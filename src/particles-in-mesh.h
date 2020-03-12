#ifndef PARTICLES_IN_MESH_H
#define PARTICLES_IN_MESH_H

#include "fmatrix.h"

void weight(fmatrix & p, int & n_active, fmatrix & wmesh, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos, const double & a_x, const double & a_y, const double & dx, const double & dy);
void electric_field_at_particles(fmatrix & efield_at_particles_x, fmatrix & efield_at_particles_y, fmatrix & efield_x, fmatrix & efield_y, fmatrix & p, const int n_active, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos, const double & a_x, const double & a_y, const double & dx, const double & dy);
void energy_field(fmatrix & kefield, fmatrix & p, int & n_active, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & wmesh, imatrix & lpos, double mass, const double & dx, const double & dt, const double & q);
void flux_field(fmatrix & ffield_x, fmatrix & ffield_y, fmatrix & p, int & n_active, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos);
double field_at_position(fmatrix & field, fmatrix & mesh_x, fmatrix & mesh_y, double x, double y, int lpos_x, int lpos_y, const double & a_x, const double & a_y, const double & dx, const double & dy);

#endif
