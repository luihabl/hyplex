#ifndef PARTICLES_IN_MESH_H
#define PARTICLES_IN_MESH_H

#include "fmatrix.h"

void weight(fmatrix & p, int & n_active, fmatrix & wmesh, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos);
void electric_field_at_particles(fmatrix & efield_at_particles_x, fmatrix & efield_at_particles_y, fmatrix & efield_x, fmatrix & efield_y, fmatrix & p, const int n_active, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos);
void energy_field(fmatrix & kefield, fmatrix & p, int & n_active, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & wmesh, imatrix & lpos, double mass);

#endif
