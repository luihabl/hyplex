#ifndef PARTICLES_IN_MESH_H
#define PARTICLES_IN_MESH_H

#include "fmatrix.h"
#include "fields.h"
#include "configuration.h"

struct pic_operations{

    double dx, dy, a_x, a_y, dt, q;

    int n_mesh_x, n_mesh_y;

    pic_operations(configuration & config);

    static void weight(fmatrix & p, int & n_active, fmatrix & wmesh, mesh_set & mesh, imatrix & lpos);
    void electric_field_at_particles(fmatrix & efield_at_particles_x, fmatrix & efield_at_particles_y, fmatrix & efield_x, fmatrix & efield_y, fmatrix & p, const int n_active, mesh_set & mesh, imatrix & lpos);
    static void kfield(fmatrix & p, int & n_active, fmatrix & kefield, mesh_set & mesh, imatrix & lpos);
    static void ufield(fmatrix & p, int & n_active, fmatrix & ffield_x, fmatrix & ffield_y, mesh_set & mesh,imatrix & lpos);
    static void weight_n(fmatrix & p, int & n_active_start, int & n_active_end, fmatrix & f, mesh_set & mesh, imatrix & lpos);
    double field_at_position(fmatrix & field, mesh_set & mesh, double x, double y, int lpos_x, int lpos_y);

};



#endif
