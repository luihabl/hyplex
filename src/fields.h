#ifndef FIELDS_H
#define FIELDS_H

#include "fmatrix.h"
#include "configuration.h"

struct mesh_set{


    fmatrix x, y, v;
    int nx, ny, nt;
    double a_x, a_y;
    double dx, dy;


    mesh_set(configuration & config);
    void init_mesh();

    double k1_x(size_t i, size_t j, int ioff=0, int joff=0);
    double k2_x(size_t i, size_t j, int ioff=0, int joff=0);
    double k3_x(size_t i, size_t j, int ioff=0, int joff=0);
    double k1_y(size_t i, size_t j, int ioff=0, int joff=0);
    double k2_y(size_t i, size_t j, int ioff=0, int joff=0);
    double k3_y(size_t i, size_t j, int ioff=0, int joff=0);

};

// struct fields{



// };



double ac_voltage_at_time(size_t i, double dt, double freq_hz, double amplitude, double phase, configuration & config);
void calculate_efield(fmatrix & efield_x, fmatrix & efield_y, fmatrix & phi, fmatrix & w_i, fmatrix & w_e, mesh_set & mesh, imatrix & electrode_mask);
double calculate_phi_zero(double sigma_old, double n_in, double q_cap, double sigma_laplace, fmatrix & phi_poisson, mesh_set mesh, fmatrix & wmesh_e, fmatrix & wmesh_i, imatrix & electrode_mask, configuration & config);
double sigma_from_phi(fmatrix & phi, mesh_set & mesh, fmatrix & wmesh_e, fmatrix & wmesh_i, imatrix & electrode_mask, configuration & config);
double calculate_sigma(double sigma_old, double phi_zero, double n_in, double cap_charge, configuration & config);
double calculate_cap_charge(double sigma_new, double sigma_old, double cap_charge_old, double n_in, configuration & config);

double physical_space(double logical_position, double a, double b, double n_mesh);
int logical_space(float physical_position, float a, float b, float n_mesh);

#endif
