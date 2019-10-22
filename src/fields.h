#ifndef FIELDS_H
#define FIELDS_H

#include "fmatrix.h"

void init_mesh(fmatrix & mesh_x, fmatrix & mesh_y, double a_x, double a_y, int n_mesh_x, int n_mesh_y);
void init_volume_mesh(fmatrix & vmesh, fmatrix & mesh_x, fmatrix & mesh_y);

double ac_voltage_at_time(size_t i, double dt, double freq_hz, double amplitude, double phase);
void calculate_efield(fmatrix & efield_x, fmatrix & efield_y, fmatrix & phi, fmatrix & w_i, fmatrix & w_e, fmatrix & mesh_x , fmatrix & mesh_y, fmatrix & vmesh);

double calculate_phi_zero(double sigma_old, double n_in, double q_cap, double sigma_laplace, fmatrix & phi_poisson, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, imatrix & electrode_mask);
double sigma_from_phi(fmatrix & phi, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, imatrix & electrode_mask);
double calculate_sigma(double sigma_old, double phi_zero, double n_in, double cap_charge);
double calculate_cap_charge(double sigma_new, double sigma_old, double cap_charge_old, double n_in);

double physical_space(double logical_position, double a, double b, double n_mesh);
int logical_space(float physical_position, float a, float b, float n_mesh);

double k1_x(fmatrix & mesh_x, size_t i, size_t j, int ioff=0, int joff=0);
double k2_x(fmatrix & mesh_x, size_t i, size_t j, int ioff=0, int joff=0);
double k3_x(fmatrix & mesh_x, size_t i, size_t j, int ioff=0, int joff=0);
double k1_y(fmatrix & mesh_y, size_t i, size_t j, int ioff=0, int joff=0);
double k2_y(fmatrix & mesh_y, size_t i, size_t j, int ioff=0, int joff=0);
double k3_y(fmatrix & mesh_y, size_t i, size_t j, int ioff=0, int joff=0);

#endif
