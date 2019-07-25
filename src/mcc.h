#ifndef MCC_H
#define MCC_H

#include "fmatrix.h"

double find_nu_prime_e(fmatrix & wmesh_n, fmatrix & vmesh);
double find_nu_prime_i(fmatrix & wmesh_n, fmatrix & vmesh);
double p_null(double nu_prime, double dt);

double kinetic_energy_ev(const fmatrix & p, const int & i, double const & mass);

fmatrix isotropic_scatter(fmatrix & p, const int & i, double chi);
void electron_elastic_collision(fmatrix & p, const int & i, const double kinetic_energy);
void electron_excitation_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double excitation_energy);
void electron_ionization_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double ionization_energy);
void ion_isotropic_collision(fmatrix & p, const int & i, const double kinetic_energy, const double ion_mass);

void collisions_e(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & p_i, int & n_active_i, imatrix & lpos_i, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n, double ion_mass, double p_null, double nu_prime);
void collisions_i(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n, double ion_mass, double p_null, double nu_prime);

#endif

