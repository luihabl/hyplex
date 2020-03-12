#ifndef MCC_H
#define MCC_H

#include "fmatrix.h"
#include "configuration.h"

double find_nu_prime_e(fmatrix & dens_n, fmatrix & vmesh, configuration & config);
double find_nu_prime_i(fmatrix & dens_n, fmatrix & vmesh, configuration & config);
double p_null(double nu_prime, double dt);

double kinetic_energy_ev(const fmatrix & p, const int & i, double const & mass, double const & dt, double const & dx, double const & q);

fmatrix isotropic_scatter(fmatrix & p, const int & i, double chi, const double & pi);
void electron_elastic_collision(fmatrix & p, const int & i, const double & kinetic_energy, const double & m_el, const double & m_i, const double & dt, const double & dx, const double & q, const double & pi);
void electron_excitation_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double excitation_energy, const double & m_el, const double & dt, const double & dx, const double & q, const double & pi);
void electron_ionization_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double ionization_energy, const double & m_el, const double & dt, const double & dx, const double & q, const double & pi);
void ion_isotropic_collision(fmatrix & p, const int & i, const double kinetic_energy, const double m_i, const double & dt, const double & dx, const double & q, const double & pi);

void collisions_e(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & p_i, int & n_active_i, imatrix & lpos_i, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n, double p_null, double nu_prime, configuration & config);
void collisions_i(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n, double p_null, double nu_prime, configuration & config);

#endif

