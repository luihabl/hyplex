#ifndef MCC_H
#define MCC_H

#include "fmatrix.h"
#include "configuration.h"
#include "fields.h"


struct mcc
{
    

    // configuration variables
    double dt, dx, dy, a_x, a_y, m_i, m_el, q, pi, e_iz, t_neutral, n_exc;
    int k_sub;
    fmatrix e_exc;

    //mcc variables
    double nu_prime_e, nu_prime_i;
    double p_null_e, p_null_i;

    mcc(configuration & config);

    double find_nu_prime_e(fmatrix & dens_n, fmatrix & vmesh);
    double find_nu_prime_i(fmatrix & dens_n, fmatrix & vmesh);
    double p_null(double nu_prime, double k_sub);
    double collision_frequency(double neutral_density, double cross_section, double kinetic_energy_ev, double mass);
    double kinetic_energy_ev(const fmatrix & p, const int & i, double const & mass);

    fmatrix isotropic_scatter(fmatrix & p, const int & i, double chi);

    void initialize_mcc(fmatrix & dens_n, fmatrix & vmesh);
    void electron_elastic_collision(fmatrix & p, const int & i, const double & kinetic_energy);
    void electron_excitation_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double excitation_energy);
    void electron_ionization_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double ionization_energy);
    void ion_isotropic_collision(fmatrix & p, const int & i, const double kinetic_energy);

    void collisions_e(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & p_i, int & n_active_i, imatrix & lpos_i, mesh_set mesh, fmatrix & dens_n);
    void collisions_i(fmatrix & p, int & n_active, imatrix & lpos, mesh_set & mesh, fmatrix & dens_n);

};





#endif

