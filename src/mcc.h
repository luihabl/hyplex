#ifndef MCC_H
#define MCC_H

#include "fmatrix.h"
#include "configuration.h"
#include "fields.h"
#include "particles.h"
#include "particles-in-mesh.h"


struct mcc
{
    
    // cross-sections
    fmatrix elastic_cs;
    fmatrix ionization_cs;
    fmatrix* excitation_cs;
    fmatrix isotropic_cs;
    fmatrix backscattering_cs;

    double e_iz;
    fmatrix e_exc;
    int n_exc;

    // configuration variables
    double dt, dx, dy, a_x, a_y, m_i, m_el, q, pi, t_neutral;
    int k_sub;

    //mcc variables
    double nu_prime_e, nu_prime_i;
    double p_null_e, p_null_i;

    particle_operations * pops;
    pic_operations * pic;

    mcc() = default;
    mcc(configuration & config, particle_operations * pops, pic_operations * pic);
    ~mcc();

    double find_nu_prime_e(fmatrix & dens_n, fmatrix & vmesh);
    double find_nu_prime_i(fmatrix & dens_n, fmatrix & vmesh);
    double p_null(double nu_prime, double k_sub);
    double collision_frequency(double neutral_density, double cross_section, double kinetic_energy_ev, double mass);
    double kinetic_energy_ev(const fmatrix & p, const int & i, double const & mass);
    double calc_total_cs(double energy, int n_exc);

    fmatrix isotropic_scatter(fmatrix & p, const int & i, double chi);

    void initialize_mcc(fmatrix & dens_n, fmatrix & vmesh);
    void electron_elastic_collision(fmatrix & p, const int & i, const double & kinetic_energy);
    void electron_excitation_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double excitation_energy);
    void electron_ionization_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double ionization_energy);
    void ion_isotropic_collision(fmatrix & p, const int & i, const double kinetic_energy);

    int collisions_e(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & p_i, int & n_active_i, imatrix & lpos_i, mesh_set mesh, fmatrix & dens_n);
    void collisions_i(fmatrix & p, int & n_active, imatrix & lpos, mesh_set & mesh, fmatrix & dens_n);

};




#endif

