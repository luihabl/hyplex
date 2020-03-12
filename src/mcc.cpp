#include "mcc.h"

#include <cmath>

#include "fmatrix.h"
#include "fmath.h"

#include "configuration.h"
#include "num-tools.h"
#include "fields.h"
#include "particles.h"
#include "random-numbers.h"
#include "particles-in-mesh.h"
#include "cross-sections.h"



mcc::mcc(configuration & config)
{

    dt = config.f("time/dt");
    k_sub = config.i("time/k_sub");
    dx = config.f("p/dx");
    dy = config.f("p/dy");
    a_x = config.f("geometry/a_x");
    a_y = config.f("geometry/a_y");
    m_i = config.f("ugas/m_i");
    m_el = config.f("electrons/m_el");
    q = config.f("physical/q");
    pi = config.f("physical/pi");
    e_iz = config.f("ugas/e_iz");
    t_neutral = config.f("neutrals/t_neutral");
    n_exc = config.i("p/n_exc");
    e_exc = config.fs("ugas/e_exc");

}


void mcc::initialize_mcc(fmatrix & dens_n, fmatrix & vmesh){
    nu_prime_e       = find_nu_prime_e(dens_n, vmesh);
    p_null_e         = p_null(nu_prime_e, 1);
    nu_prime_i       = find_nu_prime_i(dens_n, vmesh);
    p_null_i         = p_null(nu_prime_i, k_sub);
}


// --------------------------------- Initial parameters for the MCC ----------------------------------

double calc_total_cs(double energy, int n_exc){
    double total_cs;
    total_cs = interp(elastic_cs, energy);
    total_cs += interp(ionization_cs, energy);
    
    for(int i = 0; i < n_exc; i++){
        total_cs += interp(excitation_cs[i], energy);
    }
    
    return total_cs;
}

double mcc::find_nu_prime_e(fmatrix & dens_n, fmatrix & vmesh)
{
	double nu_prime = 0.0;
    double neutral_density = dens_n.max();
    
	double total_cs;
	double velocity;
	double energy;

    for (size_t i = 0; i < elastic_cs.n1; i++)
    {
        energy = elastic_cs.val[i * 2 + 0];

        total_cs = calc_total_cs(energy, n_exc);

        velocity = sqrt(2 * energy * q / m_el);
        nu_prime = nu_prime < neutral_density * total_cs * velocity ? neutral_density * total_cs * velocity : nu_prime;
    }
    
    for (size_t i = 0; i < ionization_cs.n1; i++)
    {
        energy = ionization_cs.val[i * 2 + 0];

        total_cs = calc_total_cs(energy, n_exc);

        velocity = sqrt(2 * energy * q / m_el);
        nu_prime = nu_prime < neutral_density * total_cs * velocity ? neutral_density * total_cs * velocity : nu_prime;
    }
    
    for(int n = 0; n < n_exc; n++){
        for (size_t i = 0; i < excitation_cs[n].n1; i++)
        {
            energy = excitation_cs[n].val[i * 2 + 0];

            total_cs = calc_total_cs(energy, n_exc);

            velocity = sqrt(2 * energy * q / m_el);
            nu_prime = nu_prime < neutral_density * total_cs * velocity ? neutral_density * total_cs * velocity : nu_prime;
        }
    }
    
	return nu_prime;
}

double mcc::find_nu_prime_i(fmatrix & dens_n, fmatrix & vmesh)
{
	double nu_prime = 0.0;
    double neutral_density = dens_n.max();

	double total_cs;
	double velocity;
	double energy;
       
    for (size_t i = 0; i < isotropic_cs.n1; i++)
    {
        energy = isotropic_cs.val[i * 2 + 0];
        total_cs = interp(isotropic_cs, energy) + interp(backscattering_cs, energy);
        velocity = sqrt(4 * energy * q / m_i); // x4 to ensure correct reference frame
        nu_prime = nu_prime < neutral_density * total_cs * velocity ? neutral_density * total_cs * velocity : nu_prime;
    }
    
	return nu_prime;
}

double mcc::p_null(double nu_prime, double k_sub)
{
	return 1 - exp(-nu_prime * dt * k_sub);
}

// ------------------------------ Auxiliary functions ------------------------------------

double mcc::collision_frequency(double neutral_density, double cross_section, double kinetic_energy_ev, double mass) {
	return neutral_density * cross_section * sqrt(2 * q * kinetic_energy_ev / mass);
}

double mcc::kinetic_energy_ev(const fmatrix & p, const int & i, double const & mass)
{
    return 0.5 * mass * (pow(p.val[i * 6 + 3] * dx / dt, 2) + pow(p.val[i * 6 + 4] * dx / dt, 2) + pow(p.val[i * 6 + 5] * dx / dt, 2)) / q;
}

// ------------------------------ Electron collisions ------------------------------------

void mcc::collisions_e(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & p_i, int & n_active_i, imatrix & lpos_i, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n)
{   
    
	imatrix particle_samples;
	double kinetic_energy = 0.0;
    double neutral_density = 0.0;
    
    int n_null = (int) floor(p_null_e * n_active);
    n_null = p_null_e * n_active - n_null > r_unif() ? n_null + 1 : n_null;


	// empirical formula for faster sampling algorithm (check the formula again for Release instead of Debug mode)!!! 
    particle_samples = sample_from_sequence_naive(n_null, n_active);
//    if (n_active > 0.03 * pow(n_null, 1.8067))
//        particle_samples = sample_from_sequence_naive(n_null, n_active);
//    else
//        particle_samples = sample_from_sequence_shuffle(n_null, n_active);

	int i = 0;
    double random_number_1 = 0.0;
	double freq_ratio_0 = 0.0;
	double freq_ratio_1 = 0.0;
	for (int k = 0; k < n_null; k++)
	{
		i = particle_samples.val[k];
		kinetic_energy = kinetic_energy_ev(p, i, m_el);
        neutral_density = field_at_position(dens_n, mesh_x, mesh_y, p.val[i * 6 + 0], p.val[i * 6 + 1], lpos.val[i * 2 + 0], lpos.val[i * 2 + 1], a_x, a_y, dx, dy);
		random_number_1 = r_unif();

		// Elastic collision:
		freq_ratio_0 = 0;
		freq_ratio_1 = collision_frequency(neutral_density, interp(elastic_cs, kinetic_energy), kinetic_energy, m_el) / nu_prime_e;
		if (random_number_1 <= freq_ratio_1)
		{
            electron_elastic_collision(p, i, kinetic_energy);
			continue;
		}

        for(int n = 0; n < n_exc; n++){
            // Excitation collisions:
            freq_ratio_0 = freq_ratio_1;
            freq_ratio_1 += collision_frequency(neutral_density, interp(excitation_cs[n], kinetic_energy), kinetic_energy, m_el) / nu_prime_e;
            if (random_number_1 > freq_ratio_0 && random_number_1 <= freq_ratio_1 && kinetic_energy >= e_exc.val[n])
            {
                electron_excitation_collision(p, i, kinetic_energy, e_exc.val[n]);
                continue;
            }
        }		
	
		// Ionization collision:
		freq_ratio_0 = freq_ratio_1;
        freq_ratio_1 += collision_frequency(neutral_density, interp(ionization_cs, kinetic_energy), kinetic_energy, m_el) / nu_prime_e;
		if (random_number_1 > freq_ratio_0 && random_number_1 <= freq_ratio_1 && kinetic_energy >= e_iz)
		{
            int new_electron_index = add_particle_copy(p, n_active, lpos, i);
            electron_ionization_collision(p, i, kinetic_energy, e_iz);
            electron_ionization_collision(p, new_electron_index, kinetic_energy, e_iz);
            add_maxwellian_particle_at_position(p_i, n_active_i, lpos_i, t_neutral, m_i, p.val[i * 6 + 0], p.val[i * 6 + 1], lpos.val[i * 2 + 0] , lpos.val[i * 2 + 1], dt, dx, q); // <========== particle is added here
			continue;
		}
	}
}

void mcc::electron_elastic_collision(fmatrix & p, const int & i, const double & kinetic_energy)
{
    double chi = acos(1.0 - 2.0 * r_unif());
    fmatrix v_scat = isotropic_scatter(p, i, chi);
    
    double delta_energy = (2 * m_el / m_i) * (1 - cos(chi));
    const double velocity_multiplier = (dt / dx) * sqrt(2 * q * (kinetic_energy * (1 - delta_energy)) / m_el);
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
}

void mcc::electron_excitation_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double excitation_energy)
{
    double chi = acos(1.0 - 2.0 * r_unif());
    fmatrix v_scat = isotropic_scatter(p, i, chi);
    
    const double velocity_multiplier = (dt / dx) * sqrt(2 * q * (kinetic_energy - excitation_energy) / m_el);
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
    
}

void mcc::electron_ionization_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double ionization_energy)
{
    double chi = acos(1.0 - 2.0 * r_unif());
    fmatrix v_scat = isotropic_scatter(p, i, chi);
    
    const double velocity_multiplier = (dt / dx) * sqrt(q * (kinetic_energy - ionization_energy) / m_el); // No x 2 because of ionization energy division
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
}


// ------------------------------ Ion collisions ------------------------------------

void mcc::collisions_i(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n)
{

	imatrix particle_samples;
    fmatrix v_neutral = fmatrix::zeros(3);
    double kinetic_energy_relative = 0.0;
    double neutral_density = 0.0;
	
	int n_null = (int) floor(p_null_i * n_active);
    n_null = p_null_i * n_active - n_null > r_unif() ? n_null + 1 : n_null;
    
    double random_number_1 = 0.0;

    particle_samples = sample_from_sequence_naive(n_null, n_active);

	int i = 0;
	double freq_ratio_0 = 0.0;
	double freq_ratio_1 = 0.0;
	for (int k = 0; k < n_null; k++)
	{
		i = particle_samples.val[k];
        neutral_density = field_at_position(dens_n, mesh_x, mesh_y, p.val[i * 6 + 0], p.val[i * 6 + 1], lpos.val[i * 2 + 0], lpos.val[i * 2 + 1], a_x, a_y, dx, dy);
        
        v_neutral = {r_norm() * (dt / dx) * sqrt(q * t_neutral / m_i),
                     r_norm() * (dt / dx) * sqrt(q * t_neutral / m_i),
                     r_norm() * (dt / dx) * sqrt(q * t_neutral / m_i)};

        p.val[i * 6 + 3] = p.val[i * 6 + 3] -  v_neutral.val[0];
        p.val[i * 6 + 4] = p.val[i * 6 + 4] -  v_neutral.val[1];
        p.val[i * 6 + 5] = p.val[i * 6 + 5] -  v_neutral.val[2];
        
        kinetic_energy_relative = kinetic_energy_ev(p, i, m_i);
        
		random_number_1 = r_unif();

		// Isotropic collision:
		freq_ratio_0 = 0.0;
		freq_ratio_1 = collision_frequency(neutral_density, interp(isotropic_cs, 0.5 * kinetic_energy_relative), kinetic_energy_relative, m_i) / nu_prime_i;
		if (random_number_1 <= freq_ratio_1)
		{
            ion_isotropic_collision(p, i, kinetic_energy_relative);
            p.val[i * 6 + 3] = p.val[i * 6 + 3] + v_neutral.val[0];
            p.val[i * 6 + 4] = p.val[i * 6 + 4] + v_neutral.val[1];
            p.val[i * 6 + 5] = p.val[i * 6 + 5] + v_neutral.val[2];

            continue;
		}

		// Backscattering collision:
		freq_ratio_0 = freq_ratio_1;
		freq_ratio_1 += collision_frequency(neutral_density, interp(backscattering_cs, 0.5 * kinetic_energy_relative), kinetic_energy_relative, m_i) / nu_prime_i;
		if (random_number_1 > freq_ratio_0 && random_number_1 <= freq_ratio_1)
		{
            p.val[i * 6 + 3] = v_neutral.val[0];
            p.val[i * 6 + 4] = v_neutral.val[1];
            p.val[i * 6 + 5] = v_neutral.val[2];

			continue;
		}
        
        // Null collision:
        if(random_number_1 > freq_ratio_1){
            p.val[i * 6 + 3] = p.val[i * 6 + 3] + v_neutral.val[0];
            p.val[i * 6 + 4] = p.val[i * 6 + 4] + v_neutral.val[1];
            p.val[i * 6 + 5] = p.val[i * 6 + 5] + v_neutral.val[2];
        }
	}
}

void mcc::ion_isotropic_collision(fmatrix & p, const int & i, const double kinetic_energy)
{
    
    double chi = acos(sqrt(1.0 - r_unif()));
    fmatrix v_scat = isotropic_scatter(p, i, chi);
    
    const double velocity_multiplier = (dt / dx) * sqrt(2 * q * (kinetic_energy * pow(cos(chi), 2)) / m_i);
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
    
}

// -------------------------- Scattering function -------------------------------------

fmatrix mcc::isotropic_scatter(fmatrix & p, const int & i, double chi)
{
    
    fmatrix v = { p.val[i * 6 + 3], p.val[i * 6 + 4], p.val[i * 6 + 5] };
    v = v / norm(v);
    
    double vx = v.val[0];
    double vy = v.val[1];
    double vz = v.val[2];
    
    double phi = 2 * pi * r_unif();
    double zeta = acos(vz);
    
    v = {
        vx*cos(chi) + vy*sin(chi)*sin(phi)/sin(zeta) + vx*vz*sin(chi)*cos(phi)/sin(zeta),
        vy*cos(chi) - vx*sin(chi)*sin(phi)/sin(zeta) + vy*vz*sin(chi)*cos(phi)/sin(zeta),
        vz*cos(chi) - (pow(vx, 2) + pow(vy, 2))*sin(chi)*cos(phi)/sin(zeta)
    };
    
    return v;
}
