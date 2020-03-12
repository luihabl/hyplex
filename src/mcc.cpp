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

double find_nu_prime_e(fmatrix & dens_n, fmatrix & vmesh, configuration & config)
{
	double nu_prime = 0.0;
    double neutral_density = dens_n.max();
    
	double total_cs;
	double velocity;
	double energy;

    const double q = config.f("physical/q");
    const double m_el = config.f("electrons/m_el");
    const int n_exc = config.i("p/n_exc");

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

double find_nu_prime_i(fmatrix & dens_n, fmatrix & vmesh, configuration & config)
{
	double nu_prime = 0.0;
    double neutral_density = dens_n.max();

	double total_cs;
	double velocity;
	double energy;

    const double q = config.f("physical/q");
    const double m_i = config.f("ugas/m_i");
       
    for (size_t i = 0; i < isotropic_cs.n1; i++)
    {
        energy = isotropic_cs.val[i * 2 + 0];
        total_cs = interp(isotropic_cs, energy) + interp(backscattering_cs, energy);
        velocity = sqrt(4 * energy * q / m_i); // x4 to ensure correct reference frame
        nu_prime = nu_prime < neutral_density * total_cs * velocity ? neutral_density * total_cs * velocity : nu_prime;
    }
    
	return nu_prime;
}

double p_null(double nu_prime, double dt)
{
	return 1 - exp(-nu_prime * dt);
}

// ------------------------------ Auxiliary functions ------------------------------------

double collision_frequency(double neutral_density, double cross_section, double kinetic_energy_ev, double mass, double q) {
	return neutral_density * cross_section * sqrt(2 * q * kinetic_energy_ev / mass);
}

double kinetic_energy_ev(const fmatrix & p, const int & i, double const & mass, double const & dt, double const & dx, double const & q)
{
    return 0.5 * mass * (pow(p.val[i * 6 + 3] * dx / dt, 2) + pow(p.val[i * 6 + 4] * dx / dt, 2) + pow(p.val[i * 6 + 5] * dx / dt, 2)) / q;
}

// ------------------------------ Electron collisions ------------------------------------

void collisions_e(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & p_i, int & n_active_i, imatrix & lpos_i, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n, double p_null, double nu_prime, configuration & config)
{   
    const double dt = config.f("time/dt");
    const double dx = config.f("p/dx");
    const double m_i = config.f("ugas/m_i");
    const double m_el = config.f("electrons/m_el");
    const double q = config.f("physical/q");
    const double pi = config.f("physical/pi");
    const double e_iz = config.f("ugas/e_iz");
    const double t_neutral = config.f("neutrals/t_neutral");
    const int n_exc = config.i("p/n_exc");
    const fmatrix e_exc = config.fs("ugas/e_exc");
    
	imatrix particle_samples;
	double kinetic_energy = 0.0;
    double neutral_density = 0.0;
    
    int n_null = (int) floor(p_null * n_active);
    n_null = p_null * n_active - n_null > r_unif() ? n_null + 1 : n_null;


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
		kinetic_energy = kinetic_energy_ev(p, i, m_el, dt, dx, q);
        neutral_density = field_at_position(dens_n, mesh_x, mesh_y, p.val[i * 6 + 0], p.val[i * 6 + 1], lpos.val[i * 2 + 0], lpos.val[i * 2 + 1]);
		random_number_1 = r_unif();

		// Elastic collision:
		freq_ratio_0 = 0;
		freq_ratio_1 = collision_frequency(neutral_density, interp(elastic_cs, kinetic_energy), kinetic_energy, m_el, q) / nu_prime;
		if (random_number_1 <= freq_ratio_1)
		{
            electron_elastic_collision(p, i, kinetic_energy, m_el, m_i, dt, dx, q, pi);
			continue;
		}

        for(int n = 0; n < n_exc; n++){
            // Excitation collisions:
            freq_ratio_0 = freq_ratio_1;
            freq_ratio_1 += collision_frequency(neutral_density, interp(excitation_cs[n], kinetic_energy), kinetic_energy, m_el, q) / nu_prime;
            if (random_number_1 > freq_ratio_0 && random_number_1 <= freq_ratio_1 && kinetic_energy >= e_exc.val[n])
            {
                electron_excitation_collision(p, i, kinetic_energy, e_exc.val[n], m_el, dt, dx, q, pi);
                continue;
            }
        }		
	
		// Ionization collision:
		freq_ratio_0 = freq_ratio_1;
        freq_ratio_1 += collision_frequency(neutral_density, interp(ionization_cs, kinetic_energy), kinetic_energy, m_el, q) / nu_prime;
		if (random_number_1 > freq_ratio_0 && random_number_1 <= freq_ratio_1 && kinetic_energy >= e_iz)
		{
            int new_electron_index = add_particle_copy(p, n_active, lpos, i);
            electron_ionization_collision(p, i, kinetic_energy, e_iz, m_el, dt, dx, q, pi);
            electron_ionization_collision(p, new_electron_index, kinetic_energy, e_iz, m_el, dt, dx, q, pi);
            add_maxwellian_particle_at_position(p_i, n_active_i, lpos_i, t_neutral, m_i, p.val[i * 6 + 0], p.val[i * 6 + 1], lpos.val[i * 2 + 0] , lpos.val[i * 2 + 1], dt, dx, q); // <========== particle is added here
			continue;
		}
	}
}

void electron_elastic_collision(fmatrix & p, const int & i, const double & kinetic_energy, const double & m_el, const double & m_i, const double & dt, const double & dx, const double & q, const double & pi)
{
    double chi = acos(1.0 - 2.0 * r_unif());
    fmatrix v_scat = isotropic_scatter(p, i, chi, pi);
    
    double delta_energy = (2 * m_el / m_i) * (1 - cos(chi));
    const double velocity_multiplier = (dt / dx) * sqrt(2 * q * (kinetic_energy * (1 - delta_energy)) / m_el);
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
}

void electron_excitation_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double excitation_energy, const double & m_el, const double & dt, const double & dx, const double & q, const double & pi)
{
    double chi = acos(1.0 - 2.0 * r_unif());
    fmatrix v_scat = isotropic_scatter(p, i, chi, pi);
    
    const double velocity_multiplier = (dt / dx) * sqrt(2 * q * (kinetic_energy - excitation_energy) / m_el);
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
    
}

void electron_ionization_collision(fmatrix & p, const int & i, const double kinetic_energy,  const double ionization_energy, const double & m_el, const double & dt, const double & dx, const double & q, const double & pi)
{
    double chi = acos(1.0 - 2.0 * r_unif());
    fmatrix v_scat = isotropic_scatter(p, i, chi, pi);
    
    const double velocity_multiplier = (dt / dx) * sqrt(q * (kinetic_energy - ionization_energy) / m_el); // No x 2 because of ionization energy division
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
}


// ------------------------------ Ion collisions ------------------------------------

void collisions_i(fmatrix & p, int & n_active, imatrix & lpos, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & dens_n, double p_null, double nu_prime, configuration & config)
{

    const double dt = config.f("time/dt");
    const double dx = config.f("p/dx");
    const double m_i = config.f("ugas/m_i");
    const double q = config.f("physical/q");
    const double pi = config.f("physical/pi");
    const double t_neutral = config.f("neutrals/t_neutral");

	imatrix particle_samples;
    fmatrix v_neutral = fmatrix::zeros(3);
    double kinetic_energy_relative = 0.0;
    double neutral_density = 0.0;
	
	int n_null = (int) floor(p_null * n_active);
    n_null = p_null * n_active - n_null > r_unif() ? n_null + 1 : n_null;
    
    double random_number_1 = 0.0;

    particle_samples = sample_from_sequence_naive(n_null, n_active);

	int i = 0;
	double freq_ratio_0 = 0.0;
	double freq_ratio_1 = 0.0;
	for (int k = 0; k < n_null; k++)
	{
		i = particle_samples.val[k];
        neutral_density = field_at_position(dens_n, mesh_x, mesh_y, p.val[i * 6 + 0], p.val[i * 6 + 1], lpos.val[i * 2 + 0], lpos.val[i * 2 + 1]);
        
        v_neutral = {r_norm() * (dt / dx) * sqrt(q * t_neutral / m_i),
                     r_norm() * (dt / dx) * sqrt(q * t_neutral / m_i),
                     r_norm() * (dt / dx) * sqrt(q * t_neutral / m_i)};

        p.val[i * 6 + 3] = p.val[i * 6 + 3] -  v_neutral.val[0];
        p.val[i * 6 + 4] = p.val[i * 6 + 4] -  v_neutral.val[1];
        p.val[i * 6 + 5] = p.val[i * 6 + 5] -  v_neutral.val[2];
        
        kinetic_energy_relative = kinetic_energy_ev(p, i, m_i, dt, dx, q);
        
		random_number_1 = r_unif();

		// Isotropic collision:
		freq_ratio_0 = 0.0;
		freq_ratio_1 = collision_frequency(neutral_density, interp(isotropic_cs, 0.5 * kinetic_energy_relative), kinetic_energy_relative, m_i, q) / nu_prime;
		if (random_number_1 <= freq_ratio_1)
		{
            ion_isotropic_collision(p, i, kinetic_energy_relative, m_i, dt, dx, q, pi);
            p.val[i * 6 + 3] = p.val[i * 6 + 3] + v_neutral.val[0];
            p.val[i * 6 + 4] = p.val[i * 6 + 4] + v_neutral.val[1];
            p.val[i * 6 + 5] = p.val[i * 6 + 5] + v_neutral.val[2];

            continue;
		}

		// Backscattering collision:
		freq_ratio_0 = freq_ratio_1;
		freq_ratio_1 += collision_frequency(neutral_density, interp(backscattering_cs, 0.5 * kinetic_energy_relative), kinetic_energy_relative, m_i, q) / nu_prime;
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

void ion_isotropic_collision(fmatrix & p, const int & i, const double kinetic_energy, const double m_i, const double & dt, const double & dx, const double & q, const double & pi)
{
    
    double chi = acos(sqrt(1.0 - r_unif()));
    fmatrix v_scat = isotropic_scatter(p, i, chi, pi);
    
    const double velocity_multiplier = (dt / dx) * sqrt(2 * q * (kinetic_energy * pow(cos(chi), 2)) / m_i);
    
    p.val[i * 6 + 3] = v_scat.val[0] * velocity_multiplier;
    p.val[i * 6 + 4] = v_scat.val[1] * velocity_multiplier;
    p.val[i * 6 + 5] = v_scat.val[2] * velocity_multiplier;
    
}

// -------------------------- Scattering function -------------------------------------

fmatrix isotropic_scatter(fmatrix & p, const int & i, double chi, const double & pi)
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
