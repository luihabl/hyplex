/*
    Functions to add, move and remove particles.
 */


#include "particles.h"

#include <iostream>
#include <algorithm>
#include <cmath>

#include "fmatrix.h"
#include "num-tools.h"
#include "configuration.h"
#include "random-numbers.h"
#include "particles-in-mesh.h"


using namespace std;


particle_operations::particle_operations(configuration & config, pic_operations & _pic): pic(_pic){


	n_mesh_x = config.i("geometry/n_mesh_x");
	n_mesh_y = config.i("geometry/n_mesh_y");
	n_thruster = config.i("geometry/n_thruster");
	a_x = config.f("geometry/a_x");
	a_y = config.f("geometry/a_x");

	dx = config.f("p/dx");
	dy = config.f("p/dy");
	dt = config.f("time/dt");

	temp_e = config.f("electrons/t_el");
	m_el = config.f("electrons/m_el");
	q = config.f("physical/q");

	k_inj = config.f("p/k_inj_el");
	v_sb = config.f("thruster/v_sb");
	v_rf = config.f("thruster/v_rf");
	freq = config.f("thruster/freq");
	duty_cycle = config.f("electrons/sqr_duty_cycle");
	
	n_factor = config.f("particles/n_factor");
	k_sub = config.i("time/k_sub");
	i_i = config.f("ions/i_i");
	alpha = config.f("p/alpha");
	omega_i = config.f("p/omega_i");
	k_phi = config.f("p/k_phi");

	c_cap = config.f("boundaries/c_cap");
}


//  ----------------------------- Particle Creation ---------------------------

void particle_operations::add_maxwellian_particles(fmatrix & p, int & n_active, const double temperature, const double mass, const size_t n_add)
{
	double v_temperature = sqrt(q * temperature / mass);

	for (size_t i = n_active; i < n_active + n_add; i++)
	{
		p.val[i * 6 + 0] = ((double) n_mesh_x - 1.0) * r_unif();
		p.val[i * 6 + 1] = (dy / dx) * ((double) n_mesh_y - 1.0) * r_unif();
		p.val[i * 6 + 2] = 0.0;
		p.val[i * 6 + 3] = (dt / dx) * r_norm(0.0, v_temperature);
		p.val[i * 6 + 4] = (dt / dx) * r_norm(0.0, v_temperature);
		p.val[i * 6 + 5] = (dt / dx) * r_norm(0.0, v_temperature);
	}
	n_active += n_add;
}

void particle_operations::add_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub){
	
	double f_n_add = floor(k_sub * n_add);
	int n_new = (r_unif() <= (k_sub * n_add - f_n_add) ?  f_n_add + 1 : f_n_add);

	double v_temperature = sqrt(q * temperature / mass);

	for (int i = n_active; i < n_active + n_new; i++)
	{	
		p.val[i * 6 + 3] = (dt / dx) * (v_temperature * sqrt(- 2 *  log(r_unif())) + v_drift);
		// p.val[i * 6 + 3] = (DT / DX) * sqrt(pow(v_temperature * sqrt(- log(r_unif())), 2) + (v_drift * v_drift));
		// p.val[i * 6 + 3] = (DT / DX) * (r_norm(0.0, v_temperature) + v_drift);
		p.val[i * 6 + 4] = (dt / dx) * r_norm(0.0, v_temperature);
		p.val[i * 6 + 5] = (dt / dx) * r_norm(0.0, v_temperature);

		p.val[i * 6 + 0] = k_sub * p.val[i * 6 + 3] * r_unif();
		p.val[i * 6 + 1] = (dy / dx) * ((double) n_thruster - 1.0) * r_unif();
		p.val[i * 6 + 2] = 0.0;
	}
	n_active += n_new;
}

void particle_operations::add_maxwellian_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const double n_add, double k_sub){
	
	double f_n_add = floor(k_sub * n_add);
	int n_new = (r_unif() < (k_sub * n_add - f_n_add) ?  f_n_add + 1 : f_n_add);
	
	double v_temperature = sqrt(q * temperature / mass);

	for (int i = n_active; i < n_active + n_new; i++)
	{	
		do {
			p.val[i * 6 + 3] = abs((dt / dx) * r_norm(0.0, v_temperature));
		} while (p.val[i * 6 + 3] < 0.0);

		p.val[i * 6 + 4] = (dt / dx) * r_norm(0.0, v_temperature);
		p.val[i * 6 + 5] = (dt / dx) * r_norm(0.0, v_temperature);

		p.val[i * 6 + 0] = k_sub * p.val[i * 6 + 3] * r_unif(); // maybe necessary to multiply K_SUB here!
		p.val[i * 6 + 1] = (dy / dx) * ((double) n_thruster - 1.0) * r_unif();
		p.val[i * 6 + 2] = 0.0;
	}
	n_active += n_new;
}

void particle_operations::add_maxwellian_particle_at_position(fmatrix & p, int & n_active, imatrix & lpos, const double temperature, const double mass, double x_pos, double y_pos, int lpos_x, int lpos_y)
{
    double v_temperature = sqrt(q * temperature / mass);
    
    p.val[n_active * 6 + 0] = x_pos;
    p.val[n_active * 6 + 1] = y_pos;
    p.val[n_active * 6 + 2] = 0.0;
    p.val[n_active * 6 + 3] = (dt / dx) * r_norm(0.0, v_temperature);
    p.val[n_active * 6 + 4] = (dt / dx) * r_norm(0.0, v_temperature);
    p.val[n_active * 6 + 5] = (dt / dx) * r_norm(0.0, v_temperature);
	lpos.val[n_active * 2 + 0] = lpos_x;
	lpos.val[n_active * 2 + 1] = lpos_y;

    n_active += 1;
}

int particle_operations::add_particle_copy(fmatrix & p, int & n_active, imatrix & lpos, const int & i)
{
    p.val[n_active * 6 + 0] = p.val[i * 6 + 0];
    p.val[n_active * 6 + 1] = p.val[i * 6 + 1];
    p.val[n_active * 6 + 2] = p.val[i * 6 + 2];
    p.val[n_active * 6 + 3] = p.val[i * 6 + 3];
    p.val[n_active * 6 + 4] = p.val[i * 6 + 4];
    p.val[n_active * 6 + 5] = p.val[i * 6 + 5];

	lpos.val[n_active * 2 + 0] = lpos.val[i * 2 + 0];
	lpos.val[n_active * 2 + 1] = lpos.val[i * 2 + 1];
    
    int copied_position = n_active;
    n_active += 1;
    
    return copied_position;
}

double particle_operations::balanced_injection(double old_n_inj, double rate_constant, fmatrix & wmesh_i, fmatrix & wmesh_e, int ill, int jll, int iur, int jur){

	double dw = 0;
	const int mesh_n2 = (int) wmesh_i.n2;

	for (int i = ill; i < iur + 1; i++)
	{
		for (int j = jll; j < jur + 1; j++)
		{
			dw += wmesh_i.val[i * mesh_n2 + j] -  wmesh_e.val[i * mesh_n2 + j];
		}
	}
	dw = dw/((iur - ill + 1) * (jur - jll + 1));
	double new_n_inj = old_n_inj + rate_constant * dw;

	return new_n_inj > 0 ? new_n_inj : 0;
}

double particle_operations::pulsed_injection(int i){
	return k_inj *  exp(- (v_sb + v_rf * sin(omega_i * (double) i)) / temp_e);
}

double particle_operations::square_injection(int i){

	double alpha_i_i = 4.0;

	double n_period = round(1 / (freq * dt));
	double n_0 = round((duty_cycle) / (freq * dt));
	
	if(i % (int) n_period > (int) (n_period - n_0)){ //if square injection does not work, check this
		return (alpha_i_i / duty_cycle) * (dt / (q * n_factor)) * i_i;
	}
	return 0;
}

//  ----------------------------- Boundaries ----------------------------------

void particle_operations::boundaries_n(fmatrix & p, int & n_active, imatrix & lpos){
	
	int n_remove = 0;
    imatrix tbremoved((size_t) (n_active * 0.5) + 150); 

    double x, y;
    const double x_max = ((double) n_mesh_x - 1);
    const double y_max = ((double) n_mesh_y - 1) * (dy / dx);
	const double y_thr = ((double) n_thruster - 1) * (dy / dx);
    
    for (int i = 0; i < n_active; i++)
    {
        x = p.val[i * 6 + 0];
        y = p.val[i * 6 + 1];
        
        if (x > x_max || y > y_max || (x < 0 && y < y_thr))
        {
            tbremoved.val[n_remove] = i;
            n_remove += 1;

        }
		else if (y < 0) {
            p.val[i * 6 + 1] = - p.val[i * 6 + 1];
            p.val[i * 6 + 4] = - p.val[i * 6 + 4];
        }
		else if (x < 0 && y >= y_thr) {
			p.val[i * 6 + 0] = - p.val[i * 6 + 0];
            p.val[i * 6 + 3] = - p.val[i * 6 + 3];
		}
    }
    
    for (int i = n_remove - 1; i >= 0; i--)
    {
        remove_particle(p, n_active, tbremoved.val[i], lpos);
    }
}

void particle_operations::boundaries_ob_count(fmatrix & p, int & n_active, imatrix & lpos, int & n_removed_ob, int & n_removed_thr)
{
	n_removed_ob = 0;
    n_removed_thr = 0;
	int n_remove = 0;
	static imatrix tbremoved(100000);

	double x, y;
	const double x_max = ((double) n_mesh_x - 1);
	const double y_max = ((double) n_mesh_y - 1) * (dy / dx);
	const double y_thr = ((double) n_thruster - 1) * (dy / dx);

	for (int i = 0; i < n_active; i++)
	{
		x = p.val[i * 6 + 0];
		y = p.val[i * 6 + 1];

		if (x > x_max || y > y_max || (x < 0 && y < y_thr))
		{
			tbremoved.val[n_remove] = i;
			n_remove += 1;

			if(x > x_max || y > y_max){
				n_removed_ob += 1;
			}
            
            if(y < y_thr){
                n_removed_thr += 1;
            }

		} 
		else if (y < 0) {
            p.val[i * 6 + 1] = - p.val[i * 6 + 1];
            p.val[i * 6 + 4] = - p.val[i * 6 + 4];
        }
		else if (x < 0 && y >= y_thr) {
			p.val[i * 6 + 0] = - p.val[i * 6 + 0];
            p.val[i * 6 + 3] = - p.val[i * 6 + 3];
		}

	}

	for (int i = n_remove - 1; i >= 0; i--)
	{
		remove_particle(p, n_active, tbremoved.val[i], lpos);
	}
}


void particle_operations::boundaries_e(fmatrix & p, int & n_active, imatrix & lpos, int n_out_i)
{

	int n_out = 0;
	static imatrix out(100000); 

	int n_out_ob = 0;
	static imatrix out_ob(100000); 

	bool in_thr, in_sym, is_crt;
	double energy, x, y, vx, vy, vz;
	const double x_max = ((double) n_mesh_x - 1);
	const double y_max = ((double) n_mesh_y - 1) * (dy / dx);
	const double y_thr = ((double) n_thruster - 1) * (dy / dx);
	
	for (int i = 0; i < n_active; i++)
	{
		x = p.val[i * 6 + 0];
		y = p.val[i * 6 + 1];

		if(x <= 0 || x >= x_max || y >= y_max || y <= 0){
			out.val[n_out] = i;
			n_out += 1;

			in_thr = (y <= y_thr) && (y > 0) && (x <= 0);
			in_sym = y <= 0;

			if(!in_thr && !in_sym){
				out_ob.val[n_out_ob] = i;
				n_out_ob += 1;
			}
		}
	}

	double e_crit = find_e_crit(n_out_i, out_ob, n_out_ob, p, n_active);
	for (int n = n_out - 1; n >= 0; n--){
		
		x = p.val[out.val[n] * 6 + 0];
		y = p.val[out.val[n] * 6 + 1];
		vx= p.val[out.val[n] * 6 + 3];
		vy= p.val[out.val[n] * 6 + 4];
		vz= p.val[out.val[n] * 6 + 5];
		
		energy = (vx * vx) + (vy * vy) + (vz * vz);

		is_crt = energy >= e_crit;
		in_thr = (y <= y_thr) && (y > 0) && (x <= 0);
		in_sym = (y <= 0) && (x >= 0) && (x <= x_max);

		if(in_sym || (!is_crt && !in_thr)){
			reflect_particle(p, n_active, out.val[n], x, y, vx, vy);
		} else {
			remove_particle(p, n_active, out.val[n], lpos);
		} 
	}
}

void particle_operations::boundaries_e_cap(fmatrix & p, int & n_active, imatrix & lpos, int & n_out_e, double v_cap, fmatrix & phi, mesh_set & mesh){
	n_out_e = 0;
	
	int n_out = 0;
	static imatrix out(100000); 

	static fmatrix phi_at_p(100000);

	double e_factor = 0.5 * m_el * dx * dx / (dt * dt * q);

	bool in_thr, in_sym, is_crt;
	double energy, x, y, vx, vy, vz;
	const double x_max = ((double) n_mesh_x - 1);
	const double y_max = ((double) n_mesh_y - 1) * (dy / dx);
	const double y_thr = ((double) n_thruster - 1) * (dy / dx);



	for (int i = 0; i < n_active; i++)
	{
		x = p.val[i * 6 + 0];
		y = p.val[i * 6 + 1];

		if(x <= 0 || x >= x_max || y >= y_max || y <= 0){
			out.val[n_out] = i;
			n_out += 1;
		}
	}

	find_phi_at_particles(phi_at_p, phi, mesh, out, n_out, p, n_active, lpos);
	for (int n = n_out - 1; n >= 0; n--){
		
		x = p.val[out.val[n] * 6 + 0];
		y = p.val[out.val[n] * 6 + 1];
		vx= p.val[out.val[n] * 6 + 3];
		vy= p.val[out.val[n] * 6 + 4];
		vz= p.val[out.val[n] * 6 + 5];
		
		energy =  e_factor * ((vx * vx) + (vy * vy) + (vz * vz));
		is_crt = energy >= phi_at_p.val[n] - v_cap;
		in_thr = (y <= y_thr) && (y > 0) && (x <= 0);
		in_sym = (y <= 0) && (x >= 0) && (x <= x_max);

		if(in_sym || (!is_crt && !in_thr)){
			reflect_particle(p, n_active, out.val[n], x, y, vx, vy);
		} else {
			remove_particle(p, n_active, out.val[n], lpos);
			if(!in_thr) n_out_e += 1;
		} 
	}
}

void particle_operations::find_phi_at_particles(fmatrix & phi_at_patricles, fmatrix & phi, mesh_set & mesh, imatrix & out, int n_out, fmatrix & p, int n_active, imatrix & lpos){
	
	const double x_max = ((double) n_mesh_x - 1);
	const double y_max = ((double) n_mesh_y - 1) * (dy / dx);
	
	double x, y;
	int lx, ly;

	for (int n = 0; n < n_out; n++)
	{	
		x = clamp_n(0.0, x_max, p.val[out.val[n] * 6 + 0]);
		y = clamp_n(0.0, y_max, p.val[out.val[n] * 6 + 1]);
		lx = lpos.val[out.val[n] * 2 + 0];
		ly = lpos.val[out.val[n] * 2 + 1];

		phi_at_patricles.val[n] = pic.field_at_position(phi, mesh, x, y, lx, ly) / k_phi;
	}
}

double particle_operations::cap_voltage(double voltage, int n_out_e, int n_out_i){
	return voltage + (n_factor * q / c_cap) * (n_out_i - n_out_e);
}

double particle_operations::find_e_crit(int n_out_i, imatrix & out, int n_out, fmatrix & p, int n_active){

	static fmatrix energy(10000);
	double vx, vy, vz;

	
	if(n_out_i >= n_out){
		return 0.0;
	}
	else if(n_out_i == 0){
		return 1e99;
	}
	else {
		for (int n = 0; n < n_out; n++)
		{
			vx =  p.val[out.val[n] * 6 + 3];
			vy =  p.val[out.val[n] * 6 + 4];
			vz =  p.val[out.val[n] * 6 + 5];
			energy.val[n] = (vx * vx) + (vy * vy) + (vz * vz);
		}
		nth_element(&energy.val[0], &energy.val[n_out - n_out_i], &energy.val[n_out]);
		return energy.val[n_out - n_out_i];
	}
}

void particle_operations::reflect_particle(fmatrix & p, int & n_active, int i, double x, double y, double vx, double vy)
{
	static const double x_max = ((double) n_mesh_x - 1);
	static const double y_max = ((double) n_mesh_y - 1) * (dy / dx);

	if(x <= 0){
		p.val[i * 6 + 0] = - x;
		p.val[i * 6 + 3] = - vx;
	} else if (x >= x_max) {
		p.val[i * 6 + 0] = 2*x_max - x;
		p.val[i * 6 + 3] = - vx;
	}

	if(y <= 0){
		p.val[i * 6 + 1] = - y;
		p.val[i * 6 + 4] = - vy;

	} else if(y >= y_max){
		p.val[i * 6 + 1] = 2*y_max - y;
		p.val[i * 6 + 4] = - vy;
	}
}


void particle_operations::remove_particle(fmatrix & p, int & n_active, int i, imatrix & lpos)
{
	// for (size_t j = 0; j < 6; j++)
	// {
	//   swap(p.val[i * 6 + j], p.val[(n_active - 1) * 6 + j]);
	// }

	// swap(lpos.val[i * 2 + 0], lpos.val[(n_active - 1) * 2 + 0]);
	// swap(lpos.val[i * 2 + 1], lpos.val[(n_active - 1) * 2 + 1]);

	p.val[i * 6 + 0] = p.val[(n_active - 1) * 6 + 0];
	p.val[i * 6 + 1] = p.val[(n_active - 1) * 6 + 1];
	p.val[i * 6 + 2] = p.val[(n_active - 1) * 6 + 2];
	p.val[i * 6 + 3] = p.val[(n_active - 1) * 6 + 3];
	p.val[i * 6 + 4] = p.val[(n_active - 1) * 6 + 4];
	p.val[i * 6 + 5] = p.val[(n_active - 1) * 6 + 5];

	lpos.val[i * 2 + 0] = lpos.val[(n_active - 1) * 2 + 0];
	lpos.val[i * 2 + 1] = lpos.val[(n_active - 1) * 2 + 1];

	n_active -= 1;
}

//  ----------------------------- Particle Movement ---------------------------

void particle_operations::move_e(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y)
{
    for (int i = 0; i < n_active; i++)
    {
        p.val[i * 6 + 3] = p.val[i * 6 + 3] - electric_field_at_particles_x.val[i];
		p.val[i * 6 + 4] = p.val[i * 6 + 4] - electric_field_at_particles_y.val[i];
        p.val[i * 6 + 0] = p.val[i * 6 + 0] + p.val[i * 6 + 3];
		p.val[i * 6 + 1] = p.val[i * 6 + 1] + p.val[i * 6 + 4];
    }
}

void particle_operations::move_i(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y)
{
    for (int i = 0; i < n_active; i++)
    {
        p.val[i * 6 + 3] = p.val[i * 6 + 3] + alpha * electric_field_at_particles_x.val[i];
		p.val[i * 6 + 4] = p.val[i * 6 + 4] + alpha * electric_field_at_particles_y.val[i];
        p.val[i * 6 + 0] = p.val[i * 6 + 0] + p.val[i * 6 + 3] * k_sub; // K_SUB must be multiplied here?
		p.val[i * 6 + 1] = p.val[i * 6 + 1] + p.val[i * 6 + 4] * k_sub; // K_SUB must be multiplied here?
    }
}

void particle_operations::move_n(fmatrix & p, int & n_active, double k_sub)
{
    for (int i = 0; i < n_active; i++)
    {
        p.val[i * 6 + 0] = p.val[i * 6 + 0] + k_sub * p.val[i * 6 + 3];
        p.val[i * 6 + 1] = p.val[i * 6 + 1] + k_sub * p.val[i * 6 + 4];
    }
}
