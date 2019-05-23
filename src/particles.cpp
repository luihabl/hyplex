/*
    Functions to add, move and remove particles.
 */


#include "particles.h"

#include <iostream>
#include <cmath>

#include "fmatrix.h"
#include "util.h"
#include "config.h"
#include "random-numbers.h"


using namespace std;

void remove_particle(fmatrix & p, int & n_active, int position, imatrix & lpos)
{
	for (size_t j = 0; j < 6; j++)
	{
		swap(p.val[position * 6 + j], p.val[(n_active - 1) * 6 + j]);
	}
	swap(lpos.val[position * 2 + 0], lpos.val[(n_active - 1) * 2 + 0]);
	swap(lpos.val[position * 2 + 1], lpos.val[(n_active - 1) * 2 + 1]);
	n_active -= 1;
}

void add_maxwellian_particles(fmatrix & p, int & n_active, const double temperature, const double mass, const size_t n_add)
{
	double v_temperature = sqrt(Q * temperature / mass);

	for (size_t i = n_active; i < n_active + n_add; i++)
	{
		p.val[i * 6 + 0] = ((double) N_MESH_X - 1.0) * r_unif();
		p.val[i * 6 + 1] = (DY / DX) * ((double) N_MESH_Y - 1.0) * r_unif();
		p.val[i * 6 + 2] = 0.0;
		p.val[i * 6 + 3] = (DT / DX) * r_norm(0.0, v_temperature);
		p.val[i * 6 + 4] = (DT / DX) * r_norm(0.0, v_temperature);
		p.val[i * 6 + 5] = (DT / DX) * r_norm(0.0, v_temperature);
	}
	n_active += n_add;
}

void add_flux_particles(fmatrix & p, int & n_active, const double temperature, const double v_drift, const double mass, const size_t n_add){

	double v_temperature = sqrt(Q * temperature / mass);

	for (size_t i = n_active; i < n_active + n_add; i++)
	{	
		p.val[i * 6 + 3] = (DT / DX) * (v_temperature * sqrt(- 2 * log(r_unif())) + v_drift);
		p.val[i * 6 + 4] = (DT / DX) * r_norm(0.0, v_temperature);
		p.val[i * 6 + 5] = (DT / DX) * r_norm(0.0, v_temperature);

		p.val[i * 6 + 0] = p.val[i * 6 + 3] * r_unif();
		p.val[i * 6 + 1] = (DY / DX) * ((double) N_THRUSTER - 1.0) * r_unif();
		p.val[i * 6 + 2] = 0.0;
	}
	n_active += n_add;
}

void add_maxwellian_particle_at_position(fmatrix & p, int & n_active, imatrix & lpos, const double temperature, const double mass, double x_pos, double y_pos, int lpos_x, int lpos_y)
{
    double v_temperature = sqrt(Q * temperature / mass);
    
    p.val[n_active * 6 + 0] = x_pos;
    p.val[n_active * 6 + 1] = y_pos;
    p.val[n_active * 6 + 2] = 0.0;
    p.val[n_active * 6 + 3] = (DT / DX) * r_norm(0.0, v_temperature);
    p.val[n_active * 6 + 4] = (DT / DX) * r_norm(0.0, v_temperature);
    p.val[n_active * 6 + 5] = (DT / DX) * r_norm(0.0, v_temperature);
	lpos.val[n_active * 2 + 0] = lpos_x;
	lpos.val[n_active * 2 + 1] = lpos_y;

    n_active += 1;
}

int add_particle_copy(fmatrix & p, int & n_active, imatrix & lpos, const int & i)
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

void boundaries(fmatrix & p, int & n_active, imatrix & lpos)
{
	int n_remove = 0;
	int tbremoved[100000]; // Modify fmatrix class to include int as template. Make a smart guess of the maximum removed particles.

	double pos_x = 0.0;
	double pos_y = 0.0;
	double x_max = ((double) N_MESH_X - 1);
	double y_max = ((double) N_MESH_Y - 1) * (DY / DX);

	for (int i = 0; i < n_active; i++)
	{
		pos_x = p.val[i * 6 + 0];
		pos_y = p.val[i * 6 + 1];

		if (pos_x <= 0 || pos_x >= x_max || pos_y >= y_max)
		{
			tbremoved[n_remove] = i;
			n_remove += 1;
		} else if (pos_y < 0) {
			p.val[i * 6 + 1] = - p.val[i * 6 + 1];
			p.val[i * 6 + 4] = - p.val[i * 6 + 4];
		}

		if (n_remove >= 100000)
		{
			printf("Out of bounds in tbremoved\n");
			break;
		}
	}

	for (int i = n_remove - 1; i >= 0; i--)
	{
		remove_particle(p, n_active, tbremoved[i], lpos);
	}

}

void move_e(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y)
{
    for (int i = 0; i < n_active; i++)
    {
        p.val[i * 6 + 3] = p.val[i * 6 + 3] - electric_field_at_particles_x.val[i];
		p.val[i * 6 + 4] = p.val[i * 6 + 4] - electric_field_at_particles_y.val[i];
        p.val[i * 6 + 0] = p.val[i * 6 + 0] + p.val[i * 6 + 3];
		p.val[i * 6 + 1] = p.val[i * 6 + 1] + p.val[i * 6 + 4];
    }
}

void move_i(fmatrix & p, int & n_active, fmatrix & electric_field_at_particles_x, fmatrix & electric_field_at_particles_y)
{
    for (int i = 0; i < n_active; i++)
    {
        p.val[i * 6 + 3] = p.val[i * 6 + 3] + ALPHA * electric_field_at_particles_x.val[i];
		p.val[i * 6 + 4] = p.val[i * 6 + 4] + ALPHA * electric_field_at_particles_y.val[i];
        p.val[i * 6 + 0] = p.val[i * 6 + 0] + p.val[i * 6 + 3] * K_SUB; // K_SUB must be multiplied here?
		p.val[i * 6 + 1] = p.val[i * 6 + 1] + p.val[i * 6 + 4] * K_SUB; // K_SUB must be multiplied here?
    }
}