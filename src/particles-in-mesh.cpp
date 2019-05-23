/*
    Functions to solve do the particle weighting, calculation of density and
    find the electric field at each particle.
*/



#include "particles-in-mesh.h"

#include <cmath>

#include "fmatrix.h"
#include "fields.h"
#include "config.h"

#include <iostream>

using namespace std;

void weight(fmatrix & p, int & n_active, fmatrix & wmesh, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos)
{
	wmesh.set_zero();

	double x_p = 0;
	double y_p = 0;

	double x_0_mesh = 0;
	double x_1_mesh = 0;
	double y_0_mesh = 0;
	double y_1_mesh = 0;

	int left_index_x = 0;
	int left_index_y = 0;

	const int mesh_n2 = (int) mesh_x.n2;

	for (int i = 0; i < n_active; i++)
	{
		x_p = p.val[i * 6 + 0];
		y_p = p.val[i * 6 + 1];

		left_index_x = lpos.val[i * 2 + 0];
		left_index_y = lpos.val[i * 2 + 1];

		x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
		x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
		y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
		y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];

		if(x_p < x_0_mesh || x_p > x_1_mesh)
		{
			lpos.val[i * 2 + 0] = left_index_x = logical_space(x_p, A_X, 1, N_MESH_X);
			x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
			x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
		}
		
		if(y_p < y_0_mesh || y_p > y_1_mesh)
		{
			lpos.val[i * 2 + 1] = left_index_y = logical_space(y_p, A_Y, DY/DX, N_MESH_Y);
			y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
			y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];
		}
		
		wmesh.val[left_index_x * mesh_n2 + left_index_y] += 			 (x_1_mesh - x_p) * (y_1_mesh - y_p) / ((x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh));
		wmesh.val[(left_index_x + 1) * mesh_n2 + left_index_y] += 		 (x_p - x_0_mesh) * (y_1_mesh - y_p) / ((x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh));
		wmesh.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  (x_p - x_0_mesh) * (y_p - y_0_mesh) / ((x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh));
		wmesh.val[left_index_x * mesh_n2 + (left_index_y + 1)] += 		 (x_1_mesh - x_p) * (y_p - y_0_mesh) / ((x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh));
	}
}

void electric_field_at_particles(fmatrix & efield_at_particles_x, fmatrix & efield_at_particles_y, fmatrix & efield_x, fmatrix & efield_y, fmatrix & p, const int n_active, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos)
{
	double x_p = 0;
	double y_p = 0;

	double x_0_mesh = 0;
	double x_1_mesh = 0;
	double y_0_mesh = 0;
	double y_1_mesh = 0;

	double cell_area = 0;
	double area_1 = 0;
	double area_2 = 0;
	double area_3 = 0;
	double area_4 = 0;

	int left_index_x = 0;
	int left_index_y = 0;

	double e_p_x = 0.0;
	double e_p_y = 0.0;

	for (int i = 0; i < n_active; i++)
	{
		x_p = p.val[i * 6 + 0];
		y_p = p.val[i * 6 + 1];

		left_index_x = lpos.val[i * 2 + 0];
		left_index_y = lpos.val[i * 2 + 1];

		x_0_mesh = mesh_x.val[left_index_x * mesh_x.n2 + left_index_y];
		x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_x.n2 + left_index_y];
		y_0_mesh = mesh_y.val[left_index_x * mesh_y.n2 + left_index_y];
		y_1_mesh = mesh_y.val[left_index_x * mesh_y.n2 + (left_index_y + 1)];

		area_1 = (x_1_mesh - x_p) * (y_1_mesh - y_p);
		area_2 = (x_p - x_0_mesh) * (y_1_mesh - y_p);
		area_3 = (x_p - x_0_mesh) * (y_p - y_0_mesh);
		area_4 = (x_1_mesh - x_p) * (y_p - y_0_mesh);
		cell_area = area_1 + area_2 + area_3 + area_4;
		
		e_p_x = efield_x.val[left_index_x * efield_x.n2 + left_index_y] 				* area_1;
		e_p_x += efield_x.val[(left_index_x + 1) * efield_x.n2 + left_index_y] 			* area_2;
		e_p_x += efield_x.val[(left_index_x + 1) * efield_x.n2 + (left_index_y + 1)] 	* area_3;
		e_p_x += efield_x.val[left_index_x * efield_x.n2 + (left_index_y + 1)] 			* area_4;

		efield_at_particles_x.val[i] = e_p_x / cell_area;

		e_p_y = efield_y.val[left_index_x * efield_y.n2 + left_index_y] 				* area_1;
		e_p_y += efield_y.val[(left_index_x + 1) * efield_y.n2 + left_index_y] 			* area_2;
		e_p_y += efield_y.val[(left_index_x + 1) * efield_y.n2 + (left_index_y + 1)] 	* area_3;
		e_p_y += efield_y.val[left_index_x * efield_y.n2 + (left_index_y + 1)] 			* area_4;

		efield_at_particles_y.val[i] = e_p_y / cell_area;
	}
}
