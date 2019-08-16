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

	const int mesh_n2 = (int) mesh_x.n2;;

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

		area_1 = (x_1_mesh - x_p) * (y_1_mesh - y_p);
		area_2 = (x_p - x_0_mesh) * (y_1_mesh - y_p);
		area_3 = (x_p - x_0_mesh) * (y_p - y_0_mesh);
		area_4 = (x_1_mesh - x_p) * (y_p - y_0_mesh);
		cell_area = area_1 + area_2 + area_3 + area_4;
		
		e_p_x =  efield_x.val[left_index_x * mesh_n2 + left_index_y] 				* area_1;
		e_p_x += efield_x.val[(left_index_x + 1) * mesh_n2 + left_index_y] 			* area_2;
		e_p_x += efield_x.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] 	* area_3;
		e_p_x += efield_x.val[left_index_x * mesh_n2 + (left_index_y + 1)] 			* area_4;

		efield_at_particles_x.val[i] = e_p_x / cell_area;

		e_p_y =  efield_y.val[left_index_x * mesh_n2 + left_index_y] 				* area_1;
		e_p_y += efield_y.val[(left_index_x + 1) * mesh_n2 + left_index_y] 			* area_2;
		e_p_y += efield_y.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] 	* area_3;
		e_p_y += efield_y.val[left_index_x * mesh_n2 + (left_index_y + 1)] 			* area_4;

		efield_at_particles_y.val[i] = e_p_y / cell_area;
	}
}

double field_at_position(fmatrix & field, fmatrix & mesh_x, fmatrix & mesh_y, double x, double y, int lpos_x, int lpos_y){

    const int mesh_n2 = (int) mesh_x.n2;

    double x_0_mesh = mesh_x.val[lpos_x * mesh_n2 + lpos_y];
	double x_1_mesh = mesh_x.val[(lpos_x + 1) * mesh_n2 + lpos_y];
    if(x < x_0_mesh || x > x_1_mesh)
	{
		lpos_x = logical_space(x, A_X, 1, N_MESH_X);
		x_0_mesh = mesh_x.val[lpos_x * mesh_n2 + lpos_y];
		x_1_mesh = mesh_x.val[(lpos_x + 1) * mesh_n2 + lpos_y];
	}
	
    double y_0_mesh = mesh_y.val[lpos_x * mesh_n2 + lpos_y];
	double y_1_mesh = mesh_y.val[lpos_x * mesh_n2 + (lpos_y + 1)];
	if(y < y_0_mesh || y > y_1_mesh)
	{
		lpos_y = logical_space(y, A_Y, DY/DX, N_MESH_Y);
		y_0_mesh = mesh_y.val[lpos_x * mesh_n2 + lpos_y];
		y_1_mesh = mesh_y.val[lpos_x * mesh_n2 + (lpos_y + 1)];
	}

    double field_val = field.val[lpos_x * mesh_n2 + lpos_y]             * (x_1_mesh - x) * (y_1_mesh - y);
    field_val       += field.val[(lpos_x + 1) * mesh_n2 + lpos_y]       * (x - x_0_mesh) * (y_1_mesh - y);
    field_val       += field.val[(lpos_x + 1) * mesh_n2 + (lpos_y + 1)] * (x - x_0_mesh) * (y - y_0_mesh);
    field_val       += field.val[lpos_x * mesh_n2 + (lpos_y + 1)]       * (x_1_mesh - x) * (y - y_0_mesh);

    return field_val / ((x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh));
}


void energy_field(fmatrix & kefield, fmatrix & p, int & n_active, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & wmesh, imatrix & lpos, double mass){
    
    kefield.set_zero();
    
    double x_p, y_p, vx, vy, vz;
    double kinetic_energy;
    
    double x_0_mesh = 0;
    double x_1_mesh = 0;
    double y_0_mesh = 0;
    double y_1_mesh = 0;
    
    double cell_area, area_1, area_2, area_3, area_4;
    
    int left_index_x = 0;
    int left_index_y = 0;
    
    const int mesh_n1 = (int) mesh_x.n1;
    const int mesh_n2 = (int) mesh_x.n2;
    const double factor = 0.5 * mass * (DX * DX) / (Q * DT * DT); // in eV
    
    for (int i = 0; i < n_active; i++)
    {
        x_p = p.val[i * 6 + 0];
        y_p = p.val[i * 6 + 1];
        vx =  p.val[i * 6 + 3];
        vy =  p.val[i * 6 + 4];
        vz =  p.val[i * 6 + 5];
        
        left_index_x = lpos.val[i * 2 + 0];
        left_index_y = lpos.val[i * 2 + 1];
        
        x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
        x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
        y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
        y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];
        
        area_1 = (x_1_mesh - x_p) * (y_1_mesh - y_p);
        area_2 = (x_p - x_0_mesh) * (y_1_mesh - y_p);
        area_3 = (x_p - x_0_mesh) * (y_p - y_0_mesh);
        area_4 = (x_1_mesh - x_p) * (y_p - y_0_mesh);
        cell_area = (x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh);
        
        kinetic_energy = vx*vx + vy*vy + vz*vz;
        
        kefield.val[left_index_x * mesh_n2 + left_index_y] +=              kinetic_energy * area_1 / cell_area;
        kefield.val[(left_index_x + 1) * mesh_n2 + left_index_y] +=        kinetic_energy * area_2 / cell_area;
        kefield.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  kinetic_energy * area_3 / cell_area;
        kefield.val[left_index_x * mesh_n2 + (left_index_y + 1)] +=        kinetic_energy * area_4 / cell_area;
    }
    
    for (int i = 0; i < mesh_n1 * mesh_n2; i++) {
        if (wmesh.val[i] > 0)
            kefield.val[i] = factor * kefield.val[i] / wmesh.val[i];
        else
            kefield.val[i] = 0;
    }
}

void flux_field(fmatrix & ffield_x, fmatrix & ffield_y, fmatrix & p, int & n_active, fmatrix & mesh_x, fmatrix & mesh_y, imatrix & lpos){
    
    ffield_x.set_zero();
    ffield_y.set_zero();
    
    double x_p, y_p, vx, vy;
    
    double x_0_mesh = 0;
    double x_1_mesh = 0;
    double y_0_mesh = 0;
    double y_1_mesh = 0;
    
    double cell_area, area_1, area_2, area_3, area_4;
    
    int left_index_x = 0;
    int left_index_y = 0;
    
    const int mesh_n2 = (int) mesh_x.n2;
    
    for (int i = 0; i < n_active; i++)
    {
        x_p = p.val[i * 6 + 0];
        y_p = p.val[i * 6 + 1];
        vx =  p.val[i * 6 + 3];
        vy =  p.val[i * 6 + 4];
        
        left_index_x = lpos.val[i * 2 + 0];
        left_index_y = lpos.val[i * 2 + 1];
        
        x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
        x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
        y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
        y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];
        
        area_1 = (x_1_mesh - x_p) * (y_1_mesh - y_p);
        area_2 = (x_p - x_0_mesh) * (y_1_mesh - y_p);
        area_3 = (x_p - x_0_mesh) * (y_p - y_0_mesh);
        area_4 = (x_1_mesh - x_p) * (y_p - y_0_mesh);
        cell_area = (x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh);
        
        ffield_x.val[left_index_x * mesh_n2 + left_index_y] +=              vx * area_1 / cell_area;
        ffield_x.val[(left_index_x + 1) * mesh_n2 + left_index_y] +=        vx * area_2 / cell_area;
        ffield_x.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  vx * area_3 / cell_area;
        ffield_x.val[left_index_x * mesh_n2 + (left_index_y + 1)] +=        vx * area_4 / cell_area;
        
        ffield_y.val[left_index_x * mesh_n2 + left_index_y] +=              vy * area_1 / cell_area;
        ffield_y.val[(left_index_x + 1) * mesh_n2 + left_index_y] +=        vy * area_2 / cell_area;
        ffield_y.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  vy * area_3 / cell_area;
        ffield_y.val[left_index_x * mesh_n2 + (left_index_y + 1)] +=        vy * area_4 / cell_area;
    }
}
