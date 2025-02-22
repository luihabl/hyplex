/*
    Functions to solve do the particle weighting, calculation of density and
    find the electric field at each particle.
*/



#include "particles-in-mesh.h"

#include <cmath>

#include "fmatrix.h"
#include "fields.h"

#include <iostream>

using namespace std;


pic_operations::pic_operations(configuration & config){

	dx = config.f("geometry/dx");
	dy = config.f("geometry/dy");
	dt = config.f("time/dt");
	a_x = config.f("geometry/a_x");
	a_y = config.f("geometry/a_y");
	q = config.f("physical/q");

	n_mesh_x = config.i("geometry/n_mesh_x");
	n_mesh_y = config.i("geometry/n_mesh_y");

}


void pic_operations::weight(fmatrix & p, int & n_active, fmatrix & wmesh, mesh_set & mesh, imatrix & lpos)
{
	wmesh.set_zero();

	fmatrix & mesh_x = mesh.x;
	fmatrix & mesh_y = mesh.y;

	double x_p = 0;
	double y_p = 0;

	double x_0_mesh = 0;
	double x_1_mesh = 0;
	double y_0_mesh = 0;
	double y_1_mesh = 0;
    
    double cell_area = 0;

	int left_index_x = 0;
	int left_index_y = 0;

	const int mesh_n1 = (int) mesh_x.n1;
    const int mesh_n2 = (int) mesh_x.n2;
    const double a_x = mesh.a_x;
	const double a_y = mesh.a_x;
	const double dx = mesh.dx;
	const double dy = mesh.dy;

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
			lpos.val[i * 2 + 0] = left_index_x = logical_space(x_p, a_x, 1, mesh_n1);
			x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
			x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
		}
		
		if(y_p < y_0_mesh || y_p > y_1_mesh)
		{
			lpos.val[i * 2 + 1] = left_index_y = logical_space(y_p, a_y, dy/dx, mesh_n2);
			y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
			y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];
		}
        
        cell_area = (x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh);
		wmesh.val[left_index_x * mesh_n2 + left_index_y] += 			 (x_1_mesh - x_p) * (y_1_mesh - y_p) / cell_area;
		wmesh.val[(left_index_x + 1) * mesh_n2 + left_index_y] += 		 (x_p - x_0_mesh) * (y_1_mesh - y_p) / cell_area;
		wmesh.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  (x_p - x_0_mesh) * (y_p - y_0_mesh) / cell_area;
		wmesh.val[left_index_x * mesh_n2 + (left_index_y + 1)] += 		 (x_1_mesh - x_p) * (y_p - y_0_mesh) / cell_area;
    }
}

void pic_operations::electric_field_at_particles(fmatrix & efield_at_particles_x, fmatrix & efield_at_particles_y, fmatrix & efield_x, fmatrix & efield_y, fmatrix & p, const int n_active, mesh_set & mesh, imatrix & lpos)
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

	size_t left_index_x = 0;
	size_t left_index_y = 0;

	size_t ip, jp, ipjp, ij; 

    const size_t mesh_n2 = mesh.x.n2;

	for (int i = 0; i < n_active; i++)
	{
		x_p = p.val[i * 6 + 0];
		y_p = p.val[i * 6 + 1];

		left_index_x = lpos.val[i * 2 + 0];
		left_index_y = lpos.val[i * 2 + 1];

		ij = left_index_x * mesh_n2 + left_index_y;
		ip = (left_index_x + 1) * mesh_n2 + left_index_y;
		jp = left_index_x * mesh_n2 + (left_index_y + 1);
		ipjp = (left_index_x + 1) * mesh_n2 + (left_index_y + 1);

		x_0_mesh = mesh.x.val[ij];
		x_1_mesh = mesh.x.val[ip];
		y_0_mesh = mesh.y.val[ij];
		y_1_mesh = mesh.y.val[jp];

		area_1 = (x_1_mesh - x_p) * (y_1_mesh - y_p);
		area_2 = (x_p - x_0_mesh) * (y_1_mesh - y_p);
		area_3 = (x_p - x_0_mesh) * (y_p - y_0_mesh);
		area_4 = (x_1_mesh - x_p) * (y_p - y_0_mesh);
		cell_area = area_1 + area_2 + area_3 + area_4;

		
		efield_at_particles_x.val[i]  =  ((efield_x.val[ij] 	* area_1)
										+ (efield_x.val[ip] 	* area_2)
										+ (efield_x.val[ipjp] 	* area_3)
										+ (efield_x.val[jp] 	* area_4)) / cell_area;

		efield_at_particles_y.val[i]  =  ((efield_y.val[ij] 	* area_1)
										+ (efield_y.val[ip] 	* area_2)
										+ (efield_y.val[ipjp] 	* area_3)
										+ (efield_y.val[jp] 	* area_4)) / cell_area;
	}
}

double pic_operations::field_at_position(fmatrix & field, mesh_set & mesh, double x, double y, int lpos_x, int lpos_y)
{

    const int mesh_n2 = (int) mesh.x.n2;
	const int n_mesh_x = (int) mesh.nx;
	const int n_mesh_y = (int) mesh.ny;

    double x_0_mesh = mesh.x.val[lpos_x * mesh_n2 + lpos_y];
	double x_1_mesh = mesh.x.val[(lpos_x + 1) * mesh_n2 + lpos_y];
    if(x < x_0_mesh || x > x_1_mesh)
	{
		lpos_x = logical_space(x, a_x, 1, n_mesh_x);
		x_0_mesh = mesh.x.val[lpos_x * mesh_n2 + lpos_y];
		x_1_mesh = mesh.x.val[(lpos_x + 1) * mesh_n2 + lpos_y];
	}
	
    double y_0_mesh = mesh.y.val[lpos_x * mesh_n2 + lpos_y];
	double y_1_mesh = mesh.y.val[lpos_x * mesh_n2 + (lpos_y + 1)];
	if(y < y_0_mesh || y > y_1_mesh)
	{
		lpos_y = logical_space(y, a_y, dy/dx, n_mesh_y);
		y_0_mesh = mesh.y.val[lpos_x * mesh_n2 + lpos_y];
		y_1_mesh = mesh.y.val[lpos_x * mesh_n2 + (lpos_y + 1)];
	}

    double field_val = field.val[lpos_x * mesh_n2 + lpos_y]             * (x_1_mesh - x) * (y_1_mesh - y);
    field_val       += field.val[(lpos_x + 1) * mesh_n2 + lpos_y]       * (x - x_0_mesh) * (y_1_mesh - y);
    field_val       += field.val[(lpos_x + 1) * mesh_n2 + (lpos_y + 1)] * (x - x_0_mesh) * (y - y_0_mesh);
    field_val       += field.val[lpos_x * mesh_n2 + (lpos_y + 1)]       * (x_1_mesh - x) * (y - y_0_mesh);

    return field_val / ((x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh));
}


void pic_operations::kfield(fmatrix & p, int & n_active, fmatrix & kefield, mesh_set & mesh, imatrix & lpos)
{
    
    kefield.set_zero();

	fmatrix & mesh_x = mesh.x;
	fmatrix & mesh_y = mesh.y;
    
    double x_p, y_p, vx, vy, vz;
    double ke;
    
    double x_0_mesh = 0;
    double x_1_mesh = 0;
    double y_0_mesh = 0;
    double y_1_mesh = 0;
    
    double cell_area;
    
    int left_index_x = 0;
    int left_index_y = 0;
    
    // const int mesh_n1 = (int) mesh_x.n1;
    // const double factor = 0.5 * mass * (dx * dx) / (q * dt * dt); // in eV
	const int mesh_n1 = (int) mesh_x.n1;
    const int mesh_n2 = (int) mesh_x.n2;
    const double a_x = mesh.a_x;
	const double a_y = mesh.a_x;
	const double dx = mesh.dx;
	const double dy = mesh.dy;

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

		if(x_p < x_0_mesh || x_p > x_1_mesh)
		{
			lpos.val[i * 2 + 0] = left_index_x = logical_space(x_p, a_x, 1, mesh_n1);
			x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
			x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
		}
		
		if(y_p < y_0_mesh || y_p > y_1_mesh)
		{
			lpos.val[i * 2 + 1] = left_index_y = logical_space(y_p, a_y, dy/dx, mesh_n2);
			y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
			y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];
		}
        
        cell_area = (x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh);
        
        ke = vx*vx + vy*vy + vz*vz;
        
        kefield.val[left_index_x * mesh_n2 + left_index_y] +=              ke * (x_1_mesh - x_p) * (y_1_mesh - y_p) / cell_area;
        kefield.val[(left_index_x + 1) * mesh_n2 + left_index_y] +=        ke * (x_p - x_0_mesh) * (y_1_mesh - y_p) / cell_area;
        kefield.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  ke * (x_p - x_0_mesh) * (y_p - y_0_mesh) / cell_area;
        kefield.val[left_index_x * mesh_n2 + (left_index_y + 1)] +=        ke * (x_1_mesh - x_p) * (y_p - y_0_mesh) / cell_area;
    }
}

void pic_operations::ufield(fmatrix & p, int & n_active, fmatrix & ufield_x, fmatrix & ufield_y, mesh_set & mesh,imatrix & lpos){
    
    ufield_x.set_zero();
    ufield_y.set_zero();

	fmatrix & mesh_x = mesh.x;
	fmatrix & mesh_y = mesh.y;
    
    double vx = 0; 
	double vy = 0;
    
	double x_p = 0;
	double y_p = 0;

	double x_0_mesh = 0;
	double x_1_mesh = 0;
	double y_0_mesh = 0;
	double y_1_mesh = 0;
    
    double cell_area = 0;

	int left_index_x = 0;
	int left_index_y = 0;

	double w1, w2, w3, w4;
    
	const int mesh_n1 = (int) mesh_x.n1;
    const int mesh_n2 = (int) mesh_x.n2;
    const double a_x = mesh.a_x;
	const double a_y = mesh.a_x;
	const double dx = mesh.dx;
	const double dy = mesh.dy;

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

		if(x_p < x_0_mesh || x_p > x_1_mesh)
		{
			lpos.val[i * 2 + 0] = left_index_x = logical_space(x_p, a_x, 1, mesh_n1);
			x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
			x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
		}
		
		if(y_p < y_0_mesh || y_p > y_1_mesh)
		{
			lpos.val[i * 2 + 1] = left_index_y = logical_space(y_p, a_y, dy/dx, mesh_n2);
			y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
			y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];
		}
        
		cell_area = (x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh);
        w1 = (x_1_mesh - x_p) * (y_1_mesh - y_p) / cell_area;
        w2 = (x_p - x_0_mesh) * (y_1_mesh - y_p) / cell_area;
        w3 = (x_p - x_0_mesh) * (y_p - y_0_mesh) / cell_area;
        w4 = (x_1_mesh - x_p) * (y_p - y_0_mesh) / cell_area;
        
        
        ufield_x.val[left_index_x * mesh_n2 + left_index_y] +=              vx * w1;
        ufield_x.val[(left_index_x + 1) * mesh_n2 + left_index_y] +=        vx * w2;
        ufield_x.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  vx * w3;
        ufield_x.val[left_index_x * mesh_n2 + (left_index_y + 1)] +=        vx * w4;
        
        ufield_y.val[left_index_x * mesh_n2 + left_index_y] +=              vy * w1;
        ufield_y.val[(left_index_x + 1) * mesh_n2 + left_index_y] +=        vy * w2;
        ufield_y.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  vy * w3;
        ufield_y.val[left_index_x * mesh_n2 + (left_index_y + 1)] +=        vy * w4;
    }
		
}

void pic_operations::weight_n(fmatrix & p, int & n_active_start, int & n_active_end, fmatrix & f, mesh_set & mesh, imatrix & lpos)
{
	fmatrix & mesh_x = mesh.x;
	fmatrix & mesh_y = mesh.y;
    
    double x_p, y_p;
    
    double x_0_mesh = 0;
    double x_1_mesh = 0;
    double y_0_mesh = 0;
    double y_1_mesh = 0;
    
    double cell_area;
    
    int left_index_x = 0;
    int left_index_y = 0;
    

	const int mesh_n1 = (int) mesh_x.n1;
    const int mesh_n2 = (int) mesh_x.n2;
    const double a_x = mesh.a_x;
	const double a_y = mesh.a_x;
	const double dx = mesh.dx;
	const double dy = mesh.dy;

    for (int i = n_active_start; i < n_active_end; i++)
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
			lpos.val[i * 2 + 0] = left_index_x = logical_space(x_p, a_x, 1, mesh_n1);
			x_0_mesh = mesh_x.val[left_index_x * mesh_n2 + left_index_y];
			x_1_mesh = mesh_x.val[(left_index_x + 1) * mesh_n2 + left_index_y];
		}
		
		if(y_p < y_0_mesh || y_p > y_1_mesh)
		{
			lpos.val[i * 2 + 1] = left_index_y = logical_space(y_p, a_y, dy/dx, mesh_n2);
			y_0_mesh = mesh_y.val[left_index_x * mesh_n2 + left_index_y];
			y_1_mesh = mesh_y.val[left_index_x * mesh_n2 + (left_index_y + 1)];
		}
        
        cell_area = (x_1_mesh - x_0_mesh) * (y_1_mesh - y_0_mesh);
        
        f.val[left_index_x * mesh_n2 + left_index_y] +=              (x_1_mesh - x_p) * (y_1_mesh - y_p) / cell_area;
        f.val[(left_index_x + 1) * mesh_n2 + left_index_y] +=        (x_p - x_0_mesh) * (y_1_mesh - y_p) / cell_area;
        f.val[(left_index_x + 1) * mesh_n2 + (left_index_y + 1)] +=  (x_p - x_0_mesh) * (y_p - y_0_mesh) / cell_area;
        f.val[left_index_x * mesh_n2 + (left_index_y + 1)] +=        (x_1_mesh - x_p) * (y_p - y_0_mesh) / cell_area;
    }
}