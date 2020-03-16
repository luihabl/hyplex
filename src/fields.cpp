/*
 
 Functions to solve Poisson's equation. Steps are initialize tridiagonal matrix, solve via Thomas' algorithm
 and calculate electric field as derivative of the potential.
 
 Tridiagonal matrix is defined as a vector N_MESH x 4 that receives the parameters for each diagonal of the matrix
 and the vector d.
 
 */

#include "fields.h"

#include <cmath>

#include "fmatrix.h"
#include "num-tools.h"
#include "configuration.h"



mesh_set::mesh_set(configuration & config)
{
	nx = config.i("geometry/n_mesh_x");
	ny = config.i("geometry/n_mesh_y");
	nt = config.i("geometry/n_thruster");

	a_x = config.f("geometry/a_x"); 
	a_y = config.f("geometry/a_y"); 
	
	dx = config.f("geometry/dx");
	dy = config.f("geometry/dy"); 
	
	x = fmatrix::zeros(nx, ny);
	y = fmatrix::zeros(nx, ny);
	v = fmatrix::zeros(nx, ny); 
}


void mesh_set::init_mesh()
{
	// spatial mesh
	for (size_t i = 0; i < x.n1; i++) {
		for (size_t j = 0; j < x.n2; j++) {
			x.val[i * x.n2 + j] = physical_space(i, a_x, 1, nx);
			y.val[i * y.n2 + j] = physical_space(j, a_y, dy/dx, ny);
		}	
	}

	// volume mesh
	size_t ip, im, jp, jm;
	for(size_t i = 0; i < v.n1; i++){
		for(size_t j = 0; j < v.n2; j++){
			
			ip = clamp_n(0, (int) x.n1 - 1, (int) i + 1);
			im = clamp_n(0, (int) x.n1 - 1, (int) i - 1);
			jp = clamp_n(0, (int) y.n2 - 1, (int) j + 1);
			jm = clamp_n(0, (int) y.n2 - 1, (int) j - 1);
		

			v.val[i * v.n2 + j] = (x.val[ip * x.n2 + j] - x.val[im * x.n2 + j]) * 
										  (y.val[i * y.n2 + jp] - y.val[i * y.n2 + jm]);
		}	
	}
}



double mesh_set::k1_x(size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return x.val[i * x.n2 + j] - x.val[(i - 1) * x.n2 + j];
}


double mesh_set::k2_x(size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return x.val[(i + 1) * x.n2 + j] - x.val[i * x.n2 + j];
}


double mesh_set::k3_x(size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return x.val[(i + 1) * x.n2 + j] - x.val[(i - 1) * x.n2 + j];
}


double mesh_set::k1_y(size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return y.val[i * y.n2 + j] - y.val[i * y.n2 + (j - 1)];
}


double mesh_set::k2_y(size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return y.val[i * y.n2 + (j + 1)] - y.val[i * y.n2 + j];
}


double mesh_set::k3_y(size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return y.val[i * y.n2 + (j + 1)] - y.val[i * y.n2 + (j - 1)];
}


// ------------------------------------- field functions ---------------------------------------


field_operations::field_operations(configuration & config){

	n_mesh_x = config.i("geometry/n_mesh_x");
	n_mesh_y =  config.i("geometry/n_mesh_y");
	eps_0 = config.f("physical/eps_0");
	c_cap = config.f("boundaries/c_cap");
	gamma = config.f("p/gamma");
	k_q = config.f("p/k_q");
	pi = config.f("physical/pi");
	dt = config.f("time/dt");


}


double field_operations::ac_voltage_at_time(size_t i, double freq_hz, double amplitude, double phase)
{
	return amplitude * sin(2 * pi * freq_hz * dt * i + phase);
}

void field_operations::calculate_efield(fmatrix & efield_x, fmatrix & efield_y, fmatrix & phi, fmatrix & w_i, fmatrix & w_e, mesh_set & mesh, imatrix & electrode_mask)
{
	int ip, im, jp, jm; //e_mask;

	for(int i = 0; i < n_mesh_x; i++){
        for(int j = 0; j < n_mesh_y; j++)
        {

			ip = clamp_n(0, n_mesh_x - 1, i + 1);
			im = clamp_n(0, n_mesh_x - 1, i - 1);
			jp = clamp_n(0, n_mesh_y - 1, j + 1);
			jm = clamp_n(0, n_mesh_y - 1, j - 1);

			efield_x.val[i * n_mesh_y + j] = (phi.val[im * n_mesh_y + j] - phi.val[ip * n_mesh_y + j]) / (mesh.x.val[ip * n_mesh_y + j] - mesh.x.val[im * n_mesh_y + j]);
			efield_y.val[i * n_mesh_y + j] = (phi.val[i * n_mesh_y + jm] - phi.val[i * n_mesh_y + jp]) / (mesh.y.val[i * n_mesh_y + jp] - mesh.y.val[i * n_mesh_y + jm]);

			// Problem with correction of the electric field: gives number of electrons different from xoopic
			// e_mask = electrode_mask.val[i * n_mesh_y + j];
			// if(e_mask == 1 || e_mask == 2){
			// 	signal_x = 2 * i - ip - im;
			// 	signal_y = 2 * j - jp - jm;

			// 	efield_x.val[i * n_mesh_y + j] += signal_x * GAMMA * (w_i.val[i * n_mesh_y + j] - w_e.val[i * n_mesh_y + j]) * (mesh_x.val[ip * n_mesh_y + j] - mesh_x.val[im * n_mesh_y + j]) / vmesh.val[i * n_mesh_y + j];
			// 	efield_y.val[i * n_mesh_y + j] += signal_y * GAMMA * (w_i.val[i * n_mesh_y + j] - w_e.val[i * n_mesh_y + j]) * (mesh_y.val[i * n_mesh_y + jp] - mesh_y.val[i * n_mesh_y + jm]) / vmesh.val[i * n_mesh_y + j];
			// }
        }
    }
}


double field_operations::calculate_phi_zero(double sigma_old, double n_in, double q_cap, double sigma_laplace, fmatrix & phi_poisson, mesh_set mesh, fmatrix & wmesh_e, fmatrix & wmesh_i, imatrix & electrode_mask){

	double sigma_poisson = sigma_from_phi(phi_poisson, mesh, wmesh_e, wmesh_i, electrode_mask);	
	return (sigma_old + sigma_poisson - q_cap + k_q * n_in) / (1 - sigma_laplace);
}

double field_operations::sigma_from_phi(fmatrix & phi, mesh_set & mesh, fmatrix & wmesh_e, fmatrix & wmesh_i, imatrix & electrode_mask){
	
	double sigma = 0, volume, dw, phi_xx, phi_yy;
	double dx1, dx2, dy1, dy2;
	int ip, im, jp, jm;

	for(int i = 0; i < n_mesh_x; i++){
		for(int j = 0; j < n_mesh_y; j++){
			if(electrode_mask.val[i * n_mesh_y + j] == 1){

				ip = clamp_n(0, n_mesh_x - 1, i + 1);
				im = clamp_n(0, n_mesh_x - 1, i - 1);
				jp = clamp_n(0, n_mesh_y - 1, j + 1);
				jm = clamp_n(0, n_mesh_y - 1, j - 1);

                dw = wmesh_i.val[i * n_mesh_y + j] - wmesh_e.val[i * n_mesh_y + j];

				dx1 = mesh.x.val[i * n_mesh_y + j] - mesh.x.val[im * n_mesh_y + j];
				dx2 = mesh.x.val[ip * n_mesh_y + j] - mesh.x.val[i * n_mesh_y + j];
				if (im == i || ip == i) dx1 = dx2 = dx1 + dx2;

				dy1 = mesh.y.val[i * n_mesh_y + j] - mesh.y.val[i * n_mesh_y + jm];
				dy2 = mesh.y.val[i * n_mesh_y + jp] - mesh.y.val[i * n_mesh_y + j];
				if (jm ==  j || jp == j) dy1 = dy2 = dy1 + dy2;

				phi_xx = phi.val[ip * n_mesh_y + j]/((dx1 + dx2)*dx2) - phi.val[i * n_mesh_y + j]/(dx1*dx2) + phi.val[im * n_mesh_y + j]/((dx1 + dx2)*dx1);
				phi_yy = phi.val[i * n_mesh_y + jp]/((dy1 + dy2)*dy2) - phi.val[i * n_mesh_y + j]/(dy1*dy2) + phi.val[i * n_mesh_y + jm]/((dy1 + dy2)*dy1);

				volume = mesh.v.val[i * n_mesh_y + j];

                sigma += (eps_0 / c_cap) * volume * (phi_xx + phi_yy + gamma * dw / volume);
			}
		}
	}

	return sigma;
}

double field_operations::calculate_sigma(double sigma_old, double phi_zero, double n_in, double cap_charge){
	return sigma_old - phi_zero - cap_charge + k_q * n_in;
}

double field_operations::calculate_cap_charge(double sigma_new, double sigma_old, double cap_charge_old, double n_in){
	return sigma_new - sigma_old + cap_charge_old - k_q * n_in;
}


double physical_space(double logical_position, double a, double b, double n_mesh){
	return b * ((a / (n_mesh - 1)) * pow(logical_position, 2) + (1 - a) * logical_position);
}


int logical_space(float physical_position, float a, float b, float n_mesh){
	return 2.0f * physical_position / (b * (1.0f - a + sqrtf((1.0f - a)*(1.0f - a) + (4.0f * a * physical_position / b)/(n_mesh - 1))));
}


