/*
 
 Functions to solve Poisson's equation. Steps are initialize tridiagonal matrix, solve via Thomas' algorithm
 and calculate electric field as derivative of the potential.
 
 Tridiagonal matrix is defined as a vector N_MESH x 4 that receives the parameters for each diagonal of the matrix
 and the vector d.
 
 */

#include "fields.h"

#include <cmath>

#include "fmatrix.h"
#include "config.h"
#include "num-tools.h"

/**
 Calculates a sinusoidal voltage level in time

 @param i step index
 @param dt time step
 @param freq_hz frequency in Hertz
 @param amplitude amplitude of the wave
 @param phase phase of the wave
 @return voltage at a given time
 */
double ac_voltage_at_time(size_t i, double dt, double freq_hz, double amplitude, double phase)
{
	return amplitude * sin(2 * PI * freq_hz * dt * i + phase);
}

/**
 Calculates the electric field at every node in the mesh in the X and Y directions. For the calculation it uses a CDS with charge correction on the Dirichlet boundaries.

 @param efield_x array holding the electric field in the X direction (output)
 @param efield_y array holding the electric field in the Y direction (output)
 @param phi array holding the electric potential
 @param w_i array holding the particle weights of ions
 @param w_e array holding the particle weights of electrons
 @param mesh_x mesh coordinates in the X direction
 @param mesh_y mesh coordinates in the Y direction
 @param vmesh volume mesh
 @todo still needs a correction for non-uniform meshes
 @todo correction for space charge
 */
void calculate_efield(fmatrix & efield_x, fmatrix & efield_y, fmatrix & phi, fmatrix & w_i, fmatrix & w_e, fmatrix & mesh_x , fmatrix & mesh_y, fmatrix & vmesh, imatrix & electrode_mask)
{
	int ip, im, jp, jm; //e_mask;
	int n_mesh_x =  (int) mesh_x.n1, n_mesh_y =  (int) mesh_y.n2;
	// double signal_x, signal_y;

	for(int i = 0; i < n_mesh_x; i++){
        for(int j = 0; j < n_mesh_y; j++)
        {

			ip = clamp_n(0, n_mesh_x - 1, i + 1);
			im = clamp_n(0, n_mesh_x - 1, i - 1);
			jp = clamp_n(0, n_mesh_y - 1, j + 1);
			jm = clamp_n(0, n_mesh_y - 1, j - 1);

			efield_x.val[i * n_mesh_y + j] = (phi.val[im * n_mesh_y + j] - phi.val[ip * n_mesh_y + j]) / (mesh_x.val[ip * n_mesh_y + j] - mesh_x.val[im * n_mesh_y + j]);
			efield_y.val[i * n_mesh_y + j] = (phi.val[i * n_mesh_y + jm] - phi.val[i * n_mesh_y + jp]) / (mesh_y.val[i * n_mesh_y + jp] - mesh_y.val[i * n_mesh_y + jm]);

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


double calculate_phi_zero(double sigma_old, double n_in, double q_cap, double sigma_laplace, fmatrix & phi_poisson, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, imatrix & electrode_mask){

	double sigma_poisson = sigma_from_phi(phi_poisson, mesh_x, mesh_y, wmesh_e, wmesh_i, vmesh, electrode_mask);	
	return (sigma_old + sigma_poisson - q_cap + K_Q * n_in) / (1 - sigma_laplace);
}

double sigma_from_phi(fmatrix & phi, fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & wmesh_e, fmatrix & wmesh_i, fmatrix & vmesh, imatrix & electrode_mask){
	
	double sigma = 0, volume, dw, phi_xx, phi_yy;
	double dx1, dx2, dy1, dy2;
	int ip, im, jp, jm;
	int n_mesh_x =  (int) mesh_x.n1, n_mesh_y =  (int) mesh_y.n2;
	
	for(int i = 0; i < n_mesh_x; i++){
		for(int j = 0; j < n_mesh_y; j++){
			if(electrode_mask.val[i * n_mesh_y + j] == 1){

				ip = clamp_n(0, n_mesh_x - 1, i + 1);
				im = clamp_n(0, n_mesh_x - 1, i - 1);
				jp = clamp_n(0, n_mesh_y - 1, j + 1);
				jm = clamp_n(0, n_mesh_y - 1, j - 1);

                dw = wmesh_i.val[i * n_mesh_y + j] - wmesh_e.val[i * n_mesh_y + j];

				dx1 = mesh_x.val[i * n_mesh_y + j] - mesh_x.val[im * n_mesh_y + j];
				dx2 = mesh_x.val[ip * n_mesh_y + j] - mesh_x.val[i * n_mesh_y + j];
				if (im == i || ip == i) dx1 = dx2 = dx1 + dx2;

				dy1 = mesh_y.val[i * n_mesh_y + j] - mesh_y.val[i * n_mesh_y + jm];
				dy2 = mesh_y.val[i * n_mesh_y + jp] - mesh_y.val[i * n_mesh_y + j];
				if (jm ==  j || jp == j) dy1 = dy2 = dy1 + dy2;

				phi_xx = phi.val[ip * n_mesh_y + j]/((dx1 + dx2)*dx2) - phi.val[i * n_mesh_y + j]/(dx1*dx2) + phi.val[im * n_mesh_y + j]/((dx1 + dx2)*dx1);
				phi_yy = phi.val[i * n_mesh_y + jp]/((dy1 + dy2)*dy2) - phi.val[i * n_mesh_y + j]/(dy1*dy2) + phi.val[i * n_mesh_y + jm]/((dy1 + dy2)*dy1);

				volume = vmesh.val[i * n_mesh_y + j];

                sigma += (EPS_0 / C_CAP) * volume * (phi_xx + phi_yy + GAMMA * dw / volume);
			}
		}
	}

	return sigma;
}

double calculate_sigma(double sigma_old, double phi_zero, double n_in, double cap_charge){
	return sigma_old - phi_zero - cap_charge + K_Q * n_in;
}

double calculate_cap_charge(double sigma_new, double sigma_old, double cap_charge_old, double n_in){
	return sigma_new - sigma_old + cap_charge_old - K_Q * n_in;
}


void init_mesh(fmatrix & mesh_x, fmatrix & mesh_y, double a_x, double a_y, int n_mesh_x, int n_mesh_y)
{
	for (size_t i = 0; i < mesh_x.n1; i++) {
		for (size_t j = 0; j < mesh_x.n2; j++) {
			
			mesh_x.val[i * mesh_x.n2 + j] = physical_space(i, a_x, 1, n_mesh_x);
			mesh_y.val[i * mesh_y.n2 + j] = physical_space(j, a_y, DY/DX, n_mesh_y);
		
		}	
	}
}


/**
 Initializes a mesh containing the values of the volume for each node.

 @param vmesh volume mesh array
 @param mesh_x array of x coordinates of the mesh
 @param mesh_y array of y coordinates of the mesh
 */
void init_volume_mesh(fmatrix & vmesh, fmatrix & mesh_x, fmatrix & mesh_y){
	size_t ip, im, jp, jm;
	for(size_t i = 0; i < vmesh.n1; i++){
		for(size_t j = 0; j < vmesh.n2; j++){
			
			ip = clamp_n(0, (int) mesh_x.n1 - 1, (int) i + 1);
			im = clamp_n(0, (int) mesh_x.n1 - 1, (int) i - 1);
			jp = clamp_n(0, (int) mesh_y.n2 - 1, (int) j + 1);
			jm = clamp_n(0, (int) mesh_y.n2 - 1, (int) j - 1);
		

			vmesh.val[i * vmesh.n2 + j] = (mesh_x.val[ip * mesh_x.n2 + j] - mesh_x.val[im * mesh_x.n2 + j]) * 
										  (mesh_y.val[i * mesh_y.n2 + jp] - mesh_y.val[i * mesh_y.n2 + jm]);
		}	
	}
}


double physical_space(double logical_position, double a, double b, double n_mesh){
	return b * ((a / (n_mesh - 1)) * pow(logical_position, 2) + (1 - a) * logical_position);
}


int logical_space(float physical_position, float a, float b, float n_mesh){
	return 2.0f * physical_position / (b * (1.0f - a + sqrtf((1.0f - a)*(1.0f - a) + (4.0f * a * physical_position / b)/(n_mesh - 1))));
}


double k1_x(fmatrix & mesh_x, size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return mesh_x.val[i * mesh_x.n2 + j] - mesh_x.val[(i - 1) * mesh_x.n2 + j];
}


double k2_x(fmatrix & mesh_x, size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return mesh_x.val[(i + 1) * mesh_x.n2 + j] - mesh_x.val[i * mesh_x.n2 + j];
}


double k3_x(fmatrix & mesh_x, size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return mesh_x.val[(i + 1) * mesh_x.n2 + j] - mesh_x.val[(i - 1) * mesh_x.n2 + j];
}


double k1_y(fmatrix & mesh_y, size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return mesh_y.val[i * mesh_y.n2 + j] - mesh_y.val[i * mesh_y.n2 + (j - 1)];
}


double k2_y(fmatrix & mesh_y, size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return mesh_y.val[i * mesh_y.n2 + (j + 1)] - mesh_y.val[i * mesh_y.n2 + j];
}


double k3_y(fmatrix & mesh_y, size_t i, size_t j, int ioff, int joff){
	i += ioff;
	j += joff;
	return mesh_y.val[i * mesh_y.n2 + (j + 1)] - mesh_y.val[i * mesh_y.n2 + (j - 1)];
}
