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
#include "util.h"

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
 */
void calculate_efield(fmatrix & efield_x, fmatrix & efield_y, fmatrix & phi, fmatrix & w_i, fmatrix & w_e, fmatrix & mesh_x , fmatrix & mesh_y, fmatrix & vmesh)
{
	/*
    for(size_t i = 0; i < phi.n1; i++){
        for(size_t j = 0; j < phi.n2; j++)
        {
            if (i == 0)
                efield_x.val[i * efield_x.n2 + j] = (phi.val[i * phi.n2 + j] - phi.val[(i + 1) * phi.n2 + j]) / (mesh_x.val[(i + 1) * mesh_x.n2 + j] - mesh_x.val[i * mesh_x.n2 + j]) - GAMMA * (w_i.val[i * w_i.n2 + j] - w_e.val[i * w_e.n2 + j]) * (mesh_x.val[(i + 1) * mesh_x.n2 + j] - mesh_x.val[i * mesh_x.n2 + j]) / vmesh.val[i * vmesh.n2 + j];
            else if(i == phi.n1 - 1)
                efield_x.val[i * efield_x.n2 + j] = (phi.val[(i - 1) * phi.n2 + j] - phi.val[i * phi.n2 + j]) / (mesh_x.val[i * mesh_x.n2 + j] - mesh_x.val[(i - 1) * mesh_x.n2 + j]) + GAMMA * (w_i.val[i * w_i.n2 + j] - w_e.val[i * w_e.n2 + j]) * (mesh_x.val[i * mesh_x.n2 + j] - mesh_x.val[(i - 1) * mesh_x.n2 + j]) / vmesh.val[i * vmesh.n2 + j];
            else
                efield_x.val[i * efield_x.n2 + j] = (phi.val[(i - 1) * phi.n2 + j] - phi.val[(i + 1) * phi.n2 + j]) / (mesh_x.val[(i + 1) * mesh_x.n2 + j] - mesh_x.val[(i - 1) * mesh_x.n2 + j]);
            
             if (j == 0)
                efield_y.val[i * efield_y.n2 + j] = (phi.val[i * phi.n2 + j] - phi.val[i * phi.n2 + (j + 1)]) / (mesh_y.val[i * mesh_y.n2 + (j + 1)] - mesh_y.val[i * mesh_y.n2 + j]);
            else if(j == phi.n2 - 1)
                efield_y.val[i * efield_y.n2 + j] = (phi.val[i * phi.n2 + (j - 1)] - phi.val[i * phi.n2 + j]) / (mesh_y.val[i * mesh_y.n2 + j] - mesh_y.val[i * mesh_y.n2 + (j - 1)]);
            else
                efield_y.val[i * efield_y.n2 + j] = (phi.val[i * phi.n2 + (j - 1)] - phi.val[i * phi.n2 + (j + 1)]) / (mesh_y.val[i * mesh_y.n2 + (j + 1)] - mesh_y.val[i * mesh_y.n2 + (j - 1)]);
        }
    }
	*/

	for(size_t i = 0; i < phi.n1; i++){
        for(size_t j = 0; j < phi.n2; j++)
        {
            if (i == 0)
                efield_x.val[i * efield_x.n2 + j] = (phi.val[i * phi.n2 + j] - phi.val[(i + 1) * phi.n2 + j]) / (mesh_x.val[(i + 1) * mesh_x.n2 + j] - mesh_x.val[i * mesh_x.n2 + j]);
            else if(i == phi.n1 - 1)
                efield_x.val[i * efield_x.n2 + j] = (phi.val[(i - 1) * phi.n2 + j] - phi.val[i * phi.n2 + j]) / (mesh_x.val[i * mesh_x.n2 + j] - mesh_x.val[(i - 1) * mesh_x.n2 + j]);
            else
                efield_x.val[i * efield_x.n2 + j] = (phi.val[(i - 1) * phi.n2 + j] - phi.val[(i + 1) * phi.n2 + j]) / (mesh_x.val[(i + 1) * mesh_x.n2 + j] - mesh_x.val[(i - 1) * mesh_x.n2 + j]);
            
            if (j == 0)
                efield_y.val[i * efield_y.n2 + j] = (phi.val[i * phi.n2 + j] - phi.val[i * phi.n2 + (j + 1)]) / (mesh_y.val[i * mesh_y.n2 + (j + 1)] - mesh_y.val[i * mesh_y.n2 + j]);
            else if(j == phi.n2 - 1)
                efield_y.val[i * efield_y.n2 + j] = (phi.val[i * phi.n2 + (j - 1)] - phi.val[i * phi.n2 + j]) / (mesh_y.val[i * mesh_y.n2 + j] - mesh_y.val[i * mesh_y.n2 + (j - 1)]);
            else
                efield_y.val[i * efield_y.n2 + j] = (phi.val[i * phi.n2 + (j - 1)] - phi.val[i * phi.n2 + (j + 1)]) / (mesh_y.val[i * mesh_y.n2 + (j + 1)] - mesh_y.val[i * mesh_y.n2 + (j - 1)]);
        }
    }

	
	// ADD HERE ELECTRIC FIELD CORRECTION <--------------------------



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
			
			ip = clamp(0, (int) mesh_x.n1 - 1, (int) i + 1);
			im = clamp(0, (int) mesh_x.n1 - 1, (int) i - 1);
			jp = clamp(0, (int) mesh_y.n2 - 1, (int) j + 1);
			jm = clamp(0, (int) mesh_y.n2 - 1, (int) j - 1);
		

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
