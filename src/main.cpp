
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "cross-sections.h"
#include "config.h"
#include "fields.h"
#include "fmatrix.h"
#include "fmath.h"
#include "mcc.h"
#include "mpi.h"
#include "particles.h"
#include "particles-in-mesh.h"
#include "rsolver.h"
#include "util.h"

#define CONFIG_PATH "src/config/config.ini"

using namespace std;
using namespace std::chrono;

// ----------------------------- Main function --------------------------------
int main(int argc, char* argv[])
{
    load_config_file(argv[1] == NULL ? CONFIG_PATH : argv[1]);
    load_cross_sections();

    // ----------------------------- Initialization ---------------------------

	// General field variables
    verbose_log("Initializing general variables");
	fmatrix phi             = fmatrix::zeros(N_MESH_X, N_MESH_Y);
	fmatrix efield_x        = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_y        = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix mesh_x          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix mesh_y          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix vmesh           = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix voltages        = fmatrix::zeros(3);
    
    // Average field variablesd
    verbose_log("Initializing average field variables");
    fmatrix phi_av          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_av_x     = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_av_y     = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix dens_e_av       = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix dens_i_av       = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    
	// Particle 1 - Electrons
    verbose_log("Initializing electrons variables");
	int n_active_e          = 0;
	fmatrix p_e             = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_e          = imatrix::zeros(N_MAX_PARTICLES, 2);
	fmatrix wmesh_e         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_x_at_p_e = fmatrix::zeros(N_MAX_PARTICLES);
    fmatrix efield_y_at_p_e = fmatrix::zeros(N_MAX_PARTICLES);
	double nu_prime_e       = find_nu_prime_e(N_NEUTRAL);
	double p_null_e         = p_null(nu_prime_e, DT);
    double n_inj_balanced_e = N_INJ_EL;
	
	// Particle 2 - Ions
    verbose_log("Initializing ions variables");
	int n_active_i          = 0;
	fmatrix p_i             = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_i          = imatrix::zeros(N_MAX_PARTICLES, 2);
	fmatrix wmesh_i         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
	fmatrix efield_x_at_p_i = fmatrix::zeros(N_MAX_PARTICLES);
    fmatrix efield_y_at_p_i = fmatrix::zeros(N_MAX_PARTICLES);
	double nu_prime_i       = find_nu_prime_i(N_NEUTRAL);
	double p_null_i         = p_null(nu_prime_i, DT * K_SUB);
    int n_out_i             = 0;

    //initializing mesh
    verbose_log("Initializing mesh");
    init_mesh(mesh_x,  mesh_y, A_X, A_Y, N_MESH_X, N_MESH_Y);
    init_volume_mesh(vmesh, mesh_x, mesh_y);

    // Initializing solver 
    verbose_log("Initializing solver");
    rsolver solver(mesh_x, mesh_y, vmesh, 4, 1);

    imatrix box_thruster = {0, 0, 0, N_THRUSTER - 1};
    
    solver.set_dirichlet_box(box_thruster, 0);
    
    imatrix box_1 = {0, N_THRUSTER, 0, N_MESH_Y - 2};
    imatrix box_2 = {0, N_MESH_Y - 1, N_MESH_X - 2, N_MESH_Y - 1};
    imatrix box_3 = {N_MESH_X - 1, 0, N_MESH_X - 1, N_MESH_Y - 1};
    imatrix box_4 = {1, 0, N_MESH_X - 2, 0};
    
    solver.set_neumann_box(box_1, 0);
    solver.set_neumann_box(box_2, 1);
    solver.set_neumann_box(box_3, 2);
    solver.set_neumann_box(box_4, 3);
    
    voltages = {VOLT_0_NORM};
    solver.assemble();

    // Printing initial information
    print_initial_info(p_null_e, p_null_i);
       
	// ----------------------------- Main loop --------------------------------
    verbose_log(" ---- Starting main loop ---- ");
    auto start = high_resolution_clock::now();
	for (int i = 0; i < N_STEPS; i++)
	{   
		// Step 1.0: particle weighting
        if(i % K_SUB == 0) weight(p_i, n_active_i, wmesh_i, mesh_x, mesh_y, lpos_i);
        weight(p_e, n_active_e, wmesh_e, mesh_x, mesh_y, lpos_e);
        
        // Step 2.0 integration of Poisson's equation
        solver.solve(phi, voltages, wmesh_i, wmesh_e);
        
//        if(hasnan(phi)){
//            cout << "phi has nan, step: " << i << endl;
//            break;
//        }

        // Step 2.1: calculation of electric field
        calculate_efield(efield_x, efield_y, phi, wmesh_i, wmesh_e, mesh_x, mesh_y, vmesh);

        // Step 2.2: field weighting
        if(i % K_SUB == 0) electric_field_at_particles(efield_x_at_p_i, efield_y_at_p_i, efield_x, efield_y, p_i, n_active_i, mesh_x, mesh_y, lpos_i);
        electric_field_at_particles(efield_x_at_p_e, efield_y_at_p_e, efield_x, efield_y, p_e, n_active_e, mesh_x, mesh_y, lpos_e);

        
        // Step 3: particles injection
        if(i % K_SUB == 0) add_flux_particles(p_i, n_active_i, T_I, VD_I, M_I, N_INJ_I, K_SUB);
        n_inj_balanced_e = balanced_injection(n_inj_balanced_e, 0.01, wmesh_i, wmesh_e, 0, 0, 0, N_THRUSTER - 1);
        add_flux_particles(p_e, n_active_e, T_EL, 0, M_EL, n_inj_balanced_e);
        
        // Step 4: integration of equations of motion
        if(i % K_SUB == 0) move_i(p_i, n_active_i, efield_x_at_p_i, efield_y_at_p_i);
        move_e(p_e, n_active_e, efield_x_at_p_e, efield_y_at_p_e);
        
        // Step 5: particle loss at boundaries
        // boundaries(p_e, n_active_e, lpos_e);
        if(i % K_SUB == 0) boundaries_i(p_i, n_active_i, lpos_i, n_out_i);
        boundaries_e(p_e, n_active_e, lpos_e, n_out_i);

        // Step 6: Monte-Carlo collisions
        collisions_e(p_e, n_active_e, lpos_e, p_i, n_active_i, lpos_i, M_I, N_NEUTRAL, p_null_e, nu_prime_e);
        if(i % K_SUB == 0) collisions_i(p_i, n_active_i, M_I, N_NEUTRAL, p_null_i, nu_prime_i);

        print_info(i, p_e, n_active_e, p_i, n_active_i, 100);
        // average_field(phi_av, phi, i);
        // average_field(efield_av_x, efield_x, i);
        // average_field(dens_e_av, wmesh_e, i);
        // average_field(dens_i_av, wmesh_i, i);

	}
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Total execution duration: " << (double) duration.count() / (1.0e6 * 60) << " min" << endl;

    // ----------------------------- Saving outputs ---------------------------

    fmatrix dens_e_corrected = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_e / vmesh;
    fmatrix dens_i_corrected = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_i / vmesh;
    fmatrix phi_corrected  = phi * (M_EL * pow(DX, 2))/(Q * pow(DT, 2));
    
//    print_fmatrix(phi_corrected);
//    print_fmatrix(solver.node_type);
//
    // efield_x = efield_x * (M_EL * DX/ (Q * pow(DT, 2)));
    save_to_csv(dens_e_corrected, "dens_e.csv");
    save_to_csv(dens_i_corrected, "dens_i.csv");
	save_to_csv(phi_corrected, "phi.csv");
    
    // ----------------------------- Finalizing -------------------------------
    delete_cross_sections_arrays();
	return 0;
}
