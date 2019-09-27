
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

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
    double v_cap            = -100.0;
    
    // Average field variablesd
    verbose_log("Initializing diagnostics variables");
    fmatrix phi_av          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_av_x     = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_av_y     = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix dens_e_av       = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix dens_i_av       = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix kefield_i       = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix kefield_e       = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix v_cap_diag      = fmatrix::zeros(N_STEPS);
    fmatrix n_active_e_diag = fmatrix::zeros(N_STEPS);
    fmatrix n_active_i_diag = fmatrix::zeros(N_STEPS);
    
	// Particle 1 - Electrons
    verbose_log("Initializing electrons variables");
	int n_active_e          = 0;
	fmatrix p_e             = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_e          = imatrix::zeros(N_MAX_PARTICLES, 2);
	fmatrix wmesh_e         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_x_at_p_e = fmatrix::zeros(N_MAX_PARTICLES);
    fmatrix efield_y_at_p_e = fmatrix::zeros(N_MAX_PARTICLES);
    double n_inj_e          = N_INJ_EL;
    int n_out_e             = 0;
	
	// Particle 2 - Ions
    verbose_log("Initializing ions variables");
	int n_active_i          = 0;
	fmatrix p_i             = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_i          = imatrix::zeros(N_MAX_PARTICLES, 2);
	fmatrix wmesh_i         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
	fmatrix efield_x_at_p_i = fmatrix::zeros(N_MAX_PARTICLES);
    fmatrix efield_y_at_p_i = fmatrix::zeros(N_MAX_PARTICLES);
    int n_out_i             = 0;
    
    // Particle 3 - Neutrals
    verbose_log("Initializing neutral variables");
    int n_active_n          = 0;
    fmatrix p_n             = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_n          = imatrix::zeros(N_MAX_PARTICLES, 2);
    fmatrix wmesh_n         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix dens_n          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix wmesh_n_av      = fmatrix::zeros(N_MESH_X, N_MESH_Y);

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
      
    // ----------------------------- DSMC loop --------------------------------
    if (EXPM_NEUTRAL == "dsmc"){
        verbose_log(" ---- Starting DSMC loop ---- ");
        for (int i = 0; i < N_STEPS_DSMC; i++){
            add_flux_particles(p_n, n_active_n, T_NEUTRAL, 0, M_I, N_INJ_N, K_SUB_DSMC);

            if(i > N_STEPS_DSMC - N_AVERAGE_DSMC){
                weight(p_n, n_active_n, wmesh_n, mesh_x, mesh_y, lpos_n);
                average_field(wmesh_n_av, wmesh_n, i - (N_STEPS_DSMC - N_AVERAGE_DSMC));
            }
            
            move_n(p_n, n_active_n, K_SUB_DSMC);
            boundaries_n(p_n, n_active_n, lpos_n);
            print_dsmc_info(i, n_active_n, 1000, N_STEPS_DSMC);    
        }
        wmesh_n = (N_FACTOR_DSMC / N_FACTOR) * wmesh_n_av;
        dens_n = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_n / vmesh;
        save_to_csv(dens_n, "dens_n.csv");
    }
    else if (EXPM_NEUTRAL == "constant") {
        verbose_log(" ---- Setting constant neutral density ---- ");
        dens_n.set_all(N_NEUTRAL);
    } 
    else if (EXPM_NEUTRAL == "load") {
        verbose_log(" ---- Loading neutral density field ---- ");
        dens_n = load_csv("output/dens_n.csv",',', N_MESH_Y);
    }
      
    // -----------------------------------------------------------------------
    
    // Initial MCC parameters
    double nu_prime_e       = find_nu_prime_e(dens_n, vmesh);
    double p_null_e         = p_null(nu_prime_e, DT);
    double nu_prime_i       = find_nu_prime_i(dens_n, vmesh);
    double p_null_i         = p_null(nu_prime_i, DT * K_SUB);
    
     
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
        double alpha = 2e-4;
        solver.solve(phi, voltages, wmesh_i, wmesh_e);
        phi_av = (1 - alpha) * phi_av + alpha * phi;

        // Step 2.1: calculation of electric field
        calculate_efield(efield_x, efield_y, phi, wmesh_i, wmesh_e, mesh_x, mesh_y, vmesh);

        // Step 2.2: field weighting
        if(i % K_SUB == 0) electric_field_at_particles(efield_x_at_p_i, efield_y_at_p_i, efield_x, efield_y, p_i, n_active_i, mesh_x, mesh_y, lpos_i);
        electric_field_at_particles(efield_x_at_p_e, efield_y_at_p_e, efield_x, efield_y, p_e, n_active_e, mesh_x, mesh_y, lpos_e);

        // Step 3: integration of equations of motion
        if(i % K_SUB == 0) move_i(p_i, n_active_i, efield_x_at_p_i, efield_y_at_p_i);
        move_e(p_e, n_active_e, efield_x_at_p_e, efield_y_at_p_e);
        
        // Step 4: particle loss at boundaries
        // boundaries(p_e, n_active_e, lpos_e);
        if(i % K_SUB == 0) boundaries_i(p_i, n_active_i, lpos_i, n_out_i);
        // boundaries_e(p_e, n_active_e, lpos_e, n_out_i);
        boundaries_e_cap(p_e, n_active_e, lpos_e, n_out_e, v_cap, phi_av, mesh_x, mesh_y);
        v_cap = cap_voltage(v_cap, n_out_e, n_out_i);
        
        // Step 5: particles injection
        double v_drift_e = 0;
        if(i % K_SUB == 0) add_flux_particles(p_i, n_active_i, T_I, V_DRIFT_I, M_I, N_INJ_I, K_SUB);
        if(INJ_MODEL == "constant")       n_inj_e = N_INJ_EL;
        else if(INJ_MODEL == "balanced")  n_inj_e = balanced_injection(n_inj_e, 0.01, wmesh_i, wmesh_e, 0, 0, 0, N_THRUSTER - 1);
        else if(INJ_MODEL == "pulsed"){
            n_inj_e = pulsed_injection(K_INJ_EL, V_SB, V_RF, T_EL, OMEGA_I, i);
            // double voltage_drift_e = (0.95 - 1) * (V_SB + V_RF * sin(2 * PI * FREQ * DT * i)) + 50;
            // voltage_drift_e = voltage_drift_e > 0 ? voltage_drift_e : 0;
            double voltage_drift_e = 20;
            v_drift_e = sqrt(2 * Q * voltage_drift_e / M_EL);
        }
        add_flux_particles(p_e, n_active_e, T_EL, v_drift_e, M_EL, n_inj_e);
        

        // Step 6: Monte-Carlo collisions
        if(i % K_SUB == 0) collisions_i(p_i, n_active_i, lpos_i, mesh_x, mesh_y, dens_n, M_I, p_null_i, nu_prime_i);
        collisions_e(p_e, n_active_e, lpos_e, p_i, n_active_i, lpos_i, mesh_x, mesh_y, dens_n, M_I, p_null_e, nu_prime_e);
        
        //  ----------------------------- Diagnostics -------------------------

        print_info(i, p_e, n_active_e, p_i, n_active_i, v_cap, 100);
        v_cap_diag.val[i] = v_cap;
        n_active_e_diag.val[i] = n_active_e;
        n_active_i_diag.val[i] = n_active_i;


        if(i % 10000 == 0) save_state(p_e, n_active_e, p_i, n_active_i, phi, wmesh_e, wmesh_i, vmesh, "_state");
        // if((i % 100 == 0) && i >= 595000) save_state(p_e, n_active_e, p_i, n_active_i, phi, wmesh_e, wmesh_i, vmesh, "_state_" + to_string(i));
        

	}
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Total execution duration: " << (double) duration.count() / (1.0e6 * 60) << " min" << endl;

    // ----------------------------- Saving outputs ---------------------------
    
    save_state(p_e, n_active_e, p_i, n_active_i, phi, wmesh_e, wmesh_i, vmesh, "_state");
    save_to_csv(v_cap_diag, "v_cap.csv");
    save_to_csv(n_active_e_diag, "n_active_e.csv");
    save_to_csv(n_active_i_diag, "n_active_i.csv");

   
    // ----------------------------- Finalizing -------------------------------
    delete_cross_sections_arrays();
	return 0;
}
