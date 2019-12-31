
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

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
#include "input-output.h"
#include "num-tools.h"
#include "dsmc.h"
#include "state_info.h"

#define CONFIG_PATH "input/config/config.ini"

using namespace std;
using namespace std::chrono;

// ----------------------------- Main function --------------------------------
int main(int argc, char* argv[])
{

    // ------------------- Loading configuration ------------------------------

    load_config_file(argv[1] == NULL ? CONFIG_PATH : argv[1]);
    load_cross_sections();

    // ------------------- Variable initialization ----------------------------
    
	// General field variables
    verbose_log("Initializing general variables");
	fmatrix phi             = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix phi_laplace     = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix phi_poisson     = fmatrix::zeros(N_MESH_X, N_MESH_Y);
	fmatrix efield_x        = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_y        = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix mesh_x          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix mesh_y          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix vmesh           = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    imatrix electrode_mask  = imatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix voltages        = fmatrix::zeros(4);
    double sigma_laplace;
    state_info state;
    
    // Doiagnostics variables
    verbose_log("Initializing diagnostics variables");
    fmatrix phi_av          = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix wmesh_e_av      = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix wmesh_i_av      = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    
    map<string, fmatrix> series;
    series["time"]          = fmatrix::zeros(N_STEPS);
    series["v_cap"]         = fmatrix::zeros(N_STEPS);
    series["n_active_e"]    = fmatrix::zeros(N_STEPS);
    series["n_active_i"]    = fmatrix::zeros(N_STEPS);
    int n_points_series     = 0;
    
	// Particle 1 - Electrons
    verbose_log("Initializing electrons variables");
	fmatrix p_e             = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_e          = imatrix::zeros(N_MAX_PARTICLES, 2);
	fmatrix wmesh_e         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix efield_x_at_p_e = fmatrix::zeros(N_MAX_PARTICLES);
    fmatrix efield_y_at_p_e = fmatrix::zeros(N_MAX_PARTICLES);
    double n_inj_e          = N_INJ_EL;
	
	// Particle 2 - Ions
    verbose_log("Initializing ions variables");
	fmatrix p_i             = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_i          = imatrix::zeros(N_MAX_PARTICLES, 2);
	fmatrix wmesh_i         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
	fmatrix efield_x_at_p_i = fmatrix::zeros(N_MAX_PARTICLES);
    fmatrix efield_y_at_p_i = fmatrix::zeros(N_MAX_PARTICLES);
    double n_inj_i          = N_INJ_I;
    
    // Particle 3 - Neutrals
    verbose_log("Initializing neutral variables");
    fmatrix dens_n          = fmatrix::zeros(N_MESH_X, N_MESH_Y);

    // ----------------------------- Mesh -------------------------------------
    verbose_log("Initializing mesh");
    init_mesh(mesh_x,  mesh_y, A_X, A_Y, N_MESH_X, N_MESH_Y);
    init_volume_mesh(vmesh, mesh_x, mesh_y);


    // ----------------------------- Solver -----------------------------------
    verbose_log("Initializing solver");

    rsolver solver(mesh_x, mesh_y, vmesh);
    setup_rsolver(solver, mesh_x, mesh_y, vmesh, electrode_mask);

    voltages = {1, 0, 1, 1};
    solver.solve(phi_laplace, voltages, wmesh_i, wmesh_e);
    voltages = {VOLT_1_NORM, VOLT_0_NORM, VOLT_1_NORM, VOLT_1_NORM};
    sigma_laplace = sigma_from_phi(phi_laplace, mesh_x, mesh_y, wmesh_e, wmesh_e, vmesh, electrode_mask);

    // ----------------------------- DSMC -------------------------------------
    
    if (EXPM_NEUTRAL == "dsmc"){
        verbose_log(" ---- Starting DSMC loop ---- ");
        run_dsmc(mesh_x, mesh_y, vmesh, dens_n);
    }
    else if (EXPM_NEUTRAL == "constant") {
        verbose_log(" ---- Setting constant neutral density ---- ");
        dens_n.set_all(N_NEUTRAL);
    } 
    else if (EXPM_NEUTRAL == "load") {
        verbose_log(" ---- Loading neutral density field ---- ");
        dens_n = load_csv("input/dens_n.csv",',', N_MESH_Y);
    }
      
    // ----------------------------- MCC --------------------------------------
    
    double nu_prime_e       = find_nu_prime_e(dens_n, vmesh);
    double p_null_e         = p_null(nu_prime_e, DT);
    double nu_prime_i       = find_nu_prime_i(dens_n, vmesh);
    double p_null_i         = p_null(nu_prime_i, DT * K_SUB);
         
    // ----------------------------- Loading initial state -------------------

    if(INITIAL_STATE == "load")  load_state(p_e, p_i, state);

	// ----------------------------- Main loop --------------------------------

    // Printing initial information
    print_initial_info(p_null_e, p_null_i);

    verbose_log(" ---- Starting main loop ---- ");
    auto start = high_resolution_clock::now();
	for (state.step = state.step_offset; state.step < N_STEPS + state.step_offset; state.step++)
	{ 
		// Step 1.0: particle weighting
        if(state.step % K_SUB == 0) weight(p_i, state.n_active_i, wmesh_i, mesh_x, mesh_y, lpos_i);
        weight(p_e, state.n_active_e, wmesh_e, mesh_x, mesh_y, lpos_e);

        // Step 2.0 integration of Poisson's equation
        solver.solve(phi_poisson, voltages, wmesh_i, wmesh_e);
        
        state.sigma_0 = state.sigma_1;
        state.phi_zero = calculate_phi_zero(state.sigma_1, state.n_out_i -  state.n_out_e, state.q_cap, sigma_laplace, phi_poisson, mesh_x, mesh_y, wmesh_e, wmesh_i, vmesh, electrode_mask);
        
        state.sigma_1 = calculate_sigma(state.sigma_0, state.phi_zero, state.n_out_i -  state.n_out_e, state.q_cap);
        state.q_cap = calculate_cap_charge(state.sigma_1, state.sigma_0, state.q_cap, state.n_out_i -  state.n_out_e);

        phi = phi_poisson + (state.phi_zero * phi_laplace);

        // Step 2.1: calculation of electric field
        calculate_efield(efield_x, efield_y, phi, wmesh_i, wmesh_e, mesh_x, mesh_y, vmesh, electrode_mask);

        // Step 2.2: field weighting
        if(state.step % K_SUB == 0) electric_field_at_particles(efield_x_at_p_i, efield_y_at_p_i, efield_x, efield_y, p_i, state.n_active_i, mesh_x, mesh_y, lpos_i);
        electric_field_at_particles(efield_x_at_p_e, efield_y_at_p_e, efield_x, efield_y, p_e, state.n_active_e, mesh_x, mesh_y, lpos_e);

        // Step 3: integration of equations of motion
        if(state.step % K_SUB == 0) move_i(p_i, state.n_active_i, efield_x_at_p_i, efield_y_at_p_i);
        move_e(p_e, state.n_active_e, efield_x_at_p_e, efield_y_at_p_e);

        // // Step 4: particle loss at boundaries
        if(state.step % K_SUB == 0) boundaries_ob_count(p_i, state.n_active_i, lpos_i, state.n_out_i);
        boundaries_ob_count(p_e, state.n_active_e, lpos_e, state.n_out_e);

        // Step 5: particles injection
        if(state.step % K_SUB == 0) add_flux_particles(p_i, state.n_active_i, T_I, V_DRIFT_I, M_I, n_inj_i, K_SUB);
        if(INJ_MODEL == "constant")       n_inj_e = N_INJ_EL;
        else if(INJ_MODEL == "balanced")  n_inj_e = balanced_injection(n_inj_e, 0.01, wmesh_i, wmesh_e, 0, 0, 0, N_THRUSTER - 1);
        else if(INJ_MODEL == "pulsed")    n_inj_e = pulsed_injection(K_INJ_EL, V_SB, V_RF, T_EL, OMEGA_I, state.step);
        else if(INJ_MODEL == "square")    n_inj_e = square_injection(4.0, FREQ, DT, SQR_DUTY_CYCLE, state.step);
        add_flux_particles(p_e, state.n_active_e, T_EL, V_DRIFT_EL, M_EL, n_inj_e, 1);

        // Step 6: Monte-Carlo collisions
        if(MCC_COLL == "on"){
            if(state.step % K_SUB == 0) collisions_i(p_i, state.n_active_i, lpos_i, mesh_x, mesh_y, dens_n, M_I, p_null_i, nu_prime_i);
            collisions_e(p_e, state.n_active_e, lpos_e, p_i, state.n_active_i, lpos_i, mesh_x, mesh_y, dens_n, M_I, p_null_e, nu_prime_e);
        }
        
        //  ----------------------------- Diagnostics -------------------------

        print_info(state, 100);

        if(state.step % 10000 == 0) 
        {
            save_state(p_e, p_i, state);
            save_fields_snapshot(phi, wmesh_e, wmesh_i, vmesh, state, "");
        }

        if(state.step % 5 == 0){
            series["time"].val[n_points_series] = state.step * DT;
            series["v_cap"].val[n_points_series] = state.phi_zero / K_PHI;
            series["n_active_i"].val[n_points_series] = state.n_active_i;
            series["n_active_e"].val[n_points_series] = state.n_active_e;
            n_points_series += 1;
        }

        if(state.step % 20000 == 0) {
            save_series(series, n_points_series, "");
        }

        if(state.step > (double) 0.9 * N_STEPS){
            int i_av = average_field_over_period(phi_av, phi, RF_PERIOD_I, N_STEPS, state.step - state.step_offset);
            if(i_av == RF_PERIOD_I){
                save_fields_snapshot(phi_av, wmesh_e_av, wmesh_i_av, vmesh, state, "_av");
            }
        }
        
	}

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Total execution duration: " << (double) duration.count() / (1.0e6 * 60) << " min" << endl;

    // ----------------------------- Saving outputs ---------------------------
    save_state(p_e, p_i, state);
    save_fields_snapshot(phi, wmesh_e, wmesh_i, vmesh, state, "");
    save_series(series, n_points_series, "");

    // ----------------------------- Finalizing -------------------------------
    delete_cross_sections_arrays();
	return 0;
}
