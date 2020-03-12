
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>

#include "cross-sections.h"
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
#include "state-info.h"
#include "arg-parser.h"
#include "configuration.h"

#define CONFIG_PATH "input/config/config.ini"

using namespace std;
using namespace std::chrono;

// ----------------------------- Main function --------------------------------
int main(int argc, char* argv[])
{
    // ------------------- Loading configuration ------------------------------

    argparser arg(argc, argv);

    string job_suffix = arg.get("name", "");
    string config_path = arg.get("config", CONFIG_PATH);

    configuration config("input/config/config.yaml");

    
    load_config_file(config_path);   
    load_cross_sections(config);
    

    // ------------------- Variable initialization ----------------------------
    
	
    verbose_log("Initializing general variables");

    // Loading variables from config

    const int n_steps = config.i("time/n_steps");
    const double dt = config.f("time/dt");
    const int k_sub = config.i("time/k_sub");
    
    const int n_mesh_x = config.i("geometry/n_mesh_x");
    const int n_mesh_y = config.i("geometry/n_mesh_y");
    const int n_thruster = config.i("geometry/n_thruster");

    const double freq = config.f("thruster/freq");
    const double v_sb = config.f("thruster/v_sb");
    const double v_rf = config.f("thruster/v_rf");
    
    const int n_max_particles = config.i("particles/n_max_particles");
    
    const double sqr_dt = config.f("electrons/sqr_duty_cycle");
    const double m_el = config.f("electrons/m_el");
    const double t_el = config.f("electrons/t_el");
    const double v_drift_el = config.f("electrons/v_drift_el");

    const double v_drift_i = config.f("ions/v_drift_i");
    const double t_i = config.f("ions/t_i");

    const double m_i = config.f("ugas/m_i");

    const double k_inj_el = config.f("p/k_inj_el");
    const double omega_i = config.f("p/omega_i");
    const double volt_0_norm = config.f("p/volt_0_norm");
    const double volt_1_norm = config.f("p/volt_1_norm");

 
    // General field variables

	fmatrix phi             = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix phi_laplace     = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix phi_poisson     = fmatrix::zeros(n_mesh_x, n_mesh_y);
	fmatrix efield_x        = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix efield_y        = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix mesh_x          = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix mesh_y          = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix vmesh           = fmatrix::zeros(n_mesh_x, n_mesh_y);
    imatrix electrode_mask  = imatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix voltages        = fmatrix::zeros(4);
    double sigma_laplace;
    state_info state;
    
    // Doiagnostics variables
    verbose_log("Initializing diagnostics variables");
    fmatrix phi_av          = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix wmesh_e_av      = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix wmesh_i_av      = fmatrix::zeros(n_mesh_x, n_mesh_y);
    
    unordered_map<string, fmatrix> series;
    series["time"]          = fmatrix::zeros(n_steps);
    series["v_cap"]         = fmatrix::zeros(n_steps);
    series["n_active_e"]    = fmatrix::zeros(n_steps);
    series["n_active_i"]    = fmatrix::zeros(n_steps);
    series["I_out_ob_e"]        = fmatrix::zeros(n_steps);
    series["I_out_ob_i"]        = fmatrix::zeros(n_steps);
    series["I_out_thr_e"]       = fmatrix::zeros(n_steps);
    series["I_out_thr_i"]       = fmatrix::zeros(n_steps);
    series["I_in_thr_e"]        = fmatrix::zeros(n_steps);
    series["I_in_thr_i"]        = fmatrix::zeros(n_steps);
    int n_points_series     = 0;
    
	// Particle 1 - Electrons
    verbose_log("Initializing electrons variables");
	fmatrix p_e             = fmatrix::zeros(n_max_particles, 6);
    imatrix lpos_e          = imatrix::zeros(n_max_particles, 2);
	fmatrix wmesh_e         = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix efield_x_at_p_e = fmatrix::zeros(n_max_particles);
    fmatrix efield_y_at_p_e = fmatrix::zeros(n_max_particles);
    double n_inj_el          = config.f("p/n_inj_el");
	
	// Particle 2 - Ions
    verbose_log("Initializing ions variables");
	fmatrix p_i             = fmatrix::zeros(n_max_particles, 6);
    imatrix lpos_i          = imatrix::zeros(n_max_particles, 2);
	fmatrix wmesh_i         = fmatrix::zeros(n_mesh_x, n_mesh_y);
	fmatrix efield_x_at_p_i = fmatrix::zeros(n_max_particles);
    fmatrix efield_y_at_p_i = fmatrix::zeros(n_max_particles);
    double n_inj_i          = config.f("p/n_inj_i");
    
    // Particle 3 - Neutrals
    verbose_log("Initializing neutral variables");
    fmatrix dens_n          = fmatrix::zeros(n_mesh_x, n_mesh_y);

    // ----------------------------- Mesh -------------------------------------
    verbose_log("Initializing mesh");
    init_mesh(mesh_x,  mesh_y, config);
    init_volume_mesh(vmesh, mesh_x, mesh_y);
    	
    // ----------------------------- Solver -----------------------------------
    verbose_log("Initializing solver");

    rsolver solver(mesh_x, mesh_y, vmesh, config);
    setup_rsolver(solver, mesh_x, mesh_y, vmesh, electrode_mask);

    voltages = {1, 0, 1, 1};
    solver.solve(phi_laplace, voltages, wmesh_i, wmesh_e);

    voltages = {volt_1_norm, volt_0_norm, volt_1_norm, volt_1_norm};
    sigma_laplace = sigma_from_phi(phi_laplace, mesh_x, mesh_y, wmesh_e, wmesh_e, vmesh, electrode_mask, config);

    // ----------------------------- DSMC -------------------------------------
    
    if (EXPM_NEUTRAL == "dsmc"){
        verbose_log(" ---- Starting DSMC loop ---- ");
        run_dsmc(mesh_x, mesh_y, vmesh, dens_n, config, "dens_n" + job_suffix + ".exdir");
    }
    else if (EXPM_NEUTRAL == "constant") {
        verbose_log(" ---- Setting constant neutral density ---- ");
        dens_n.set_all(N_NEUTRAL);
    } 
    else if (EXPM_NEUTRAL == "load") {
        verbose_log(" ---- Loading neutral density field ---- ");
        load_fmatrix(dens_n, config.s("project/input_path") + "dens_n.exdir", "dens_n");
    }

    // ----------------------------- MCC --------------------------------------
    
    double nu_prime_e       = find_nu_prime_e(dens_n, vmesh, config);
    double p_null_e         = p_null(nu_prime_e, dt);
    double nu_prime_i       = find_nu_prime_i(dens_n, vmesh, config);
    double p_null_i         = p_null(nu_prime_i, dt * (double) k_sub);
         
    // ----------------------------- Loading initial state -------------------

    if(INITIAL_STATE == "load")  load_state(p_e, p_i, state);

	// ----------------------------- Main loop --------------------------------

    // Printing initial information
    print_initial_info(p_null_e, p_null_i);

    verbose_log(" ---- Starting main loop ---- ");
    auto start = high_resolution_clock::now();
	for (state.step = state.step_offset; state.step < n_steps + state.step_offset; state.step++)
	{ 
		// Step 1.0: particle weighting
        if(state.step % k_sub == 0) weight(p_i, state.n_active_i, wmesh_i, mesh_x, mesh_y, lpos_i);
        weight(p_e, state.n_active_e, wmesh_e, mesh_x, mesh_y, lpos_e);

        // Step 2.0 integration of Poisson's equation
        solver.solve(phi_poisson, voltages, wmesh_i, wmesh_e);
        
        state.sigma_0 = state.sigma_1;
        state.phi_zero = calculate_phi_zero(state.sigma_1, state.n_out_ob_i -  state.n_out_ob_e, state.q_cap, sigma_laplace, phi_poisson, mesh_x, mesh_y, wmesh_e, wmesh_i, vmesh, electrode_mask, config);
        
        state.sigma_1 = calculate_sigma(state.sigma_0, state.phi_zero, state.n_out_ob_i -  state.n_out_ob_e, state.q_cap, config);
        state.q_cap = calculate_cap_charge(state.sigma_1, state.sigma_0, state.q_cap, state.n_out_ob_i -  state.n_out_ob_e, config);

        phi = phi_poisson + (state.phi_zero * phi_laplace);

        // Step 2.1: calculation of electric field
        calculate_efield(efield_x, efield_y, phi, wmesh_i, wmesh_e, mesh_x, mesh_y, vmesh, electrode_mask);

        // Step 2.2: field weighting
        if(state.step % k_sub == 0) electric_field_at_particles(efield_x_at_p_i, efield_y_at_p_i, efield_x, efield_y, p_i, state.n_active_i, mesh_x, mesh_y, lpos_i);
        electric_field_at_particles(efield_x_at_p_e, efield_y_at_p_e, efield_x, efield_y, p_e, state.n_active_e, mesh_x, mesh_y, lpos_e);

        // Step 3: integration of equations of motion
        if(state.step % k_sub == 0) move_i(p_i, state.n_active_i, efield_x_at_p_i, efield_y_at_p_i);
        move_e(p_e, state.n_active_e, efield_x_at_p_e, efield_y_at_p_e);

        // // Step 4: particle loss at boundaries
        if(state.step % k_sub == 0) boundaries_ob_count(p_i, state.n_active_i, lpos_i, state.n_out_ob_i, state.n_out_thr_i);
        boundaries_ob_count(p_e, state.n_active_e, lpos_e, state.n_out_ob_e, state.n_out_thr_e);

        // Step 5: particles injection
        if(state.step % k_sub == 0) add_flux_particles(p_i, state.n_active_i, t_i, v_drift_i, m_i, n_inj_i, k_sub);
        if(INJ_MODEL == "constant")       n_inj_el = n_inj_el;
        else if(INJ_MODEL == "balanced")  n_inj_el = balanced_injection(n_inj_el, 0.01, wmesh_i, wmesh_e, 0, 0, 0, n_thruster - 1);
        else if(INJ_MODEL == "pulsed")    n_inj_el = pulsed_injection(k_inj_el, v_sb, v_rf, t_el, omega_i, state.step);
        else if(INJ_MODEL == "square")    n_inj_el = square_injection(4.0, freq, dt, sqr_dt, state.step);
        add_flux_particles(p_e, state.n_active_e, t_el, v_drift_el, m_el, n_inj_el, 1);

        // Step 6: Monte-Carlo collisions
        if(MCC_COLL == "on"){
            if(state.step % K_SUB == 0) collisions_i(p_i, state.n_active_i, lpos_i, mesh_x, mesh_y, dens_n, p_null_i, nu_prime_i, config);
            collisions_e(p_e, state.n_active_e, lpos_e, p_i, state.n_active_i, lpos_i, mesh_x, mesh_y, dens_n, p_null_e, nu_prime_e, config);
        }
        
        //  ----------------------------- Diagnostics -------------------------

        print_info(state, 100);

        if(state.step % 10000 == 0) 
        {
            save_state(p_e, p_i, state, job_suffix);
            save_fields_snapshot(phi, wmesh_e, wmesh_i, vmesh, state, job_suffix);
        }

        // if(state.step % 1 == 0){
        //     series["time"].val[n_points_series] = state.step * DT;
        //     series["v_cap"].val[n_points_series] = state.phi_zero / K_PHI;
        //     series["n_active_i"].val[n_points_series] = state.n_active_i;
        //     series["n_active_e"].val[n_points_series] = state.n_active_e;
            
        //     series["I_out_ob_e"].val[n_points_series] = state.n_out_ob_e * N_FACTOR * Q / DT;
        //     series["I_out_ob_i"].val[n_points_series] = state.n_out_ob_i * N_FACTOR * Q / DT;
        //     series["I_out_thr_e"].val[n_points_series] = state.n_out_thr_e * N_FACTOR * Q / DT;
        //     series["I_out_thr_i"].val[n_points_series] = state.n_out_thr_i * N_FACTOR * Q / DT;
        //     series["I_in_thr_e"].val[n_points_series] = n_inj_e * N_FACTOR * Q / DT;
        //     series["I_in_thr_i"].val[n_points_series] = n_inj_i * N_FACTOR * Q / DT;
        //     n_points_series += 1;
        // }

        // if(state.step % 50000 == 0) {
        //     save_series(series, n_points_series, state, "");
        // }
        
        //  if(state.step % 1 == 0) {
        //     fmatrix dens_i = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_i / vmesh;
        //     save_field_series(dens_i, state, 1, "_dens_i");
        //     fmatrix dens_e = (4 / pow(DX, 2)) *  N_FACTOR * wmesh_e / vmesh;
        //     save_field_series(dens_e, state, 1, "_dens_e");
            
        //      fmatrix kefield = fmatrix::zeros(N_MESH_X, N_MESH_Y);
        //      energy_field(kefield, p_i, state.n_active_i, mesh_x, mesh_y, vmesh, lpos_i, M_I);
        //      save_field_series(kefield, state, 1, "_ke_i");
            
        //      energy_field(kefield, p_e, state.n_active_e, mesh_x, mesh_y, vmesh, lpos_e, M_EL);
        //      save_field_series(kefield, state, 1, "_ke_e");
            
        //      fmatrix phi_corrected = phi / K_PHI;
        //      save_field_series(phi_corrected, state, 1, "_phi");
        //  }

        // if(state.step - state.step_offset > 0.9 * (double) n_steps){
        //     int i_av = average_field_over_period(phi_av, phi, RF_PERIOD_I, n_steps, state.step - state.step_offset);
        //     average_field(wmesh_e_av, wmesh_e, state.step - state.step_offset);
        //     average_field(wmesh_i_av, wmesh_i, state.step - state.step_offset);
        //     if(i_av == RF_PERIOD_I){
        //         save_fields_snapshot(phi_av, wmesh_e_av, wmesh_i_av, vmesh, state, "_av");
        //     }
        // }
        
	}

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Total execution duration: " << (double) duration.count() / (1.0e6 * 60) << " min" << endl;

    // ----------------------------- Saving outputs ---------------------------
    save_state(p_e, p_i, state, job_suffix);
    save_fields_snapshot(phi, wmesh_e, wmesh_i, vmesh, state, job_suffix);
    save_series(series, n_points_series, state, job_suffix);

    // ----------------------------- Finalizing -------------------------------
    delete_cross_sections_arrays();
	return 0;
}
