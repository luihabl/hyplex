
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>

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
#include "clock.h"

#define CONFIG_PATH "input/config/config.yaml"

using namespace std;
using namespace std::chrono;

// ----------------------------- Main function --------------------------------
int main(int argc, char* argv[])
{

    // ------------------- Loading configuration ------------------------------

    argparser arg(argc, argv);

    string job_suffix = arg.get("name", "");
    string config_path = arg.get("config", CONFIG_PATH);

    configuration config(config_path);

    // ------------------- Variable initialization ----------------------------
    
	
    verbose_log("Initializing general variables");

    // Loading variables from config
    
    const int n_steps           = config.i("time/n_steps");
    const int rf_period_i       = config.i("p/rf_period_i");
    const int k_sub             = config.i("time/k_sub");
    const int n_mesh_x          = config.i("geometry/n_mesh_x");
    const int n_mesh_y          = config.i("geometry/n_mesh_y");
    const int n_thruster        = config.i("geometry/n_thruster");
    const int n_max_particles   = config.i("particles/n_max_particles");
    const double dt             = config.f("time/dt");
    const double q              = config.f("physical/q");
    const double m_el           = config.f("electrons/m_el");
    const double t_el           = config.f("electrons/t_el");
    const double v_drift_el     = config.f("electrons/v_drift_el");
    const double v_drift_i      = config.f("ions/v_drift_i");
    const double t_i            = config.f("ions/t_i");
    const double m_i            = config.f("ugas/m_i");
    const double volt_0_norm    = config.f("p/volt_0_norm");
    const double volt_1_norm    = config.f("p/volt_1_norm");
    const double k_phi          = config.f("p/k_phi");
    const double n_factor       = config.f("particles/n_factor");
    const bool mcc_coll         = config.b("neutrals/mcc_coll");
    const string inj_model      = config.s("electrons/inj_model");
 
    // General field variables

	fmatrix phi             = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix phi_laplace     = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix phi_poisson     = fmatrix::zeros(n_mesh_x, n_mesh_y);
	fmatrix efield_x        = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix efield_y        = fmatrix::zeros(n_mesh_x, n_mesh_y);
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
    series["time"]              = fmatrix::zeros(n_steps);
    series["v_cap"]             = fmatrix::zeros(n_steps);
    series["n_active_e"]        = fmatrix::zeros(n_steps);
    series["n_active_i"]        = fmatrix::zeros(n_steps);
    series["I_out_ob_e"]        = fmatrix::zeros(n_steps);
    series["I_out_ob_i"]        = fmatrix::zeros(n_steps);
    series["I_out_thr_e"]       = fmatrix::zeros(n_steps);
    series["I_out_thr_i"]       = fmatrix::zeros(n_steps);
    series["I_in_thr_e"]        = fmatrix::zeros(n_steps);
    series["I_in_thr_i"]        = fmatrix::zeros(n_steps);
    int n_points_series         = 0;
    
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

    // ----------------------------- Operation objects -------------------------
    verbose_log("Initializing mesh");

    mesh_set mesh(config);
    mesh.init_mesh();
    
    field_operations fields(config);
    pic_operations pic(config);
    particle_operations pops(config, pic);

    // ----------------------------- Solver -----------------------------------
    verbose_log("Initializing solver");

    rsolver solver(mesh, config);
    setup_rsolver(solver, mesh, electrode_mask);

    voltages = {1, 0, 1, 1};
    solver.solve(phi_laplace, voltages, wmesh_i, wmesh_e);

    voltages = {volt_1_norm, volt_0_norm, volt_1_norm, volt_1_norm};
    sigma_laplace = fields.sigma_from_phi(phi_laplace, mesh, wmesh_e, wmesh_i, electrode_mask);

    // ----------------------------- DSMC -------------------------------------
    
    string expm_neutral = config.s("neutrals/expm_neutral");

    if (expm_neutral == "dsmc"){
        verbose_log(" ---- Starting DSMC loop ---- ");
        run_dsmc(mesh, dens_n, config, "dens_n" + job_suffix + ".exdir");
    }
    else if (expm_neutral == "constant") {
        verbose_log(" ---- Setting constant neutral density ---- ");
        dens_n.set_all(config.f("neutrals/n_neutral"));
    } 
    else if (expm_neutral == "load") {
        verbose_log(" ---- Loading neutral density field ---- ");
        load_fmatrix(dens_n, config.s("project/input_path") + "dens_n.exdir", "dens_n");
    }

    // ----------------------------- MCC --------------------------------------
    
    mcc coll(config, pops, pic);
    coll.initialize_mcc(dens_n, mesh.v);
    
    // ----------------------------- Loading initial state -------------------

    if(config.s("simulation/initial_state") == "load")  load_state(p_e, p_i, state, config);

	// ----------------------------- Main loop --------------------------------

    // Printing initial information
    print_initial_info(coll.p_null_e, coll.p_null_i, config);

    verbose_log(" ---- Starting main loop ---- ");
    auto start = now();
	for (state.step = state.step_offset; state.step < n_steps + state.step_offset; state.step++)
	{
		// Step 1.0: particle weighting
        if(state.step % k_sub == 0) pic.weight(p_i, state.n_active_i, wmesh_i, mesh, lpos_i);
        pic.weight(p_e, state.n_active_e, wmesh_e, mesh, lpos_e);

        // Step 2.0 integration of Poisson's equation
        solver.solve(phi_poisson, voltages, wmesh_i, wmesh_e);
        
        state.sigma_0 = state.sigma_1;
        state.phi_zero = fields.calculate_phi_zero(state.sigma_1, state.n_out_ob_i -  state.n_out_ob_e, state.q_cap, sigma_laplace, phi_poisson, mesh, wmesh_e, wmesh_i, electrode_mask);
        
        state.sigma_1 = fields.calculate_sigma(state.sigma_0, state.phi_zero, state.n_out_ob_i -  state.n_out_ob_e, state.q_cap);
        state.q_cap = fields.calculate_cap_charge(state.sigma_1, state.sigma_0, state.q_cap, state.n_out_ob_i -  state.n_out_ob_e);

        phi = phi_poisson + (state.phi_zero * phi_laplace);

        // Step 2.1: calculation of electric field
        fields.calculate_efield(efield_x, efield_y, phi, wmesh_i, wmesh_e, mesh, electrode_mask);

        // Step 2.2: field weighting
        if(state.step % k_sub == 0) pic.electric_field_at_particles(efield_x_at_p_i, efield_y_at_p_i, efield_x, efield_y, p_i, state.n_active_i, mesh, lpos_i);
        pic.electric_field_at_particles(efield_x_at_p_e, efield_y_at_p_e, efield_x, efield_y, p_e, state.n_active_e, mesh, lpos_e);

        // Step 3: integration of equations of motion
        if(state.step % k_sub == 0) pops.move_i(p_i, state.n_active_i, efield_x_at_p_i, efield_y_at_p_i);
        pops.move_e(p_e, state.n_active_e, efield_x_at_p_e, efield_y_at_p_e);

        // // Step 4: particle loss at boundaries
        if(state.step % k_sub == 0) pops.boundaries_ob_count(p_i, state.n_active_i, lpos_i, state.n_out_ob_i, state.n_out_thr_i);
        pops.boundaries_ob_count(p_e, state.n_active_e, lpos_e, state.n_out_ob_e, state.n_out_thr_e);

        // Step 5: particles injection
        if(state.step % k_sub == 0) pops.add_flux_particles(p_i, state.n_active_i, t_i, v_drift_i, m_i, n_inj_i, k_sub);
        if(inj_model == "constant");  //n_inj_el = n_inj_el;
        else if(inj_model == "balanced")  n_inj_el = pops.balanced_injection(n_inj_el, 0.01, wmesh_i, wmesh_e, 0, 0, 0, n_thruster - 1);
        else if(inj_model == "pulsed")    n_inj_el = pops.pulsed_injection(state.step);
        else if(inj_model == "square")    n_inj_el = pops.square_injection(state.step);
        pops.add_flux_particles(p_e, state.n_active_e, t_el, v_drift_el, m_el, n_inj_el, 1);

        // Step 6: Monte-Carlo collisions
        if(mcc_coll){
            if(state.step % k_sub == 0) coll.collisions_i(p_i, state.n_active_i, lpos_i, mesh, dens_n);
            coll.collisions_e(p_e, state.n_active_e, lpos_e, p_i, state.n_active_i, lpos_i, mesh, dens_n);
        }
        
        //  ----------------------------- Diagnostics -------------------------

        print_info(state, 100, config);

         if(state.step % 1 == 0){
             series["time"].val[n_points_series] = state.step * dt;
             series["v_cap"].val[n_points_series] = state.phi_zero / k_phi;
             series["n_active_i"].val[n_points_series] = state.n_active_i;
             series["n_active_e"].val[n_points_series] = state.n_active_e;
            
             series["I_out_ob_e"].val[n_points_series] = state.n_out_ob_e * n_factor * q / dt;
             series["I_out_ob_i"].val[n_points_series] = state.n_out_ob_i * n_factor * q / dt;
             series["I_out_thr_e"].val[n_points_series] = state.n_out_thr_e * n_factor * q / dt;
             series["I_out_thr_i"].val[n_points_series] = state.n_out_thr_i * n_factor * q / dt;
             series["I_in_thr_e"].val[n_points_series] = n_inj_el * n_factor * q / dt;
             series["I_in_thr_i"].val[n_points_series] = n_inj_i * n_factor * q / dt;
             n_points_series += 1;
         }
        
        if(state.step % 10000 == 0)
        {
            save_state(p_e, p_i, state, config, job_suffix);
            save_fields_snapshot(phi, wmesh_e, wmesh_i, mesh, state, config, job_suffix);
        }

         if(state.step % 50000 == 0) {
             save_series(series, n_points_series, state, config, job_suffix);
         }
        
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

         if(state.step - state.step_offset > 0.0 * (double) n_steps){
             int i_av = average_field_over_period(phi_av, phi, rf_period_i, n_steps, state.step - state.step_offset);
             average_field(wmesh_e_av, wmesh_e, state.step - state.step_offset);
             average_field(wmesh_i_av, wmesh_i, state.step - state.step_offset);
             if(i_av == rf_period_i){
                 save_fields_snapshot(phi_av, wmesh_e_av, wmesh_i_av, mesh, state, config, "_av" + job_suffix);
             }
         }
        
//        print_fmatrix(phi, 10, 25);
//        system("read -p 'Press Enter to continue...' var");
        
	}

    auto stop = now();
    std::cout << "Total execution duration: " << tdiff(start, stop) / 60 << " min" << endl;

    // ----------------------------- Saving outputs ---------------------------
    save_state(p_e, p_i, state, config, job_suffix);
    save_fields_snapshot(phi, wmesh_e, wmesh_i, mesh, state, config, job_suffix);
    save_series(series, n_points_series, state, config, job_suffix);

    // ----------------------------- Finalizing -------------------------------
    // delete_cross_sections_arrays();
	return 0;
}
