
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
#include "diagnostics.h"
#include "clock.h"

#define DEFAULT_CONFIG_PATH "input/config/config.yaml"

using namespace std;
using namespace std::chrono;

// ----------------------------- Main function --------------------------------
int main(int argc, char* argv[])
{
    MPI_Init(&argc,&argv);
    
    int mpi_rank, mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    double mpi_size_double = (double) mpi_size;

    auto start_utc = sys_now();
    cout << "Hyplex " << GIT_VERSION << endl;
    cout << "Starting at " <<  time_to_string(start_utc) << " UTC" << endl << endl;

    // ------------------- Loading configuration ------------------------------

    argparser arg(argc, argv);

    string job_name = arg.get("name", "");
    string config_path = arg.get("config", DEFAULT_CONFIG_PATH);

    configuration config(arg.get("config", DEFAULT_CONFIG_PATH));
    config.set_job_name(job_name);

    state_info state;

    // ------------------- Variable initialization ----------------------------

    // Loading variables from config
    
    const int n_steps           = config.i("time/n_steps");
    const int k_sub             = config.i("time/k_sub");
    const int n_mesh_x          = config.i("geometry/n_mesh_x");
    const int n_mesh_y          = config.i("geometry/n_mesh_y");
    const int n_thruster        = config.i("geometry/n_thruster");
    const int n_max_particles   = config.i("particles/n_max_particles");
    const int verbosity         = config.i("simulation/verbosity");
    const double m_el           = config.f("electrons/m_el");
    const double t_el           = config.f("electrons/t_el");
    const double v_drift_el     = config.f("electrons/v_drift_el");
    const double v_drift_i      = config.f("ions/v_drift_i");
    const double t_i            = config.f("ions/t_i");
    const double m_i            = config.f("ugas/m_i");
    const double volt_0_norm    = config.f("p/volt_0_norm");
    const double volt_1_norm    = config.f("p/volt_1_norm");
    const bool mcc_coll         = config.b("neutrals/mcc_coll");
    const string inj_model      = config.s("electrons/inj_model");
    const int n_mesh_total      = config.i("geometry/n_mesh_x") * config.i("geometry/n_mesh_y");

    verbose_log("Initializing general variables", verbosity >= 1);
 
    // General field variables

	fmatrix phi             = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix phi_laplace     = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix phi_poisson     = fmatrix::zeros(n_mesh_x, n_mesh_y);
	fmatrix efield_x        = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix efield_y        = fmatrix::zeros(n_mesh_x, n_mesh_y);
    imatrix electrode_mask  = imatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix voltages        = fmatrix::zeros(4);
    double sigma_laplace;
    
    // Doiagnostics variables
    verbose_log("Initializing diagnostics variables", verbosity >= 1);
    fmatrix phi_av          = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix wmesh_e_av      = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix wmesh_i_av      = fmatrix::zeros(n_mesh_x, n_mesh_y);
    
	// Particle 1 - Electrons
    verbose_log("Initializing electrons variables", verbosity >= 1);
	fmatrix p_e             = fmatrix::zeros(n_max_particles / mpi_size, 6);
    imatrix lpos_e          = imatrix::zeros(n_max_particles / mpi_size, 2);
	fmatrix wmesh_e         = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix wmesh_e_global   = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix efield_x_at_p_e = fmatrix::zeros(n_max_particles / mpi_size);
    fmatrix efield_y_at_p_e = fmatrix::zeros(n_max_particles / mpi_size);
    double n_inj_el          = config.f("p/n_inj_el");
	
	// Particle 2 - Ions
    verbose_log("Initializing ions variables", verbosity >= 1);
	fmatrix p_i             = fmatrix::zeros(n_max_particles / mpi_size, 6);
    imatrix lpos_i          = imatrix::zeros(n_max_particles / mpi_size, 2);
	fmatrix wmesh_i         = fmatrix::zeros(n_mesh_x, n_mesh_y);
    fmatrix wmesh_i_global   = fmatrix::zeros(n_mesh_x, n_mesh_y);
	fmatrix efield_x_at_p_i = fmatrix::zeros(n_max_particles / mpi_size);
    fmatrix efield_y_at_p_i = fmatrix::zeros(n_max_particles / mpi_size);
    double n_inj_i          = config.f("p/n_inj_i");
    
    int n_out_ob_i_global = 0, n_out_ob_e_global = 0;
    
    // Particle 3 - Neutrals
    verbose_log("Initializing neutral variables", verbosity >= 1);
    fmatrix dens_n          = fmatrix::zeros(n_mesh_x, n_mesh_y);

    // ----------------------------- Operation objects -------------------------
    verbose_log("Initializing mesh", verbosity >= 1);

    mesh_set mesh(config);
    mesh.init_mesh();

    field_operations fields(config);
    pic_operations pic(config);
    particle_operations pops(config, pic);
    diagnostics diag(config, state);

    // ----------------------------- Solver -----------------------------------
    verbose_log("Initializing solver", verbosity >= 1);

    rsolver solver(mesh, config);
    setup_rsolver(solver, mesh, electrode_mask);

    voltages = {1, 0, 1, 1};
    solver.solve(phi_laplace, voltages, wmesh_i, wmesh_e);

    voltages = {volt_1_norm, volt_0_norm, volt_1_norm, volt_1_norm};
    sigma_laplace = fields.sigma_from_phi(phi_laplace, mesh, wmesh_e, wmesh_i, electrode_mask);

    // ----------------------------- DSMC -------------------------------------
    
    string expm_neutral = config.s("neutrals/expm_neutral");

    if (expm_neutral == "dsmc"){
        verbose_log(" ---- Starting DSMC loop ---- ", verbosity >= 1);
        // run_dsmc(mesh, dens_n, config, "dens_n" + job_name + ".exdir");
        run_dsmc_pump(mesh, dens_n, config, "dens_n" + job_name + ".exdir");
    }
    else if (expm_neutral == "constant") {
        verbose_log(" ---- Setting constant neutral density ---- ", verbosity >= 1);
        dens_n.set_all(config.f("neutrals/n_neutral"));
    } 
    else if (expm_neutral == "load") {
        verbose_log(" ---- Loading neutral density field ---- ", verbosity >= 1);
        load_fmatrix(dens_n, config.s("project/input_path") + "dens_n.exdir", "dens_n");
    }

    // ----------------------------- MCC --------------------------------------
    
    mcc coll(config, pops, pic);
    coll.initialize_mcc(dens_n, mesh.v);
    
    // ----------------------------- Loading initial state -------------------

    if(config.s("simulation/initial_state") == "load")  load_state(p_e, p_i, state, config);

    // ----------------------------- Output manager ---------------------------
    output_manager output(start_utc, state, config, mesh);
    output.save_initial_data();

	// ----------------------------- Main loop --------------------------------

    // Printing initial information
    print_initial_info(coll.p_null_e, coll.p_null_i, config);

    tmatrix<system_clock::time_point> tp(10);

    verbose_log(" ---- Starting main loop ---- ", verbosity >= 1);
    auto start = sys_now();
	for (state.step = state.step_offset; state.step < n_steps + state.step_offset; state.step++)
	{
        tp.val[0] = sys_now();

		// Step 1.0: particle weighting
        if(state.step % k_sub == 0) pic.weight(p_i, state.n_active_i, wmesh_i, mesh, lpos_i);
        pic.weight(p_e, state.n_active_e, wmesh_e, mesh, lpos_e);

        tp.val[1] = sys_now();

        // Step 2.0 integration of Poisson's equation
        
        MPI_Reduce(wmesh_i.val, wmesh_i_global.val, n_mesh_total, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(wmesh_e.val, wmesh_e_global.val, n_mesh_total, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&state.n_out_ob_i, &n_out_ob_i_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&state.n_out_ob_e, &n_out_ob_e_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if(mpi_rank == 0){

            solver.solve(phi_poisson, voltages, wmesh_i_global, wmesh_e_global);

            tp.val[2] = sys_now();

            state.sigma_0 = state.sigma_1;

            state.phi_zero = fields.calculate_phi_zero(state.sigma_1, n_out_ob_i_global - n_out_ob_e_global, state.q_cap, sigma_laplace, phi_poisson, mesh, wmesh_e_global, wmesh_i_global, electrode_mask);

            state.sigma_1 = fields.calculate_sigma(state.sigma_0, state.phi_zero, n_out_ob_i_global - n_out_ob_e_global, state.q_cap);
            state.q_cap = fields.calculate_cap_charge(state.sigma_1, state.sigma_0, state.q_cap, n_out_ob_i_global - n_out_ob_e_global);

            phi = phi_poisson + (state.phi_zero * phi_laplace);

            tp.val[3] = sys_now();

            // Step 2.1: calculation of electric field
            fields.calculate_efield(efield_x, efield_y, phi, wmesh_i_global, wmesh_e_global, mesh, electrode_mask);

        }
        MPI_Bcast(phi.val, n_mesh_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(efield_x.val, n_mesh_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(efield_y.val, n_mesh_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        tp.val[4] = sys_now();
        
        // Step 2.2: field weighting
        if(state.step % k_sub == 0) pic.electric_field_at_particles(efield_x_at_p_i, efield_y_at_p_i, efield_x, efield_y, p_i, state.n_active_i, mesh, lpos_i);
        pic.electric_field_at_particles(efield_x_at_p_e, efield_y_at_p_e, efield_x, efield_y, p_e, state.n_active_e, mesh, lpos_e);

        tp.val[5] = sys_now();

        // Step 3: integration of equations of motion
        if(state.step % k_sub == 0) pops.move_i(p_i, state.n_active_i, efield_x_at_p_i, efield_y_at_p_i);
        pops.move_e(p_e, state.n_active_e, efield_x_at_p_e, efield_y_at_p_e);

        tp.val[6] = sys_now();

        // // Step 4: particle loss at boundaries
        if(state.step % k_sub == 0) pops.boundaries_ob_count(p_i, state.n_active_i, lpos_i, state.n_out_ob_i, state.n_out_thr_i);
        pops.boundaries_ob_count(p_e, state.n_active_e, lpos_e, state.n_out_ob_e, state.n_out_thr_e);

        tp.val[7] = sys_now();

        // Step 5: particles injection
        if(state.step % k_sub == 0) state.n_in_thr_i = pops.add_flux_particles(p_i, state.n_active_i, t_i, v_drift_i, m_i, n_inj_i / mpi_size_double, k_sub);
        else state.n_in_thr_i = 0;
        
        if(inj_model == "constant");  //n_inj_el = n_inj_el;
        else if(inj_model == "balanced")  n_inj_el = pops.balanced_injection(n_inj_el, 0.01, wmesh_i, wmesh_e, 0, 0, 0, n_thruster - 1);
        else if(inj_model == "pulsed")    n_inj_el = pops.pulsed_injection(state.step);
        else if(inj_model == "square")    n_inj_el = pops.square_injection(state.step);
        state.n_in_thr_e = pops.add_flux_particles(p_e, state.n_active_e, t_el, v_drift_el, m_el, n_inj_el / mpi_size_double, 1);
        
        // Step 6: Monte-Carlo collisions
        if(mcc_coll){
            if(state.step % k_sub == 0) coll.collisions_i(p_i, state.n_active_i, lpos_i, mesh, dens_n);
            coll.collisions_e(p_e, state.n_active_e, lpos_e, p_i, state.n_active_i, lpos_i, mesh, dens_n);
        }

        tp.val[8] = sys_now();

        //  ----------------------------- Diagnostics -------------------------

        output.print_info();

        diag.update_series(n_inj_el, n_inj_i);
        output.save_state(p_e, p_i);
        output.save_fields_snapshot(phi, wmesh_e_global, wmesh_i_global, mesh, "");
        output.save_series(diag);
        output.save_distributions(diag, p_e, p_i);
        output.update_metadata();
        output.fields_rf_average(phi, wmesh_e_global, wmesh_i_global, mesh);

        tp.val[9] = sys_now();

        output.print_loop_timing(tp);
    }

    auto stop = sys_now();
    if(mpi_rank==0) std::cout << "Total execution duration: " << tdiff_h(start, stop) << " hours" << endl;

    // ----------------------------- Saving outputs ---------------------------
    output.save_state(p_e, p_i, config.b("diagnostics/state/end_save"));
    output.save_fields_snapshot(phi, wmesh_e_global, wmesh_i_global, mesh, "", config.b("diagnostics/fields_snapshot/end_save"));
    output.save_series(diag, config.b("diagnostics/series/end_save"));
    output.save_distributions(diag, p_e, p_i, config.b("diagnostics/vdist/end_save"));
    output.update_metadata("completed", config.b("diagnostics/metadata/end_save"));
    
   
	return 0;
}
