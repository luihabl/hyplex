#include "simulation.h"
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


simulation::simulation(mpi_info & _mpi, configuration & _config) : config(_config), mpi(_mpi) {
    setup();
}

void simulation::setup() {

    auto start_utc = clk::sys_now();
    cout << "Hyplex " << GIT_VERSION << endl;
    cout << "Starting at " <<  clk::time_to_string(start_utc) << " UTC" << endl << endl;

    is_benchmark = config.s("boundaries/ob_type") == "benchmark";

    n_steps           = config.i("time/n_steps");
    k_sub             = config.i("time/k_sub");
    n_mesh_x          = config.i("geometry/n_mesh_x");
    n_mesh_y          = config.i("geometry/n_mesh_y");
    n_thruster        = config.i("geometry/n_thruster");
    n_max_particles   = config.i("particles/n_max_particles");
    verbosity         = config.i("simulation/verbosity");
    m_el              = config.f("electrons/m_el");
    t_el              = config.f("electrons/t_el");
    v_drift_el        = config.f("electrons/v_drift_el");
    v_drift_i         = config.f("ions/v_drift_i");
    t_i               = config.f("ions/t_i");
    m_i               = config.f("ugas/m_i");
    volt_0_norm       = config.f("p/volt_0_norm");
    volt_1_norm       = config.f("p/volt_1_norm");
    mcc_coll          = config.b("neutrals/mcc_coll");
    inj_model         = config.s("electrons/inj_model");
    n_mesh_total      = config.i("geometry/n_mesh_x") * config.i("geometry/n_mesh_y");
    dt                = config.f("time/dt");
    freq              = config.f("thruster/freq");

    io::verbose_log("Initializing general variables", verbosity >= 1);
 
    // General field variables

    phi             = fmatrix::zeros(n_mesh_x, n_mesh_y);
    phi_laplace     = fmatrix::zeros(n_mesh_x, n_mesh_y);
    phi_poisson     = fmatrix::zeros(n_mesh_x, n_mesh_y);
    efield_x        = fmatrix::zeros(n_mesh_x, n_mesh_y);
    efield_y        = fmatrix::zeros(n_mesh_x, n_mesh_y);
    electrode_mask  = imatrix::zeros(n_mesh_x, n_mesh_y);
    voltages        = fmatrix::zeros(4);
    
    // Doiagnostics variables
    io::verbose_log("Initializing diagnostics variables", verbosity >= 1);

    phi_av          = fmatrix::zeros(n_mesh_x, n_mesh_y);
    wmesh_e_av      = fmatrix::zeros(n_mesh_x, n_mesh_y);
    wmesh_i_av      = fmatrix::zeros(n_mesh_x, n_mesh_y);
    
    // Particle 1 - Electrons
    io::verbose_log("Initializing electrons variables", verbosity >= 1);

    p_e             = fmatrix::zeros(n_max_particles / mpi.size, 6);
    lpos_e          = imatrix::zeros(n_max_particles / mpi.size, 2);
    wmesh_e         = fmatrix::zeros(n_mesh_x, n_mesh_y);
    wmesh_e_global  = fmatrix::zeros(n_mesh_x, n_mesh_y);
    efield_x_at_p_e = fmatrix::zeros(n_max_particles / mpi.size);
    efield_y_at_p_e = fmatrix::zeros(n_max_particles / mpi.size);
    n_inj_el        = config.f("p/n_inj_el");
        
    // Particle 2 - Ions
    io::verbose_log("Initializing ions variables", verbosity >= 1);

    p_i             = fmatrix::zeros(n_max_particles / mpi.size, 6);
    lpos_i          = imatrix::zeros(n_max_particles / mpi.size, 2);
    wmesh_i         = fmatrix::zeros(n_mesh_x, n_mesh_y);
    wmesh_i_global  = fmatrix::zeros(n_mesh_x, n_mesh_y);
    efield_x_at_p_i = fmatrix::zeros(n_max_particles / mpi.size);
    efield_y_at_p_i = fmatrix::zeros(n_max_particles / mpi.size);
    n_inj_i         = config.f("p/n_inj_i");
    
    n_out_ob_i_global = 0;
    n_out_ob_e_global = 0;
    
    // Particle 3 - Neutrals
    io::verbose_log("Initializing neutral variables", verbosity >= 1);
    dens_n          = fmatrix::zeros(n_mesh_x, n_mesh_y);

    // ----------------------------- Operation objects -------------------------
    io::verbose_log("Initializing mesh", verbosity >= 1);

    mesh = mesh_set(config);
    mesh.init_mesh();

    fields = field_operations(config);
    pic = pic_operations(config);
    pops = particle_operations(config, &pic);
    diag = diagnostics(&config, &state);

    // ----------------------------- Solver -----------------------------------


    output = io::output_manager(start_utc, &state, &config, &mesh);
}

void simulation::gen_neutral_field() {
    
    expm_neutral = config.s("neutrals/expm_neutral");

    if (expm_neutral == "dsmc"){
        io::verbose_log(" ---- Starting DSMC loop ---- ", verbosity >= 1);
        cout << "Process " << mpi.rank << " starting DSMC loop" << endl;
        dsmc::run_dsmc(mesh, dens_n, config);
    }
    else if (expm_neutral == "constant") {
        io::verbose_log(" ---- Setting constant neutral density ---- ", verbosity >= 1);
        dens_n.set_all(config.f("neutrals/n_neutral"));
    } 
    else if (expm_neutral == "load") {
        io::verbose_log(" ---- Loading neutral density field ---- ", verbosity >= 1);
        io::load_fmatrix(dens_n, config.s("project/input_path") + "dens_n.exdir", "dens_n");
    }

    coll = mcc(config, &pops, &pic);
    coll.initialize_mcc(dens_n, mesh.v);
}

void simulation::load_state() {
    if(config.s("simulation/initial_state") == "load")  {
        io::load_state(p_e, p_i, state, config);
    }
}

void simulation::run() {

    io::verbose_log("Initializing field solver", verbosity >= 1);

    rsolver solver(&mesh, &config);
    solver.setup(mesh, electrode_mask);

    if(is_benchmark)
        voltages = {0, 1};
    else
        voltages = {1, 0, 1, 1};

    solver.solve(phi_laplace, voltages, wmesh_i, wmesh_e);

    if(is_benchmark)
        voltages = {0, 0};
    else
        voltages = {volt_1_norm, volt_0_norm, volt_1_norm, volt_1_norm};
    
    sigma_laplace = fields.sigma_from_phi(phi_laplace, mesh, wmesh_e, wmesh_i, electrode_mask);

    if(is_benchmark) { //load initial particles
        pops.add_maxwellian_particles(p_e, state.n_active_e, t_el, m_el, (size_t) ceil(config.f("benchmark/n_total") / (double) mpi.size));
        pops.add_maxwellian_particles(p_i, state.n_active_i, t_i, m_i, (size_t) ceil(config.f("benchmark/n_total") / (double) mpi.size));
    }

    // Printing initial information
    io::print_initial_info(coll.p_null_e, coll.p_null_i, config);

    tmatrix<system_clock::time_point> tp(10);

    io::verbose_log(" ---- Starting main loop ---- ", verbosity >= 1);
    auto start = clk::sys_now();
	for (state.step = state.step_offset; state.step < n_steps + state.step_offset; state.step++)
	{
        tp.val[0] = clk::sys_now();

		// Step 1.0: particle weighting
        if(state.step % k_sub == 0) pic_operations::weight(p_i, state.n_active_i, wmesh_i, mesh, lpos_i);
        pic_operations::weight(p_e, state.n_active_e, wmesh_e, mesh, lpos_e);

        tp.val[1] = clk::sys_now();

        // Step 2.0 integration of Poisson's equation
        
        MPI_Reduce(wmesh_i.val, wmesh_i_global.val, n_mesh_total, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(wmesh_e.val, wmesh_e_global.val, n_mesh_total, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&state.n_out_ob_i, &n_out_ob_i_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&state.n_out_ob_e, &n_out_ob_e_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if(mpi.rank == 0){

            solver.solve(phi_poisson, voltages, wmesh_i_global, wmesh_e_global);

            tp.val[2] = clk::sys_now();

            if(is_benchmark)
            {
                state.phi_zero = volt_1_norm * std::sin(2 * M_PI * freq * dt * (double) state.step);
            }
            else
            {      
                state.sigma_0 = state.sigma_1;

                state.phi_zero = fields.calculate_phi_zero(state.sigma_1, n_out_ob_i_global - n_out_ob_e_global, state.q_cap, sigma_laplace, phi_poisson, mesh, wmesh_e_global, wmesh_i_global, electrode_mask);

                state.sigma_1 = fields.calculate_sigma(state.sigma_0, state.phi_zero, n_out_ob_i_global - n_out_ob_e_global, state.q_cap);
                state.q_cap = fields.calculate_cap_charge(state.sigma_1, state.sigma_0, state.q_cap, n_out_ob_i_global - n_out_ob_e_global);
            }

            phi = phi_poisson + (state.phi_zero * phi_laplace);

            tp.val[3] = clk::sys_now();

            // Step 2.1: calculation of electric field
            fields.calculate_efield(efield_x, efield_y, phi, wmesh_i_global, wmesh_e_global, mesh, electrode_mask);

        }
        MPI_Bcast(phi.val, n_mesh_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(efield_x.val, n_mesh_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(efield_y.val, n_mesh_total, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        tp.val[4] = clk::sys_now();
        
        // Step 2.2: field weighting
        if(state.step % k_sub == 0) pic.electric_field_at_particles(efield_x_at_p_i, efield_y_at_p_i, efield_x, efield_y, p_i, state.n_active_i, mesh, lpos_i);
        pic.electric_field_at_particles(efield_x_at_p_e, efield_y_at_p_e, efield_x, efield_y, p_e, state.n_active_e, mesh, lpos_e);

        tp.val[5] = clk::sys_now();

        // Step 3: integration of equations of motion
        if(state.step % k_sub == 0) pops.move_i(p_i, state.n_active_i, efield_x_at_p_i, efield_y_at_p_i);
        pops.move_e(p_e, state.n_active_e, efield_x_at_p_e, efield_y_at_p_e);

        tp.val[6] = clk::sys_now();

        // // Step 4: particle loss at boundaries
        if(is_benchmark)
        {
            if(state.step % k_sub == 0) pops.boundaries_benckmark(p_i, state.n_active_i, lpos_i);
            pops.boundaries_benckmark(p_e, state.n_active_e, lpos_e);
        }
        else
        {
            if(state.step % k_sub == 0) pops.boundaries_ob_count(p_i, state.n_active_i, lpos_i, state.n_out_ob_i, state.n_out_thr_i, diag.p_i_removed, diag.n_removed_i, true);
            pops.boundaries_ob_count(p_e, state.n_active_e, lpos_e, state.n_out_ob_e, state.n_out_thr_e, diag.p_e_removed, diag.n_removed_e, true);
        }
        
        tp.val[7] = clk::sys_now();

        // Step 5: particles injection

        if(!is_benchmark)
        {
            if(state.step % k_sub == 0) state.n_in_thr_i = pops.add_flux_particles(p_i, state.n_active_i, t_i, v_drift_i, m_i, n_inj_i / mpi.size_double, k_sub);
            else state.n_in_thr_i = 0;
            
            if(inj_model == "constant");  //n_inj_el = n_inj_el;
            else if(inj_model == "balanced")  n_inj_el = pops.balanced_injection(n_inj_el, 0.01, wmesh_i, wmesh_e, 0, 0, 0, n_thruster - 1);
            else if(inj_model == "pulsed")    n_inj_el = pops.pulsed_injection(state.step);
            else if(inj_model == "square")    n_inj_el = pops.square_injection(state.step);
            state.n_in_thr_e = pops.add_flux_particles(p_e, state.n_active_e, t_el, v_drift_el, m_el, n_inj_el / mpi.size_double, 1);
        }

        // Step 6: Monte-Carlo collisions
        if(mcc_coll){
            if(state.step % k_sub == 0) coll.collisions_i(p_i, state.n_active_i, lpos_i, mesh, dens_n);
            diag.n_e_iz = coll.collisions_e(p_e, state.n_active_e, lpos_e, p_i, state.n_active_i, lpos_i, mesh, dens_n);
        }

        tp.val[8] = clk::sys_now();

        //  ----------------------------- Diagnostics -------------------------

        output.refresh_file();
        output.print_info();
        diag.update_all(mesh, p_e, p_i, lpos_e, lpos_i); // <- add starting step and saving step for boundary current distribution.
                                                         // <- remove steps from the specific functions in update_all and put them inside update_all, this will make it easier to call specific functions individually for field series
        output.save_state(p_e, p_i);
        output.save_fields_snapshot(phi, wmesh_e_global, wmesh_i_global, diag, mesh, ""); 
        output.save_series(diag);
        output.save_distributions(diag);
        output.update_metadata();
        output.fields_average(phi, wmesh_e_global, wmesh_i_global, mesh);
        
        tp.val[9] = clk::sys_now();
        output.print_loop_timing(tp);
    }

    auto stop = clk::sys_now();
    if(mpi.rank==0) cout << "Total execution duration: " << clk::tdiff_h(start, stop) << " hours" << endl;

    // ----------------------------- Saving outputs ---------------------------
    diag.update_internal_wmesh(mesh, p_e, p_i, lpos_e, lpos_i, config.b("diagnostics/fields_snapshot/end_save")); // <- it is bad that this is not included in "update_all"
    diag.update_ufield(mesh, p_e, p_i, lpos_e, lpos_i, config.b("diagnostics/fields_snapshot/end_save")); // <- it is bad that this is not included in "update_all"
    diag.update_kfield(mesh, p_e, p_i, lpos_e, lpos_i, config.b("diagnostics/fields_snapshot/end_save")); // <- it is bad that this is not included in "update_all"
    
    output.save_state(p_e, p_i, config.b("diagnostics/state/end_save"));
    output.save_fields_snapshot(phi, wmesh_e_global, wmesh_i_global, diag, mesh, "",  config.b("diagnostics/fields_snapshot/end_save"));
    output.save_series(diag, config.b("diagnostics/series/end_save"));
    output.save_distributions(diag, config.b("diagnostics/vdist/end_save"));
    output.update_metadata("completed", config.b("diagnostics/metadata/end_save"));

}