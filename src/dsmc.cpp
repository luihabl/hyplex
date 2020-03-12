#include "dsmc.h"

#include "fmatrix.h"
#include "fmath.h"
#include "particles.h"
#include "particles-in-mesh.h"
#include "num-tools.h"
#include "fields.h"
#include "configuration.h"
#include "input-output.h"


void run_dsmc(mesh_set & mesh, fmatrix & dens_n, configuration & config, string output_name){
    
    int n_active_n = 0;
    particle_operations pops(config);

    fmatrix p_n             = fmatrix::zeros(config.i("particles/n_max_particles"), 6);
    imatrix lpos_n          = imatrix::zeros(config.i("particles/n_max_particles"), 2);
    fmatrix wmesh_n         = fmatrix::zeros(config.i("geometry/n_mesh_x"), config.i("geometry/n_mesh_y"));
    fmatrix wmesh_n_av      = fmatrix::zeros(config.i("geometry/n_mesh_x"), config.i("geometry/n_mesh_y"));

    const int n_steps_dsmc = config.i("time/n_steps_dsmc");
    const int n_average_dsmc = config.i("time/n_average_dsmc");
    const int k_sub_dsmc = config.i("time/k_sub_dsmc");
    const double t_neutral = config.f("neutrals/t_neutral");
    const double m_i = config.f("ugas/m_i");
    const double n_inj_n = config.f("p/n_inj_n");
    const double a_x = config.f("geometry/a_x");
    const double a_y = config.f("geometry/a_y");
    const double dx = config.f("p/dx");
    const double dy = config.f("p/dy");
    
    for (int i = 0; i < n_steps_dsmc; i++){
        pops.add_flux_particles(p_n, n_active_n, t_neutral, 0, m_i, n_inj_n, k_sub_dsmc);

        if(i > n_steps_dsmc - n_average_dsmc){
            weight(p_n, n_active_n, wmesh_n, mesh, lpos_n, a_x, a_y, dx, dy);
            average_field(wmesh_n_av, wmesh_n, i - (n_steps_dsmc - n_average_dsmc));
        }
        
        pops.move_n(p_n, n_active_n, k_sub_dsmc);
        pops.boundaries_n(p_n, n_active_n, lpos_n);
        print_dsmc_info(i, n_active_n, 1000, n_steps_dsmc);    
    }

    wmesh_n = (config.f("particles/n_factor_dsmc") / config.f("particles/n_factor")) * wmesh_n_av;
    dens_n = (4.0 / pow(config.f("p/dx"), 2)) *  config.f("particles/n_factor") * wmesh_n / mesh.v;

    save_fmatrix(dens_n, config.s("project/output_path") + output_name, "dens_n");
}

