#include "dsmc.h"

#include "fmatrix.h"
#include "fmath.h"
#include "particles.h"
#include "particles-in-mesh.h"
#include "num-tools.h"
#include "fields.h"
#include "configuration.h"
#include "input-output.h"

#include "mpi.h"


void run_dsmc(mesh_set & mesh, fmatrix & dens_n, configuration & config, string output_name){
    
    int n_active_n = 0;
    pic_operations pic(config);
    particle_operations pops(config, pic);

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

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    for (int i = 0; i < n_steps_dsmc; i++){
        pops.add_flux_particles(p_n, n_active_n, t_neutral, 0, m_i, n_inj_n / (double) size, k_sub_dsmc);

        if(i > n_steps_dsmc - n_average_dsmc){
            pic.weight(p_n, n_active_n, wmesh_n, mesh, lpos_n);
            average_field(wmesh_n_av, wmesh_n, i - (n_steps_dsmc - n_average_dsmc));
        }
        
        pops.move_n(p_n, n_active_n, k_sub_dsmc);
        pops.boundaries_n(p_n, n_active_n, lpos_n);
        if(rank == 0) print_dsmc_info(i, n_active_n * size, 1000, n_steps_dsmc);
    }

    MPI_Allreduce(wmesh_n_av.val, wmesh_n.val, config.i("geometry/n_mesh_x") * config.i("geometry/n_mesh_y"), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dens_n = config.f("particles/n_factor_dsmc") *  (4.0 / pow(config.f("geometry/dx"), 2)) * wmesh_n / mesh.v;

    if(rank == 0){
        save_fmatrix(dens_n, config.s("project/output_path") + output_name, "dens_n");
    }
}

