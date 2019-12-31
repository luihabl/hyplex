#include "dsmc.h"

#include "fmatrix.h"
#include "fmath.h"
#include "particles.h"
#include "particles-in-mesh.h"
#include "num-tools.h"
#include "fields.h"
#include "config.h"
#include "input-output.h"


void run_dsmc(fmatrix & mesh_x, fmatrix & mesh_y, fmatrix & vmesh, fmatrix & dens_n){
    
    int n_active_n;
    fmatrix p_n = fmatrix::zeros(N_MAX_PARTICLES, 6);
    imatrix lpos_n          = imatrix::zeros(N_MAX_PARTICLES, 2);
    fmatrix wmesh_n         = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    fmatrix wmesh_n_av      = fmatrix::zeros(N_MESH_X, N_MESH_Y);
    
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
    dens_n = (4.0 / pow(DX, 2)) *  N_FACTOR * wmesh_n / vmesh;
        
    save_to_csv(dens_n, "dens_n.csv");
}

