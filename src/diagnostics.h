#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "fmatrix.h"
#include "configuration.h"
#include "state-info.h"
#include "fields.h"
#include <string>

using namespace std;
class diagnostics{

    private:
        double k_v, dt, k_phi, n_factor, q;
        int series_measure_step, n_mesh_x, n_mesh_y;
        configuration & config;
        state_info & state;

        fmatrix dist_e_x, dist_e_y;
        fmatrix dist_i_x, dist_i_y;
        fmatrix kefield_e, kefield_i;
        fmatrix vfield_e_x, vfield_e_y;
        fmatrix vfield_i_x, vfield_i_y;


        void initialize_series();
    
    public:
        int series_size;
        unordered_map<string, fmatrix> gseries; //global series
        unordered_map<string, fmatrix> lseries; //local series - MPI_Reduce before saving
        fmatrix tmp_array;
        tmatrix<string> gseries_keys;
        tmatrix<string> lseries_keys;

        int n_v_e, n_v_i;
        fmatrix dist_e_global_x, dist_e_global_y;
        fmatrix dist_i_global_x, dist_i_global_y;
        fmatrix vlim_e, vlim_i;

        fmatrix kefield_e_global, kefield_i_global;
        fmatrix vfield_e_x_global, vfield_e_y_global;
        fmatrix vfield_i_x_global, vfield_i_y_global;

        int n_points_series;
        diagnostics(configuration & _config, state_info & _state);
        void velocity_distribution(fmatrix & p, int & n_active, int vcol, double v_0, double v_1,  fmatrix & dmesh);
        void update_distributions(fmatrix & p_e, fmatrix & p_i);
        void update_series(double n_inj_el, double n_inj_i);
        void update_velocity_field(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, fmatrix & wmesh_e_global, fmatrix & wmesh_i_global, imatrix & lpos_e, imatrix & lpos_i);
        void update_energy_field(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, fmatrix & wmesh_e_global, fmatrix & wmesh_i_global, imatrix & lpos_e, imatrix & lpos_i);
};


#endif
