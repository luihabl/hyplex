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
        int series_measure_step, n_mesh_x, n_mesh_y, step_save_vdist, step_save_fields;
        configuration & config;
        state_info & state;

        fmatrix dist_e_x, dist_e_y;
        fmatrix dist_i_x, dist_i_y;
        fmatrix pfield_e, pfield_i;
        fmatrix ffield_e_x, ffield_e_y;
        fmatrix ffield_i_x, ffield_i_y;
        fmatrix izfield;


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

        fmatrix pfield_e_global, pfield_i_global;
        fmatrix ffield_e_x_global, ffield_e_y_global;
        fmatrix ffield_i_x_global, ffield_i_y_global;
        fmatrix izfield_global;

        int n_points_series;
        diagnostics(configuration & _config, state_info & _state);
        void velocity_distribution(fmatrix & p, int & n_active, int vcol, double v_0, double v_1,  fmatrix & dmesh);
        void update_distributions(fmatrix & p_e, fmatrix & p_i);
        void update_series(double n_inj_el, double n_inj_i);
        void update_ffield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, fmatrix & wmesh_e_global, fmatrix & wmesh_i_global, imatrix & lpos_e, imatrix & lpos_i);
        void update_pfield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, fmatrix & wmesh_e_global, fmatrix & wmesh_i_global, imatrix & lpos_e, imatrix & lpos_i);
        void update_izfield(mesh_set & mesh, fmatrix & p_i, imatrix & lpos_i, int n_iz);
};


#endif
