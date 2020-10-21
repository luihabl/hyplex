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
        double k_v, dt, k_phi, n_factor, q, izfield_start_progress;
        int series_measure_step, n_mesh_x, n_mesh_y, step_save_vdist, step_save_fields, steps_since_last_izfield_save, n_steps;
        configuration & config;
        state_info & state;

        fmatrix wmesh_e, wmesh_i;
        fmatrix dist_e_x, dist_e_y;
        fmatrix dist_i_x, dist_i_y;
        fmatrix kfield_e, kfield_i;
        fmatrix ufield_e_x, ufield_e_y;
        fmatrix ufield_i_x, ufield_i_y;
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

        fmatrix wmesh_e_global, wmesh_i_global;
        fmatrix kfield_e_global, kfield_i_global;
        fmatrix ufield_e_x_global, ufield_e_y_global;
        fmatrix ufield_i_x_global, ufield_i_y_global;
        fmatrix izfield_global;

        int n_points_series;
        diagnostics(configuration & _config, state_info & _state);
        void velocity_distribution(fmatrix & p, int & n_active, int vcol, double v_0, double v_1,  fmatrix & dmesh);
        void update_distributions(fmatrix & p_e, fmatrix & p_i);
        void update_series(double n_inj_el, double n_inj_i);
        void update_ufield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, imatrix & lpos_e, imatrix & lpos_i, bool force = false);
        void update_kfield(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, imatrix & lpos_e, imatrix & lpos_i, bool force = false);
        void update_izfield(mesh_set & mesh, fmatrix & p_i, imatrix & lpos_i, int n_iz);
        void izfield_set_zero();
        void update_internal_wmesh(mesh_set & mesh, fmatrix & p_e, fmatrix & p_i, imatrix & lpos_e, imatrix & lpos_i, bool force = false);
        void set_internal_wmesh(fmatrix & _wmesh_e_global, fmatrix & _wmesh_i_global);
};


#endif
