#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "fmatrix.h"
#include "configuration.h"
#include "state-info.h"
#include <string>

using namespace std;
class diagnostics{

    private:
        double k_v, dt, k_phi, n_factor, q;
        int series_measure_step;
        configuration & config;
        state_info & state;
        void initialize_series();
    
    public:
        int series_size;
        unordered_map<string, fmatrix> gseries;
        unordered_map<string, fmatrix> lseries;
        fmatrix tmp_array;
        tmatrix<string> gseries_keys;
        tmatrix<string> lseries_keys;
        int n_points_series;
        diagnostics(configuration & _config, state_info & _state);
        void velocity_distribution(fmatrix & p, int & n_active, int vcol, double v_0, double v_1,  fmatrix & dmesh);
        void update_series(double n_inj_el, double n_inj_i);

};


#endif
