#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "fmatrix.h"
#include "configuration.h"


class diagnostics{

    private:
        double k_v;
        configuration & config;



    public:
        diagnostics(configuration & _config);
        void velocity_distribution(fmatrix & p, int & n_active, int vcol, int n_v, double v_0, double v_1,  fmatrix & dmesh);
        

};


#endif