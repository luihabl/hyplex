#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "fmatrix.h"



class diagnostics{

    



    public:
        diagnostics();
        void velocity_distribution(fmatrix & p, int & n_active, int & vcol, int n_v, double v_0, double v_1,  fmatrix & dmesh);

};








#endif