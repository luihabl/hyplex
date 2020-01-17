#ifndef STATEINFO_H
#define STATEINFO_H


struct state_info{
    
    // Circuit variables
    double q_cap            = 0;
    double sigma_0          = 0;
    double sigma_1          = 0;
    double phi_zero         = 0;

    // Particle variables
    int n_active_e          = 0;
    int n_active_i          = 0;
    int n_out_e             = 0;
    int n_out_i             = 0;

    // State variables
    int step                = 0;
    int step_offset         = 0;

};

#endif