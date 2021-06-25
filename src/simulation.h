#ifndef SIMULATION_H
#define SIMULATION_H

#include "configuration.h"
#include "fmatrix.h"
#include "fields.h"
#include "rsolver.h"
#include "state-info.h"
#include "particles-in-mesh.h"
#include "mcc.h"
#include "input-output.h"

struct mpi_info {
    int size;
    double size_double;
    int rank;
};

class simulation{

    public:
        // ---- Constructors / destructors
        simulation(mpi_info & mpi, configuration & config);

        // ---- Methods
        void run();
        void gen_neutral_field();
        void setup();
        void load_state();

        // --- Auxiliary objects

        field_operations fields;
        pic_operations pic;
        particle_operations pops;
        diagnostics diag;
        mcc coll;

        configuration & config;
        state_info state;
        mpi_info mpi;

        io::output_manager output;

        // ---- Simulation variables

        int n_steps;        
        int k_sub;          
        int n_mesh_x;       
        int n_mesh_y;       
        int n_thruster;     
        int n_max_particles;
        int verbosity;      
        double m_el;        
        double t_el;        
        double v_drift_el;  
        double v_drift_i;   
        double t_i;         
        double m_i;         
        double volt_0_norm; 
        double volt_1_norm; 
        bool mcc_coll;      
        string inj_model;   
        int n_mesh_total;   

        fmatrix phi;           
        fmatrix phi_laplace;   
        fmatrix phi_poisson;   
        fmatrix efield_x;      
        fmatrix efield_y;      
        imatrix electrode_mask;
        fmatrix voltages;      
        double sigma_laplace;

        fmatrix phi_av;
        fmatrix wmesh_e_av;
        fmatrix wmesh_i_av;

        fmatrix p_e;
        imatrix lpos_e;
        fmatrix wmesh_e;
        fmatrix wmesh_e_global;
        fmatrix efield_x_at_p_e;
        fmatrix efield_y_at_p_e;
        double n_inj_el;
            
        fmatrix p_i;
        imatrix lpos_i;
        fmatrix wmesh_i;
        fmatrix wmesh_i_global;
        fmatrix efield_x_at_p_i;
        fmatrix efield_y_at_p_i;
        double n_inj_i;

        int n_out_ob_i_global;
        int n_out_ob_e_global;

        fmatrix dens_n;

        mesh_set mesh;

        string expm_neutral;


};





#endif