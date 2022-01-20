#ifndef RSOLVER_H
#define RSOLVER_H

#include <iostream>
#include "mpi.h"
#include "HYPRE_struct_ls.h"
#include "fmatrix.h"
#include "fields.h"
#include "configuration.h"

class rsolver{
    
private:
    // private variables
    HYPRE_StructGrid     hypre_grid;
    HYPRE_StructStencil  hypre_stencil;
    HYPRE_StructMatrix   hypre_A;
    HYPRE_StructVector   hypre_b;
    HYPRE_StructVector   hypre_x;
    HYPRE_StructSolver   hypre_solver;
    imatrix              dirichlet_boxes, neumann_boxes, node_type;
    imatrix              ll, ur, inner_ll, inner_ur;
    fmatrix              dirichlet_input;
    mesh_set *           mesh;
    configuration *      config;
    int                  n_mesh_x, n_mesh_y, n_solve;
    int                  n_dirichlet, n_input_dirichlet, n_neumann;
    
    int stencil_indices[5] = {0, 1, 2, 3, 4};
    int stencil_offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

    // private methods
    void init_grid();
    void init_matrix();
    void init_vectors();
    void find_node_types();
    void set_stencils();
    void set_inner_nodes();
    void set_dirichlet_nodes();
    void set_neumann_nodes();
    void add_dirichlet_input(int i, int j, int stencil, int n_bc);
    bool in_box(int i, int j, int ill, int jll, int iur, int jur);
    int dirichlet_boundary_number(int i, int j);
    
public:
    // constructors and destructors
    rsolver() = default;
    rsolver(mesh_set * mesh, int n_neumann, int n_dirichlet, configuration * config);
    rsolver(mesh_set * mesh, configuration * config);
    ~rsolver();
    
    // public methods
    void late_init(int n_neumann, int n_dirichlet);
    void assemble();
    void solve(fmatrix & solution, fmatrix & voltages, fmatrix & w_i, fmatrix & w_e);
    void set_dirichlet_box(imatrix & box, int counter);
    void set_neumann_box(imatrix & box, int counter);
    int get_node_type(int i, int j, int ioff=0, int joff=0);

    void setup(mesh_set & mesh, imatrix & electrode_mask);
    void setup_benchmark(mesh_set & mesh, imatrix & electrode_mask);

};

#endif
