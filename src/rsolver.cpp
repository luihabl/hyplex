
#include "rsolver.h"

#include <iostream>
#include "mpi.h"
#include "HYPRE_struct_ls.h"
#include "fmatrix.h"
#include "fields.h"
#include "config.h"
#include "util.h"

#include <chrono>
using namespace std::chrono;

using namespace std;
// ------------------- constructors / destructor -------------------------

rsolver::rsolver(fmatrix & _mesh_x, fmatrix & _mesh_y, fmatrix & _vmesh, int _n_neumann, int _n_dirichlet): mesh_x(_mesh_x), mesh_y(_mesh_y), vmesh(_vmesh){

    MPI_Init(NULL, NULL);

    n_neumann = _n_neumann;
    n_dirichlet = _n_dirichlet;
    n_mesh_x = (int) _mesh_x.n1;
    n_mesh_y = (int) _mesh_y.n2;
    n_input_dirichlet = 0;
    n_solve = 0;

    ll = ur = inner_ll = inner_ur = imatrix::zeros(2);
    dirichlet_boxes = imatrix::zeros(n_dirichlet, 4);
    neumann_boxes = imatrix::zeros(n_neumann, 4);
    node_type = imatrix::zeros(n_mesh_x, n_mesh_y);
    dirichlet_input = fmatrix::zeros(2 * n_mesh_x + 2 * n_mesh_y, 4);

    inner_ll = {1, 1};
    inner_ur = {n_mesh_x - 2, n_mesh_y - 2};
}

rsolver::~rsolver(){
    HYPRE_StructGridDestroy(hypre_grid);
    HYPRE_StructStencilDestroy(hypre_stencil);
    HYPRE_StructMatrixDestroy(hypre_A);
    HYPRE_StructVectorDestroy(hypre_b);
    HYPRE_StructVectorDestroy(hypre_x);
    MPI_Finalize();
}


// ------------------- initialization ----------------------------

void rsolver::init_grid(){
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &hypre_grid);
    HYPRE_StructGridSetExtents(hypre_grid, inner_ll.val, inner_ur.val);
    n_solve += (inner_ur.val[0] - inner_ll.val[0] + 1) * (inner_ur.val[1] - inner_ll.val[1] + 1);

    // Neumann boundary points
    for (int i = 0; i < n_neumann; i++){
        
        ll = {neumann_boxes.val[i * 4 + 0], neumann_boxes.val[i * 4 + 1]};
        ur = {neumann_boxes.val[i * 4 + 2], neumann_boxes.val[i * 4 + 3]};
        
        HYPRE_StructGridSetExtents(hypre_grid, ll.val, ur.val);
        n_solve += (ur.val[0] - ll.val[0] + 1) * (ur.val[1] - ll.val[1] + 1);
    }
    HYPRE_StructGridAssemble(hypre_grid);
}

void rsolver::init_matrix(){
    HYPRE_StructStencilCreate(2, 5, &hypre_stencil);
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, hypre_grid, hypre_stencil, &hypre_A);
    HYPRE_StructMatrixInitialize(hypre_A);
}

void rsolver::init_vectors(){
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, &hypre_b);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, hypre_grid, &hypre_x);
    HYPRE_StructVectorInitialize(hypre_b);
    HYPRE_StructVectorInitialize(hypre_x);
}

void rsolver::assemble(){
    init_grid();
    find_node_types();
    init_matrix();
    set_stencils();
    init_vectors();
}

// ------------------- setting boundaries -------------------------

void rsolver::set_dirichlet_box(imatrix & box, int counter){
    for (int j = 0; j < 4; j++)
        dirichlet_boxes.val[counter * 4 + j] = box.val[j];
}

void rsolver::set_neumann_box(imatrix & box, int counter){
    for (int j = 0; j < 4; j++)
        neumann_boxes.val[counter * 4 + j] = box.val[j];
}

void rsolver::find_node_types(){
    for (int n = 0; n < n_neumann; n++)
        node_type.setbox_value(2, neumann_boxes.val[n * 4 + 0], neumann_boxes.val[n * 4 + 1], neumann_boxes.val[n * 4 + 2], neumann_boxes.val[n * 4 + 3]);
    for (int n = 0; n < n_dirichlet; n++)
        node_type.setbox_value(1, dirichlet_boxes.val[n * 4 + 0], dirichlet_boxes.val[n * 4 + 1], dirichlet_boxes.val[n * 4 + 2], dirichlet_boxes.val[n * 4 + 3]);
}

int rsolver::get_node_type(int i, int j, int ioff, int joff){
    size_t ipos = i + ioff;
    size_t jpos = j + joff;

    if((ipos < 0 || ipos > node_type.n1 - 1) || (jpos < 0 || jpos > node_type.n2 - 1)) return -1;
    else return node_type.val[ipos * node_type.n2 + jpos];
}

// --------------------- setting stencils --------------------

void rsolver::set_stencils(){

    for (int n = 0; n < 5; n++)
        HYPRE_StructStencilSetElement(hypre_stencil, n, stencil_offsets[n]);

    set_inner_nodes();
    set_neumann_nodes();
    set_dirichlet_nodes();
    
    HYPRE_StructMatrixAssemble(hypre_A);
}

void rsolver::set_inner_nodes(){
    for (int i = 0; i < n_mesh_x; i++)
        for (int j = 0; j < n_mesh_y; j++) { 
            if(get_node_type(i, j) == 0){

                double vol = vmesh.val[i * vmesh.n2 + j];
                double d0 =   1/ (k1_x(mesh_x, i, j) * k2_x(mesh_x, i, j)) + 1/(k1_y(mesh_y, i, j) * k2_y(mesh_y, i, j));
                double d1 = - 1/ (k1_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
                double d2 = - 1/ (k2_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
                double d3 = - 1/ (k1_y(mesh_y, i, j) * k3_y(mesh_y, i, j));
                double d4 = - 1/ (k2_y(mesh_y, i, j) * k3_y(mesh_y, i, j)); 
                double val[5] = {vol * d0, vol * d1, vol * d2, vol * d3, vol * d4};               

                ll = {i, j};
                HYPRE_StructMatrixSetValues(hypre_A, ll.val, 5, stencil_indices, val);
            }            
        }
}

void rsolver::set_dirichlet_nodes(){
    double val[] = {0.0};
    for (int i = 0; i < n_mesh_x; i++)
        for (int j = 0; j < n_mesh_y; j++)  
            if(get_node_type(i, j) == 1) 
                for (int k = 1; k < 5; k++) 
                    if (get_node_type(i, j, stencil_offsets[k][0], stencil_offsets[k][1]) == 0 || get_node_type(i, j, stencil_offsets[k][0], stencil_offsets[k][1]) == 2) {
                      
                        int stencil_index[] = {k - (stencil_offsets[k][0] + stencil_offsets[k][1])};
                        ll = ur = {i + stencil_offsets[k][0], j + stencil_offsets[k][1]};
                        
                        add_dirichlet_input(i + stencil_offsets[k][0], j + stencil_offsets[k][1], stencil_index[0], dirichlet_boundary_number(i, j));
                      
                        HYPRE_StructMatrixSetBoxValues(hypre_A, ll.val, ur.val, 1, stencil_index, val);
                    }
}


void rsolver::set_neumann_nodes(){
    int k = 0;
    for (int i = 0; i < n_mesh_x; i++)
        for (int j = 0; j < n_mesh_y; j++)
            
            if(get_node_type(i, j) == 2){

                // Set neumann for each possible direction
                k = 1;
                if (get_node_type(i, j, stencil_offsets[k][0], stencil_offsets[k][1]) == 0) {

                    double vol = vmesh.val[i * vmesh.n2 + j];

                    double dx0 =   1 / (k1_x(mesh_x, i, j) * k1_x(mesh_x, i, j)) + 1/(k1_y(mesh_y, i, j) * k2_y(mesh_y, i, j));
                    double dx1 = - 1 / (k1_x(mesh_x, i, j) * k1_x(mesh_x, i, j));
                    double dx2 =   0.0; 
                    double dx3 = - 1 / (k1_y(mesh_y, i, j) * k3_y(mesh_y, i, j));
                    double dx4 = - 1 / (k2_y(mesh_y, i, j) * k3_y(mesh_y, i, j)); 

                    double val[5] = {vol * dx0, vol * dx1, vol * dx2, vol * dx3, vol * dx4};

                    ll = {i, j};
                    HYPRE_StructMatrixSetValues(hypre_A, ll.val, 5, stencil_indices, val);
                }

                k = 2;
                if (get_node_type(i, j, stencil_offsets[k][0], stencil_offsets[k][1]) == 0) {

                    double vol = vmesh.val[i * vmesh.n2 + j];

                    double dx0 =   1 / (k2_x(mesh_x, i, j) * k2_x(mesh_x, i, j)) + 1/(k1_y(mesh_y, i, j) * k2_y(mesh_y, i, j));
                    double dx1 =   0.0;
                    double dx2 = - 1 / (k2_x(mesh_x, i, j) * k2_x(mesh_x, i, j));
                    double dx3 = - 1 / (k1_y(mesh_y, i, j) * k3_y(mesh_y, i, j));
                    double dx4 = - 1 / (k2_y(mesh_y, i, j) * k3_y(mesh_y, i, j)); 

                    double val[5] = {vol * dx0, vol * dx1, vol * dx2, vol * dx3, vol * dx4};

                    ll = {i, j};
                    HYPRE_StructMatrixSetValues(hypre_A, ll.val, 5, stencil_indices, val);
                }

                k = 3;
                if (get_node_type(i, j, stencil_offsets[k][0], stencil_offsets[k][1]) == 0) {

                    double vol = vmesh.val[i * vmesh.n2 + j];

                    double dx0 =   1 / (k1_x(mesh_x, i, j) * k2_x(mesh_x, i, j)) + 1/(k1_y(mesh_y, i, j) * k1_y(mesh_y, i, j));
                    double dx1 = - 1 / (k1_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
                    double dx2 = - 1 / (k2_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
                    double dx3 = - 1 / (k1_y(mesh_y, i, j) * k1_y(mesh_y, i, j));
                    double dx4 =   0.0;

                    double val[5] = {vol * dx0, vol * dx1, vol * dx2, vol * dx3, vol * dx4};

                    ll = {i, j};
                    HYPRE_StructMatrixSetValues(hypre_A, ll.val, 5, stencil_indices, val);
                }

                k = 4;
                if (get_node_type(i, j, stencil_offsets[k][0], stencil_offsets[k][1]) == 0) {

                    double vol = vmesh.val[i * vmesh.n2 + j];

                    double dx0 =   1 / (k1_x(mesh_x, i, j) * k2_x(mesh_x, i, j)) + 1/(k2_y(mesh_y, i, j) * k2_y(mesh_y, i, j));
                    double dx1 = - 1 / (k1_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
                    double dx2 = - 1 / (k2_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
                    double dx3 =   0.0;
                    double dx4 = - 1 / (k2_y(mesh_y, i, j) * k2_y(mesh_y, i, j));

                    double val[5] = {vol * dx0, vol * dx1, vol * dx2, vol * dx3, vol * dx4};

                    ll = {i, j};
                    HYPRE_StructMatrixSetValues(hypre_A, ll.val, 5, stencil_indices, val);
                }
            }
}

void rsolver::add_dirichlet_input(int i, int j, int stencil, int n_bc){
    
    dirichlet_input.val[n_input_dirichlet * 4 + 0] =  i;
    dirichlet_input.val[n_input_dirichlet * 4 + 1] =  j;
    dirichlet_input.val[n_input_dirichlet * 4 + 2] =  n_bc;

    switch(stencil){
        case 1: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = vmesh.val[i * vmesh.n2 + j] / (k1_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
            break;
        case 2: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = vmesh.val[i * vmesh.n2 + j] / (k2_x(mesh_x, i, j) * k3_x(mesh_x, i, j));
            break;
        case 3: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = vmesh.val[i * vmesh.n2 + j] / (k1_y(mesh_y, i, j) * k3_y(mesh_y, i, j));
            break;
        case 4: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = vmesh.val[i * vmesh.n2 + j] / (k2_y(mesh_y, i, j) * k3_y(mesh_y, i, j));
            break;
        default:
            cout << "stencil error in dirichlet bc" << endl;
            break;
    }
    
    n_input_dirichlet++;
}

int rsolver::dirichlet_boundary_number(int i, int j){
    for (int n = 0; n < n_dirichlet; n++)
    {
        if(in_box(i, j, dirichlet_boxes.val[n * 4 + 0], dirichlet_boxes.val[n * 4 + 1], dirichlet_boxes.val[n * 4 + 2], dirichlet_boxes.val[n * 4 + 3])){
            return n;
        }
    }

    return -1;   
}

// ------------------- solving --------------------------------

void rsolver::solve(fmatrix & solution, fmatrix & voltages, fmatrix & w_i, fmatrix & w_e){

    for(int i = 0; i < n_mesh_x; i++)
        for(int j = 0; j < n_mesh_y; j++){
            if(get_node_type(i, j) == 0 || get_node_type(i, j) == 2)
            {   
                ll = {i, j};
                HYPRE_StructVectorSetValues(hypre_b, ll.val, GAMMA * (w_i.val[i * w_i.n2 + j] - w_e.val[i * w_e.n2 + j]));
                HYPRE_StructVectorSetValues(hypre_x, ll.val, 0.0);
            } 
        }

    for (int n = 0; n < n_input_dirichlet; n++)
    {
        ll = {(int) dirichlet_input.val[n * 4 + 0], (int) dirichlet_input.val[n * 4 + 1]};
        HYPRE_StructVectorAddToValues(hypre_b, ll.val, voltages.val[(int) dirichlet_input.val[n * 4 + 2]] * dirichlet_input.val[n * 4 + 3]);
    }

    HYPRE_StructVectorAssemble(hypre_b);
    HYPRE_StructVectorAssemble(hypre_x);

    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &hypre_solver);
//    HYPRE_StructSMGSetTol(hypre_solver, 1.0e-6);
    HYPRE_StructSMGSetup(hypre_solver, hypre_A, hypre_b, hypre_x);
    HYPRE_StructSMGSolve(hypre_solver, hypre_A, hypre_b, hypre_x);

    double val[1];
    for(int i = 0; i < n_mesh_x; i++)
        for(int j = 0; j < n_mesh_y; j++){
            if(get_node_type(i, j) == 0 || get_node_type(i, j) == 2)
            {   
                ll = {i, j};
                HYPRE_StructVectorGetValues(hypre_x, ll.val, val);
                solution.val[i * solution.n2 + j] = val[0];
            }
        }
    for (int n = 0; n < n_dirichlet; n++)
    {
        solution.setbox_value(voltages.val[n], dirichlet_boxes.val[n * 4 + 0], dirichlet_boxes.val[n * 4 + 1], dirichlet_boxes.val[n * 4 + 2], dirichlet_boxes.val[n * 4 + 3]);
    }
    
    HYPRE_StructSMGDestroy(hypre_solver);
}
// ------------------- utilities ------------------------------

bool rsolver::in_box(int i, int j, int ill, int jll, int iur, int jur){
    return (i <= iur) && (i >= ill) && (j <= jur) && (j >= jll);
}





