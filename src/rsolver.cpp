
#include "rsolver.h"

#include <chrono>
#include <iostream>

#include "mpi.h"
#include "HYPRE_struct_ls.h"
#include "fmatrix.h"
#include "fields.h"
#include "configuration.h"

#define SOLVER_PRINT_LEVEL 0
#define SOLVER_TOLERANCE 1e-6
#define SOLVER_F(A, B) A ## SMG ## B // Solvers: SMG, PCG

using namespace std::chrono;

using namespace std;
// ------------------- constructors / destructor -------------------------

rsolver::rsolver(mesh_set * _mesh, int _n_neumann, int _n_dirichlet, configuration * _config): mesh(_mesh), config(_config){

    n_neumann = _n_neumann;
    n_dirichlet = _n_dirichlet;
    n_mesh_x = (int) _mesh->nx;
    n_mesh_y = (int) _mesh->ny;
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

rsolver::rsolver(mesh_set * _mesh, configuration * _config): mesh(_mesh), config(_config){

    n_mesh_x = (int) _mesh->nx;
    n_mesh_y = (int) _mesh->ny;
    n_input_dirichlet = 0;
    n_solve = 0;

    ll = ur = inner_ll = inner_ur = imatrix::zeros(2);
    node_type = imatrix::zeros(n_mesh_x, n_mesh_y);
    dirichlet_input = fmatrix::zeros(2 * n_mesh_x + 2 * n_mesh_y, 4);

    inner_ll = {1, 1};
    inner_ur = {n_mesh_x - 2, n_mesh_y - 2};
}

void rsolver::late_init(int _n_neumann, int _n_dirichlet){
    n_neumann = _n_neumann;
    n_dirichlet = _n_dirichlet;

    dirichlet_boxes = imatrix::zeros(n_dirichlet, 4);
    neumann_boxes = imatrix::zeros(n_neumann, 4);
}

rsolver::~rsolver(){
    HYPRE_StructGridDestroy(hypre_grid);
    HYPRE_StructStencilDestroy(hypre_stencil);
    HYPRE_StructMatrixDestroy(hypre_A);
    HYPRE_StructVectorDestroy(hypre_b);
    HYPRE_StructVectorDestroy(hypre_x);
}


// ------------------- initialization ----------------------------

void rsolver::init_grid(){
    HYPRE_StructGridCreate(MPI_COMM_SELF, 2, &hypre_grid);
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
    HYPRE_StructMatrixCreate(MPI_COMM_SELF, hypre_grid, hypre_stencil, &hypre_A);
    HYPRE_StructMatrixInitialize(hypre_A);
}

void rsolver::init_vectors(){
    HYPRE_StructVectorCreate(MPI_COMM_SELF, hypre_grid, &hypre_b);
    HYPRE_StructVectorCreate(MPI_COMM_SELF, hypre_grid, &hypre_x);
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

                double vol = mesh->v.val[i * mesh->v.n2 + j];
                double d0 =   1/ (mesh->k1_x(i, j) * mesh->k2_x(i, j)) + 1/(mesh->k1_y(i, j) * mesh->k2_y(i, j));
                double d1 = - 1/ (mesh->k1_x(i, j) * mesh->k3_x(i, j));
                double d2 = - 1/ (mesh->k2_x(i, j) * mesh->k3_x(i, j));
                double d3 = - 1/ (mesh->k1_y(i, j) * mesh->k3_y(i, j));
                double d4 = - 1/ (mesh->k2_y(i, j) * mesh->k3_y(i, j)); 
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
    for (int i = 0; i < n_mesh_x; i++)
        for (int j = 0; j < n_mesh_y; j++)
            
            if(get_node_type(i, j) == 2){
                
                double d0 = 0;
                double d1 = 0;
                double d2 = 0;
                double d3 = 0;
                double d4 = 0;

                // Neumann point in the x direction
                if (get_node_type(i, j, stencil_offsets[2][0], stencil_offsets[2][1]) == -1) {

                    d0 +=  1 / (mesh->k1_x(i, j) * mesh->k1_x(i, j));
                    d1 = - 1 / (mesh->k1_x(i, j) * mesh->k1_x(i, j));
                    d2 =   0.0;

                }
                else if (get_node_type(i, j, stencil_offsets[1][0], stencil_offsets[1][1]) == -1) {

                    d0 +=  1 / (mesh->k2_x(i, j) * mesh->k2_x(i, j));
                    d1 =   0.0;
                    d2 = - 1 / (mesh->k2_x(i, j) * mesh->k2_x(i, j));

                }
                else {
                    d0 +=  1 / (mesh->k1_x(i, j) * mesh->k2_x(i, j));
                    d1 = - 1 / (mesh->k1_x(i, j) * mesh->k3_x(i, j));
                    d2 = - 1 / (mesh->k2_x(i, j) * mesh->k3_x(i, j));
                }
                
                // Neumann point in the y direction
                if (get_node_type(i, j, stencil_offsets[4][0], stencil_offsets[4][1]) == -1) {

                    d0 +=  1 / (mesh->k1_y(i, j) * mesh->k1_y(i, j));
                    d3 = - 1 / (mesh->k1_y(i, j) * mesh->k1_y(i, j));
                    d4 =   0.0;

                }
                else if (get_node_type(i, j, stencil_offsets[3][0], stencil_offsets[3][1]) == -1) {

                    d0 +=  1 / (mesh->k2_y(i, j) * mesh->k2_y(i, j));
                    d3 =   0.0;
                    d4 = - 1 / (mesh->k2_y(i, j) * mesh->k2_y(i, j));
                    
                }
                else {
                    d0 +=  1 / (mesh->k1_y(i, j) * mesh->k2_y(i, j));
                    d3 = - 1 / (mesh->k1_y(i, j) * mesh->k3_y(i, j));
                    d4 = - 1 / (mesh->k2_y(i, j) * mesh->k3_y(i, j));
                }
                
                double vol = mesh->v.val[i * mesh->v.n2 + j];
                double val[5] = {vol * d0, vol * d1, vol * d2, vol * d3, vol * d4};
                
                ll = {i, j};
                HYPRE_StructMatrixSetValues(hypre_A, ll.val, 5, stencil_indices, val);
            }
}

void rsolver::add_dirichlet_input(int i, int j, int stencil, int n_bc){
    
    dirichlet_input.val[n_input_dirichlet * 4 + 0] =  i;
    dirichlet_input.val[n_input_dirichlet * 4 + 1] =  j;
    dirichlet_input.val[n_input_dirichlet * 4 + 2] =  n_bc;

    switch(stencil){
        case 1: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = mesh->v.val[i * mesh->v.n2 + j] / (mesh->k1_x(i, j) * mesh->k3_x(i, j));
            break;
        case 2: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = mesh->v.val[i * mesh->v.n2 + j] / (mesh->k2_x(i, j) * mesh->k3_x(i, j));
            break;
        case 3: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = mesh->v.val[i * mesh->v.n2 + j] / (mesh->k1_y(i, j) * mesh->k3_y(i, j));
            break;
        case 4: 
            dirichlet_input.val[n_input_dirichlet * 4 + 3] = mesh->v.val[i * mesh->v.n2 + j] / (mesh->k2_y(i, j) * mesh->k3_y(i, j));
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
    double gamma = config->f("p/gamma");
    for(int i = 0; i < n_mesh_x; i++)
        for(int j = 0; j < n_mesh_y; j++){
            if(get_node_type(i, j) == 0 || get_node_type(i, j) == 2)
            {   
                ll = {i, j};
                HYPRE_StructVectorSetValues(hypre_b, ll.val, gamma * (w_i.val[i * w_i.n2 + j] - w_e.val[i * w_e.n2 + j]));
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

    SOLVER_F(HYPRE_Struct, Create)(MPI_COMM_SELF, &hypre_solver);
    SOLVER_F(HYPRE_Struct, SetTol)(hypre_solver, SOLVER_TOLERANCE);
    // SOLVER_F(HYPRE_Struct, SetPrintLevel)(hypre_solver, SOLVER_PRINT_LEVEL);
    SOLVER_F(HYPRE_Struct, Setup)(hypre_solver, hypre_A, hypre_b, hypre_x);
    SOLVER_F(HYPRE_Struct, Solve)(hypre_solver, hypre_A, hypre_b, hypre_x);

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
    
    SOLVER_F(HYPRE_Struct, Destroy)(hypre_solver);
}
// ------------------- utilities ------------------------------

bool rsolver::in_box(int i, int j, int ill, int jll, int iur, int jur){
    return (i <= iur) && (i >= ill) && (j <= jur) && (j >= jll);
}


void rsolver::setup(mesh_set & mesh, imatrix & electrode_mask){
    
    int n_dirichlet = 0, n_neumann = 0;

    if(config->s("boundaries/ob_type") == "dirichlet"){
        n_dirichlet = 3;
        n_neumann   = 2;
    }
    if(config->s("boundaries/ob_type") == "neumann"){
        n_dirichlet = 1;
        n_neumann   = 4;
    }

    late_init(n_neumann, n_dirichlet);
    
    imatrix box_thruster        = {0, 0, 0, config->i("geometry/n_thruster") - 1};
    imatrix box_top_thruster    = {0, config->i("geometry/n_thruster"), 0, config->i("geometry/n_mesh_y") - 2};
    imatrix box_ob_top          = {0, config->i("geometry/n_mesh_y") - 1, config->i("geometry/n_mesh_x") - 2, config->i("geometry/n_mesh_y") - 1};
    imatrix box_ob_right        = {config->i("geometry/n_mesh_x") - 1, 0, config->i("geometry/n_mesh_x") - 1, config->i("geometry/n_mesh_y") - 1};
    imatrix box_sym             = {1, 0, config->i("geometry/n_mesh_x") - 2, 0};

    if(config->s("boundaries/ob_type") == "neumann"){
        set_dirichlet_box(box_thruster, 0);
        set_neumann_box(box_top_thruster, 0);
        set_neumann_box(box_ob_top, 1);
        set_neumann_box(box_ob_right, 2);
        set_neumann_box(box_sym, 3);
    }

    if(config->s("boundaries/ob_type") == "dirichlet"){
        
        set_neumann_box(box_sym, 0);
        set_neumann_box(box_top_thruster, 1);

        set_dirichlet_box(box_ob_top, 0);
        set_dirichlet_box(box_thruster, 1);
        set_dirichlet_box(box_ob_right, 2);
            
        electrode_mask.setbox_value(1, box_ob_top.val[0], box_ob_top.val[1], box_ob_top.val[2], box_ob_top.val[3]);
        electrode_mask.setbox_value(1, box_ob_right.val[0], box_ob_right.val[1], box_ob_right.val[2], box_ob_right.val[3]);
    }

    electrode_mask.setbox_value(2, box_thruster.val[0], box_thruster.val[1], box_thruster.val[2], box_thruster.val[3]);

    assemble();

}   


void rsolver::setup_benchmark(mesh_set & mesh, imatrix & electrode_mask){
    
    int n_dirichlet = 2;
    int n_neumann   = 2;

    late_init(n_neumann, n_dirichlet);
    
    imatrix box_left_electrode  = {0, 0, 0, config->i("geometry/n_mesh_y") - 1};
    imatrix box_right_electrode = {config->i("geometry/n_mesh_x") - 1, 0, config->i("geometry/n_mesh_x") - 1, config->i("geometry/n_mesh_y") - 1};
    imatrix box_top_sym         = {1, config->i("geometry/n_mesh_y") - 1, config->i("geometry/n_mesh_x") - 2, config->i("geometry/n_mesh_y") - 1};
    imatrix box_bot_sym         = {1, 0, config->i("geometry/n_mesh_x") - 2, 0};

    set_neumann_box(box_top_sym, 0);
    set_neumann_box(box_bot_sym, 1);
    set_dirichlet_box(box_left_electrode, 0);
    set_dirichlet_box(box_right_electrode, 1);

    electrode_mask.setbox_value(1, box_right_electrode.val[0], box_right_electrode.val[1], box_right_electrode.val[2], box_right_electrode.val[3]);
    electrode_mask.setbox_value(2, box_left_electrode.val[0], box_left_electrode.val[1], box_left_electrode.val[2], box_left_electrode.val[3]);

    assemble();
}   


