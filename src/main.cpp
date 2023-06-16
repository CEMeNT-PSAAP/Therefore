/*brief: assembly functions for matrices to compute therefore commands
date: May 23rd 2023
auth: J Piper Morgan (morgjack@oregonstate.edu)*/

#include <iostream>
#include <vector>
#include "legendre.h"
#include "util.h"
#include "builders.h"
#include "H5Cpp.h"
#include "mkl_lapacke.h"
//#include <Eigen/Dense>
//#include <cusparse_v2.h>
//#include <cuda.h>

void eosPrint(ts_solutions state);

// row major!!!!!!
//extern "C" void dgesv_( int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *lbd, int *info  );
//extern "C" void LAPACKE_dgesv_( LAPACK_ROW_MAJOR, int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *lbd, int *info  );
//-llapacke
std::vector<double> row2colSq(std::vector<double> row);

// i space, m is angle, k is time, g is energy group

int main(void){

    print_title();

    using namespace std;
    
    // problem definition
    // eventually from an input deck
    double dx = 0.5;
    double dt = 0.5;
    vector<double> v = {1};
    vector<double> xsec_total = {1};
    vector<double> xsec_scatter = {0};
    vector<double> Q = {1};
    double Length = 1;
    double IC_homo = 0;
    
    int N_cells = 2; 
    int N_angles = 2; 
    int N_time = 1;
    int N_groups = 1;

    // 4 = N_subspace (2) * N_subtime (2)
    int N_mat = 4 * N_cells * N_angles * N_groups;

    // N_cm is the size of the row major vector
    int N_rm = N_mat*N_mat;

    // homogeneous initial condition vector
    // will be stored as the solution at time=0
    vector<double> IC(N_mat, 1.0);
    for (int p=0; p<N_cells*2; p++){IC[p] = IC[p]*IC_homo;}

    // actual computation below here

    // generate g-l quadrature angles and weights
    double weights[N_angles];
    double angles[N_angles];
    legendre_compute_glr(N_angles, angles, weights);

    // problem space class construction
    problem_space ps;
    ps.dt = dt;
    ps.dx = dx;
    ps.N_angles = N_angles;
    ps.N_cells = N_cells;
    ps.N_groups = N_groups;
    ps.N_time = N_time;
    ps.angles = angles;
    ps.weights = weights;
    ps.convergence_tolerance = 1e-9;
    ps.initialize_from_previous = true;
    ps.max_iteration = int(1);

    // cell construction;
    vector<cell> cells;

    for (int i=0; i<N_cells; i++){
        cell cellCon;
        cellCon.cell_id = i;
        cellCon.x_left = i*dx;
        cellCon.xsec_scatter = xsec_scatter;
        cellCon.xsec_total = xsec_total;
        cellCon.dx = dx;
        cellCon.v = v;
        cellCon.dt = dt;
        cellCon.Q = Q;

        cells.push_back(cellCon);
    }
    
    // initial condition stored as first element of solution vector
    vector<ts_solutions> solutions;

    ts_solutions sol_con;
    sol_con.aflux = IC;
    sol_con.time = 0.0;
    sol_con.spectral_radius = 0.0;
    sol_con.N_step = 0;
    sol_con.number_iteration = 0;

    solutions.push_back(sol_con);

    // allocation of the whole ass mat
    vector<double> A(N_rm);
    
    // vector org angular flux from last iteration
    vector<double> aflux_last(N_mat, 0.0);
    // vector org converged angular flux from previous time step
    vector<double> aflux_previous(N_mat, 0.0);
    // initializing the inital previous as the IC
    aflux_previous = IC;

    int nrhs = 1; // one column in b
    int lda = N_mat;
    int ldb = N_mat;
    std::vector<int> i_piv(N_mat, 0);  // pivot column vector
    int info;


    // generation of the whole ass mat
    A_gen(A, cells, ps);

    print_rm(A);

    //A = row2colSq(A);

    //print_vec_sd(A);
    //print_rm(A);

    //cout <<"thru" << endl;

    //print_cm(A);

    vector<double> b(N_mat);


    // time step loop
    for(int t=0; t<N_time; ++t){

        
        if (ps.initialize_from_previous){
            // all the angular fluxes start from the previous converged time step
            aflux_last = solutions[t].aflux;
        } else {
            // all angular fluxes start this time step iteration from 0
            fill(aflux_last.begin(), aflux_last.end(), 0.0);
        }

        vector<double> b(N_mat, 0.0);
        
        int itter = 0;          // iteration counter
        double error = 1;       // error from current iteration
        double error_n1 = 1;    // error back one iteration (from last)
        double error_n2 = 1;    // error back two iterations
        bool converged = true; // converged boolean
        double spec_rad;

        while (converged){
            // build b
            b_gen(b, aflux_previous, aflux_last, cells, ps);
            // reminder: last refers to iteration, previous refers to time step

            print_vec_sd(b);
            // solve Ax=b
            info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, N_mat, nrhs, &A[0], lda, &i_piv[0], &b[0], ldb );
            //Lapack solver dgesv_( &N_mat, &nrhs, &*A.begin(), &lda, &*i_piv.begin(), &*b.begin(), &ldb, &info );
            //print_vec_sd(b);

            print_rm(A);
            print_vec_sd(b);

            // compute the relative error between the last and current iteration
            error = infNorm_error(aflux_last, b);

            // move errors back
            error_n1 = error;
            error_n2 = error_n1;

            // compute spectral radius
            spec_rad = abs(error-error_n1) / abs(error_n1 - error_n2);

            // too allow for a error computation we need at least three cycles
            if (itter > 3){
                // if relative error between the last and just down iteration end the time step
                if ( error < ps.convergence_tolerance ){ converged = false; }
            }

            if (itter > ps.max_iteration){
                cout << "WARNING: Computation did not converge" << endl;
                cout << "       itter: " << itter << endl;
                cout << "       error: " << error << endl;
                cout << "" << endl;
                converged = false;
            }

        itter++;

        } // end convergence loop

        // store solution vector org
        ts_solutions save_timestep;
        save_timestep.time = (t+1)*ps.dt;
        save_timestep.spectral_radius = spec_rad;
        save_timestep.N_step = t+1;
        save_timestep.number_iteration = itter;

        // print end of step information
        eosPrint(save_timestep);
        //print_vec_sd(b);
        save_timestep.aflux = b;
        //print_vec_sd(save_timestep.aflux);
        solutions.push_back(save_timestep);


    } // end of time step loop

    WholeProblem wp = WholeProblem(cells, ps, solutions);
    wp.PublishUnsorted();
    
    return(0);
} // end of main



std::vector<double> row2colSq(std::vector<double> row){
    /*brief */

    int SIZE = sqrt(row.size());

    std::vector<double> col(SIZE);

    for (int i = 0; i < SIZE; ++i){
        for (int j = 0; j < SIZE; ++j){
            col[ i * SIZE + j ] = row[ j * SIZE + i ];
        }
    }

    return(col);
}

void std2cuda_bsr(problem_space ps){
    int block_size = 4*ps.N_angles*ps.N_groups;
    int number_row_blocks = ps.N_cells;
    int number_col_blocks = ps.N_cells;
    int number_nonzero_blocks = ps.N_cells;

}
/*
void A_gen_c_g(){
    
    breif: assmebles a coeeficant matrix within a given group and cell for all angles
    NOTE: ROW MAJOR FORMAT
    
}
*/

void eosPrint(ts_solutions state){
    using namespace std;
    /* brief: end of time step printer*/

    if (state.N_step == 0){
        // print header
        printf("Time    Step Table\n");
        printf("step    time    error   spectral\n");
        printf("                        radius\n");
        printf("============================================================================\n");
    } 
    
    // print cycle information table
    printf("%3d     %2.3e    %1.4e     %5f\n", state.N_step, state.time, state.final_error, state.spectral_radius);
}