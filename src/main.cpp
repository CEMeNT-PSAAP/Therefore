/*brief: assembly functions for matrices to compute therefore commands
date: May 23rd 2023
auth: J Piper Morgan (morgjack@oregonstate.edu)*/

#include <iostream>
#include <vector>
#include "util.h"
#include "builders.h"
//#include "H5Cpp.h"
#include "lapacke.h"
//#include <Eigen/Dense>
//#include <cusparse_v2.h>
//#include <cuda.h>

/* compile notes and prerecs

    In ubuntu yyou need these libarires:
        sudo apt-get install libblas-dev checkinstall
        sudo apt-get install libblas-doc checkinstall
        sudo apt-get install liblapacke-dev checkinstall
        sudo apt-get install liblapack-doc checkinstall

    Should be able to configure with: 
        g++ main.cpp -std=c++20 -llapacke
*/

void eosPrint(ts_solutions state);

// row major!!!!!!
//extern "C" void dgesv_( int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *lbd, int *info  );
//extern "C" void LAPACKE_dgesv_( LAPACK_ROW_MAJOR, int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *lbd, int *info  );
//-llapacke
std::vector<double> row2colSq(std::vector<double> row);

// i space, m is angle, k is time, g is energy group

const bool print_mats = false;
const bool cycle_print = true;

int main(void){

    print_title();

    using namespace std;
    
    // problem definition
    // eventually from an input deck
    double dx = 0.1;
    double dt = 0.5;
    vector<double> v = {1, .5};
    vector<double> xsec_total = {1, 0.5};
    vector<double> xsec_scatter = {0.25, 0.1};
    vector<double> Q = {1, 0};
    double Length = 1;
    double IC_homo = 0;
    
    int N_cells = 10; 
    int N_angles = 4; 
    int N_time = 1;
    int N_groups = 2;

    // 4 = N_subspace (2) * N_subtime (2)
    int N_mat = 4 * N_cells * N_angles * N_groups;

    // N_cm is the size of the row major vector
    int N_rm = N_mat*N_mat;

    // homogeneous initial condition vector
    // will be stored as the solution at time=0
    vector<double> IC(N_mat, 0.0);
    for (int p=0; p<N_cells*2; p++){IC[p] = IC[p]*IC_homo;}

    // actual computation below here

    // generate g-l quadrature angles and weights
    vector<double> weights(N_angles, 0.0);
    vector<double> angles(N_angles, 0.0);

    quadrature(angles, weights);


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
    ps.initialize_from_previous = true;
    ps.max_iteration = int(1e4);

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
    int ldb = 1; // 
    std::vector<int> i_piv(N_mat, 0);  // pivot column vector
    int info;

    // generation of the whole ass mat
    A_gen(A, cells, ps);

    if (print_mats){
        print_rm(A);
    }

    print_rm(A);

    vector<double> b(N_mat);

    cout << "howdy" << endl;

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
        bool converged = true;  // converged boolean
        double spec_rad;

        //vector<double> A_copy;
        vector<double> A_copy(N_mat);

        while (converged){

            // lapack requires a copy of data that it uses for row piviot (A after _dgesv != A)
            A_copy = A;

            b_gen(b, aflux_previous, aflux_last, cells, ps);
            // reminder: last refers to iteration, previous refers to time step
            
            if (print_mats){
                cout << "Cycle: " << itter << endl;
                cout << "RHS" << endl;
                print_vec_sd(b);
            }
            // solve Ax=b
            info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, N_mat, nrhs, &A_copy[0], lda, &i_piv[0], &b[0], ldb );

            //info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );
            //Lapack solver dgesv_( &N_mat, &nrhs, &*A.begin(), &lda, &*i_piv.begin(), &*b.begin(), &ldb, &info );

            if (print_mats){
                cout << "x" <<endl;
                print_vec_sd(b);
            }

            if( info > 0 ) {
                printf( "The diagonal element of the triangular factor of A,\n" );
                printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
                printf( "the solution could not be computed.\n" );
                exit( 1 );
            }

            // compute the relative error between the last and current iteration
            error = infNorm_error(aflux_last, b);

            // compute spectral radius
            spec_rad = abs(error-error_n1) / abs(error_n1 - error_n2);

            // too allow for a error computation we need at least three cycles
            if (itter > 3){
                // if relative error between the last and just down iteration end the time step
                if ( error < ps.convergence_tolerance ){ converged = false; }
            }

            if (itter > ps.max_iteration){
                cout << "WARNING: Computation did not converge after " << ps.max_iteration << "iterations" << endl;
                cout << "       itter: " << itter << endl;
                cout << "       error: " << error << endl;
                cout << "" << endl;
                converged = false;
            }

            aflux_last = b;
            itter++;

            if (cycle_print){
                cout << "cycle: " <<itter<<" of time step: "<<t<<endl;
                cout << "error: " <<error<<" error-1: "<<error_n1<<" error-2: "<<error_n2<<endl; 
                cout << "spec rad: " << spec_rad <<endl;
                cout << "" << endl;
            }

            error_n2 = error_n1;
            error_n1 = error;

        } // end convergence loop

        string ext = ".csv";
        string file_name = "afluxUnsorted";
        string dt = to_string(t);

        file_name = file_name + dt + ext;

        std::ofstream output(file_name);
        output << "TIME STEP: " << t << "Unsorted solution vector" << endl;
        output << "N_space: " << ps.N_cells << " N_groups: " << ps.N_groups << " N_angles: " << ps.N_angles << endl;
        for (int i=0; i<b.size(); i++){
            output << b[i] << "," << endl;
        }

        cout << "file saved under: " << file_name << endl;

        // store solution vector org
        //ts_solutions save_timestep;
        //save_timestep.time = (t+1)*ps.dt;
        //save_timestep.spectral_radius = spec_rad;
        //save_timestep.N_step = t+1;
        //save_timestep.number_iteration = itter;

        // print end of step information
        //eosPrint(save_timestep);
        //print_vec_sd(b);
        //save_timestep.aflux = b;
        //print_vec_sd(save_timestep.aflux);
        //solutions.push_back(save_timestep);
        print_vec_sd(b);

    } // end of time step loop

    //WholeProblem wp = WholeProblem(cells, ps, solutions);
    //wp.PublishUnsorted();
    
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