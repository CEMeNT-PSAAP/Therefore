/*brief: assembly functions for matricies to compute therefore commands
date: May 23rd 2023
auth: J Piper Morgan (morgjack@oregonstate.edu)*/

#include <iostream>
#include <vector>
#include "legendre.h"
#include "util.h"
#include "base_mats.h"
#include <Eigen/Dense>
//#include <cusparse_v2.h>
//#include <cuda.h>

// row major!!!!!!
void b_gen(std::vector<double> &b, std::vector<double> aflux_previous, std::vector<double> aflux_last, std::vector<cell> cells, problem_space ps);
void A_c_gen(int i, std::vector<double> &A_c, std::vector<cell> cells, problem_space ps);
void A_gen(std::vector<double> &A, std::vector<cell> cells, problem_space ps);
Eigen::VectorXd std2eigen(std::vector<double> sv);



// i space, m is angle, k is time, g is energy group

int main(void){

    print_title();

    using namespace std;
    
    // problem deffinition
    // eventually from an input deck
    double dx = 0.5;
    double dt = 0.5;
    vector<double> v = {1, 0.5};
    vector<double> xsec_total = {2, 1};
    vector<double> xsec_scatter = {1, 0.5};
    vector<double> Q = {0, 1};
    double Length = 1;
    double IC_homo = 0;
    
    int N_cells = 2; 
    int N_angles = 2; 
    int N_time = 5;
    int N_groups = 1;

    // homogenious initil condition vector
    // will be stored as the soultion at time=0
    vector<double> IC(N_cells*2, 1.0);
    for (int p=0; p<N_cells*2; p++){IC[p] = IC[p]*IC_homo;}

    // actaul computation below here

    // generate g-l quadrature angles and weigths
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
    
    // initial condition stored as first element of soultion vector
    vector<ts_soultion> soultions;

    ts_soultion sol_con;
    sol_con.aflux = IC;
    sol_con.aflux_h = IC;
    sol_con.time = 0.0;
    sol_con.specral_radius = 0.0;
    sol_con.N_step = 0;
    sol_con.number_itteration = 0;

    soultions.push_back(sol_con);


    // 4 = N_subspace (2) * N_subtime (2)
    int N_mat = 4 * N_cells * N_angles * N_groups;

    // N_cm is the size of the row major vetor
    int N_rm = N_mat*N_mat;

    // allocation of the whole ass mat
    vector<double> A(N_rm);
    
    // generation of the whole ass mat
    A_gen(A, cells, ps);

    print_rm(A);

    bool converged = false;
    int itter = 0;
    double error = 1;

    // vector org angular flux from last itteration
    vector<double> aflux_last(N_mat);
    // vector org converged angular flux from previous time step
    vector<double> aflux_previous(N_mat);
    // initilizing the inital previous as the IC
    aflux_previous = IC;

    vector<double> b(N_mat);

    b_gen(b, aflux_previous, aflux_last, cells, ps);

    print_vec_sd(b);

    cout << "" << endl;
    Eigen::VectorXd b_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());

    //Eigen::Vector3d b_eig = std2eigen(b);

    cout << b_eig << endl;


    Eigen::MatrixXd A_eig(N_mat, N_mat, Eigen::RowMajor); //A_eig = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>(A.data(), A.size());
    //A_eig[0,0] = A[0];
    cout << A_eig <<endl;

    /*
    for(int t=0; t<N_time; ++t){

        aflux_last = soultions[t].aflux;
        //aflux_last_h = soultions[t].aflux_h;

        vector<double> b(N_mat, 0.0);
        
        while (converged){

            //build b
            b_gen(b, aflux_previous, aflux_last, cells, ps);

            //solve Ax=b

            //compute error

            if (itter > 3){
                if (error < ps.convergence_tolarance){
                    converged == true;
                }
            }

        itter++;
        }

        // store soultion vector org
    }*/

    // convert vetor org to more friendly format

    // save state

    return(0);
}


void b_gen(std::vector<double> &b, std::vector<double> aflux_previous, std::vector<double> aflux_last, std::vector<cell> cells, problem_space ps){
    //breif: builds b
    vector<double> b_small;

    // size of the cell blocks in all groups and angle
    int size_cellBlocks = ps.N_angles*ps.N_groups*4;
    // size of the group blocks in all angle within a cell
    int size_groupBlocks = ps.N_angles*4;
    // size of the angle blocks within a group and angle
    int size_angleBlocks = 4;
    // helper index
    int index_start;

    for (int i=0; i<ps.N_cells; i++){
        

        for (int g=0; g<ps.N_groups; g++){

            // angular fluxes from the right bound (lhs of cell at right) last itteration
            double af_rb;
            // angular fluxes from the left bound (rhs if cell at left) last itteration
            double af_lb;
            // angular fluxes from the right bound (lhs of cell at right) k+1/2 last itteration
            double af_hn_rb;
            // angular fluxes from the left bound (rhs of cell at left) k+1/2 last itteration
            double af_hn_lb;
            // k is time step index
            

            for (int j=0; j<ps.N_angles; j++){
                // the first index in the smallest chunk of 4
                index_start = (i*(size_cellBlocks) + g*(size_groupBlocks) + 4*j);
                // 4 blocks orginized af_l, af_r, af_hn_l, af_hn_r

                // angular flux from the k-1+1/2 from within the cell
                double af_hl_l = aflux_previous[index_start+2];
                double af_hl_r = aflux_previous[index_start+3];

                // negative angle
                if (ps.angles[g] < 0){
                    if (i == ps.N_cells-1){ // right boundary condition
                        af_rb = ps.boundary_condition();
                        af_hn_rb = ps.boundary_condition();
                    } else { // pulling information from right to left
                        af_rb = aflux_last[index_start+1+4];
                        af_hn_rb = aflux_last[index_start+3+4];
                    }

                    b_small = b_neg(cells[i], g, ps.angles[j], af_hl_l, af_hl_r, af_rb, af_hn_rb);

                // positive angles
                } else {
                    if (i == 0){ // left boundary condition
                        af_lb    = ps.boundary_condition();
                        af_hn_lb = ps.boundary_condition();
                    } else { // pulling information from left to right
                        af_lb = aflux_last[index_start-2];
                        af_hn_lb = aflux_last[index_start+2-2];
                    }

                    b_small = b_pos(cells[i], g, ps.angles[j], af_hl_l, af_hl_r, af_rb, af_hn_rb);
                }
                b[index_start] = b_small[0];
                b[index_start+1] = b_small[1];
                b[index_start+2] = b_small[2];
                b[index_start+3] = b_small[3];
            }
        }
    }
}


void A_gen(std::vector<double> &A, std::vector<cell> cells, problem_space ps){ 

    int dimA_c = 4 * ps.N_groups * ps.N_angles;

    for (int i=0; i<ps.N_cells; i++){
        
        vector<double> A_c(dimA_c*dimA_c, 0.0);

        A_c_gen(i, A_c, cells, ps);

        int A_id_start = dimA_c*ps.N_cells * dimA_c*i + dimA_c*i;

        for (int r=0; r<dimA_c; r++){
            for (int c=0; c<dimA_c; c++){
                int id_wp = A_id_start + r * (dimA_c*ps.N_cells) + c ;
                int id_c = dimA_c*r + c;
                A[id_wp] = A_c[id_c];
            }
        }
    }
}


void A_c_gen(int i, std::vector<double> &A_c, std::vector<cell> cells, problem_space ps){
    /*
    breif: assembles a coefficant matrix within all groups and angles in a cell
    NOTE: ROW MAJOR FORMAT
    */

   for (int g=0; g<ps.N_groups; g++){

        vector<double> A_c_g(4*ps.N_angles * 4*ps.N_angles);
        vector<double> A_c_g_a(4*4);
        vector<double> S(4*ps.N_angles * 4*ps.N_angles);

        for (int j=0; j<ps.N_angles; j++){
            if (ps.angles[j] > 0){
                A_c_g_a = A_pos_rm(cells[i], ps.angles[j], g);
            } else {
                A_c_g_a = A_neg_rm(cells[i], ps.angles[j], g);
            }

            // push it into an all angle cellwise fuck me
            for (int r=0; r<4; r++){
                for (int c=0; c<4; c++){
                    // awfully confusing I know
                    // id_acell = (moves us betten angle blocks) + (moves us to diagonal) + 
                    //            (moves us in rows w/in an angle) + (moves us in col w/in an angle)
                    int id_acell  = ((4*ps.N_angles * 4*j) + (4*j) + (4*ps.N_angles)*r + (c));
                    int id_ancell = 4*r + c;
                    A_c_g[id_acell] = A_c_g_a[id_ancell];
                }
            }
        }

        S = scatter(cells[i], ps.weights, ps.N_angles, g);

        int index_start = 4*g*ps.N_angles * 4*ps.N_angles*ps.N_groups + 4*g*ps.N_angles;
        int Adim_angle = 4*ps.N_angles; 

        for (int r=0; r<Adim_angle; r++){
            for (int c=0; c<Adim_angle; c++){

                int id_group = Adim_angle*r + c;
                int id_c_g = index_start + r*(Adim_angle*ps.N_groups) + c;

                A_c[id_c_g] = A_c_g[id_group] - S[id_group];
            }
        }
    
    }
}


Eigen::VectorXd std2eigen(std::vector<double> sv){
    Eigen::VectorXd ev = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(sv.data(), sv.size());
    return(ev);
}


void std2cuda_bsr(problem_space ps){
    int block_size = 4*ps.N_angles*ps.N_groups;
    int number_row_blocks = ps.N_cells;
    int number_col_blocks = ps.N_cells;
    int number_nonzero_blocks = ps.N_cells;

}
/*
void A_gen_c_g(){
    /*
    breif: assmebles a coeeficant matrix within a given group and cell for all angles
    NOTE: ROW MAJOR FORMAT
    
}
*/