/*brief: assembly functions for matricies to compute therefore commands
date: May 23rd 2023
auth: J Piper Morgan (morgjack@oregonstate.edu)*/

#include <iostream>
#include <vector>
#include "legendre.h"
#include "util.h"
#include "base_mats.h"
//#include <cusparse_v2.h>
//#include <cuda.h>

// row major!!!!!!



// i space, m is angle, k is time, g is energy group

int main(void){

    print_title();

    using namespace std;
    
    double dx = 1;
    double dt = 0.5;
    vector<double> v = {1, 0.5};
    vector<double> xsec_total = {2, 1};
    vector<double> xsec_scatter = {1, 0.5};
    double Length = 1;

    int N_cells = 1; 
    int N_angles = 2; 
    int N_time = 5;
    int N_groups = 2;

    double weights[N_angles];
    double angles[N_angles];
    legendre_compute_glr(N_angles, angles, weights);

    //cell construction;
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

        cells.push_back(cellCon);
    }

    // 4 = N_subspace (2) * N_subtime (2)
    int N_mat = 4 * N_cells * N_angles * N_groups;

    // N_cm is the size of the row major vetor
    int N_rm = N_mat*N_mat;

    vector<double> A(N_rm);

    for (int i=0; i<N_cells; i++){

        vector<double> A_group(4 * N_groups * N_angles * 4 * N_groups * N_angles);

        for (int g=0; g<N_groups; g++){

            vector<double> A_cell_group(4*N_angles * 4*N_angles);
            vector<double> A_an_cell(4*4);
            vector<double> S(4*N_angles * 4*N_angles);

            for (int j=0; j<N_angles; j++){
                if (angles[j] > 0){
                    A_an_cell = A_pos_rm(cells[i], angles[j], g);
                } else {
                    A_an_cell = A_neg_rm(cells[i], angles[j], g);
                }

                // push it into an all angle cellwise fuck me
                for (int r=0; r<4; r++){
                    for (int c=0; c<4; c++){
                        // awfully confusing I know
                        // id_acell = (moves us betten angle blocks) + (moves us to diagonal) + 
                        //            (moves us in rows w/in an angle) + (moves us in col w/in an angle)
                        int id_acell  = ((4*N_angles * 4*j) + (4*j) + (4*N_angles)*r + (c));
                        int id_ancell = 4*r + c;
                        A_cell_group[id_acell] = A_an_cell[id_ancell];
                    }
                }
            }

            S = scatter(cells[i], weights, N_angles, g);


            int index_start = 4*g*N_angles * 4*N_angles*N_groups + 4*g*N_angles;
            int Adim_angle = 4*N_angles; 

            for (int r=0; r<Adim_angle; r++){
                for (int c=0; c<Adim_angle; c++){

                    int id_group = Adim_angle*r + c;
                    int id_cell_group = index_start + r*(Adim_angle*N_groups) + c;
                    
                    A_group[id_cell_group] = A_cell_group[id_group] - S[id_group];
                }
            }

        }

        print_rm(A_group);


    }



    return(0);
}
