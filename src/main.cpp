/*brief: assembly functions for matricies to compute therefore commands
date: May 23rd 2023
auth: J Piper Morgan (morgjack@oregonstate.edu)*/

#include <iostream>
#include <vector>

// row major!!!!!!

class cell{
    public:
        int cell_id;
        float x_cell_center;
        float xsec_scatter;
        float xsec;
};

class problem_space{
    public:
        // physical space
        int N_cells; 
        int N_angles; 
        int N_time;
        float Length;
        float dt;
        float dx;
        float t_max;
        float material_source;
        float velocity;

        // comptuational
        int hardware_precision;
        float convergence_tolarance;
};

class boundary_condition{
    public:
        char side;
        int cell_id;
        int type;
        float magnitude;
};


// i space, m is angle, k is time, g is energy group

int main(){
    using namespace std;
    

    int N_cells = 4; 
    int N_angles = 2; 
    int N_time = 5;

    // 4 = N_subspace (2) * N_subtime (2)
    int N_mat = 4 * N_cells * N_angles;

    // N_cm is the size of the row major vetor
    int N_rm = N_mat*N_mat;

    vector<double> A_cm(N_rm);
    
    for (int i=0; i<N_cells; i++){
        for (int j=0; j<N_angles; j++){
            
        }

    }


}
