#include <iostream>
#include <vector>
#include "class.h"

// ROW MAJOR!!!!!!!!!


std::vector<double> A_neg_rm(cell cell, double mu, int group){
    double gamma = (cell.dx*cell.xsec_total[group])/2;
    double timer = cell.dx/(cell.v[group]*cell.dt);
    double timer2 = cell.dx/(2*cell.v[group]*cell.dt);
    double a = mu/2;

    std::vector<double> A_n = {-a + gamma, a,          timer2,            0,
                               -a,        -a + gamma,  0,                 timer2,
                               -timer,     0,          timer - a + gamma, a,
                                0,        -timer,     -a,                 timer -a + gamma};
    
    return(A_n);
    }

std::vector<double> A_pos_rm(cell cell, double mu, int group){
    double gamma = (cell.dx*cell.xsec_total[group])/2;
    double timer = cell.dx/(cell.v[group]*cell.dt);
    double timer2 = cell.dx/(2*cell.v[group]*cell.dt);
    double a = mu/2;

    std::vector<double> A_p = {a + gamma, a,         timer2,            0,
                              -a,         a + gamma, 0,                 timer2,
                              -timer,     0,         timer + a + gamma, a,
                               0,        -timer,    -a,                 timer +a + gamma};

    return(A_p);
    }

std::vector<double> scatter(cell cell, double *w, int N, int group){
    std::vector<double> S ((4*N*4*N));
    double beta = cell.dx*cell.xsec_scatter[group]/4;

    for (int ca=0; ca<N; ca++){
        for (int ra=0; ra<N; ra++){

            S[4*ra+0 + 4*4*ca*N + 0*N*4]   = beta*w[ra];
            S[4*ra+1 + 4*4*ca*N + 1*N*4] = beta*w[ra];
            S[4*ra+2 + 4*4*ca*N + 2*N*4] = beta*w[ra];
            S[4*ra+3 + 4*4*ca*N + 3*N*4] = beta*w[ra];
        }
    }

    return(S);
}