#include <iostream>
#include <vector>
#include "class.h"

// ROW MAJOR!!!!!!!!!


std::vector<double> A_neg(cell cell, float mu){
    float gamma = (cell.dx*cell.xsec_total)/2;
    float timer = cell.dx/(cell.v*cell.dt);
    float timer2 = cell.dx/(2*cell.v*cell.dt);
    float a = mu/2;

    std::vector<double> A_n = {-a + gamma, a,          timer2,            0,
                               -a,        -a + gamma,  0,                 timer2,
                               -timer,     0,          timer - a + gamma, a,
                                0,        -timer,     -a,                 timer -a + gamma};
    
    return(A_n);
    }

std::vector<double> A_pos(cell cell, float mu){
    float gamma = (cell.dx*cell.xsec_total)/2;
    float timer = cell.dx/(cell.v*cell.dt);
    float timer2 = cell.dx/(2*cell.v*cell.dt);
    float a = mu/2;

    std::vector<double> A_p = {a + gamma, a,         timer2,            0,
                              -a,         a + gamma, 0,                 timer2,
                              -timer,     0,         timer + a + gamma, a,
                               0,        -timer,    -a,                 timer +a + gamma};

    return(A_p);
    }