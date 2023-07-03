#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include<algorithm>



void rm2md(){

}

void print_rm(std::vector<double> vec){
    using namespace std;

    // Assuming square
    int N = sqrt(vec.size());

    if (N > 100){
        cout << ">>>Warning: Matrix is over 100x100<<<" << endl;
    }

    cout << "Matrix is of size ["<<N<<","<<N<<"]"<<endl;

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            printf("%5.2f ", vec[i*N+j]) ;
        }
        printf("\n");
    }
}



void print_cm(std::vector<double> vec){
    using namespace std;
    cout << "in print cm" << endl;
    // Assuming square
    int N = sqrt(vec.size());

    cout << "in print cm" << endl;

    if (N > 100){
        cout << ">>>Warning: Matrix is over 100x100<<<" << endl;
    }

    cout << "Matrix is of size ["<<N<<","<<N<<"]"<<endl;

    for (int row=0; row<N; row++){
        for (int col=0; col<N; col++){
            printf("%5.2f ", vec[row+N*col]);
        }
        printf("\n");
    }
}



void print_vec(int N, double *vec){
    using namespace std;

    cout << "" << endl;

    cout << "Size of Vector: "<< N << endl;

    for (int i=0; i<N; i++){
        cout << vec[i] << endl;
    }

    cout << "" << endl;
}



void print_vec_sd(std::vector<double> vec){
    using namespace std;

    int N = vec.size();

    cout << "" << endl;

    cout << "Size of Vector: "<< N << endl;

    for (int i=0; i<N; i++){
        cout << vec[i] << endl;
    }

    cout << "" << endl;
}

/*
int main(){
    std::vector<double> vec_dub = {0,1,2,3,4,5,6,7,8};
    std::vector<int> vec_int = {0,1,2,3,4,5,6,7,8};

    print_rm(vec_dub);

    return(0);
}*/


void print_title(){
    std::ifstream f("title.txt");

    if (f.is_open())
        std::cout << f.rdbuf();
}

double relativeError(double a, double b){
    return (abs(a - b) / b);
}

std::vector<double> relativeError_ebe ( std::vector<double> v1, std::vector<double> v2 ){
    /*brei: computing the relative error between two vectors element by element*/

    if (v1.size() != v2.size()){
        std::cout << "WARNING: vectors are not of same size --relativeError_ebe" << std::endl;
    }

    std::vector<double> error(v1.size());

    for (int i=0; i<v1.size(); i++){
        error[i] = relativeError(v1[i], v2[i]);
    }

    return(error);
}


double infNorm_error ( std::vector<double> v1, std::vector<double> v2 ){
    /*brief: computing the max relative error between two vectors*/

    using namespace std;

    std::vector<double> error = relativeError_ebe(v1, v2);

    double relError = *max_element(error.begin(), error.end());

    return(relError);

}

inline void outofbounds_check(int index, std::vector<double> &vec){
    using namespace std;

    if ( index < 0 ) {
        cout<<">>>>>>>>>>>>ERROR<<<<<<<<<<<<"<<endl;
        cout<<"sometihng was indexed under 0"<<endl;
    } else if ( index >= vec.size() ) {
        cout<<">>>>>>>>>>>>>>>>>>>>ERROR<<<<<<<<<<<<<<<<<<<<"<<endl;
        cout<<"sometihng was indexed over a vectors max size"<<endl;
    }
}