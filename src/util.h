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

    // Assuming square
    int N = sqrt(vec.size());

    if (N > 100){
        cout << ">>>Warning: Matrix is over 100x100<<<" << endl;
    }

    cout << "Matrix is of size ["<<N<<","<<N<<"]"<<endl;

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            printf("%5.2f ", vec[j*N+i]);
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

void print_vec_sd_int(std::vector<int> vec){
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

double relativeError(double a, double b, double max){
    return ( abs(a - b) / max);
}

std::vector<double> relativeError_ebe ( std::vector<double> v1, std::vector<double> v2){
    using namespace std;
    /*brei: computing the relative error between two vectors element by element*/

    if (v1.size() != v2.size()){
        std::cout << "WARNING: vectors are not of same size --relativeError_ebe" << std::endl;
    }

    // taking the absoulte value of all vectors
    //for (int i=0; i<v1.size(); i++){
    //    v1[i] = abs(v1[i]);
    //    v2[i] = abs(v2[i]);
    //}

    // find the maximum element between all vectors
    vector<double> max_vec = {*max_element(v2.begin(), v2.end()), *max_element(v1.begin(), v1.end())};
    double max = *max_element( max_vec.begin(), max_vec.end() ) ;

    // allocate an error vector
    vector<double> error(v1.size());

    // compute the relative error
    for (int i=0; i<v1.size(); i++){
        error[i] = relativeError(v1[i], v2[i], max);
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