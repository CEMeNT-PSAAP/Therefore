#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <fstream>


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


void print_vec(int N, double* vec){
    using namespace std;

    cout << "" << endl;

    cout << N << endl;

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