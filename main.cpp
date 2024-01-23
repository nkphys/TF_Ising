#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;
#include "Matrix.h"
#include "Hamiltonian.h"
#include "Parameters.h"

#include "random"


int main(int argc, char *argv[]) {


    string ex_string_original =argv[0];
    cout<<"'"<<ex_string_original<<"'"<<endl;
//    string ex_string;
//    ex_string = ex_string_original.substr(ex_string_original.length()-5);
//    cout<<"'"<<ex_string<<"'"<<endl;


//        string ModelType = argv[1];
        string inputfile_ = argv[1];

//        if (argc<3) { throw std::invalid_argument("USE:: executable inputfile"); }



        Parameters Parameters_;
        Parameters_.Initialize(inputfile_);


        mt19937_64 Generator_(Parameters_.DisorderSeed); //for random fields


        Hamiltonian Hamiltonian_(Parameters_, Generator_);
        Hamiltonian_.Save_Static_Corr();


        Hamiltonian_.Save_B_matrices();

        Hamiltonian_.Get_Energy_and_kinkdensity();
        Hamiltonian_.Get_Entropy();
        Hamiltonian_.Get_MI();
        //Hamiltonian_.Bob_Sx_diff();
        Hamiltonian_.Get_2pointSzSzcorrelations();
        //Hamiltonian_.Get_4pointSxSxcorrelations();
        //Hamiltonian_.Get_4pointSxtSz0correlations();
        Hamiltonian_.Get_4pointSzSzcorrelations();


    cout << "--------THE END--------" << endl;
} // main
