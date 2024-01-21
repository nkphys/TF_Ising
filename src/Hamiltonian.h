#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters.h"
#include "pfapack.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);
//zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);



//dgesdd	(character JOBZ, integer M, integer N, double precision  dimension( lda, * ) A,
//             integer LDA, double precision dimension( * ) S, double precision dimension( ldu, * ) U,
//             integer LDU, double precision dimension( ldvt, * ) VT, integer LDVT,
//             double precision dimension( * ) WORK, integer LWORK, integer dimension( * ) IWORK,
//             integer INFO)
extern "C" void dgesdd_ (char *, int *, int *, double *, int *, double *, double *, int *, double *, int*,
                         double *, int *, int *, int *);

extern "C" void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*,
                         std::complex<double> *, int *, double * , int *, int *);


extern "C" void zgesvd_ (char *, char *, int *, int *, std::complex<double> *,
                         int *, double *, std::complex<double> *, int *, std::complex<double> *,
                         int*, std::complex<double> *, int *, double * , int *);

class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__)
        :Parameters_(Parameters__)
    {
        Initialize();
        Create_Scheduler();
        Create_M_mat();
    }


    void Initialize();    //::DONE 
    void Create_M_mat();   //::DONE
    void Diagonalize(char option);
    complex<double> Pfaffian(Matrix<complex<double>> & A);
    void Create_Scheduler();
    void Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_);

    void Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_);

    void Perform_SVD2_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_);

    complex<double> Get_Xi_2pointCorr(int site_1, int time_slice1, int site_2, int time_slice2);

    void Get_2pointSzSzcorrelations();
    void Get_4pointSzSzcorrelations();
    void Get_4pointSxSxcorrelations();
    void Get_4pointSxtSz0correlations();
    void Save_B_matrices();
    void Bob_Sx_diff();

    void Get_A_matrix(int time_slice, Mat_2_Complex_doub &A_mat);
    void Get_B_matrix(int t_slice);
    complex<double> Abs_Get_n_point_Xi_Corr(Mat_1_int t_array, Mat_1_int index_array);

    void Save_Static_Corr();
    double Fermi(double E);

    void Get_Energy_and_kinkdensity();

    void Get_Entropy();
    void Get_MI();
    void Diagonalize_for_eigvals(Matrix<complex<double>> &MatTemp, Mat_1_doub &EigvalTemp);
    Parameters &Parameters_;
    int ns_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Mat_1_doub eigs_;
    Matrix<double> M_mat;

    Mat_3_Complex_doub B_mat_slices;
    Mat_1_real Time_;
    Mat_1_real Gamma_;
    Mat_1_real Js_;

    int No_TimeSlices;
    double Time_max;
    double dt_;
    double t_butterfly;
    int t_butterfly_slice;

    Matrix<complex<double>> SzSzMat;

    Mat_2_real Xi_2point_corr;

    Mat_2_Complex_doub Static_Corr;


};


complex<double> Hamiltonian::Pfaffian(Matrix<complex<double>> & A){
      complex<double> pfaffian;
      int info=0;
      int N = A.n_col();
      vector<complex<double> > WORK(N);
      vector<double> RWORK(N-1);
      vector<int> IWORK(N);
      int LWORK;

      //workspace query for optimal workspace
      LWORK=-1;
      zskpfa_("L", "P", &N, &A(0,0), &N, &pfaffian, &IWORK[0],
              &WORK[0], &LWORK, &RWORK[0], &info);

      LWORK=static_cast<int>(real(WORK[0]));
      WORK.resize(LWORK);

      //now do the real calculation
      zskpfa_("L", "P", &N, &A(0,0), &N, &pfaffian, &IWORK[0],
              &WORK[0], &LWORK, &RWORK[0], &info);

      //cout << "The pfaffian is " << pfaffian << endl;

      return pfaffian;

}

void Hamiltonian::Get_A_matrix(int time_slice, Mat_2_Complex_doub &A_mat){


assert(time_slice>0);

A_mat.resize(2*ns_);
for(int i=0;i<2*ns_;i++){
A_mat[i].resize(2*ns_);
}


Ham_.resize(2*ns_, 2*ns_);

for(int m=0;m<2*ns_;m++){

    if(m%2==0){
    if((m+1)<2*ns_){
    Ham_(m,m+1) = -0.5*iota_complex*Parameters_.Hx*Gamma_[time_slice-1];
    }
    if((m-1)>=0){
    Ham_(m,m-1)=  0.5*iota_complex*Parameters_.Jzz*Js_[time_slice-1];
    }
    }
    else{ //m is odd
        if((m-1)>=0){
        Ham_(m,m-1) = 0.5*iota_complex*Parameters_.Hx*Gamma_[time_slice-1];
        }
        if((m+1)<2*ns_){
        Ham_(m,m+1)=  -0.5*iota_complex*Parameters_.Jzz*Js_[time_slice-1];
        }
    }

}

char Dflag = 'V';
Diagonalize(Dflag);



//Ham_(row_,n)
for(int n=0;n<2*ns_;n++){
for(int m=0;m<2*ns_;m++){
A_mat[n][m]=0.0;
for(int np=0;np<2*ns_;np++){
A_mat[n][m] += Ham_(n,np)*
               exp(-4.0*iota_complex*PI*dt_*(eigs_[np]))*
               conj(Ham_(m,np));
}
}
}

}



void Hamiltonian::Save_B_matrices(){

B_mat_slices.resize(Parameters_.t_measure_slice + 1);

for(int t_slice=0;t_slice<=Parameters_.t_measure_slice;t_slice++){
  Get_B_matrix(t_slice);
  cout<<"t_slice = "<<t_slice<<endl;
}

}

void Hamiltonian::Get_B_matrix(int t_slice){

B_mat_slices[t_slice].resize(2*ns_);
for(int i=0;i<2*ns_;i++){
B_mat_slices[t_slice][i].resize(2*ns_);
}

if(t_slice==0){
for(int i=0;i<2*ns_;i++){
for(int j=0;j<2*ns_;j++){
    if(i==j){
     B_mat_slices[t_slice][i][j]=1.0;
    }
    else{
     B_mat_slices[t_slice][i][j]=0.0;
    }
}
}
}
else{//time_slice>0
Mat_2_Complex_doub A_mat;
//Mat_2_Complex_doub B_mat_old;
//Get_B_matrix(time_slice-1);

Get_A_matrix(t_slice, A_mat);

for(int i=0;i<2*ns_;i++){
for(int j=0;j<2*ns_;j++){
    B_mat_slices[t_slice][i][j]=0.0;
    for(int k=0;k<2*ns_;k++){
    //B_mat_slices[t_slice][i][j]+=B_mat_slices[t_slice-1][i][k]*A_mat[k][j];
    B_mat_slices[t_slice][i][j]+=B_mat_slices[t_slice-1][k][j]*A_mat[i][k];
    }
}}


}


}


void Hamiltonian::Save_Static_Corr(){


Matrix<double> VT_, U_;
Mat_1_doub Sigma_;
Perform_SVD(M_mat, VT_, U_, Sigma_);


Matrix<double> UTU, UUT, VTV, VVT;

UTU.resize(ns_,ns_);UUT.resize(ns_,ns_);
VTV.resize(ns_,ns_);VVT.resize(ns_,ns_);


for(int i=0;i<ns_;i++){
    for(int j=0;j<ns_;j++){
        UTU(i,j)=0;
        for(int n=0;n<ns_;n++){
        UTU(i,j) += U_(n,i)*U_(n,j);
        UUT(i,j) += U_(i,n)*U_(j,n);

        VTV(i,j) += VT_(i,n)*VT_(j,n);
        VVT(i,j) += VT_(n,i)*VT_(n,j);
        }
    }
}

string UTU_file_str="UTU.txt";
string UUT_file_str="UUT.txt";
string VTV_file_str="VTV.txt";
string VVT_file_str="VVT.txt";

ofstream UTU_file(UTU_file_str.c_str());
ofstream UUT_file(UUT_file_str.c_str());
ofstream VTV_file(VTV_file_str.c_str());
ofstream VVT_file(VVT_file_str.c_str());

for(int i=0;i<ns_;i++){
    for(int j=0;j<ns_;j++){
        UTU_file<<UTU(i,j)<<"  ";
        UUT_file<<UUT(i,j)<<"  ";
        VTV_file<<VTV(i,j)<<"  ";
        VVT_file<<VVT(i,j)<<"  ";
    }
    UTU_file<<endl;
    UUT_file<<endl;
    VTV_file<<endl;
    VVT_file<<endl;
}


Static_Corr.resize(2*ns_);
for(int i=0;i<2*ns_;i++){
Static_Corr[i].resize(2*ns_);
}


for(int mp=0;mp<2*ns_;mp++){
    for(int np=0;np<2*ns_;np++){

        int site_1 = int((1.0*mp+0.5)/2.0);
        int type_1 = mp - 2*site_1;

        int site_2 = int((1.0*np+0.5)/2.0);
        int type_2 = np - 2*site_2;

        if(type_1==type_2){
            if(site_1==site_2){
              Static_Corr[mp][np]=1.0;
            }
            else{
              Static_Corr[mp][np]=0.0;
            }
        }
        else{
          Static_Corr[mp][np]=0;

          if(type_1==0 && type_2==1){
          for(int n=0;n<ns_;n++){
           Static_Corr[mp][np] += iota_complex*VT_(n,site_1)*U_(site_2,n)*(1.0-2.0*Fermi(2.0*Sigma_[n]));
          }
          }
          else{
          for(int n=0;n<ns_;n++){
           Static_Corr[mp][np] -= iota_complex*VT_(n,site_2)*U_(site_1,n)*(1.0-2.0*Fermi(2.0*Sigma_[n]));
          }
          }
        }

    }
}

string fileout_str="Static_Corr.txt";
ofstream fileout_(fileout_str.c_str());
for(int i=0;i<2*ns_;i++){
    for(int j=0;j<2*ns_;j++){
    fileout_<<i<< "  "<<j<<"  "<<Static_Corr[i][j].real()<<"  "<<Static_Corr[i][j].imag()<<endl;
    }
    fileout_<<endl;
}


}

double Hamiltonian::Fermi(double E){

    return ( 1.0 - (1.0/( exp((E)*Parameters_.beta) + 1.0)));

}

complex<double> Hamiltonian::Get_Xi_2pointCorr(int index_1, int time_slice1, int index_2, int time_slice2){

complex<double> val;
//Mat_2_Complex_doub B1_mat, B2_mat;
//Xi majorans's --> a (type=0) ,b(type=1) majorana
//a0    b0    a1    b1    a2    b2 ----------- a_{ns-1}    b_{ns-1}
//Xi0   Xi1   Xi2   Xi3   Xi4   Xi5----------- Xi_{2ns-2}  Xi_{2ns -1}
//type + site*2 = index

int TYPE_a=0;
int TYPE_b=1;

//Get_B_matrix(time_slice1, B1_mat);
//Get_B_matrix(time_slice2, B2_mat);

val=0.0;
int np;
for(int mp=0;mp<2*ns_;mp++){
    //for(int np=0;np<2*ns_;np++){
    np=mp;
       val += B_mat_slices[time_slice1][index_1][mp]*B_mat_slices[time_slice2][index_2][np]*Static_Corr[mp][np];

   np=mp+1;
   if(np<2*ns_){
   val += B_mat_slices[time_slice1][index_1][mp]*B_mat_slices[time_slice2][index_2][np]*Static_Corr[mp][np];
    }
    np=mp-1;
    if(np>=0){
    val += B_mat_slices[time_slice1][index_1][mp]*B_mat_slices[time_slice2][index_2][np]*Static_Corr[mp][np];
    }
       //}
}

return val;
}


void Hamiltonian::Initialize(){

    ns_=Parameters_.NSites;
    M_mat.resize(ns_,ns_);

    No_TimeSlices = Parameters_.No_TimeSlices;
    Time_max = Parameters_.t_anneal;
    dt_ = Parameters_.dt_;

    t_butterfly=Parameters_.t_butterfly;
    t_butterfly_slice= int((t_butterfly/dt_)+0.5);
} // ----------




void Hamiltonian::Diagonalize_for_eigvals(Matrix<complex<double>> &MatTemp, Mat_1_doub &EigvalTemp){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz='N';
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=MatTemp.n_row();
    int lda=MatTemp.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    EigvalTemp.resize(MatTemp.n_row());
    fill(EigvalTemp.begin(),EigvalTemp.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(MatTemp(0,0)),&lda,&(EigvalTemp[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(MatTemp(0,0)),&lda,&(EigvalTemp[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}



void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}


void Hamiltonian::Create_M_mat(){


    for(int i=0;i<ns_;i++){
     M_mat(i,i)=Parameters_.Hx*Gamma_[0];
     if(i<ns_-1){//OBC
     M_mat(i,i+1)=-1.0*Parameters_.Jzz*(Js_[0]);
     }
    }


}




void Hamiltonian::Create_Scheduler(){

    //assert(abs(Restart_Time)<0.0000000001);

    cout<<"No of Time slices = "<<Parameters_.No_TimeSlices<<endl;

    Time_.resize(No_TimeSlices);
    Gamma_.resize(No_TimeSlices);
    Js_.resize(No_TimeSlices);

    double time_normalized;
    double gamma_temp, js_temp;
    for(int time_ind=0;time_ind<No_TimeSlices;time_ind++){

        time_normalized=(time_ind*dt_)/(Time_max);

       // cout<<time_ind<<"  "<<time_normalized<<endl;
        //assert(time_normalized>=0 && time_normalized<=1.0);

        for(int time_ind2=0;time_ind2<Parameters_.Time_bare.size()-1;time_ind2++){
            if(time_normalized>=Parameters_.Time_bare[time_ind2] && time_normalized<=Parameters_.Time_bare[time_ind2+1] ){
                gamma_temp = (abs(Parameters_.Time_bare[time_ind2+1]-time_normalized)*Parameters_.Gamma_bare[time_ind2]
                        +abs(Parameters_.Time_bare[time_ind2]-time_normalized)*Parameters_.Gamma_bare[time_ind2+1])*
                        (1.0/(Parameters_.Time_bare[time_ind2+1]-Parameters_.Time_bare[time_ind2]));

                js_temp = (abs(Parameters_.Time_bare[time_ind2+1]-time_normalized)*Parameters_.Js_bare[time_ind2]
                        +abs(Parameters_.Time_bare[time_ind2]-time_normalized)*Parameters_.Js_bare[time_ind2+1])*
                        (1.0/(Parameters_.Time_bare[time_ind2+1]-Parameters_.Time_bare[time_ind2]));

                break;
            }
        }

        Time_[time_ind]=time_ind*dt_;
        Gamma_[time_ind]=gamma_temp;
        Js_[time_ind]=js_temp;
    }


    string created_schd_str = "Created_Scheduler.txt";
    ofstream created_schd_stream(created_schd_str.c_str());

    created_schd_stream<<"# time   Gamma    Js"<<endl;
    for(int time_ind=0;time_ind<Time_.size();time_ind++){
        created_schd_stream<<Time_[time_ind]<<"  "<<Gamma_[time_ind]<<"  "<<Js_[time_ind]<<endl;
    }

}


void Hamiltonian::Get_Energy_and_kinkdensity(){


   // Get_Xi_2pointCorr(int index_1, int time_slice1, int index_2, int time_slice2);

    complex<double> kinkdensity;
    string fileout = "Energy_and_kinkdensity.txt";
    ofstream fileout_str(fileout.c_str());
    fileout_str<<"#time Energy(t)  kinkdensity(t)"<<endl;

    complex<double> energy;
    for(int time_slice=0;time_slice<Parameters_.t_measure_slice;time_slice++){

    energy=0.0;
    kinkdensity=0.0;
    for(int n=0;n<ns_-1;n++){
        energy += -1.0*iota_complex*Parameters_.Jzz*Js_[time_slice]*
                   Get_Xi_2pointCorr(2*n+1, time_slice, 2*n+2,time_slice);
        kinkdensity +=(1.0/ns_)*((1.0*one_complex) + iota_complex*Get_Xi_2pointCorr(2*n+1, time_slice, 2*n+2,time_slice))*0.5;
    }

    for(int n=0;n<ns_;n++){
        energy += 1.0*iota_complex*Parameters_.Hx*Gamma_[time_slice]*
                   Get_Xi_2pointCorr(2*n+1, time_slice, 2*n,time_slice);
    }



    fileout_str<<time_slice*dt_<<"   "<<energy.real()<<"  "<<energy.imag()<<"   "<<kinkdensity.real()<<"  "<<kinkdensity.imag()<<endl;

    }

}


void Hamiltonian::Bob_Sx_diff(){

int time_slice_1=t_butterfly_slice;
int site_1=(ns_/2);

//<Sx_site2(t2)> - <Sz_site1(0)Sx_site2(t2)Sz_site1(0)>

complex<double> No_AliceFlip, With_AliceFlip;
complex<double> Diff_;
//NEEDS PFAFFIAN

Mat_1_int time_array, index_array;

int fraction = int ((Parameters_.dt_measure/Parameters_.dt_)+0.5);


string fileout = "Bob_Sx_measurement.txt";
ofstream fileout_str(fileout.c_str());
fileout_str<<"#time1  site1   time2  site2  Diff.r Diff.i  No_Aflip.r No_Aflip.i  With_Aflip.r  With_Aflip.i"<<endl;


for(int time_slice_2=t_butterfly_slice;time_slice_2<=Parameters_.t_measure_slice;time_slice_2=time_slice_2+fraction){

for(int site_2=0;site_2<ns_;site_2++){
time_array.clear();index_array.clear();


time_array.push_back(0);index_array.push_back(2*site_1);//a(site_1)
for(int m=0;m<site_1;m++){
time_array.push_back(0);
index_array.push_back(2*m + 1); //b(m)
time_array.push_back(0);
index_array.push_back(2*m); //a(m)
}

//time_array.push_back(time_slice_2);index_array.push_back(2*site_2);//a(site_1)
//for(int m=0;m<site_2;m++){
//time_array.push_back(time_slice_2);
//index_array.push_back(2*m + 1); //b(m)
//time_array.push_back(time_slice_2);
//index_array.push_back(2*m); //a(m)
//}

time_array.push_back(time_slice_2);
index_array.push_back(2*site_2 + 1); //b(m)
time_array.push_back(time_slice_2);
index_array.push_back(2*site_2); //a(m)



time_array.push_back(0);index_array.push_back(2*site_1);//a(site_1)
for(int m=0;m<site_1;m++){
time_array.push_back(0);
index_array.push_back(2*m + 1); //b(m)
time_array.push_back(0);
index_array.push_back(2*m); //a(m)
}


With_AliceFlip=iota_complex*Abs_Get_n_point_Xi_Corr(time_array, index_array);

time_array.clear();index_array.clear();
time_array.push_back(time_slice_2);
index_array.push_back(2*site_2 + 1); //b(m)
time_array.push_back(time_slice_2);
index_array.push_back(2*site_2); //a(m)

No_AliceFlip=iota_complex*Abs_Get_n_point_Xi_Corr(time_array, index_array);
Diff_= No_AliceFlip -With_AliceFlip;

fileout_str<<setprecision(10)<<time_slice_1*dt_<<"   "<<site_1<<"   "<<time_slice_2*dt_<<"   "<<site_2<<"   "
          <<Diff_.real()<<"  "<<Diff_.imag()<<"   "
          <<No_AliceFlip.real()<<"   "<<No_AliceFlip.imag()<<"   "
          <<With_AliceFlip.real()<<"   "<<With_AliceFlip.imag()<<endl;

}
fileout_str<<endl;
}

}


/*
void Hamiltonian::Get_static_2pointSzSzCorrelationMatrix(){

    SzSzMat.resize(ns_,ns_);
    Mat_1_int time_array, index_array;

    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){

        }
    }

}
*/



void Hamiltonian::Get_Entropy(){

    int n_SitesA;
    double Entropy_A;
    double eps=1e-10;
   // cout<<"eps = "<<eps<<endl;
   // cout<<"log2(8) = "<<log2(8.0)<<endl;

    Mat_1_int time_array, index_array;

    Matrix<complex<double>> Cmat;
    Mat_1_doub vm_;

    int fraction = int ((Parameters_.dt_measure/Parameters_.dt_)+0.5);

    string fileout = "Entropy_l_and_Nminusl.txt";
    ofstream fileout_str(fileout.c_str());
    fileout_str<<"#time   n_sites[0-l]  Entropy[0-l]"<<endl;

    for(int time_slice_1=0;time_slice_1<=Parameters_.t_measure_slice;time_slice_1=time_slice_1+fraction){

        for(int site_l=0;site_l<=(ns_/2);site_l++){
        n_SitesA=site_l+1;

        Cmat.resize(n_SitesA*2,n_SitesA*2);
        Cmat.fill(0.0);

        for(int i=0;i<2*n_SitesA;i++){
            for(int j=0;j<2*n_SitesA;j++){
                Cmat(i,j)=Get_Xi_2pointCorr(i,time_slice_1,j, time_slice_1);
            }
        }

        Diagonalize_for_eigvals(Cmat, vm_);

        Entropy_A=0;


        for(int n=0;n<vm_.size();n++){
            if(vm_[n]>0){
            double pm=(vm_[n])/2.0;
            Entropy_A += (-1.0*pm)*log(pm);
            }
        }


        fileout_str<<time_slice_1*dt_<<"   "<< n_SitesA<<"  "<<Entropy_A<<endl;


        }
        fileout_str<<endl;

    }

}



void Hamiltonian::Get_MI(){

    assert(ns_%2==0);

    int n_SitesB;
    int n_Sites_ApB;
    double Entropy_A; //[N/2-l .... N/2-1]
    double Entropy_B; //[N/2 .... N/2+l-1]
    double Entropy_ApB;
    double MI_bw_A_and_B;
    double eps=1e-10;
    int ip,jp;
   // cout<<"eps = "<<eps<<endl;
   // cout<<"log2(8) = "<<log2(8.0)<<endl;

    Mat_1_int time_array, index_array;

    Matrix<complex<double>> CmatA, CmatB, CmatApB;
    Mat_1_doub vm_;
    double pm;

    int fraction = int((Parameters_.dt_measure/Parameters_.dt_)+0.5);

    string fileout = "MI_bw_A_and_B.txt";
    ofstream fileout_str(fileout.c_str());
    fileout_str<<"#time   size(A)=l  MI_A_and_B    SA    SB   SAB"<<endl;

    for(int time_slice_1=0;time_slice_1<=Parameters_.t_measure_slice;time_slice_1=time_slice_1+fraction){

        for(int n_SitesA=1;n_SitesA<=(ns_/2);n_SitesA++){

        CmatA.resize(n_SitesA*2,n_SitesA*2);
        CmatA.fill(0.0);
        for(int i=(ns_/2)-n_SitesA;i<=(ns_/2)-1;i++){
            for(int j=(ns_/2)-n_SitesA;j<=(ns_/2)-1;j++){

                ip=i-((ns_/2)-n_SitesA);
                jp=j-((ns_/2)-n_SitesA);

                for(int type1=0;type1<2;type1++){
                for(int type2=0;type2<2;type2++){
                CmatA(2*ip+type1,2*jp+type2)=Get_Xi_2pointCorr(2*i+type1,time_slice_1,2*j+type2, time_slice_1);
                }}

           }
        }


        CmatB.resize(n_SitesA*2,n_SitesA*2);
        CmatB.fill(0.0);
        for(int i=(ns_/2);i<=(ns_/2)+n_SitesA-1;i++){
            for(int j=(ns_/2);j<=(ns_/2)+n_SitesA-1;j++){
                ip = i- (ns_/2);
                jp = j- (ns_/2);
                for(int type1=0;type1<2;type1++){
                for(int type2=0;type2<2;type2++){
                CmatB(2*ip+type1,2*jp+type2)=Get_Xi_2pointCorr(2*i+type1,time_slice_1,2*j+type2, time_slice_1);
                }}
           }
        }


        CmatApB.resize(n_SitesA*4,n_SitesA*4);
        CmatApB.fill(0.0);
        for(int i=(ns_/2)-n_SitesA;i<=(ns_/2)+n_SitesA-1;i++){
            for(int j=(ns_/2)-n_SitesA;j<=(ns_/2)+n_SitesA-1;j++){
                ip = i - ((ns_/2)-n_SitesA);
                jp = j - ((ns_/2)-n_SitesA);
                for(int type1=0;type1<2;type1++){
                for(int type2=0;type2<2;type2++){
                CmatApB(2*ip+type1,2*jp+type2)=Get_Xi_2pointCorr(2*i+type1,time_slice_1,2*j+type2, time_slice_1);
                }}
           }
        }

        Diagonalize_for_eigvals(CmatA, vm_);
        Entropy_A=0;
        for(int n=0;n<vm_.size();n++){
            if(vm_[n]>0){
            pm=(vm_[n])/2.0;
            Entropy_A += (-1.0*pm)*log(pm);
            }
        }


        Diagonalize_for_eigvals(CmatB, vm_);
        Entropy_B=0;
        for(int n=0;n<vm_.size();n++){
            if(vm_[n]>0){
            pm=(vm_[n])/2.0;
            Entropy_B += (-1.0*pm)*log(pm);
            }
        }

        Diagonalize_for_eigvals(CmatApB, vm_);
        Entropy_ApB=0;
        for(int n=0;n<vm_.size();n++){
            if(vm_[n]>0){
            pm=(vm_[n])/2.0;
            Entropy_ApB += (-1.0*pm)*log(pm);
            }
        }


        MI_bw_A_and_B = Entropy_A + Entropy_B -Entropy_ApB;
        fileout_str<<time_slice_1*dt_<<"   "<< n_SitesA<<"  "<<MI_bw_A_and_B<<"   "<<Entropy_A<<"  "<<Entropy_B<<"  "<<Entropy_ApB<<endl;


        }
        fileout_str<<endl;

    }

}

void Hamiltonian::Get_2pointSzSzcorrelations(){

int site_1=(ns_/2);
int time_slice_2=t_butterfly_slice;

Mat_1_int time_array, index_array;

int fraction = int ((Parameters_.dt_measure/Parameters_.dt_)+0.5);


string fileout = "Sz0_Szt_Corr.txt";
ofstream fileout_str(fileout.c_str());
fileout_str<<"#time1 site1  time2  site2  <Sz1(t1)Sz2(t2)>"<<endl;

for(int time_slice_1=t_butterfly_slice;time_slice_1<=Parameters_.t_measure_slice;time_slice_1=time_slice_1+fraction){

for(int site_2=0;site_2<ns_;site_2++){
time_array.clear();index_array.clear();


time_array.push_back(time_slice_1);index_array.push_back(2*site_1);//a(site_1)
for(int m=0;m<site_1;m++){
time_array.push_back(time_slice_1);
index_array.push_back(2*m + 1); //b(m)
time_array.push_back(time_slice_1);
index_array.push_back(2*m); //a(m)
}

time_array.push_back(time_slice_2);index_array.push_back(2*site_2); //a(site_2)
for(int n=0;n<site_2;n++){
time_array.push_back(time_slice_2);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_2);
index_array.push_back(2*n); //a(n)
}


//cout<<time_slice_1<<"   "<<site_1<<"   "<<0.0<<"   "<<site_2<<endl;
fileout_str<<setprecision(10)<<time_slice_1*dt_<<"   "<<site_1<<"   "<<0.0<<"   "<<site_2<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).real()<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).imag()<<endl;
}

fileout_str<<endl;

cout<<time_slice_1<<" done SzSz"<<endl;
}



string fileout2 = "Szt_Szt_Corr.txt";
ofstream fileout2_str(fileout2.c_str());
fileout2_str<<"#time site1  time  site2  <Sz1(t)Sz2(t)>"<<endl;


for(int time_slice_1=0;time_slice_1<=Parameters_.t_measure_slice;time_slice_1=time_slice_1+fraction){

for(int site_2=0;site_2<ns_;site_2++){
time_array.clear();index_array.clear();


time_array.push_back(time_slice_1);index_array.push_back(2*site_1);//a(site_1)
for(int m=0;m<site_1;m++){
time_array.push_back(time_slice_1);
index_array.push_back(2*m + 1); //b(m)
time_array.push_back(time_slice_1);
index_array.push_back(2*m); //a(m)
}

time_array.push_back(time_slice_1);index_array.push_back(2*site_2); //a(site_2)
for(int n=0;n<site_2;n++){
time_array.push_back(time_slice_1);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_1);
index_array.push_back(2*n); //a(n)
}


//cout<<time_slice_1<<"   "<<site_1<<"   "<<0.0<<"   "<<site_2<<endl;
fileout2_str<<setprecision(10)<<time_slice_1*dt_<<"   "<<site_1<<"   "<<0.0<<"   "<<site_2<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).real()<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).imag()<<endl;
}

fileout2_str<<endl;

cout<<time_slice_1<<" done SztSzt"<<endl;
}




}





void Hamiltonian::Get_4pointSzSzcorrelations(){

int site_1=(ns_/2);
int site_3=site_1;

int time_slice_2=t_butterfly_slice;
int time_slice_4=t_butterfly_slice;

int fraction = int ((Parameters_.dt_measure/Parameters_.dt_)+0.5);
Mat_1_int time_array, index_array;


string fileout = "Szt_Sz0_Szt_Sz0Corr.txt";
ofstream fileout_str(fileout.c_str());
fileout_str<<"#time1 site1  time2  site2 time3  site3   time4   site4 <Sz1(t1)Sz2(t2)Sz3(t3)Sz4(t4)>"<<endl;

for(int time_slice_1=t_butterfly_slice;time_slice_1<=Parameters_.t_measure_slice;time_slice_1=time_slice_1+fraction){
int time_slice_3=time_slice_1;

for(int site_2=0;site_2<ns_;site_2++){
int site_4=site_2;
time_array.clear();index_array.clear();


time_array.push_back(time_slice_1);index_array.push_back(2*site_1);//a(site_1)
for(int m=0;m<site_1;m++){
time_array.push_back(time_slice_1);
index_array.push_back(2*m + 1); //b(m)
time_array.push_back(time_slice_1);
index_array.push_back(2*m); //a(m)
}

time_array.push_back(time_slice_2);index_array.push_back(2*site_2); //a(site_2)
for(int n=0;n<site_2;n++){
time_array.push_back(time_slice_2);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_2);
index_array.push_back(2*n); //a(n)
}

time_array.push_back(time_slice_3);index_array.push_back(2*site_3); //a(site_2)
for(int n=0;n<site_3;n++){
time_array.push_back(time_slice_3);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_3);
index_array.push_back(2*n); //a(n)
}

time_array.push_back(time_slice_4);index_array.push_back(2*site_4); //a(site_2)
for(int n=0;n<site_4;n++){
time_array.push_back(time_slice_4);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_4);
index_array.push_back(2*n); //a(n)
}



fileout_str<<setprecision(10)<<time_slice_1*dt_<<"   "<<site_1<<"   "<<0.0<<"   "<<site_2<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).real()<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).imag()<<endl;
}

fileout_str<<endl;
cout<<time_slice_1<<" done OTOC"<<endl;
}

}




void Hamiltonian::Get_4pointSxtSz0correlations(){

int site_1=(ns_/2);
int site_3=site_1;

int time_slice_2=t_butterfly_slice;
int time_slice_4=t_butterfly_slice;

int fraction = int ((Parameters_.dt_measure/Parameters_.dt_)+0.5);
Mat_1_int time_array, index_array;


string fileout = "Sxt_Sz0_Sxt_Sz0Corr.txt";
ofstream fileout_str(fileout.c_str());
fileout_str<<"#time1 site1  time2  site2 time3  site3   time4   site4 <Sz1(t1)Sz2(t2)Sz3(t3)Sz4(t4)>"<<endl;
int m;
for(int time_slice_1=t_butterfly_slice;time_slice_1<=Parameters_.t_measure_slice;time_slice_1=time_slice_1+fraction){
int time_slice_3=time_slice_1;

for(int site_2=0;site_2<ns_;site_2++){
int site_4=site_2;
time_array.clear();index_array.clear();

m=site_1;
time_array.push_back(time_slice_1);
index_array.push_back(2*m + 1); //b(m)
time_array.push_back(time_slice_1);
index_array.push_back(2*m); //a(m)


time_array.push_back(time_slice_2);index_array.push_back(2*site_2); //a(site_2)
for(int n=0;n<site_2;n++){
time_array.push_back(time_slice_2);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_2);
index_array.push_back(2*n); //a(n)
}

m=site_3;
time_array.push_back(time_slice_3);
index_array.push_back(2*m + 1); //b(n)
time_array.push_back(time_slice_3);
index_array.push_back(2*m); //a(n)


time_array.push_back(time_slice_4);index_array.push_back(2*site_4); //a(site_2)
for(int n=0;n<site_4;n++){
time_array.push_back(time_slice_4);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_4);
index_array.push_back(2*n); //a(n)
}




fileout_str<<setprecision(10)<<time_slice_1*dt_<<"   "<<site_1<<"   "<<0.0<<"   "<<site_2<<"  "<<-1.0*Abs_Get_n_point_Xi_Corr(time_array, index_array).real()<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).imag()<<endl;
}

fileout_str<<endl;
cout<<time_slice_1<<" done OTOC Sxt_Sz0_Sxt_Sz0"<<endl;
}

}



void Hamiltonian::Get_4pointSxSxcorrelations(){

int site_1=(ns_/2);
int site_3=site_1;

int time_slice_2=t_butterfly_slice;
int time_slice_4=t_butterfly_slice;

int fraction = int ((Parameters_.dt_measure/Parameters_.dt_)+0.5);
Mat_1_int time_array, index_array;


string fileout = "Sxt_Sx0_Sxt_Sx0Corr.txt";
ofstream fileout_str(fileout.c_str());
fileout_str<<"#time1 site1  time2  site2 time3  site3   time4   site4 <Sz1(t1)Sz2(t2)Sz3(t3)Sz4(t4)>"<<endl;
int m,n;
for(int time_slice_1=t_butterfly_slice;time_slice_1<=Parameters_.t_measure_slice;time_slice_1=time_slice_1+fraction){
int time_slice_3=time_slice_1;

for(int site_2=0;site_2<ns_;site_2++){
int site_4=site_2;
time_array.clear();index_array.clear();

m=site_1;
time_array.push_back(time_slice_1);
index_array.push_back(2*m + 1); //b(m)
time_array.push_back(time_slice_1);
index_array.push_back(2*m); //a(m)


n=site_2;
time_array.push_back(time_slice_2);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_2);
index_array.push_back(2*n); //a(n)

n=site_3;
time_array.push_back(time_slice_3);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_3);
index_array.push_back(2*n); //a(n)


n=site_4;
time_array.push_back(time_slice_4);
index_array.push_back(2*n + 1); //b(n)
time_array.push_back(time_slice_4);
index_array.push_back(2*n); //a(n)




fileout_str<<setprecision(10)<<time_slice_1*dt_<<"   "<<site_1<<"   "<<0.0<<"   "<<site_2<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).real()<<"  "<<Abs_Get_n_point_Xi_Corr(time_array, index_array).imag()<<endl;
}

fileout_str<<endl;
cout<<time_slice_1<<" done OTOCxxxx"<<endl;
}

}






complex<double> Hamiltonian::Abs_Get_n_point_Xi_Corr(Mat_1_int t_array, Mat_1_int index_array){



assert(t_array.size()==index_array.size());

Matrix<complex<double>> A_mat;
A_mat.resize(t_array.size(),t_array.size());
for(int i=0;i<A_mat.n_row();i++){
    for(int j=0;j<A_mat.n_row();j++){
    A_mat(i,j)=0.0;
    }
}

for(int n=0;n<A_mat.n_row();n++){
    for(int m=0;m<A_mat.n_col();m++){
        if(n<m){
         A_mat(n,m)=Get_Xi_2pointCorr(index_array[n], t_array[n], index_array[m], t_array[m]);
         A_mat(m,n)=-1.0*A_mat(n,m);
        }
    }
}


//string fileout_str="A_mat.txt";
//ofstream fileout_(fileout_str.c_str());
//for(int i=0;i<A_mat.n_row();i++){
//    for(int j=0;j<A_mat.n_row();j++){
//    fileout_<<i<< "  "<<j<<"  "<<A_mat(i,j).real()<<"  "<<A_mat(i,j).imag()<<endl;
//    }
//    fileout_<<endl;
//}


/*
Mat_1_doub Sigma_;
Matrix<complex<double>> VT_, U_;
Perform_SVD2_complex(A_mat, VT_, U_, Sigma_);

double Det_=1.0;
for(int i=0;i<Sigma_.size();i++){
    Det_ = Det_*Sigma_[i];
}
Det_=abs(Det_);

return sqrt(Det_);
*/
complex<double> pfaffian_;
pfaffian_=Pfaffian(A_mat);

return pfaffian_;
}

void Hamiltonian::Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_){


    char jobz='A'; //A,S,O,N

    int m=A_.n_row();
    int n=A_.n_col();
    int lda=A_.n_row();
    int ldu=A_.n_row();
    int ldvt=n;

    Sigma_.clear();
    Sigma_.resize(min(m,n));

    U_.resize(ldu,m);

    VT_.resize(ldvt,n);

    vector<double> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min(m,n));

    // query:
    dgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(iwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
    // real work:
    dgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(iwork[0]), &info);
    if (info!=0) {
        if(info>0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info>0.\n");}
        if(info<0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info<0.\n");
        }
    }


}



void Hamiltonian::Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_){


    char jobz='A'; //A,S,O,N

    int m=A_.n_row();
    int n=A_.n_col();
    int lda=A_.n_row();
    int ldu=A_.n_row();
    int ldvt=n;

    Sigma_.clear();
    Sigma_.resize(min(m,n));

    U_.resize(ldu,m);

    VT_.resize(ldvt,n);


    vector<complex<double>> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min(m,n));
    int lrwork = max( (5*min(m,n)*min(m,n)) + 7*min(m,n), (2*max(m,n)*min(m,n)) + (2*min(m,n)*min(m,n)) + min(m,n) )+10;
    lrwork=2*lrwork;
    vector<double> rwork(lrwork);

    // query:
    zgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(rwork[0]), &(iwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real())+10;
    work.resize(lwork);
    // real work:
    zgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(rwork[0]), &(iwork[0]), &info);
    if (info!=0) {
        if(info>0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zgesdd: failed with info>0.\n");}
        if(info<0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zgesdd: failed with info<0.\n");
        }
    }

    // Ham_.print();



}



void Hamiltonian::Perform_SVD2_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_){


    char jobu='A'; //A,S,O,N
    char jobv='A';
    int m=A_.n_row();
    int n=A_.n_col();
    int lda=A_.n_row();
    int ldu=A_.n_row();
    int ldvt=n;

    Sigma_.clear();
    Sigma_.resize(min(m,n));

    U_.resize(ldu,m);

    VT_.resize(ldvt,n);


    vector<complex<double>> work(3);
    int info;
    int lwork= -1;
    int lrwork = (5*min(m,n));
    vector<double> rwork(lrwork);

    // query://here
    zgesvd_(&jobu, &jobv, &m, &n, &(A_(0,0)),
            &lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)),
            &ldvt, &(work[0]), &lwork, &(rwork[0]), &info);

    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real())+10;
    work.resize(lwork);
    // real work:
    zgesvd_(&jobu, &jobv, &m, &n, &(A_(0,0)),
            &lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)),
            &ldvt, &(work[0]), &lwork, &(rwork[0]), &info);
    if (info!=0) {
        if(info>0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zgesdd: failed with info>0.\n");}
        if(info<0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zgesdd: failed with info<0.\n");
        }
    }

    // Ham_.print();



}


// ----------


#endif
