#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:
    int NSites;


    string Scheduler_File_;
    string use_scheduler_;
    bool Use_Scheduler;

    Mat_1_real Time_bare;
    Mat_1_real Gamma_bare;
    Mat_1_real Js_bare;

    double Hx;
    double Jzz;

    double Disorder_Strength, RandomDisorderSeed;

    double dw_dos, eta_dos;
    double w_min, w_max;

    int No_TimeSlices;
    double dt_;
    double dt_measure;
    double t_anneal;
    double t_measure;
    double t_butterfly;
    int t_measure_slice;
    double Temperature,beta;

    bool ReadDisorder;
    string ReadDisorderString;
    string DisorderSeedFile;

    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){


    cout << "____________________________________" << endl;
    cout << "Reading the inputfile name: " << inputfile_ << endl;
    cout << "____________________________________" << endl;



    NSites = int(matchstring(inputfile_,"Total_sites"));

    RandomDisorderSeed = int(matchstring(inputfile_,"RandomDisorderSeed"));
    ReadDisorderString = matchstring2(inputfile_,"ReadDisorderConf");
    Disorder_Strength = matchstring(inputfile_,"Disorder_Strength");
    Temperature = matchstring(inputfile_,"Temperature");

    Hx = matchstring(inputfile_,"Hx");
    Jzz = matchstring(inputfile_,"Jzz");

    Scheduler_File_ = matchstring2(inputfile_,"Scheduler_File");
    use_scheduler_ = matchstring2(inputfile_,"Use_Scheduler");


    dt_ = matchstring(inputfile_,"dt");
    dt_measure = matchstring(inputfile_,"dt_measure");
    t_anneal = matchstring(inputfile_,"t_anneal");
    t_measure = matchstring(inputfile_,"t_measure");
    t_butterfly =matchstring(inputfile_,"t_butterfly");

    assert(dt_measure>dt_);

    t_measure_slice = int((t_measure/dt_)+0.5) + 1;
    No_TimeSlices = int((t_anneal/dt_)+0.5) + 1;

    //dw_dos, eta_dos
    //dw_dos = matchstring(inputfile_, "dw_dos");
    //eta_dos = matchstring(inputfile_, "eta_dos");

    Dflag = 'N';


    if(ReadDisorderString=="true"){
        ReadDisorder = true;
        DisorderSeedFile = matchstring2(inputfile_,"DisorderSeedFile");
    }
    else{
        ReadDisorder = false;
    }


    //beta=(11605.0/Temperature);
    beta=(1.0/Temperature);



    if(use_scheduler_ == "true"){
        Use_Scheduler=true;
        cout<<"Scheduler is used i.e. time dependent Hamiltonian"<<endl;
    }
    else{
        Use_Scheduler=false;
    }


    if(Use_Scheduler){
    string line2;
    double temp_t, temp_h, temp_J;
    ifstream scheduler_stream(Scheduler_File_.c_str());
    while(getline(scheduler_stream,line2)){
    stringstream line_ss(line2);
    line_ss>>temp_t>>temp_h>>temp_J;
    Time_bare.push_back(temp_t);
    Gamma_bare.push_back(temp_h);
    Js_bare.push_back(temp_J);
    }
    }





    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



