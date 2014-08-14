#include <iostream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "molecule.h"
int main()
{
    int nat,nmol,npar;
    vector<Molecule> formal;
    vector<double> params,rcc,rco,rcha,rchb,ydat;
    cin>>nat>>nmol>>npar;
    Read_params(params,npar);
    Read_file(formal,nat,nmol);
    Read_data(rcc,rco,rcha,rchb,ydat,nmol);
    sim_anneal(formal,params,ydat);
    /*PRINTING POTENTIAL TO VERIFY THE PROGRAMS ACCURACY*/
    ofstream logfile ("outlog");
    if (logfile.is_open())
    {
        for (int k=0; k<nmol; k++)
        {
            //   cout.fixed;
            logfile.precision(5);
            logfile<<rco.at(k)<<setw(15)<<potent_df(formal[k],params)<<setw(15)<<ydat.at(k)<<endl;
        }
        logfile.close();
    }
    ofstream parfile ("outpar");
    if (parfile.is_open())
    {
        for (int p=0; p<npar; p++)
        {
            parfile.precision(5);
            parfile<<p+1<<setw(15)<<params.at(p)<<endl;
        }
        parfile.close();
    }
    /*CLEARING VECTORS FROM MEMORY (DESTROYING THEM)*/
    formal.clear();
    params.clear();
    // cout<<"Successfully read this thing" <<endl;
    return 0;
}
