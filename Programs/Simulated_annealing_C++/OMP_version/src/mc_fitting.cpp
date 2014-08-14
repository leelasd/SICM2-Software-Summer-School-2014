#include<iostream>
#include<iomanip>
#include <algorithm>    // std::copy
#include <cmath>
#include <vector>
#include <ctime>
#include <omp.h>
#include "molecule.h"
using namespace std;

double Chisq(vector<Molecule>& formal,vector<double>& params,vector<double>& ydat)
{
    double chi=0.0;
    if(formal.size()==ydat.size())
    {
        for(unsigned i=0; i<ydat.size(); i++)
        {
            chi+=(potent_df(formal[i],params)-ydat.at(i))*(potent_df(formal[i],params)-ydat.at(i));
        }
    }
    else
    {
        cout<<"Size of data is not consistent"<<endl;
    }
    return chi;
}
void sim_anneal(vector<Molecule>& formal,vector<double>& pinit,vector<double>& ydat,int moves,double temp)
{
    // cout<<"Simulated Annealing is called"<<endl;
    double old_ene,new_ene,boltz,old;
    srand(time(NULL));
    //const double temp = 2;
    const double alpha=0.9;
    const double delta = 0.005;
    unsigned int item,accept=0,reject=0;
    for (double T = temp; T > temp*1e-6; T *= alpha) //T = T * alpha which used as a cooling schedule
    {
        double beta = 1/T;
    	double start = omp_get_wtime();
        for (unsigned int i=0; i<moves; i++)
        {
            //cout<<"Doing Iteration "<<i <<"  At Temp "<<T<<endl;
            old_ene = Chisq(formal,pinit,ydat);
            item = int(((double) rand() / (RAND_MAX))*pinit.size());
            old = pinit.at(item);
            pinit.at(item) = old + delta*(((double) rand() / (RAND_MAX))-0.5);
            new_ene = Chisq(formal,pinit,ydat);
            boltz = exp(-beta*(new_ene-old_ene));
            if(((double) rand() / (RAND_MAX))<boltz || new_ene<old_ene)
            {
                accept=accept+1;

            }
            else
            {
                reject=reject+1;
                pinit.at(item)=old;
            }
        }
    	double end = omp_get_wtime();
        cout<<"Chisq is " <<setprecision(3)<<new_ene<<setw(8)<<"  At T = "<<setprecision(6)<<T<<setw(8)<<" Time(sec) = "<<end-start<<endl;
    }
}
