#include<iostream>
#include <algorithm>    // std::copy
#include <cmath>
#include <vector>
#include <ctime>
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
void sim_anneal(vector<Molecule>& formal,vector<double>& pinit,vector<double>& ydat)
{
    // cout<<"Simulated Annealing is called"<<endl;
    double old_ene,new_ene,boltz,old;
    srand(time(NULL));
    //const double temp = 2;
    const double alpha=0.9;
    const double delta = 0.005;
    unsigned int item,accept=0,reject=0,moves=1000;
    for (double T = 8; T > 8e-5; T *= alpha) //T = T * alpha which used as a cooling schedule
    {
        double beta = 1/T;
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
        cout<<"Chisq is " <<new_ene<<"  At T = "<< T <<endl;
    }
}
