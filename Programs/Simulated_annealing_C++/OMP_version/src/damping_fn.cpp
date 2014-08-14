#include<iostream>
#include <algorithm>    // std::copy
#include <cmath>
#include <vector>       // std::vector
#include "molecule.h"
using namespace std;
double potent_df(Molecule& Mol,vector<double>& params)
{
    double potent = 0.0;
    vector<double>tbpar(7);
    omp_set_num_threads(2)
        #pragma omp parallel for reduction(+:potent)  schedule(static)
        {
            for (int nat=0; nat<16; nat++)
            {
                for (int nbt=26; nbt<30; nbt++)
                {
                    //cout<<nat<<"\t"<<nbt<<"\n";
                    if(nbt==26)
                    {
                        copy (params.begin(), params.begin()+7, tbpar.begin() );
                        potent=potent+tbd_damping(tbpar,Mol.Atoms[nat],Mol.Atoms[nbt]);
                    }
                    else if(nbt==27)
                    {
                        copy (params.begin()+7, params.begin()+14, tbpar.begin() );
                        potent=potent+tbd_damping(tbpar,Mol.Atoms[nat],Mol.Atoms[nbt]);
                    }
                    else if(nbt==28 || nbt==29)
                    {
                        copy (params.begin()+14, params.begin()+21, tbpar.begin() );
                        potent=potent+tbd_damping(tbpar,Mol.Atoms[nat],Mol.Atoms[nbt]);
                    }
                }
            }
            return potent;
        }
}
double tbd_damping(vector<double>&pars,Atom& A,Atom& B)
{
    double tbd,fh,rh;
    double dist=atom_dist(A,B);
    rh = pars[6];
    if(dist < rh)
    {
        fh = exp(-1*pow((rh/dist)-1,2));
    }
    else
    {
        fh = 1.0;
    }
    tbd=pars[0]*exp(-pars[1]*dist)+((pars[2]*pow(dist,-2))+(pars[3]*pow(dist,-4))+(pars[4]*pow(dist,-6))+(pars[5]*pow(dist,-8)))*fh;
    return tbd;
}
