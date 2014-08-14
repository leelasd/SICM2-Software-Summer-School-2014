#include<iostream>
#include <algorithm>    // std::copy
#include <cmath>
#include <vector>       // std::vector
#include "molecule.h"
using namespace std;
double potent(Molecule& Mol,vector<double>& params)
{
    double potent = 0.0;
//cout<<"Called Potent"<<"\n";
    vector<double>tbpar(6);
    for (int nat=0; nat<16; nat++)
    {
        for (int nbt=26; nbt<30; nbt++)
        {
            //cout<<nat<<"\t"<<nbt<<"\n";
            if(nbt==26)
            {
                copy (params.begin(), params.begin()+6, tbpar.begin() );
                potent=potent+twobody(tbpar,Mol.Atoms[nat],Mol.Atoms[nbt]);
            }
            else if(nbt==27)
            {
                copy (params.begin()+6, params.begin()+12, tbpar.begin() );
                potent=potent+twobody(tbpar,Mol.Atoms[nat],Mol.Atoms[nbt]);
            }
            else if(nbt==28 || nbt==29)
            {
                copy (params.begin()+12, params.begin()+18, tbpar.begin() );
                potent=potent+twobody(tbpar,Mol.Atoms[nat],Mol.Atoms[nbt]);
            }
        }

    }
    return potent;
}
double twobody(vector<double>&pars,Atom& A,Atom& B)
{
    double tbd;
    double dist=atom_dist(A,B);
    tbd=pars[0]*exp(-pars[1]*dist)+(pars[2]*pow(dist,-2))+(pars[3]*pow(dist,-4))+(pars[4]*pow(dist,-6))+(pars[5]*pow(dist,-8));
    return tbd;
}
void print_vec(vector<double>& tes)
{
    for (unsigned int n=0; n<tes.size(); n++) cout<<n<<"\t"<<tes[n]<<endl;
}
