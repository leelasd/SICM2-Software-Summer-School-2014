#include<iostream>
#include<cmath>
#include <fstream>
#include "molecule.h"
using namespace std;
void Read_params(vector<double>& pars,int npar)
{
    for (int l=0; l<npar; l++)
    {
        double part;
        cin>>part;
        pars.push_back(part);
    }
}
void Read_data(vector<double>& rcc,vector<double>& rco,vector<double>& rcha,vector<double>& rchb,vector<double>& ydat,int nat)
{
    ifstream myfile ("ydat");
    if (myfile.is_open())
    {
        for (int l=0; l<nat; l++)
        {
            double p1,p2,p3,p4,p5;
            myfile>>p1>>p2>>p3>>p4>>p5;
            rcc.push_back(p1);
            rco.push_back(p2);
            rcha.push_back(p3);
            rchb.push_back(p4);
            ydat.push_back(p5);
        }
        myfile.close();
    }
    // cout<<"DATA READING IS DONE"<<endl;
}
void Read_file(vector<Molecule>& formal,int nat,int nmol)
{
    formal.resize(nmol);
    for (int k=0; k<nmol; k++)
    {
        formal[k].Nat=nat;
        formal[k].Atoms.resize(nat);
        Read_Mol(formal[k]);
    }
}
void Print_Mol_Inter(Molecule& Mol,int Na,int Nb)
{
    for (int i=0; i<Na; i++)
    {
        for (int k=Na; k<Mol.Nat; k++)
        {
            cout<<i+1<<k+1<<"\t"<<atom_dist(Mol.Atoms[i],Mol.Atoms[k])<<"\t";
        }
        cout <<"\t"<<endl;
    }
}
double atom_dist(Atom& A,Atom& B)
{
    double dist = sqrt(pow(A.x-B.x,2)+pow(A.y-B.y,2)+pow(A.z-B.z,2));
    return dist;
}
void Read_Mol(Molecule& Mol)
{
    for (int i=0; i<Mol.Nat; i++)
    {
        read_atom(Mol.Atoms[i]);
    }
}
void Print_Mol(Molecule& Mol)
{
    for (int i=0; i<Mol.Nat; i++)
    {
        print_atom(Mol.Atoms[i]);
    }
}

void read_atom(Atom& Atoms)
{
    cin>>Atoms.At>>Atoms.x>>Atoms.y>>Atoms.z;
}
void print_atom(Atom& Atoms)
{
    cout<<Atoms.At<<"\t"<<Atoms.x<<"\t"<<Atoms.y<<"\t"<<Atoms.z<<endl;
}
