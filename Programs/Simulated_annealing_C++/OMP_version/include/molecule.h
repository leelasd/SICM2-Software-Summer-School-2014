#ifndef MOLECULE_H_INCLUDED
#define MOLECULE_H_INCLUDED
#include <vector>
using namespace std;
struct Atom
{
    char   At;
    double x,y,z;
};
struct Molecule
{
    int Nat;
    vector<Atom> Atoms;
};
void Read_file(vector<Molecule>&,int,int);
void Read_data(vector<double>&,vector<double>&,vector<double>&,vector<double>&,vector<double>&,int);
void sim_anneal(vector<Molecule>&,vector<double>&,vector<double>&,int,double);
void Read_params(vector<double>&,int);
void read_atom(Atom&);
void print_atom(Atom&);
void Read_Mol(Molecule&);
void Print_Mol(Molecule&);
double atom_dist(Atom&,Atom&);
void Print_Mol_Inter(Molecule&,int,int);
double potent(Molecule&,vector<double>&);
double potent_df(Molecule&,vector<double>&);
double twobody(vector<double>&,Atom&,Atom&);
double tbd_damping(vector<double>&,Atom&,Atom&);
double totpot(vector<Molecule>&,int,vector<double>&);
void print_vec(vector<double>&);
#endif // MOLECULE_H_INCLUDED
