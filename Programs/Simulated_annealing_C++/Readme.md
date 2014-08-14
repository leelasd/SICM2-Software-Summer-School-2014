Fitting with Simulated Annealing Procedure
=========
## *Leela S. Dodda*
### Yale University
#### leela.dodda@yale.edu


This code is written with respect to specific system i.e Formaldehyde - Pyrene system. It is used for obtaining parameters to represent the non-covalent interactions between two molecules. This code is written in C++ and effeorts are underway to parallelize using OpenMP

This Program takes the interaction energies, corresponding geometries and initial parameters as their input and gives the fitted values of parameters as the output. 

List of Files
----
###include/ 
* **molecule.h**  : New data types such as Molecule type and Atom type are declared here and also the list of all functions in this program. 
 
 
### src/


* **main.cpp**              : Need to change this to setup Temperature and number of moves at each temperature 
 
* **damping_fn.cpp**        : Damping functions are used along with Buckingham potential to represent two-body potential
 
* **mc_fitting.cpp**        : Simulated Annealing Algorithm and Chisq function
 
* **molecule_funcs.cpp**    : Functions to read and write Molecule and Atom type data are defined here 
 
* **potential_funcs.cpp**   : Buckingham potential to represent 2-body potential are defined here

### makefile
* Made using *Make* functionality available in Linux  easy to understand and anyone can edit


Requirements
-----------
* g++ with C++11 standard 
* make 

License
-----------
MIT

