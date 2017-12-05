#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include "main.h"

struct Config {
	std::string job_name = "ez-nma";
	std::string nma_coor = "";
    std::string ref_coor = "";
    integer n_modes = 10;
	real r_cutoff = 10.0;
    real tol = 1e-6;
};

bool read_config(const std::string& inp_name, Config& config);

bool is_good_config(const Config& config);

void print_config(const Config&);

bool read_coor(const std::string& pdb_name, std::vector<Coor>& coor); 

/*************** sample configuration file ********************
# Configuration file for EZNMA, comments lines start with "#"
      
#1 [job_name] as name prefix for output files.
# not mandatory, default is <ez-nma>
job_name	clpp-14mer

#2 [nma_coor] specifies coordinates perform normal mode analysis
# madatory, need this to build Hessian array.
#nma_coor	clpp-monomer-newstate.charmm.pdb
#nma_coor	clpp-allatom-zincbinding.pdb
nma_coor	spastin-hexamer-without-mtbd-ordered.pdb

#3 [r_cutoff] specifies cutoff distance to build spring network.
# For GNM, r_cutoff is recommended to be set between 7 and 8A (Kundu 2002);
# For ENM, r_cutoff is recommended to be set between 13 to 15A (Eyal 2006).
# not mandatory, default is 10.0
r_cutoff	10.0

#4 [n_modes] specifies numer of normal modes with lowest 
# frequencies (energies) needed to write.
# not mandatory, default is 10, maxium is 3*Natom for
# ENM and Natom for GNM.
n_modes		10

#5 [ref_coor] specifies coordinates file to calculate overlap.
# not required if overlap calculation is not needed. (currently not supported.)
# ref_coor referce_coordiante.pdb

#6 [tol] specifies the error tolerance to diagonalize Hessian
# matrix, default is 1e-6. It's not recommended for user to
# modify this value thus not shown here.
**************************************************************/

#endif