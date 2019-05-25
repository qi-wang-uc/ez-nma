#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include "main.h"

struct Config {
    std::string job_name = "ez-nma";
    std::string nma_coor = "";
    std::string ref_coor = "";
    std::string vmd_file = "";
    integer n_modes = 10;
    real r_cutoff = 10.0;
    real tol = 1e-6;
};

bool read_config(const std::string& inp_name, Config& config);

bool is_good_config(const Config& config);

void print_config(const Config&);

bool read_coor(const std::string& pdb_name, std::vector<Coor>& coor); 

#endif