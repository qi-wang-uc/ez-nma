#include <iostream>
#include <string>
#include <ctime>
#include "../include/read.h"
#include "../include/nma.h"
#include "../include/write.h"
#include "../include/util.h"

int main(int argc, char* argv[]) {
    clock_t start_time = clock();
    std::ios_base::sync_with_stdio(false);
    std::cout << "EZNMA> ***********************************************" << std::endl;
    std::cout << "EZNMA> Normal Mode Analysis with Elastic Network Model" << std::endl
              << "EZNMA> version 1.0, Feburary 2017." << std::endl;
    std::cout << "EZNMA> ***********************************************" << std::endl;
    if(argc!=2) {
        std::cout << "ERROR> Missing input file!" << std::endl;
        return 1;
    }
    /* 1. read input */
    std::string inp_name = argv[1];
    Config config;
    if(!read_config(inp_name, config)) { return 1; }
    print_config(config);
    /* 2. run nma */
    build_hessian(config.nma_coor, config.r_cutoff);
    diag_hessian(config.tol);
    // calc_overlap(config.job_name, config.ref_coor);
    /* 3. write output */
    std::string out_name = config.job_name + std::string("_nma_data_cutoff") + 
                           real2str(config.r_cutoff) + std::string(".dat");
    write_data(out_name, config.n_modes, config.tol);
    std::cout << "EZNMA> Normal termination of the program." << std::endl;
    clock_t end_time = clock();
    time_stat(start_time, end_time);
    return 0;
}