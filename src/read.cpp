#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include "../include/read.h"

bool read_config(const std::string& inp_name, Config& config) {
    std::ifstream inp_file(inp_name);
    /* check if input file exists. */
    if(!inp_file.is_open()) {
        std::cout << "ERROR> Cannot open file [" << inp_name << "]." << std::endl;
        return false;
    } else {
        std::cout << "ReadConfig> Reading input from file [" << inp_name << "]" << std::endl;
    }
    /* now read configuration */
    std::string each_line;
    std::stringstream each_stream;
    std::string buff;
    auto is_space = [](char c) {return std::isspace(c);};
    while(getline(inp_file, each_line)) {
        /* filter out comment lines, empty lines and lines of spaces */
        if(0==each_line.substr(0,1).compare("#") ||
            std::all_of(each_line.begin(),each_line.end(),is_space))
            continue;
        each_stream.str(each_line);
        each_stream >> buff;
        if(0==buff.compare("job_name")) {
            each_stream >> config.job_name;
        } else if (0==buff.compare("nma_coor")) {
            each_stream >> config.nma_coor;
        } else if (0==buff.compare("ref_coor")) {
            each_stream >> config.ref_coor;
        } else if (0==buff.compare("r_cutoff")) {
            each_stream >> config.r_cutoff;
        } else if (0==buff.compare("n_modes") ) {
            each_stream >> config.n_modes;
        } else if (0==buff.compare("vmd_file")) {
            each_stream >> config.vmd_file;
        } else {
            std::cout << "ERROR> Unrecognized parameter: "
                    << buff << std::endl;
        }
        each_stream.clear();
    }
    std::cout << "ReadConfig> Done." << std::endl;
    return true;
}

void print_config(const Config& config) {
    std::cout << "After reading, the following parameters will be used:" << std::endl;
    std::cout << "********************************" << std::endl
              << "job_name is: " << config.job_name << std::endl
              << "nma_coor is: " << config.nma_coor << std::endl
              << "ref_coor is: " << config.ref_coor << std::endl
              << "vmd_file is: " << config.vmd_file << std::endl
              << "n_modes  is: " << config.n_modes  << std::endl
              << "r_cutoff is: " << config.r_cutoff << std::endl
              << "********************************" << std::endl;
}

bool read_coor(const std::string& pdb_name, std::vector<Coor>& coor) {
	std::cout << "ReadCoor> Reading coordinates from file [" << pdb_name
	          << "]" << std::endl;
	std::ifstream inp_file(pdb_name);
	if(!inp_file.is_open()) {
        std::cout << "ERROR> Coordinate file not found!" << std::endl;
        return false;
    }
	std::string each_line;
	unsigned int i_beads = 0;
	while(std::getline(inp_file, each_line)){
		if(each_line.substr(0,4)!="ATOM" || each_line.substr(13,2)!="CA"){
			continue;
		}
        coor.push_back(Coor(std::stof(each_line.substr(31,8)),
                            std::stof(each_line.substr(39,8)),
                            std::stof(each_line.substr(47,8))));
		i_beads++;
	}
    inp_file.close();
	std::cout << "ReadCoor> After reading, (" << i_beads 
	          << ") CA atoms are recorded." << std::endl;
	if(i_beads<2) {
        std::cout << "ERROR> Not enough atoms for NMA!" << std::endl;
        return false;
    } 
    std::cout << "ReadCoor> Done." << std::endl;
    return true;
}