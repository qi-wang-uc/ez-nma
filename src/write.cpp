#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <utility>
#include "../include/nma.h"
#include "../include/write.h"


void write_data(const std::string& out_name, integer& n_modes, const real& tol) {
    std::cout << "WriteData> Wrting data to " << out_name << std::endl;
	const auto dim3 = 3*get_natom();
	if (n_modes > dim3) {
		std::cout << "WriteData> Requested output more than maximum eigenvectors." 
				  << std::endl
				  << "WriteData> Writing all the eigenvectors ..." << std::endl;
		n_modes = dim3;
	}
	/* Sort the eigenvalues */
	std::vector<std::pair<integer, real> > X; 
	for(integer i=0; i<dim3; ++i) {
		real val = get_H_elem(i*dim3+i) > tol ? get_H_elem(i*dim3+i) : 0.0;
		X.push_back(std::make_pair(i, val));
	}
	auto cmp_lt = [](const std::pair<integer, real>& a,
					 const std::pair<integer, real>& b) { 
					 return a.second < b.second;};
	std::sort(X.begin(), X.end(), cmp_lt);
	/* 3. write data into output file */
	auto width_num = std::to_string(n_modes).size()+1;
	std::ofstream out_file(out_name);
	for(integer i=0; i<n_modes; ++i) {
		// Write mode # and eigenvalue
		out_file << std::right << std::setw(width_num) << i+1 << " " 
				 << std::left << std::setprecision(6) << std::setiosflags(std::ios::fixed) 
				 << X[i].second << " : ";
		// Write eigenvectors.
		for(integer j=0; j<dim3; ++j) {
			out_file << get_E_elem(dim3*j+X[i].first) << " ";
		}
		out_file << std::endl;
	}
	out_file.close();
	std::cout << "WriteData> Done." << std::endl;
}
