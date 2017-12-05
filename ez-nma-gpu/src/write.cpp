#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <utility>
#include "../include/write.h"

void write_data(Leading_Dim leading_dim, integer n_modes, const std::string& out_name, const real tol) {
    std::cout << "WriteData> Wrting data to " << out_name << std::endl;
	if (n_modes > leading_dim.LD_data) {
		std::cout << "WriteData> Requested output more than maximum eigenvectors." 
				  << std::endl
				  << "WriteData> Writing all the eigenvectors ..." << std::endl;
		n_modes = leading_dim.LD_data;
    }
    
	/* Sort the eigenvalues in ascending order */
	std::vector<std::pair<integer, real> > X; 
	for(integer i=0; i<leading_dim.LD_matrix; ++i) {
		real val = h_H[i*leading_dim.LD_matrix+i] > tol ? h_H[i*leading_dim.LD_matrix+i] : tol;
		val = h_H[i*leading_dim.LD_matrix+i];
		X.push_back(std::make_pair(i, val));
	}
	auto cmp_lt = [](const std::pair<integer, real>& a,
					 const std::pair<integer, real>& b) { 
					 return a.second < b.second;};
    std::sort(X.begin(), X.end(), cmp_lt);
    
	/* Write data into output file */
	auto width_num = std::to_string(n_modes).size()+1;
	std::ofstream out_file(out_name);
	for(integer mode_counter=0, i=leading_dim.LD_diff(); mode_counter<n_modes; mode_counter++,i++) {
		// Write mode # and eigenvalue
		out_file << std::right << std::setw(width_num) << mode_counter+1 << " " 
				 << std::left << std::setprecision(6) << std::setiosflags(std::ios::fixed) 
				 << X[i].second << " : ";
		// Write eigenvectors.
		for(integer j=0; j<leading_dim.LD_data; j++) {
			out_file << h_E[leading_dim.LD_matrix*X[i].first +j] << " ";
		}
		out_file << std::endl;
	}
	out_file.close();
	std::cout << "WriteData> Done." << std::endl << std::endl;
}
