#ifndef NMA_H
#define NMA_H

#include <string>
#include "main.h"


void build_hessian(const std::string& inp_name, const real& r_cutoff);

void diag_hessian(const real& tol);

void sort_eigen(const real& tol);

void calc_overlap(const std::string& job_name, const std::string& ref_name);

const real& get_H_elem(unsigned int query);

const real& get_E_elem(unsigned int query);

unsigned int get_natom ();

#endif