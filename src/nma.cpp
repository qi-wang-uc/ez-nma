#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <omp.h>
#include "../include/main.h"
#include "../include/read.h"
#include "../include/util.h"

static std::vector<Coor> nma_coor;	// coordinates to construct H (3N)
static std::vector<Coor> ref_coor;	// coordintes to calculate overlap (3N)
static std::vector<real> H;	// Hessian matrix (3N x 3N);
static std::vector<real> E;	// Eigenvectors of H(3N x 3N);

void build_hessian(const std::string& inp_name, const real& r_cutoff) {
	/* 1. read coordinates */
	if(!read_coor(inp_name, nma_coor)) return;
	unsigned int natom = nma_coor.size();
	unsigned int dim3 = 3 * natom;
	std::cout << "BuildHessian> Constructing Hessian Matrix of (" 
			  << dim3 * dim3 << ") elements ..." << std::endl;
	H.resize(dim3*dim3);
	std::fill(H.begin(), H.end(), 0.0);
	/* 2. build Hessian array, write VMD file if req'd */
	std::string vmd_name = std::string("contact-map-cutoff") + real2str(r_cutoff) + std::string(".tcl");
	std::ofstream vmd_file(vmd_name);
	std::cout << "BuildHessian> Building off-diagonal blocks ..." << std::endl;
	unsigned int pair_counter = 0;
	real cutoff2 = r_cutoff * r_cutoff;
	for(unsigned int i=0; i<dim3; i+=3) {
		auto a1 = i/3; auto coor1 = nma_coor[a1];
		for(unsigned int j=0; j<dim3; j+=3) {
			auto a2 = j/3; auto coor2 = nma_coor[a2];
			if(i!=j) {
				real dist2 = (coor1 - coor2).norm2();
				if(dist2 < cutoff2) {
					/* write contact map */
					if(i<j && vmd_file.is_open()) {
						vmd_file << "draw line { "
							 << coor1.x << " " << coor1.y << " " << coor1.z << " } { "
							 << coor2.x << " " << coor2.y << " " << coor2.z << " } width 3\n";
					}
					/* calculate off-diagonal elements */
					pair_counter++;
					dist2 = -1.0/dist2;
					H[(i+0)*dim3 + j+0] = (coor1.x-coor2.x)*(coor1.x-coor2.x)*dist2;
					H[(i+0)*dim3 + j+1] = (coor1.x-coor2.x)*(coor1.y-coor2.y)*dist2;
					H[(i+0)*dim3 + j+2] = (coor1.x-coor2.x)*(coor1.z-coor2.z)*dist2;
					H[(i+1)*dim3 + j+0] = (coor1.y-coor2.y)*(coor1.x-coor2.x)*dist2;
					H[(i+1)*dim3 + j+1] = (coor1.y-coor2.y)*(coor1.y-coor2.y)*dist2;
					H[(i+1)*dim3 + j+2] = (coor1.y-coor2.y)*(coor1.z-coor2.z)*dist2;
					H[(i+2)*dim3 + j+0] = (coor1.z-coor2.z)*(coor1.x-coor2.x)*dist2;
					H[(i+2)*dim3 + j+1] = (coor1.z-coor2.z)*(coor1.y-coor2.y)*dist2;
					H[(i+2)*dim3 + j+2] = (coor1.z-coor2.z)*(coor1.z-coor2.z)*dist2;
				}
			}			
		}
	}
	if(vmd_file.is_open()) vmd_file.close();
	/* calculate on-diagonal elements */
	for(unsigned int i=0; i<dim3; i+=3) {
		for(unsigned int j=0; j<dim3; j+=3) {
			if(i!=j){
				H[(i+0)*dim3 + i+0] -= H[(i+0)*dim3 + j+0];
				H[(i+0)*dim3 + i+1] -= H[(i+0)*dim3 + j+1];
				H[(i+0)*dim3 + i+2] -= H[(i+0)*dim3 + j+2];
				H[(i+1)*dim3 + i+0] -= H[(i+1)*dim3 + j+0];
				H[(i+1)*dim3 + i+1] -= H[(i+1)*dim3 + j+1];
				H[(i+1)*dim3 + i+2] -= H[(i+1)*dim3 + j+2];
				H[(i+2)*dim3 + i+0] -= H[(i+2)*dim3 + j+0];
				H[(i+2)*dim3 + i+1] -= H[(i+2)*dim3 + j+1];
				H[(i+2)*dim3 + i+2] -= H[(i+2)*dim3 + j+2];
			}
		}
	}
	pair_counter /= 2;
	std::cout <<"BuildHessian> Done. (" << pair_counter << ") pairs of contacts found." << std::endl;
}

void diag_hessian(const real& tol) {
	unsigned int natom = nma_coor.size();
	unsigned int dim3 = 3* natom; 
	auto fdim = static_cast<real>(dim3*dim3);
	std::cout << "DiagHessian> Diagonalizing Hessian matrix ..." << std::endl;
	E.resize(dim3*dim3);
	std::fill(E.begin(), E.end(), 0.0);
	for(unsigned int i=0; i<dim3; i++)
		E[i*dim3+i] = 1.0;
	real sum_offd = 0.0;
	for(unsigned int i=0; i<dim3; i++) {		
		for(unsigned int j=0; j<dim3; j++) {
			if(i!=j) sum_offd += H[i*dim3+j]*H[i*dim3+j];
		}
	}
	if(sum_offd < tol) return;
	real avg_offd = 0.5*sum_offd/fdim;
	unsigned int itr_counter = 0;
	while(sum_offd > tol) {
		itr_counter++;
		std::cerr << std::fixed
				  <<"DiagHessian> Iteration ("<< std::setw(8) << itr_counter <<") | "
				  << "Converging index: " << std::setw(12) << std::setprecision(6) << sum_offd << " | "
				  << "Target value: " << tol << "\n";
				  // << std::endl;
		for(unsigned int i=0; i<dim3-1; i++) {
			for(unsigned int j=i+1; j<dim3; j++) {
				if(H[j*dim3+i]*H[j*dim3+i] < avg_offd) continue;
				sum_offd -= 2.0*H[j*dim3+i]*H[j*dim3+i];
				avg_offd = 0.5*sum_offd/fdim;
// step3. Calculate coefficients [c] and [s] for Givens matrix.
				real beta = (H[j*dim3+j]-H[i*dim3+i])/(2.0*H[j*dim3+i]);
				real coeff = 0.5*beta/sqrt(1.0+beta*beta);
				real s = sqrt(std::max(0.5+coeff, 0.0));
				real c = sqrt(std::max(0.5-coeff, 0.0));
// step4. Update rows [i] and [j] of Hessian matrix (Givens matrix pre-multiplies Hessian).
			//	#pragma omp parallel for
				for(unsigned int k=0; k<dim3; k++) {
					real cs =  c*H[i*dim3+k] + s*H[j*dim3+k];
					real sc = -s*H[i*dim3+k] + c*H[j*dim3+k];
					H[i*dim3+k] = cs;
					H[j*dim3+k] = sc;
				}
// step5. Givens matrix post-multiplies Hessian. Also calculate eigenvectors. 
			//	#pragma omp parallel for
				for(unsigned int k=0; k<dim3; k++) {
					real cs =  c*H[k*dim3+i] + s*H[k*dim3+j];
					real sc = -s*H[k*dim3+i] + c*H[k*dim3+j];
					H[k*dim3+i] = cs;
					H[k*dim3+j] = sc;
					cs =  c*E[k*dim3+i] + s*E[k*dim3+j];
					sc = -s*E[k*dim3+i] + c*E[k*dim3+j];
					E[k*dim3+i] = cs;
					E[k*dim3+j] = sc;
				}
			}
		}
	}
	std::cout << "DiagHessian> Done." << std::endl;
}

void calc_overlap(const std::string& job_name, const std::string& ref_name) {
	if(ref_name.empty()) return;
}

const real& get_H_elem(unsigned int query) {
	return H[query];
}

const real& get_E_elem(unsigned int query) {
	return E[query];
}

unsigned int get_natom () {
	return nma_coor.size();
}

