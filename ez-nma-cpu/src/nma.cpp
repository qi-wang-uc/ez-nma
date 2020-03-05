#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <<numeric>
#include <omp.h>
#include "../include/main.h"
#include "../include/read.h"
#include "../include/util.h"

static std::vector<Coor> nma_coor;	// coordinates to construct H (3N)
static std::vector<Coor> ref_coor;	// coordintes to calculate overlap (3N)
static std::vector<real> H;	// Hessian matrix (3N x 3N);
static std::vector<real> E;	// Eigenvectors of H(3N x 3N);
/* This will decrease the complexity, but the running time actually increased. */
// static std::vector<real> P;	// Temporary matrix (3N x 3N) to hold off-diagonal elements;
// static std::vector<real> Q;	// Temporary matrix (3N x 3N) to transform H to P;

void build_hessian(const std::string& inp_name, const real& r_cutoff) {
    /* 1. read coordinates */
    if(!read_coor(inp_name, nma_coor)) return;
    integer natom = nma_coor.size();
    integer LD = 3 * natom;
    std::cout << "BuildHessian> Constructing Hessian Matrix of (" 
              << LD * LD << ") elements ..." << std::endl;
    H.resize(LD*LD); std::fill(H.begin(), H.end(), 0.0);
    // P.resize(LD*LD); std::fill(P.begin(), P.end(), 0.0);
    // Q.resize(LD*LD); std::fill(Q.begin(), Q.end(), 1.0); 
    // for(integer i=0; i<LD; i++)
    // 	Q[i*LD+i] = 0.0;
    /* 2. build Hessian array, write VMD file if req'd */
    std::string vmd_name = std::string("contact-map-cutoff") + real2str(r_cutoff) + std::string(".tcl");
    std::ofstream vmd_file(vmd_name);
    std::cout << "BuildHessian> Building off-diagonal blocks ..." << std::endl;
    integer pair_counter = 0;
    real cutoff2 = r_cutoff * r_cutoff;
    for(integer i=0; i<LD; i+=3) {
        auto a1 = i/3; auto coor1 = nma_coor[a1];
        for(integer j=0; j<LD; j+=3) {
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
                    H[(i+0)*LD + j+0] = (coor1.x-coor2.x)*(coor1.x-coor2.x)*dist2;
                    H[(i+0)*LD + j+1] = (coor1.x-coor2.x)*(coor1.y-coor2.y)*dist2;
                    H[(i+0)*LD + j+2] = (coor1.x-coor2.x)*(coor1.z-coor2.z)*dist2;
                    H[(i+1)*LD + j+0] = (coor1.y-coor2.y)*(coor1.x-coor2.x)*dist2;
                    H[(i+1)*LD + j+1] = (coor1.y-coor2.y)*(coor1.y-coor2.y)*dist2;
                    H[(i+1)*LD + j+2] = (coor1.y-coor2.y)*(coor1.z-coor2.z)*dist2;
                    H[(i+2)*LD + j+0] = (coor1.z-coor2.z)*(coor1.x-coor2.x)*dist2;
                    H[(i+2)*LD + j+1] = (coor1.z-coor2.z)*(coor1.y-coor2.y)*dist2;
                    H[(i+2)*LD + j+2] = (coor1.z-coor2.z)*(coor1.z-coor2.z)*dist2;
                }
            }			
        }
    }
    if(vmd_file.is_open()) vmd_file.close();
    /* calculate on-diagonal elements */
    for(integer i=0; i<LD; i+=3) {
        for(integer j=0; j<LD; j+=3) {
            if(i!=j){
                H[(i+0)*LD + i+0] -= H[(i+0)*LD + j+0];
                H[(i+0)*LD + i+1] -= H[(i+0)*LD + j+1];
                H[(i+0)*LD + i+2] -= H[(i+0)*LD + j+2];
                H[(i+1)*LD + i+0] -= H[(i+1)*LD + j+0];
                H[(i+1)*LD + i+1] -= H[(i+1)*LD + j+1];
                H[(i+1)*LD + i+2] -= H[(i+1)*LD + j+2];
                H[(i+2)*LD + i+0] -= H[(i+2)*LD + j+0];
                H[(i+2)*LD + i+1] -= H[(i+2)*LD + j+1];
                H[(i+2)*LD + i+2] -= H[(i+2)*LD + j+2];
            }
        }
    }
    pair_counter /= 2;
    std::cout <<"BuildHessian> Done. (" << pair_counter << ") pairs of contacts found." << std::endl;
}

void diag_hessian(const real& tol) {
    integer natom = nma_coor.size();
    integer LD = 3* natom; 
    auto fdim = static_cast<real>(LD*LD);
    std::cout << "DiagHessian> Diagonalizing Hessian matrix ..." << std::endl;
    E.resize(LD*LD);
    std::fill(E.begin(), E.end(), 0.0);
    for(integer i=0; i<LD; i++)
        E[i*LD+i] = 1.0;
    real sum_offd = 0.0;
    for(integer i=0; i<LD; i++) {		
        for(integer j=0; j<LD; j++) {
            if(i!=j) sum_offd += H[i*LD+j]*H[i*LD+j];
        }
    }
    if(sum_offd < tol) return;
    real avg_offd = 0.5*sum_offd/fdim;
    integer itr_counter = 0;
    integer sweep_counter = 1;
    while(sum_offd > tol) {
        itr_counter++;
        std::cerr << std::fixed
                  <<"DiagHessian> Sweep ("<< std::setw(8) << itr_counter <<") | "
                  << "Converging index: " << std::setw(12) << std::setprecision(6) << sum_offd << " | "
                  << "Target value: " << tol << "\n";
                  // << std::endl;
        for(integer i=0; i<LD-1; i++) {
            for(integer j=i+1; j<LD; j++) {
                if(H[j*LD+i]*H[j*LD+i] < avg_offd) continue;
                // std::transform(H.cbegin(), H.cend(), Q.cbegin(), P.begin(), std::multiplies<real>());
                // sum_offd = std::inner_product(P.cbegin(), P.cend(), P.cbegin(), 0.0);
                sum_offd = std::inner_product(H.cbegin(), H.cend(), H.cbegin(), 0.0);
                for (integer i=0; i<LD; i++) {
                    sum_offd -= H[i*LD+i]*H[i*LD+i];
                }
                avg_offd = sum_offd/fdim;
// step3. Calculate coefficients [c] and [s] for Givens matrix.
                real beta = (H[j*LD+j]-H[i*LD+i])/(2.0*H[j*LD+i]);
                real coeff = 0.5*beta/sqrt(1.0+beta*beta);
                real s = sqrt(std::max(0.5+coeff, 0.0));
                real c = sqrt(std::max(0.5-coeff, 0.0));
// step4. Update rows [i] and [j] of Hessian matrix (Givens matrix pre-multiplies Hessian).
            //	#pragma omp parallel for
                for(integer k=0; k<LD; k++) {
                    real cs =  c*H[i*LD+k] + s*H[j*LD+k];
                    real sc = -s*H[i*LD+k] + c*H[j*LD+k];
                    H[i*LD+k] = cs;
                    H[j*LD+k] = sc;
                }
// step5. Givens matrix post-multiplies Hessian. Also calculate eigenvectors. 
            //	#pragma omp parallel for
                for(integer k=0; k<LD; k++) {
                    real cs =  c*H[k*LD+i] + s*H[k*LD+j];
                    real sc = -s*H[k*LD+i] + c*H[k*LD+j];
                    H[k*LD+i] = cs;
                    H[k*LD+j] = sc;
                    cs =  c*E[k*LD+i] + s*E[k*LD+j];
                    sc = -s*E[k*LD+i] + c*E[k*LD+j];
                    E[k*LD+i] = cs;
                    E[k*LD+j] = sc;
                }
            }
        }
        sweep_counter++;
    }
    std::cout << "DiagHessian> Done." << std::endl;
}

void calc_overlap(const std::string& job_name, const std::string& ref_name) {
    if(ref_name.empty()) return;
}

const real& get_H_elem(integer query) {
    return H[query];
}

const real& get_E_elem(integer query) {
    return E[query];
}

integer get_natom () {
    return nma_coor.size();
}
