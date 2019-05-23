#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <thrust/fill.h>
#include "../include/read.h"
#include "../include/util.h"
#include "../include/build.h"

std::vector<Coor> nma_coor;     // coordinates to construct H (3N)
std::vector<Coor> ref_coor;     // coordintes to calculate overlap (3N)
thrust::host_vector<real> h_H;  // Hessian matrix (3N x 3N);
thrust::host_vector<real> h_E;  // Eigenvectors matrix (3N x 3N);

/* variable LD is abbr. for Leading Dimension of the matrix, which is 3N */
Leading_Dim build_hessian(const std::string& inp_name, const real& r_cutoff) {
    /* 1. read coordinates */
    auto result = Leading_Dim(0, 0);
    if(!read_coor(inp_name, nma_coor)) return result;
    integer natom = nma_coor.size();
    integer LD_data = 3*natom;
    integer LD_matrix = (0==LD_data%2) ? LD_data : (LD_data+1);
    std::cout << "BuildHessian> Constructing Hessian Matrix of ("
              << LD_matrix * LD_matrix << ") elements with leading dimension ("
              << LD_matrix << ") ..." << std::endl;
    std::cout << "BuildHessian> *Note: If Hessian matrix had an odd leading dimension,"
              << " one more row and column with zeroes will be padded." << std::endl;
    /* init hessian array */
    h_H.resize(LD_matrix*LD_matrix);
    std::fill(h_H.begin(), h_H.end(), 0.0);
    h_E.resize(LD_matrix*LD_matrix);
    std::fill(h_E.begin(), h_E.end(), 0.0);
    for(integer i=0; i<LD_data; i++) h_E[i*LD_matrix+i] = 1.0;
    /* 2. build Hessian array */
    std::cout << "BuildHessian> Building off-diagonal blocks ..." << std::endl;
    integer pair_counter = 0;
    real cutoff2 = r_cutoff * r_cutoff;
    for(integer i=0; i<LD_data; i+=3) {
        auto a1 = i/3; auto coor1 = nma_coor[a1];
        for(integer j=0; j<LD_data; j+=3) {
            auto a2 = j/3; auto coor2 = nma_coor[a2];
            if(i!=j) {
                real dist2 = (coor1 - coor2).norm2();
                if(dist2 < cutoff2) {
                    /* calculate off-diagonal elements */
                    pair_counter++;
                    dist2 = -1.0/dist2;
                    h_H[(i+0)*LD_matrix + j+0] = (coor1.x-coor2.x)*(coor1.x-coor2.x)*dist2;
                    h_H[(i+0)*LD_matrix + j+1] = (coor1.x-coor2.x)*(coor1.y-coor2.y)*dist2;
                    h_H[(i+0)*LD_matrix + j+2] = (coor1.x-coor2.x)*(coor1.z-coor2.z)*dist2;
                    h_H[(i+1)*LD_matrix + j+0] = (coor1.y-coor2.y)*(coor1.x-coor2.x)*dist2;
                    h_H[(i+1)*LD_matrix + j+1] = (coor1.y-coor2.y)*(coor1.y-coor2.y)*dist2;
                    h_H[(i+1)*LD_matrix + j+2] = (coor1.y-coor2.y)*(coor1.z-coor2.z)*dist2;
                    h_H[(i+2)*LD_matrix + j+0] = (coor1.z-coor2.z)*(coor1.x-coor2.x)*dist2;
                    h_H[(i+2)*LD_matrix + j+1] = (coor1.z-coor2.z)*(coor1.y-coor2.y)*dist2;
                    h_H[(i+2)*LD_matrix + j+2] = (coor1.z-coor2.z)*(coor1.z-coor2.z)*dist2;
                }
            }
        }
    }
    /* calculate on-diagonal elements */
    std::cout << "BuildHessian> Building diagonal blocks ..." << std::endl;
    for(integer i=0; i<LD_data; i+=3) {
        for(integer j=0; j<LD_data; j+=3) {
            if(i!=j){
                h_H[(i+0)*LD_matrix + i+0] -= h_H[(i+0)*LD_matrix + j+0];
                h_H[(i+0)*LD_matrix + i+1] -= h_H[(i+0)*LD_matrix + j+1];
                h_H[(i+0)*LD_matrix + i+2] -= h_H[(i+0)*LD_matrix + j+2];
                h_H[(i+1)*LD_matrix + i+0] -= h_H[(i+1)*LD_matrix + j+0];
                h_H[(i+1)*LD_matrix + i+1] -= h_H[(i+1)*LD_matrix + j+1];
                h_H[(i+1)*LD_matrix + i+2] -= h_H[(i+1)*LD_matrix + j+2];
                h_H[(i+2)*LD_matrix + i+0] -= h_H[(i+2)*LD_matrix + j+0];
                h_H[(i+2)*LD_matrix + i+1] -= h_H[(i+2)*LD_matrix + j+1];
                h_H[(i+2)*LD_matrix + i+2] -= h_H[(i+2)*LD_matrix + j+2];
            }
        }
    }
    pair_counter /= 2;
    std::cout <<"BuildHessian> Done. (" << pair_counter << ") pairs of contacts found." 
              << std::endl << std::endl;
    result.LD_data = LD_data;
    result.LD_matrix = LD_matrix;
    return result;
}
