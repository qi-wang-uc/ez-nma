#ifndef NMA_H
#define NMA_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/system/cuda/experimental/pinned_allocator.h>
#include <string>
#include <vector>
#include "main.h"

struct Leading_Dim {
    integer LD_data;    // Raw leading dimension of Hessian, odd or even number.
    integer LD_matrix;  // Actual leading dimension of Hessian after building, must be even number for parallel Jacobi.
    Leading_Dim() {}
    Leading_Dim(integer a, integer b) {
        LD_data   = a;
        LD_matrix = b;
    }
    integer LD_diff() {
        return LD_matrix-LD_data;
    }
};

Leading_Dim build_hessian(const std::string& inp_name, const real& r_cutoff);

extern std::vector<Coor> nma_coor;	    // coordinates to construct H (3N)
extern std::vector<Coor> ref_coor;	    // coordintes to calculate overlap (3N)
extern thrust::host_vector<real> h_H;	// Hessian matrix (3N x 3N);
extern thrust::host_vector<real> h_E;	// Eigenvectors matrix (3N x 3N);

#endif