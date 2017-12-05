#include <iostream>
#include <iomanip>
#include <cuda.h>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/copy.h>
#include "../include/diag.h"
#include "../include/build.h"
#include "../include/main.h"

/* Calculate index for each processing unit using Chess Tournament (CT) algorithm. */
__device__ void calc_index_ij(integer bid, integer iter, integer LD, integer& i, integer& j) {
	integer LD_1 = LD - 1;
	integer index1 = (bid+iter) % LD_1;
	integer index2 = (0==bid) ? LD_1 : (LD_1-bid+iter) % LD_1;
	i = min(index1, index2);
	j = max(index1, index2);
}

/* Givens matrix pre-multiplies Hessian. */
__global__ void jacobi_kernel_step1 (real* H, real* M, real* S, real* C, integer LD, integer iter, real sweep_thresh, integer stride) {
	integer bid = blockIdx.x;
	real s = S[bid];
	real c = C[bid];
	integer i, j;
	calc_index_ij(bid, iter, LD, i, j);
	real H_ji = H[j*LD+i];
	for(integer tid=threadIdx.x; tid<LD; tid+=stride) {
		real H_ik = H[i*LD+tid];
		real H_jk = H[j*LD+tid];
		if(H_ji*H_ji < sweep_thresh) {
			// inefficient global memory access pattern. should be improved.
			M[tid*LD+i] = H_ik;
			M[tid*LD+j] = H_jk;
		} else {
			// Update rows [i] and [j] for Hessian matrix (pre-multiply).
			M[tid*LD+i] =  c*H_ik + s*H_jk;
			M[tid*LD+j] = -s*H_ik + c*H_jk;
		}
	}
}

/* Givens matrix post-multiplies Hessian. Also calculate eigenvectors. */
__global__ void jacobi_kernel_step2 (real* H, real* M, real* S, real* C, integer LD, integer iter, real sweep_thresh, integer stride, real* E) {
	integer bid = blockIdx.x;
	real s = S[bid];
	real c = C[bid];
	integer i, j;
	calc_index_ij(bid, iter, LD, i, j);
	real H_ji = H[j*LD+i];
	__syncthreads();
	for(integer tid=threadIdx.x; tid<LD; tid+=stride) {
		// M elements is stored column-wise to improve memory coalesce
		real M_ki = M[i*LD+tid];
		real M_kj = M[j*LD+tid];
		if(H_ji*H_ji < sweep_thresh) {
			// H is symmetric, so use H_ik(Hjk) instead of H_ki(H_kj) for memory coalesce
			H[i*LD+tid] = M_ki;
			H[j*LD+tid] = M_kj;
		} else {
			// H is symmetric, so use H_ik(H_jk) instead of H_ki(H_kj) for memory coalesce
			H[i*LD+tid] =  c*M_ki + s*M_kj;
			H[j*LD+tid] = -s*M_ki + c*M_kj;
			// Apply rotations to eigenvectors, storing in row-major.
			real E_ik = E[i*LD+tid];
			real E_jk = E[j*LD+tid];
			E[i*LD+tid] =  c*E_ik + s*E_jk;
			E[j*LD+tid] = -s*E_ik + c*E_jk;
		}
	}
}

/* Calculate S and C for next Givens rotation, quick and dirty */
__global__ void calc_param_kernel (real* H, real *S, real *C, integer LD, integer iter) {
	integer bid = blockIdx.x;
	integer i, j;
	calc_index_ij(bid, iter, LD, i, j);
	real beta = (H[j*LD+j]-H[i*LD+i]) / (2.0*H[j*LD+i]);
	real coeff = 0.5*beta / sqrt(1.0+beta*beta);
	real s = sqrt(fmax(0.5+coeff, 0.0));
	real c = sqrt(fmax(0.5-coeff, 0.0));
	S[bid] = s; 
	C[bid] = c;
}

void diag_hessian_gpu(integer LD, real tol) {
	/********************** MEMORY SPACE ALLOCATION **********************/
	const integer dev_id = 0;
	cudaDeviceProp dev_prop;
	cudaGetDeviceProperties(&dev_prop, dev_id);
	std::cout << "DiagHessian> Diagonalizing Hessian matrix on device [" << dev_prop.name 
			  << "] with computability (" << dev_prop.major << "." << dev_prop.minor << ")" 
			  << std::endl;
	std::cout << "DiagHessian> If encountering CUDA kernel launch failure, please change the computability in Makefile accordingly." << std::endl;
	// These 2 numbers will be used many times
	const integer LDSQ = LD*LD;
	const integer LDHF = LD/2;
	thrust::device_vector<real> d_H = h_H; // Hessian matrix
	thrust::device_vector<real> d_E = h_E; // Eigenvector matrix
	thrust::device_vector<real> d_M(LDSQ, 0.0); // (Auxiliary) interMediate matrix
	thrust::device_vector<real> d_S(LDHF, 0.0); // S elements array of Givens rotation matrix
	thrust::device_vector<real> d_C(LDHF, 0.0); // C elements array of Givens rotation matrix
	real* pd_H = thrust::raw_pointer_cast(d_H.data());
	real* pd_E = thrust::raw_pointer_cast(d_E.data());
	real* pd_M = thrust::raw_pointer_cast(d_M.data());
	real* pd_S = thrust::raw_pointer_cast(d_S.data());
	real* pd_C = thrust::raw_pointer_cast(d_C.data());

	/********************** KERNEL LAUNCHING CONFIGURATION **********************/
	const integer stride = dev_prop.maxThreadsPerBlock;
	// kernel launch parameters for calculating Givens rotation coefficients
	integer num_threads_param = 1;
	integer num_blocks_param  = LDHF;
	dim3 dim_block_param(num_threads_param, 1, 1);
	dim3 dim_grid_param(num_blocks_param, 1, 1);
	// kernel launch parameters for jacobi-sweep 
	integer num_threads_jacobi = stride;
	integer num_blocks_jacobi = LDHF;
	dim3 dim_block_jacobi(num_threads_jacobi, 1, 1);
	dim3 dim_grid_jacobi(num_blocks_jacobi, 1, 1);
	
	/********************** DATA INITIALIZATION **********************/
	// h_H and h_E are currently idle and can be resued as temporary data containers
	std::fill(h_E.begin(), h_E.end(), 1.0);
	for(integer i=0; i<LD; i++) h_E[i*LD+i] = 0.0;	// h_E will zero out diagonal elements of its element-wise multiplier.
	thrust::device_vector<real> d_N = h_E; // h_N as auxiliary matrix on device works like h_E.
	std::transform(h_H.cbegin(), h_H.cend(), h_E.cbegin(), h_H.begin(), std::multiplies<real>());
	real offd_sumsq = std::inner_product(h_H.cbegin(), h_H.cend(), h_H.cbegin(), 0.0);
	if (offd_sumsq < tol) return;
	std::cout << "DiagHessian> Converging index before Jacobi sweep: " << std::fixed << offd_sumsq << std::endl;
	real delta_sum = offd_sumsq;
	real sweep_thresh = 0.5*offd_sumsq/LDSQ;

	/********************** RUN JACOBI SWEEPS **********************/
	integer sweep_counter = 0;
	while(offd_sumsq > tol) {
		integer n_iter = LD-1;
		for(integer iter=0; iter<n_iter; iter++) {
			calc_param_kernel<<<dim_grid_param, dim_block_param>>>(pd_H, pd_S, pd_C, LD, iter);
			jacobi_kernel_step1<<<dim_grid_jacobi, dim_block_jacobi>>>(pd_H, pd_M, pd_S, pd_C, LD, iter, sweep_thresh, stride);
			jacobi_kernel_step2<<<dim_grid_jacobi, dim_block_jacobi>>>(pd_H, pd_M, pd_S, pd_C, LD, iter, sweep_thresh, stride, pd_E);
		}
		sweep_counter++;
		// Transform d_M into off-diagonal matrix of H and calculate the sum of squares.
		thrust::transform(d_H.cbegin(), d_H.cend(), d_N.cbegin(), d_M.begin(), thrust::multiplies<real>());
		offd_sumsq = thrust::inner_product(d_M.cbegin(), d_M.cend(), d_M.cbegin(), 0.0);
		sweep_thresh = 0.5*offd_sumsq/LDSQ;
		delta_sum -= offd_sumsq;
		std::cout << std::fixed
				  << "DiagHessian> Sweep#" << std::setw(12) << sweep_counter
				  << " | ConvergeIndex:"   << std::setw(16) << offd_sumsq
				  << " | TargetValue:"     << std::setw(12) << tol
				  << " | DeltaOffdSumsq:"  << std::setw(12) << delta_sum
				  << std::endl;
		delta_sum = offd_sumsq;
 	}
	// copy diagonal matrix and eigenvectors back to host.
	thrust::copy(d_H.begin(), d_H.end(), h_H.begin());
	thrust::copy(d_E.begin(), d_E.end(), h_E.begin());
	std::cout << "DiagHessian> Done." << std::endl << std::endl;
}
