#ifndef DIAG_CUH
#define DIAG_CUH

#include "main.h"
#include <thrust/host_vector.h>

void diag_hessian_gpu(integer LD, real tol);

#endif