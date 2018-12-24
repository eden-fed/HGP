#pragma once

#include "gmm/gmm.h"


typedef gmm::wsvector<double> GMMSparseVector; //sparse vector which is optimized for write operations
typedef gmm::row_matrix<GMMSparseVector> GMMSparseRowMatrix; //sparse row matrix such that each row is a sparse vector
typedef gmm::csr_matrix<double, 1> GMMCompressedRowMatrix; //sparse row matrix in compressed form with 1-based indexing (like in Fortran/Matlab)
typedef gmm::dense_matrix<double> GMMDenseColMatrix; //dense column major matrix - all memory is contiguous

typedef gmm::wsvector<std::complex<double> > GMMSparseComplexVector; //sparse complex vector which is optimized for write operations
typedef gmm::row_matrix<GMMSparseComplexVector> GMMSparseComplexRowMatrix; //sparse complex row matrix such that each row is a sparse vector


typedef gmm::csr_matrix<double, 0> GMMCompressed0RowMatrix; //sparse real row matrix in compressed form with 0-based indexing (like in C)
typedef gmm::csr_matrix<double, 1> GMMCompressed1RowMatrix; //sparse real row matrix in compressed form with 1-based indexing (like in Fortran/Matlab)

typedef gmm::csr_matrix<std::complex<double>, 0> GMMCompressed0ComplexRowMatrix; //sparse complex row matrix in compressed form with 0-based indexing (like in C)
typedef gmm::csr_matrix<std::complex<double>, 1> GMMCompressed1ComplexRowMatrix; //sparse complex row matrix in compressed form with 1-based indexing (like in Fortran/Matlab)


typedef gmm::dense_matrix<std::complex<double> > GMMDenseComplexColMatrix; //dense complex column major matrix - all memory is contiguous
