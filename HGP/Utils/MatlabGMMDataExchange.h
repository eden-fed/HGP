#pragma once

#include "Utils/GMM_Macros.h"
#include "Utils/Limit_Util.h"


namespace MatlabGMMDataExchange
{


	int SetEngineSparseMatrix(const char* name, GMMSparseRowMatrix& A);
	int SetEngineSparseMatrix(const char* name, GMMSparseComplexRowMatrix& A);
	int SetEngineDenseMatrix(const char* name, GMMDenseColMatrix& A);
	int SetEngineDenseMatrix(const char* name, GMMDenseComplexColMatrix& A);


	int GetEngineDenseMatrix(const char* name, GMMDenseComplexColMatrix& A);
	int GetEngineDenseMatrix(const char* name, GMMDenseColMatrix& A);
	int GetEngineSparseMatrix(const char* name, GMMSparseRowMatrix& A);
	int GetEngineSparseMatrix(const char* name, GMMSparseComplexRowMatrix& A);
	int GetEngineCompressedSparseMatrix(const char* name, GMMCompressed0ComplexRowMatrix& A);
	int GetEngineCompressedSparseMatrix(const char* name, GMMCompressed0RowMatrix& A);

	template<class SparseRowMatrixType>
	bool isMatrixValid(const SparseRowMatrixType& M);
}


template<class SparseRowMatrixType>
bool MatlabGMMDataExchange::isMatrixValid(const SparseRowMatrixType& M)
{
	gmm::linalg_traits<SparseRowMatrixType>::const_row_iterator rowIter = mat_row_const_begin(M);
	gmm::linalg_traits<SparseRowMatrixType>::const_row_iterator rowIterEnd = mat_row_const_end(M);

	for(; rowIter != rowIterEnd; rowIter++)
	{
		gmm::linalg_traits<SparseRowMatrixType>::const_sub_row_type row = gmm::linalg_traits<SparseRowMatrixType>::row(rowIter);

		gmm::linalg_traits<gmm::linalg_traits<SparseRowMatrixType>::const_sub_row_type>::const_iterator elementIter = vect_const_begin(row);
		gmm::linalg_traits<gmm::linalg_traits<SparseRowMatrixType>::const_sub_row_type>::const_iterator elementIterEnd = vect_const_end(row);

		for(; elementIter != elementIterEnd; elementIter++)
		{
			if(!LIMIT::isFinite(*elementIter)) //this will be true if val is NAN (not a valid number)
			{
				return false;
			}
		}
	}
	return true;
}
