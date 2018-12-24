
#include "MatlabGMMDataExchange.h"
#include "MatlabInterface.h"

#include <string>
#include <vector>
#include <algorithm>




int MatlabGMMDataExchange::SetEngineDenseMatrix(const char* name, GMMDenseColMatrix& A)
{
	if(A.ncols() == 0 || A.nrows() == 0)
	{
		return -1;
	}
	MatlabInterface& matlab = MatlabInterface::GetEngine();
	int res = matlab.SetEngineRealMatrix(name, A.nrows(), A.ncols(), &A.front(), true);
	return res;
}



int MatlabGMMDataExchange::SetEngineDenseMatrix(const char* name, GMMDenseComplexColMatrix& A)
{
	if(A.ncols() == 0 || A.nrows() == 0)
	{
		return -1;
	}
	MatlabInterface& matlab = MatlabInterface::GetEngine();
	int res = matlab.SetEngineComplexMatrix(name, A.nrows(), A.ncols(), &A.front(), true);
	return res;
}




int MatlabGMMDataExchange::SetEngineSparseMatrix(const char* name, GMMSparseRowMatrix& A)
{
	int nRows = A.nrows();
	int nCols = A.ncols();

	if(nCols == 0 || nRows == 0)
	{
		return -1;
	}
	std::vector<unsigned> iv;
	std::vector<unsigned> jv;
	std::vector<double> dv;

	for (int i = 0; i < gmm::mat_nrows(A); ++i)
	{
		gmm::wsvector<double>::iterator iter = A[i].begin();

		for ( ; iter != A[i].end(); ++iter)
		{
			int j = iter->first;
			double v = iter->second;

			iv.push_back(i);
			jv.push_back(j);
			dv.push_back(v);
		}
	}

	MatlabInterface& matlab = MatlabInterface::GetEngine();

	if(dv.size() == 0)
	{
		int res = matlab.CreateAllZerosSparseMatrix(name, nRows, nCols);
		return res;
	}

	int res = matlab.SetEngineSparseRealMatrix(name, iv.size(), &iv[0], &jv[0], &dv[0], nRows, nCols);
	return res;
}




int MatlabGMMDataExchange::SetEngineSparseMatrix(const char* name, GMMSparseComplexRowMatrix& A)
{
	int nRows = A.nrows();
	int nCols = A.ncols();

	if(nCols == 0 || nRows == 0)
	{
		return -1;
	}
	std::vector<unsigned> iv;
	std::vector<unsigned> jv;
	std::vector<std::complex<double> > dv;

	for (int i = 0; i < gmm::mat_nrows(A); ++i)
	{
		gmm::wsvector<std::complex<double> >::iterator iter = A[i].begin();

		for ( ; iter != A[i].end(); ++iter)
		{
			int j = iter->first;
			std::complex<double> v = iter->second;

			iv.push_back(i);
			jv.push_back(j);
			dv.push_back(v);
		}
	}

	MatlabInterface& matlab = MatlabInterface::GetEngine();

	if(dv.size() == 0)
	{
		int res = matlab.CreateAllZerosSparseMatrix(name, nRows, nCols);
		return res;
	}

	int res = matlab.SetEngineSparseComplexMatrix(name, iv.size(), &iv[0], &jv[0], &dv[0], nRows, nCols);
	return res;
}


int MatlabGMMDataExchange::GetEngineDenseMatrix(const char* name, GMMDenseComplexColMatrix& A)
{
	unsigned int m = 0;
	unsigned int n = 0;

	MatlabInterface& matlab = MatlabInterface::GetEngine();

	bool matrixExists = matlab.GetMatrixDimensions(name, m, n);
	if(!matrixExists)
	{
		return -1;
	}
	
	A.resize(m, n);

	int res = matlab.GetEngineComplexMatrix(name, m, n, &A.front(), true);

	return res;
}


int MatlabGMMDataExchange::GetEngineDenseMatrix(const char* name, GMMDenseColMatrix& A)
{
	unsigned int m = 0;
	unsigned int n = 0;

	MatlabInterface& matlab = MatlabInterface::GetEngine();

	bool matrixExists = matlab.GetMatrixDimensions(name, m, n);
	if(!matrixExists)
	{
		return -1;
	}

	A.resize(m, n);

	int res = matlab.GetEngineRealMatrix(name, m, n, &A.front(), true);

	return res;
}




int MatlabGMMDataExchange::GetEngineSparseMatrix(const char* name, GMMSparseRowMatrix& A)
{
	std::vector<unsigned int> rowind;
	std::vector<unsigned int> colind;
	std::vector<double> vals;
	unsigned int m, n, nentries;

	MatlabInterface& matlab = MatlabInterface::GetEngine();

	int res = matlab.GetSparseRealMatrix(name, rowind, colind, vals, nentries, m, n);

	A.resize(m, n);

	for(unsigned int i = 0; i < nentries; i++)
	{
		A(rowind[i], colind[i]) = vals[i];
	}

	return 0;
}


int MatlabGMMDataExchange::GetEngineSparseMatrix(const char* name, GMMSparseComplexRowMatrix& A)
{
	std::vector<unsigned int> rowind;
	std::vector<unsigned int> colind;
	std::vector<std::complex<double> > vals;
	unsigned int m, n, nentries;

	MatlabInterface& matlab = MatlabInterface::GetEngine();
	
	int res = matlab.GetSparseComplexMatrix(name, rowind, colind, vals, nentries, m, n);

	A.resize(m, n);

	for(unsigned int i = 0; i < nentries; i++)
	{
		A(rowind[i], colind[i]) = vals[i];
	}

	return 0;
}




int MatlabGMMDataExchange::GetEngineCompressedSparseMatrix(const char* name, GMMCompressed0RowMatrix& A)
{
	MatlabInterface& matlab = MatlabInterface::GetEngine();

	unsigned int m, n;
	int res = matlab.GetEngineEncodedSparseRealMatrix(name, A.ir, A.jc, A.pr, m, n);
	if(res == 0)
	{
		A.nr = m;
		A.nc = n;
		return 0;
	}
	else
	{
		return -1;
	}
}



int MatlabGMMDataExchange::GetEngineCompressedSparseMatrix(const char* name, GMMCompressed0ComplexRowMatrix& A)
{
	MatlabInterface& matlab = MatlabInterface::GetEngine();

	unsigned int m, n;
	int res = matlab.GetEngineEncodedSparseComplexMatrix(name, A.ir, A.jc, A.pr, m, n);
	if(res == 0)
	{
		A.nr = m;
		A.nc = n;
		return 0;
	}
	else
	{
		return -1;
	}
}


