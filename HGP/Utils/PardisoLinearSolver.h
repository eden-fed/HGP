#pragma once


#include "UTILS/GMM_Macros.h"
#include "Eigen\SparseCore"
#include "pardiso.h"


class PardisoLinearSolver
{
public:

	PardisoLinearSolver(int matrixType=11,  int solver_in = 0);
	~PardisoLinearSolver();

	void setMatrixType(int matrixType);
	bool init(bool isHighlyIndefinite=false);
	bool preprocess();
	bool reorder();
	bool numericalFactorization();
	bool selectiveInverse();
	bool getMatrixInGMMformat(GMMSparseRowMatrix& M);
	bool getMatrixInGMMformat(GMMCompressed1RowMatrix& M);
	bool createPardisoFormatMatrix(const GMMCompressed1RowMatrix& M, std::vector<int>& rowsOfElementsToSet, std::vector<int>& colsOfElementsToSet);
	bool createPardisoFormatMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& M);//here we assume that the matrix is already upper triangular and that the diagonal elements are set
	bool createPardisoFormatMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& M, std::vector<int>& rowsOfElementsToSet, std::vector<int>& colsOfElementsToSet);//here we assume that we need to set elements

	template<class DenseColMatrixType>
	bool solve(DenseColMatrixType& RHS, DenseColMatrixType& X);
	
	template<class DenseColMatrixType>
	bool oneTimeSolve(DenseColMatrixType& RHS, DenseColMatrixType& X);

private:

	void *pt[64];			/* Internal solver memory pointer pt */
	int mtype;				/* matrix type */
	int solver;				/* solver method */
	int error;				/* error output */
	int iparm[64];			/* Pardiso control parameters */
	double dparm[64];		/* Pardiso control parameters */
	int maxfct;				/* Maximum number of numerical factorizations.  */
	int mnum;				/* Which factorization to use. */
	int msglvl;				/* Print statistical information  */
	/* Sparse csr-compressed sparse row format matrix */
	std::vector<double> a;	/* Nonzero values of the coefficient matrix */
	std::vector<int> ja;	/* column indices of the sparse matrix */
	std::vector<int> ia;	/* ia(i) points to the first column index of row i in the array ja */
	/**/
	std::vector<double> b;	/* Right hand side */
	std::vector<double> x;	/* Solution */
	int n;					/* Number of equations in the sparse linear systems */
	int phase;				/* Which phase to perform. */
	double   ddum;			/* Double dummy */
	int      idum;			/* Integer dummy. */
	int nrhs;				/* Number of right-hand sides */
	int mNumProcessors;		/* Number of processors */
	bool mIsSymmetric;
	int mNonzeros;
	bool mPardisoInitialized;

	template<class DenseColMatrixType>
	bool createRHSvector(DenseColMatrixType& RHS);

	template<class DenseColMatrixType>
	bool updateX(DenseColMatrixType& X);
};



template<class DenseColMatrixType>
bool PardisoLinearSolver::solve(DenseColMatrixType& RHS, DenseColMatrixType& X)
{

	if (!createRHSvector(RHS))
		return false;

	/* -------------------------------------------------------------------- */
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */
	phase = 33;

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], &idum, &nrhs, iparm, &msglvl, &b[0], &x[0], &error, dparm);
	if (error != 0) {
		cout << "\nERROR during solution: " << error << endl;
		return false;
	}

	//Number of performed iterative refinement steps=iparm[6]

	if (!updateX(X))
		return false;

	return true;
}

template<class DenseColMatrixType>
bool PardisoLinearSolver::oneTimeSolve(DenseColMatrixType& RHS, DenseColMatrixType& X)
{
	phase = 13;	/*Reordering And Symbolic Factorization*/

	if (!createRHSvector(RHS))
		return false;

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], &idum, &nrhs, iparm, &msglvl, &b[0], &x[0], &error, dparm);

	if (error != 0) {
		cout << "ERROR during oneTimeSolve: " << error << endl;
		return false;
	}

	//Number of performed iterative refinement steps=iparm[6]

	if (!updateX(X))
		return false;

	return true;
}

template<class DenseColMatrixType>
bool PardisoLinearSolver::createRHSvector(DenseColMatrixType& RHS)
{
	nrhs = RHS.ncols();

	b = RHS;
	x.resize(n*nrhs);

	//check vector
	pardiso_chkvec(&n, &nrhs, &b[0], &error);
	if (error != 0) {
		cout << "\nERROR right hand side: " << error << endl;
		return false;
	}
	return true;
}

template<class DenseColMatrixType>
bool PardisoLinearSolver::updateX(DenseColMatrixType& X)
{
	X.assign(x.begin(), x.end());
	return true;
}

