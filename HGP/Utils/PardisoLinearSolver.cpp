#include <string>
#include "PardisoLinearSolver.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"


PardisoLinearSolver::PardisoLinearSolver(int matrixType, int solver_in) : mtype(matrixType), solver(solver_in)
{
	mIsSymmetric = (mtype == 2) || (mtype == -2);
	mPardisoInitialized = false;
	idum = 0;
	ddum = 0;
}
bool PardisoLinearSolver::init(bool isHighlyIndefinite)
{
	error = 0;
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);		/*Setup Pardiso control parameters*/
	if (error != 0)
	{
		if (error == -10)
			std::cout << "No license file found \n" << error;
		if (error == -11)
			std::cout << "License is expired \n" << error;
		if (error == -12)
			std::cout << "Wrong username or hostname \n" << error;

		return false;
	}

	/* Numbers of processors, value of OMP_NUM_THREADS */
	char *var = getenv("OMP_NUM_THREADS");
	if (var != NULL)
		sscanf(var, "%d", &mNumProcessors);
	else{
		std::cout << "OMP_NUM_THREADS is not set \n" << error;
		return false;
	}

	mPardisoInitialized = true;

	//*****set parameters for PARDISO:*****
	iparm[2] = mNumProcessors;
	maxfct = 1;
	mnum = 1;
	msglvl = 0; //don't print statistical information
	iparm[7] = 0;//Max numbers of iterative refinement steps (with the default value 0 it may perform 1 or 2 steps)

	/* parallelization of reordering and Factorization steps*/
	iparm[27] = 1;//parallel Reordering METIS
	iparm[23] = 0;//Parallel Numerical Factorization.

	/* selective inverse flags */
	iparm[35] = 1; // not overwrite internal factor L 

	/* improve accuracy for Symmetric indefinite Matrices */
	if (mtype == -2 && isHighlyIndefinite){
		iparm[10] = 1;// symmetric weighted matchings algorithm for pivoting
		iparm[12] = 1;// symmetric weighted matchings algorithm for pivoting
	}

	return true;
}

PardisoLinearSolver::~PardisoLinearSolver()
{
	if (mPardisoInitialized){
		/* -------------------------------------------------------------------- */
		/* ..  Termination and release of memory.                               */
		/* -------------------------------------------------------------------- */
		phase = -1;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, &ia[0], &ja[0], &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
	}

}


void PardisoLinearSolver::setMatrixType(int matrixType)
{
	mtype = matrixType;
	mIsSymmetric = (mtype == 2) || (mtype == -2);
}

//Here we assume that the matrix is already upper triangular (if symmetric) and that the diagonal elements are set.
bool PardisoLinearSolver::createPardisoFormatMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& M)
{
	n = M.rows();

	a.clear();
	ja.clear();
	ia.clear();

	mNonzeros = M.nonZeros();
	a = std::vector<double>(M.valuePtr(), M.valuePtr() + mNonzeros);
	ja = std::vector<int>(M.innerIndexPtr(), M.innerIndexPtr() + mNonzeros);
	ia = std::vector<int>(M.outerIndexPtr(), M.outerIndexPtr() + n + 1);


	//Convert matrix from 0-based C-notation to Fortran 1-based notation.
	for (int i = 0; i < n + 1; i++) {
		ia[i] += 1;
	}
	for (int i = 0; i < mNonzeros; i++) {
		ja[i] += 1;
	}

	//check matrix
	pardiso_chkmatrix(&mtype, &n, &a[0], &ia[0], &ja[0], &error);
	if (error != 0) {
		std::cout << "ERROR in consistency of matrix: " << error << std::endl;
		return false;
	}

	return true;
}

/*Reordering and symbolic factorization (phase 1) 
			+
numerical factorization (phase 2)*/
bool PardisoLinearSolver::preprocess()
{
	phase = 12;	

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], &idum, &idum, iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0) {
		std::cout << "ERROR during factorization: " << error << std::endl;
		return false;
	}

	return true;
}

/*Reordering and symbolic factorization (phase 1)*/
bool PardisoLinearSolver::reorder()
{
	phase = 11;											/*Reordering And Symbolic Factorization*/

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], &idum, &idum, iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0) {
		std::cout << "ERROR during Reordering and symbolic factorization: " << error << std::endl;
		return false;
	}

	return true;
}

/*numerical factorization (phase 2)*/
bool PardisoLinearSolver::numericalFactorization()
{
	phase = 22;	

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], &idum, &idum, iparm, &msglvl, &ddum, &ddum, &error, dparm);

	if (error != 0) {
		std::cout << "\nERROR during numerical factorization: " << error;
		return false;
	}

	return true;
}

bool PardisoLinearSolver::selectiveInverse()
{
	//Inverse factorization 
	if (solver == 0)
	{
		phase = -22;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &a[0], &ia[0], &ja[0], &idum, &idum, iparm, &msglvl, &ddum, &ddum, &error, dparm);

		if (error != 0) {
			std::cout << "\nERROR during inverse: %d", error;
			return false;
		}

		return true;
	}
	else return false;

}

//return full symmetric matrix
bool PardisoLinearSolver::getMatrixInGMMformat(GMMSparseRowMatrix& M)
{

	M.resize(n, n);
	int a_index = 0;
	for (int i = 0; i < n; i++){
		for (int j = ia[i] - 1; j < ia[i + 1] - 1; j++){
			int col = ja[j] - 1;
			int row = i;
			M(row, col) = a[a_index];
			if (row != col)
				M(col, row) = a[a_index];
			a_index++;
		}
	}
	return true;
}

//return an upper triangular matrix
bool PardisoLinearSolver::getMatrixInGMMformat(GMMCompressed1RowMatrix& M)
{
	if (M.nrows() != n || M.ncols() != n)
		return false;

	M.ir.assign(ja.begin(), ja.end());
	M.pr.assign(a.begin(), a.end());
	M.jc.assign(ia.begin(), ia.end());

	return true;
}
