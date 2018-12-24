#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <windows.h>

#include "MatlabInterface.h"

#include "engine.h"  // Matlab engine header
#define NO_LOGGER
// 
// #ifndef TARGET_API_VERSION
// #error TARGET_API_VERSION must be defined
// #endif
// 
// #if TARGET_API_VERSION <= 700 //R2017b
// #error Compiling the code in this file on Matlab version R2017b or older will cause run time errors due to the different encoding of complex valued arrays
// #endif
// 


#ifndef NO_LOGGER
    #include "logger.h"
    extern Logger Log; // the log object named Log is instantiated elsewhere
    #define ERROR(msg)   Log.error() << msg << std::endl
    #define WARNING(msg) Log.warn() << msg << std::endl
    #define MESSAGE(msg) Log.print() << msg << std::endl
    #define MESSAGE2(msg) Log.print() << msg << std::flush
#else // ifndef NO_LOGGER
    #define ERROR(msg)   std::cerr << "ERROR: "   << msg << std::endl
    #define WARNING(msg) std::cerr << "WARNING: " << msg << std::endl
    #define MESSAGE(msg) std::cerr << msg << std::endl
	#define MESSAGE2(msg) std::cout << msg << std::flush
#endif // ifndef NO_LOGGER

using namespace std;

////////////////////////////////////////////////////////////
// Explicit instantiation of exported template methods

// Creates a matrix in matlab.
template int MatlabInterface::SetEngineIndexMatrix         <int>(const char *name, unsigned int m, unsigned int n, const                  int *vals, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template int MatlabInterface::SetEngineIndexMatrix<unsigned int>(const char *name, unsigned int m, unsigned int n, const         unsigned int *vals, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template int MatlabInterface::SetEngineRealMatrix       <double>(const char *name, unsigned int m, unsigned int n, const               double *vals, bool colmaj);
template int MatlabInterface::SetEngineComplexMatrix    <double>(const char *name, unsigned int m, unsigned int n, const std::complex<double> *vals, bool colmaj);
template int MatlabInterface::SetEngineRealMatrix       <float>(const char *name, unsigned int m, unsigned int n, const               float *vals, bool colmaj);
template int MatlabInterface::SetEngineComplexMatrix    <float>(const char *name, unsigned int m, unsigned int n, const std::complex<float> *vals, bool colmaj);

template
int MatlabInterface::SetEngineEncodedSparseRealMatrix<unsigned int, double>(const char *name, unsigned int n,
        const unsigned int *rowind, const unsigned int *colind, const double *vals);

template
int MatlabInterface::SetEngineSparseRealMatrix<int, double>(const char *name, unsigned int n,
        const int *rowind, const int *colind, const double *vals, unsigned int nrows, unsigned int ncols);
template
int MatlabInterface::SetEngineSparseRealMatrix<unsigned int, double>(const char *name, unsigned int n,
        const unsigned int *rowind, const unsigned int *colind, const double *vals, unsigned int nrows, unsigned int ncols);

template
int MatlabInterface::SetEngineEncodedSparseComplexMatrix<unsigned int, double>(const char *name, unsigned int n,
        const unsigned int *rowind, const unsigned int *colind, const std::complex<double> *vals);

template
int MatlabInterface::SetEngineSparseComplexMatrix<unsigned int, double>(const char *name, unsigned int n,
        const unsigned int *rowind, const unsigned int *colind, const std::complex<double> *vals, unsigned int nrows, unsigned int ncols);

// Reads a matrix from matlab.
template int MatlabInterface::GetEngineIndexMatrix<int>(const char *name, unsigned int m, unsigned int n,         int *vals, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template int MatlabInterface::GetEngineIndexMatrix<unsigned int>(const char *name, unsigned int m, unsigned int n,         unsigned int *vals, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template int MatlabInterface::GetEngineRealMatrix       <double>(const char *name, unsigned int m, unsigned int n,               double *vals, bool colmaj);
template int MatlabInterface::GetEngineComplexMatrix    <double>(const char *name, unsigned int m, unsigned int n, std::complex<double> *vals, bool colmaj);
template int MatlabInterface::GetEngineRealMatrix       <float>(const char *name, unsigned int m, unsigned int n,               float *vals, bool colmaj);
template int MatlabInterface::GetEngineComplexMatrix    <float>(const char *name, unsigned int m, unsigned int n, std::complex<float> *vals, bool colmaj);


// Reads a variable-sized matrix from matlab.
template int MatlabInterface::GetEngineIndexMatrix<int>(const char *name, unsigned int& m, unsigned int& n, std::vector<         int>& vals, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template int MatlabInterface::GetEngineIndexMatrix<unsigned int>(const char *name, unsigned int& m, unsigned int& n, std::vector<         unsigned int>& vals, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template int MatlabInterface::GetEngineRealMatrix       <double>(const char *name, unsigned int& m, unsigned int& n, std::vector<         double      >& vals, bool colmaj);
template int MatlabInterface::GetEngineComplexMatrix    <double>(const char *name, unsigned int& m, unsigned int& n, std::vector<std::complex<double> >& vals, bool colmaj);


// matrix creation helper functions
template mxArray* MatlabInterface::CreateIndexMatrix<unsigned int>(unsigned int m, unsigned int n, const         unsigned int *vals, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template mxArray* MatlabInterface::CreateRealMatrix       <double>(unsigned int m, unsigned int n, const               double *vals, bool colmaj);
template mxArray* MatlabInterface::CreateComplexMatrix    <double>(unsigned int m, unsigned int n, const std::complex<double> *vals, bool colmaj);

template
mxArray* MatlabInterface::CreateEncodedSparseRealMatrix<unsigned int, double>(unsigned int n,
        const unsigned int *rowind, const unsigned int *colind, const double *vals);

template
mxArray* MatlabInterface::CreateEncodedSparseComplexMatrix<unsigned int, double>(unsigned int n,
        const unsigned int *rowind, const unsigned int *colind, const std::complex<double> *vals);

// matrix copy-back helper functions
template int MatlabInterface::CopyFromIndexMatrix<unsigned int>(mxArray *M, unsigned int m, unsigned int n,         unsigned int *dest, bool colmaj); // shifts indices by 1 (C/C++ vs matlab)
template int MatlabInterface::CopyFromRealMatrix       <double>(mxArray *M, unsigned int m, unsigned int n,               double *dest, bool colmaj);
template int MatlabInterface::CopyFromComplexMatrix    <double>(mxArray *M, unsigned int m, unsigned int n, std::complex<double> *dest, bool colmaj);

////////////////////////////////////////////////////////////

DWORD MatlabInterface::mThreadIdUsedToOpenEngine = 0;

MatlabInterface::MatlabInterface() : m_ep(NULL), mOutputStringBuffer(NULL)
{

	mOutputStringBuffer = new char[MAX_OUTPUT_BUFFER_SIZE];
	memset(mOutputStringBuffer, 0, MAX_OUTPUT_BUFFER_SIZE);

	cout << "Starting Matlab engine" << endl;

	EngineOpen();
	EngineSetupUI();
}



MatlabInterface::~MatlabInterface()
{
//  EngineClose();
	delete[] mOutputStringBuffer;
	mOutputStringBuffer = NULL;

	//make sure that Matlab output buffer is nullified to prevent Matlab from accessing memory that is deallocated.
	engOutputBuffer(m_ep, NULL, 0);
}


void MatlabInterface::EngineOpen()
{
	// Start the MATLAB engine locally by executing the string
	// "matlab"
	//
	// To start the session on a remote host, use the name of
	// the host as the string rather than \0
	//
	// For more complicated cases, use any string with whitespace,
	// and that string will be executed literally to start MATLAB
	//
	if (!(m_ep = engOpen("\0"))) {
		ERROR("Can't start MATLAB engine");
		return;
	}
	//let Matlab know where to issue to textual output
	int res = engOutputBuffer(m_ep, mOutputStringBuffer, MAX_OUTPUT_BUFFER_SIZE - 1);

	if (res != 0)
	{
		ERROR("Unable to set Matlab output buffer");
	}

	mThreadIdUsedToOpenEngine = GetCurrentThreadId();
}


void MatlabInterface::EngineSetupUI()
{
	int res = engSetVisible(m_ep, false); //turn the white matlab command window OFF
	if (res != 0)
	{
		ERROR("Unable to turn off Matlab engine's white console window");
	}
	res = engEvalString(m_ep, "desktop");  //turn the entire GUI of matlab ON
	if (res != 0)
	{
		ERROR("Unable to turn on Matlab full GUI");
	}
}


void MatlabInterface::EngineClose()
{
    engClose(m_ep);
    m_ep = NULL;
}


void
MatlabInterface::DestroyMatrix(mxArray *&M)
{
    if (M) {
        mxDestroyArray(M);
        M = NULL;
    }
}



// This also shifts all the values by 1 to account for the difference
// between matlab and C/C++ indexing
template <typename T>
mxArray*
MatlabInterface::CreateIndexMatrix(unsigned int m, unsigned int n, const T *vals, bool colmaj)
{
    mxArray *M = mxCreateDoubleMatrix(m, n, mxREAL);
    if (vals != NULL) {
        mxDouble *pM = mxGetDoubles(M);
        // note that matlab expects the data in column-major order
        for (unsigned int j = 0; j < n; ++j)
            for (unsigned int i = 0; i < m; ++i)
            {
                unsigned int idxM = j*m+i;
                unsigned int idx = colmaj ? idxM : i*n+j;
                pM[idxM] = double(vals[idx]+1);
            }
    }
    return M;
}


//this function can be applied to single or double precision floating point values (T), but either way, the created Matlab matrix will be double precision
template <typename T>
mxArray*
MatlabInterface::CreateRealMatrix(unsigned int m, unsigned int n, const T *vals, bool colmaj)
{
    mxArray *M = mxCreateDoubleMatrix(m, n, mxREAL);
    if (vals != NULL) {
		mxDouble *pM = mxGetDoubles(M);
        // note that matlab expects the data in column-major order

		if(colmaj) //column-major
		{
			int numElements = m*n;

			for(int i = 0; i < numElements; i++)
			{
				pM[i] = (double)vals[i];
			}
		}
		else //row-major
		{
			for (unsigned int j = 0; j < n; ++j)
			{
				for (unsigned int i = 0; i < m; ++i)
				{
					unsigned int idxM = j*m+i;
					unsigned int idx = i*n+j;
					pM[idxM] = double(vals[idx]);
				}
			}
		}
    }
    return M;
}


//this function can be applied to single or double precision floating point values (T), but either way, the created Matlab matrix will be double precision
template <typename T>
mxArray*
MatlabInterface::CreateComplexMatrix(unsigned int m, unsigned int n, const std::complex<T> *vals, bool colmaj)
{
	mxArray *M = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
	if (vals != NULL) {
		mxComplexDouble *p = mxGetComplexDoubles(M);

		// note that matlab expects the data in column-major order

		if(colmaj) //column-major
		{
			int numElements = m*n;
			memcpy(p, vals, numElements * mxGetElementSize(M));
		}
		else //row-major
		{
			for (unsigned int j = 0; j < n; ++j)
			{
				for (unsigned int i = 0; i < m; ++i)
				{
					unsigned int idxM = j*m+i;
					unsigned int idx = i*n+j;
					p[idxM].real = vals[idx].real();
					p[idxM].imag = vals[idx].imag();
				}
			}
		}
	}
	return M;
}



template <typename IndexType, typename ValueType>
mxArray*
MatlabInterface::CreateEncodedSparseRealMatrix(unsigned int n,
        const IndexType *rowind, const IndexType *colind, const ValueType *vals)
{
    mxArray *M = mxCreateDoubleMatrix(n, 3, mxREAL);
    if (vals != NULL) {
        assert(rowind != NULL);
        assert(colind != NULL);
        mxDouble *pM = mxGetDoubles(M);
        // note that matlab expects the data in column-major order
        for (unsigned int i = 0; i < n; ++i)
        {
            pM[0*n+i] = double(rowind[i]+1);
            pM[1*n+i] = double(colind[i]+1);
            pM[2*n+i] = double(  vals[i]  );
        }
    }
    return M;
}



template <typename IndexType, typename ValueType>
mxArray* MatlabInterface::CreateEncodedSparseComplexMatrix(unsigned int n, const IndexType *rowind, const IndexType *colind, const std::complex<ValueType> *vals)
{
    mxArray *M = mxCreateDoubleMatrix(n, 3, mxCOMPLEX);
    if (vals != NULL) {
        assert(rowind != NULL);
        assert(colind != NULL);
		mxComplexDouble* pM = mxGetComplexDoubles(M);

		// note that matlab expects the data in column-major order
		for (unsigned int i = 0; i < n; ++i)
		{
			pM[0 * n + i].real = double(rowind[i] + 1);
			pM[0 * n + i].imag = double(0);

			pM[1 * n + i].real = double(colind[i] + 1);
			pM[1 * n + i].imag = double(0);

			pM[2 * n + i].real = double(vals[i].real());
			pM[2 * n + i].imag = double(vals[i].imag());
		}
    }
    return M;
}



mxArray*
MatlabInterface::CreateStringArray(unsigned int n, const char **val)
{
    assert(val != NULL);
    mxArray *M = mxCreateCharMatrixFromStrings(n, val);
    assert(M);
    return M;
}




// This also shifts all the values by 1 to account for the difference
// between matlab and C/C++ indexing
template <typename T>
int
MatlabInterface::CopyFromIndexMatrix(mxArray *M, unsigned int m, unsigned int n, T *dest, bool colmaj)
{
    assert(dest != NULL);
    if ( mxGetM(M) != m || mxGetN(M) != n ) {
        ERROR("CopyFromIndexMatrix: expected size " << m << "x" << n << ". Got " << mxGetM(M) << "x" << mxGetN(M) << ".");
        return 1;
    }
    assert(mxGetM(M) == m);
    assert(mxGetN(M) == n);
	mxDouble *pM = mxGetDoubles(M);
    // note that matlab expects the data in column-major order
    for (unsigned int j = 0; j < n; ++j)
        for (unsigned int i = 0; i < m; ++i)
        {
            unsigned int idxM = j*m+i;
            unsigned int idx = colmaj ? idxM : i*n+j;
            dest[idx] = T(pM[idxM]-1);
        }

    return 0;
}

template <typename T>
int
MatlabInterface::CopyFromRealMatrix(mxArray *M, unsigned int m, unsigned int n, T *dest, bool colmaj)
{
    assert(dest != NULL);
    if ( mxGetM(M) != m || mxGetN(M) != n ) {
        ERROR("CopyFromRealMatrix: expected size " << m << "x" << n << ". Got " << mxGetM(M) << "x" << mxGetN(M) << ".");
        return 1;
    }
    assert(mxGetM(M) == m);
    assert(mxGetN(M) == n);
	mxDouble *pM = mxGetDoubles(M);
    // note that matlab expects the data in column-major order
    for (unsigned int j = 0; j < n; ++j)
        for (unsigned int i = 0; i < m; ++i)
        {
            unsigned int idxM = j*m+i;
            unsigned int idx = colmaj ? idxM : i*n+j;
            dest[idx] = T(pM[idxM]);
        }

    return 0;
}


template <typename T>
int
MatlabInterface::CopyFromComplexMatrix(mxArray *M, unsigned int m, unsigned int n, std::complex<T> *dest, bool colmaj)
{
    assert(dest != NULL);
    if ( mxGetM(M) != m || mxGetN(M) != n ) {
        ERROR("CopyFromComplexMatrix: expected size " << m << "x" << n << ". Got " << mxGetM(M) << "x" << mxGetN(M) << ".");
        return 1;
    }
    assert(mxGetM(M) == m);
    assert(mxGetN(M) == n);
	assert(mxIsComplex(M));

	mxComplexDouble *pM = mxGetComplexDoubles(M);
//	bool pure_real = (pMi == NULL); // seems to be NULL when the matrix doesn't have complex values

    // note that matlab expects the data in column-major order
    for (unsigned int j = 0; j < n; ++j)
        for (unsigned int i = 0; i < m; ++i)
        {
            unsigned int idxM = j*m+i;
            unsigned int idx = colmaj ? idxM : i*n+j;
//			dest[idx] = std::complex<T>(T(pMr[idxM]), pure_real ? T(0) : T(pMi[idxM]));
			dest[idx] = std::complex<T>(T(pM[idxM].real), T(pM[idxM].imag));
        }

    return 0;
}

// Creates a matrix in the Matlab engine.
template <typename T>
int
MatlabInterface::SetEngineIndexMatrix(const char *name, unsigned int m, unsigned int n, const T *vals, bool colmaj) // shifts indices by 1 (C/C++ vs matlab)
{
    mxArray *ary = CreateIndexMatrix(m, n, vals, colmaj);
    assert(ary);
    int res = engPutVariable(name, ary);
    mxDestroyArray(ary);
    return res;
}

template <typename T>
int
MatlabInterface::SetEngineRealMatrix(const char *name, unsigned int m, unsigned int n, const T *vals, bool colmaj)
{
    mxArray *ary = CreateRealMatrix(m, n, vals, colmaj);
    assert(ary);
    int res = engPutVariable(name, ary);
    mxDestroyArray(ary);
    return res;
}

template <typename T>
int
MatlabInterface::SetEngineComplexMatrix(const char *name, unsigned int m, unsigned int n, const std::complex<T> *vals, bool colmaj)
{
    mxArray *ary = CreateComplexMatrix(m, n, vals, colmaj);
    assert(ary);
    int res = engPutVariable(name, ary);
    mxDestroyArray(ary);
    return res;
}

template <typename IndexType, typename ValueType>
int
MatlabInterface::SetEngineEncodedSparseRealMatrix(const char *name, unsigned int n,
        const IndexType *rowind, const IndexType *colind, const ValueType *vals)
{
    mxArray *ary = CreateEncodedSparseRealMatrix(n, rowind, colind, vals);
    assert(ary);
    int res = engPutVariable(name, ary);
    mxDestroyArray(ary);
    return res;
}

template <typename IndexType, typename ValueType>
int
MatlabInterface::SetEngineSparseRealMatrix(const char *name, unsigned int n,
        const IndexType *rowind, const IndexType *colind, const ValueType *vals, unsigned int nrows, unsigned int ncols)
{
    int res = SetEngineEncodedSparseRealMatrix(name, n, rowind, colind, vals);
    if (res != 0) return res;
    char cmd[1024];
    if (ncols == 0 || nrows == 0) {
        sprintf(cmd, "%s = sparse(%s(:,1), %s(:,2), %s(:,3));", name, name, name, name);
    }
    else {
        assert(ncols > 0 && nrows > 0);
        sprintf(cmd, "%s = sparse(%s(:,1), %s(:,2), %s(:,3), %d, %d);", name, name, name, name, nrows, ncols);
    }
    res = engEvalString(m_ep, cmd);
    return res;
}

template <typename IndexType, typename ValueType>
int
MatlabInterface::SetEngineEncodedSparseComplexMatrix(const char *name, unsigned int n,
        const IndexType *rowind, const IndexType *colind, const std::complex<ValueType> *vals)
{
    mxArray *ary = CreateEncodedSparseComplexMatrix(n, rowind, colind, vals);
    assert(ary);
    int res = engPutVariable(name, ary);
    mxDestroyArray(ary);
    return res;
}

template <typename IndexType, typename ValueType>
int
MatlabInterface::SetEngineSparseComplexMatrix(const char *name, unsigned int n,
        const IndexType *rowind, const IndexType *colind, const std::complex<ValueType> *vals, unsigned int nrows, unsigned int ncols)
{
    int res = SetEngineEncodedSparseComplexMatrix(name, n, rowind, colind, vals);
    if (res != 0) return res;
    char cmd[1024];
    if (ncols == 0 || nrows == 0) {
        sprintf(cmd, "%s = sparse(%s(:,1), %s(:,2), %s(:,3));", name, name, name, name);
    }
    else {
        assert(ncols > 0 && nrows > 0);
        sprintf(cmd, "%s = sparse(%s(:,1), %s(:,2), %s(:,3), %d, %d);", name, name, name, name, nrows, ncols);
    }
    res = engEvalString(m_ep, cmd);
    return res;
}

int MatlabInterface::CreateAllZerosSparseMatrix(const char* name, int nRows, int nCols)
{
	char cmd[1024];
	sprintf(cmd, "%s = sparse(%d, %d);", name, nRows, nCols);
	int res = engEvalString(m_ep, cmd);
	return res;
}

// Create a string in matlab.
int
MatlabInterface::SetEngineStringArray(const char *name, unsigned int n, const char **val)
{
    assert(val != NULL);
    mxArray *ary = CreateStringArray(n, val);
    assert(ary);
    int res = engPutVariable(name, ary);
    mxDestroyArray(ary);
    return res;
}


// Reads a matrix from the matlab engine.
template <typename T>
int
MatlabInterface::GetEngineIndexMatrix(const char *name, unsigned int m, unsigned int n, T *dest, bool colmaj)
{
    assert(dest != NULL);
    mxArray *ary = engGetVariable(name);
    if (ary == NULL) return 2;
    assert(ary != NULL);
    int res = CopyFromIndexMatrix(ary, m, n, dest, colmaj);
    mxDestroyArray(ary);
    return res;
}

template <typename T>
int
MatlabInterface::GetEngineRealMatrix(const char *name, unsigned int m, unsigned int n, T *dest, bool colmaj)
{
    assert(dest != NULL);
    mxArray *ary = engGetVariable(name);
    if (ary == NULL) return 2;
    assert(ary != NULL);
    int res = CopyFromRealMatrix(ary, m, n, dest, colmaj);
    mxDestroyArray(ary);
    return res;
}



template <typename T>
int
MatlabInterface::GetEngineComplexMatrix(const char *name, unsigned int m, unsigned int n, std::complex<T> *dest, bool colmaj)
{
    assert(dest != NULL);
    mxArray *ary = engGetVariable(name);
    if (ary == NULL) return 2;
    assert(ary != NULL);
    int res = CopyFromComplexMatrix(ary, m, n, dest, colmaj);
    mxDestroyArray(ary);
    return res;
}

// Reads a matrix from the matlab engine.
template <typename T>
int
MatlabInterface::GetEngineIndexMatrix(const char* name, unsigned int& m, unsigned int& n, std::vector<T>& dest, bool colmaj)
{
    mxArray *ary = engGetVariable(name);
    if (ary == NULL)
    {
    	m = 0;
    	n = 0;
    	dest.clear();
    	return 2;
    }
    assert(ary != NULL);
    m = mxGetM(ary);
    n = mxGetN(ary);
    dest.resize(m*n);
    
    int res = CopyFromIndexMatrix(ary, m, n, &dest[0], colmaj);
    mxDestroyArray(ary);
    return res;
}

template <typename T>
int
MatlabInterface::GetEngineRealMatrix(const char* name, unsigned int& m, unsigned int& n, std::vector<T>& dest, bool colmaj)
{
    mxArray *ary = engGetVariable(name);
    if (ary == NULL)
    {
    	m = 0;
    	n = 0;
    	dest.clear();
    	return 2;
    }
    assert(ary != NULL);
    m = mxGetM(ary);
    n = mxGetN(ary);
    dest.resize(m*n);
    
    int res = CopyFromRealMatrix(ary, m, n, &dest[0], colmaj);
    mxDestroyArray(ary);
    return res;
}



template <typename T>
int
MatlabInterface::GetEngineComplexMatrix(const char *name, unsigned int& m, unsigned int& n, std::vector<std::complex<T> >& dest, bool colmaj)
{
    mxArray *ary = engGetVariable(name);
    if (ary == NULL)
    {
    	m = 0;
    	n = 0;
    	dest.clear();
    	return 2;
    }
    assert(ary != NULL);
    m = mxGetM(ary);
    n = mxGetN(ary);
    dest.resize(m*n);
    
    int res = CopyFromComplexMatrix(ary, m, n, &dest[0], colmaj);
    mxDestroyArray(ary);
    return res;
}


//evaluate a single line. returns non-zero on error.
int MatlabInterface::Eval(const char *matlab_code)
{
	int res = engEvalString(m_ep, matlab_code);
	//if ( 1&& res != 0 )	{
	int retries = 5; // doesn't help
	while ( res != 0 )	{
		ERROR("Error " << res << " running matlab command \"" << matlab_code << "\"");
		if ( --retries < 0 )
			break;
		ERROR("Resetting the engine and trying again " << retries);
		//assert(0);
		GetEngine(true);
		//res = engEvalString(m_ep, "1;"); // first eval after a problem usually fails
		res = engEvalString(m_ep, matlab_code);
	}
	return res;
}


void MatlabInterface::EvalToCout(const char *matlab_code, bool saveCommandString)
{
	//save the command string in Matlab
	if (saveCommandString) {
//		Log.println(matlab_code);
		SetEngineStringArray("LAST_EVAL_TO_COUT", 1, &matlab_code);
		Eval("if ~exist('EVAL_TO_COUT', 'var'); EVAL_TO_COUT={LAST_EVAL_TO_COUT}; end; EVAL_TO_COUT = [LAST_EVAL_TO_COUT; EVAL_TO_COUT(1:min(99,end))];");
	}

	Eval(matlab_code);

	if(mOutputStringBuffer[0] != 0)
	{
		MESSAGE2(mOutputStringBuffer);
	}
}

std::string MatlabInterface::EvalToString(const char *matlab_code, bool saveCommandString)
{
	int res;
	return EvalToString(matlab_code, res, saveCommandString);
}

std::string MatlabInterface::EvalToString(const char *matlab_code, int &res, bool saveCommandString)
{
	//save the command string in Matlab
	if (saveCommandString) {
		//		Log.println(matlab_code);
		SetEngineStringArray("LAST_EVAL_TO_COUT", 1, &matlab_code);
		Eval("if ~exist('EVAL_TO_COUT', 'var'); EVAL_TO_COUT={LAST_EVAL_TO_COUT}; end; EVAL_TO_COUT = [LAST_EVAL_TO_COUT; EVAL_TO_COUT(1:min(99,end))];");
	}

	res = Eval(matlab_code); // This returns an API runtime error and not a script error.
	if (res != 0) {
		std::ostringstream oss;
		oss << "ERROR: Matlab command failed with error code " << res << ".\n";
		return oss.str();
	}

// 	if (buf[0] == '>' && buf[1] == '>' && buf[2] == ' ')
// 		buf += 3;
// 	if (buf[0] == '\n') ++buf;

	return std::string(mOutputStringBuffer);
}


int
MatlabInterface::AddScriptPath(const char *path)
{
    std::string matlab_addpath = std::string("addpath('") + path + "')";
    int res = this->Eval(matlab_addpath.c_str());
    return res;
}


// Get a singleton MatlabInterface object
MatlabInterface &MatlabInterface::GetEngine(bool restart)
{
	static MatlabInterface matlab;
	
	DWORD currentThreadID = GetCurrentThreadId();
	//cout << "Current thread id: " << currentThreadID << endl;

	if (restart || (currentThreadID != mThreadIdUsedToOpenEngine))
	{
		matlab.EngineClose();
		matlab.EngineOpen();
	}

	return matlab;
}


int MatlabInterface::GetSparseRealMatrix(const char* name, std::vector<unsigned int>& rowind, std::vector<unsigned int>& colind, std::vector<double>& vals, unsigned int& nentries, unsigned int& m, unsigned int& n)
{ 
	assert(name != NULL && name[0] != 0);

	bool res = GetMatrixDimensions(name, m, n);

	if(!res || m <= 0 || n <= 0)
	{
		return -1;
	}

	char cmd[1024];
	std::string tempName("temp_");
	tempName = tempName + name;

	sprintf(cmd, "clear %s; [%s(:, 1), %s(:, 2), %s(:, 3)] = find(%s);", tempName.c_str(), tempName.c_str(), tempName.c_str(), tempName.c_str(), name);

	int result = engEvalString(m_ep, cmd);
	
    mxArray *M = engGetVariable(tempName.c_str());
  
	if(M == NULL)
	{
        WARNING("matrix " << name << " could not be loaded");
        return - 1;
    }
    
	nentries =  mxGetM(M);
    int ncols =  mxGetN(M);
    assert(ncols == 3); 
    
	mxDouble *pM = mxGetDoubles(M);
	if(pM == NULL) return -1;

    rowind.resize(nentries);
    colind.resize(nentries);
    vals.resize(nentries);

    for(unsigned int i = 0; i < nentries; ++i)
    {
        rowind[i] = (unsigned int)(pM[           i]) - 1; assert(rowind[i] >= 0); 
        colind[i] = (unsigned int)(pM[  nentries+i]) - 1; assert(colind[i] >= 0); 
        vals  [i] = (double      )(pM[2*nentries+i]);
    }
	mxDestroyArray(M);

	sprintf(cmd, "clear %s;", tempName.c_str());
	result = engEvalString(m_ep, cmd);

	return result;
}


int MatlabInterface::GetSparseComplexMatrix(const char* name, std::vector<unsigned int>& rowind, std::vector<unsigned int>& colind, std::vector<std::complex<double> >& vals, unsigned int& nentries, unsigned int& m, unsigned int& n)
{ 
	assert(name != NULL && name[0] != 0);

	bool res = GetMatrixDimensions(name, m, n);

	if(!res || m <= 0 || n <= 0)
	{
		return -1;
	}

	char cmd[1024];
	std::string tempName("temp_");
	tempName = tempName + name;

	sprintf(cmd, "clear %s; [%s(:, 1), %s(:, 2), %s(:, 3)] = find(%s);", tempName.c_str(), tempName.c_str(), tempName.c_str(), tempName.c_str(), name);

	int result = engEvalString(m_ep, cmd);

	mxArray *M = engGetVariable(tempName.c_str());

	if(M == NULL)
	{
		WARNING("matrix " << name << " could not be loaded");
		return - 1;
	}

	nentries =  mxGetM(M);
	int ncols =  mxGetN(M);
	assert(ncols == 3); 

	mxComplexDouble *pM = mxGetComplexDoubles(M);

	if(pM == NULL) return -1;

	rowind.resize(nentries);
	colind.resize(nentries);
	vals.resize(nentries);

	for(unsigned int i = 0; i < nentries; ++i)
	{
		rowind[i] = (unsigned int)(pM[i].real) - 1;
		assert(rowind[i] >= 0); 
		colind[i] = (unsigned int)(pM[nentries+i].real) - 1;
		assert(colind[i] >= 0); 
		vals[i] = std::complex<double>(pM[2*nentries+i].real, pM[2*nentries+i].imag);
	}
	mxDestroyArray(M);

	sprintf(cmd, "clear %s;", tempName.c_str());
	result = engEvalString(m_ep, cmd);

	return result;
}



int MatlabInterface::GetEngineEncodedSparseRealMatrix(const char* name, std::vector<unsigned int>& Ir, std::vector<unsigned int>& Jc, std::vector<double>& vals, unsigned int& m, unsigned int& n)
{ 
	assert(name != NULL && name[0] != 0);

	mxArray *M = engGetVariable(name);

	if(M == NULL || !mxIsSparse(M))
	{ 
		WARNING("matrix " << name << " could not be loaded");
		return -1;
	}
	m =  mxGetM(M);
	n =  mxGetN(M);

	mwIndex nzmax =  mxGetNzmax(M);

	mxDouble *pM = mxGetDoubles(M);
	mwIndex *ir = mxGetIr(M);
	mwIndex *jc = mxGetJc(M);
	if(pM == NULL || ir == NULL || jc == NULL) return -1;

	Ir.resize(nzmax);
	Jc.resize(n + 1);
	vals.resize(nzmax);

	memcpy(pM, &vals.front(), vals.size()*sizeof(vals.front()));
	memcpy(ir, &Ir.front(), Ir.size()*sizeof(Ir.front()));
	memcpy(jc, &Jc.front(), Jc.size()*sizeof(Jc.front()));

	return 0;
}

int MatlabInterface::GetEngineEncodedSparseComplexMatrix(const char* name, std::vector<unsigned int>& Ir, std::vector<unsigned int>& Jc, std::vector<std::complex<double> >& vals, unsigned int& m, unsigned int& n)
{ 
	assert(name != NULL && name[0] != 0);

	mxArray *M = engGetVariable(name);

	if(M == NULL || !mxIsSparse(M))
	{ 
		WARNING("matrix " << name << " could not be loaded");
		return -1;
	}
	m =  mxGetM(M);
	n =  mxGetN(M);

	mwIndex nzmax =  mxGetNzmax(M);

	mxComplexDouble *pM = mxGetComplexDoubles(M);
	mwIndex *ir = mxGetIr(M);
	mwIndex *jc = mxGetJc(M);
	if(pM == NULL || ir == NULL || jc == NULL) return -1;

	Ir.resize(nzmax);
	Jc.resize(n + 1);
	vals.resize(nzmax);

	for(int i = 0; i < nzmax; i++)
	{
		vals[i] = std::complex<double>(pM[i].real, pM[i].imag);
	}
	memcpy(&Ir.front(), ir, Ir.size()*sizeof(Ir.front()));
	memcpy(&Jc.front(), jc, Jc.size()*sizeof(Jc.front()));

	return 0;
}


bool MatlabInterface::GetMatrixDimensions(const char* variableName, unsigned int& m, unsigned int& n)
{
	m = 0;
	n = 0;

	std::string command("size(");
	command = command + variableName;
	command = command + ")";

	int res = Eval(command.c_str());
	if(res != 0)
	{
		return false;
	}
	std::stringstream stream(std::stringstream::in | std::stringstream::out);

	stream << mOutputStringBuffer;

	std::string str1; //to capture the "ans"
	std::string str2; //to capture the "="
	stream >> str1 >> str2;

	if((stream >> m >> n) && m > 0 && n > 0)
	{
		return true;
	}
	else
	{
		return false;
	}

	return true;
}

int MatlabInterface::engPutVariable(Engine *ep, const char *name, const mxArray *arr)
{
	int ret = 0;

	ret = ::engPutVariable(m_ep, name, arr);

	int retries = 5; // doesn't help
	while ( ret != 0 ) {
		ERROR("!! Error in engPutVariable() " << name);
		if ( --retries < 0 )
			break;
		ERROR("Resetting the engine and trying again " << retries);
		GetEngine(true);
		ret = ::engPutVariable(m_ep, name, arr);
	}

	return ret;
}

// Patched engPutVariable() to support struct members, e.g.
// C.setAsMatlabVariable("HarmonicConvex.C");
int MatlabInterface::engPutVariable(const char *name, const mxArray *arr)
{
	int ret = 0;

	string s = name;
	int pos = s.find_first_of('.');
	if ( pos > 0 ) {
		// set a struct member
		char struct_name[255], member_name[255];
		int len;
		len = s.copy(struct_name, pos, 0);
		struct_name[len] = 0;
		len = s.copy(member_name, s.length() - pos, pos+1);
		member_name[len] = 0;
		//cout << struct_name << "->" << member_name << endl;
		ret = engPutVariable(m_ep, "TEMPORARY_VAR", arr);
		char cmd[1024];
		sprintf(cmd, "%s.%s = TEMPORARY_VAR; clear TEMPORARY_VAR;", struct_name, member_name);
		//cout << cmd << endl;
		EvalToCout(cmd, false);
	} else
		ret = engPutVariable(m_ep, name, arr);

	return ret;
}

mxArray *MatlabInterface::engGetVariable(Engine *ep, const char *name)
{
	mxArray *arr = 0;

	arr = ::engGetVariable(m_ep, name);

	int retries = 5; // doesn't help
	while ( arr == NULL ) {
		ERROR("!! Error in engGetVariable()" << name);
		if ( --retries < 0 )
			break;
		ERROR("Resetting the engine and trying again " << retries);
		GetEngine(true);
		arr = ::engGetVariable(m_ep, name);
	}

	return arr;
}

mxArray *MatlabInterface::engGetVariable(const char *name)
{
	mxArray *arr = 0;

	string s = name;
	int pos = s.find_first_of('.');
	if ( pos > 0 ) {
		// set a struct member
		char struct_name[255], member_name[255];
		int len;
		len = s.copy(struct_name, pos, 0);
		struct_name[len] = 0;
		len = s.copy(member_name, s.length() - pos, pos+1);
		member_name[len] = 0;
		//cout << struct_name << "->" << member_name << endl;
		char cmd[1024];
		sprintf(cmd, "TEMPORARY_VAR = %s.%s;", struct_name, member_name);
		//cout << cmd << endl;
		EvalToCout(cmd, false);
		arr = ::engGetVariable(m_ep, "TEMPORARY_VAR");
		EvalToCout("clear TEMPORARY_VAR;", false);
	} else
		arr = ::engGetVariable(m_ep, name);

	return arr;
}
