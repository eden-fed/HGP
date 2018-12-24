/**
 * A higher-level interface to the matlab engine.
 * It also supports the use of a matlab singleton object
 */

#pragma once

#include <complex>
#include <cassert>
#include <map>
#include <string>
#include <vector>


// from Matlab
struct engine;
struct mxArray_tag;
typedef struct engine Engine;
typedef struct mxArray_tag mxArray;
typedef unsigned long	DWORD;


class MatlabInterface
{
public:

	//#define MAX_OUTPUT_BUFFER_SIZE (65536)
	#define MAX_OUTPUT_BUFFER_SIZE (262144)

	~MatlabInterface();


	// Using a SINGLETON object throughout the application
	// Get a singleton MatlabInterface object
	static MatlabInterface &GetEngine(bool restart = false);
	void EngineOpen();
	void EngineSetupUI();
	void EngineClose();

    static void DestroyMatrix(mxArray *&M);

    // Creates a matrix in matlab.
    template <typename T>
    int SetEngineIndexMatrix(const char *name, unsigned int m, unsigned int n, const T *vals, bool colmaj=false); // shifts indices by 1 (C/C++ vs matlab)
    template <typename T>
    int SetEngineRealMatrix(const char *name, unsigned int m, unsigned int n, const T *vals, bool colmaj=false);
    template <typename T>
    int SetEngineComplexMatrix(const char *name, unsigned int m, unsigned int n, const std::complex<T> *vals, bool colmaj=false);

    template <typename IndexType, typename ValueType>
    int SetEngineEncodedSparseRealMatrix(const char *name, unsigned int n,
            const IndexType *rowind, const IndexType *colind, const ValueType *vals);

    template <typename IndexType, typename ValueType>
    int SetEngineSparseRealMatrix(const char *name, unsigned int n,
            const IndexType *rowind, const IndexType *colind, const ValueType *vals, unsigned int nrows=0, unsigned int ncols=0);

    template <typename IndexType, typename ValueType>
    int SetEngineEncodedSparseComplexMatrix(const char *name, unsigned int n,
            const IndexType *rowind, const IndexType *colind, const std::complex<ValueType> *vals);

    template <typename IndexType, typename ValueType>
    int SetEngineSparseComplexMatrix(const char *name, unsigned int n,
            const IndexType *rowind, const IndexType *colind, const std::complex<ValueType> *vals, unsigned int nrows=0, unsigned int ncols=0);

	int CreateAllZerosSparseMatrix(const char* name, int nRows, int nCols);


    // Create a string in matlab.
    int SetEngineStringArray(const char *name, unsigned int n, const char **val);

    // Reads a matrix from matlab.
    template <typename T>
    int GetEngineIndexMatrix(const char *name, unsigned int m, unsigned int n, T *vals, bool colmaj=false); // shifts indices by 1 (C/C++ vs matlab)
    template <typename T>
    int GetEngineRealMatrix(const char *name, unsigned int m, unsigned int n, T *vals, bool colmaj=false);
    template <typename T>
    int GetEngineComplexMatrix(const char *name, unsigned int m, unsigned int n, std::complex<T> *vals, bool colmaj=false);

    template <typename T>
    int GetEngineIndexMatrix(const char *name, unsigned int &m, unsigned int &n, std::vector< T > &vals, bool colmaj=false); // shifts indices by 1 (C/C++ vs matlab)
    template <typename T>
    int GetEngineRealMatrix(const char *name, unsigned int &m, unsigned int &n, std::vector< T > &vals, bool colmaj=false);
    template <typename T>
    int GetEngineComplexMatrix(const char *name, unsigned int &m, unsigned int &n, std::vector< std::complex<T> > &vals, bool colmaj=false);

    void GetEncodedSparseRealMatrix(const char* name, unsigned int*& rowind, unsigned int*& colind, double*& vals, unsigned int& nentries);
	
	int GetEngineEncodedSparseRealMatrix(const char* name, std::vector<unsigned int>& Ir, std::vector<unsigned int>& Jc, std::vector<double>& vals, unsigned int& m, unsigned int& n);
	int GetEngineEncodedSparseComplexMatrix(const char* name, std::vector<unsigned int>& Ir, std::vector<unsigned int>& Jc, std::vector<std::complex<double> >& vals, unsigned int& m, unsigned int& n);

	int GetSparseRealMatrix(const char* name, std::vector<unsigned int>& rowind, std::vector<unsigned int>& colind, std::vector<double>& vals, unsigned int& nentries, unsigned int& m, unsigned int& n);
	int GetSparseComplexMatrix(const char* name, std::vector<unsigned int>& rowind, std::vector<unsigned int>& colind, std::vector<std::complex<double> >& vals, unsigned int& nentries, unsigned int& m, unsigned int& n);

	// Eval in-place string
	int Eval(const char *matlab_code);
    //int Eval(const char *matlab_code, char *output_buffer = NULL, int buffer_size = 0);
	std::string EvalToString(const char *matlab_code, bool saveCommandString = true);
	std::string EvalToString(const char *matlab_code, int &res, bool saveCommandString = true);

	//Eval and send the Matlab output to std::cout
	void EvalToCout(const char *matlab_code, bool saveCommandString = true);

    // Add a path to the matlab directory
    int AddScriptPath(const char* path);

	bool GetMatrixDimensions(const char* variableName, unsigned int& m, unsigned int& n);

	int engPutVariable(const char *name, const mxArray *arr);
	int engPutVariable(Engine *ep, const char *name, const mxArray *arr);
	mxArray *engGetVariable(const char *name);
	mxArray *engGetVariable(Engine *ep, const char *name);

private:
	MatlabInterface();
    // copying a MatlabInterface object is disallowed
    MatlabInterface(const MatlabInterface &);
    MatlabInterface &operator=(const MatlabInterface &);

	void quickRestart();

    // Note that the arrays created by these matrices *must be destroyed*
    // afterward.
    // matrix creation helper functions
    template <typename T>
    static mxArray* CreateIndexMatrix(unsigned int m, unsigned int n, const T *vals, bool colmaj=false); // shifts indices by 1 (C/C++ vs matlab)
    template <typename T>
    static mxArray* CreateRealMatrix(unsigned int m, unsigned int n, const T *vals, bool colmaj=false);
    template <typename T>
    static mxArray* CreateComplexMatrix(unsigned int m, unsigned int n, const std::complex<T> *vals, bool colmaj=false);

    template <typename IndexType, typename ValueType>
    static mxArray* CreateEncodedSparseRealMatrix(unsigned int n,
            const IndexType *rowind, const IndexType *colind, const ValueType *vals);

    template <typename IndexType, typename ValueType>
    static mxArray* CreateEncodedSparseComplexMatrix(unsigned int n,
            const IndexType *rowind, const IndexType *colind, const std::complex<ValueType> *vals);

    // string creation helper functions
    static mxArray* CreateStringArray(unsigned int n, const char **val);

    // matrix copy-back helper functions
    template <typename T>
    static int CopyFromIndexMatrix(mxArray *M, unsigned int m, unsigned int n, T *dest, bool colmaj=false); // shifts indices by 1 (C/C++ vs matlab)
    template <typename T>
    static int CopyFromRealMatrix(mxArray *M, unsigned int m, unsigned int n, T *dest, bool colmaj=false);
    template <typename T>
    static int CopyFromComplexMatrix(mxArray *M, unsigned int m, unsigned int n, std::complex<T> *dest, bool colmaj=false);


private:
    Engine *m_ep;
	char* mOutputStringBuffer;
	static DWORD mThreadIdUsedToOpenEngine;
};
