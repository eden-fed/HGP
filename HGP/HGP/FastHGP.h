#pragma once

#include "HarmonicParametrization.h"
#include "CGAL\CGAL_Mesh.h"
#include "Utils/GMM_Macros.h"
#include "Borders.h"
#include "Parser/Parser.h"
#include "Eigen/SparseCore"
#include <unordered_set>
#include <complex>


class FastHGP:HarmonicParametrization
{
public:
	FastHGP(){ mMethodName = "FastHGP"; }
	~FastHGP(){}

	bool run(std::string& objpath, std::string& vfPath);


private:
	//-------------Variables-------------------
	bool mCalcFramesFromVecField, mATPandNewtonPassed, mFixCot;
	int mNumDOF, mSizeOfMatrix, mSegSize, mNumOfBorderVertices;
	std::vector<Vertex_handle> mConesAndMetaMap;//all the cones and meta vertices 
	Eigen::SparseMatrix<double, Eigen::RowMajor> mKKtEigen;//KKT matrix
	GMMSparseComplexRowMatrix mHarmonicBasisInAllV;//Harmonic basis matrix in the size of mSizeOfSystemVar X mNumDOF - but only rows that are cones or near cones are filled (the others are just 0 rows)
	std::vector<int> mIndicesOfMetaVerticesInUVbyGeneralIndex;//Used for Tutte initial value
	std::vector<std::vector<Halfedge_handle> > mMetaVerticesInBorderByHalfEdge;//for each border - the halfedges corresponding to meta vertices in this border
	std::vector<bool> mIsMeta; //array in the size of all vertices that indicates if a vertex is meta 
	std::unordered_set<int> mConesAndNearConesMapOfRowsInKKT;//all the rows to calculate in the harmonic basis (use hash map to prevent duplications)
	double mTotalTime;


	//--------------Methods--------------
	bool initialize(Borders& borders);
	bool loadMesh(std::string& objpath, std::string& vfPath);
	void getSettings();
	void getConesMap();
	void getBordersMapAndSetMetaVertices(Borders& borders);
	void getVerticesThatMustBeMeta(Borders& borders);
	bool TransferBorderFacesAndVerticesToMatlab();
	bool runAlgorithm(Borders& borders);
	void constructKKTmatrix(Borders& borders, int& conesConstraintsStartRow);
	void FillLaplacianInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues);
	void updateTermInLaplacianByHalfEdge(Halfedge_handle h, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues);
	void updateTermInLaplacianByIndices(int i, int j, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues);
	void FillRotationConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT);
	void FillMetaVerticesConstraintsInKKT(Borders& borders, std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT);
	void FillConesConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow);
	void SetElementsForPARDISO(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow);
	bool constructHarmonicBasisAndSendToMATLAB(int conesConstraintsStartRow);
	bool calculateHarmonicBasisInPARDISO(int conesConstraintsStartRow);
	bool getATPInitialValue();
	bool getTutteInitialValue(Borders& borders);
	void solveSystemToFindInitialValueEntriesOnNonMainBorder(Borders& borders, int& numMetaVerticesMainBorder, GMMDenseColMatrix& initialValue);
	bool runNewton(int conesConstraintsStartRow, GMMDenseColMatrix& RHS);
	bool getAllUVsByRHS(GMMDenseColMatrix& RHS);

	bool testResult();
	void checkLocationOfFoldoversAndPrint(std::vector<Facet_handle>& flippedTriangles, int& numFoldovers, int& numFoldsNearCones, int& numFoldsNearBorder);
	void fixCotFoldovers(std::vector<Facet_handle>& flippedTriangles);
	void putVertexInKernelUsingCVX(Vertex_handle v);
	void sendValuesToMatlabReport();

};

typedef complex<double> Complex;
