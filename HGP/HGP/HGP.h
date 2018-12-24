#pragma once

#include "HarmonicParametrization.h"
#include "CGAL\CGAL_Mesh.h"
#include "Utils/GMM_Macros.h"
#include "Parser/Parser.h"
#include <complex>


class HGP:HarmonicParametrization
{
public:
	HGP(){ mMethodName = "HGP"; mHasCones = true; }
	~HGP(){}

	bool run(std::string& objPath, std::string& vfPath);


private:
	//-------------Variables-------------------
	int N;
	GMMDenseColMatrix F, halfEdges;
	std::vector<int> seamHalfEdgesIndices;
	bool useFrameFixing,frameStatus;
	bool calcFramesFromVecField;

	//-----------------------------------------

	//~~~~~~~~~~~~~~~~Methods~~~~~~~~~~~~~~~~~~
	bool loadMesh(std::string& objPath, std::string& vfPath);
	void transferMeshToMatlab();
	void setHalfEdgesMap();
	void setRotationsConstraints();
	void setBoundaryFaces();
	void updateHalfEdgesMetric(bool firstTime);
	void setHarmonicInternalVerticesConstraints();
	void setHarmonicSeamVerticesConstraints();
	void setFramesInMatlab(bool firstTime);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
};

typedef complex<double> Complex;
