#pragma once

#include "CGAL\CGAL_Mesh.h"
#include "Utils/GMM_Macros.h"
#include "Parser/Parser.h"

class HarmonicParametrization
{
public:
	HarmonicParametrization(){}
	~HarmonicParametrization(){}

	virtual bool run(std::string& objpath, std::string& vfPath) = 0;


protected:
	//-------------Variables-------------
	std::string mMethodName;
	Mesh mCgalMesh;
	MeshBuffer mMeshBuffer;
	std::vector<int> mIndexOfHEinSystem;//this array maps from halfedges indices to the index in the system (the variables in the system are the vertices + duplications of vertices on the cut)
	bool mHasCones;
	int mSizeOfSystemVar;
	GMMDenseColMatrix mRotationConstraints, mUVs;


	//-------------Methods-------------
	void buildTriangleMesh(const std::vector<Point_3>& vertices, const std::vector<unsigned int>& faceIndices, Mesh& mesh);
	void setHalfEdgesInCgal();
	void updateIndexOfHEinSystem();
	void updatHalfedgeUVs();
	void checkForFoldovers(std::vector<Facet_handle>& flippedTriangles);
	void visualize();
	int nextOnSeam(int h);
	int degreeOnSeam(int h);
	void calcDistortion();
	double calcK(Facet_const_handle& face);
	void coneAngleDetection(int& numWrongAngles, int& numWrongConeAngles);
	double oneRingAngle(Vertex_iterator& v);
};

typedef complex<double> Complex;
