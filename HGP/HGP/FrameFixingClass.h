#pragma once
#include "CGAL/CGAL_Mesh.h"
#include "Utils/GMM_Macros.h"
#include "Parser/Parser.h"


class FrameFixingClass
{
private:
	Mesh& cgalMesh;
	MeshBuffer& meshAdditionalData;
	GMMDenseColMatrix& F;
	std::vector<int> badCones;
	std::vector<int> oneRingFaces;
	std::vector<double> newAngles;
	Halfedge_around_vertex_circulator start;

	//private methods
	void searchForBadCone();
	void extractOneRingAngles();
	void setLocalEmbedding();

public:
	FrameFixingClass(Mesh& cgal_Mesh, MeshBuffer& m1, GMMDenseColMatrix& Faces) : cgalMesh(cgal_Mesh), meshAdditionalData(m1), F(Faces) {}
	~FrameFixingClass() {}
	bool runFrameFixingProcedure(bool useFrameFixing);
};

typedef complex<double> Complex;