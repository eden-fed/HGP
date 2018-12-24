#pragma once
#include "CGAL\CGAL_Mesh.h"

class Borders
{
public:
	
	Borders(Mesh& mesh);
	int numBorders() const;
	int numVertices(unsigned int index) const;
	Vertex_handle vertex(unsigned int borderIndex, unsigned int vertexIndex) const;
	double borderCurvature(unsigned int index, bool useTargetMetric);
	static double getBorderHalfEdgePath(const Vertex_handle& source, const Vertex_handle& target, std::vector<Halfedge_handle>& path);
	static double getBorderHalfEdges(const Vertex_handle& vertex, std::vector<Halfedge_handle>& path);
	static double getPath3DLength(const std::vector<Halfedge_handle>& path);
	int genus();

private:

	MeshAdaptor mAdaptor;
	MeshFeatureExtractor mFeatureExtractor;
	std::vector<std::vector<Vertex_handle> > mBorders;
};
