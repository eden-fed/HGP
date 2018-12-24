
#include "Borders.h"




Borders::Borders(Mesh& mesh) : mAdaptor(mesh), mFeatureExtractor(mAdaptor)
{
	const Skeleton& skeleton = mFeatureExtractor.get_borders();

	int numBorders = skeleton.size();

	mBorders.resize(numBorders);

	for(int i = 0; i < numBorders; i++)
	{
		Border* border = skeleton[i];
		int borderSize = border->size();

		std::vector<Vertex_handle>& currBorder = mBorders[i];
		currBorder.reserve(borderSize);

		for(Border::iterator iter = border->begin(); iter != border->end(); iter++)
		{
			currBorder.push_back(*iter);
		}
	}
}




int Borders::numBorders() const
{
	return mBorders.size();
}

int Borders::genus()
{
	return mFeatureExtractor.get_genus();
}


int Borders::numVertices(unsigned int index) const
{
	assert(index < mBorders.size());
	return mBorders[index].size();
}



Vertex_handle Borders::vertex(unsigned int borderIndex, unsigned int vertexIndex) const
{
	assert(borderIndex < mBorders.size());
	assert(vertexIndex < mBorders[borderIndex].size());

	return mBorders[borderIndex][vertexIndex];
}



///////////////////////////////////////////////////////////
// computes the Geodesic curvature of the border curve //
///////////////////////////////////////////////////////////
// index - which border to use                         //
// useTargetMetric - whether to use the metric from the  //
//                   3D embedding or from the prescribed //
//                   edge lengths                        //
///////////////////////////////////////////////////////////
double Borders::borderCurvature(unsigned int index, bool useTargetMetric)
{
	assert(index < mBorders.size());

	double totalCurvature = 0.0;

	int numV = numVertices(index);

	for(int i = 0; i < numV; i++)
	{
		Vertex_handle v = vertex(index, i);
		assert(v->is_border());

		double kappa = 0.0;
		
		if(useTargetMetric)
		{
			kappa = v->targetMetricGaussianCurvature();
		}
		else
		{
			kappa = v->gaussianCurvature();
		}
		
		totalCurvature += kappa;
	}
	return totalCurvature;
}




//finds a list of border halfedges connecting "source" to "target"
//returns the uv length of the polyline defined by this path and 0.0 if there is no path connecting the two vertices
double Borders::getBorderHalfEdgePath(const Vertex_handle& source, const Vertex_handle& target, std::vector<Halfedge_handle>& path)
{
	double length = 0.0;

	if(!source->is_border() || !target->is_border())
	{
		return 0.0;
	}

	Halfedge_around_vertex_circulator hec = source->vertex_begin();
	Halfedge_around_vertex_circulator hec_end = hec;

	CGAL_For_all(hec, hec_end)
	{
		if(hec->opposite()->is_border()) break;
	}
	
	path.clear();

	bool connected = false;

	for(Halfedge_handle he = hec->opposite(); he->vertex() != source; he = he->next())
	{
		path.push_back(he);
		length += he->uvLength();

		if(he->vertex() == target)
		{
			connected = true;
			break;
		}
	}
	
	if(connected)
	{
		return length;
	}
	else
	{
		return 0.0; //there is no border path from source to target
	}
}



double Borders::getPath3DLength(const std::vector<Halfedge_handle>& path)
{
	double pathLength = 0.0;
	for(int i = 0; i < path.size(); i++)
	{
		double edgeLength = path[i]->length(); //edge length from 3D embedding
		pathLength += edgeLength;
	}
	return pathLength;
}



//find the halfedges of the border that the vertex belongs to
//returns the uv length of the polyline border or 0.0 if such a border is not found
double Borders::getBorderHalfEdges(const Vertex_handle& vertex, std::vector<Halfedge_handle>& path)
{
	double length = 0.0;

	if(!vertex->is_border())
	{
		return 0.0;
	}

	Halfedge_around_vertex_circulator hec = vertex->vertex_begin();
	Halfedge_around_vertex_circulator hec_end = hec;

	CGAL_For_all(hec, hec_end)
	{
		if(hec->opposite()->is_border()) break;
	}

	path.clear();

	Halfedge_handle he = hec->opposite();

	for(; he->vertex() != vertex; he = he->next())
	{
		path.push_back(he);
		length += he->uvLength();
	}
	path.push_back(he); //add the last one
	length += he->uvLength();

	return length;
}


