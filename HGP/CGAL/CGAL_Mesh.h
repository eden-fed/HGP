#pragma once

#include "CGAL_Vertex.h"
#include "CGAL_Halfedge.h"
#include "CGAL_Face.h"
#include "CGAL\Polyhedron_3.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL\Parameterization_polyhedron_adaptor_3.h>
#include <CGAL\Parameterization_mesh_feature_extractor.h>

template<class kernel, class items>
class Polyhedron : public CGAL::Polyhedron_3<kernel, items>
{
public:

	Polyhedron() {}
	virtual ~Polyhedron() {}

	void updateAllGlobalIndices();
	bool updateVerticesLocalIndices();
	bool updateVerticesXYLocalIndices(int& numFree, int& numDependent);
	bool updateFacesLocalIndices();
	int countVertexStates(int& numFreeVertices, int& numFixedVertices) const;
	int fixBorderVertices(); //set all border vertices to be fixed and all the rest to be free. return the number of border vertices
	int fixBorderFaces(); //set all border faces to be fixed and all the rest to be free. return the number of border faces
	bool getBorderVerticesCCW(std::vector<Vertex_handle>& borderVertices, const Vertex_const_handle* firstVertex = NULL);
	bool getBorderHalfEdges(std::vector<Halfedge_handle>& borderHalfedges);
	void getFixedVertexIndexMap(std::vector<Vertex_handle>& map);
	void getFixedVertexIndexMap(std::vector<Vertex_const_handle>& map) const;
	void getFreeVertexIndexMap(std::vector<Vertex_handle>& map);
	void getFreeVertexIndexMap(std::vector<Vertex_const_handle>& map) const;
	bool hasValidMetric() const;
	void setFreeAllVertices();
	void updateTargetMetricFromConformalFactors(bool useTargetMetric);
	void updateTargetMetricFromUVs();
	void updateTargetMetricFrom3dEmbedding();
	double averageEdgeLength();
	double averageTargetMetricEdgeLength();
	double averageUVEdgeLength();
	Facet_handle face(unsigned int index);
	Facet_const_handle face(unsigned int index) const;
	Vertex_handle vertex(unsigned int index);
	Vertex_const_handle vertex(unsigned int index) const;
	double area() const;
	double eulerNumber(bool useTargetMetric = false) const;
	void setName(const std::string& name);
	std::string name() const;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Halfedge_handle halfEdge(unsigned int index);
	Halfedge_const_handle halfEdge(unsigned int index) const;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

private:

	std::vector<Facet_handle> mFaceIndexMap;
	std::vector<Vertex_handle> mVertexIndexMap;
	std::string mName; //name for the mesh
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	std::vector<Halfedge_handle> mHalfEdgeIndexMap;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
};




template<class kernel, class items>
void Polyhedron<kernel, items>::setName(const std::string& name)
{
	mName = name;
}


template<class kernel, class items>
std::string Polyhedron<kernel, items>::name() const
{
	return mName;
}


template<class kernel, class items>
typename Polyhedron<kernel, items>::Facet_handle
Polyhedron<kernel, items>::face(unsigned int index)
{
	assert(index >= 0 && index < mFaceIndexMap.size());
	return mFaceIndexMap[index];
}



template<class kernel, class items>
typename Polyhedron<kernel, items>::Facet_const_handle
Polyhedron<kernel, items>::face(unsigned int index) const
{
	assert(index >= 0 && index < mFaceIndexMap.size());
	return mFaceIndexMap[index];
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template<class kernel, class items>
typename Polyhedron<kernel, items>::Halfedge_handle
Polyhedron<kernel, items>::halfEdge(unsigned int index)
{
	assert(index >= 0 && index < mHalfEdgeIndexMap.size());
	return mHalfEdgeIndexMap[index];
}



template<class kernel, class items>
typename Polyhedron<kernel, items>::Halfedge_const_handle
Polyhedron<kernel, items>::halfEdge(unsigned int index) const
{
	assert(index >= 0 && index < mHalfEdgeIndexMap.size());
	return mHalfEdgeIndexMap[index];
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<class kernel, class items>
typename Polyhedron<kernel, items>::Vertex_handle
Polyhedron<kernel, items>::vertex(unsigned int index)
{
	assert(index >= 0 && index < mVertexIndexMap.size());
	return mVertexIndexMap[index];
}



template<class kernel, class items>
typename Polyhedron<kernel, items>::Vertex_const_handle
Polyhedron<kernel, items>::vertex(unsigned int index) const
{
	assert(index >= 0 && index < mVertexIndexMap.size());
	return mVertexIndexMap[index];
}


template<class kernel, class items>
bool Polyhedron<kernel, items>::hasValidMetric() const
{
	for(Face_const_iterator f = facets_begin(); f != facets_end(); f++)
	{
		if(!f->hasValidMetric())
		{
			return false;
		}
	}
	return true;
}



template<class kernel, class items>
void Polyhedron<kernel, items>::getFixedVertexIndexMap(std::vector<Vertex_handle>& map)
{
	int n = size_of_vertices();

	std::vector<Vertex_handle> fixedVertices;
	fixedVertices.reserve(n);

	for(Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi)
	{
		if(!vi->isFree())
		{
			fixedVertices.push_back(vi);
		}
	}
	int numFixedVertices = fixedVertices.size();
	map.resize(numFixedVertices, NULL);

	for(int i = 0; i < numFixedVertices; i++)
	{
		map[fixedVertices[i]->localIndex()] = fixedVertices[i];
	}
}



template<class kernel, class items>
void Polyhedron<kernel, items>::getFixedVertexIndexMap(std::vector<Vertex_const_handle>& map) const
{
	int n = size_of_vertices();

	std::vector<Vertex_const_handle> fixedVertices;
	fixedVertices.reserve(n);

	for(Vertex_const_iterator vi = vertices_begin(); vi != vertices_end(); ++vi)
	{
		if(vi->state() == Vertex::Base::FIXED_VERTEX)
		{
			fixedVertices.push_back(vi);
		}
	}
	int numFixedVertices = fixedVertices.size();
	map.resize(numFixedVertices, NULL);

	for(int i = 0; i < numFixedVertices; i++)
	{
		map[fixedVertices[i]->localIndex()] = fixedVertices[i];
	}
}



template<class kernel, class items>
void Polyhedron<kernel, items>::getFreeVertexIndexMap(std::vector<Vertex_handle>& map)
{
	int n = size_of_vertices();

	std::vector<Vertex_handle> freeVertices;
	freeVertices.reserve(n);

	for(Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi)
	{
		if(vi->isFree())
		{
			freeVertices.push_back(vi);
		}
	}
	int numFreeVertices = freeVertices.size();
	map.resize(numFreeVertices, NULL);

	for(int i = 0; i < numFreeVertices; i++)
	{
		map[freeVertices[i]->localIndex()] = freeVertices[i];
	}
}



template<class kernel, class items>
void Polyhedron<kernel, items>::getFreeVertexIndexMap(std::vector<Vertex_const_handle>& map) const
{
	int n = size_of_vertices();

	std::vector<Vertex_const_handle> freeVertices;
	freeVertices.reserve(n);

	for(Vertex_const_iterator vi = vertices_begin(); vi != vertices_end(); ++vi)
	{
		if(vi->state() == Vertex::Base::FREE_VERTEX)
		{
			freeVertices.push_back(vi);
		}
	}
	int numFreeVertices = freeVertices.size();
	map.resize(numFreeVertices, NULL);

	for(int i = 0; i < numFreeVertices; i++)
	{
		map[freeVertices[i]->localIndex()] = freeVertices[i];
	}
}



template<class kernel, class items>
int Polyhedron<kernel, items>::fixBorderVertices()
{
	int numBorderVertices = 0;

	for(Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi)
	{
		if(vi->is_border())
		{
			vi->isFree() = false;
			numBorderVertices++;
		}
		else
		{
			vi->isFree() = true;
		}
	}
	bool res = updateVerticesLocalIndices();
	if(!res)
	{
		return -1;
	}
	return numBorderVertices;
}


//set all border faces to be fixed and all the rest to be free. return the number of border faces
//note that a border face is a face which has edge on the border (vertex on the border alone is not enough)
template<class kernel, class items>
int Polyhedron<kernel, items>::fixBorderFaces()
{
	int numBorderFaces = 0;

	for(Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi)
	{
		if(fi->is_border())
		{
			fi->state() = Face::Base::FIXED_FACE;
			numBorderVertices++;
		}
		else
		{
			fi->state() = Face::Base::FREE_FACE;
		}
	}
	bool res = updateVerticesLocalIndices();
	if(!res)
	{
		return -1;
	}
	return numBorderFaces;
}



//get a list of all the border vertices ordered counter clockwise. If firstVertex is provided (not NULL) then the function
//will arrange the border vertices such that the first element in the list will be that vertex.
//first vertex is optional. however, if it is provided it must point to a border vertex
template<class kernel, class items>
bool Polyhedron<kernel, items>::getBorderVerticesCCW(std::vector<Vertex_handle>& borderVertices, const Vertex_const_handle* firstVertex)
{
	int k = size_of_border_halfedges();

	if(k <= 3)
	{
		return false;
	}

	borderVertices.clear();
	borderVertices.resize(k);

	Halfedge_iterator borderEdge = border_edges_begin();
	if(!borderEdge->is_border()) borderEdge =  borderEdge->opposite();

	int firstVertexIndex = 0;
	for(int i = 0; i < k; i++, borderEdge = borderEdge->prev())
	{
		Vertex_handle v = borderEdge->vertex();
		borderVertices[i] = v;
		if(firstVertex && (Vertex_const_handle)v == (*firstVertex))
		{
			firstVertexIndex = i;
		}
	}

	if(firstVertex == NULL)
	{
		return true;
	}
	if(!(*firstVertex)->is_border())
	{
		return false;
	}
	else
	{
		std::rotate(borderVertices.begin(), borderVertices.begin() + firstVertexIndex, borderVertices.end());
		return true;
	}
}


//This function assumes that the mesh has exactly one border.
//The halfedges are ordered in clockwise manner
template<class kernel, class items>
bool Polyhedron<kernel, items>::getBorderHalfEdges(std::vector<Halfedge_handle>& borderHalfedges)
{
	int k = size_of_border_halfedges();

	if(k <= 3)
	{
		assert(0);
		return false;
	}

	borderHalfedges.clear();
	borderHalfedges.resize(k);

	Halfedge_iterator borderHalfedge = border_edges_begin();
	if(!borderHalfedge->is_border()) borderHalfedge =  borderHalfedge->opposite();

	Halfedge_around_facet_circulator hc = borderHalfedge->facet_begin();
	const Halfedge_around_facet_circulator hcEnd = hc;

	int i = 0;

	CGAL_For_all(hc, hcEnd)
	{
		borderHalfedges[i] = hc;
		i++;
	}
	if(i != k)
	{
		assert(0); //this should mean that the mesh has more than one boundary
		return false;
	}
	
	return true;
}


template<class kernel, class items>
double Polyhedron<kernel, items>::area() const
{
	double sumOfAllFaceAreas = 0.0;

	Facet_const_handle f = facets_begin();
	const Facet_const_handle fEnd = facets_end();

	CGAL_For_all(f, fEnd)
	{
		sumOfAllFaceAreas += f->area();
	}

	return sumOfAllFaceAreas;
}


//uses Gauss–Bonnet theorem to compute the Euler characteristic number of a mesh (with or without boundary).
//the function sums the integrated Gaussian curvature (angle deficit without area normalization).
//the input parameter useTargetMetric determines whether the curvature is being computed based on the 3D embedding or based on the prescribed target metric (edge lengths).
template<class kernel, class items>
double Polyhedron<kernel, items>::eulerNumber(bool useTargetMetric = false) const
{
	double sumOfAngles = 0.0;

	for(Vertex_const_iterator vi = vertices_begin(); vi != vertices_end(); ++vi)
	{
		double angleDeficit = 0.0;
		
		if(useTargetMetric)
		{
			angleDeficit = vi->targetMetricGaussianCurvature(false);
		}
		else
		{
			angleDeficit = vi->gaussianCurvature(false);
		}
		sumOfAngles += angleDeficit;
	}

	double eulerNumber = sumOfAngles / (2.0*CGAL_PI);
	return eulerNumber;
}



template<class kernel, class items>
double Polyhedron<kernel, items>::averageEdgeLength()
{
	int numEdges = size_of_halfedges() / 2;

	double average = 0.0;

	for(Edge_iterator e = edges_begin(); e != edges_end(); e++)
	{
		average += e->length();
	}

	average /= numEdges;

	return average;
}



template<class kernel, class items>
double Polyhedron<kernel, items>::averageTargetMetricEdgeLength()
{
	int numEdges = size_of_halfedges() / 2;

	double average = 0.0;

	for(Edge_iterator e = edges_begin(); e != edges_end(); e++)
	{
		average += e->targetMetric();
	}

	average /= numEdges;

	return average;
}



template<class kernel, class items>
double Polyhedron<kernel, items>::averageUVEdgeLength()
{
	int numEdges = size_of_halfedges() / 2;

	double average = 0.0;

	for(Edge_iterator e = edges_begin(); e != edges_end(); e++)
	{
		average += e->uvLength();
	}

	average /= numEdges;

	return average;
}




template<class kernel, class items>
void Polyhedron<kernel, items>::setFreeAllVertices()
{
	Vertex_iterator vi = vertices_begin();

	for(int i = 0; vi != vertices_end(); vi++, i++)
	{
		vi->isFree() = true;
		vi->isXfree() = true;
		vi->isYfree() = true;
	}
}

template<class kernel, class items>
void Polyhedron<kernel, items>::updateTargetMetricFromConformalFactors(bool useTargetMetric)
{
	for(Halfedge_handle h = halfedges_begin(); h != halfedges_end(); h++)
	{
		Halfedge_handle he = h;
		if(h->is_border())
		{
			he = h->opposite();
		}
		double oldEdgeLength = 0.0;
		if(useTargetMetric)
		{
			oldEdgeLength = he->targetMetric();
		}
		else
		{
			oldEdgeLength = he->length();
		}
		double newEdgeLength = he->scaleFactor()*oldEdgeLength;
		assert(newEdgeLength > 0.0);

		h->targetMetric() = newEdgeLength;
	}
}


//update the target metric with edge lengths taken from the parameterization
template<class kernel, class items>
void Polyhedron<kernel, items>::updateTargetMetricFromUVs()
{
	for(Halfedge_handle h = halfedges_begin(); h != halfedges_end(); h++)
	{
		std::complex<double> v1 = h->vertex()->uvC();
		std::complex<double> v2 = h->opposite()->vertex()->uvC();
		double edgeLength = std::abs(v1-v2);
		h->targetMetric() = edgeLength;
	}
}


//update the target metric with edge lengths taken from the 3D embedding
template<class kernel, class items>
void Polyhedron<kernel, items>::updateTargetMetricFrom3dEmbedding()
{
	for(Halfedge_handle h = halfedges_begin(); h != halfedges_end(); h++)
	{
		h->targetMetric() = h->length();
	}
}


template<class kernel, class items>
void Polyhedron<kernel, items>::updateAllGlobalIndices()
{
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	int numHalfEdges = size_of_halfedges();
	assert(numHalfEdges > 0);
	mHalfEdgeIndexMap.resize(numHalfEdges);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Halfedge_handle he = halfedges_begin();

	for(int j = 0; he != halfedges_end(); he++, j++)
	{
		he->index() = j;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		mHalfEdgeIndexMap[j] = he;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}

	int numVertices = size_of_vertices();
	assert(numVertices > 0);

	mVertexIndexMap.resize(numVertices);

	Vertex_iterator v = vertices_begin();

	for(int j = 0; v != vertices_end(); v++, j++)
	{
		v->index() = j;
		mVertexIndexMap[j] = v;
	}


	int numFaces = size_of_facets();
	assert(numFaces > 0);

	mFaceIndexMap.resize(numFaces);
	
	Facet_handle f = facets_begin();

	for(int j = 0; f != facets_end(); f++, j++)
	{
		f->index() = j;
		mFaceIndexMap[j] = f;
	}
}


template<class kernel, class items>
bool Polyhedron<kernel, items>::updateVerticesLocalIndices()
{
	int numFree = 0;
	int numFixed = 0;
	for(Vertex_iterator vi = vertices_begin(); vi != vertices_end(); vi++)
	{
		if(vi->isFree())
		{
			vi->localIndex() = numFree;
			numFree++;
		}
		else if(!vi->isFree())
		{
			vi->localIndex() = numFixed;
			numFixed++;
		}
	}

	if(numFixed + numFree == size_of_vertices())
	{
		return true;
	}
	else
	{
		return false;
	}
}

//numFree and numDependent are both input and output
//the counting always starts from the given values and the value of the variables is then further updated
//this is useful when you want to work with multiple meshes and you want to have only one pool of indices
template<class kernel, class items>
bool Polyhedron<kernel, items>::updateVerticesXYLocalIndices(int& numFree, int& numDependent)
{
	for(Vertex_iterator vi = vertices_begin(); vi != vertices_end(); vi++)
	{
		if(vi->isXfree())
		{
			vi->xLocalIndex() = numFree;
			numFree++;
		}
		else
		{
			vi->xLocalIndex() = numDependent;
			numDependent++;
		}
		if(vi->isYfree())
		{
			vi->yLocalIndex() = numFree;
			numFree++;
		}
		else
		{
			vi->yLocalIndex() = numDependent;
			numDependent++;
		}
	}

	if(numDependent + numFree == 2*size_of_vertices())
	{
		return true;
	}
	else
	{
		return false;
	}
}

template<class kernel, class items>
bool Polyhedron<kernel, items>::updateFacesLocalIndices()
{
	int numFree = 0;
	int numFixed = 0;
	for(Facet_iterator fi = facets_begin(); fi != facets_end(); fi++)
	{
		if(fi->state() == FBase::FREE_FACE)
		{
			fi->localIndex() = numFree;
			numFree++;
		}
		else if(fi->state() == FBase::FIXED_FACE)
		{
			fi->localIndex() = numFixed;
			numFixed++;
		}
	}

	if(numFixed + numFree == size_of_faces())
	{
		return true;
	}
	else
	{
		return false;
	}
}


template<class kernel, class items>
int Polyhedron<kernel, items>::countVertexStates(int& numFreeVertices, int& numFixedVertices) const
{
	numFreeVertices = 0;
	numFixedVertices = 0;

	for(Vertex_const_iterator vi = vertices_begin(); vi != vertices_end(); ++vi)
	{
		if(vi->isFree())
		{
			numFreeVertices++;
		}
		else if(!vi->isFree())
		{
			numFixedVertices++;
		}
	}
	int totalNumVertices = size_of_vertices();
	if(totalNumVertices != numFixedVertices + numFreeVertices)
	{
		cerr << "The number of fixed and free vertices is not consistent with the total number of vertices\n";
	}

	return totalNumVertices;
}


template<class H, class K>
class MeshBuilder : public CGAL::Modifier_base<H>
{
	typedef typename K::Point_3 Point_3;
	typedef typename H::Vertex_handle Vertex_handle;

public:

	MeshBuilder(const std::vector<Point_3>* vertices, const std::vector<int>* faceIndices, const std::vector<int>* vertexCount = NULL); //set vertexCount to NULL for triangle meshes
	void operator()(H& hds);

protected:

	MeshBuilder(){}; //shouldn't be used

protected:

	const std::vector<Point_3>* mVertices;
	const std::vector<int>* mFaceIndices;
	const std::vector<int>* mVertexCount;
};  



template<class H, class K>
MeshBuilder<H, K>::MeshBuilder(const std::vector<Point_3>* vertices, const std::vector<int>* faceIndices, const std::vector<int>* vertexCount) :
mVertices(vertices), mFaceIndices(faceIndices), mVertexCount(vertexCount)
{
}


template<class H, class K>
void MeshBuilder<H, K>::operator()(H& hds)
{
	int numVertices = mVertices->size();
	int numFaces = 0;
	
	if(mVertexCount != NULL) //polygon mesh
	{
		numFaces = mVertexCount->size();
	}
	else //triangle mesh
	{
		assert(mFaceIndices->size() % 3 == 0);
		numFaces = mFaceIndices->size() / 3;
	}

	if(numFaces < 1 || numVertices < 3)
	{
		std::cerr << "Can't create CGAL HDS from face's list!\n";
		return;
	}

	CGAL::Polyhedron_incremental_builder_3<H> builder(hds, true);

	builder.begin_surface(numVertices, numFaces);

	for(int i = 0; i < numVertices; i++)
	{
		Vertex_handle v = builder.add_vertex((*mVertices)[i]);
	}

	int index = 0;

	for(int i = 0; i < numFaces; i++)
	{
		int numVerticesInCurrFace = 0;
		if(mVertexCount != NULL) //polygon mesh
		{
			numVerticesInCurrFace = (*mVertexCount)[i];
		}
		else
		{
			numVerticesInCurrFace = 3; //triangle mesh
		}

		builder.begin_facet();

		for(int j = 0; j < numVerticesInCurrFace; j++)
		{
			assert((*mFaceIndices)[index] >= 0);
			builder.add_vertex_to_facet((*mFaceIndices)[index]);
			index++;
		}
		builder.end_facet();
	}

	builder.end_surface();
	hds.normalize_border();
}



struct PolyhedronItems : public CGAL::Polyhedron_items_3
{
	template<class Refs, typename Traits>
	struct Vertex_wrapper
	{
		typedef Vertex_Base<Refs, Traits> Vertex;
	};

	template<class Refs, typename Traits>
	struct Halfedge_wrapper
	{
		typedef HalfedgeBase<Refs, Traits> Halfedge;
	};

	template<class Refs, typename Traits>
	struct Face_wrapper
	{
		typedef FaceBase<Refs, Traits> Face;
	};
};




typedef Polyhedron<CGAL::Simple_cartesian<double>, PolyhedronItems> Mesh;



typedef  Mesh::Traits::Kernel Kernel;
typedef  Kernel::Point_3 Point;
typedef  Kernel::Point_3 Point_3;
typedef  Kernel::Point_2 Point_2;
typedef  Kernel::Vector_3 Vector;
typedef  Kernel::Vector_3 Vector_3;

typedef  Mesh::HalfedgeDS HDS;
typedef  Mesh::Vertex_handle Vertex_handle;
typedef  Mesh::Vertex_const_handle Vertex_const_handle;
typedef  Mesh::Halfedge_handle Halfedge_handle;
typedef  Mesh::Halfedge_const_handle Halfedge_const_handle;
typedef  Mesh::Facet_handle Facet_handle;
typedef  Mesh::Facet_const_handle Facet_const_handle;

typedef  Mesh::Vertex_iterator Vertex_iterator;
typedef  Mesh::Vertex_const_iterator Vertex_const_iterator;
typedef  Mesh::Facet_iterator Facet_iterator;
typedef  Mesh::Facet_const_iterator Facet_const_iterator;
typedef  Mesh::Halfedge_iterator Halfedge_iterator;
typedef  Mesh::Halfedge_const_iterator Halfedge_const_iterator;
typedef  Mesh::Edge_iterator Edge_iterator;
typedef  Mesh::Edge_const_iterator Edge_const_iterator;

typedef  Mesh::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
typedef  Mesh::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
typedef  Mesh::Halfedge_around_facet_const_circulator Halfedge_around_facet_const_circulator;
typedef  Mesh::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;

typedef  CGAL::Parameterization_polyhedron_adaptor_3<Mesh> MeshAdaptor;
typedef  CGAL::Parameterization_mesh_feature_extractor<MeshAdaptor> MeshFeatureExtractor;

typedef  MeshFeatureExtractor::Skeleton Skeleton;
typedef  MeshFeatureExtractor::Border Border;