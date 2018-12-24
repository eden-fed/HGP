#pragma once
#include "CGAL\HalfedgeDS_face_base.h"
#include "CGAL\circulator.h"

template <class Refs, typename Traits>
class FaceBase : public CGAL::HalfedgeDS_face_base<Refs, CGAL::Tag_true>
{
	typedef typename Traits::Kernel::Vector_3 Vector_3;
	typedef typename Traits::Kernel::Point_3 Point_3;


public:

	typedef enum FACE_STATE
	{
		FREE_FACE,
		FIXED_FACE
	};

public:

	FaceBase();

	void getHalfedges(Halfedge_handle h[3]);
	void getVertices(Vertex_handle v[3]);
	void getVertices(Vertex_const_handle v[3]) const;
	void getPoints(Point_3 p[3]) const;
	Vector_3 grad_u() const;
	Vector_3 grad_v() const;
	Vector_3& grad_u();
	Vector_3& grad_v();
	Vector_3 normal() const;
	typename FaceBase<Refs, Traits>::FACE_STATE state() const;
	typename FaceBase<Refs, Traits>::FACE_STATE& state();
	const std::complex<double>& prescribedMu() const;
	std::complex<double>& prescribedMu();
	const std::complex<double>& userComplex() const;
	std::complex<double>& userComplex();
	const double& user() const;
	double& user();
	const std::complex<double>& e1() const; //1st halfedge in local coordinate system
	const std::complex<double>& e2() const; //2nd halfedge in local coordinate system
	const std::complex<double>& e3() const; //3rd halfedge in local coordinate system
	double doubleArea() const;
	double area() const;
	double doubleUVArea() const;
	double uvArea() const;
	bool hasValidMetric() const;
	void updateHalfedgesInLocalCoords(bool useTargetMetric);
	bool computeHalfedgesInLocalCoords(std::complex<double>& e1, std::complex<double>& e2, std::complex<double>& e3, bool useTargetMetric) const;
	int& index();
	int index() const;
	int& localIndex();
	int localIndex() const;
	bool is_border_face() const; //checks if this face has either border edge or border vertex 
	bool is_border() const; //checks if this face has an edge on the border
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	double k1() const;
	double& k1();
	double k2() const;
	double& k2();
	Vector_3 kv1() const;
	Vector_3& kv1();
	Vector_3 kv2() const;
	Vector_3& kv2();
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private:

	//all the faces which have FREE_FACE state are numbered from 0 to numFreeFaces - 1
	//and all the faces which have FIXED_FACE state are numbered from 0 to numFixedFaces - 1
	int mLocalIndex;
	int mIndex;
	FACE_STATE mState;
	Vector_3 mGrad_u, mGrad_v;
	std::complex<double> mPrescribedMu;
	std::complex<double> m_e1, m_e2, m_e3;
	double mUser; //for debugging, checking and temporary computations
	std::complex<double> mUserComplex; //for debugging, checking and temporary computations
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	double _k1;
	double _k2;
	Vector_3 _kv1;
	Vector_3 _kv2;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
};




template <class Refs, typename Traits>
FaceBase<Refs, Traits>::FaceBase() :
mIndex(-1),
mGrad_u(0.0, 0.0, 0.0),
mGrad_v(0.0, 0.0, 0.0),
mPrescribedMu(0.0, 0.0),
m_e1(0.0, 0.0),
m_e2(0.0, 0.0),
m_e3(0.0, 0.0),
mUser(0.0),
mUserComplex(0.0, 0.0),
mState(FREE_FACE),
mLocalIndex(-1),
_k1(0),
_k2(0),
_kv1(0.0, 0.0, 0.0),
_kv2(0.0, 0.0, 0.0)
{
}



template<class Refs, typename Traits>
int& FaceBase<Refs, Traits>::index()
{
	return mIndex;
}



template<class Refs, typename Traits>
int FaceBase<Refs, Traits>::index() const
{
	return mIndex;
}

template<class Refs, typename Traits>
int& FaceBase<Refs, Traits>::localIndex()
{
	return mLocalIndex;
}


template<class Refs, typename Traits>
int FaceBase<Refs, Traits>::localIndex() const
{
	return mLocalIndex;
}


template<class Refs, typename Traits>
typename FaceBase<Refs, Traits>::FACE_STATE FaceBase<Refs, Traits>::state() const
{
	return mState;
}

template<class Refs, typename Traits>
typename FaceBase<Refs, Traits>::FACE_STATE& FaceBase<Refs, Traits>::state()
{
	return mState;
}


template <class Refs, typename Traits>
const std::complex<double>& FaceBase<Refs, Traits>::e1() const
{
	return m_e1;
}


template <class Refs, typename Traits>
const std::complex<double>& FaceBase<Refs, Traits>::e2() const
{
	return m_e2;
}


template <class Refs, typename Traits>
const std::complex<double>& FaceBase<Refs, Traits>::e3() const
{
	return m_e3;
}




template <class Refs, typename Traits>
std::complex<double>& FaceBase<Refs, Traits>::prescribedMu()
{
	return mPrescribedMu;
}



template <class Refs, typename Traits>
const std::complex<double>& FaceBase<Refs, Traits>::prescribedMu() const
{
	return mPrescribedMu;
}


template <class Refs, typename Traits>
std::complex<double>& FaceBase<Refs, Traits>::userComplex()
{
	return mUserComplex;
}



template <class Refs, typename Traits>
const std::complex<double>& FaceBase<Refs, Traits>::userComplex() const
{
	return mUserComplex;
}


template <class Refs, typename Traits>
double& FaceBase<Refs, Traits>::user()
{
	return mUser;
}



template <class Refs, typename Traits>
const double& FaceBase<Refs, Traits>::user() const
{
	return mUser;
}


template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3 FaceBase<Refs, Traits>::normal() const
{
	Point_3 p[3];
	getPoints(p);
	Vector_3 triangleNormal = cross_product(p[1] - p[0], p[2] - p[0]);

	double sqLength = triangleNormal.squared_length();

	if(sqLength != 0.0)
	{
		triangleNormal = triangleNormal/sqrt(sqLength);
	}
	return triangleNormal;
}


//returns the unsigned area. Assumes that the polygon is simple and planar
template <class Refs, typename Traits>
double FaceBase<Refs, Traits>::doubleArea() const
{
	Point_3 p[3];
	getPoints(p);
	Vector_3 u = p[0]-p[2];
	Vector_3 v = p[1]-p[2];

	Vector_3 normal = cross_product(u, v);

	Halfedge_const_handle h = halfedge();

	double doublePolygonArea = sqrt(normal.squared_length());

	if(!h->is_triangle()) //this face is not a triangle. here we assume that the polygon is simple and lies on a plane. otherwise the output area will be wrong
	{
		//otherwise the direction of the normal is not well defined and when we compute the dot product later on we can get wrong sign.
		//this can happen if we have a polygon where several vertices are collinear.
		//it is possible to solve this issue by searching for using another triplets of vertices of this polygon to compute the normal or by using SVD to fit the best plane
		//assert(doublePolygonArea > 1e-6);
		//if ( !(doublePolygonArea > 1e-6) )
			//Log.warn("FaceBase::doubleArea(), assert(doublePolygonArea > 1e-6) failed");

		Point_3 referecePoint = p[0]; //can be any point on the plane of the polygon

		Halfedge_around_facet_const_circulator hc = h->facet_begin();
		Halfedge_around_facet_const_circulator hc_end = hc;

		doublePolygonArea = 0.0;

		CGAL_For_all(hc, hc_end)
		{
			u = hc->vertex()->point() - referecePoint;
			v = hc->prev()->vertex()->point() - referecePoint;

			Vector_3 cross = cross_product(u, v);

			double unsignedDoubleArea = sqrt(cross.squared_length());

			double dot = cross*normal;

			if(dot > 0.0)
			{
				doublePolygonArea += unsignedDoubleArea;
			}
			else
			{
				doublePolygonArea -= unsignedDoubleArea;
			}
		}
		if(doublePolygonArea < 0.0)
		{
			doublePolygonArea = -doublePolygonArea;
		}
	}
	return doublePolygonArea;
}


template <class Refs, typename Traits>
double FaceBase<Refs, Traits>::area() const
{
	return 0.5*doubleArea();
}


template <class Refs, typename Traits>
double FaceBase<Refs, Traits>::doubleUVArea() const
{
	assert(halfedge()->is_triangle());

	Halfedge_const_handle h1 = halfedge();
	Halfedge_const_handle h2 = h1->next();
	Halfedge_const_handle h3 = h2->next();

	Point_3 uv1 = h1->vertex()->uv();
	Point_3 uv2 = h2->vertex()->uv();
	Point_3 uv3 = h3->vertex()->uv();

	Vector_3 cross = cross_product(uv1 - uv3, uv2 - uv3);
	return sqrt(cross.squared_length());
}

template <class Refs, typename Traits>
double FaceBase<Refs, Traits>::uvArea() const
{
	return 0.5*doubleUVArea();
}


template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3 FaceBase<Refs, Traits>::grad_u() const
{
	return mGrad_u;
}

template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3 FaceBase<Refs, Traits>::grad_v() const
{
	return mGrad_v;
}

template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3& FaceBase<Refs, Traits>::grad_u()
{
	return mGrad_u;
}

template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3& FaceBase<Refs, Traits>::grad_v()
{
	return mGrad_v;
}


template <class Refs, typename Traits>
bool FaceBase<Refs, Traits>::is_border() const
{
	Halfedge_const_handle h1 = halfedge();
	Halfedge_const_handle h2 = h1->next();
	Halfedge_const_handle h3 = h2->next();
	
	if(h1->opposite()->is_border() || h2->opposite()->is_border() || h3->opposite()->is_border())
	{
		return true;
	}
	else
	{
		return false;
	}
}


template <class Refs, typename Traits>
bool FaceBase<Refs, Traits>::is_border_face() const
{
	Halfedge_const_handle h1 = halfedge();
	Halfedge_const_handle h2 = h1->next();
	Halfedge_const_handle h3 = h2->next();

	if(h1->vertex()->is_border() || h2->vertex()->is_border() || h3->vertex()->is_border())
	{
		return true;
	}
	else
	{
		return false;
	}
}


template <class Refs, typename Traits>
void FaceBase<Refs, Traits>::getPoints(Point_3 p[3]) const
{
	Halfedge_const_handle h = halfedge();
	p[0] = h->opposite()->vertex()->point();
	p[1] = h->vertex()->point();
	p[2] = h->next()->vertex()->point();
}

template <class Refs, typename Traits>
void FaceBase<Refs, Traits>::getVertices(Vertex_handle v[3])
{
	Halfedge_handle h = halfedge();
	v[0] = h->opposite()->vertex();
	v[1] = h->vertex();
	v[2] = h->next()->vertex();
}

template <class Refs, typename Traits>
void FaceBase<Refs, Traits>::getVertices(Vertex_const_handle v[3]) const
{
	Halfedge_const_handle h = halfedge();
	v[0] = h->opposite()->vertex();
	v[1] = h->vertex();
	v[2] = h->next()->vertex();
}


template <class Refs, typename Traits>
void FaceBase<Refs, Traits>::getHalfedges(Halfedge_handle h[3])
{
	h[0] = halfedge();
	h[1] = h[0]->next();
	h[2] = h[1]->next();
}

//checks if the target metric satisfies the triangle inequality
template <class Refs, typename Traits>
bool FaceBase<Refs, Traits>::hasValidMetric() const
{
	double a = halfedge()->targetMetric();
	double b = halfedge()->opposite()->targetMetric();
	double c = halfedge()->next()->targetMetric();

	if(a>0.0 && b>0.0 && c>0.0 && a+b>c && a+c>b && b+c>a)
	{
		return true;
	}

	return false;
}


//for each triangle we define a *local* orthonormal coordinate system such that the first vector is a unit vector in the direction of the first halfedge.
//the origin is at p2 and the second vector is rotation by 90 of the first.
//e1, e2, e3 are not the coordinate system itself but rather the coordinates of the original vector edges in that coordinate system.
//when useTargetMetric is true, we use the prescribed edge lengths rather than the embedding of the face in 3D
//
//note that this function doesn't update the local variables m_e1, m_e2, m_e3.
//for computing and updating the variables use updateHalfedgesInLocalCoords instead
template <class Refs, typename Traits>
bool FaceBase<Refs, Traits>::computeHalfedgesInLocalCoords(std::complex<double>& e1, std::complex<double>& e2, std::complex<double>& e3, bool useTargetMetric) const
{
	Halfedge_const_handle h1 = halfedge();
	Halfedge_const_handle h2 = h1->next();
	Halfedge_const_handle h3 = h1->prev();

	if(useTargetMetric)
		//not tested yet!! in order to test compare with the result when useTargetMetric is false
	{
		double l1 = h1->targetMetric();
		double l2 = h2->targetMetric();
		double l3 = h3->targetMetric();

		if(l1 <= 0.0 || l2 <= 0.0 || l3 <= 0.0)
		{
			return false;
		}

		double cosine_13 = (l1*l1 + l3*l3 - l2*l2)/(2.0*l1*l3);
		double sine_13 = sqrt(1.0 - cosine_13*cosine_13);

		e1 = std::complex<double>(l1, 0.0);
		e3 = -l3*std::complex<double>(cosine_13, sine_13);
		e2 = -e1-e3;
	}
	else
	{
		Point_3 p2 = h3->vertex()->point();
		Point_3 p3 = h1->vertex()->point();
		Point_3 p1 = h2->vertex()->point();

		Vector_3 U = p3-p2;
		Vector_3 V = p1-p2;
		
		if(U == Vector_3(0.0, 0.0, 0.0) || V == Vector_3(0.0, 0.0, 0.0))
		{
			return false;
		}

		double u_length = sqrt(U.squared_length());
		double dot = U*V;
		double cross = sqrt(cross_product(U, V).squared_length());

		e1 = std::complex<double>(u_length, 0.0);
		e3 = -std::complex<double>(dot/u_length, cross/u_length);
		e2 = -e1-e3;
	}
	return true;
}

template <class Refs, typename Traits>
void FaceBase<Refs, Traits>::updateHalfedgesInLocalCoords(bool useTargetMetric)
{
	bool res = computeHalfedgesInLocalCoords(m_e1, m_e2, m_e3, useTargetMetric);
	assert(res);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template <class Refs, typename Traits>
double FaceBase<Refs, Traits>::k1() const
{
	return _k1;
}

template <class Refs, typename Traits>
double& FaceBase<Refs, Traits>::k1()
{
	return _k1;
}

template <class Refs, typename Traits>
double FaceBase<Refs, Traits>::k2() const
{
	return _k2;
}

template <class Refs, typename Traits>
double& FaceBase<Refs, Traits>::k2()
{
	return _k2;
}


template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3 FaceBase<Refs, Traits>::kv1() const
{
	return _kv1;
}

template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3& FaceBase<Refs, Traits>::kv1()
{
	return _kv1;
}


template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3 FaceBase<Refs, Traits>::kv2() const
{
	return _kv2;
}

template <class Refs, typename Traits>
typename Traits::Kernel::Vector_3& FaceBase<Refs, Traits>::kv2()
{
	return _kv2;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%