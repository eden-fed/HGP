#pragma once
#include <complex>
#include "CGAL\HalfedgeDS_vertex_base.h"
//#include "CGAL_Macros.h"

template<class Refs, typename Traits>
class Vertex_Base : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, typename Traits::Kernel::Point_3>
{
	typedef typename Traits::Kernel::Vector_3 Vector_3;
	typedef typename Traits::Kernel::Point_3 Point_3;

public:

	#define GAUSS_CURVATURE_ZERO_TOLERANCE (1e-10)

	Vertex_Base();
	Vertex_Base(const Point_3& p);
	Vertex_Base(const Vertex_Base& v);

	int& index();
	int index() const;
	int& localIndex();
	int localIndex() const;
	int& xLocalIndex();
	int xLocalIndex() const;
	int& yLocalIndex();
	int yLocalIndex() const;
	int& userIndex();
	int userIndex() const;
	int& fullRotationIndex();
	int fullRotationIndex() const;
	bool is_border() const;
	bool isFree() const;
	bool isXfree() const;
	bool isYfree() const;
	bool& isFree();
	bool& isXfree();
	bool& isYfree();
	Point_3& uv();
	const Point_3& uv() const;
	std::complex<double> uvC() const;
	std::complex<double> uvC(const std::complex<double>& uv);
	Vector_3& stretchDirOnUnitSphere();
	const Vector_3& stretchDirOnUnitSphere() const;
	std::complex<double>& stretchDirOnExtendedPlane();
	const std::complex<double>& stretchDirOnExtendedPlane() const;
	std::complex<double>& stretchDirOnUnitDisk();
	const std::complex<double>& stretchDirOnUnitDisk() const;
	double coneAngle() const;
	double gaussianCurvature(bool mixedCellAreaNormalization = false) const;
	double uvGeodesicCurvature() const;
	double phiMetricGaussianCurvature() const;
	double& phiMetricGaussianCurvature();
	double targetMetricConeAngle() const;
	double targetMetricGaussianCurvature(bool mixedCellAreaNormalization = false) const;
	double prescribedConeAngle() const;
	double prescribedGaussianCurvature() const;
	double& prescribedGaussianCurvature();
	double uvConeAngle(bool UVOnHalfEdge = false) const;
	double uvGaussianCurvature(bool mixedCellAreaNormalization = false) const;
	double conformalFactor() const;
	double& conformalFactor();
	std::complex<double> quadraticDifferential() const;
	std::complex<double>& quadraticDifferential();
	std::complex<double> projectedPoint() const;
	std::complex<double>& projectedPoint();
	double mixedCellArea(bool useTargetMetric = false) const;
	double uvMixedCellArea() const;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	bool& isCone();
	bool isCone() const;
	bool& onCut();
	bool onCut() const;
	double getConeAngle() const;
	void setConeAngle(double ang);

	// for Dijkstra:
	int& dPrevIndex();
	int dPrevIndex() const;
	double& distanceFromSource();
	double distanceFromSource() const;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


private:

	int mIndex;

	//all the free vertices are numbered from 0 to numFreeVertices - 1
	//and all the fixed vertices are numbered from 0 to numFixedVertices - 1
	int mLocalIndex;
	int mXLocalIndex;
	int mYLocalIndex;
	int mUserIndex;
	int mFullRotationIndex;
	bool mIsFree; //is the whole vertex free
	bool mIsXfree; //is the x component of the vertex free
	bool mIsYfree; //is the y component of the vertex free
	Point_3 mUV;
	Vector_3 mStretchDirOnUnitSphere;
	std::complex<double> mStretchDirOnExtendedPlane;
	std::complex<double> mStretchDirOnUnitDisk;
	double mConformalFactor;
	double mPrescribedGaussianCurvature;
	double mPhiMetricGaussianCurvature;
	std::complex<double> mQuadraticDifferential;
	std::complex<double> mProjectedPoint;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	bool _isCone;
	bool _onCut;
	double _coneAngle;

	int dPrev;	//previous halfedge index in Dijkstra algorithem
	double _distanceFromSource;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
};







template<class Refs, typename Traits>
Vertex_Base<Refs, Traits>::Vertex_Base() :
mIndex(-1),
mLocalIndex(-1),
mXLocalIndex(-1),
mYLocalIndex(-1),
mUserIndex(-1),
mFullRotationIndex(0),
mIsFree(true),
mIsXfree(true),
mIsYfree(true),
mConformalFactor(0.0),
mPrescribedGaussianCurvature(0.0),
mPhiMetricGaussianCurvature(0.0),
mUV(0.0, 0.0, 0.0),
mQuadraticDifferential(0.0, 0.0),
mProjectedPoint(0.0, 0.0),
mStretchDirOnUnitSphere(0.0, 0.0, 0.0),
mStretchDirOnExtendedPlane(0.0, 0.0),
mStretchDirOnUnitDisk(0.0, 0.0),
_isCone(false),
_onCut(false),
_coneAngle(0.0),
dPrev(-1)
{
}


template<class Refs, typename Traits>
Vertex_Base<Refs, Traits>::Vertex_Base(const Vertex_Base& v) :
mIndex(-1),
mLocalIndex(-1),
mXLocalIndex(-1),
mYLocalIndex(-1),
mUserIndex(-1),
mFullRotationIndex(0),
mIsFree(true),
mIsXfree(true),
mIsYfree(true),
mConformalFactor(0.0),
mPrescribedGaussianCurvature(0.0),
mPhiMetricGaussianCurvature(0.0),
mUV(0.0, 0.0, 0.0),
mQuadraticDifferential(0.0, 0.0),
mProjectedPoint(0.0, 0.0),
mStretchDirOnUnitSphere(0.0, 0.0, 0.0),
mStretchDirOnExtendedPlane(0.0, 0.0),
mStretchDirOnUnitDisk(0.0, 0.0),
_isCone(false),
_onCut(false),
_coneAngle(0.0),
dPrev(-1),
CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point_3>(v)
{
}


template<class Refs, typename Traits>
Vertex_Base<Refs, Traits>::Vertex_Base(const Point_3& p) :
mIndex(-1),
mLocalIndex(-1),
mXLocalIndex(-1),
mYLocalIndex(-1),
mUserIndex(-1),
mFullRotationIndex(0),
mIsFree(true),
mIsXfree(true),
mIsYfree(true),
mConformalFactor(0.0),
mPrescribedGaussianCurvature(0.0),
mPhiMetricGaussianCurvature(0.0),
mUV(0.0, 0.0, 0.0),
mQuadraticDifferential(0.0, 0.0),
mProjectedPoint(0.0, 0.0),
mStretchDirOnUnitSphere(0.0, 0.0, 0.0),
mStretchDirOnExtendedPlane(0.0, 0.0),
mStretchDirOnUnitDisk(0.0, 0.0),
_isCone(false),
_onCut(false),
_coneAngle(0.0),
dPrev(-1),
CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point_3>(p)
{
}



template<class Refs, typename Traits>
std::complex<double>& Vertex_Base<Refs, Traits>::quadraticDifferential()
{
	return mQuadraticDifferential;
}

template<class Refs, typename Traits>
std::complex<double> Vertex_Base<Refs, Traits>::quadraticDifferential() const
{
	return mQuadraticDifferential;
}


template<class Refs, typename Traits>
std::complex<double>& Vertex_Base<Refs, Traits>::projectedPoint()
{
	return mProjectedPoint;
}

template<class Refs, typename Traits>
std::complex<double> Vertex_Base<Refs, Traits>::projectedPoint() const
{
	return mProjectedPoint;
}


template<class Refs, typename Traits>
double& Vertex_Base<Refs, Traits>::conformalFactor()
{
	return mConformalFactor;
}

template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::conformalFactor() const
{
	return mConformalFactor;
}


template<class Refs, typename Traits>
typename Traits::Kernel::Vector_3& Vertex_Base<Refs, Traits>::stretchDirOnUnitSphere()
{
	return mStretchDirOnUnitSphere;
}

template<class Refs, typename Traits>
const typename Traits::Kernel::Vector_3& Vertex_Base<Refs, Traits>::stretchDirOnUnitSphere() const
{
	return mStretchDirOnUnitSphere;
}


template<class Refs, typename Traits>
std::complex<double>& Vertex_Base<Refs, Traits>::stretchDirOnExtendedPlane()
{
	return mStretchDirOnExtendedPlane;
}

template<class Refs, typename Traits>
const std::complex<double>& Vertex_Base<Refs, Traits>::stretchDirOnExtendedPlane() const
{
	return mStretchDirOnExtendedPlane;
}


template<class Refs, typename Traits>
std::complex<double>& Vertex_Base<Refs, Traits>::stretchDirOnUnitDisk()
{
	return mStretchDirOnUnitDisk;
}

template<class Refs, typename Traits>
const std::complex<double>& Vertex_Base<Refs, Traits>::stretchDirOnUnitDisk() const
{
	return mStretchDirOnUnitDisk;
}



template<class Refs, typename Traits>
typename Traits::Kernel::Point_3& Vertex_Base<Refs, Traits>::uv()
{
	return mUV;
}

template<class Refs, typename Traits>
const typename Traits::Kernel::Point_3& Vertex_Base<Refs, Traits>::uv() const
{
	return mUV;
}


template<class Refs, typename Traits>
std::complex<double> Vertex_Base<Refs, Traits>::uvC(const std::complex<double>& uv)
{
	double u = mUV.x();
	double v = mUV.y();

	mUV = Point_3(uv.real(), uv.imag(), 0.0);

	return std::complex<double>(u, v);
}

template<class Refs, typename Traits>
std::complex<double> Vertex_Base<Refs, Traits>::uvC() const
{
	return std::complex<double>(mUV.x(), mUV.y());
}



template<class Refs, typename Traits>
bool Vertex_Base<Refs, Traits>::isFree() const
{
	return mIsFree;
}

template<class Refs, typename Traits>
bool& Vertex_Base<Refs, Traits>::isFree()
{
	return mIsFree;
}

template<class Refs, typename Traits>
bool Vertex_Base<Refs, Traits>::isXfree() const
{
	return mIsXfree;
}

template<class Refs, typename Traits>
bool& Vertex_Base<Refs, Traits>::isXfree()
{
	return mIsXfree;
}

template<class Refs, typename Traits>
bool Vertex_Base<Refs, Traits>::isYfree() const
{
	return mIsYfree;
}

template<class Refs, typename Traits>
bool& Vertex_Base<Refs, Traits>::isYfree()
{
	return mIsYfree;
}

template<class Refs, typename Traits>
int& Vertex_Base<Refs, Traits>::index()
{
	return mIndex;
}

template<class Refs, typename Traits>
int Vertex_Base<Refs, Traits>::index() const
{
	return mIndex;
}

template<class Refs, typename Traits>
int& Vertex_Base<Refs, Traits>::localIndex()
{
	return mLocalIndex;
}


template<class Refs, typename Traits>
int Vertex_Base<Refs, Traits>::localIndex() const
{
	return mLocalIndex;
}

template<class Refs, typename Traits>
int& Vertex_Base<Refs, Traits>::xLocalIndex()
{
	return mXLocalIndex;
}

template<class Refs, typename Traits>
int Vertex_Base<Refs, Traits>::xLocalIndex() const
{
	return mXLocalIndex;
}

template<class Refs, typename Traits>
int& Vertex_Base<Refs, Traits>::yLocalIndex()
{
	return mYLocalIndex;
}

template<class Refs, typename Traits>
int Vertex_Base<Refs, Traits>::yLocalIndex() const
{
	return mYLocalIndex;
}

template<class Refs, typename Traits>
int& Vertex_Base<Refs, Traits>::userIndex()
{
	return mUserIndex;
}

template<class Refs, typename Traits>
int Vertex_Base<Refs, Traits>::userIndex() const
{
	return mUserIndex;
}

template<class Refs, typename Traits>
int& Vertex_Base<Refs, Traits>::fullRotationIndex()
{
	return mFullRotationIndex;
}

template<class Refs, typename Traits>
int Vertex_Base<Refs, Traits>::fullRotationIndex() const
{
	return mFullRotationIndex;
}


//compute the discrete Gaussian curvature at the vertex.
//if mixedCellAreaNormalization is false the angle deficit is being computed (which is an integrated quantity).
//otherwise, the angle deficit is further divided by the area of the mixed cell element to provide a pointwise approximation of the curvature.
//computations are done based on the 3D embedding of the mesh.
template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::gaussianCurvature(bool mixedCellAreaNormalization = false) const
{
	double k = 0.0;

	if(is_border())
	{
		k = CGAL_PI-coneAngle();
	}
	else //internal vertex
	{
		k = 2.0*CGAL_PI-coneAngle();
	}
	if(mixedCellAreaNormalization)
	{
		k = k / mixedCellArea(false);
	}

	if(std::abs(k) < GAUSS_CURVATURE_ZERO_TOLERANCE) k = 0.0; //just to get nicer values for the common case where the metric is flat
	return k;
}


//compute the discrete Gaussian curvature at the vertex.
//if mixedCellAreaNormalization is false the angle deficit is being computed (which is an integrated quantity).
//otherwise, the angle deficit is further divided by the area of the mixed cell element to provide a pointwise approximation of the curvature.
//computations are done based on the uv domain which is treated as a 3D mesh (embedded in 2D).
template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::uvGaussianCurvature(bool mixedCellAreaNormalization = false) const
{
	double k = 0.0;

	if(is_border())
	{
		k = CGAL_PI-uvConeAngle();
	}
	else //internal vertex
	{
		k = 2.0*CGAL_PI-uvConeAngle();
	}
	if(mixedCellAreaNormalization)
	{
		k = k / uvMixedCellArea();
	}

	if(std::abs(k) < GAUSS_CURVATURE_ZERO_TOLERANCE) k = 0.0; //just to get nicer values for the common case where the metric is flat
	return k;
}


//computes the geodesic curvature of a boundary curve at a boundary vertex.
//use the uv layout rather than the 3D mesh.
//note that this is different from uvGaussianCurvature that uses the sum of angles taken from adjacent triangles
template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::uvGeodesicCurvature() const
{
	if(!is_border())
	{
		assert(0);
		return 0.0;
	}

	Halfedge_const_handle h = halfedge();
	const Halfedge_const_handle h_end = h;
	bool foundIt = false;
	do
	{
		if(h->is_border())
		{
			foundIt = true;
			break;
		}

		h = h->next_on_vertex();
	} while(h != h_end);

	assert(foundIt);

	std::complex<double> pa = h->next()->vertex()->uvC();
	std::complex<double> pb = uvC();
	std::complex<double> pc = h->opposite()->vertex()->uvC();

	double k = 0.0;

	//we use the orient test in order to get a robust decision about the sign of the curvature
	double cross = CounterClockWise(pa, pb, pc);
	if(cross == 0.0)
	{
		k = 0.0;
	}
	else
	{
		std::complex<double> U = pc - pb;
		std::complex<double> V = pa - pb;

		double dot = U.real()*V.real() + U.imag()*V.imag();
		double lengthSQ = abs(U*conj(U)*V*conj(V));
		assert(lengthSQ > 0.0);

		double num = dot / sqrt(lengthSQ);

		if(num <= -1.0)
		{
			k = 0.0;
		}
		else if(num >= 1.0) //the angle between U and V is zero
		{
			if(cross > 0.0)
			{
				k = CGAL_PI;
			}
			else //cross < 0.0
			{
				k = -CGAL_PI;
			}
		}
		else
		{
			double angle = acos(num);

			if(cross > 0.0)
			{
				k = CGAL_PI - abs(angle);
			}
			else //cross < 0.0
			{
				k = abs(angle) - CGAL_PI;
			}
		}
	}

	return k;
}



template<class Refs, typename Traits>
bool Vertex_Base<Refs, Traits>::is_border() const
{
	Halfedge_const_handle h = halfedge();
	if(h == NULL) return true;
	const Halfedge_const_handle h_end = h;
	do
	{
		if(h->is_border())
		{
			return true;
		}
		h = h->next_on_vertex();
	} while (h != h_end);
	return false;
}


template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::coneAngle() const
{
	double anglesSum = 0.0;

	Halfedge_const_handle h = halfedge();
	const Halfedge_const_handle h_end = h;

	do
	{
		if(!h->is_border())
		{
			anglesSum += h->angle();
		}
		h = h->next_on_vertex();
	} while(h != h_end);

	return anglesSum;
}



template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::prescribedConeAngle() const
{
	if(is_border())
	{
		return CGAL_PI-mPrescribedGaussianCurvature;
	}
	else //internal vertex
	{
		return 2.0*CGAL_PI-mPrescribedGaussianCurvature;
	}
}



template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::uvConeAngle(bool UVOnHalfEdge = false) const
{
	double anglesSum = 0.0;

	Halfedge_const_handle h = halfedge();
	const Halfedge_const_handle h_end = h;

	do
	{
		if(!h->is_border())
		{
			anglesSum += h->uvAngle(UVOnHalfEdge);
		}
		h = h->next_on_vertex();
	} while(h != h_end);

	return anglesSum;
}



template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::prescribedGaussianCurvature() const
{
	return mPrescribedGaussianCurvature;
}


template<class Refs, typename Traits>
double& Vertex_Base<Refs, Traits>::prescribedGaussianCurvature()
{
	return mPrescribedGaussianCurvature;
}



template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::phiMetricGaussianCurvature() const
{
	return mPhiMetricGaussianCurvature;
}


template<class Refs, typename Traits>
double& Vertex_Base<Refs, Traits>::phiMetricGaussianCurvature()
{
	return mPhiMetricGaussianCurvature;
}



//compute the cone angle of the abstract (i.e., we don't have embedding for it) manifold based on the target metric
template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::targetMetricConeAngle() const
{
	double cone = 0.0;

	Halfedge_const_handle h = halfedge();
	const Halfedge_const_handle h_end = h;

	do{
		if(!h->is_border())
		{
			cone += h->angleFromTargetMetric();
		}
		h = h->next_on_vertex();
	} while(h != h_end);

	return cone;
}



//compute the discrete Gaussian curvature at the vertex.
//if mixedCellAreaNormalization is false the angle deficit is being computed (which is an integrated quantity).
//otherwise, the angle deficit is further divided by the area of the mixed cell element to provide a pointwise approximation of the curvature.
//computations are done based on the prescribed target metric (edge lengths).
template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::targetMetricGaussianCurvature(bool mixedCellAreaNormalization = false) const
{
	double k = 0.0;

	if(is_border())
	{
		k = CGAL_PI-targetMetricConeAngle();
	}
	else
	{
		k = 2.0*CGAL_PI-targetMetricConeAngle();
	}
	if(mixedCellAreaNormalization)
	{
		k = k / mixedCellArea(true);
	}

	if(std::abs(k) < GAUSS_CURVATURE_ZERO_TOLERANCE) k = 0.0; //just to get nicer values for the common case where the metric is flat

	return k;
}

//this function computes the mixed cell area (combination of Voronoi and barycentric area) for the given vertex.
//area computation are either based on 3D embedding or the target metric.
template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::mixedCellArea(bool useTargetMetric = false) const
{
	double area = 0.0;

	Halfedge_const_handle h = halfedge();
	const Halfedge_const_handle h_end = h;

	do
	{
		if(!h->is_border())
		{
			area += h->prev()->mixedCellArea(useTargetMetric);
		}
		h = h->next_on_vertex();
	} while(h != h_end);

	return area;
}


//this function computes the mixed cell area (combination of Voronoi and barycentric area) for the given vertex.
//area computation is based on the uv domain (which is treated as a mesh).
template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::uvMixedCellArea() const
{
	double area = 0.0;

	Halfedge_const_handle h = halfedge();
	const Halfedge_const_handle h_end = h;

	do
	{
		if(!h->is_border())
		{
			area += h->prev()->uvMixedCellArea();
		}
		h = h->next_on_vertex();
	} while(h != h_end);

	return area;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template<class Refs, typename Traits>
bool& Vertex_Base<Refs, Traits>::isCone()
{
	return _isCone;
}

template<class Refs, typename Traits>
bool Vertex_Base<Refs, Traits>::isCone() const
{
	return _isCone;
}

template<class Refs, typename Traits>
bool& Vertex_Base<Refs, Traits>::onCut()
{
	return _onCut;
}

template<class Refs, typename Traits>
bool Vertex_Base<Refs, Traits>::onCut() const
{
	return _onCut;
}

template<class Refs, typename Traits>
void Vertex_Base<Refs, Traits>::setConeAngle(double ang)
{
	_coneAngle = ang;
}

template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::getConeAngle() const
{
	return _coneAngle;
}


template<class Refs, typename Traits>
int& Vertex_Base<Refs, Traits>::dPrevIndex()
{
	return dPrev;
}

template<class Refs, typename Traits>
int Vertex_Base<Refs, Traits>::dPrevIndex() const
{
	return dPrev;
}

template<class Refs, typename Traits>
double& Vertex_Base<Refs, Traits>::distanceFromSource()
{
	return _distanceFromSource;
}

template<class Refs, typename Traits>
double Vertex_Base<Refs, Traits>::distanceFromSource() const
{
	return _distanceFromSource;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%