#pragma once

#include "CGAL\HalfedgeDS_halfedge_base.h"
#include "Utils\Limit_Util.h"

using namespace std;

template<class Refs, typename Traits>
class HalfedgeBase : public CGAL::HalfedgeDS_halfedge_base<Refs, CGAL::Tag_true, CGAL::Tag_true, CGAL::Tag_true>
{
	typedef typename Traits::Kernel::Vector_3 Vector_3;
	typedef typename Traits::Kernel::Point_3 Point_3;
	#define CONFORMAL_FACTOR_EPSILON (1e-8)
public:

	HalfedgeBase();
	Point_3& uv();
	const Point_3& uv() const;
	double& numeric1();
	double numeric1() const;
	double cot(bool useTargetMetric = false) const;
	double uvCot() const;
	double meanValueWeight(bool useTargetMetric = false) const;
	int& index();
	int index() const;
	int& userIndex();
	int userIndex() const;
	bool is_border_edge() const;
	double angle() const;
	double angleFromTargetMetric() const;
	double uvAngle(bool UVOnHalfEdge = false) const;
	double scaleFactor() const;
	double squared_length() const;
	double length() const;
	double uvLength() const;
	double targetMetric() const;
	double& targetMetric();
	bool is_separating_edge() const;
	double mixedCellArea(bool useTargetMetric = false) const;
	double uvMixedCellArea() const;
	double triangleAreaFromEdgeLengths(double a, double b, double c) const;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	bool isCut() const;
	bool& isCut();
	int rotation() const;
	int& rotation();
	int matching() const;
	int& matching();

	double get_dMetric() const;
	void set_dMetric(double val);

	double rotationAngle() const;
	double& rotationAngle();
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

private:

	Point_3 mUV;
	int mIndex;
	int mUserIndex;
	double mNumeric1; //uninitialized variable used as an auxiliary for external functions
	double mTargetMetric;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	bool _isCut;
	int rot;
	int _matching;

	double dMetric; //for Dijkstra's algoritem
	double _rotAng;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
};



template<class Refs, typename Traits>
HalfedgeBase<Refs, Traits>::HalfedgeBase() : mUV(0.0, 0.0, 0.0), mIndex(-1), mTargetMetric(0.0), mNumeric1(0.0), mUserIndex(-1), _isCut(false), rot(-1), _matching(0)
{

}

template<class Refs, typename Traits>
double& HalfedgeBase<Refs, Traits>::numeric1()
{
	return mNumeric1;
}

template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::numeric1() const
{
	return mNumeric1;
}

template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::scaleFactor() const
{
	double phi1 = vertex()->conformalFactor();
	double phi2 = opposite()->vertex()->conformalFactor();

	double scale = 0.0;

	if(std::abs(phi1-phi2) < CONFORMAL_FACTOR_EPSILON)
	{
		scale = exp(phi1);
	}
	else
	{
		scale = (exp(phi1) - exp(phi2)) / (phi1 - phi2);
	}
	return scale;
}

template<class Refs, typename Traits>
int& HalfedgeBase<Refs, Traits>::index()
{
	return mIndex;
}

template<class Refs, typename Traits>
int HalfedgeBase<Refs, Traits>::index() const
{
	return mIndex;
}

template<class Refs, typename Traits>
int& HalfedgeBase<Refs, Traits>::userIndex()
{
	return mUserIndex;
}

template<class Refs, typename Traits>
int HalfedgeBase<Refs, Traits>::userIndex() const
{
	return mUserIndex;
}

template<class Refs, typename Traits>
typename Traits::Kernel::Point_3& HalfedgeBase<Refs, Traits>::uv()
{
	return mUV;
}

template<class Refs, typename Traits>
const typename Traits::Kernel::Point_3& HalfedgeBase<Refs, Traits>::uv() const
{
	return mUV;
}


template<class Refs, typename Traits>
bool HalfedgeBase<Refs, Traits>::is_border_edge() const
{
	return this->is_border() || this->opposite()->is_border();
}

//this function computes the two parts of the weight
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::meanValueWeight(bool useTargetMetric = false) const
{
	assert(face()->is_triangle()); //making sure we don't have an n-gon or border edge
	assert(!is_border_edge());

	if(useTargetMetric) //compute the mean value weights based on the assigned edge lengths
	{
		double a = next()->mTargetMetric;
		double b = prev()->mTargetMetric;
		double c = mTargetMetric;
		double d = opposite()->prev()->mTargetMetric;
		double e = opposite()->next()->mTargetMetric;
		assert(a > 0.0 && b > 0.0 && c > 0.0 && d > 0.0 && e > 0.0);

		double alphaA_half = 0.0;
		double alphaB_half = 0.0;

		double nominatorA = (c-a-b)*(c+a-b);
		double denominatorA = (a+b+c)*(a-b-c);
		double nominatorB = (c-d-e)*(c+d-e);
		double denominatorB = (d+e+c)*(d-e-c);

		if((b+c > a) && (a+c > b) && (a+b > c)) //check triangle inequality
		{
			assert(denominatorA < 0.0);
			alphaA_half = sqrt(nominatorA / denominatorA);
		}
		else //triangle inequality doesn't hold
		{
			cerr << "Error at face: " << face()->index() << ". Triangle inequality doesn't hold" << endl;
			return 1.0; //this leads to less numerical issues
		}

		if((e+c > d) && (d+c > e) && (d+e > c)) //check triangle inequality
		{
			assert(denominatorB < 0.0);
			alphaB_half = sqrt(nominatorB / denominatorB);
		}
		else //triangle inequality doesn't hold
		{
			cerr << "Error at face: " << face()->index() << ". Triangle inequality doesn't hold" << endl;
			return 1.0; //this leads to less numerical issues
		}
		assert(c > 0.0); //if c is 0.0 we should have returned earlier since triangle inequality wouldn't hold
		
		double weight = (alphaA_half + alphaB_half)/c;
		assert(weight > 0.0); //mean value weights are supposed to be positive

		return weight;
	}
	else
	{
		Point_3 p0 = prev()->vertex()->point();
		Point_3 p1 = vertex()->point();
		Point_3 p2 = next()->vertex()->point();
		Point_3 p3 = opposite()->next()->vertex()->point();

		Vector_3 u = p1-p0;
		Vector_3 v = p2-p0;
		Vector_3 w = p3-p0;

		double u_length = sqrt(u.squared_length());
		double v_length = sqrt(v.squared_length());
		double w_length = sqrt(w.squared_length());

		double cross_uv = sqrt(cross_product(u, v).squared_length());
		double cross_wu = sqrt(cross_product(w, u).squared_length());

		double alphaUV_half = 0.0;
		double alphaWU_half = 0.0;

		if(cross_uv <= 0.0)
		{
			cerr << "Face: " << face()->index() << " is singular!" << endl;
			return 1.0; //this leads to less numerical issues
		}
		else
		{
			alphaUV_half = (u_length*v_length - u*v) / cross_uv;
		}

		if(cross_wu <= 0.0)
		{
			cerr << "Face: " << face()->index() << " is singular!" << endl;
			return 1.0; //this leads to less numerical issues
		}
		else
		{
			alphaWU_half = (w_length*u_length - w*u) / cross_wu;
		}

		double weight = (alphaUV_half + alphaWU_half)/u_length;
		assert(weight > 0.0); //mean value weights are supposed to be positive

		return weight;
	}
}


//this function computes only the part of the weight that belongs to the specific halfedge
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::cot(bool useTargetMetric) const
{
	if(is_border())
	{
		return 0.0;
	}

	if(useTargetMetric) //compute the cot weights based on the assigned edge lengths
	{
		double a = next()->mTargetMetric;
		double b = prev()->mTargetMetric;
		double c = mTargetMetric;
		assert(a > 0.0 && b > 0.0 && c > 0.0);

		double nominator = a*a + b*b - c*c;
		double denominatorSq = (b+c-a)*(a+b-c)*(a+c-b)*(a+b+c);

		if((b+c > a) && (a+c > b) && (a+b > c))
		{
			assert(denominatorSq > 0.0);
			double cot = nominator / sqrt(denominatorSq);
			return cot;
		}
		else //triangle inequality doesn't hold
		{
			cerr << "Error at face: " << face()->index() << ". Triangle inequality doesn't hold" << endl;
			
			//return 0.0; //that's what was recommended by Springborn Schroder
			return 1.0; //but this leads to less numerical issues
		}
	}
	else //use the mesh embedding
	{
		Point_3 p0 = opposite()->vertex()->point();
		Point_3 p1 = vertex()->point();
		Point_3 p2 = prev()->opposite()->vertex()->point();

		Vector_3 u = p0-p2;
		Vector_3 v = p1-p2;
		Vector_3 cross = cross_product(u, v);
		double doubleArea = sqrt(cross.squared_length());

		if(doubleArea <= 0.0)
		{
			cerr << "Face: " << face()->index() << " is singular!" << endl;
			return 1.0; //this leads to less numerical issues
		}
		double cotan = u*v / doubleArea;
		return cotan;
	}
}


//this function computes only the part of the weight that belongs to the specific halfedge.
//the computation of the cotangent is based on the uv domain which is treated here as a 3D mesh (embedded in 2D).
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::uvCot() const
{
	if(is_border())
	{
		return 0.0;
	}

	Point_3 p0 = prev()->vertex()->uv();
	Point_3 p1 = vertex()->uv();
	Point_3 p2 = next()->vertex()->uv();

	Vector_3 V20 = p0-p2;
	Vector_3 V21 = p1-p2;
	Vector_3 cross = cross_product(V20, V21);
	double doubleArea = sqrt(cross.squared_length());

	if(doubleArea <= 0.0)
	{
		cerr << "Face: " << face()->index() << " is singular!" << endl;
		return 1.0; //this leads to less numerical issues
	}

	double cotan = V20*V21 / doubleArea;
	return cotan;
}


template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::squared_length() const
{
	return (vertex()->point() - opposite()->vertex()->point()).squared_length();
}


template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::length() const
{
	return sqrt(squared_length());
}


template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::uvLength() const
{
	return sqrt((vertex()->uv() - opposite()->vertex()->uv()).squared_length());
}




//computes the angle that the halfedge points to
//returns a number in the range 0..PI
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::angle() const
{
	if(is_border())
	{
		return 0.0;
	}

	double u_sq = next()->squared_length();
	double v_sq = squared_length();
	double w_sq = prev()->squared_length();

	double cosine = (u_sq + v_sq - w_sq) / sqrt(4.0*u_sq*v_sq);
	double alpha = acos(cosine);

	assert(LIMIT::isFinite(alpha));
	assert(alpha > 0.0);

	return alpha;
}


//computes the angle that the halfedge points to, based on the uv layout
//returns a number in the range 0..PI
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::uvAngle(bool UVOnHalfEdge = false) const
{
	if(is_border())
	{
		return 0.0;
	}

	Point_3 p_w = vertex()->uv();
	Point_3 p_u = prev()->vertex()->uv();
	Point_3 p_v = next()->vertex()->uv();

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (UVOnHalfEdge)
	{
		p_w = uv();
		p_u = prev()->uv();
		p_v = next()->uv();
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	double u_sq = (p_v - p_w).squared_length();
	double v_sq = (p_w - p_u).squared_length();
	double w_sq = (p_u - p_v).squared_length();

	double cosine = (u_sq + v_sq - w_sq) / sqrt(4.0*u_sq*v_sq);
	double alpha = acos(cosine);

	assert(LIMIT::isFinite(alpha));
	assert(alpha > 0.0);
	
	return alpha;
}



template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::targetMetric() const
{
	return mTargetMetric;
}

template<class Refs, typename Traits>
double& HalfedgeBase<Refs, Traits>::targetMetric()
{
	return mTargetMetric;
}


//computes the angle that the halfedge points to by using the target metric (assigned edgeLength)
//returns a number in the range 0..PI
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::angleFromTargetMetric() const
{
	if(is_border())
	{
		return 0.0;
	}

	double a = next()->mTargetMetric;
	double b = mTargetMetric;
	double c = prev()->mTargetMetric;

	assert(a > 0.0);
	assert(b > 0.0);
	assert(c > 0.0);

	double cosine = (a*a + b*b - c*c) / (2.0*a*b);
	if(cosine >= 1.0)
	{
		return 0.0;
	}
	if(cosine <= -1.0)
	{
		return CGAL_PI;
	}
	double alpha = acos(cosine);
	assert(alpha == alpha); //this will be false if alpha is NAN (not a valid number)
	assert(alpha > 0.0);
	return alpha;
}

//a separating edge is an edge that has two vertices on the border but the edge itself is not a border edge.
//for meshes with disk topology, the existence of such an edge indicates that the underlying graph is not three-connected
template<class Refs, typename Traits>
bool HalfedgeBase<Refs, Traits>::is_separating_edge() const
{
	bool borderEdge = this->is_border_edge();

	if(borderEdge) return false;

	if(this->vertex()->is_border() && this->opposite()->vertex()->is_border())
	{
		return true;
	}
	return false;
}


//for triangle meshes an halfedge uniquely defines a vertex-triangle pair.
//this function computes the part of the mixed cell area (combination of Voronoi and barycentric area) that belongs to the vertex-triangle pair.
//based on the Caltech paper: "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds". See Fig 4 and Section 3.3.
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::mixedCellArea(bool useTargetMetric = false) const
{
	assert(face()->is_triangle()); //making sure we don't have an n-gon
	assert(!is_border()); //border halfedges belong to the outer face and are not associated with a vertex-triangle pair.

	double area = 0.0;

	if(useTargetMetric) //compute the area from the prescribed edge lengths
	{
		double a = mTargetMetric;
		double b = next()->mTargetMetric;
		double c = prev()->mTargetMetric;

		double a_2 = a*a;
		double b_2 = b*b;
		double c_2 = c*c;

		if(b_2 + c_2 <= a_2) //the triangle angle of the vertex-triangle pair is >= 90 degrees (easy to see from the cosine law)
		{
			area = 0.5*triangleAreaFromEdgeLengths(a, b, c);
		}
		else if((a_2 + c_2 <= b_2) || (a_2 + b_2 <= c_2)) //one of the other angles is >= 90 degrees
		{
			area = 0.25*triangleAreaFromEdgeLengths(a, b, c);
		}
		else //the triangle is acute and we can use the Voronoi area as the circumcenter of the triangle is located inside the triangle
		{
			area = 0.125*(b_2*next()->cot(true) + c_2*prev()->cot(true));
		}
	}
	else //compute the area from the 3D embedding
	{
		double a_2 = squared_length();
		double b_2 = next()->squared_length();
		double c_2 = prev()->squared_length();

		if(b_2 + c_2 <= a_2) //the triangle angle of the vertex-triangle pair is >= 90 degrees (easy to see from the cosine law)
		{
			area = 0.5*face()->area();
		}
		else if((a_2 + c_2 <= b_2) || (a_2 + b_2 <= c_2)) //one of the other angles is >= 90 degrees
		{
 			area = 0.25*face()->area();
		}
		else //the triangle is acute and we can use the Voronoi area as the circumcenter of the triangle is located inside the triangle
		{
			area = 0.125*(b_2*next()->cot(false) + c_2*prev()->cot(false));
		}
	}
	assert(area > 0.0);

	return area;
}


//for triangle meshes an halfedge uniquely defines a vertex-triangle pair.
//this function computes the part of the mixed cell area (combination of Voronoi and barycentric area) that belongs to the vertex-triangle pair.
//based on the Caltech paper: "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds". See Fig 4 and Section 3.3.
//all computations are based on the uv domain which is treated as a mesh. the 3d embedding is ignored.
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::uvMixedCellArea() const
{
	assert(face()->is_triangle()); //making sure we don't have an n-gon
	assert(!is_border()); //border halfedges belong to the outer face and are not associated with a vertex-triangle pair.

	double area = 0.0;


	Point_3 p0 = prev()->vertex()->uv();
	Point_3 p1 = vertex()->uv();
	Point_3 p2 = next()->vertex()->uv();

	double a_2 = (p1-p0).squared_length();
	double b_2 = (p2-p1).squared_length();
	double c_2 = (p0-p2).squared_length();

	if(b_2 + c_2 <= a_2) //the triangle angle of the vertex-triangle pair is >= 90 degrees (easy to see from the cosine law)
	{
		area = 0.5*face()->uvArea();
	}
	else if((a_2 + c_2 <= b_2) || (a_2 + b_2 <= c_2)) //one of the other angles is >= 90 degrees
	{
		area = 0.25*face()->uvArea();
	}
	else //the triangle is acute and we can use the Voronoi area as the circumcenter of the triangle is located inside the triangle
	{
		area = 0.125*(b_2*next()->uvCot() + c_2*prev()->uvCot());
	}

	assert(area > 0.0);

	return area;
}


//a utility function to compute triangle area from edge lengths based on Heron's formula.
template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::triangleAreaFromEdgeLengths(double a, double b, double c) const
{
	double s = 0.5*(a + b + c); //s is half the triangle's perimeter

	double triangleArea = sqrt(s*(s-a)*(s-b)*(s-c));
	
	return triangleArea;
}

//-------------------Alon addition--------------------

template<class Refs, typename Traits>
bool HalfedgeBase<Refs, Traits>::isCut() const
{
	return (_isCut);
}

template<class Refs, typename Traits>
bool& HalfedgeBase<Refs, Traits>::isCut()
{
	return (_isCut);
}


template<class Refs, typename Traits>
int HalfedgeBase<Refs, Traits>::rotation() const
{
	return (rot);
}

template<class Refs, typename Traits>
int& HalfedgeBase<Refs, Traits>::rotation()
{
	return (rot);
}

template<class Refs, typename Traits>
int HalfedgeBase<Refs, Traits>::matching() const
{
	return (_matching);
}

template<class Refs, typename Traits>
int& HalfedgeBase<Refs, Traits>::matching()
{
	return (_matching);
}


template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::get_dMetric() const
{
	return (dMetric);
}

template<class Refs, typename Traits>
void HalfedgeBase<Refs, Traits>::set_dMetric(double val)
{
	dMetric = val;
}


template<class Refs, typename Traits>
double HalfedgeBase<Refs, Traits>::rotationAngle() const
{
	return (_rotAng);
}

template<class Refs, typename Traits>
double& HalfedgeBase<Refs, Traits>::rotationAngle()
{
	return (_rotAng);
}
//-----------------------------------------------