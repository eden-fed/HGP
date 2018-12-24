
#include "Vec2.h"

Vec2::Vec2(double x, double y)
{
	v[0] = x;
	v[1] = y;
}

double& Vec2::operator[](int index)
{
	return v[index];
}

double Vec2::operator[](int index) const
{
	return v[index];
}
