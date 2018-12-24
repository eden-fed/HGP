
#include "Vec3.h"

Vec3::Vec3(double x, double y, double z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

double& Vec3::operator[](int index)
{
	return v[index];
}

double Vec3::operator[](int index) const
{
	return v[index];
}
