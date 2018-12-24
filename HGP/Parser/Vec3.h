#pragma once 
#include <iostream>
using namespace std;

class Vec3
{
public:
	double v[3];
	Vec3(){}
	Vec3(double x, double y, double z);
	~Vec3(){}
	double& operator[](int index);//{ return v[index]; }
	double operator[](int index) const;// { return v[index]; }
};