#pragma once 
#include <iostream>
using namespace std;

class Vec2
{
public:
	double v[2];
	Vec2(){}
	Vec2(double x, double y);//{ v[0] = x; v[1] = y;}
	~Vec2(){}
	double& operator[](int index);//{ return v[index]; }
	double operator[](int index) const; //{ return v[index]; }
};