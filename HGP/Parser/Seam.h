#pragma once 
#include <iostream>
using namespace std;

class Seam
{
public:
	int faceidx, posidx, rot;
	Seam(){}
	Seam(int f, int p, int r);//{ faceidx = f; posidx = p; rot = r; }
	~Seam(){}
};