#pragma once 

#include <vector>
#include "Cone.h"
#include "Seam.h"
#include "Vec2.h"
#include "Vec3.h"


class MeshBuffer
{
public:
	MeshBuffer();
	~MeshBuffer();

	vector<int> idx_sizes, sample_sizes, sample_row_sizes;
	vector<int> idx_uv, idx_pos, idx_nor;
	vector<Vec3> positions, normals;
	vector<Vec2> uvs;
	vector<Seam> seams;
	vector<Cone> cones;
};