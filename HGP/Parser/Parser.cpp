
#include "Parser.h"
#include <iostream>
using namespace std;
#include <fstream>

bool Parser::loadOBJ(const char *filename, MeshBuffer *mb, MeshBuffer *smb)
{
	bool has_normals = false, has_uvs = false;

	double x, y, z, l;
	//Vec8 cv;
	//Vec9 tn;
	unsigned int posidx, noridx, uvidx, faceidx;
	unsigned int n;
	int rot;

	string line, token;

	ifstream ifs(filename);
	if (!ifs.is_open())
	{
		cerr << "Failed to open file: " << filename << endl;
		return false;
	}
	istringstream iss;

	while (getline(ifs, line))
	{
		iss.clear();
		iss.str(line);

		iss >> token;

		if (token.compare("f") == 0)
		{
			mb->idx_sizes.push_back(0);
			if (smb) {
				smb->idx_sizes.push_back(0);
				smb->sample_sizes.push_back(0);
				smb->sample_row_sizes.push_back(0);
			}

			while (!iss.eof())
			{
				bool valid_entry = true;

				iss >> token; if (iss.fail()) break;
				if (sscanf(token.c_str(), "%d/%d/%d", &posidx, &uvidx, &noridx) == 3)
				{
					if (!has_uvs) mb->idx_uv.reserve(mb->idx_pos.capacity());
					if (!has_normals) mb->idx_nor.reserve(mb->idx_pos.capacity());

					has_uvs = has_normals = true;
				}
				else if (sscanf(token.c_str(), "%d//%d", &posidx, &noridx) == 2)
				{
					if (!has_normals) mb->idx_nor.reserve(mb->idx_pos.capacity());

					has_normals = true;
				}
				else if (sscanf(token.c_str(), "%d/%d", &posidx, &uvidx) == 2)
				{
					if (!has_uvs) mb->idx_uv.reserve(mb->idx_pos.capacity());

					has_uvs = true;
				}
				else if (sscanf(token.c_str(), "%d/", &posidx) == 1)
				{
					// nothing to change.
				}
				else
				{
					valid_entry = false;

				}
				if (valid_entry)
				{
					posidx--; noridx--; uvidx--;
					mb->idx_sizes.back()++;
					mb->idx_pos.push_back(posidx);
					if (has_normals) mb->idx_nor.push_back(noridx);
					if (has_uvs)     mb->idx_uv.push_back(uvidx);
				}
			}
		}
		else if (token.compare("v") == 0)
		{
			iss >> x >> y >> z;
			if (!iss.fail())
			{
				mb->positions.push_back(Vec3(x, y, z));
			}
		}
		else if (token.compare("vn") == 0)
		{
			iss >> x >> y >> z;
			l = sqrt(x*x + y*y + z*z); if (l == 0) l = 1;
			if (!iss.fail())
			{
				mb->normals.push_back(Vec3(x / l, y / l, z / l));
			}
		}
		else if (token.compare("vt") == 0)
		{
			iss >> x >> y;
			if (!iss.fail())
			{
				mb->uvs.push_back(Vec2(x, y));
			}
		}
		else if (token.compare("sm") == 0)
		{
			iss >> faceidx >> posidx; faceidx--; posidx--;
			iss >> rot;
			mb->seams.push_back(Seam(faceidx, posidx, rot));
		}
		else if (token.compare("c") == 0)
		{
			iss >> posidx; posidx--;
			iss >> x;
			if (!iss.fail())
			{
				mb->cones.push_back(Cone(posidx, x));
			}
		}
		else if (token.compare("#num_v") == 0)
		{
			iss >> n; if (!iss.fail()) mb->positions.reserve(n);
		}
		else if (token.compare("#num_vt") == 0)
		{
			iss >> n; if (!iss.fail()) mb->uvs.reserve(n);
		}
		else if (token.compare("#num_vn") == 0)
		{
			iss >> n; if (!iss.fail()) mb->normals.reserve(n);
		}
		else if ((smb) && (token.compare("#num_$b") == 0))
		{
			iss >> n; if (!iss.fail()) smb->idx_pos.reserve(n);
		}
		else if ((smb) && (token.compare("#num_$v") == 0))
		{
			iss >> n;
			if (!iss.fail()) smb->positions.reserve(n);
		}
		else if ((smb) && (token.compare("#num_$vt") == 0))
		{
			iss >> n; if (!iss.fail()) smb->uvs.reserve(n);
		}
		else if ((smb) && (token.compare("#num_$vn") == 0))
		{
			iss >> n; if (!iss.fail()) smb->normals.reserve(n);
		}
		else if (token.compare("#num_f") == 0)
		{
			iss >> n;
			if (!iss.fail())
			{
				mb->idx_sizes.reserve(n);
				if (smb)
				{
					smb->sample_sizes.reserve(n);
					smb->sample_row_sizes.reserve(n);
				}
			}
		}
		else if (token.compare("#num_f_v") == 0)
		{
			iss >> n; if (!iss.fail()) mb->idx_pos.reserve(n);
		}
	
	}

	ifs.close();

	if (!has_normals) mb->idx_nor.clear();
	if (!has_uvs) mb->idx_uv.clear();

	cout << "Loaded " << filename << endl;

	return true;
}

bool Parser::loadVectorField(const char *filename, Mesh& cgalMesh)
{
	std::ifstream fin(filename);

	if (!fin)	// problem reading the file
		return false;

	std::string line, token;
	istringstream iss, temp;

	while (getline(fin, line))	// search for k1
	{
		iss.clear();
		iss.str(line);
		iss >> token;
		if (token.compare("k1") == 0)
			break;
	}

	double k1, k2, kv1x, kv1y, kv1z, kv2x, kv2y, kv2z;
	int faceIndex = 0;

	while (getline(fin, line))	// scan data of k1,k2,k1V,k2V
	{
		iss.clear();
		iss.str(line);
		temp.clear();
		temp.str(line);
		temp >> token;
		if (token.compare("crossfield_angles") == 0)
			break;

		iss >> k1 >> k2 >> kv1x >> kv1y >> kv1z >> kv2x >> kv2y >> kv2z;

		if (cgalMesh.size_of_facets() == faceIndex)
			return false;	//error

		auto f = cgalMesh.face(faceIndex);
		f->k1() = k1;
		f->k2() = k2;
		f->kv1() = Vector_3(kv1x, kv1y, kv1z);
		f->kv2() = Vector_3(kv2x, kv2y, kv2z);
		faceIndex++;
	}

	bool foundFlag = false;
	while (getline(fin, line))	// search for MI_matchings_and_sharp
	{
		temp.clear();
		temp.str(line);
		temp >> token;
		if (token.compare("MI_matchings_and_sharp") == 0)
		{
			foundFlag = true;
			break;
		}
	}

	if (!foundFlag)
		return false;	//error

	int v1, v2, v3, h1, h2, h3;
	faceIndex = 0;
	while (getline(fin, line))	// scan data of matching
	{
		iss.clear();
		iss.str(line);
		iss >> v1 >> v2 >> v3 >> h1 >> h2 >> h3;

		auto f = cgalMesh.face(faceIndex);
		auto h = f->halfedge();

#ifdef DOUBLE_CHECK
		Vertex_handle checkV[3];
		f->getVertices(checkV);
		bool check0 = checkV[0]->index() == v1 || checkV[0]->index() == v2 || checkV[0]->index() == v3;
		bool check1 = checkV[1]->index() == v1 || checkV[1]->index() == v2 || checkV[1]->index() == v3;
		bool check2 = checkV[2]->index() == v1 || checkV[2]->index() == v2 || checkV[2]->index() == v3;
		if (!(check0 && check1 && check2))
			return false;	//problem!!
#endif

		while (h->vertex()->index() != v1)
		{
			h = h->next();
			if (h == f->halfedge())
				return false;	//error
		}
		h->matching() = h1;
		h = h->next();
		h->matching() = h2;
		h = h->next();
		h->matching() = h3;

		faceIndex++;
	}


	return true;
}

void Parser::setMeshAdditionalData(Mesh& cgalMesh, MeshBuffer& m1)
{
	for (int i = 0; i < (int)m1.cones.size(); ++i)
	{
		Vertex_handle v = cgalMesh.vertex(m1.cones[i].posidx);
		v->isCone() = true;
		v->setConeAngle(m1.cones[i].x);
		v->onCut() = true;
	}

	for (int i = 0; i < (int)m1.seams.size(); ++i)
	{
		int k = 0;
		Facet_handle f = cgalMesh.face(m1.seams[i].faceidx);
		Halfedge_handle h = f->halfedge();

		while (h->vertex()->index() != m1.seams[i].posidx)
		{
			h = h->next();
			++k;
		}
		h->isCut() = true;
		h->vertex()->onCut() = true;
		h->rotation() = m1.seams[i].rot;
	}
}
