#include "HarmonicParametrization.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "CGAL\CGAL_Macros.h"

#define M_2PI		6.28318530717958647693

void HarmonicParametrization::buildTriangleMesh(const std::vector<Point_3>& vertices, const std::vector<unsigned int>& faceIndices, Mesh& mesh)
{
	//we assume here that the indices are not negative!
	MeshBuilder<HDS, Kernel> meshBuilder(&vertices, (std::vector<int>*)(&faceIndices), NULL);
	mesh.clear();
	mesh.delegate(meshBuilder);
	mesh.updateAllGlobalIndices();
}

void HarmonicParametrization::setHalfEdgesInCgal()
{
	mRotationConstraints.resize((int)mMeshBuffer.seams.size(), 5);
	std::set<int> visitedHe;
	int row = 0;
	for (int i = 0; i < (int)mMeshBuffer.seams.size(); ++i)
	{
		Facet_handle f = mCgalMesh.face(mMeshBuffer.seams[i].faceidx);
		Halfedge_handle h = f->halfedge();

		while (h->vertex()->index() != mMeshBuffer.seams[i].posidx)
			h = h->next();

		h->isCut() = true;
		h->vertex()->onCut() = true;
		h->rotation() = mMeshBuffer.seams[i].rot;

		if (visitedHe.count(h->index()) == 0)
		{

			mRotationConstraints(row, 0) = h->index();
			mRotationConstraints(row, 1) = h->prev()->index();
			mRotationConstraints(row, 2) = mMeshBuffer.seams[i].rot;
			mRotationConstraints(row, 3) = h->opposite()->index();
			mRotationConstraints(row, 4) = h->opposite()->prev()->index();
			visitedHe.insert(h->opposite()->index());
			row++;
		}
	}
}

void HarmonicParametrization::updateIndexOfHEinSystem()
{
	mIndexOfHEinSystem.resize(mCgalMesh.size_of_halfedges());
	auto heIt = mCgalMesh.halfedges_begin();
	int index1, index2;
	std::vector<bool> visit;
	visit.resize(mCgalMesh.size_of_halfedges());
	for (int i = 0; i < (int)visit.size(); ++i)
		visit[i] = false;
	//int num = 0;
	int num = 1;
	bool flag = false;
	while (heIt != mCgalMesh.halfedges_end())
	{
		if (visit[heIt->index()])
		{
			heIt++;
			continue;
		}

		if (!heIt->vertex()->onCut())//if the vertex is not on the cut, then all the halfedges are mapped to the same variable
		{
			auto oneRing = heIt->vertex()->vertex_begin();
			auto start = oneRing;
			index1 = oneRing->index();
			mIndexOfHEinSystem[index1] = num;
			visit[index1] = true;
			oneRing++;
			while (oneRing != start)
			{
				index2 = oneRing->index();
				mIndexOfHEinSystem[index2] = num;
				visit[index2] = true;
				oneRing++;
			}
			num++;
		}
		else//if it is on the cut, split the vertex
		{
			auto oneRing = heIt->vertex()->vertex_begin();
			while ((!oneRing->isCut()) && (!oneRing->opposite()->is_border()))
				oneRing++;
			auto start = oneRing;
			index1 = start->index();
			mIndexOfHEinSystem[index1] = num;
			visit[index1] = true;
			oneRing++;

			while (1)	//its o.k, the break condition have to be true
			{
				bool endFlag = false;
				while (!oneRing->isCut())
				{
					flag = true;
					index2 = oneRing->index();
					mIndexOfHEinSystem[index2] = num;
					visit[index2] = true;
					if (oneRing->is_border())
					{
						if (oneRing == start)
							endFlag = true;
						oneRing++;
						break;
					}
					oneRing++;
				}
				num++;
				if ((oneRing == start) || (endFlag))
					break;
				index1 = oneRing->index();
				mIndexOfHEinSystem[index1] = num;
				visit[index1] = true;
				oneRing++;
			}
		}
		heIt++;
	}
	mSizeOfSystemVar = num - 1;

}

//Update the UV's of the halfedges when the input array UVs is ordered such there is a vector for y values and vector for x values ([x , y])
void HarmonicParametrization::updatHalfedgeUVs()
{
	Facet_iterator faceIt = mCgalMesh.facets_begin();
	while (faceIt != mCgalMesh.facets_end())
	{
		Halfedge_handle h = faceIt->halfedge();
		h->uv() = Point_3(mUVs(mIndexOfHEinSystem[h->index()] - 1, 0), mUVs(mIndexOfHEinSystem[h->index()] - 1, 1), 0);
		h = h->next();
		h->uv() = Point_3(mUVs(mIndexOfHEinSystem[h->index()] - 1, 0), mUVs(mIndexOfHEinSystem[h->index()] - 1, 1), 0);
		h = h->next();
		h->uv() = Point_3(mUVs(mIndexOfHEinSystem[h->index()] - 1, 0), mUVs(mIndexOfHEinSystem[h->index()] - 1, 1), 0);
		faceIt++;
	}
}

void HarmonicParametrization::checkForFoldovers(std::vector<Facet_handle>& flippedTriangles)
{
	flippedTriangles.clear();

	Point_2 a, b, c;
	Facet_iterator faceIt = mCgalMesh.facets_begin();
	while (faceIt != mCgalMesh.facets_end())
	{
		Halfedge_handle h = faceIt->halfedge();
		a = Point_2(h->uv().x(), h->uv().y());
		h = h->next();
		b = Point_2(h->uv().x(), h->uv().y());
		h = h->next();
		c = Point_2(h->uv().x(), h->uv().y());
		Kernel::Triangle_2 t(a, b, c);
		if (t.orientation() <= 0){
			flippedTriangles.push_back(faceIt);
		}
		faceIt++;
	}
}
void HarmonicParametrization::visualize()
{
	std::string visMatlab = MatlabInterface::GetEngine().EvalToString("fprintf(num2str(visMatlab))");

	if (visMatlab.c_str()[0] == '0')	//save some time...
		return;
	//now we transfer information to Matlab

	if (mHasCones){
		int startConeIndex = mMeshBuffer.cones[0].posidx;
		auto startCone = mCgalMesh.vertex(startConeIndex);
		auto halfEdge = startCone->vertex_begin();
		while (!halfEdge->isCut())
			halfEdge++;

		int currentHe = halfEdge->index();
		int color = 1;;
		std::vector<pair<int, int>> seamTraverse3d, seamTraverse2d;

		int start = currentHe;

		std::map<int, bool> visitedHe;
		std::map<std::pair<int, int>, int> colorsMap3d, colorsMap2d;
		auto he = mCgalMesh.halfEdge(currentHe);
		do
		{
			currentHe = this->nextOnSeam(he->index());
			he = mCgalMesh.halfEdge(currentHe);

			if (visitedHe.count(currentHe) == 0)
			{
				seamTraverse3d.push_back(std::pair<int, int>(he->prev()->vertex()->index(), he->vertex()->index()));
				seamTraverse2d.push_back(std::pair<int, int>(this->mIndexOfHEinSystem[he->prev()->index()], this->mIndexOfHEinSystem[he->index()]));

				visitedHe[currentHe] = true;
				visitedHe[he->opposite()->index()] = true;

				colorsMap3d[std::pair<int, int>(he->prev()->vertex()->index(), he->vertex()->index())] = color;
				colorsMap2d[std::pair<int, int>(this->mIndexOfHEinSystem[he->prev()->index()], this->mIndexOfHEinSystem[he->index()])] = color;

				if ((he->vertex()->isCone()) || (degreeOnSeam(currentHe) > 2))
					color++;
			}
		} while (currentHe != start);

		do
		{
			currentHe = this->nextOnSeam(he->index());
			he = mCgalMesh.halfEdge(currentHe);

			if ((colorsMap2d.count(std::pair<int, int>(this->mIndexOfHEinSystem[he->opposite()->index()], this->mIndexOfHEinSystem[he->opposite()->prev()->index()])) == 0) && (colorsMap3d.count(std::pair<int, int>(he->prev()->vertex()->index(), he->vertex()->index())) == 1))
			{
				seamTraverse2d.push_back(std::pair<int, int>(this->mIndexOfHEinSystem[he->opposite()->index()], this->mIndexOfHEinSystem[he->opposite()->prev()->index()]));
				colorsMap2d[std::pair<int, int>(this->mIndexOfHEinSystem[he->opposite()->index()], this->mIndexOfHEinSystem[he->opposite()->prev()->index()])] = colorsMap2d[std::pair<int, int>(this->mIndexOfHEinSystem[he->prev()->index()], this->mIndexOfHEinSystem[he->index()])];
			}
		} while (currentHe != start);


		GMMDenseColMatrix seam3d(seamTraverse3d.size(), 3);
		GMMDenseColMatrix seam2d(seamTraverse2d.size(), 3);
		for (int i = 0; i < seamTraverse3d.size(); ++i)
		{
			seam3d(i, 0) = seamTraverse3d[i].first + 1;
			seam3d(i, 1) = seamTraverse3d[i].second + 1;
			seam3d(i, 2) = colorsMap3d[seamTraverse3d[i]];
		}
		for (int i = 0; i < seamTraverse2d.size(); ++i)
		{
			seam2d(i, 0) = seamTraverse2d[i].first;
			seam2d(i, 1) = seamTraverse2d[i].second;
			seam2d(i, 2) = colorsMap2d[seamTraverse2d[i]];
		}

		std::string maltabVarName = mMethodName + ".Visualize.seam3d";
		MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), seam3d);
		maltabVarName = mMethodName + ".Visualize.seam2d";
		MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), seam2d);

		GMMDenseColMatrix conesMapMatrix(20, mMeshBuffer.cones.size());
		for (int i = 0; i < (int)mMeshBuffer.cones.size(); ++i)
		{
			auto cc = mCgalMesh.vertex(mMeshBuffer.cones[i].posidx);
			auto oneRing = cc->vertex_begin();
			int row = 0;
			int col = i;
			do
			{
				conesMapMatrix(row, col) = this->mIndexOfHEinSystem[oneRing->index()];
				row++;
				oneRing++;
			} while (oneRing != cc->vertex_begin());
		}
		maltabVarName = mMethodName + ".Visualize.conesMapMatrix";
		MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), conesMapMatrix);
		/*GMMDenseColMatrix cones(mMeshBuffer.cones.size(), 1);

		for (int i = 0; i < (int)mMeshBuffer.cones.size(); ++i)
			cones(i, 0) = mMeshBuffer.cones[i].posidx + 1;	//matlab indexing
		maltabVarName = mMethodName + ".Visualize.cones";
		MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), cones);*/
	}
}

int HarmonicParametrization::nextOnSeam(int h)
{
	Halfedge_around_vertex_circulator res = mCgalMesh.halfEdge(h)->vertex_begin();
	while (res->index() != h)
		res++;
	res++;
	while (!res->isCut())
		res++;
	return(res->opposite()->index());
}

int HarmonicParametrization::degreeOnSeam(int h)
{
	int deg = 0;
	auto he = mCgalMesh.halfEdge(h);
	auto start = he->vertex_begin();
	auto hec = start;
	do
	{
		if (hec->isCut())
			++deg;
		hec++;
	} while (hec->index() != start->index());

	return (deg);
}
void HarmonicParametrization::calcDistortion()
{
	double average_k = 0.0;
	std::vector<double> k(mCgalMesh.size_of_facets(), 0.0);
	int i = 0;

	for_each_const_facet(f, mCgalMesh)
	{
		double curr_k = 0.0;
		curr_k = this->calcK(f);
		double currArea = f->area();
		average_k += currArea*curr_k;
		k[i] = curr_k;
		i++;
	}

	double min_k = *(std::min_element(k.begin(), k.end()));
	double max_k = *(std::max_element(k.begin(), k.end()));

	GMMDenseColMatrix distortionVector(k.size(), 1);
	for (int i = 0; i < (int)k.size(); ++i)
		distortionVector(i, 0) = k[i];
	std::string maltabVarName = mMethodName + ".Result.k";
	MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), distortionVector);
}

double HarmonicParametrization::calcK(Facet_const_handle& face)
{
	Complex e1, e2, e3;
	bool res = face->computeHalfedgesInLocalCoords(e1, e2, e3, false);

	if (!res) //this means that the triangle is degenerated
	{
		return 1.0;
	}

	Complex t1(imag(e1), -real(e1));
	Complex t2(imag(e2), -real(e2));
	Complex t3(imag(e3), -real(e3));

	Complex f1 = Complex(face->halfedge()->next()->uv().x(), face->halfedge()->next()->uv().y());
	Complex f2 = Complex(face->halfedge()->prev()->uv().x(), face->halfedge()->prev()->uv().y());
	Complex f3 = Complex(face->halfedge()->uv().x(), face->halfedge()->uv().y());

	Complex numerator = (f1*t1 + f2*t2 + f3*t3);
	Complex denominator = (f1*conj(t1) + f2*conj(t2) + f3*conj(t3));

	if (denominator == 0.0)
	{
		return 1.0;
	}

	double k = abs(numerator / denominator);
	return k;
}

void HarmonicParametrization::coneAngleDetection(int& numWrongAngles, int& numWrongConeAngles)
{
	std::vector<double> problem, desire;
	std::vector<int> problemIndices;
	double sum = 0;
	numWrongConeAngles = 0;

	for (Vertex_iterator v = mCgalMesh.vertices_begin(); v != mCgalMesh.vertices_end(); v++)
	{
		sum = oneRingAngle(v);
		if (v->isCone())
		{
			if (std::abs(v->getConeAngle()*M_PI - sum) > 0.01)
			{
				problemIndices.push_back(v->index());
				problem.push_back(sum);
				desire.push_back(v->getConeAngle()*M_PI);
				++numWrongConeAngles;
			}
		}
		else
		{
			if (!v->is_border())
			{
				if (std::abs(sum - M_2PI) > 0.01)
				{
					problemIndices.push_back(v->index());
					problem.push_back(sum);
					desire.push_back(M_2PI);
				}
			}
		}
	}
	if (problem.size() == 0)
	{
		std::string maltabCommand = mMethodName + ".Result.problemVerticesIndices=[];" + mMethodName + ".Result.problemVerticesAngles=[];" + mMethodName + ".Result.desiredVerticesAngles = [];";
		MatlabInterface::GetEngine().EvalToCout(maltabCommand.c_str());
	}
	else
	{
		GMMDenseColMatrix problemVerticesIndices(problem.size(), 1), problemVerticesAngles(problem.size(), 1), desiredAngles(problem.size(), 1);
		for (int i = 0; i < (int)problem.size(); ++i)
		{
			problemVerticesIndices(i, 0) = problemIndices[i];
			problemVerticesAngles(i, 0) = problem[i];
			desiredAngles(i, 0) = desire[i];
		}
		std::string maltabVarName = mMethodName + ".Result.problemVerticesIndices";
		MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), problemVerticesIndices);
		maltabVarName = mMethodName + ".Result.problemVerticesAngles";
		MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), problemVerticesAngles);
		maltabVarName = mMethodName + ".Result.desiredVerticesAngles";
		MatlabGMMDataExchange::SetEngineDenseMatrix(maltabVarName.c_str(), desiredAngles);
	}
	numWrongAngles = problem.size();

}

double HarmonicParametrization::oneRingAngle(Vertex_iterator& v)
{
	double sum = 0;
	auto oneRing = v->vertex_begin();
	do
	{
		if (oneRing->is_border())
		{
			oneRing++;
			continue;
		}

		Point_3 p0 = oneRing->uv();
		Point_3 p1 = oneRing->prev()->uv();
		Point_3 p2 = oneRing->next()->uv();

		Vector_3 u = p1 - p0;
		Vector_3 v = p2 - p0;

		double u_length = sqrt(u.squared_length());
		double v_length = sqrt(v.squared_length());

		double cosTeta = (u.x() * v.x() + u.y() * v.y()) / (u_length*v_length);
		sum += std::acos(cosTeta);

		oneRing++;

	} while (oneRing != v->vertex_begin());
	return (sum);
}
