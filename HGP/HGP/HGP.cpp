#include "HGP.h"
#include "Utils/MatlabInterface.h"

#include "Borders.h"
#include "Utils/MatlabGMMDataExchange.h"

#include "FrameFixingClass.h"
#include "CGAL\CGAL_Macros.h"

#include <CGAL\Timer.h>

bool HGP::run(std::string& objPath, std::string& vfPath)
{

	bool loadFlag = loadMesh(objPath, vfPath);
	if (!loadFlag){
		return false;
	}

	Parser::setMeshAdditionalData(mCgalMesh, mMeshBuffer);
	transferMeshToMatlab();
	setHalfEdgesMap();
	setRotationsConstraints();
	setBoundaryFaces();

	MatlabInterface::GetEngine().Eval("HGP_settings");
	GMMDenseColMatrix gmmItNum(1, 1), gmmFrameFix(1, 1);
	MatlabGMMDataExchange::GetEngineDenseMatrix("maxIt", gmmItNum);
	MatlabGMMDataExchange::GetEngineDenseMatrix("useFrameFixing", gmmFrameFix);
	MatlabInterface::GetEngine().Eval("clear maxIt useFrameFixing");

	int maxIt = gmmItNum(0, 0);
	useFrameFixing = gmmFrameFix(0, 0) == 1;

	bool firstIteration = true;
	mUVs.resize(mSizeOfSystemVar, 2);
	MatlabInterface::GetEngine().EvalToCout("HGP.FrameFix.use = 0;");
	std::vector<Facet_handle> flippedTriangles;
	while (maxIt > 0)
	{
		maxIt--;
		updateHalfEdgesMetric(firstIteration);
		setHarmonicInternalVerticesConstraints();
		setHarmonicSeamVerticesConstraints();
		setFramesInMatlab(firstIteration); 
		MatlabInterface::GetEngine().EvalToCout("HGP_iteration"); 
		MatlabGMMDataExchange::GetEngineDenseMatrix("HGP.UV", mUVs);
		updatHalfedgeUVs();
		checkForFoldovers(flippedTriangles);
		if (maxIt == 0)
			break;
		FrameFixingClass frameFix(mCgalMesh,mMeshBuffer,F);
		frameStatus = frameFix.runFrameFixingProcedure(useFrameFixing);
		firstIteration = false;
		if (!frameStatus)
			MatlabInterface::GetEngine().EvalToCout("HGP.FrameFix.use = 1;");
		if ((flippedTriangles.size() == 0) && (frameStatus))
			break;	
	}
	
	if (flippedTriangles.size() == 0)
	{
		MatlabInterface::GetEngine().EvalToCout("HGP.Result.flips = [];"); 
	}
	else
	{
		GMMDenseColMatrix flips(flippedTriangles.size(), 1);
		for (int i = 0; i < flippedTriangles.size(); ++i)
			flips(i, 0) = flippedTriangles[i]->index();
		MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.Result.flips", flips);
	}

	int dummy;
	coneAngleDetection(dummy, dummy);
	calcDistortion();
	
	Borders borderClass(mCgalMesh);
	GMMDenseColMatrix genus(1, 1), numBorders(1, 1);
	genus(0,0) = borderClass.genus();
	numBorders(0, 0) = borderClass.numBorders();
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.Result.genus", genus);
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.Result.numBorders", numBorders);

	visualize();
	
	MatlabInterface::GetEngine().EvalToCout("HGP_report");

	return true;
}

bool HGP::loadMesh(std::string& objPath, std::string& vfPath)
{
	MatlabInterface::GetEngine().EvalToCout("HGP.meshName=meshName;");
	MatlabInterface::GetEngine().EvalToCout("clear meshName");

	MeshBuffer m2;
	CGAL::Timer loadObjTimer;
	loadObjTimer.start();
	bool objReadFlag = Parser::loadOBJ(objPath.c_str(), &mMeshBuffer, &m2);
	loadObjTimer.stop();
	if (!objReadFlag)
		return false;
	std::cout << "Total time to load the mesh: " << loadObjTimer.time() << " seconds.\n";

	// -----create cgal mesh-------
	std::vector<Point_3> pVec;
	pVec.resize(mMeshBuffer.positions.size());
	for (int i = 0; i < (int)pVec.size(); ++i)
		pVec[i] = Point_3(mMeshBuffer.positions[i][0], mMeshBuffer.positions[i][1], mMeshBuffer.positions[i][2]);
	std::vector<unsigned int> faces;
	faces.resize((int)mMeshBuffer.idx_pos.size());
	for (int i = 0; i < (int)faces.size(); ++i)
		faces[i] = mMeshBuffer.idx_pos[i];
	buildTriangleMesh(pVec, faces, mCgalMesh);
	// -----------------------------

	// -----load vector field or mat file with frames-------
	const std::string ffieldPostfix = ".ffield";
	const std::string matPostfix = ".mat";

	if (vfPath.size() >= ffieldPostfix.size() && vfPath.compare(vfPath.size() - ffieldPostfix.size(), ffieldPostfix.size(), ffieldPostfix) == 0){//check if this is a vector field file
		calcFramesFromVecField = true;
		MatlabInterface::GetEngine().EvalToCout("HGP.calcFramesFromVecField=true;");
		bool res = Parser::loadVectorField(vfPath.c_str(), mCgalMesh);
		if (!res){
			cout << "problem reading the file vector field file" << endl;
			return false;
		}
	}
	else if (vfPath.size() >= matPostfix.size() && vfPath.compare(vfPath.size() - matPostfix.size(), matPostfix.size(), matPostfix) == 0){//check if this is a mat file with frames
		calcFramesFromVecField = false;
		const char *temp = vfPath.c_str();
		MatlabInterface::GetEngine().SetEngineStringArray("matLocation", 1, &temp);

		MatlabInterface::GetEngine().EvalToCout("matFileExist=exist(matLocation, 'file');");
		GMMDenseColMatrix matFileExist(1, 1);
		MatlabGMMDataExchange::GetEngineDenseMatrix("matFileExist", matFileExist);
		if (!(bool)matFileExist(0, 0)){
			cout << "The frames file does not exist" << endl;
			return false;
		}
		MatlabInterface::GetEngine().EvalToCout("load(matLocation);");
		MatlabInterface::GetEngine().EvalToCout("HGP.frames=frames;");
		MatlabInterface::GetEngine().EvalToCout("HGP.calcFramesFromVecField=false;");
		MatlabInterface::GetEngine().Eval("clear matLocation matFileExist frames");

	}
	else{
		cout << "ffiled or frames are not provided" << endl;
		return false;
	}

	return true;
}

void HGP::transferMeshToMatlab()
{
	//******* Tranfer Vertices ************
	this->N = this->mCgalMesh.size_of_vertices();
	GMMDenseColMatrix V(N, 3);
	int fN = mCgalMesh.size_of_facets(); 
	GMMDenseColMatrix vectorFieldKv1(fN, 3);
	this->halfEdges.resize(fN, 3);
	this->F.resize(mCgalMesh.size_of_facets(), 3);

	auto vIt = mCgalMesh.vertices_begin();
	while (vIt != mCgalMesh.vertices_end())
	{
		int index = vIt->index();
		auto p = vIt->point();
		V(index, 0) = p[0];
		V(index, 1) = p[1];
		V(index, 2) = p[2];
		vIt++;
	}
	//**************************************

	//********* Transfer faces and VF ****************
	auto fIt = mCgalMesh.facets_begin();
	while (fIt != mCgalMesh.facets_end())
	{
		Vertex_const_handle v[3];
		int index = fIt->index();
		fIt->getVertices(v);

		int first = this->mMeshBuffer.idx_pos[3 * index];
		int count = 0;
		for (count = 0; count < 3; ++count)
			if (v[count]->index() == first)
				break;
		F(index, 0) = v[count]->index();
		F(index, 1) = v[(count+1)%3]->index();
		F(index, 2) = v[(count+2)%3]->index();

		auto he = fIt->halfedge();
		while (this->F(index, 0) != he->vertex()->index())
			he = he->next();

		halfEdges(index, 0) = he->index();
		he = he->next();
		halfEdges(index, 1) = he->index();
		he = he->next();
		halfEdges(index, 2) = he->index();

		if (calcFramesFromVecField){
			double l = fIt->kv1().squared_length();
			vectorFieldKv1(fIt->index(), 0) = fIt->kv1().x() / l;
			vectorFieldKv1(fIt->index(), 1) = fIt->kv1().y() / l;
			vectorFieldKv1(fIt->index(), 2) = fIt->kv1().z() / l;
		}

		fIt++;
	}
	//**************************************************

	GMMDenseColMatrix conesIndices(mMeshBuffer.cones.size(), 1);
	for (int i = 0; i < mMeshBuffer.cones.size(); ++i)
		conesIndices(i, 0) = mMeshBuffer.cones[i].posidx + 1;


	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.V", V);
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.F", F);
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.conesIndices", conesIndices);
	if (calcFramesFromVecField)
		MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.vectorFieldKv1", vectorFieldKv1);

	MatlabInterface::GetEngine().Eval("HGP.F=HGP.F+1;");
}

void HGP::setHalfEdgesMap()
{
	this->setHalfEdgesInCgal();

	//update seamHalfEdgesIndices
	this->seamHalfEdgesIndices.clear();
	for (int i = 0; i < (int)mMeshBuffer.seams.size(); ++i)
	{
		Facet_handle f = mCgalMesh.face(mMeshBuffer.seams[i].faceidx);
		Halfedge_handle h = f->halfedge();
		while (h->vertex()->index() != mMeshBuffer.seams[i].posidx)
			h = h->next();
		seamHalfEdgesIndices.push_back(h->index());
	}

	this->updateIndexOfHEinSystem();

	for (int i = 0; i < (int)mCgalMesh.size_of_facets(); ++i)
	{
		halfEdges(i, 0) = this->mIndexOfHEinSystem[(int)halfEdges(i, 0)];
		halfEdges(i, 1) = this->mIndexOfHEinSystem[(int)halfEdges(i, 1)];
		halfEdges(i, 2) = this->mIndexOfHEinSystem[(int)halfEdges(i, 2)];
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.halfEdges", halfEdges);
}

void HGP::setRotationsConstraints()
{
	GMMSparseRowMatrix rotMatrix(mMeshBuffer.seams.size(), 2 * mSizeOfSystemVar);

	double cosAngle[4], sinAngle[4];	//right now the code support only rotation of pi/2. you can extend below (instead of using this array) 
	cosAngle[0] = 1; sinAngle[0] = 0;
	cosAngle[1] = 0; sinAngle[1] = -1;
	cosAngle[2] = -1; sinAngle[2] = 0;
	cosAngle[3] = 0; sinAngle[3] = 1;

	/*cosAngle[4] = -0.5; sinAngle[4] = -1 * sqrt(0.75);
	cosAngle[5] = -0.5; sinAngle[5] = sqrt(0.75);*/

	for (int i = 0; i < (int)mMeshBuffer.seams.size()/2; ++i)
	{
		mRotationConstraints(i, 0) = this->mIndexOfHEinSystem[mRotationConstraints(i, 0)];
		mRotationConstraints(i, 1) = this->mIndexOfHEinSystem[mRotationConstraints(i, 1)];
		mRotationConstraints(i, 3) = this->mIndexOfHEinSystem[mRotationConstraints(i, 3)];
		mRotationConstraints(i, 4) = this->mIndexOfHEinSystem[mRotationConstraints(i, 4)];

		//if ((i % 2 == 0))
		//{
		rotMatrix(2*i, (int)mRotationConstraints(i, 0) - 1) += 1;
		rotMatrix(2*i, (int)mRotationConstraints(i, 1) - 1) += -1;
		rotMatrix(2*i, (int)mRotationConstraints(i, 4) - 1) += -1 * cosAngle[(int)mRotationConstraints(i, 2)];
		rotMatrix(2*i, (int)mRotationConstraints(i, 3) - 1) += cosAngle[(int)mRotationConstraints(i, 2)];
		rotMatrix(2*i, (int)mRotationConstraints(i, 4) - 1 + mSizeOfSystemVar) += sinAngle[(int)mRotationConstraints(i, 2)];
		rotMatrix(2*i, (int)mRotationConstraints(i, 3) - 1 + mSizeOfSystemVar) += -1 * sinAngle[(int)mRotationConstraints(i, 2)];

		rotMatrix(2*i + 1, (int)mRotationConstraints(i, 0) - 1 + mSizeOfSystemVar) += 1;
		rotMatrix(2*i + 1, (int)mRotationConstraints(i, 1) - 1 + mSizeOfSystemVar) += -1;
		rotMatrix(2*i + 1, (int)mRotationConstraints(i, 4) - 1) += -1 * sinAngle[(int)mRotationConstraints(i, 2)];
		rotMatrix(2*i + 1, (int)mRotationConstraints(i, 3) - 1) += sinAngle[(int)mRotationConstraints(i, 2)];
		rotMatrix(2*i + 1, (int)mRotationConstraints(i, 4) - 1 + mSizeOfSystemVar) += -1 * cosAngle[(int)mRotationConstraints(i, 2)];
		rotMatrix(2*i + 1, (int)mRotationConstraints(i, 3) - 1 + mSizeOfSystemVar) += cosAngle[(int)mRotationConstraints(i, 2)];
		//}
		
	}

	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.rotMatrix", rotMatrix);

}

void HGP::setBoundaryFaces()
{
	std::vector<int> seamFaces, coneFaces;
	std::map<int, bool> visitedSeamFaces;
	std::vector<int> importantVertices;
	if (!mCgalMesh.is_closed())
	{
		std::vector<bool> coneList;
		coneList.resize(mCgalMesh.size_of_vertices());
		for (int i = 0; i < (int)coneList.size(); ++i)
			coneList[i] = false;
		for (int i = 0; i < (int)mMeshBuffer.cones.size(); ++i)
		{
			importantVertices.push_back(mMeshBuffer.cones[i].posidx);
			coneList[mMeshBuffer.cones[i].posidx] = true;
		}
		auto vIt = mCgalMesh.vertices_begin();
		while (vIt != mCgalMesh.vertices_end())
		{
			if ((vIt->is_border()) && (!coneList[vIt->index()]))
			{
				coneList[vIt->index()] = true;
				importantVertices.push_back(vIt->index());
			}
			vIt++;
		}
	}
	else
	{
		for (int i = 0; i < (int)mMeshBuffer.cones.size(); ++i)
			importantVertices.push_back(mMeshBuffer.cones[i].posidx);
	}

	for (int i = 0; i < (int)importantVertices.size(); ++i)
	{
		auto v = mCgalMesh.vertex(importantVertices[i]);
		auto oneRing = v->vertex_begin();
		do
		{
			if (!oneRing->is_border())
				coneFaces.push_back(oneRing->face()->index());
			oneRing++;
		} while (oneRing->index() != v->vertex_begin()->index());
	}

	GMMDenseColMatrix BVFaces(coneFaces.size(), 1);
	for (int i = 0; i < (int)coneFaces.size(); ++i)
		BVFaces(i, 0) = coneFaces[i];

	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.BVFaces", BVFaces);
	MatlabInterface::GetEngine().Eval("HGP.BVFaces=HGP.BVFaces+1;HGP.BVFaces=unique(HGP.BVFaces);");
}

void HGP::updateHalfEdgesMetric(bool firstTime)
{
	if (firstTime)
	{
		for (Halfedge_handle h = mCgalMesh.halfedges_begin(); h != mCgalMesh.halfedges_end(); h++)
			h->targetMetric() = h->length();
	}
	else
	{
		for (Halfedge_handle h = mCgalMesh.halfedges_begin(); h != mCgalMesh.halfedges_end(); h++)
		{
			std::complex<double> v1 = Complex(h->uv().x(), h->uv().y());
			std::complex<double> v2 = Complex(h->prev()->uv().x(), h->prev()->uv().y());
			double edgeLength = std::abs(v1 - v2);
			h->targetMetric() = edgeLength;
		}
	}
}

void HGP::setHarmonicInternalVerticesConstraints()
{
	GMMSparseRowMatrix W(mSizeOfSystemVar, mSizeOfSystemVar);
	std::vector<int> goodRows;
	Mesh::Halfedge_around_vertex_circulator hdgAroundV;
	std::vector<bool> vertexList;
	double sum = 0;

	vertexList.resize(mCgalMesh.size_of_vertices());
	for (int i = 0; i < (int)vertexList.size(); ++i)
		vertexList[i] = false;

	auto hIt = mCgalMesh.halfedges_begin();
	while (hIt != mCgalMesh.halfedges_end())
	{
		if ((hIt->isCut()) || (hIt->vertex()->onCut()) || (hIt->vertex()->isCone()) || (hIt->vertex()->is_border()))
		{
			hIt++;
			continue;
		}

		if (vertexList[hIt->vertex()->index()])
		{
			hIt++;
			continue;
		}

		hdgAroundV = hIt->vertex_begin();
		int index = this->mIndexOfHEinSystem[hdgAroundV->index()] - 1;
		goodRows.push_back(index + 1);
		std::vector<int> activeRows;

		do // traverse through connected edges
		{
			double weight = hdgAroundV->opposite()->meanValueWeight(true);	//this is oposite because the way i implement and because its not symetric
			W(index, this->mIndexOfHEinSystem[hdgAroundV->prev()->index()] - 1) = weight;
			activeRows.push_back(this->mIndexOfHEinSystem[hdgAroundV->prev()->index()] - 1);
			sum += weight;
			hdgAroundV++;
		} while (hdgAroundV != hIt->vertex_begin());

		vertexList[hIt->vertex()->index()] = true;
		W(index, index) = -1;

		for (int qq = 0; qq < (int)activeRows.size(); ++qq)
			W(index, activeRows[qq]) = W(index, activeRows[qq]) / sum;

		sum = 0;
		activeRows.clear();
		hIt++;
	}

	GMMDenseColMatrix rowsToTake(goodRows.size(), 1);
	for (int i = 0; i < (int)goodRows.size(); ++i)
		rowsToTake(i, 0) = goodRows[i];

	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.W", W);
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.rowsToTake", rowsToTake);
}

void HGP::setHarmonicSeamVerticesConstraints()
{
	double sum = 0;
	GMMSparseComplexRowMatrix WSeamX(mSizeOfSystemVar, mSizeOfSystemVar), WSeamY(mSizeOfSystemVar, mSizeOfSystemVar);
	GMMSparseComplexRowMatrix WSeamX2(mSizeOfSystemVar, mSizeOfSystemVar), WSeamY2(mSizeOfSystemVar, mSizeOfSystemVar);
	GMMSparseComplexRowMatrix WSeamU1(mSizeOfSystemVar, mSizeOfSystemVar), WSeamU2(mSizeOfSystemVar, mSizeOfSystemVar);
	GMMSparseComplexRowMatrix WSeamV1(mSizeOfSystemVar, mSizeOfSystemVar), WSeamV2(mSizeOfSystemVar, mSizeOfSystemVar);

	std::vector<int> dontTake;
	std::map<int, bool> visitedHE, conesVisit;

	int rowNumber = 0;

	for (int i = 0; i < (int)seamHalfEdgesIndices.size(); ++i)
	{
		auto curH = mCgalMesh.halfEdge(seamHalfEdgesIndices[i]);
		if ((curH->vertex()->isCone()) || (visitedHE.count(curH->index()) == 1) || (curH->vertex()->is_border()))
		{
			if ((curH->vertex()->isCone() || (curH->vertex()->is_border())) && (conesVisit.count(curH->vertex()->index()) == 0))
			{
				dontTake.push_back(mIndexOfHEinSystem[curH->index()]);
				conesVisit[curH->vertex()->index()] = true;
			}
			continue;
		}

		std::vector<int> activeInside, activeRotated;

		int index = this->mIndexOfHEinSystem[curH->index()] - 1;
		auto oneRing = curH->vertex_begin();
		while (oneRing->index() != curH->index())
			oneRing++;

		auto start = oneRing;
		bool endFlag = false;
		sum = 0;
		while (!endFlag)
		{
			double weight = 1;
			weight = oneRing->opposite()->meanValueWeight(true);	//this is oposite because the way i implement and because its not symetric
			WSeamX(rowNumber, this->mIndexOfHEinSystem[oneRing->prev()->index()] - 1) = weight;
			WSeamY(rowNumber, this->mIndexOfHEinSystem[oneRing->prev()->index()] - 1) = weight;
			sum += weight;
			activeInside.push_back(this->mIndexOfHEinSystem[oneRing->prev()->index()] - 1);

			oneRing++;
			if (oneRing->isCut())
				endFlag = true;
		}

		visitedHE[oneRing->index()] = true;

		WSeamX(rowNumber, index) = -1 * sum;
		WSeamY(rowNumber, index) = -1 * sum;
		activeInside.push_back(index);

		int last_rx = 1, last_ry = 0;
		while (oneRing->index() != start->index())
		{
			int r = oneRing->rotation();
			int rx, ry;
			switch (r)
			{
			case 0: rx = 1; ry = 0;
				break;
			case 1: rx = 0; ry = 1;
				break;
			case 2: rx = -1; ry = 0;
				break;
			case 3: rx = 0; ry = -1;
				break;
			}
			int temprx = rx, tempry = ry;
			rx = temprx*last_rx - tempry*last_ry;
			ry = temprx*last_ry + tempry*last_rx;
			last_rx = rx;
			last_ry = ry;

			bool flag1 = true;

			while (flag1)
			{
				double weight = oneRing->opposite()->meanValueWeight(true);	//this is oposite because the way i implement and because its not symetric
				int i1 = this->mIndexOfHEinSystem[oneRing->prev()->index()] - 1;
				int i2 = this->mIndexOfHEinSystem[oneRing->index()] - 1;

				WSeamU1(rowNumber, i1) += weight*rx;
				WSeamU1(rowNumber, i2) += -1 * weight*rx;
				WSeamV1(rowNumber, i1) += -1 * weight*ry;
				WSeamV1(rowNumber, i2) += weight * ry;

				WSeamU2(rowNumber, i1) += weight*ry;
				WSeamU2(rowNumber, i2) += -1 * weight*ry;
				WSeamV2(rowNumber, i1) += weight*rx;
				WSeamV2(rowNumber, i2) += -1 * weight*rx;

				sum += weight;
				activeRotated.push_back(i1);
				activeRotated.push_back(i2);

				oneRing++;
				if (oneRing->isCut())
					flag1 = false;
			}

			visitedHE[oneRing->index()] = true;
		}

		sort(activeRotated.begin(), activeRotated.end());
		activeRotated.erase(unique(activeRotated.begin(), activeRotated.end()), activeRotated.end());

		for (int k = 0; k < (int)activeInside.size(); ++k)
		{
			WSeamX(rowNumber, activeInside[k]) = WSeamX(rowNumber, activeInside[k]) / sum;
			WSeamY(rowNumber, activeInside[k]) = WSeamY(rowNumber, activeInside[k]) / sum;
		}

		for (int k = 0; k < (int)activeRotated.size(); ++k)
		{
			WSeamU1(rowNumber, activeRotated[k]) = WSeamU1(rowNumber, activeRotated[k]) / sum;
			WSeamV1(rowNumber, activeRotated[k]) = WSeamV1(rowNumber, activeRotated[k]) / sum;
			WSeamU2(rowNumber, activeRotated[k]) = WSeamU2(rowNumber, activeRotated[k]) / sum;
			WSeamV2(rowNumber, activeRotated[k]) = WSeamV2(rowNumber, activeRotated[k]) / sum;
		}
		rowNumber++;
	}


	GMMDenseColMatrix coneIndices(dontTake.size(), 1), rNum(1, 1);
	bool take = true;
	int cur = 0;

	for (int i = 0; i < (int)dontTake.size(); ++i)
		coneIndices(i, 0) = dontTake[i];

	rNum(0, 0) = rowNumber;

	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.WSeamX", WSeamX);
	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.WSeamY", WSeamY);

	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.WSeamU1", WSeamU1);
	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.WSeamU2", WSeamU2);
	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.WSeamV1", WSeamV1);
	MatlabGMMDataExchange::SetEngineSparseMatrix("HGP.WSeamV2", WSeamV2);

	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.coneIndices", coneIndices);
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.rNum", rNum);
}

void HGP::setFramesInMatlab(bool firstTime)
{
	if (firstTime)
		MatlabInterface::GetEngine().EvalToCout("setFrames1");
	else
	{
		GMMDenseColMatrix frameFixingStatus(1, 1);
		frameFixingStatus(0, 0) = frameStatus;
		MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.FrameFix.status", frameFixingStatus);
		MatlabInterface::GetEngine().EvalToCout("setFrames2");
	}
}

