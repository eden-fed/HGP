#include "FastHGP.h"
#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/PardisoLinearSolver.h"
#include "CGAL\CGAL_Macros.h"
#include <CGAL\Timer.h>

#define M_2PI		6.28318530717958647693

bool FastHGP::run(std::string& objpath, std::string& vfPath)
{
	//Load mesh
	bool res = loadMesh(objpath, vfPath);
	if (!res){
		cout << "Failed to load mesh" << endl;
		return false;
	}

	//Get borders of mesh
	Borders borders(mCgalMesh);
	if (borders.genus() != 0){
		cout << "The code only supports genus 0 for now" << endl;
		return false;
	}
	if (mHasCones && borders.numBorders() > 1){
		cout << "The code does not support more than 1 border for models with cones for now" << endl;
		return false;
	}
	mNumOfBorderVertices = 0;
	for (int i = 0; i < borders.numBorders(); i++)
		mNumOfBorderVertices += borders.numVertices(i);

	//get settings from MATLAB
	getSettings();

	//Initialize mesh and arrays, and sent to MATLAB
	res = initialize(borders);
	if (!res) return false;

	//Create KKT, calculate basis, run ATP+Newton, and find all UV's
	res = runAlgorithm(borders);
	if (!res) return false;

	//Update UV's for halfedges in the mesh
	updatHalfedgeUVs();

	//Decide if failed, succeeded or partially succeeded
	res = testResult();

	sendValuesToMatlabReport();
	calcDistortion();
	visualize();

	MatlabInterface::GetEngine().EvalToCout("FastHGP_report");

	return true;
}

bool FastHGP::initialize(Borders& borders)
{
	if (mHasCones)
		setHalfEdgesInCgal();
	updateIndexOfHEinSystem();
	getBordersMapAndSetMetaVertices(borders);
	bool hasBorder = mConesAndMetaMap.size() > 0;
	if (!hasBorder && !mHasCones){
		cout << "Model must have cones or border" << endl;
		return false;
	}
	if (mHasCones)
		getConesMap();

	//number of degrees of freedom
	mNumDOF = hasBorder ? mConesAndMetaMap.size() : mConesAndMetaMap.size() - 1;

	bool res = TransferBorderFacesAndVerticesToMatlab();
	if (!res){
		cout << "Failed to transfer to MATLAB" << endl;
		return false;
	}

	return true;
}


bool FastHGP::loadMesh(std::string& objpath, std::string& vfPath)
{
	MatlabInterface::GetEngine().EvalToCout("FastHGP.meshName=meshName;");
	MatlabInterface::GetEngine().EvalToCout("clear meshName");

	//get from gui
	MeshBuffer m;
	CGAL::Timer loadObjTimer;
	loadObjTimer.start();
	bool res = Parser::loadOBJ(objpath.c_str(), &mMeshBuffer, &m);
	loadObjTimer.stop();
	if (!res) return false;
	std::cout << "Time to load the mesh: " << loadObjTimer.time() << " seconds.\n";

	// -----create cgal mesh-------
	CGAL::Timer createCGALmeshTimer;
	createCGALmeshTimer.start();
	std::vector<Point_3> pVec;
	pVec.resize(mMeshBuffer.positions.size());
	for (int i = 0; i < (int)pVec.size(); ++i)
		pVec[i] = Point_3(mMeshBuffer.positions[i][0], mMeshBuffer.positions[i][1], mMeshBuffer.positions[i][2]);
	std::vector<unsigned int> faces;
	faces.resize((int)mMeshBuffer.idx_pos.size());
	for (int i = 0; i < (int)faces.size(); ++i)
		faces[i] = mMeshBuffer.idx_pos[i];
	buildTriangleMesh(pVec, faces, mCgalMesh);
	createCGALmeshTimer.stop();
	std::cout << "Time to create CGAL mesh: " << createCGALmeshTimer.time() << " seconds.\n";

	mHasCones = mMeshBuffer.cones.size() > 0;

	if (mHasCones) {
		// -----load vector field or mat file with frames-------
		const std::string ffieldPostfix = ".ffield";
		const std::string matPostfix = ".mat";
		bool hasFrames = false;

		if (vfPath.size() >= ffieldPostfix.size() && vfPath.compare(vfPath.size() - ffieldPostfix.size(), ffieldPostfix.size(), ffieldPostfix) == 0) {//check if this is a vector field file
			mCalcFramesFromVecField = true;
			hasFrames = Parser::loadVectorField(vfPath.c_str(), mCgalMesh);
		}
		else if (vfPath.size() >= matPostfix.size() && vfPath.compare(vfPath.size() - matPostfix.size(), matPostfix.size(), matPostfix) == 0) {//check if this is a mat file with frames
			mCalcFramesFromVecField = false;
			const char *temp = vfPath.c_str();
			MatlabInterface::GetEngine().SetEngineStringArray("FastHGP.matLocation", 1, &temp);
			hasFrames = true;
		}
		if (!hasFrames) {
			cout << "ffiled or frames are not provided, ignoring the cones... " << endl;
			mMeshBuffer.cones.clear();
			mHasCones = false;
		}
		else {
			// -----update cones and seam-------
			Parser::setMeshAdditionalData(mCgalMesh, mMeshBuffer);
		}
	}

	return true;
}


void FastHGP::getSettings()
{
	//send statistics to MATLAB
	GMMDenseColMatrix numBorderVertices(1, 1), numOfVertices(1, 1), numOfFaces(1, 1), numOfCones(1, 1);

	numOfVertices(0, 0) = mCgalMesh.size_of_vertices();
	MatlabGMMDataExchange::SetEngineDenseMatrix("numV", numOfVertices);
	numOfFaces(0, 0) = mCgalMesh.size_of_facets();
	MatlabGMMDataExchange::SetEngineDenseMatrix("numF", numOfFaces);
	numOfCones(0, 0) = mMeshBuffer.cones.size();
	MatlabGMMDataExchange::SetEngineDenseMatrix("numC", numOfCones);
	numBorderVertices(0, 0) = mNumOfBorderVertices;
	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.numB", numBorderVertices);

	//call gui
	MatlabInterface::GetEngine().Eval("FastHGP_settings");

	//get setting from MATLAB
	GMMDenseColMatrix segSize(1, 1);
	MatlabGMMDataExchange::GetEngineDenseMatrix("segSize", segSize);
	mSegSize = (int)segSize(0, 0);
	GMMDenseColMatrix fixCot(1, 1);
	MatlabGMMDataExchange::GetEngineDenseMatrix("fixCot", fixCot);
	mFixCot = (fixCot(0, 0) == 1);
	MatlabInterface::GetEngine().Eval("clear fixCot segSize numV numF numC");

}

void FastHGP::getConesMap()
{
	for (int i = 0; i < (int)mMeshBuffer.cones.size(); ++i)
		mConesAndMetaMap.push_back(mCgalMesh.vertex(mMeshBuffer.cones[i].posidx));
}

//set segSize to 1 to set all boundary vertices (don't use meta vertices)
void FastHGP::getBordersMapAndSetMetaVertices(Borders& borders)
{
	mIndicesOfMetaVerticesInUVbyGeneralIndex.resize(mCgalMesh.size_of_vertices());

	mIsMeta.resize(mCgalMesh.size_of_vertices(), false);
	getVerticesThatMustBeMeta(borders);

	mMetaVerticesInBorderByHalfEdge.resize(borders.numBorders());

	//find meta vertices
	for (int i = 0; i < borders.numBorders(); i++){
		int borderSegmentSize = mSegSize;
		std::vector<Halfedge_handle>& currBorder = mMetaVerticesInBorderByHalfEdge[i];

		std::vector<Halfedge_handle> borderPath; //halfedges on boundary
		borders.getBorderHalfEdges(borders.vertex(i, 0), borderPath);

		if (borderPath.size() < 3 * mSegSize){//make sure we have at least 3 meta vertices in this border
			borderSegmentSize = borderPath.size() / 3;
		}

		int verticesInSeg = borderSegmentSize;
		for (int j = borderPath.size() - 1; j >= 0; j--){

			Halfedge_handle HD = borderPath[j];

			if (verticesInSeg >= borderSegmentSize || mIsMeta[HD->vertex()->index()]){//meta vertex
				currBorder.push_back(HD);
				verticesInSeg = 1;
				mIndicesOfMetaVerticesInUVbyGeneralIndex[HD->vertex()->index()] = mConesAndMetaMap.size();
				mConesAndMetaMap.push_back(HD->vertex());//update borders map
				mIsMeta[HD->vertex()->index()] = true;
			}
			else
				verticesInSeg++;
		}
		currBorder.push_back(currBorder[0]);
	}

}


void FastHGP::getVerticesThatMustBeMeta(Borders& borders)
{

	// mark vertices on bridges
	GMMDenseColMatrix isVertexOnBridge(mCgalMesh.size_of_vertices(), 1); //indicate if a vertex in the mesh is on a bridge
	gmm::clear(isVertexOnBridge);

	int numBorders = borders.numBorders();
	std::vector<Vertex_handle> firstBridgeVertexInBorder(numBorders, NULL);


	for (int i = 0; i < borders.numBorders(); i++){
		for (int j = 0; j < borders.numVertices(i); j++){
			Vertex_handle v = borders.vertex(i, j);

			//if a vertex is on cut and on border, is should be meta vertex
			if (v->onCut()){
				mIsMeta[v->index()] = true;
				continue;
			}

			//find bridges (edges whose deletion disconnects the graph) - and set the 2 vertices on the bridge to be meta vertices
			Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
			const Mesh::Halfedge_around_vertex_circulator hEnd = h;

			double sumWeights = 0.0;
			CGAL_For_all(h, hEnd)
			{
				if (!h->is_border_edge()){
					Vertex_handle v1 = h->vertex();
					Vertex_handle v2 = h->opposite()->vertex();
					if (v1->is_border() && v2->is_border())
					{
						isVertexOnBridge(v1->index(), 0) = 1;
						isVertexOnBridge(v2->index(), 0) = 1;

						//add bridge vertices to meta
						mIsMeta[v1->index()] = true;
						mIsMeta[v2->index()] = true;

						if (firstBridgeVertexInBorder[i] == NULL)
							firstBridgeVertexInBorder[i] = v;
					}
				}
			}
		}
	}

	//select a vertex between every 2 bridge vertices.
	for (int i = 0; i < borders.numBorders(); i++){

		if (firstBridgeVertexInBorder[i] == NULL)
			continue;

		Halfedge_around_vertex_circulator hec = firstBridgeVertexInBorder[i]->vertex_begin();
		Halfedge_around_vertex_circulator hec_end = hec;

		CGAL_For_all(hec, hec_end)
		{
			if (hec->is_border()) break;
		}
		Halfedge_handle he = hec;

		do
		{
			Halfedge_handle he_next = he->next();
			while (!mIsMeta[he_next->vertex()->index()])
				he_next = he_next->next();
			Halfedge_handle he2 = he_next;
			while (he->vertex() != he2->opposite()->vertex() && he->next()->vertex() != he2->opposite()->vertex()){
				he = he->next();
				he2 = he2->prev();
			}
			mIsMeta[he2->opposite()->vertex()->index()] = true;

			he = he_next;
		} while ((he->vertex() != firstBridgeVertexInBorder[i]));
	}
}

bool FastHGP::TransferBorderFacesAndVerticesToMatlab()
{
	//******* Tranfer duplicated Vertices by halfedges************
	GMMDenseColMatrix verticesByHalfedges(mSizeOfSystemVar, 3);
	auto heIt = mCgalMesh.halfedges_begin();
	while (heIt != mCgalMesh.halfedges_end())
	{
		int index = mIndexOfHEinSystem[heIt->index()] - 1;
		auto p = heIt->vertex()->point();
		verticesByHalfedges(index, 0) = p[0];
		verticesByHalfedges(index, 1) = p[1];
		verticesByHalfedges(index, 2) = p[2];
		heIt++;
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.V", verticesByHalfedges);

	//********* Transfer border faces corresponding to Halfedges and vector field****************
	std::unordered_set<Facet_handle> facesNearCones;//use hash to prevent duplications

	for (int i = 0; i < mConesAndMetaMap.size(); i++){
		Vertex_handle v = mConesAndMetaMap[i];
		Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
		const Mesh::Halfedge_around_vertex_circulator hEnd = h;
		CGAL_For_all(h, hEnd)
			if (!h->is_border()){
				facesNearCones.insert(h->face());
			}
	}

	int numNearConesFaces = facesNearCones.size();
	GMMDenseColMatrix facesByHalfedges(numNearConesFaces, 3);
	GMMDenseColMatrix vectorFieldKv1(numNearConesFaces, 3);

	GMMDenseColMatrix nearConesFaces(numNearConesFaces, 1);

	unordered_set<Facet_handle>::const_iterator itr;
	int i = 0;
	for (itr = facesNearCones.begin(); itr != facesNearCones.end(); ++itr){
		Facet_handle f = (*itr);
		Halfedge_handle he[3];
		f->getHalfedges(he);
		facesByHalfedges(i, 0) = mIndexOfHEinSystem[he[0]->index()];
		facesByHalfedges(i, 1) = mIndexOfHEinSystem[he[1]->index()];
		facesByHalfedges(i, 2) = mIndexOfHEinSystem[he[2]->index()];

		mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[0]->index()]-1);
		mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[1]->index()]-1);
		mConesAndNearConesMapOfRowsInKKT.insert(mIndexOfHEinSystem[he[2]->index()]-1);

		if (mHasCones){
			if (mCalcFramesFromVecField){
				double l = f->kv1().squared_length();
				vectorFieldKv1(i, 0) = f->kv1().x() / l;
				vectorFieldKv1(i, 1) = f->kv1().y() / l;
				vectorFieldKv1(i, 2) = f->kv1().z() / l;
			}
			else{
				nearConesFaces(i, 0) = f->index() + 1;
			}
		}
		++i;
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.F", facesByHalfedges);

	if (mHasCones){
		if (mCalcFramesFromVecField){
			MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.vectorFieldKv1", vectorFieldKv1);
		}
		else{//get the mat file if exists, otherwise there are no frames
			MatlabInterface::GetEngine().EvalToCout("matFileExist=exist(FastHGP.matLocation, 'file');");
			GMMDenseColMatrix matFileExist(1, 1);
			MatlabGMMDataExchange::GetEngineDenseMatrix("matFileExist", matFileExist);
			if (!(bool)matFileExist(0, 0)){
				cout << "The frames file does not exist" << endl;
				return false;
			}
			MatlabInterface::GetEngine().EvalToCout("load(FastHGP.matLocation);");
			MatlabGMMDataExchange::SetEngineDenseMatrix("nearConesFaces", nearConesFaces);
			MatlabInterface::GetEngine().EvalToCout("FastHGP.frames=frames(nearConesFaces);");
			MatlabInterface::GetEngine().Eval("clear matLocation matFileExist frames nearConesFaces");
		}
	}
	return true;
}

bool FastHGP::runAlgorithm(Borders& borders)
{
	//*** run algorithm to calculate UV's ***
	GMMDenseColMatrix RHS;
	int conesConstraintsStartRow;

	CGAL::Timer totalTime;
	totalTime.start();

	constructKKTmatrix(borders, conesConstraintsStartRow);

	bool res = constructHarmonicBasisAndSendToMATLAB(conesConstraintsStartRow);
	if (!res){
		cout << "Failed to construct Harmonic Basis" << endl;
		return false;
	}

	bool ATPSucsess = mHasCones ? getATPInitialValue() : getTutteInitialValue(borders);
	bool NewtonSucsess = runNewton(conesConstraintsStartRow, RHS);
	mATPandNewtonPassed = ATPSucsess&&NewtonSucsess;

	totalTime.stop();
	mTotalTime = totalTime.time();
	std::cout << "Total time: " << mTotalTime << "\n";

	//There is a bug in PARDISO - solve is not working after selective inverse in some cases (for example - if there are zeros in the diagonal)
	//And even if it works it is not accurate.
	//So we calculate it in MATLAB, but we don't count this time, as it takes much more time in MATLAB.
	//If PARDISO would worked as expected, we only need to perform one back substitution (the factorization is already calculated before the selective inverse phase), and this time is negligible.
	res = getAllUVsByRHS(RHS);
	if (!res){
		cout << "Failed to Calculate UVs" << endl;
		return false;
	}

	return true;
}


void FastHGP::constructKKTmatrix(Borders& borders, int& conesConstraintsStartRow)
{
	std::vector<Eigen::Triplet<double>> tripletListValues;

	tripletListValues.reserve(6 * mCgalMesh.size_of_vertices() + mNumDOF);

	FillLaplacianInKKT(tripletListValues);

	conesConstraintsStartRow = mSizeOfSystemVar;
	if (mHasCones){
		conesConstraintsStartRow = 2 * mSizeOfSystemVar;
		FillRotationConstraintsInKKT(tripletListValues, conesConstraintsStartRow);
	}

	FillMetaVerticesConstraintsInKKT(borders, tripletListValues, conesConstraintsStartRow);

	mSizeOfMatrix = mHasCones ? conesConstraintsStartRow + 2 * mNumDOF : conesConstraintsStartRow + mNumDOF;

	FillConesConstraintsInKKT(tripletListValues, conesConstraintsStartRow);
	SetElementsForPARDISO(tripletListValues, conesConstraintsStartRow);


	mKKtEigen.resize(mSizeOfMatrix, mSizeOfMatrix);
	mKKtEigen.setFromTriplets(tripletListValues.begin(), tripletListValues.end());
	mKKtEigen.makeCompressed();
}


void FastHGP::FillLaplacianInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues)
{
	//fill the matrix by going on the triangles one by one
	for_each_facet(f, mCgalMesh)
	{
		Halfedge_handle h1 = f->halfedge();
		Halfedge_handle h2 = h1->next();
		Halfedge_handle h3 = h2->next();

		double term1 = -0.5*h1->cot(false);
		double term2 = -0.5*h2->cot(false);
		double term3 = -0.5*h3->cot(false);

		//v1
		updateTermInLaplacianByHalfEdge(h3, term3, tripletListValues); //v1->v2 and v2->v1

		//v2
		updateTermInLaplacianByHalfEdge(h1, term1, tripletListValues); //v2->v3 and v3->v2

		//v3
		updateTermInLaplacianByHalfEdge(h2, term2, tripletListValues); //v3->v1 and v1->v3
	}
}


void FastHGP::updateTermInLaplacianByHalfEdge(Halfedge_handle h, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues)
{
	const int i = mIndexOfHEinSystem[h->prev()->index()] - 1;
	const int j = mIndexOfHEinSystem[h->index()] - 1;

	if (term != 0){
		updateTermInLaplacianByIndices(i, j, term, tripletListValues);
		if (mHasCones){
			updateTermInLaplacianByIndices(mSizeOfSystemVar + i, mSizeOfSystemVar + j, term, tripletListValues);
		}
	}
}


void FastHGP::updateTermInLaplacianByIndices(int i, int j, const double& term, std::vector<Eigen::Triplet<double>>& tripletListValues)
{
	tripletListValues.push_back(Eigen::Triplet<double>(i, i, -term));
	tripletListValues.push_back(Eigen::Triplet<double>(j, j, -term));
	// only store upper triangular part
	if (j >= i){ //v1->v2
		tripletListValues.push_back(Eigen::Triplet<double>(i, j, term));
	}
	else {//v2->v1
		tripletListValues.push_back(Eigen::Triplet<double>(j, i, term));
	}
}


void FastHGP::FillRotationConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT)
{
	int seamSize = mRotationConstraints.nrows();
	double cosAngle[4], sinAngle[4];	//right now the code support only rotation of pi/2. you can extend below (instead of using this array) 
	cosAngle[0] = 1; sinAngle[0] = 0;
	cosAngle[1] = 0; sinAngle[1] = -1;
	cosAngle[2] = -1; sinAngle[2] = 0;
	cosAngle[3] = 0; sinAngle[3] = 1;

	for (int i = 0; i < seamSize / 2; ++i)
	{
		int v0 = mIndexOfHEinSystem[mRotationConstraints(i, 0)] - 1;
		int v1 = mIndexOfHEinSystem[mRotationConstraints(i, 1)] - 1;
		int v3 = mIndexOfHEinSystem[mRotationConstraints(i, 3)] - 1;
		int v4 = mIndexOfHEinSystem[mRotationConstraints(i, 4)] - 1;
		int rot = mRotationConstraints(i, 2);

		// if we use Eigen, then we only fill upper triangular part, then only A^t is relevant
		int curRow = rowInKKT + 2 * i;

		tripletListValues.push_back(Eigen::Triplet<double>(v0, curRow, 1));
		tripletListValues.push_back(Eigen::Triplet<double>(v1, curRow, -1));
		tripletListValues.push_back(Eigen::Triplet<double>(v4, curRow, -1 * cosAngle[rot]));
		tripletListValues.push_back(Eigen::Triplet<double>(v3, curRow, cosAngle[rot]));
		tripletListValues.push_back(Eigen::Triplet<double>(v4 + mSizeOfSystemVar, curRow, sinAngle[rot]));
		tripletListValues.push_back(Eigen::Triplet<double>(v3 + mSizeOfSystemVar, curRow, -1 * sinAngle[rot]));

		curRow = curRow + 1;

		tripletListValues.push_back(Eigen::Triplet<double>(v0 + mSizeOfSystemVar, curRow, 1));
		tripletListValues.push_back(Eigen::Triplet<double>(v1 + mSizeOfSystemVar, curRow, -1));
		tripletListValues.push_back(Eigen::Triplet<double>(v4, curRow, -1 * sinAngle[rot]));
		tripletListValues.push_back(Eigen::Triplet<double>(v3, curRow, sinAngle[rot]));
		tripletListValues.push_back(Eigen::Triplet<double>(v4 + mSizeOfSystemVar, curRow, -1 * cosAngle[rot]));
		tripletListValues.push_back(Eigen::Triplet<double>(v3 + mSizeOfSystemVar, curRow, cosAngle[rot]));
	}

	rowInKKT = rowInKKT + seamSize;
}


void FastHGP::FillMetaVerticesConstraintsInKKT(Borders& borders, std::vector<Eigen::Triplet<double>>& tripletListValues, int& rowInKKT)
{
	//calculate weights to meta vertices
	for (int k = 0; k < borders.numBorders(); k++){
		std::vector<Halfedge_handle>& currBorder = mMetaVerticesInBorderByHalfEdge[k];
		for (int i = 0; i < currBorder.size() - 1; i++){
			Halfedge_handle metaVertexHD = currBorder[i];
			Halfedge_handle nextMetaVertexHD = currBorder[i + 1];

			//find meta Edge Length
			Halfedge_handle h = metaVertexHD;
			double metaEdgeLength = 0.0;
			while (h != nextMetaVertexHD){
				metaEdgeLength += h->length();
				h = h->prev();
			}

			//update meta vertices and t
			h = metaVertexHD->prev();
			double lengthFromFirstMeta = 0.0;
			while (h != nextMetaVertexHD){
				const int firstMetaVertex = mIndexOfHEinSystem[metaVertexHD->index()] - 1;
				const int secondMetaVertex = mIndexOfHEinSystem[nextMetaVertexHD->next()->opposite()->index()] - 1;//if the vertex is on cut then this gives the correct halfedge, and if not that it is the same.

				lengthFromFirstMeta += h->next()->length();
				const double t = (1 - lengthFromFirstMeta / metaEdgeLength);

				// if we use Eigen, then we only fill upper triangular part, then only A^t is relevant
				tripletListValues.push_back(Eigen::Triplet<double>(firstMetaVertex, rowInKKT, t));
				tripletListValues.push_back(Eigen::Triplet<double>(secondMetaVertex, rowInKKT, 1 - t));
				tripletListValues.push_back(Eigen::Triplet<double>(mIndexOfHEinSystem[h->index()] - 1, rowInKKT, -1));

				if (mHasCones){
					++rowInKKT;
					tripletListValues.push_back(Eigen::Triplet<double>(firstMetaVertex + mSizeOfSystemVar, rowInKKT, t));
					tripletListValues.push_back(Eigen::Triplet<double>(secondMetaVertex + mSizeOfSystemVar, rowInKKT, 1 - t));
					tripletListValues.push_back(Eigen::Triplet<double>(mIndexOfHEinSystem[h->index()] - 1 + mSizeOfSystemVar, rowInKKT, -1));
				}

				++rowInKKT;
				h = h->prev();
			}
		}
	}
}


void FastHGP::FillConesConstraintsInKKT(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow)
{
	for (int i = 0; i < mNumDOF; i++){
		int rowInKKT = conesConstraintsStartRow + i;
		int indexOfConeInSystem = mIndexOfHEinSystem[mConesAndMetaMap[i]->halfedge()->index()] - 1; // we only need to set one appearance of the cone, so set the one that correspond to cone->halfedge

		// if we use Eigen, then we only fill upper triangular part, then only B^t is relevant (this is why the rows and cols are reversed)

		//for the specific cone we want the equation x+i*y=1, convertining it to real matrix: [re(1) ,-im(1) ; im(1) re(1)] = [1, 0; 0 , 1]
		//all the other real and imaginary values of cones are 0
		tripletListValues.push_back(Eigen::Triplet<double>(indexOfConeInSystem, rowInKKT, 1));
		if (mHasCones)
			tripletListValues.push_back(Eigen::Triplet<double>(mSizeOfSystemVar + indexOfConeInSystem, mNumDOF + rowInKKT, 1));

	}
}


void FastHGP::SetElementsForPARDISO(std::vector<Eigen::Triplet<double>>& tripletListValues, int conesConstraintsStartRow)
{
	//set diagonal elements
	for (int i = 0; i < mSizeOfMatrix; i++)
		tripletListValues.push_back(Eigen::Triplet<double>(i, i, 0));

	//set rows by mConesAndNearConesMapOfRowsInKKT
	unordered_set<int>::const_iterator itr;
	for (itr = mConesAndNearConesMapOfRowsInKKT.begin(); itr != mConesAndNearConesMapOfRowsInKKT.end(); itr++){
		int row = (*itr);
		int row2 = row + mSizeOfSystemVar;
		for (int j = 0; j < mNumDOF; j++){
			int col = j + conesConstraintsStartRow;// because we want the cols in the upper right block of the kkt matrix
			tripletListValues.push_back(Eigen::Triplet<double>(row, col, 0));
			if (mHasCones)
				tripletListValues.push_back(Eigen::Triplet<double>(row2, col, 0));
		}
	}
}


bool FastHGP::constructHarmonicBasisAndSendToMATLAB(int conesConstraintsStartRow)
{
	//**********calculate harmonic basis in PARDISO************
	bool res = calculateHarmonicBasisInPARDISO(conesConstraintsStartRow);
	if (!res) return false;

	//******* Calculate frames and initialize matrices in MATLAB **********
	if (mHasCones && mCalcFramesFromVecField){
		MatlabInterface::GetEngine().EvalToCout("FastHGP.frames = getFramesFromVectorField(FastHGP);");
	}
	MatlabInterface::GetEngine().EvalToCout("createJmatrix");
	return true;
}

bool FastHGP::calculateHarmonicBasisInPARDISO(int conesConstraintsStartRow)
{
	//**********calculate harmonic basis in PARDISO************

	int	mtype = -2;//real and symmetric indefinite
	PardisoLinearSolver pardisoSolver(mtype);

	bool res = pardisoSolver.init(true);
	if (!res) return false;
	res = pardisoSolver.createPardisoFormatMatrix(mKKtEigen);
	if (!res) return false;
	res = pardisoSolver.preprocess();
	if (!res) return false;

	GMMCompressed1RowMatrix selectiveInvSol(mKKtEigen.rows(), mKKtEigen.cols());
	res = pardisoSolver.selectiveInverse();
	if (!res) return false;
	res = pardisoSolver.getMatrixInGMMformat(selectiveInvSol);
	if (!res) return false;

	//********** Create Harmonic basis matrix for all vertices - **********
	mHarmonicBasisInAllV.resize(mSizeOfSystemVar, mNumDOF);

	// h_i if near boundary
	//set elements
	unordered_set<int>::const_iterator itr;
	for (itr = mConesAndNearConesMapOfRowsInKKT.begin(); itr != mConesAndNearConesMapOfRowsInKKT.end(); itr++){
		int row = (*itr);
		for (int j = 0; j < mNumDOF; j++){
			int col = j + conesConstraintsStartRow;

			double realPart = selectiveInvSol(row, col);
			double imagPart = mHasCones ? selectiveInvSol(row + mSizeOfSystemVar, col) : 0;

			mHarmonicBasisInAllV(row, j) = Complex(realPart, imagPart);
		}
	}

	MatlabGMMDataExchange::SetEngineSparseMatrix("HarmonicBasis", mHarmonicBasisInAllV);
	return true;
}


bool FastHGP::getATPInitialValue()
{
	MatlabInterface::GetEngine().EvalToCout("FastHGP = ATPForInitialValue( FastHGP );");
	GMMDenseColMatrix ATPSucsess(1, 1);
	MatlabGMMDataExchange::GetEngineDenseMatrix("FastHGP.ATPSucsess", ATPSucsess);
	MatlabInterface::GetEngine().EvalToCout("[ FastHGP.UVonCones, FastHGP.FixedIndices, FastHGP.FixedValues ] = fixFirstCone( FastHGP.UVonCones );");
	return (bool)ATPSucsess(0,0);
}


bool FastHGP::getTutteInitialValue(Borders& borders)
{
	GMMDenseColMatrix initialValue(mNumDOF * 2, 1);

	int mainBorderIndex = 0;
	Vertex_handle mainBorderVertex = borders.vertex(0, 0);
	int mainBorderVertexIndex = mainBorderVertex->index();

	/******* calculate disk for main border *******/
	std::vector<Halfedge_handle> pathMainBorder;
	borders.getBorderHalfEdges(mainBorderVertex, pathMainBorder);
	int numVerticesMainBorder = pathMainBorder.size();

	double totalBorderLength = 0.0;
	for (int i = 0; i < numVerticesMainBorder; i++)
	{
		double edgeLength = pathMainBorder[i]->length();
		totalBorderLength += edgeLength;
	}
	assert(totalBorderLength > 0.0);

	double totalAngle = 0.0;
	int numMetaVerticesMainBorder = 0;
	double angle = 0.0;
	for (int i = numVerticesMainBorder - 1; i >= 0; --i)
	{
		totalAngle += angle;
		if (mIsMeta[pathMainBorder[i]->vertex()->index()]){
			double x = cos(totalAngle);
			double y = sin(totalAngle);
			int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[pathMainBorder[i]->vertex()->index()];
			initialValue(indexInUV, 0) = x;
			initialValue(indexInUV + mNumDOF, 0) = y;
			numMetaVerticesMainBorder++;
		}
		double edgeLength = pathMainBorder[i]->length();
		angle = 2.0*M_PI*edgeLength / totalBorderLength;
	}

	if (borders.numBorders() != 1){
		solveSystemToFindInitialValueEntriesOnNonMainBorder(borders, numMetaVerticesMainBorder, initialValue);
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.UVonCones", initialValue);
	MatlabInterface::GetEngine().EvalToCout("[ FastHGP.UVonCones, FastHGP.FixedIndices, FastHGP.FixedValues ] = fixFirstCone( FastHGP.UVonCones );");
	MatlabInterface::GetEngine().EvalToCout("FastHGP.ATPSucsess=1;");

	return true;
}


void FastHGP::solveSystemToFindInitialValueEntriesOnNonMainBorder(Borders& borders, int& numMetaVerticesMainBorder, GMMDenseColMatrix& initialValue)
{
	/******** construct matrix and RHS to find entries in the initial value that correspond to other borders********/
	int numEquations = mNumDOF - numMetaVerticesMainBorder;
	GMMDenseColMatrix TutteLaplacianMatrix(numEquations, numEquations);
	GMMDenseColMatrix TutteRHS(numEquations, 2);
	gmm::clear(TutteLaplacianMatrix);
	gmm::clear(TutteRHS);

	std::vector<int> indicesOfNonMainBorderInAllBorder(mNumDOF, -1);// save indices in eliminated system
	Vertex_handle v;
	int rowInMat = 0;
	for (int i = 1; i < borders.numBorders(); i++){
		for (int j = 0; j < borders.numVertices(i); j++){
			v = borders.vertex(i, j);
			if (mIsMeta[v->index()]){
				int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[v->index()];
				indicesOfNonMainBorderInAllBorder[indexInUV] = rowInMat;
				rowInMat++;
			}
		}
	}

	rowInMat = 0;
	for (int i = 1; i < borders.numBorders(); i++){
		for (int j = 0; j < borders.numVertices(i); j++){
			v = borders.vertex(i, j);
			if (mIsMeta[v->index()]){
				Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();
				const Mesh::Halfedge_around_vertex_circulator hEnd = h;
				int valance = 0;
				CGAL_For_all(h, hEnd)
				{
					++valance;
				}
				double sumWeights = 0.0;
				CGAL_For_all(h, hEnd)
				{
					double w;
					w = 1.0 / (double)valance;
					Vertex_handle neighborVertex = h->opposite()->vertex();
					if (neighborVertex->is_border() && mIsMeta[neighborVertex->index()]){//meta vertex - another variable of the system
						int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[neighborVertex->index()];
						TutteLaplacianMatrix(rowInMat, indicesOfNonMainBorderInAllBorder[indexInUV]) -= w;
					}
					else{//near border or border non meta - need to use harmonic basis
						for (int k = 0; k < mNumDOF; k++){
							double hb = real(mHarmonicBasisInAllV(mIndexOfHEinSystem[neighborVertex->halfedge()->index()] - 1, k));//there are no cones and cut, so the basis is real and not complex.
							if (indicesOfNonMainBorderInAllBorder[k] != -1){
								TutteLaplacianMatrix(rowInMat, indicesOfNonMainBorderInAllBorder[k]) -= w*hb;
							}
							else{
								TutteRHS(rowInMat, 0) += w*hb*initialValue(k, 0);
								TutteRHS(rowInMat, 1) += w*hb*initialValue(k + mNumDOF, 0);
							}
						}
					}
					sumWeights += w;
				}
				TutteLaplacianMatrix(rowInMat, rowInMat) += 1.0;
				rowInMat++;
			}
		}
	}

	/******** solve linear system ********/
	MatlabGMMDataExchange::SetEngineDenseMatrix("LaplacianMatrixForInitialVaue", TutteLaplacianMatrix);
	MatlabGMMDataExchange::SetEngineDenseMatrix("RHSforInitialValue", TutteRHS);

	MatlabInterface::GetEngine().EvalToCout("x=LaplacianMatrixForInitialVaue\\RHSforInitialValue;");

	GMMDenseColMatrix initialValueNonMain(numEquations, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("x", initialValueNonMain);;

	for (int i = 1; i < borders.numBorders(); i++){
		for (int j = 0; j < borders.numVertices(i); j++){
			v = borders.vertex(i, j);
			if (mIsMeta[v->index()]){ // only calculate for meta vertex
				int indexInUV = mIndicesOfMetaVerticesInUVbyGeneralIndex[v->index()];
				double x = initialValueNonMain(indicesOfNonMainBorderInAllBorder[indexInUV], 0);
				double y = initialValueNonMain(indicesOfNonMainBorderInAllBorder[indexInUV], 1);
				initialValue(indexInUV, 0) = x;
				initialValue(indexInUV + mNumDOF, 0) = y;
			}
		}
	}

	MatlabInterface::GetEngine().EvalToCout("clear x LaplacianMatrixForInitialVaue RHSforInitialValue");
}


bool FastHGP::runNewton(int conesConstraintsStartRow, GMMDenseColMatrix& RHS)
{
	MatlabInterface::GetEngine().EvalToCout("FastHGP = Newton( FastHGP );");

	GMMDenseColMatrix NewtonSucsess(1, 1);
	MatlabGMMDataExchange::GetEngineDenseMatrix("FastHGP.NewtonSucsess", NewtonSucsess);

	GMMDenseColMatrix UVonCones(mNumDOF, 2);
	MatlabGMMDataExchange::GetEngineDenseMatrix("FastHGP.UVonCones", UVonCones);

	gmm::clear(RHS);
	if (mHasCones)
		gmm::resize(RHS, mSizeOfMatrix, 1);// [zeros(2 * sizeOfSystemVar + seamSize,1);uvsOnCones(2*numCones,1)]
	else
		gmm::resize(RHS, mSizeOfMatrix, 2);// [zeros(sizeOfSystemVar,1);uOnCones(numCones,1) , zeros(sizeOfSystemVar,1);vOnCones(numCones,1)]

	//update RHS 
	// update RHS by UVonCones in corresponding elements
	for (int i = 0; i < mNumDOF; i++){
		int row = conesConstraintsStartRow + i;
		RHS(row, 0) = UVonCones(i, 0);
		if (mHasCones)
			RHS(row + mNumDOF, 0) = UVonCones(i, 1);
		else
			RHS(row, 1) = UVonCones(i, 1);
	}

	return (bool)NewtonSucsess(0, 0);
}


bool FastHGP::getAllUVsByRHS(GMMDenseColMatrix& RHS)
{
	GMMDenseColMatrix UVs(RHS.nrows(), RHS.ncols());
	int res;

	//res = mPardisoSolver.solve(RHS, mUVs); //This is not working for some reason after selective inverse (but only if there are 0 in diagonal of matrix) - bug in PARDISO

	//***PARIDSO is not accurate so we use MATLAB***

	//create GMM matrix from Eigen:
	const Eigen::SparseMatrix<double> M = mKKtEigen;
	GMMSparseRowMatrix kktGMM(M.rows(), M.cols());

	for (unsigned k = 0; k < (unsigned)M.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
		{
			kktGMM(it.row(), it.col()) = it.value();
			kktGMM(it.col(), it.row()) = it.value();
		}
	}

	//solve in matlab:
	MatlabGMMDataExchange::SetEngineSparseMatrix("M", kktGMM);
	MatlabGMMDataExchange::SetEngineDenseMatrix("RHS", RHS);
	MatlabInterface::GetEngine().EvalToCout("[X, ~] = general_solve(M, RHS);");
	MatlabGMMDataExchange::GetEngineDenseMatrix("X", UVs);
	MatlabInterface::GetEngine().Eval("clear M RHS X");


	mUVs.resize(mSizeOfSystemVar, 2);
	for (int i = 0; i < mSizeOfSystemVar; ++i){
		mUVs(i, 0) = UVs(i, 0);
		if (mHasCones)
			mUVs(i, 1) = UVs(i + mSizeOfSystemVar, 0);
		else
			mUVs(i, 1) = UVs(i, 1);
	}

	return true;
}

bool FastHGP::testResult()
{
	int numFoldovers, numFoldsNearCones, numFoldsNearBorder, numWrongAngles, numWrongConeAngles;
	std::vector<Facet_handle> flippedTriangles;
	checkForFoldovers(flippedTriangles);

	if (mFixCot && mATPandNewtonPassed){
		fixCotFoldovers(flippedTriangles);
		checkForFoldovers(flippedTriangles);
	}
	
	checkLocationOfFoldoversAndPrint(flippedTriangles, numFoldovers, numFoldsNearCones, numFoldsNearBorder);
	coneAngleDetection(numWrongAngles, numWrongConeAngles);

	if (mATPandNewtonPassed && numFoldovers == 0 && numWrongAngles == 0)
		cout << "result: success" << endl;
	else if (mATPandNewtonPassed && numFoldsNearCones == 0 && numFoldsNearBorder == 0 && numWrongConeAngles == 0)
		cout << "result: partial success" << endl;
	else{
		cout << "result: fail" << endl;
		return false;
	}
	return true;
}

void FastHGP::checkLocationOfFoldoversAndPrint(std::vector<Facet_handle>& flippedTriangles, int& numFoldovers, int& numFoldsNearCones, int& numFoldsNearBorder)
{
	std::vector<int> flippedTrianglesNearCones;
	std::vector<int> flippedTrianglesNearCut;
	std::vector<int> flippedTrianglesNearBorder;
	std::vector<int> flippedTrianglesinterior;

	if (flippedTriangles.size() == 0)
	{
		MatlabInterface::GetEngine().EvalToCout("FastHGP.Result.flips = [];");
	}
	else
	{
		GMMDenseColMatrix flips(flippedTriangles.size(), 1);

		for (int i = 0; i < flippedTriangles.size(); ++i){
			Facet_handle faceIt = flippedTriangles[i];
			flips(i, 0) = faceIt->index();
			Halfedge_handle h = faceIt->halfedge();
			if (faceIt->is_border_face())
				flippedTrianglesNearBorder.push_back(faceIt->index());
			if (h->vertex()->onCut() || h->next()->vertex()->onCut() || h->prev()->vertex()->onCut())
				if (h->vertex()->isCone() || h->next()->vertex()->isCone() || h->prev()->vertex()->isCone())
					flippedTrianglesNearCones.push_back(faceIt->index());
				else
					flippedTrianglesNearCut.push_back(faceIt->index());
			else if (!faceIt->is_border_face())
				flippedTrianglesinterior.push_back(faceIt->index());
		}

		MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.Result.flips", flips);

		std::stringstream out;
		for (int i = 0; i < flippedTriangles.size(); ++i){
			out << flippedTriangles[i]->index() << " ";
		}
		std::string infoString(out.str());
		cout << infoString << endl;
	}

	numFoldovers = flippedTriangles.size();
	numFoldsNearCones = flippedTrianglesNearCones.size();
	numFoldsNearBorder = flippedTrianglesNearBorder.size();

}

void FastHGP::fixCotFoldovers(std::vector<Facet_handle>& flippedTriangles)
{
	if (flippedTriangles.size() == 0)
		return;
	cout << "There are " << flippedTriangles.size() << " local foldovers resulting from cot weights, fixing..." << endl;
	CGAL::Timer fixCotTime;
	fixCotTime.start();

	for (int i = 0; i < flippedTriangles.size(); ++i){
		Facet_handle faceIt = flippedTriangles[i];
		Vertex_handle v[3];
		faceIt->getVertices(v);
		for (int j = 0; j < 3; ++j){
			if (!v[j]->isCone() && !v[j]->is_border() && !v[j]->onCut())//only fix interior vertices for simplicity
			{
				double sum = oneRingAngle(v[j]);
				if (std::abs(sum - M_2PI) > 0.01)
				{
					putVertexInKernelUsingCVX(v[j]);
				}
			}
		}
	}

	fixCotTime.stop();
	std::cout << "Fix cot foldovers time: " << fixCotTime.time() << "\n";

}

void FastHGP::putVertexInKernelUsingCVX(Vertex_handle v)
{
	//get values for MATLAB function
	Mesh::Halfedge_around_vertex_circulator h = v->vertex_begin();//clockwise
	const Mesh::Halfedge_around_vertex_circulator hEnd = h;
	int numNeighbors = 0;
	CGAL_For_all(h, hEnd)
	{
		++numNeighbors;
	}

	GMMDenseColMatrix UV(2, 1);
	UV(0, 0) = v->halfedge()->uv().x();
	UV(1, 0) = v->halfedge()->uv().y();

	GMMDenseColMatrix oneRing(2, numNeighbors);
	int index = numNeighbors - 1;//I start from the end since I need it ordered counterclockwise
	CGAL_For_all(h, hEnd)
	{
		oneRing(0, index) = h->opposite()->uv().x();
		oneRing(1, index) = h->opposite()->uv().y();
		--index;
	}

	//calculate in MATLAB
	MatlabGMMDataExchange::SetEngineDenseMatrix("oneRing", oneRing);
	MatlabGMMDataExchange::SetEngineDenseMatrix("UV", UV);
	MatlabInterface::GetEngine().EvalToCout("[ newUV ] = putVertexInKernel( UV, oneRing );");
	MatlabGMMDataExchange::GetEngineDenseMatrix("newUV", UV);

	//Update halfedge UVs
	Point_3 UVp = Point_3(UV(0, 0), UV(1, 0), 0);
	CGAL_For_all(h, hEnd)
	{
		h->uv() = UVp;
	}
	//update mUVs array
	mUVs(mIndexOfHEinSystem[h->index()] - 1, 0) = UV(0, 0);
	mUVs(mIndexOfHEinSystem[h->index()] - 1, 1) = UV(1, 0);
}

void FastHGP::sendValuesToMatlabReport()
{

	//******* Tranfer Vertices ************
	GMMDenseColMatrix V(mCgalMesh.size_of_vertices(), 3);
	GMMDenseColMatrix F(mCgalMesh.size_of_facets(), 3);
	GMMDenseColMatrix UVsIndices(mCgalMesh.size_of_facets(), 3);

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

	//********* Transfer faces and UVsIndices ****************
	auto fIt = mCgalMesh.facets_begin();
	while (fIt != mCgalMesh.facets_end())
	{
		Vertex_const_handle v[3];
		int index = fIt->index();
		fIt->getVertices(v);

		int first = mMeshBuffer.idx_pos[3 * index];
		int count = 0;
		for (count = 0; count < 3; ++count)
			if (v[count]->index() == first)
				break;
		F(index, 0) = v[count]->index() + 1;
		F(index, 1) = v[(count + 1) % 3]->index() + 1;
		F(index, 2) = v[(count + 2) % 3]->index() + 1;

		auto he = fIt->halfedge();
		while (F(index, 0) - 1 != he->vertex()->index())
			he = he->next();

		UVsIndices(index, 0) = mIndexOfHEinSystem[he->index()];
		he = he->next();
		UVsIndices(index, 1) = mIndexOfHEinSystem[he->index()];
		he = he->next();
		UVsIndices(index, 2) = mIndexOfHEinSystem[he->index()];

		fIt++;
	}

	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.V", V);
	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.F", F);
	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.halfEdges", UVsIndices);
	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.UV", mUVs);

	//*********** Transfer cones indices **************
	if (mHasCones){
		GMMDenseColMatrix conesIndices(mMeshBuffer.cones.size(), 1);
		for (int i = 0; i < mMeshBuffer.cones.size(); ++i)
			conesIndices(i, 0) = mMeshBuffer.cones[i].posidx + 1;
		MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.conesIndices", conesIndices);
	}
	else{
		MatlabInterface::GetEngine().EvalToCout("FastHGP.conesIndices=[];");
	}

	GMMDenseColMatrix numDOF(1, 1);
	numDOF(0, 0) = mNumDOF;
	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.numDOF", numDOF);

	/*** Total time ***/
	GMMDenseColMatrix timeGMM(1, 1);
	timeGMM(0, 0) = mTotalTime;
	MatlabGMMDataExchange::SetEngineDenseMatrix("FastHGP.Result.totalTime", timeGMM);
}
