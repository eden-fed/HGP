#include "FrameFixingClass.h"


#include "Utils/MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"


void FrameFixingClass::searchForBadCone()
{
	for (int i = 0; i < (int)meshAdditionalData.cones.size(); ++i)
	{
		auto v = cgalMesh.vertex(meshAdditionalData.cones[i].posidx);
		double ang = v->uvConeAngle(true);
		ang = ang / M_PI;
		if (std::abs(ang - v->getConeAngle()) > 0.1)
			badCones.push_back(v->index());
	}
}

void FrameFixingClass::extractOneRingAngles()
{
	auto v = cgalMesh.vertex(badCones[0]);
	auto h = v->vertex_begin();
	while (!h->isCut())
		h++;
	start = h;

	std::vector<double> oneRingAngles, n1, n2;
	std::vector<Complex> rotOffset;
	Complex rr = Complex(1, 0);
	do
	{
		if (h->is_border())
		{
			h++;
			continue;
		}
		oneRingAngles.push_back(h->uvAngle(true));
		n1.push_back(h->next()->uvAngle(true));
		n2.push_back(h->prev()->uvAngle(true));
		oneRingFaces.push_back(h->face()->index());
		rotOffset.push_back(rr);
		h++;
		if (h->isCut())
		{
			int rx, ry;
			switch (h->rotation())
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
			rr = rr*Complex(rx, ry);
		}
	} while (h != start);


	GMMDenseColMatrix badOneRing(oneRingFaces.size(), 1);
	for (int i = 0; i < (int)oneRingFaces.size(); ++i)
		badOneRing(i, 0) = oneRingFaces[i] + 1;
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.FrameFix.badOneRing", badOneRing);

	GMMDenseColMatrix originalOneRingAngle(oneRingAngles.size(), 3), oneRingConeAngle(1, 1);
	GMMDenseComplexColMatrix rotationOffset(rotOffset.size(), 1);

	oneRingConeAngle(0, 0) = v->getConeAngle()*M_PI;
	for (int i = 0; i < (int)oneRingAngles.size(); ++i)
	{
		originalOneRingAngle(i, 0) = oneRingAngles[i];//(oneRingAngles[i] / sum) * v->getConeAngle()*M_PI;
		originalOneRingAngle(i, 1) = n1[i];
		originalOneRingAngle(i, 2) = n2[i];
	}
	for (int i = 0; i < (int)oneRingFaces.size(); ++i)
		rotationOffset(i, 0) = rotOffset[i];
	 
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.FrameFix.rotationOffset", rotationOffset);
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.FrameFix.oneRingConeAngle", oneRingConeAngle);
	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.FrameFix.originalOneRingAngle", originalOneRingAngle);



	MatlabInterface::GetEngine().EvalToCout("HGP.FrameFix.newOneRingAngles = frameFixABF(HGP);");

	GMMDenseColMatrix newOneRingAngles(oneRingAngles.size(), 1);
	newAngles.resize(oneRingAngles.size());
	MatlabGMMDataExchange::GetEngineDenseMatrix("HGP.FrameFix.newOneRingAngles", newOneRingAngles);
	for (int i = 0; i < (int)oneRingAngles.size(); ++i)
		newAngles[i] = newOneRingAngles(i, 0);

}

void FrameFixingClass::setLocalEmbedding()
{
	auto h = start;
	GMMDenseComplexColMatrix embeddedOneRing(oneRingFaces.size(), 3);
	Complex I(0, 1), rotation(1, 0);
	int row = 0;
	do
	{
		if (h->is_border())
		{
			h++;
			continue;
		}
		double A = sqrt((h->uv() - h->prev()->uv()).squared_length());
		double B = sqrt((h->uv() - h->next()->uv()).squared_length());
		double C = sqrt(A*A + B*B - 2 * A*B*std::cos(newAngles[row]));
		double y = (1 / (2 * A))*(C*C - A*A - B*B);
		double x = -1 * sqrt(B*B - y*y);

		int index = h->vertex()->index();
		if (F(oneRingFaces[row], 0) == index)
		{
			embeddedOneRing(row, 0) = 0;
			embeddedOneRing(row, 1) = Complex(x, y)*rotation;
			embeddedOneRing(row, 2) = I * (double)-1 * A*rotation;
		}
		else if (F(oneRingFaces[row], 1) == index)
		{
			embeddedOneRing(row, 0) = I * (double)-1 * A*rotation;
			embeddedOneRing(row, 1) = 0;
			embeddedOneRing(row, 2) = Complex(x, y)*rotation;
		}
		else if (F(oneRingFaces[row], 2) == index)
		{
			embeddedOneRing(row, 0) = Complex(x, y)*rotation;
			embeddedOneRing(row, 1) = I * (double)-1 * A*rotation;
			embeddedOneRing(row, 2) = 0;
		}

		rotation *= Complex(cos(newAngles[row]), -1 * sin(newAngles[row]));
		row++;
		h++;
	} while (h != start);

	MatlabGMMDataExchange::SetEngineDenseMatrix("HGP.FrameFix.embeddedOneRing", embeddedOneRing);
}

bool FrameFixingClass::runFrameFixingProcedure(bool useFrameFixing)
{
	if (!useFrameFixing)	//dont do anything
		return true;
	searchForBadCone();
	if (badCones.size() != 1)
		return true;
	extractOneRingAngles();
	setLocalEmbedding();
	return false;
}