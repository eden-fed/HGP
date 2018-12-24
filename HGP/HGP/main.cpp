#include "HGP.h"
#include <direct.h>
#include "Utils\MatlabInterface.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "FastHGP.h"

#define BUFFER_SIZE 4000
void setMatlabPath();
void getMeshFromGui(int& methodIndex, std::string& objPath, std::string& vfPath);


void main()
{
	setMatlabPath();

	int methodIndex;
	std::string objPath, vfPath;
	getMeshFromGui(methodIndex, objPath, vfPath);

	if (methodIndex == 1){
		HGP program;
		program.run(objPath, vfPath);
	}
	else{
		FastHGP program;
		program.run(objPath, vfPath);
	}
}

void setMatlabPath()
{
	char path[BUFFER_SIZE];
	_fullpath(path, "..\\", BUFFER_SIZE);
	int index = strlen(path);
	//path[index - 4] = '\0';
	std::string matlabPath(path);
	matlabPath += "MatlabScripts\\";
	MatlabInterface::GetEngine().AddScriptPath(matlabPath.c_str());

	std::string matlabPath2 =	matlabPath + "FastHGP\\";
	MatlabInterface::GetEngine().AddScriptPath(matlabPath2.c_str());
	matlabPath2 = matlabPath + "HGP\\";
	MatlabInterface::GetEngine().AddScriptPath(matlabPath2.c_str());
}

void getMeshFromGui(int& methodIndex, std::string& objPath, std::string& vfPath)
{
	//get from gui
	MatlabInterface::GetEngine().Eval("LoadMesh");
	objPath = MatlabInterface::GetEngine().EvalToString("fprintf(objLocation)");
	vfPath = MatlabInterface::GetEngine().EvalToString("fprintf(vfLocation)");
	GMMDenseColMatrix methodIndexGMM(1, 1);
	MatlabGMMDataExchange::GetEngineDenseMatrix("methodIndex", methodIndexGMM);
	methodIndex = methodIndexGMM(0, 0);

	MatlabInterface::GetEngine().EvalToCout("temp=strsplit(objLocation,'\\');meshName=temp{end};meshName=meshName(1:end-4);");
	MatlabInterface::GetEngine().Eval("clear objLocation vfLocation methodIndex temp");
}