#include "regH2D.h"
using namespace std;

regH2D::regH2D(std::string D, std::string bM){
	Prog=D;
	baseMesh=bM;
	numTest=2;
	Geometry[0]=D+"_low.xml";
	Geometry[1]=D+"_high.xml";
	Regress=D+".ok";
	Output=D+".err";
	// DEMO and .OK PATH
	Path=regH2D::demoPath();
	okPath.assign(REG_PATH);
	okPath+="/Demo/";
	meshPath=regH2D::mPath();
}
std::string regH2D::makeCommand(int i){
	std::string prog="",mesh="",comm="";
	prog=regDemo::Path+regDemo::Prog;
	mesh=regDemo::Geometry[i];
	comm=prog+" "+mesh;//+" "+mesh;
	comm+=" > "+regDemo::Output;
	return comm;
};
