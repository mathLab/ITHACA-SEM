#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sys/dir.h>
#include <cmath>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion.h>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

#include <SolverUtils/SolverUtilsDeclspec.h>


using namespace Nektar;

int main(int argc, char *argv[])
{
	// Variable initialisation
	
	cout << "N1  N2  N3  lenght  Re  Tw  Tinf  Mach  beta  x_hump  x_humpin  x_humpout  h_0" << endl;
	cout << " x_hump et tt entre [0,1] " <<endl;
	cout << "h_0: vraie valeur" << endl;

	int N1=int(atof(argv[1]));
	int N2=int(atof(argv[2]));
	int N3=int(atof(argv[3]));
	
	NekDouble lenght= atof(argv[4]);
	NekDouble Re= atof(argv[5]);
	NekDouble Twall =atof(argv[6]);
	NekDouble Tinf =atof(argv[7]);
	NekDouble Mach =atof(argv[8]);
	NekDouble dx, dy ;
	
	NekDouble beta=atof(argv[9]);
	NekDouble x_hump=atof(argv[10]);
	NekDouble x_humpin=atof(argv[11]);
	NekDouble x_humpout=atof(argv[12]);
	NekDouble h_0 =atof(argv[13]);
	
	
	
	ofstream file;
	file.open("MEANBL.DAT");
		
	file << "  Goertler Swept Flat Plate Test Case 1 " << endl;
	file << "       0       0       1       0       0       0 "<< endl;
	file << "		Re		" <<   "  " << "		Mach		" <<
	 "  " << "		Sweep" <<  " " << 
	"	Points	" <<  "  " << "		Temp" <<
	 "  " << "		Chord		" <<  "  " <<
	 "		Twall		" << endl;
	 
	file <<  scientific << setprecision(10) << Re << "  " << Mach << "  " << 0.0 << "  " << int(N1+N2+N3-2) << "  " <<
	Tinf << "  " << lenght << "  " << Twall << endl;
	file << "		x/c		" <<   "  " << "		y/c		" <<
	 "  " << "		s/c		" <<  "  " << "				Cp" <<
	 "  " << "				Cq	" << endl;
	
	for (int i=0; i<N1; i++)
	{
		dx= ((x_humpin)/(N1-1))*i;
		dy=0;
		
		file << "  " << dx << "  " << dy << "  " << dx << "  " << 0.0 << "  " << 0.0 << endl;
	}
	
	for (int i=1; i<N2; i++)
	{
		dx= x_humpin + ((x_humpout-x_humpin)/(N2-1))*i;
		dy= (h_0/lenght)*exp(-beta*(dx-x_hump)*(dx-x_hump));;
		
		file << "  " << dx << "  " << dy << "  " << dx << "  " << 0.0 << "  " << 0.0 << endl;
	}
		
	for (int i=1; i<N3; i++)
	{
		dx= x_humpout + ((1-x_humpout)/(N3-1))*i;
		dy=0;
		
		file << "  " << dx << "  " << dy << "  " << dx << "  " << 0.0 << "  " << 0.0 << endl;
	}	

    return 0;
}




