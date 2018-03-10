#include <MultiRegions/ExpList1D.h>
#include <LibUtilities/Kernel/kernel.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    int j,e;
	if (argc != 2)
	{
		fprintf(stderr,"Usage: Post-processor meshfile \n"); 
		exit(1);
	}
	/************************* Mesh File ***********************/
	// Read the mesh
    SpatialDomains::MeshGraphSharedPtr graph1D = SpatialDomains::MeshGraph::Read(vSession);

	/***********************************************************/

	// Construct an object from the class ExpList1D.
        // This is the class which represents a multi-elemental expansion
	// This object can be constructed based on the input mesh
        MultiRegions::ExpList1DSharedPtr u =
            MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vSession,graph1D);

	/***********************************************************/
	// Get the number of elements
	int nelements = u->GetExpSize();
	NekDouble h = 1.0/u->GetExpSize(); // mesh spacing

	// Get the number of modes
	Array<OneD, LibUtilities::BasisSharedPtr>  base = (u->GetExp(0))->GetBase();
	int order = base[0]->GetNumModes();

	
	/****** Read in the data file and update the coefficients ***/
	// Contains the local polynomial modes for each element
	Array<OneD,NekDouble> uhat_pre = u->UpdateCoeffs();
	
	ifstream inFile;
	inFile.open("uhatPre.txt");
	double x;
	int index = 0;
	while (inFile >> x)
	{
		uhat_pre[index] = x;
		index++;
	}
	inFile.close();
	

	/************** Define the kernel object *******************/
		
	LibUtilities::KernelSharedPtr post_kernel =
		MemoryManager<LibUtilities::Kernel>::AllocateSharedPtr(order);
	post_kernel->UpdateKernelBreaks(h);
	
	/******************* Post-Processing ***********************/
	// Construct the appropriate post-processed multi-element expansion object
	int post_order = 2*order;
	LibUtilities::PointsType pType = base[0]->GetPointsType();
	LibUtilities::PointsKey  pKey(post_order+1,pType);
	LibUtilities::BasisKey upost_bkey(base[0]->GetBasisType(),post_order,pKey);
	MultiRegions::ExpList1DSharedPtr u_post =
        MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vSession,upost_bkey,graph1D);

	// The post-processing for the first element
	// Step1: Define the element ID on which the post-processing is done
	e = 0;

	// Step2: Calculate the points used for evaluating the post-processed solution
	StdRegions::StdExpansionSharedPtr elmExp = u_post->GetExp(e);
	int eval_npoints = elmExp->GetTotPoints();
	Array<OneD,NekDouble> eval_points(eval_npoints);
	elmExp->GetCoords(eval_points);

	// Step3: Call the appropriate post-processing function
	Array<OneD,NekDouble> ustar_elm(eval_npoints);
	u->PostProcess(post_kernel,eval_points,ustar_elm,h,e);
	
	// The Post-Processing for the entire domain
	// Step1: Calculate the points used for evaluating the post-processed solution
	int eval_npoints2 = u_post->GetNpoints();
	Array<OneD,NekDouble> eval_points2(eval_npoints2);
	u_post->GetCoords(eval_points2);
	
	// Step2: Call the appropriate post-processing function
	Array<OneD,NekDouble> ustar(eval_npoints2);
	u->PostProcess(post_kernel,eval_points2,ustar,h);
	
	// Calculate the post-processed coefficients for the entire domain
	Array<OneD,NekDouble> uhat_post = u_post->UpdateCoeffs();
	u_post->FwdTrans(ustar,uhat_post);
	
	ofstream postCoeff;
	postCoeff.open("postCoeff.txt");
	cout << "Printing the uhat_post: " << endl;  // These coefficients are different than my previous result, given the ustars are the same
	for(e = 0; e < nelements; e++)
	{
		for(j = 0; j < post_order; j++)
		{
			cout << uhat_post[e*post_order+j] << "  ";
			postCoeff << uhat_post[e*post_order+j] << endl;
		}
		cout << endl;
	}

}
