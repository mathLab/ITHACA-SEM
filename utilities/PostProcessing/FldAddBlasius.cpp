/* ===============================================================================
 * Generation of an .fld file for the Blasius boundary layer. 
 * Requirements: 
 *                a) Session file with the mesh in the physical space and some 
 *                   data to define the BL properly
 * 
 *                b) Blasius similarity solution consistent with the dimensions
 *                   of the mesh file (in the physical space)
=============================================================================== */

/* =====================================
 * Author: Gianmarco Mengaldo 
 * Generation: dd/mm/aa = 08/03/12
===================================== */ 

//! Loading cc libraries
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <iomanip>

//! Loading Nektar++ libraries
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList0D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList1DHomogeneous2D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>



//! Nektar++ namespace
using namespace Nektar;

//! Main
int main(int argc, char *argv[])
{
    //! Setting up the decimal precision to machine precision
    setprecision (16);
    
    //! Auxiliary counters for the x and y directions
    int  i, j, k, m, yLevel, xElement;
    NekDouble tmp;
    
    //! Auxiliary variables
    char    LocalString[1000];
    string  GlobalString[10000];
    bool    inspection = 1;

    //! Check for the command line
    if(argc != 3)
    {
        fprintf(stderr,"Usage: FldAddBlasius  meshfile fieldfile\n");
        exit(1);
    }

    //! Reading the session file
    LibUtilities::SessionReaderSharedPtr vSession = LibUtilities::SessionReader::CreateInstance(argc, argv);
    //! Reading the mesh from session file
    //string meshfile(argv[argc-2]);
    //SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);
    
    //! Loading the parameters to define the BL
    NekDouble Re;
    NekDouble L;
    NekDouble U_inf;
    NekDouble x;
    NekDouble nu;
    NekDouble C;
    
    vSession->LoadParameter("Re",       Re,     1.0);
    vSession->LoadParameter("L",        L,      1.0);
    vSession->LoadParameter("U_inf",    U_inf,  1.0);
    vSession->LoadParameter("x",        x,      1.0);

    
    /*  Read in mesh from input file and create an object of class 
     *  MeshGraph1D to encaplusate the mesh
     */
    SpatialDomains::MeshGraphSharedPtr graphShPt; 
    graphShPt = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
    
    /*  Feed our spatial discretisation object with the information
     *  coming from the session file i.e. initialise all the memory
     */
    MultiRegions::ContField2DSharedPtr Domain;
    Domain = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(vSession,graphShPt,vSession->GetVariable(0));
    
    //! Get the total number of elements
    int nElements;
    nElements = Domain->GetExpSize();
    std::cout << "Number of elements          = " << nElements << std::endl;
    
    /*  Get the total number of quadrature points (depends on n. modes)
     *  For nodal expansion the number of physical points is equal 
     *  to the number of coefficients
     */
    int nQuadraturePts;
    nQuadraturePts = Domain->GetTotPoints();
    std::cout << "Number of quadrature points = " << nQuadraturePts << std::endl;

    
    /*  Create 3 arrays "x_QuadraturePts", "y_QuadraturePts", "z_QuadraturePts" 
     *  to hold the coordinates values and fill it up with the the function GetCoords 
     */
    Array<OneD,NekDouble> x_QuadraturePts;
    Array<OneD,NekDouble> y_QuadraturePts;
    Array<OneD,NekDouble> z_QuadraturePts;
    x_QuadraturePts = Array<OneD,NekDouble>(nQuadraturePts);
    y_QuadraturePts = Array<OneD,NekDouble>(nQuadraturePts);
    z_QuadraturePts = Array<OneD,NekDouble>(nQuadraturePts);
    Domain->GetCoords(x_QuadraturePts,y_QuadraturePts,z_QuadraturePts);
    
    FILE *qdpoints_data; 
    qdpoints_data = fopen("qdpoints_data.txt","w+"); 
    for(j=0; j<=nQuadraturePts-1;j++)
    {
        fprintf(qdpoints_data,"%f %f %f\n", x_QuadraturePts[j], y_QuadraturePts[j], z_QuadraturePts[j]);
    }
    fclose(qdpoints_data);
    
    //! Temporary check of the coordinates 
    if(inspection == 0)
    {
        std::cout << "xqPts0 = " << x_QuadraturePts[0] << std::endl;
        std::cout << "yqPts0 = " << y_QuadraturePts[0] << std::endl;
        std::cout << "zqPts0 = " << z_QuadraturePts[0] << std::endl;
    
        std::cout << "xqPts1 = " << x_QuadraturePts[1] << std::endl;
        std::cout << "yqPts1 = " << y_QuadraturePts[1] << std::endl;
        std::cout << "zqPts1 = " << z_QuadraturePts[1] << std::endl;
    
        std::cout << "xqPts2 = " << x_QuadraturePts[2] << std::endl;
        std::cout << "yqPts2 = " << y_QuadraturePts[2] << std::endl;
        std::cout << "zqPts2 = " << z_QuadraturePts[2] << std::endl;
    
        std::cout << "xqPts4 = " << x_QuadraturePts[4] << std::endl;
        std::cout << "yqPts4 = " << y_QuadraturePts[4] << std::endl;
        std::cout << "zqPts4 = " << z_QuadraturePts[4] << std::endl;
    }
    
    //! Interface points
    int nInterfacePts;
    nInterfacePts = graphShPt->GetNvertices();
    
    //! Coordinates of the interface points
    Array<OneD,NekDouble> x_InterfacePts;
    Array<OneD,NekDouble> y_InterfacePts;
    Array<OneD,NekDouble> z_InterfacePts;
    x_InterfacePts = Array<OneD,NekDouble>(nInterfacePts);
    y_InterfacePts = Array<OneD,NekDouble>(nInterfacePts);
    z_InterfacePts = Array<OneD,NekDouble>(nInterfacePts);
    
    for (int i = 0; i < nInterfacePts; i++)
    {
        graphShPt->GetVertex(i)->GetCoords(x_InterfacePts[i],y_InterfacePts[i],z_InterfacePts[i]);
    }
    
    /* ----------------------------------------------------------------------------------------
     * Getting the number of quadrature points and the coordinates along x and y directions
     * !!! NOTE: Valid only for structured meshes !!!
     --------------------------------------------------------------------------------------- */ 
    
    //! Number of quadrature points along x
    int nQuadraturePts_x = 0;
    for (i = 0; i < nQuadraturePts; i++)
    {
        if(y_QuadraturePts[i] == 0)
        {
            nQuadraturePts_x++;
        }
    }
    std::cout<< "Quadrature Points along x    = "<< nQuadraturePts_x <<std::endl;

    //! Coordinates of the quadrature points along x (NOTE: Valid only for structured meshes)
    Array<OneD,NekDouble> x_QuadraturePts_x(nQuadraturePts_x);
    Array<OneD,NekDouble> y_QuadraturePts_x(nQuadraturePts_x);
    Array<OneD,NekDouble> z_QuadraturePts_x(nQuadraturePts_x);
    j = 0;
    for (i = 0; i < nQuadraturePts; i++)
    {
        if(y_QuadraturePts[i] == 0)
        {
            x_QuadraturePts_x[j] = x_QuadraturePts[i];
            y_QuadraturePts_x[j] = y_QuadraturePts[i];
            z_QuadraturePts_x[j] = z_QuadraturePts[i];
            j++;
        }
    }
    
    //! Number of quadrature points along y
    int nQuadraturePts_y = 0;
    for (i = 0; i < nQuadraturePts; i++)
    {
        if(x_QuadraturePts[i] == 0)
        {
            nQuadraturePts_y++;
        }    
    }
    std::cout<< "Quadrature Points along y    = "<< nQuadraturePts_y <<std::endl;

    
    //! Coordinates of the quadrature points along y (NOTE: Valid only for structured meshes)
    Array<OneD,NekDouble> x_QuadraturePts_y(nQuadraturePts_y);
    Array<OneD,NekDouble> y_QuadraturePts_y(nQuadraturePts_y);
    Array<OneD,NekDouble> z_QuadraturePts_y(nQuadraturePts_y);
    j = 0;
    for (i = 0; i < nQuadraturePts; i++)
    {
        if(x_QuadraturePts[i] == 0)
        {
            x_QuadraturePts_y[j] = x_QuadraturePts[i];
            y_QuadraturePts_y[j] = y_QuadraturePts[i];
            z_QuadraturePts_y[j] = z_QuadraturePts[i];
            j++;
        }    
    }   
    /* ----------------------------------------------------------------------------------------
     * ------------------------------------------------------------------------------------- */

    
    /*  Reading the .txt file with the similarity variable eta, and the similarity solutions:
     *  f(eta) and f'(eta) 
     */
    const char *txtfile_char;
    string txtfile(argv[argc-1]);
    txtfile_char = txtfile.c_str();
    
    ifstream pFile(txtfile_char);
    int numLines = 120;
    NekDouble d;
    NekDouble GlobalArray[numLines][3];

    for (j = 0; j <= 2; j++)
    {
        for (i=0; i<=numLines-1; i++) 
        {
            pFile >> d;
            GlobalArray[i][j] = d;
        }
    }

    //! Temporary check of the data extraction
    if(inspection == 0)
    {
        for (j = 0; j <= 2; j++)
        {
            for (i=0; i<=numLines-1; i++) 
            {
                std::cout<<"GlobalArray:   " << setprecision (16) << GlobalArray[i][j]<< std::endl;            
            }
        }
    }
    
    //! Definition of the arrays for BL computations
    Array<OneD,NekDouble> eta;
    Array<OneD,NekDouble> f;
    Array<OneD,NekDouble> df;
    
    Array<OneD,NekDouble> y;
    Array<OneD,NekDouble> u;
    Array<OneD,NekDouble> v;
    
    eta = Array<OneD,NekDouble>(numLines);
    f   = Array<OneD,NekDouble>(numLines);
    df  = Array<OneD,NekDouble>(numLines);

    y = Array<OneD,NekDouble>(numLines);
    u = Array<OneD,NekDouble>(numLines);
    v = Array<OneD,NekDouble>(numLines);

    nu = U_inf*L/Re;
    C = sqrt(2*nu*x/U_inf);
    
    /*  Saving eta, f and df in separate arrays and computing the physical variables using 
     *  BL data provided by the session file
     */
    for (i=0; i<=numLines-1; i++) 
    {
        eta[i] = GlobalArray[i][0];
        f[i]   = GlobalArray[i][1];
        df[i]  = GlobalArray[i][2];
        
        y[i]   = C*eta[i];
        u[i]   = U_inf*df[i];
        v[i]   = nu*sqrt(U_inf/(2*nu*x))*(y[i]*sqrt(U_inf/(2*nu*x))*df[i] - f[i]);
    }
    
    //! Temporary check of the physical variables
    if(inspection == 0)
    {
        for (i=0; i<=numLines-1; i++) 
        {
            std::cout<<"y:   " << setprecision (16) << y[i] << std::endl;  
            std::cout<<"u:   " << setprecision (16) << u[i] << std::endl;            
            std::cout<<"v:   " << setprecision (16) << v[i] << std::endl;            
        }
    }

    //! Reordering y grid points
    for(j=0; j<=nQuadraturePts_y-1; j++)
    {
        for(i=j+1; i<=nQuadraturePts_y-1; i++)
        {
            if(y_QuadraturePts_y[j] > y_QuadraturePts_y[i])
            {
                tmp = y_QuadraturePts_y[j];
                y_QuadraturePts_y[j] = y_QuadraturePts_y[i];
                y_QuadraturePts_y[i] = tmp;
            }
        }
    }
    
    //! Temporary check of the physical variables
    if(inspection == 0)
    {
        for(j=0; j<=nQuadraturePts_y-1; j++)
        {
            std::cout << "y_QuadPts = " << y_QuadraturePts_y[j] << std::endl;
        }
    }
    
    //! Piecewise linear interpolation of the blasius profile along y lagrangian points ------------
    k = 0;
    Array<OneD,NekDouble> uk;
    Array<OneD,NekDouble> vk;
    uk = Array<OneD,NekDouble>(nQuadraturePts_y);
    vk = Array<OneD,NekDouble>(nQuadraturePts_y);
    
    for(i=0; i<=numLines-1; i++)
    {
        while((y_QuadraturePts_y[k] >= y[i]) & (y_QuadraturePts_y[k] <= y[i+1])) 
        {
            uk[k] = (y_QuadraturePts_y[k] - y[i])*(u[i+1] - u[i])/(y[i+1] - y[i]) + u[i];
            k = k+1;
            if (k == nQuadraturePts_y + 1) break;
        }
    }
    
    
    FILE *original_data; 
    original_data = fopen("original_data.txt","w+"); 
    for(j=0; j<=numLines-1;j++)
    {
        fprintf(original_data,"%f %f\n", y[j], u[j]);
    }
    fclose(original_data);
    
    FILE *interpolation_data; 
    interpolation_data = fopen("interpolation_data.txt","w+"); 
    for(j=0; j<=nQuadraturePts_y-1;j++)
    {
        fprintf(interpolation_data,"%f %f\n", y_QuadraturePts_y[j], uk[j]);
    }
    fclose(interpolation_data);
    //! --------------------------------------------------------------------------------------------

    
    
    //! Definition of the 2D expansion using the mesh data specified on the session file ----------  
    MultiRegions::ExpList2DSharedPtr Exp2D;
    Exp2D = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vSession,graphShPt);
    //! --------------------------------------------------------------------------------------------

    
    
    //! Filling the 2D expansion using a recursive algorithm based on the mesh ordering ------------
    /*  NOTE:   1) The mesh must have a specific ordering.
     *          2) The mesh must have a BL refinement block defined using a specific criterion. 
     *          3) At the moment parallel BL only.
     *          4) No recursive algorithm on the farfield block (implicit assumption that the 
     *             entire BL (and possibly part of the farfield) will be in the BL refinement 
     *             block. 
     *          5) Only field u at the moment --> v and p are coming.
     */
    
    Array<OneD,NekDouble> ukGlobal;
    ukGlobal = Array<OneD,NekDouble>(nQuadraturePts);
    m = 0;
    int numModes = 3;
    int nElement_x = nQuadraturePts_x/(numModes+1);
    std::cout<< "nElement_x = " << nElement_x << std::endl;
    
    //! Loops on the BL refinement block
    for(xElement=0; xElement<=nElement_x-1; xElement++)
    {
        for(yLevel=0; yLevel<=55; yLevel++)
        {
            for(j=0; j<=3; j++)
            {
                ukGlobal[m] = uk[yLevel];
                m = m+1;
            }
        }
    }
    
    //! Loops on the farfield block
    for(i=m; i<=nQuadraturePts-1; i++)
    {
        ukGlobal[i] = uk[nQuadraturePts_y-1];
    }
  
    //! Temporary check of the physical variables
    if(inspection == 0)
    {
    std::cout << "i = " << i << std::endl;
    std::cout << "m = " << m << std::endl;
    std::cout << "yLevel = " << yLevel << std::endl;

    std::cout << "ukGlobal = " << ukGlobal[0] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[1] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[2] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[3] << std::endl;
    
    std::cout << "ukGlobal = " << ukGlobal[4] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[5] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[6] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[7] << std::endl;
    
    std::cout << "ukGlobal = " << ukGlobal[8] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[9] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[10] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[11] << std::endl;
    
    std::cout << "ukGlobal = " << ukGlobal[56*4-4] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4-3] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4-2] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4-1] << std::endl;
    
    std::cout << "ukGlobal = " << ukGlobal[56*4] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+1] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+2] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+3] << std::endl;
    
    std::cout << "ukGlobal = " << ukGlobal[56*4+4] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+5] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+6] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+7] << std::endl;
    
    std::cout << "ukGlobal = " << ukGlobal[56*4+8] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+9] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+10] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[56*4+11] << std::endl;
    
    std::cout << "ukGlobal = " << ukGlobal[136*88-1] << std::endl;
    std::cout << "ukGlobal = " << ukGlobal[136*88] << std::endl;
    }
    
    //! Copying the ukGlobal vector (with the same pattern of m_phys) in m_phys 
    Vmath::Vcopy(nQuadraturePts, ukGlobal, 1, Exp2D->UpdatePhys(), 1);
    
    //! Expansion coefficient extraction (necessary to write the .fld file)
    Exp2D->FwdTrans(Exp2D->GetPhys(), Exp2D->UpdateCoeffs());
    //! --------------------------------------------------------------------------------------------


    
    //! Generation .FLD file with one field only (at the moment) -----------------------------------
    string blasius = "blasius.fld";
    int nFields = 1;
    
    //! Only 1 field at the moment ===> Getting the field definition
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef = Exp2D->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    //! Writing the .fld file
    for(j = 0; j <= nFields - 1; j++)
    {
		for(i = 0; i < FieldDef.size(); i++)
		{
            FieldDef[i]->m_fields.push_back("w");
			Exp2D->AppendFieldData(FieldDef[i], FieldData[i]);
		}
    }
    graphShPt->Write(blasius, FieldDef, FieldData);
    //! --------------------------------------------------------------------------------------------

    
    
    //! WORK IN PROGRESS. Gianmarco +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    //----------------------------------------------
    // Define Expansion
    /*
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int addfields = 1;
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + addfields);
	
    MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

    // Define Homogeneous expansion
    int nplanes = fielddef[0]->m_numModes[1];

    // Choose points to be at evenly spaced points
    const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
    const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[1],nplanes,Pkey);
    NekDouble ly = fielddef[0]->m_homogeneousLengths[0];

    Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,ly,useFFT,graphShPt);
    Exp[0] = Exp2DH1;

    for(i = 1; i < nfields; ++i)
    {
        Exp[i] = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(*Exp2DH1);
    }
    //----------------------------------------------

    
    //----------------------------------------------
    // Compute gradients of fields and compute reentricity
    ASSERTL0(nfields >= 2, "Need two fields (u,v) to add reentricity");
    int nq = Exp[0]->GetNpoints();
    Array<OneD, NekDouble> grad_u[2], grad_v[2];
    Array<OneD, NekDouble> mag_cross(nq);
    for (i = 0; i < 2; ++i)
    {
    	grad_u[i] = Array<OneD,NekDouble>(nq);
    	grad_v[i] = Array<OneD,NekDouble>(nq);
    }
    Exp[0]->PhysDeriv(Exp[0]->GetPhys(), grad_u[0], grad_u[1]);
    Exp[1]->PhysDeriv(Exp[1]->GetPhys(), grad_v[0], grad_v[1]);

    // Compute cross product magnitude
    for (i = 0; i < nq; ++i)
    {
    	mag_cross[i] = grad_u[0][i] * grad_v[1][i] - grad_u[1][i] * grad_v[0][i];
    }
    Exp[nfields]->FwdTrans(mag_cross, Exp[nfields]->UpdateCoeffs());
    //----------------------------------------------

    

    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string  fldfilename(argv[2]);
    string  out = fldfilename.substr(0, fldfilename.find_last_of("."));
    string  endfile("_add.fld");
    out +=  endfile;
     
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(j = 0; j < nfields + addfields; ++j)
    {
		for(i = 0; i < FieldDef.size(); ++i)
		{
			if (j >= nfields)
			{
				FieldDef[i]->m_fields.push_back("w");
			}
			else
			{
				FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
			}
			Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
		}
    }
    graphShPt->Write(out, FieldDef, FieldData);
     */
    //-----------------------------------------------

    return 0;
}

