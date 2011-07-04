#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList0D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList1DHomogeneous2D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: ./FldCalcBCs  meshfile fieldfile\n");
        exit(1);
    }

    LibUtilities::CommSharedPtr vComm
            = LibUtilities::GetCommFactory().CreateInstance("Serial",argc,argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraph graph;
    SpatialDomains::MeshGraphSharedPtr graphShPt = graph.Read(meshfile);
    //----------------------------------------------
    // Also read and store the boundary conditions
    SpatialDomains::MeshGraph *meshptr = graphShPt.get();
    SpatialDomains::BoundaryConditionsSharedPtr boundaryConditions;        
    boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                                        ::AllocateSharedPtr(meshptr);
    boundaryConditions->Read(meshfile);    

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
	bool useFFT = false;
    //----------------------------------------------


    //----------------------------------------------
    // Define Expansion
    int expdim   = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields);

    switch(expdim)
    {
    case 1:
        {
            SpatialDomains::MeshGraph1DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph1D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }

            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 2,"NumHomogeneousDir is only set up for 1 or 2");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

                // Define Homogeneous expansion
                int nplanes = fielddef[0]->m_numModes[1];

                // choose points to be at evenly spaced points at
                const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[1],nplanes,Pkey);
                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];

                Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(vComm,Bkey,ly,useFFT,*mesh);
                Exp[0] = Exp2DH1;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(*Exp2DH1);
                }
            }
			else if(fielddef[0]->m_numHomogeneousDir == 2)
            {
                MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;
				
                // Define Homogeneous expansion
                int nylines = fielddef[0]->m_numModes[1];
				int nzlines = fielddef[0]->m_numModes[2];
				
                // choose points to be at evenly spaced points at
                const LibUtilities::PointsKey PkeyY(nylines,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(fielddef[0]->m_basis[1],nylines,PkeyY);
				
				const LibUtilities::PointsKey PkeyZ(nzlines,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);
                
				NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
				NekDouble lz = fielddef[0]->m_homogeneousLengths[1];
				
                Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(vComm,BkeyY,BkeyZ,ly,lz,useFFT,*mesh);
                Exp[0] = Exp3DH2;
				
                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(*Exp3DH2);
                }
            }
            else
            {
                MultiRegions::ExpList1DSharedPtr Exp1D;
                Exp1D = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vComm,*mesh);
                Exp[0] = Exp1D;
                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*Exp1D);
                }
            }
        }
        break;
    case 2:
        {
            SpatialDomains::MeshGraph2DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }

            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                // Define Homogeneous expansion
                int nplanes = fielddef[0]->m_numModes[2];

                // choose points to be at evenly spaced points at
                // nplanes + 1 points
                const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

				Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vComm,Bkey,lz,useFFT,*mesh);
                Exp[0] = Exp3DH1;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(*Exp3DH1);
                }
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vComm,*mesh);
                Exp[0] =  Exp2D;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*Exp2D);
                }
            }
        }
        break;
    case 3:
        {
            SpatialDomains::MeshGraph3DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph3D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamic cast failed");
            }

            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(vComm,*mesh);
            Exp[0] =  Exp3D;

            for(i = 1; i < nfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(*Exp3D);
            }
        }
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Copy data from file
    for(j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j]);
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------
    
    
    
    
    
    
    
}    