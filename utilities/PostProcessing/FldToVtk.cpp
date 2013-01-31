#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
using namespace Nektar;

#include <sys/stat.h>

int fexist( const char *filename ) {
  struct stat buffer ;
  if ( stat( filename, &buffer ) ) return 0 ;
  return 1 ;
}

int main(int argc, char *argv[])
{
    unsigned int i,j;

    if(argc < 3)
    {
        fprintf(stderr,"Usage: FldToVtk  meshfile  fieldfile(s)\n");
        exit(1);
    }


    int nExtraPoints;
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    vSession->LoadParameter("OutputExtraPoints",nExtraPoints,0);

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graphShPt
                                    = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    for (int n = 1; n < argc; ++n)
    {
        string fname = std::string(argv[n]);
        int fdot = fname.find_last_of('.');
        if (fdot != std::string::npos)
        {
            string ending = fname.substr(fdot);

            // If .chk or .fld we exchange the extension in the output file.
            // For all other files (e.g. .bse) we append the extension to avoid
            // conflicts.
            if (ending == ".chk" || ending == ".fld")
            {
                fname = fname.substr(0,fdot);
            }
            else if (ending == ".xml")
            {
                continue;
            }
        }
        fname = fname + ".vtu";

        if (argc > 3)
        {
            if (fexist(fname.c_str()))
            {
                cout << "Skipping converted file: " << argv[n] << endl;
                continue;
            }
            cout << "Processing " << argv[n] << endl;
        }

        //----------------------------------------------
        // Import field file.
        string fieldfile(argv[n]);
        vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
        vector<vector<NekDouble> > fielddata;
        graphShPt->Import(fieldfile,fielddef,fielddata);
        bool useFFT = false;
		bool dealiasing = false;
        //----------------------------------------------

        //----------------------------------------------
        // Set up Expansion information
        for(i = 0; i < fielddef.size(); ++i)
        {
            vector<LibUtilities::PointsType> ptype;
            for(j = 0; j < 3; ++j)
            {
                ptype.push_back(LibUtilities::ePolyEvenlySpaced);
            }
            
            fielddef[i]->m_pointsDef = true;
            fielddef[i]->m_points    = ptype; 
            
            vector<unsigned int> porder;
            if(fielddef[i]->m_numPointsDef == false)
            {
                for(j = 0; j < fielddef[i]->m_numModes.size(); ++j)
                {
                    porder.push_back(fielddef[i]->m_numModes[j]+nExtraPoints);
                }
                
                fielddef[i]->m_numPointsDef = true;
            }
            else
            {
                for(j = 0; j < fielddef[i]->m_numPoints.size(); ++j)
                {
                    porder.push_back(fielddef[i]->m_numPoints[j]+nExtraPoints);
                }
            }
            fielddef[i]->m_numPoints = porder;
            
        }
        graphShPt->SetExpansions(fielddef);
        //----------------------------------------------


        //----------------------------------------------
        // Define Expansion
        int expdim  = graphShPt->GetMeshDimension();
        int nfields = fielddef[0]->m_fields.size();
        Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields);

        switch(expdim)
        {
        case 1:
            {
                ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 2,"NumHomogeneousDir is only set up for 1 or 2");

                if(fielddef[0]->m_numHomogeneousDir == 1)
                {
                    MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

                    // Define Homogeneous expansion
                    //int nplanes = fielddef[0]->m_numModes[1];
					int nplanes; 
					vSession->LoadParameter("HomModesZ",nplanes,fielddef[0]->m_numModes[1]);

                    // choose points to be at evenly spaced points at
                    const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                    const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[1],nplanes,Pkey);
                    NekDouble ly = fielddef[0]->m_homogeneousLengths[0];

                    Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,ly,useFFT,dealiasing,graphShPt);

                    for(i = 1; i < nfields; ++i)
                    {
                        Exp[i] = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(*Exp2DH1);
                    }
                }
				else if(fielddef[0]->m_numHomogeneousDir == 2)
				{
					MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;
					
					// Define Homogeneous expansion
					//int nylines = fielddef[0]->m_numModes[1];
					//int nzlines = fielddef[0]->m_numModes[2];
					
					int nylines;
					int nzlines;
					vSession->LoadParameter("HomModesY",nylines,fielddef[0]->m_numModes[1]);
					vSession->LoadParameter("HomModesZ",nzlines,fielddef[0]->m_numModes[2]);
					
					
					// choose points to be at evenly spaced points at
					const LibUtilities::PointsKey PkeyY(nylines+1,LibUtilities::ePolyEvenlySpaced);
					const LibUtilities::BasisKey  BkeyY(fielddef[0]->m_basis[1],nylines,PkeyY);
					
					const LibUtilities::PointsKey PkeyZ(nzlines+1,LibUtilities::ePolyEvenlySpaced);
					const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);
					
					NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
					NekDouble lz = fielddef[0]->m_homogeneousLengths[1];
					
					Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(vSession,BkeyY,BkeyZ,ly,lz,useFFT,dealiasing,graphShPt);
					Exp[0] = Exp3DH2;
					
					for(i = 1; i < nfields; ++i)
					{
						Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(*Exp3DH2);
					}
				}
                else
                {
                    MultiRegions::ExpList1DSharedPtr Exp1D;
                    Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                    ::AllocateSharedPtr(vSession,graphShPt);
                    Exp[0] = Exp1D;
                    for(i = 1; i < nfields; ++i)
                    {
                        Exp[i] = MemoryManager<MultiRegions::ExpList1D>
                                                        ::AllocateSharedPtr(*Exp1D);
                    }
                }
            }
            break;
        case 2:
            {
                ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

                if(fielddef[0]->m_numHomogeneousDir == 1)
                {
                    MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                    // Define Homogeneous expansion
                    //int nplanes = fielddef[0]->m_numModes[2];
					
					int nplanes; 
					vSession->LoadParameter("HomModesZ",nplanes,fielddef[0]->m_numModes[2]);

                    // choose points to be at evenly spaced points at
                    // nplanes + 1 points
                    const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                    const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                    NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                    Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,dealiasing,graphShPt);
                    Exp[0] = Exp3DH1;

                    for(i = 1; i < nfields; ++i)
                    {
                        Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(*Exp3DH1);
                    }
                }
                else
                {
                    MultiRegions::ExpList2DSharedPtr Exp2D;
                    Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                            ::AllocateSharedPtr(vSession,graphShPt);
                    Exp[0] =  Exp2D;

                    for(i = 1; i < nfields; ++i)
                    {
                        Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                                                            ::AllocateSharedPtr(*Exp2D);
                    }
                }
            }
            break;
        case 3:
            {
                MultiRegions::ExpList3DSharedPtr Exp3D;
                Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                        ::AllocateSharedPtr(vSession,graphShPt);
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
        // Copy data from field file
        for(j = 0; j < nfields; ++j)
        {
            for(int i = 0; i < fielddata.size(); ++i)
            {
                Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                            fielddata[i],
                                            fielddef [i]->m_fields[j],
                                            Exp[j]->UpdateCoeffs());
            }
            Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
        }
        //----------------------------------------------

        //----------------------------------------------
        // Write solution
        //string   outname(strtok(argv[n],"."));
        //outname += ".vtu";
        ofstream outfile(fname.c_str());
        Exp[0]->WriteVtkHeader(outfile);

        // For each field write out field data for each expansion.
        for(i = 0; i < Exp[0]->GetNumElmts(); ++i)
        {
            Exp[0]->WriteVtkPieceHeader(outfile,i);
            // For this expansion, write out each field.
            for(j = 0; j < Exp.num_elements(); ++j)
            {
                Exp[j]->WriteVtkPieceData(outfile,i, fielddef[0]->m_fields[j]);
            }
            Exp[0]->WriteVtkPieceFooter(outfile,i);
        }
        Exp[0]->WriteVtkFooter(outfile);
        cout << "Written file: " << fname << endl;
        //----------------------------------------------
    }
    return 0;
}

