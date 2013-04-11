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

#include <sys/stat.h>

int fexist( const char *filename ) {
  struct stat buffer ;
  if ( stat( filename, &buffer ) ) return 0 ;
  return 1 ;
}

int main(int argc, char *argv[])
{
    unsigned int     i, j;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;

    if(argc < 3)
    {
        fprintf(stderr,"Usage: %s meshfile fieldfile\n",argv[0]);
        exit(1);
    }

    bool Extrude2DWithHomogeneous = false;
	
	bool SingleModePlot=false;
	bool HalfModePlot=false;

	
    int nExtraPoints, nExtraPlanes;
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);


    vSession->LoadParameter("OutputExtraPoints",nExtraPoints,0);
    vSession->LoadParameter("OutputExtraPlanes",nExtraPlanes,0);

    vSession->MatchSolverInfo("Extrude2DWithHomogeneous","True",Extrude2DWithHomogeneous,false);
    vSession->MatchSolverInfo("ModeType","SingleMode",SingleModePlot,false);
    vSession->MatchSolverInfo("ModeType","HalfMode",HalfModePlot,false);


    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);
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
            else if (ending == ".gz")
            {
                fname = fname.substr(0,fdot);
                fdot = fname.find_last_of('.');
                ASSERTL0(fdot != std::string::npos,
                         "Error: expected file extension before .gz.");
                ending = fname.substr(fdot);
                ASSERTL0(ending == ".xml",
                         "Compressed non-xml files are not supported.");
                continue;
            }
            else if (ending == ".xml")
            {
                continue;
            }
        }
#ifdef TECPLOT
        fname = fname + ".dat";
#else
        fname = fname + ".pos";
#endif

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
        vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
        vector<vector<NekDouble> > fielddata;
        LibUtilities::Import(fieldfile,fielddef,fielddata);
        //----------------------------------------------

        if(Extrude2DWithHomogeneous) // Set up Homogeneous information
        {
            NekDouble length;
            vSession->LoadParameter("LZ",length,1);
            fielddef[0]->m_numHomogeneousDir = 1;
            fielddef[0]->m_numModes.push_back(2); // Have to set this to 2 as default
            fielddef[0]->m_homogeneousZIDs.push_back(0);
            fielddef[0]->m_homogeneousLengths.push_back(length);
            fielddef[0]->m_basis.push_back(LibUtilities::eFourier);
        }

        if(SingleModePlot) // Set Up printing of perturbation
        {
            fielddef[0]->m_numModes.push_back(4); // Have to set this to 4 as default
            fielddef[0]->m_basis.push_back(LibUtilities::eFourier); //Initialisation of a standard Fourier Expansion
        }

        if(HalfModePlot)
        {
            fielddef[0]->m_numModes.push_back(4); // Have to set this to 4 as default
            fielddef[0]->m_basis.push_back(LibUtilities::eFourier); //Initialisation of a standard Fourier Expansion
            fielddef[0]->m_basis.push_back(LibUtilities::eFourierHalfModeIm);//Initialisation of a HalfModeFourierIm Expansion

        }

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
        bool useFFT = false;
        bool dealiasing = false;
        //----------------------------------------------


        //----------------------------------------------
        // Define Expansion
        int expdim   = graphShPt->GetMeshDimension();
        int nfields = fielddef[0]->m_fields.size();
        Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields);

        //auxiliary expansion for plotting perturbations
        Array<OneD, MultiRegions::ExpListSharedPtr> Exp1(nfields);

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
                    Exp1D = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vSession,graphShPt);
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
                ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

                if(fielddef[0]->m_numHomogeneousDir == 1)
                {
                    MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
                    MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1_aux;
                    MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1_Im;

                    if(SingleModePlot)
                    {
                        int nplanes = fielddef[0]->m_numModes[2];

                        const LibUtilities::PointsKey Pkey(nplanes+nExtraPlanes,LibUtilities::eFourierSingleModeSpaced);
                        //for plotting perturbations (4planes)
                        const LibUtilities::PointsKey Pkey1(nplanes+nExtraPlanes+2+1,LibUtilities::ePolyEvenlySpaced);

                        //SingleMode Basis
                        const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                        //Fourier expansion
                        const LibUtilities::BasisKey  Bkey1(fielddef[0]->m_basis[3],nplanes+2,Pkey1);

                        NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                       //Fourier SingleMode Expansion with two points
                        Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,dealiasing,graphShPt,fielddef[0]->m_fields[0]);

                        //Fourier 4 modes expansion
                        Exp3DH1_aux = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey1,lz,useFFT,dealiasing,graphShPt,fielddef[0]->m_fields[0]);
                        Exp1[0]= Exp3DH1_aux;


                        //Define Homogeneous standard 4 plane Fourier Expansion
                        for(i = 1; i < nfields; ++i)
                        {
                            Exp1[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                            ::AllocateSharedPtr(*Exp3DH1_aux);
                        }


                    }
                    else if(HalfModePlot)
                    {
                        int nplanes = fielddef[0]->m_numModes[2];

                        const LibUtilities::PointsKey Pkey(nplanes+nExtraPlanes,LibUtilities::eFourierSingleModeSpaced);
                        //for plotting perturbations (4planes)
                        const LibUtilities::PointsKey Pkey1(nplanes+nExtraPlanes+3+1,LibUtilities::ePolyEvenlySpaced);

                        //FourierHalfModeRe Basis
                        const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                        //Fourier expansion
                        const LibUtilities::BasisKey  Bkey1(fielddef[0]->m_basis[3],nplanes+3,Pkey1);
                        //FourierHalfModeIm Expansion
                        const LibUtilities::BasisKey  Bkey2(fielddef[0]->m_basis[4],nplanes,Pkey);


                        NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                        //FourierHalfModeRe Expansion
                        Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,dealiasing,graphShPt,fielddef[0]->m_fields[0]);
                        //FourierHalfModeIm Expansion
                        Exp3DH1_Im = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey2,lz,useFFT,dealiasing,graphShPt,fielddef[0]->m_fields[2]);

                        //Fourier 4 modes expansion
                        Exp3DH1_aux = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey1,lz,useFFT,dealiasing,graphShPt,fielddef[0]->m_fields[0]);
                        Exp1[0]= Exp3DH1_aux;


                        //Define Homogeneous standard 4 plane Fourier Expansion
                        for(i = 1; i < nfields; ++i)
                        {
                            Exp1[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                            ::AllocateSharedPtr(*Exp3DH1_aux);
                        }


                    }
                    else
                    {
                        //int nplanes = fielddef[0]->m_numModes[2];
                        int nplanes;
                        vSession->LoadParameter("HomModesZ",nplanes,fielddef[0]->m_numModes[2]);

                        // choose points to be at evenly spaced points at
                        // nplanes + 1 points
                        const LibUtilities::PointsKey Pkey(nplanes+nExtraPlanes+1,LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                        NekDouble lz = fielddef[0]->m_homogeneousLengths[0];


                        Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,dealiasing,graphShPt,fielddef[0]->m_fields[0]);
                    }


                    //it is a FourierSingleMode or HalfMode in case
                    if(HalfModePlot)
                    {
                        Exp[0] = Exp3DH1;
                        for(i = 1; i < nfields; ++i)
                        {
                            //w must have imaginary basis
                            if(i==2)
                            {
                                Exp[2] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                                ::AllocateSharedPtr(*Exp3DH1_Im);
                            }
                            else
                            {
                                Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                                        ::AllocateSharedPtr(*Exp3DH1);
                            }
                        }


                    }
                    else
                    {
                        Exp[0] = Exp3DH1;
                        for(i = 1; i < nfields; ++i)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                            ::AllocateSharedPtr(*Exp3DH1);

                        }
                    }
                }
                else
                {
                    MultiRegions::ExpList2DSharedPtr Exp2D;
                    Exp2D = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vSession,graphShPt,true,fielddef[0]->m_fields[0]);
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

       if(Extrude2DWithHomogeneous)
       {
            // Need to set this back to 1 to read 2D field
            // Perhaps could set up Extra parameters?
            fielddef[0]->m_numModes[2] = 1;
        }

        //----------------------------------------------
        // Copy data to file
        for(j = 0; j < nfields; ++j)
        {
            for(int i = 0; i < fielddata.size(); ++i)
            {

                Exp[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],
                                            fielddef[i]->m_fields[j],
                                            Exp[j]->UpdateCoeffs());
            }

            if(SingleModePlot)
            {
                //it is on two planes for single mode
                int dim=Exp[j]->GetNcoeffs();
                //copy the single mode on the 4planes expansion
                Vmath::Vcopy(dim,&Exp[j]->GetCoeffs()[0],1,&Exp1[j]->UpdateCoeffs()[dim],1);

                Exp1[j]->BwdTrans(Exp1[j]->GetCoeffs(),Exp1[j]->UpdatePhys());

            }
            else if(HalfModePlot)
            {
                //it is one planes for single mode
                int dim=Exp[j]->GetNcoeffs();
                //copy the  mode on the 4planes expansion
                if(j==2)
                {
                    //copy on the 4th plane
                    //Vmath::Vcopy(dim,&Exp[j]->GetCoeffs()[0],1,&Exp1[j]->UpdateCoeffs()[3*dim],1);
                    Vmath::Vcopy(dim,&Exp[j]->GetCoeffs()[0],1,&Exp1[j]->UpdateCoeffs()[2*dim],1);

                }
                else
                {
                    Vmath::Vcopy(dim,&Exp[j]->GetCoeffs()[0],1,&Exp1[j]->UpdateCoeffs()[2*dim],1);
                }

                Exp1[j]->BwdTrans(Exp1[j]->GetCoeffs(),Exp1[j]->UpdatePhys());

            }
            else
            {
                Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
            }
        }
        //----------------------------------------------

        //----------------------------------------------
        // Write solution  depending on #define
    #ifdef TECPLOT

        if(SingleModePlot || HalfModePlot)
        {
            std::string var = "";

            for(int j = 0; j < Exp1.num_elements(); ++j)
            {
                var = var + ", " + fielddef[0]->m_fields[j];
            }

            ofstream outfile(fname.c_str());

            Exp1[0]->WriteTecplotHeader(outfile,var);
            for(int i = 0; i < Exp1[0]->GetNumElmts(); ++i)
            {
                Exp1[0]->WriteTecplotZone(outfile,i);
                for(int j = 0; j < Exp1.num_elements(); ++j)
                {
                    Exp1[j]->WriteTecplotField(outfile,i);
                }
            }
        }
        else{

            std::string var = "";

            for(int j = 0; j < Exp.num_elements(); ++j)
            {
                var = var + ", " + fielddef[0]->m_fields[j];
            }

            ofstream outfile(fname.c_str());

            Exp[0]->WriteTecplotHeader(outfile,var);
            for(int i = 0; i < Exp[0]->GetNumElmts(); ++i)
            {
                Exp[0]->WriteTecplotZone(outfile,i);
                for(int j = 0; j < Exp.num_elements(); ++j)
                {
                    Exp[j]->WriteTecplotField(outfile,i);
                }
            }

        }
    #else
        for(i = 0; i < nfields; ++i)
        {
            ofstream outstrm(fname.c_str());
    
            Exp[i]->WriteToFile(outstrm,eGmsh);
            outstrm.close();
        }
    #endif
        //----------------------------------------------
    }
    return 0;
}

