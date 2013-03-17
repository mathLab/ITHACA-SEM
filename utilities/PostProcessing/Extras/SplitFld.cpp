#include <cstdio>
#include <cstdlib>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>


using namespace Nektar;

int main(int argc, char *argv[])
{

    void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
        SpatialDomains::BoundaryConditionsSharedPtr& boundaryConditions,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables);
         Array<OneD, int> GetReflectionIndex(Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,
         int Ireg);
         
    if(argc != 3)
    {
        fprintf(stderr,"Usage: SplitFld  meshfile fieldfile\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);


    //----------------------------------------------

    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //----------------------------------------------

    // Also read and store the boundary conditions
    SpatialDomains::BoundaryConditionsSharedPtr boundaryConditions;        
    boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                                        ::AllocateSharedPtr(vSession, graphShPt);
    SpatialDomains::BoundaryConditions bcs(vSession, graphShPt);                                        
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-1]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    // Define Expansion    
    int nfields; 
    nfields = fielddef[0]->m_fields.size(); 
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp; 
    Exp = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);    

    std::string solvtype = vSession->GetSolverInfo("SOLVERTYPE");
    if(solvtype == "CoupledLinearisedNS" && vSession->DefinesSolverInfo("HOMOGENEOUS") )
    {
      	    
         SetFields(graphShPt,boundaryConditions,vSession,Exp,nfields-1);
 	 //decomment
         //nfields = nfields-1;
//start
         int lastfield = nfields-1;
         cout<<"Set pressure: "<<lastfield<<endl;           
         int nplanes = fielddef[0]->m_numModes[2];
         const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::ePolyEvenlySpaced);
         const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
         NekDouble lz = fielddef[0]->m_homogeneousLengths[0];
         MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
         Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,false,false,graphShPt,fielddef[0]->m_fields[0]);
         Exp[lastfield] = Exp3DH1;        
//end       
    }
    else
    {

    	 SetFields(graphShPt,boundaryConditions,vSession,Exp,nfields);
    }
    //----------------------------------------------    
	
    //----------------------------------------------
    // Copy data from field file
    for(int j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddef.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                        fielddef [i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans_IterPerExp(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
 
    //----------------------------------------------


    // Write solution to file with additional computed fields
    string   fldfilename(argv[2]);
    string   out = fldfilename.substr(0, fldfilename.find_last_of("."));
    string   endfile("split.fld");


    //Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(1);

    string outfile; 
    string var;


    Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(Exp.num_elements());
    //NB in case of homo fields you CANNOT use the BwdTrans 
    //because the Im comp is set to 0
    for(int i = 0; i < Exp.num_elements(); ++i)
    {
        fieldcoeffs[i] = Exp[i]->UpdateCoeffs();
    }

        
     // copy Data into FieldData and set variable
    /*
      for(int g=0; g<Exp[0]->GetPlane(1)->GetNcoeffs(); g++)
      {
      cout<<"g="<<g<<"  coeff f0="<<Exp[lastfield]->GetPlane(0)->GetCoeff(g)<<" f1="<<Exp[lastfield]->GetPlane(1)->GetCoeff(g)<<endl;
      }  
    */
    
     for(int j =0; j<nfields; j++)
     {
          outfile = out;
          std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
              = Exp[j]->GetFieldDefinitions(); 
          std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
          //fieldcoeffs[j] = fields[j]->UpdateCoeffs();             	    
          for(int i = 0; i < FieldDef.size(); ++i)
          {
               var = fielddef[i]->m_fields[j];
               // Could do a search here to find correct variable
               FieldDef[i]->m_fields.push_back(var);  	     
               Exp[j]->AppendFieldData(FieldDef[i], FieldData[i],fieldcoeffs[j]);

          }
          outfile += "_"+var+"_"+endfile;  
          LibUtilities::Write(outfile,FieldDef,FieldData); 

      }

      //-----------------------------------------------

    return 0;
}




	// Define Expansion       		
        void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
	    SpatialDomains::BoundaryConditionsSharedPtr &boundaryConditions,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables)
	{		
		// Setting parameteres for homogenous problems
		NekDouble LhomX;           ///< physical length in X direction (if homogeneous) 
		NekDouble LhomY;           ///< physical length in Y direction (if homogeneous)
		NekDouble LhomZ;           ///< physical length in Z direction (if homogeneous)
		
		bool DeclareCoeffPhysArrays = true;		
		int npointsX;              ///< number of points in X direction (if homogeneous)
		int npointsY;              ///< number of points in Y direction (if homogeneous)
                int npointsZ;              ///< number of points in Z direction (if homogeneous)	
		int HomoDirec       = 0;
		bool useFFT = false;	
		bool deal = false;        
		///Parameter for homogeneous expansions		
		enum HomogeneousType
		{
			eHomogeneous1D,
			eHomogeneous2D,
			eHomogeneous3D,
			eNotHomogeneous
		};

		enum HomogeneousType HomogeneousType = eNotHomogeneous;

		if(session->DefinesSolverInfo("HOMOGENEOUS"))
		{ 			
			std::string HomoStr = session->GetSolverInfo("HOMOGENEOUS");
			//m_spacedim          = 3;
			
			if((HomoStr == "HOMOGENEOUS1D")||(HomoStr == "Homogeneous1D")||
			   (HomoStr == "1D")||(HomoStr == "Homo1D"))
			{
				HomogeneousType = eHomogeneous1D;
				npointsZ        = session->GetParameter("HomModesZ");
				LhomZ           = session->GetParameter("LZ");
				HomoDirec       = 1;				
			}
			
			if((HomoStr == "HOMOGENEOUS2D")||(HomoStr == "Homogeneous2D")||
			   (HomoStr == "2D")||(HomoStr == "Homo2D"))
			{
				HomogeneousType = eHomogeneous2D;
				npointsY        = session->GetParameter("HomModesY");
				LhomY           = session->GetParameter("LY");
				npointsZ        = session->GetParameter("HomModesZ");
				LhomZ           = session->GetParameter("LZ");
				HomoDirec       = 2;
			}
			
			if((HomoStr == "HOMOGENEOUS3D")||(HomoStr == "Homogeneous3D")||
			   (HomoStr == "3D")||(HomoStr == "Homo3D"))
			{
				HomogeneousType = eHomogeneous3D;
				npointsX        = session->GetParameter("HomModesX");
				LhomX           = session->GetParameter("LX");
				npointsY        = session->GetParameter("HomModesY");
				LhomY           = session->GetParameter("LY");
				npointsZ        = session->GetParameter("HomModesZ");
				LhomZ           = session->GetParameter("LZ");
				HomoDirec       = 3;
			}
			
			if(session->DefinesSolverInfo("USEFFT"))
			{
				useFFT = true;
			}
		}		
			
	    int i;		
	    int expdim   = mesh->GetMeshDimension();
	    //Exp= Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);    
  	    // I can always have 3 variables in a 2D mesh (oech vel component i a function which can depend on 1-3 var)
        // Continuous Galerkin projection

            switch(expdim)
            {
                case 1:
                {
                    if(HomogeneousType == eHomogeneous2D)
                    {
                        const LibUtilities::PointsKey PkeyY(npointsY,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyY(LibUtilities::eFourier,npointsY,PkeyY);
                        const LibUtilities::PointsKey PkeyZ(npointsZ,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,npointsZ,PkeyZ);

                        for(i = 0 ; i < nvariables; i++)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous2D>
                                ::AllocateSharedPtr(session,BkeyY,BkeyZ,LhomY,LhomZ,useFFT,deal,mesh,session->GetVariable(i));
                        }
                    }
                    else
                    {
                        for(i = 0 ; i < nvariables; i++)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(session,mesh,session->GetVariable(i));
                        }
                    }

                    break;
                }
            case 2:
                {   
                    if(HomogeneousType == eHomogeneous1D)
                    {
                        const LibUtilities::PointsKey PkeyZ(npointsZ,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,npointsZ,PkeyZ);
                        for(i = 0 ; i < nvariables; i++)
                        {                        	
                            Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
                                ::AllocateSharedPtr(session,BkeyZ,LhomZ,useFFT,deal,mesh,session->GetVariable(i));                                    
                        }
                    }
                    else
                    {                   	    
                        i = 0;
                        MultiRegions::ContField2DSharedPtr firstfield;
                        firstfield = MemoryManager<MultiRegions::ContField2D>
                            ::AllocateSharedPtr(session,mesh,session->GetVariable(i),DeclareCoeffPhysArrays);

                        Exp[0] = firstfield;
                        for(i = 1 ; i < nvariables; i++)
                        {                        	
                            Exp[i] = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(*firstfield,mesh,session->GetVariable(i),DeclareCoeffPhysArrays);
                        }
                    }

                    break;
                }
                case 3:
                    {
                        if(HomogeneousType == eHomogeneous3D)
                        {
                            ASSERTL0(false,"3D fully periodic problems not implemented yet");
                        }
                        else
                        {
                            i = 0;
                            MultiRegions::ContField3DSharedPtr firstfield =
                                MemoryManager<MultiRegions::ContField3D>
                                ::AllocateSharedPtr(session,mesh,session->GetVariable(i));

                            Exp[0] = firstfield;
                            for(i = 1 ; i < nvariables; i++)
                            {
                                Exp[i] = MemoryManager<MultiRegions::ContField3D>
                                    ::AllocateSharedPtr(*firstfield,mesh,session->GetVariable(i));
                            }
                        }
                        break;
                    }
            default:
                ASSERTL0(false,"Expansion dimension not recognised");
                break;
            }   
        }   	   

