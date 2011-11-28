#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList0D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList1DHomogeneous2D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <LocalRegions/MatrixKey.h>


using namespace Nektar;

int main(int argc, char *argv[])
{
    void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
        SpatialDomains::BoundaryConditionsSharedPtr& boundaryConditions,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables);
    void Extractlayerdata(Array<OneD, int> Iregions, int coordim, SpatialDomains::MeshGraphSharedPtr &mesh,   
    	    SpatialDomains::BoundaryConditions &bcs,
    	    Array<OneD,MultiRegions::ExpListSharedPtr> &infields,  
    	    MultiRegions::ContField1DSharedPtr &outfieldx,
       	    MultiRegions::ContField1DSharedPtr &outfieldy,
       	    MultiRegions::ExpListSharedPtr &streak);
    void Manipulate(Array<OneD, int> Iregions, int coordim, SpatialDomains::MeshGraphSharedPtr &mesh,   
    	    SpatialDomains::BoundaryConditions &bcs,
    	    Array<OneD,MultiRegions::ExpList1DSharedPtr> &infields,  
    	    Array<OneD,MultiRegions::ExpList1DSharedPtr> &outfieldx,
       	       Array<OneD,MultiRegions::ExpList1DSharedPtr> &outfieldy,
       	       MultiRegions::ExpListSharedPtr &streak);
    void WriteBcs(string variable, int region, string fieldfile, SpatialDomains::MeshGraphSharedPtr &mesh, 
    	    MultiRegions::ContField1DSharedPtr &outregionfield);
    void WriteFld(string outfile, SpatialDomains::MeshGraphSharedPtr &mesh, Array<OneD,
    	                       MultiRegions::ExpListSharedPtr> &fields );
    
    int i,j;
    if(argc != 4)
    {
        fprintf(stderr,"Usage: ./FldCalcBCs  meshfile fieldfile  streakfile\n");
        exit(1);
    }
 
    //----------------------------------------------
    string meshfile(argv[argc-3]);
  
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);
    //----------------------------------------------
  
    // Also read and store the boundary conditions
    SpatialDomains::MeshGraph *meshptr = graphShPt.get();
    SpatialDomains::BoundaryConditionsSharedPtr boundaryConditions;        
    boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                                        ::AllocateSharedPtr(vSession, graphShPt);
    SpatialDomains::BoundaryConditions bcs(vSession, graphShPt);                                        
    //----------------------------------------------
      
    // Import field file.
    string fieldfile(argv[argc-2]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    // Define Expansion    
    int nfields; 
    nfields = fielddef[0]->m_fields.size(); 
    Array<OneD, MultiRegions::ExpListSharedPtr> fields; 
    fields= Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);    
cout<<nfields<<endl;    

    std::string solvtype = vSession->GetSolverInfo("SOLVERTYPE");
    if(solvtype == "CoupledLinearisedNS")
    {
      	    
         SetFields(graphShPt,boundaryConditions,vSession,fields,nfields-1);
         int lastfield = nfields-1;
         cout<<"Set pressure: "<<lastfield<<endl;           
         int nplanes = fielddef[0]->m_numModes[2];
         const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::ePolyEvenlySpaced);
         const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
         NekDouble lz = fielddef[0]->m_homogeneousLengths[0];
         MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
         Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,false,graphShPt,fielddef[0]->m_fields[0]);
         fields[lastfield] = Exp3DH1;        
       
    }
    else
    {
    	 SetFields(graphShPt,boundaryConditions,vSession,fields,nfields);
    }
    //----------------------------------------------       
 
    // Copy data from file:fill fields with the fielddata
    for(j = 0; j < nfields; ++j)
    {  	       	    
        for(int i = 0; i < fielddata.size(); ++i)
        {
            fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j]);
        }             
        fields[j]->BwdTrans_IterPerExp(fields[j]->GetCoeffs(),fields[j]->UpdatePhys());               
    }
    //----------------------------------------------
                            
    //--------------------------------------------------------------    



    // determine the I regions
    //hypothesis: the number of I regions is the same for all the variables
    //hypothesis: all the I regions have the same nq points
    int nIregions, lastIregion; 
    const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConditions  = fields[0]->GetBndConditions();      
    Array<OneD, int> Iregions =Array<OneD, int>(bndConditions.num_elements(),-1);    
    //3 internal layers:
    Array<OneD, int> Ilayers =Array<OneD, int>(3,-1);     
    nIregions=0;
    int nbnd= bndConditions.num_elements();    
    for(int r=0; r<nbnd; r++)
    {
    	  if(bndConditions[r]->GetUserDefined().GetEquation()=="CalcBC")
    	  {
    	  	  lastIregion=r;
    	  	  Iregions[r]=r;
    	  	  Ilayers[nIregions]=r;
    	  	  nIregions++;
    	  }    	  
    } 
    ASSERTL0(nIregions>0,"there is any boundary region with the tag USERDEFINEDTYPE=""CalcBC"" specified");



    // import the streak 

cout<<"streak"<<endl; 
    //import the streak field
    //add y to obtain the streak w=w'+y and write fld file           
    MultiRegions::ExpListSharedPtr streak; 
    streak = MemoryManager<MultiRegions::ContField2D>
          ::AllocateSharedPtr(vSession, graphShPt, "w",true);
    //streak = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr
    //          (vSession, graphShPt, true, "w");
    
    string streakfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> streakdef;
    vector<vector<NekDouble> > streakdata;
    graphShPt->Import(streakfile,streakdef,streakdata); 
    //attention: streakdata.size==2 because the session file is 3DHomo1D but in
    //reality the streak is real quantity
    

    // Copy data from file:fill streak with the streakdata
    	       	    
    for(int i = 0; i < streakdata.size(); ++i)
    {        	
            streak->ExtractDataToCoeffs(streakdef[i],streakdata[i],streakdef[i]->m_fields[0]);
    }             
    streak->BwdTrans(streak->GetCoeffs(),streak->UpdatePhys());
    int totpoints = fields[0]->GetPlane(0)->GetTotPoints();

    
    
//    Vmath::Vcopy(totpoints,&(streak->GetPhys()[0]),1,&(fields[1]->GetPlane(0)->UpdatePhys()[0]),1);

/*
    for(int q=0; q<10; q++)
    {
cout<<"streak="<<streak->GetPhys()[q]<<"  field="<<fields[1]->GetPlane(0)->GetPhys()[q]<<endl;
cout<<"   1Dstreak="<<Istreak[Ilayers[1]]->GetPhys()[q]<<endl;    	    
    	    
    }
*/
                                         
/*   
    //SetFields(graphShPt,boundaryConditions,vSession,streak,nvar);
    int nq=streak[0]->GetTotPoints();  
cout<<"nq 3D="<<nq<<"  nq 2D="<<nq2D<<"fields[0] nq2D+1="<<(fields[0]->GetPhys())[nq2D+1]<<endl;    
    Array<OneD, NekDouble> x(nq,0.0);
    Array<OneD, NekDouble> y(nq,0.0);
    //fields[0]->GetCoords(x,y); 
  
    Array<OneD, NekDouble> z(nq,0.0);    
    // z point are always 0 in the mesh!!! something related to Lz number of planes??    
    streak[0]->GetCoords(x,y,z);   
cout<<z[2*nq2D+1]<<endl;
    int nplanes = vSession->GetParameter("HomModesZ");
    for(int j=0; j<nplanes; j++)
    {
        for(int i=0; i<nq2D; i++)
        {
             (streak[0]->UpdatePhys())[i+j*nq2D]=0;
             (streak[1]->UpdatePhys())[i+j*nq2D]=0;
             //w=v+y    
             (streak[2]->UpdatePhys())[i+j*nq2D]=(fields[1]->GetPhys())[i]+y[i]; 
cout<<"x="<<x[i+j*nq2D]<<"  y="<<y[i+j*nq2D]<<"  z="<<z[i+j*nq2D]<<"  v="<<(fields[1]->GetPhys())[i]
<<"  streak="<<(streak[2]->GetPhys())[i+j*nq2D]<<endl;	  
cout<<"plane="<<j<<endl;
	}
    }	
*/     
/*
    //write fld file for the streak
    string   file = fieldfile.substr(0,fieldfile.find_last_of("."));
    file= file+"_streak.fld";
    
    Array<OneD, MultiRegions::ExpListSharedPtr> dump;
    dump = Array<OneD, MultiRegions::ExpListSharedPtr>(2);
    dump[0] =streak;
    dump[1] =streak;
    
 
   WriteFld(file,graphShPt,dump);
*/   	
    //---------------------------------------------------------------   

    //set 1D output fields:
    //initialise fields
    const SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
/*
    MultiRegions::ExpList1DSharedPtr outfieldx;    		
    MultiRegions::ExpList1DSharedPtr outfieldy;
    outfieldx = MemoryManager<MultiRegions::ExpList1D>
                ::AllocateSharedPtr(*(bregions[lastIregion]), graphShPt, true);    
    outfieldy = MemoryManager<MultiRegions::ExpList1D>
    		::AllocateSharedPtr(*(bregions[lastIregion]), graphShPt, true);
*/

    MultiRegions::ExpList1DSharedPtr outfieldxtmp;    		
    MultiRegions::ExpList1DSharedPtr outfieldytmp;
    MultiRegions::ContField1DSharedPtr outfieldx;    		
    MultiRegions::ContField1DSharedPtr outfieldy;
    //initialie fields
    //const SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
    //declarecoeffsphysarray?????????!!!! true?!      

    outfieldxtmp = MemoryManager<MultiRegions::ExpList1D>
                ::AllocateSharedPtr(*(bregions[lastIregion]), graphShPt, true);    
    outfieldytmp = MemoryManager<MultiRegions::ExpList1D>
    		::AllocateSharedPtr(*(bregions[lastIregion]), graphShPt, true);

    outfieldx = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(vSession, *outfieldxtmp);

    outfieldy = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(vSession, *outfieldytmp);    
   int a = outfieldy->GetNcoeffs();
   int a1 = outfieldytmp->GetNcoeffs();
cout<<"kokok  a="<<a<<"  a1="<<a1<<endl;  
//YOU NEED TO DEFINE THE OUTFIELDS AS ContField2D to use them in the 2D NS!!!!!!!!!!!
//(not 3DHomogeneous1D!!!!!!!

/*  
    //define the outfields:
    int nplanes = fielddef[0]->m_numModes[1];
    int nzlines = fielddef[0]->m_numModes[2];    
    //ATTENTO Polyevenlyspaced NON Fourier!!!
    const LibUtilities::PointsKey PkeyZ(nzlines,LibUtilities::ePolyEvenlySpaced);
    const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);
    NekDouble lz = fielddef[0]->m_homogeneousLengths[1];
    for(int r=0; r<nlayers; r++)
    {
    	    
         bndfieldx[r]= MemoryManager<MultiRegions::ExpList2DHomogeneous1D>
                             ::AllocateSharedPtr(vSession,BkeyZ,lz,false,graphShPt);
cout<<"OOOK"<<endl;                               
         bndfieldy[r]= MemoryManager<MultiRegions::ExpList2DHomogeneous1D>
                             ::AllocateSharedPtr(vSession,BkeyZ,lz,false,graphShPt);
    }
    //--------------------------------------------------------
*/    
    
    //manipulate data
   
    //for 2 variables(u,v) only:
    int coordim = graphShPt->GetMeshDimension();              	   
    Extractlayerdata(Ilayers,coordim, graphShPt, bcs, fields, outfieldx,outfieldy,streak);       	       

    //--------------------------------------------------------------------------------------


  

    //write bcs files: one for each I region and each variable
//    for(int s=0; s<nbnd; s++)
//    {    	 
//          if(Iregions[s]!=-1)
//      	  {
      	        string var="u";
      	        WriteBcs(var,lastIregion, fieldfile,graphShPt,outfieldx);
      	        var="v";
      	        WriteBcs(var,lastIregion, fieldfile,graphShPt,outfieldy);      	      	      
//      	  }
//     }
    
}
				
	// Define Expansion       		
	void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
	    SpatialDomains::BoundaryConditionsSharedPtr &boundaryConditions,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables)
	{		
		// Setting parameteres for homogenous problems
        	MultiRegions::GlobalSysSolnType solnType;
		NekDouble LhomX;           ///< physical length in X direction (if homogeneous) 
		NekDouble LhomY;           ///< physical length in Y direction (if homogeneous)
		NekDouble LhomZ;           ///< physical length in Z direction (if homogeneous)
		
		bool DeclareCoeffPhysArrays = true;		
		int npointsX;              ///< number of points in X direction (if homogeneous)
		int npointsY;              ///< number of points in Y direction (if homogeneous)
		int npointsZ;              ///< number of points in Z direction (if homogeneous)		
		int HomoDirec       = 0;
		bool useFFT = false;	        
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
                                ::AllocateSharedPtr(session,BkeyY,BkeyZ,LhomY,LhomZ,useFFT,mesh,session->GetVariable(i));
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
                                ::AllocateSharedPtr(session,BkeyZ,LhomZ,useFFT,mesh,session->GetVariable(i));                                
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

       void  Extractlayerdata(Array<OneD, int> Iregions,int coordim, SpatialDomains::MeshGraphSharedPtr &mesh,
       	       SpatialDomains::BoundaryConditions &bcs,
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &infields, 
       	       MultiRegions::ContField1DSharedPtr &outfieldx,
       	       MultiRegions::ContField1DSharedPtr &outfieldy,
       	       MultiRegions::ExpListSharedPtr &streak)
        {
            //3 I regions are expected: (the layer is the last region)
            ASSERTL0(Iregions.num_elements()==3, "something wrong with the number of I layers");
            //int Iup =Iregions[0];
            int Ireg =Iregions[2];
cout<<"layer region="<<Ireg<<endl;          
            //take the pressure from infields
            MultiRegions::ExpListSharedPtr pressure =infields[3];
            Array<OneD, MultiRegions::ExpListSharedPtr> Iexp =infields[0]->GetBndCondExpansions();
            //MultiRegions::ExpListSharedPtr Ilayerup  = Iexp[Iup];
            MultiRegions::ExpListSharedPtr Ilayer = Iexp[Ireg];
            //take the streak data along the layers (getbondcondexpansions does not work!!!)          
            //Array<OneD, MultiRegions::ExpListSharedPtr> Istreak = streak->GetBndCondExpansions();
            //Array<OneD, MultiRegions::ExpListSharedPtr> Istreak =
              //                                     infields[1]->GetPlane(0)->GetBndCondExpansions();
            //MultiRegions::ExpListSharedPtr Istreakup = Istreak[Iup];
            //MultiRegions::ExpListSharedPtr Istreakdown = Istreak[Idown];

            
            
            
            
            int np = streak->GetTotPoints();
            Array<OneD, NekDouble> gradx(np);
            Array<OneD, NekDouble> grady(np); 

            
/******************************************************************************************/        
/*
//test gradx,grady
            Array<OneD,NekDouble> x(np);
            Array<OneD,NekDouble> y(np);  
            Array<OneD,NekDouble> z(np); 
            streak->GetCoords(x,y,z);
                        
            
            for(int a=0; a<np; a++)
            {
                streak->UpdatePhys()[a] = x[a]*( 1- x[a] );             	    
            }
*/
/*****************************************************************************************/
            
            int totcoeffs = streak->GetNcoeffs();
            Array<OneD, NekDouble> gradxcoeffs(totcoeffs);
            Array<OneD, NekDouble> gradycoeffs(totcoeffs);

            streak->PhysDeriv(streak->GetPhys(), gradx, grady);    
            streak->FwdTrans_IterPerExp(gradx, gradxcoeffs);
            streak->FwdTrans_IterPerExp(grady, gradycoeffs);
            Array<OneD, Array<OneD, NekDouble> > normals;
            normals = Array<OneD, Array<OneD, NekDouble> >(2);            
            Array<OneD, Array<OneD, NekDouble> > tangents;
	    tangents = Array<OneD, Array<OneD, NekDouble> >(2);		     

            int nq1D = Ilayer->GetPlane(0)->GetTotPoints();
           
            Array<OneD,NekDouble> x0d(nq1D);
            Array<OneD,NekDouble> x1d(nq1D);  
            Array<OneD,NekDouble> x2d(nq1D); 
            Ilayer->GetPlane(0)->GetCoords(x0d,x1d,x2d);
            Array<OneD, NekDouble> nx(nq1D);
            Array<OneD, NekDouble> ny(nq1D);
            Array<OneD, NekDouble> tx(nq1D);
            Array<OneD, NekDouble> ty(nq1D);            
            

            
            
cout<<"region="<<Iregions[0]<<" nq1D="<<nq1D<<endl;

            //initialise the curvature array:
            Array<OneD, NekDouble> curv=Array<OneD, NekDouble>(nq1D,0.0);
            

 


            const SpatialDomains::BoundaryRegionCollection &bregions
                = bcs.GetBoundaryRegions();
            SpatialDomains::Composite Icompreg;

            SpatialDomains::SegGeomSharedPtr segmentGeomreg;
            SpatialDomains::ElementEdgeVectorSharedPtr elementreg;
            StdRegions::EdgeOrientation orientreg;            
            SpatialDomains::BoundaryRegion::iterator regIt;
            Array<OneD, unsigned int> bmapreg;
            //fields to store the Re Im parts in physical space:  
       	    Array<OneD,NekDouble> Rephysreg (nq1D);
  	    Array<OneD,NekDouble> Imphysreg (nq1D);
  	    
  	    //Array<OneD, NekDouble> stphysup (nq1D);
  	    Array<OneD, NekDouble> stphysreg (nq1D);
  	    Array<OneD, NekDouble> stphysgradxreg (nq1D);
  	    Array<OneD, NekDouble> stphysgradyreg (nq1D);
  	    
  	    
  	    
  	    
  	    
            //nummodes*nedges_layerdown
            int Nregcoeffs;  

            for(
       	        regIt= bregions[Ireg]->begin();
       	        regIt !=bregions[Ireg]->end();
       	        ++regIt)
       	    {
                Icompreg=regIt->second;





		//assuming that ncoeffs is equal for all the edges (k=0)
		int nummodes = Ilayer->GetExp(0)->GetNcoeffs();
		Array<OneD, int> offsetup(Icompreg->size());
		offsetup[0]=0;
                int offsetregIExp;

            	//array coeffs:             	
                Nregcoeffs =   Icompreg->size()*nummodes;
            	Array<OneD, NekDouble> Recoeffsreg (Nregcoeffs);
            	Array<OneD, NekDouble> Imcoeffsreg (Nregcoeffs);

            	Array<OneD, NekDouble> stcoeffsreg (Nregcoeffs);

            	//grad coeffs
            	Array<OneD, NekDouble> stgradxcoeffsreg (Nregcoeffs);
            	Array<OneD, NekDouble> stgradycoeffsreg (Nregcoeffs);            
            
            
            
            
                        	
            	
		//define nq points for each edge
            	int nqedge =nq1D/Icompreg->size();
            	
		Array<OneD,NekDouble> x0edge(nqedge);
		Array<OneD,NekDouble> x1edge(nqedge);    
                Array<OneD, NekDouble> Rephysedgereg(nqedge);
                Array<OneD, NekDouble> Imphysedgereg(nqedge);
                

                for(int f=0; f<2; ++f)
                {
                     normals[f]  = Array<OneD, NekDouble>(nqedge); 
                     tangents[f] = Array<OneD, NekDouble>(nqedge);
                }                
        	
                for(int k=0; k< Icompreg->size(); k++)//loop over segments of each layer
                {
                    
                     if(
                             !(segmentGeomreg = boost::dynamic_pointer_cast<
                             	     SpatialDomains::SegGeom>((*Icompreg)[k]))
                     )
		     { ASSERTL0(false, "dynamic cast to a SegGeom failed");  }
		     
		     int EIDreg = segmentGeomreg->GetEid();	     
                     offsetregIExp = Ilayer->GetPlane(0)->GetCoeff_Offset(k);
cout<<"edge="<<EIDreg<<"  CHECK OFFSET="<<offsetregIExp<<endl;                                   
                     elementreg = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(mesh)
                                    ->GetElementsFromEdge(segmentGeomreg);
	     	     int elmtidreg = ((*elementreg)[0]->m_Element)->GetGlobalID();   
                     int dimelmtreg = ((*elementreg)[0]->m_Element)->GetNumEdges();
                    
                     orientreg =(boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>
                     	           ((*elementreg)[0]->m_Element))->GetEorient((*elementreg)[0]
                     	           	   ->m_EdgeIndx);
                     
                     
		     //	setup map between global and local ids
                     map<int, int>EdgeGIDreg;
		     int id,cnt,f;
                     for(cnt=f=0; f<dimelmtreg; f++)
		     {
			   id = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>
                     	           ((*elementreg)[0]->m_Element))->GetEid(f);
			   EdgeGIDreg[id]=cnt++;		     	     
		     }
cout<<" GIDreg="<<elmtidreg<<"  num edges in this element="<<dimelmtreg<<endl;
	             int localidedge = EdgeGIDreg.find(EIDreg)->second;
cout<<"locaidedge from m_EdgeInx ="<<localidedge<<endl;	             
	             localidedge = (*elementreg)[0]->m_EdgeIndx;
cout<<"check map local id of Eid="<<EIDreg<<"  is="<<localidedge<<endl;	             


		     //set bmaps -------------
		     Array<OneD, unsigned int > bmapedgereg;
		     Array<OneD, int > signreg;
		     pressure->GetExp(elmtidreg)->GetEdgeToElementMap(EdgeGIDreg.find(EIDreg)->second,
		     	     orientreg,bmapedgereg,signreg);		     		     

cout<<"SIGN down="<<signreg[0]<<endl;
		     //Extracting the coeffs...
		     Array<OneD, NekDouble> Recoeffs = pressure->GetPlane(0)->GetCoeffs();
		     Array<OneD, NekDouble> Imcoeffs = pressure->GetPlane(1)->GetCoeffs();  
		     Array<OneD, NekDouble> stcoeffs = streak->GetCoeffs();
                     //OFFSET EDGE REFERING TO ELEMENT EXPLIST MISSING		     
		     int Reoffsetreg = pressure->GetPlane(0)->GetCoeff_Offset(elmtidreg);
		     int Imoffsetreg = pressure->GetPlane(1)->GetCoeff_Offset(elmtidreg);
		     int stoffsetreg = streak->GetCoeff_Offset(elmtidreg);
		     
cout<<" Re offsetdown="<<Reoffsetreg<<"   Im offset="<<Imoffsetreg<<endl;		     
		     

/*********************************************************************************/

		     int ReNumoffset = Ilayer->GetPlane(0)->GetNcoeffs();
                     int ImNumoffset = Ilayer->GetPlane(1)->GetNcoeffs();	
//cout<<" OFFSET TEst Re="<<ReNumoffset<<"  Im="<<ImNumoffset<<endl;                     
		     int Reoffsetedgereg = Ilayer->GetPlane(0)->GetCoeff_Offset(k);
cout<<"k="<<k<<"  Offsetedge="<<Reoffsetedgereg<<endl;		     
		     
		     //hypothesis: 2 planes (0 Re, 1 Im)
		     int nplanecoeffs = pressure->GetExp(elmtidreg)->GetNcoeffs();
		     Array<OneD, NekDouble> Recoeffsedgereg(bmapedgereg.num_elements());	
		     Array<OneD, NekDouble> Imcoeffsedgereg(bmapedgereg.num_elements());
		     Array<OneD, NekDouble> stcoeffsedgereg(bmapedgereg.num_elements());
		     
		     
		     //grad edge coeffs
		     Array<OneD, NekDouble> stgradxcoeffsedge (bmapedgereg.num_elements());
		     Array<OneD, NekDouble> stgradycoeffsedge (bmapedgereg.num_elements());		     
cout<<" Recoeffsedgeup num elements="<<Recoeffsedgereg.num_elements()
<<"  n bmapedge="<<bmapedgereg.num_elements()<<endl;	



/*************************************************************************************/


            	     int nqed = Ilayer->GetPlane(0)->GetTotPoints();
            	     Array<OneD,NekDouble> x0ed(nqed);
            	     Array<OneD,NekDouble> x1ed(nqed);  
            	     Array<OneD,NekDouble> x2ed(nqed); 
            	     Ilayer->GetPlane(0)->GetExp(k)->GetCoords(x0ed,x1ed,x2ed);



		     for(int d=0; d<bmapedgereg.num_elements(); d++)
		     {	 
		     	   Recoeffsedgereg[d]=Recoeffs[Reoffsetreg+bmapedgereg[d]];		
                           Imcoeffsedgereg[d]=Imcoeffs[Imoffsetreg+bmapedgereg[d]];
                           
                           //streak edge coeffs
                           stcoeffsedgereg[d] = stcoeffs[stoffsetreg + bmapedgereg[d]];
                           stgradxcoeffsedge[d] = gradxcoeffs[stoffsetreg + bmapedgereg[d]];
                           stgradycoeffsedge[d] = gradycoeffs[stoffsetreg + bmapedgereg[d]];
                           
		     	   Recoeffsreg[offsetregIExp +d]  = Recoeffsedgereg[d];
		     	   Imcoeffsreg[offsetregIExp +d]  = Imcoeffsedgereg[d];
		     	   
		     	   //streak coeffs
		     	   stcoeffsreg[offsetregIExp +d] = stcoeffsedgereg[d];
		     	   stgradxcoeffsreg[offsetregIExp +d] = stgradxcoeffsedge[d];
		     	   stgradycoeffsreg[offsetregIExp +d] = stgradycoeffsedge[d];
		     	   
//cout<<"x="<<x0ed[d]<<"  y="<<x1ed[d]<<"  Re coeff="<<Recoeffsdown[offsetdownIExp +d]<<"  Im coeff="<<Imcoeffsdown[offsetdownIExp +d]<<endl; 		     	   
		     }

cout<<k<<"ok  dimcoeffarray="<<Recoeffsedgereg.num_elements()<<"  dimedgeup="<<Rephysedgereg.num_elements()<<endl;
cout<<"Nlayercoeffs="<<Recoeffsreg.num_elements()<<"   Nlayerphys="<<Rephysreg.num_elements()<<endl;
		     
		     
		     
cout<<"extract normal for edge="<<k<<endl;	
	     
cout<<"extract tangent for edge="<<k<<endl;
              	     //k is the number of the edge according to the composite list
		     LocalRegions::Expansion1DSharedPtr edgeexp = 
		     		boost::dynamic_pointer_cast<LocalRegions::Expansion1D>
		       (Ilayer->GetPlane(0)->GetExp(k) );
		     int localEid = edgeexp->GetLeftAdjacentElementEdge();
                     normals = edgeexp->
                            GetLeftAdjacentElementExp()->GetEdgeNormal(localEid);
		     LocalRegions::SegExpSharedPtr  bndSegExp = 
		          boost::dynamic_pointer_cast<LocalRegions::SegExp>(Ilayer->GetPlane(0)->GetExp(k));                             
	             tangents = (bndSegExp)->GetMetricInfo()->GetEdgeTangent();
                     for(int e=0; e< nqedge; e++)
                     {
cout<<" nx="<<normals[0][e]<<"    ny="<<normals[1][e]<<"   e="<<e<<endl; 
//cout<<"tannn tx="<<tangents[0][e]<<"   ty="<<tangents[1][e]<<endl;
                          nx[k*nqedge +e] = normals[0][e];
                          ny[k*nqedge +e] = normals[1][e];
                          tx[k*nqedge +e] = tangents[0][e];
                          ty[k*nqedge +e] = tangents[1][e];
                     }
		
 

                }
cout<<"comp closed"<<endl;         
                //bwd transform:
		Ilayer->GetPlane(0)->BwdTrans(Recoeffsreg, Rephysreg);
		Ilayer->GetPlane(1)->BwdTrans(Imcoeffsreg, Imphysreg);
		//streak phys values:
		Ilayer->GetPlane(0)->BwdTrans(stcoeffsreg, stphysreg);
		Ilayer->GetPlane(0)->BwdTrans_IterPerExp(stgradxcoeffsreg, stphysgradxreg);
		Ilayer->GetPlane(0)->BwdTrans_IterPerExp(stgradycoeffsreg, stphysgradyreg);
            }  
cout<<"Dim Rephysreg="<<Rephysreg.num_elements()<<"  nq1D="<<nq1D<<endl;            








             //TANGENTS WRONG WITH 3DHOMOGENEOUS 1D...	     
             //calculate the curvature:
             Array<OneD, NekDouble> dtx(nq1D);
             Array<OneD, NekDouble> dty(nq1D);
             Array<OneD, NekDouble> txcoeffs (Nregcoeffs);
             Array<OneD, NekDouble> tycoeffs (Nregcoeffs);
             Ilayer->GetPlane(0)->FwdTrans(tx, txcoeffs);
             Ilayer->GetPlane(0)->FwdTrans(ty, tycoeffs);
             Ilayer->GetPlane(0)->BwdTrans(txcoeffs, tx);
             Ilayer->GetPlane(0)->BwdTrans(tycoeffs, ty);
             Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, tx, dtx);
             Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, ty, dty);
cout<<"x     y      tx      ty      dtx        dty       curv"<<endl;             
             for(int t=0; t<nq1D; t++)
	     {
	            curv[t]=sqrt(dtx[t]*dtx[t] +dty[t]*dty[t]);     
cout<<setw(13)<<x0d[t]<<"     "<<setw(13)<<x1d[t]<<"      "<<setw(13)<<tx[t]<<setw(13)<<"     "<<ty[t]<<"     "
<<setw(13)<<dtx[t]<<"      "<<dty[t]<<"     "<<curv[t]<<endl;
             }		     
		    



           //normalise the pressure  norm*sqrt[ (\int Psquare)/Area]=1 =>
           // norm = sqrt[ Area /  (\int Pquare)]
           Array<OneD, NekDouble> P_square(np, 0.0);
           Array<OneD, NekDouble> P_tmp(np);
           Array<OneD, NekDouble> I (np, 1.0);             
           NekDouble Area = pressure->GetPlane(0)->PhysIntegral(I);
           //NekDouble IntPsquare = pressure->PhysIntegral(
        
           Vmath::Vmul(np, pressure->GetPlane(0)->GetPhys(),1,pressure->GetPlane(0)->GetPhys(),1,P_square,1);
           Vmath::Vmul(np, pressure->GetPlane(0)->GetPhys(),1,pressure->GetPlane(0)->GetPhys(),1,P_tmp,1);
           Vmath::Vadd(np, P_tmp,1,P_square,1, P_square,1);
           NekDouble IntPsquare = pressure->L2();
           NekDouble Int = pressure->PhysIntegral(P_square);
           NekDouble norm = sqrt(Area/Int);

cout<<" norm="<<norm<<"   IntPsquare="<<IntPsquare<<"   Area="<<Area<<"  Int="<<Int<<endl;







            //fill output fields:
            for(int g=0; g<nq1D; g++)
            {
            	(outfieldx->UpdatePhys())[g] = 
            	   Rephysreg[g];
            	(outfieldy->UpdatePhys())[g] =
            	   Imphysreg[g];	   
//cout<<"curvature x="<<x0d[g]<<"   k="<<curv[g]<<endl;		
            }

cout<<"Down: x"<<"   y"<<"  z"<<"    Re p"<<"    Im p"<<endl;                      	    
	    //print out the fields
            

cout<<"derivative of the pressure"<<endl;
            Array<OneD, NekDouble> dP_re = Array<OneD, NekDouble>(nq1D,0.0);
            Array<OneD, NekDouble> dP_im = Array<OneD, NekDouble>(nq1D,0.0);

            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS,Rephysreg,dP_re);   
            Ilayer->GetPlane(1)->PhysDeriv(MultiRegions::eS,Imphysreg,dP_im);   
            //attempt to smooth the derivative:
            Array<OneD, NekDouble> dP_re_coeffs (Nregcoeffs); 
            Array<OneD, NekDouble> dP_im_coeffs (Nregcoeffs);   
            Vmath::Vcopy(nq1D,dP_re,1,Ilayer->GetPlane(0)->UpdatePhys(),1);
            
	    //outfieldx->FwdTrans(Ilayer->GetPlane(0)->UpdatePhys(), dP_re_coeffs);
	    //Ilayer->GetPlane(1)->FwdTrans(dP_im, dP_im_coeffs);
	    //Ilayer->GetPlane(0)->BwdTrans(dP_re_coeffs, dP_re);
	    //Ilayer->GetPlane(1)->BwdTrans(dP_im_coeffs, dP_im);
            
            Array<OneD, NekDouble> dP_square = Array<OneD, NekDouble>(nq1D, 0.0);
cout<<"x"<<"  P_re"<<"  dP_re"<<"   streak"<<"   dstreak"<<"   pjump"<<endl;
	    for(int s=0; s<nq1D; s++)
            {
                dP_square[s] = dP_re[s]*dP_re[s] +dP_im[s]*dP_im[s];
//cout<<x0d[s]<<"    "<<x1d[s]<<"    "<<x2d[s]<<"   "<<Rephysdown[s]<<"   "<<dP_re[s]<<"  "<<dP_im[s]<<endl;
	    }            
cout<<"dim streak="<<stphysreg.num_elements()<<endl;
	    Array<OneD, NekDouble> dUreg (nq1D);
	    //Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eN, stphysreg,dUreg);
	    for(int w=0; w<nq1D; w++)
	    {
	         dUreg[w] = nx[w]*stphysgradxreg[w] +ny[w]*stphysgradyreg[w]; 	    
	    }
/*
	    //attempting to smooth mu field...
cout<<"ncoeffs="<<Ndowncoeffs<<endl;	    
	    Array<OneD, NekDouble> dU_coeffs (Nregcoeffs);
	    outfieldx->->FwdTrans(dUdown, dU_coeffs);
	    outfieldx->->BwdTrans(dU_coeffs, dUdown);
*/	    
	    Array<OneD, NekDouble> mu53  (nq1D);
	    Array<OneD, NekDouble> d2v  (nq1D);	    
	    double pow = 1.0/3.0;
	    double base;
            for(int y=0; y<nq1D; y++)
            {
                base =dUreg[y];
                mu53[y] = std::pow ((base*base),pow);
                mu53[y] = 1/(base*mu53[y]);
                d2v[y] = mu53[y]*dP_square[y];
                
            }

            
	    //attempting to smooth field...	    
	    Array<OneD, NekDouble> prod_coeffs (Nregcoeffs);
	    //outfieldx->FwdTrans(d2v, prod_coeffs);
	    //outfieldx->BwdTrans(prod_coeffs, d2v);            

	    Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, d2v ,d2v);
	    //attempting t smooth the vjump
	    outfieldx->FwdTrans(d2v, prod_coeffs);
	    outfieldx->BwdTrans(prod_coeffs, d2v);  	    
	    NekDouble n0 = 2.6789385347077476337;

	    Array<OneD, NekDouble> pjump(nq1D);
	    Array<OneD, NekDouble> vjump(nq1D);
	    for(int g=0; g<nq1D; g++)
	    {
                //double po = (2.09)**0.5;
                pjump[g] = n0*mu53[g]*dP_square[g] ;
                vjump[g] = n0*d2v[g];
cout<<setw(14)<<x0d[g]<<"       "<<Rephysreg[g]<<"      "<<dP_re[g]<<"     "
<<setw(13)<<"      "<<stphysreg[g]<<"     "<<dUreg[g]<<"    "<<
mu53[g]<<"       "<<setw(13)<<curv[g]<<"    "<<dP_square[g]<<"     "<<
pjump[g]<<setw(13)<<"        "<<vjump[g]<<endl;
//cout<<setw(14)<<x0d[g]<<"     "<<stphysgradxreg[g]<<"     "<<stphysgradyreg[g]
//<<"     "<<dUreg[g]<<endl;                                          
	    }

//PAY ATTENTION to the sign (vjump -/+ pjump)*t!!!!!!!!!!
            for(int j=0; j<nq1D; j++)
            {
/*            	              	    
            	(outfieldx->UpdatePhys())[j] = 
            	  (vjump[j]-pjump[j])*tx[j];
            	(outfieldy->UpdatePhys())[j] =
            	   (vjump[j]-pjump[j])*ty[j];		   
*/		
            }
          

// gamma(1/3)= 2.6789385347077476337            
            // need the tangents related to the expList1D outfieldx[region]
/*
           LocalRegions::SegExpSharedPtr  bndSegExp =  boost::dynamic_pointer_cast<LocalRegions::SegExp>(outfieldx[region]->GetExp(0)); 
           tangents = (bndSegExp)->GetMetricInfo()->GetEdgeTangent();
           for(int w=0; w<nq1D; w++)
           {       
cout<<"tangent tx="<<tangents[0][w]<<" ty="<<tangents[1][w]<<"  x="<<x0[w]<<"  y="<<x1[w]<<endl
           }
*/           
	}
	
	
	
	
	
	
       void  Manipulate(Array<OneD, int> Iregions,int coordim, SpatialDomains::MeshGraphSharedPtr &mesh,
       	       SpatialDomains::BoundaryConditions &bcs,
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &infields, 
       	       Array<OneD,MultiRegions::ExpList1DSharedPtr> &outfieldx,
       	       Array<OneD,MultiRegions::ExpList1DSharedPtr> &outfieldy,
       	       MultiRegions::ExpListSharedPtr &streak)
       
       {
       	       
       	       
       }       	       
	
	
	
	
	
	
	
	
	
	
	
	
	void WriteBcs(string variable, int region, string fieldfile,SpatialDomains::MeshGraphSharedPtr &mesh,
		MultiRegions::ContField1DSharedPtr &outregionfield)
	{			
		string   outfile = fieldfile.substr(0,fieldfile.find_last_of("."));
		outfile +="_"+variable+"_";
    		char ibnd[16]="";
    		sprintf(ibnd,"%d",region);
    		outfile +=ibnd;		
    		string   endfile(".bc");    		
    		outfile += endfile;
    		Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(1);   
                outregionfield->FwdTrans(outregionfield->GetPhys(),outregionfield->UpdateCoeffs()); 
                fieldcoeffs[0] = outregionfield->UpdateCoeffs();		
		std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                = outregionfield->GetFieldDefinitions();               
                std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
                // copy Data into FieldData and set variable             
            	//for(int i = 0; i < FieldDef.size(); ++i)
          	//{                
            	    //FieldDef[i]->m_fields.push_back(variable);
            	    FieldDef[0]->m_fields.push_back(variable);            	    
            	    //outregionfield->AppendFieldData(FieldDef[i], FieldData[i]);
            	    outregionfield->AppendFieldData(FieldDef[0], FieldData[0]);            	    
            	    //outregionfield->AppendFieldData(FieldDef[i], FieldData[i],fieldcoeffs[0]);
            	//}
cout<<"fielddata"<<FieldDef.size()<<endl;             	
            	mesh->Write(outfile,FieldDef,FieldData);            		
	}
    
        void WriteFld(string outfile,SpatialDomains::MeshGraphSharedPtr &mesh, Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
        {
            //variables array dimension is defined by fields!!!	
            Array<OneD, std::string> variables(fields.num_elements());   
            if( variables.num_elements()==2)
            {            	    
               variables[0]="u";            
               variables[1]="v";
            }
            else if( variables.num_elements()==3)
            {
                variables[0]="u";
                variables[1]="v";
                variables[2]="w";
            }
            else
            {
            	    ASSERTL0(false, " something goes wrong... ");
            }
            Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(fields.num_elements());
           
            std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
            = fields[0]->GetFieldDefinitions();
         
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
            //this cycle is needed in order to use the the appendfielddata with coeffs function!!!
            for(int i = 0; i < fields.num_elements(); ++i)
            {
            	if (fields[i]->GetPhysState()==true)
            	{	
                    //fields[i]->FwdTrans_IterPerExp(fields[i]->GetPhys(),fields[i]->UpdateCoeffs());
                    fields[i]->FwdTrans(fields[i]->GetPhys(),fields[i]->UpdateCoeffs());
                    
                }
                fieldcoeffs[i] = fields[i]->UpdateCoeffs();
            }
            
            // copy Data into FieldData and set variable
            for(int j = 0; j < fieldcoeffs.num_elements(); ++j)
            {
            	//fieldcoeffs[j] = fields[j]->UpdateCoeffs();             	    
            	for(int i = 0; i < FieldDef.size(); ++i)
          	{
          	     // Could do a search here to find correct variable
          	     FieldDef[i]->m_fields.push_back(variables[j]);
          	     //remark: appendfielddata without fieldcoeffs input gives always the first field!!!
          	     //fields[0]->AppendFieldData(FieldDef[i], FieldData[i]);          	     
          	     fields[0]->AppendFieldData(FieldDef[i], FieldData[i],fieldcoeffs[j]);
          	}
            }
            mesh->Write(outfile,FieldDef,FieldData);            
            
        }
    
