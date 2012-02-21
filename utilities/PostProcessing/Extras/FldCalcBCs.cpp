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
    void Extractlayerdata(Array<OneD, int> Iregions, int coordim, 
            SpatialDomains::MeshGraphSharedPtr &mesh,   
            LibUtilities::SessionReaderSharedPtr &session,
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

//REMARK:
/**
To create another bcs file comment everything from '//start' to '//end'
decomment the lines which follow '//decomment'
 and run with the command: ./FldCalcBCs  meshfile fieldfile  fieldfile
**/


    
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
   

    std::string solvtype = vSession->GetSolverInfo("SOLVERTYPE");
    if(solvtype == "CoupledLinearisedNS")
    {
      	    
         SetFields(graphShPt,boundaryConditions,vSession,fields,nfields-1);
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
         fields[lastfield] = Exp3DH1;        
//end       
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
    	  if(bndConditions[r]->GetUserDefined()==SpatialDomains::eCalcBC)
    	  {
    	  	  lastIregion=r;
    	  	  Iregions[r]=r;
    	  	  Ilayers[nIregions]=r;
    	  	  nIregions++;
//cout<<"Iregion="<<r<<endl;    	  	  
    	  }    	  
    } 
    ASSERTL0(nIregions>0,"there is any boundary region with the tag USERDEFINEDTYPE=""CalcBC"" specified");



    // import the streak 
    //import the streak field
    //add y to obtain the streak w=w'+y and write fld file           
    MultiRegions::ExpListSharedPtr streak; 
//start
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
//end   	
    //---------------------------------------------------------------   

    //set 1D output fields:
    //initialise fields
    const SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();

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
    Extractlayerdata(Ilayers,coordim, graphShPt,vSession, bcs, fields, outfieldx,outfieldy,streak);       	       

    //--------------------------------------------------------------------------------------


  

    //write bcs files: one for each I region and each variable
    string var="u";
    WriteBcs(var,lastIregion, fieldfile,graphShPt,outfieldx);
    var="v";
    WriteBcs(var,lastIregion, fieldfile,graphShPt,outfieldy);      	      	      
    
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

       void  Extractlayerdata(Array<OneD, int> Iregions, int coordim, 
            SpatialDomains::MeshGraphSharedPtr &mesh,   
            LibUtilities::SessionReaderSharedPtr &session,
    	    SpatialDomains::BoundaryConditions &bcs,
    	    Array<OneD,MultiRegions::ExpListSharedPtr> &infields,  
    	    MultiRegions::ContField1DSharedPtr &outfieldx,
       	    MultiRegions::ContField1DSharedPtr &outfieldy,
       	    MultiRegions::ExpListSharedPtr &streak)
        {
            //1 I regions is expected: (the layer is the last region)
            ASSERTL0(Iregions.num_elements()==3, "something wrong with the number of I layers");
            //int Iup =Iregions[0];
            int Ireg =Iregions[2];
cout<<"layer region="<<Ireg<<endl;   

    
            //take the pressure from infields
            MultiRegions::ExpListSharedPtr pressure =infields[3];
            Array<OneD, MultiRegions::ExpListSharedPtr> Iexp =infields[0]->GetBndCondExpansions();
            MultiRegions::ExpListSharedPtr Ilayer = Iexp[Ireg];

            int nq1D = Ilayer->GetPlane(0)->GetTotPoints();  
//decomment            int nq1D = Ilayer->GetTotPoints();  
            Array<OneD,NekDouble> x0d(nq1D);
            Array<OneD,NekDouble> x1d(nq1D);  
            Array<OneD,NekDouble> x2d(nq1D); 

            Ilayer->GetPlane(0)->GetCoords(x0d,x1d,x2d); 
//decomment            Ilayer->GetCoords(x0d,x1d,x2d); 
//start
            //take the streak data along the layers (getbondcondexpansions does not work!!!)          
            //Array<OneD, MultiRegions::ExpListSharedPtr> Istreak = streak->GetBndCondExpansions();
            //Array<OneD, MultiRegions::ExpListSharedPtr> Istreak =
              //                                     infields[1]->GetPlane(0)->GetBndCondExpansions();
            //MultiRegions::ExpListSharedPtr Istreakup = Istreak[Iup];
            //MultiRegions::ExpListSharedPtr Istreakdown = Istreak[Idown];

            
          
            
            
            int np = streak->GetTotPoints();
            Array<OneD, NekDouble> gradx(np);
            Array<OneD, NekDouble> grady(np); 
                     
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


           

            Array<OneD, NekDouble> nx(nq1D);
            Array<OneD, NekDouble> ny(nq1D);
            Array<OneD, NekDouble> tx(nq1D);
            Array<OneD, NekDouble> ty(nq1D);            
    
            
            
            


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
  	    
  	    //jacobian for the layer, std_deriv  
  	    Array<OneD, NekDouble> Jac (nq1D);
  	    Array<OneD, NekDouble> dPrestd (nq1D);
  	    
  	    //Array<OneD, NekDouble> stphysup (nq1D);
  	    Array<OneD, NekDouble> stphysreg (nq1D);
  	    Array<OneD, NekDouble> stphysgradxreg (nq1D);
  	    Array<OneD, NekDouble> stphysgradyreg (nq1D);
 
            //metrics definitions
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > deriv;   
            Array<OneD, NekDouble> dxde1 (nq1D);
	    Array<OneD, NekDouble> dyde1 (nq1D);
	    Array<OneD, NekDouble> dxde2 (nq1D);
	    Array<OneD, NekDouble> dyde2 (nq1D);
       	    Array<OneD, NekDouble> gmat0 (nq1D);
	    Array<OneD, NekDouble> gmat3 (nq1D);
            Array<OneD, NekDouble> gmat1 (nq1D);
	    Array<OneD, NekDouble> gmat2 (nq1D);
            
            Array<OneD, int> Elmtid(Ilayer->GetPlane(0)->GetExpSize());
            Array<OneD, int> Edgeid(Ilayer->GetPlane(0)->GetExpSize());
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
            	Array<OneD, NekDouble> stgradxcoeffsreg (Nregcoeffs);  
            	Array<OneD, NekDouble> stgradycoeffsreg (Nregcoeffs);  

                        	
            	
		//define nq points for each edge
            	int nqedge =nq1D/Icompreg->size();
            	
		Array<OneD,NekDouble> x0edge(nqedge);
		Array<OneD,NekDouble> x1edge(nqedge);    
                Array<OneD, NekDouble> Rephysedgereg(nqedge);
                Array<OneD, NekDouble> Imphysedgereg(nqedge);
		Array<OneD, NekDouble> stgradxcoeffsedge(nqedge);
		Array<OneD, NekDouble> stgradycoeffsedge(nqedge);

        	
             	Array<TwoD, const NekDouble>  gmat;
                Array<OneD, NekDouble>  gmat0_2D (nqedge*nqedge);
                Array<OneD, NekDouble>  gmat3_2D (nqedge*nqedge);
                Array<OneD, NekDouble>  gmat1_2D (nqedge*nqedge);
                Array<OneD, NekDouble>  gmat2_2D (nqedge*nqedge);   	             	
           
                
                
                
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
//cout<<"edge="<<EIDreg<<"  CHECK OFFSET="<<offsetregIExp<<endl;                                   
                     elementreg = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(mesh)
                                    ->GetElementsFromEdge(segmentGeomreg);
	     	     int elmtidreg = ((*elementreg)[0]->m_Element)->GetGlobalID();   
                     int dimelmtreg = ((*elementreg)[0]->m_Element)->GetNumEdges();
                     Elmtid[k] = elmtidreg;
		     Edgeid[k] = EIDreg;
                     orientreg =(boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>
                     	           ((*elementreg)[0]->m_Element))->GetEorient((*elementreg)[0]
                     	           	   ->m_EdgeIndx);

                     gmat = pressure->GetPlane(0)->GetExp(elmtidreg)->GetMetricInfo()->GetGmat();
                     int nq2D = nqedge*nqedge;                     	           
                     
		     //	setup map between global and local ids
                     map<int, int>EdgeGIDreg;
		     int id,cnt,f;
                     for(cnt=f=0; f<dimelmtreg; f++)
		     {
			   id = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>
                     	           ((*elementreg)[0]->m_Element))->GetEid(f);
			   EdgeGIDreg[id]=cnt++;		     	     
		     }
//cout<<" GIDreg="<<elmtidreg<<"  num edges in this element="<<dimelmtreg<<endl;
	             int localidedge = EdgeGIDreg.find(EIDreg)->second;
//cout<<"locaidedge from m_EdgeInx ="<<localidedge<<endl;	             
	             localidedge = (*elementreg)[0]->m_EdgeIndx;
//cout<<"check map local id of Eid="<<EIDreg<<"  is="<<localidedge<<endl;	             


		     //set bmaps -------------
		     Array<OneD, unsigned int > bmapedgereg;
		     Array<OneD, int > signreg;
		     pressure->GetExp(elmtidreg)->GetEdgeToElementMap(EdgeGIDreg.find(EIDreg)->second,
		     	     orientreg,bmapedgereg,signreg);		     		     

//cout<<"SIGN down="<<signreg[0]<<endl;
		     //Extracting the coeffs...
		     Array<OneD, NekDouble> Recoeffs = pressure->GetPlane(0)->GetCoeffs();
		     Array<OneD, NekDouble> Imcoeffs = pressure->GetPlane(1)->GetCoeffs();  
		     Array<OneD, NekDouble> stcoeffs = streak->GetCoeffs();
                     //OFFSET EDGE REFERING TO ELEMENT EXPLIST MISSING		     
		     int Reoffsetreg = pressure->GetPlane(0)->GetCoeff_Offset(elmtidreg);
		     int Imoffsetreg = pressure->GetPlane(1)->GetCoeff_Offset(elmtidreg);
		     int stoffsetreg = streak->GetCoeff_Offset(elmtidreg);
		     
//cout<<" Re offsetdown="<<Reoffsetreg<<"   Im offset="<<Imoffsetreg<<endl;		     
		     



		     int ReNumoffset = Ilayer->GetPlane(0)->GetNcoeffs();
                     int ImNumoffset = Ilayer->GetPlane(1)->GetNcoeffs();	
//cout<<" OFFSET TEst Re="<<ReNumoffset<<"  Im="<<ImNumoffset<<endl;                     
		     int Reoffsetedgereg = Ilayer->GetPlane(0)->GetCoeff_Offset(k);
//cout<<"k="<<k<<"  Offsetedge="<<Reoffsetedgereg<<endl;		     
		     
		     //hypothesis: 2 planes (0 Re, 1 Im)
		     int nplanecoeffs = pressure->GetExp(elmtidreg)->GetNcoeffs();
		     Array<OneD, NekDouble> Recoeffsedgereg(bmapedgereg.num_elements());	
		     Array<OneD, NekDouble> Imcoeffsedgereg(bmapedgereg.num_elements());
		     Array<OneD, NekDouble> stcoeffsedgereg(bmapedgereg.num_elements());
		     
		     
		     //grad edge coeffs		     
//cout<<" Recoeffsedgeup num elements="<<Recoeffsedgereg.num_elements()
//<<"  n bmapedge="<<bmapedgereg.num_elements()<<endl;	



		     //ncoeffs per edge fwd transf..
		     int ncoeffs2Delmt = pressure->GetPlane(0)->GetExp(elmtidreg)->GetNcoeffs();	     
                  
            



            	     int nqed = Ilayer->GetPlane(0)->GetTotPoints();

		     for(int d=0; d<bmapedgereg.num_elements(); d++)
		     {
			   //pressure	 
		     	   Recoeffsedgereg[d]=Recoeffs[Reoffsetreg+bmapedgereg[d]];		
                           Imcoeffsedgereg[d]=Imcoeffs[Imoffsetreg+bmapedgereg[d]];
 		     	   Recoeffsreg[offsetregIExp +d]  = Recoeffsedgereg[d];
		     	   Imcoeffsreg[offsetregIExp +d]  = Imcoeffsedgereg[d];
                          
                           //streak coeffs
                           stcoeffsedgereg[d] = stcoeffs[stoffsetreg + bmapedgereg[d]];
                           stgradxcoeffsedge[d] = gradxcoeffs[stoffsetreg + bmapedgereg[d]];
                           stgradycoeffsedge[d] = gradycoeffs[stoffsetreg + bmapedgereg[d]];                           
		     	   stcoeffsreg[offsetregIExp +d] = stcoeffsedgereg[d];
		     	   stgradxcoeffsreg[offsetregIExp +d] = stgradxcoeffsedge[d];
		     	   stgradycoeffsreg[offsetregIExp +d] = stgradycoeffsedge[d];      
		     	   
		     }
		     
		     
		     
		     
		     
		     
		     
		     //store the jacobian, std_deriv                      	           
                     Array<OneD, NekDouble> jacedge (nqedge);
                     Array<OneD, NekDouble> Pre_edge (nqedge);                     
                     Array<OneD, NekDouble> dPre_edge (nqedge);
                     jacedge = Ilayer->GetPlane(0)->GetExp(k)->GetMetricInfo()->GetJac();


                     //check if the metrix is the same for streak and Ilayer obj
                     StdRegions::StdExpansion1DSharedPtr edgestdexp;
                     edgestdexp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D>  (
                     	     Ilayer->GetPlane(0)->GetExp(k)  );
                     edgestdexp->BwdTrans(Recoeffsedgereg, Pre_edge);                     
//cout<<"call PhysTensorDERIV......"<<endl;                     
                     edgestdexp->StdPhysDeriv(Pre_edge, 
                     	      dPre_edge);
                     for(int q=0; q<nqedge; q++)
                     {
                     	  Jac[k*nqedge +q]= jacedge[q];
                     	  dPrestd[k*nqedge +q] = dPre_edge[q];
                     	  
                     }		     
		     
		     
		     
		     
		     
		     
		     
		     
//cout<<"Nlayercoeffs="<<Recoeffsreg.num_elements()<<"   Nlayerphys="<<Rephysreg.num_elements()<<endl;
		     
		     
		     
//cout<<"extract normal for edge="<<k<<endl;	
	     
//cout<<"extract tangent for edge="<<k<<endl;
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

                          nx[k*nqedge +e] = normals[0][e];
                          ny[k*nqedge +e] = normals[1][e];
                          tx[k*nqedge +e] = tangents[0][e];
                          ty[k*nqedge +e] = tangents[1][e];
//cout<<" nx="<<normals[0][e]<<"    ny="<<normals[1][e]<<"   e="<<e<<endl; 
//cout<<"tannn tx="<<tangents[0][e]<<"   ty="<<tangents[1][e]<<endl;
                     }
		


		     //extract metrics factors...
                     //int nq2D= nqedge*nqedge;
                     deriv = pressure->GetPlane(0)->GetExp(elmtidreg)->GetMetricInfo()->GetDeriv();

                     Array<OneD, NekDouble> derivelmt(nq2D);
                     int offsetregphys = Ilayer->GetPlane(0)->GetPhys_Offset(k);   
                     Array<OneD, NekDouble> deriv_i(nqedge);
                     Vmath::Vcopy(nq2D, &(deriv[0][0][0]),1, &(derivelmt[0]),1);                
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge, derivelmt, deriv_i);
                     Vmath::Vcopy(nqedge, &(deriv_i[0]),1, &(dxde1[offsetregphys]),1);                
                     Vmath::Zero(nqedge, deriv_i,1);
                     Vmath::Zero(nq2D, derivelmt,1);
                     Vmath::Vcopy(nq2D, &(deriv[1][1][0]),1, &(derivelmt[0]),1);                
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge, derivelmt, deriv_i);
                     Vmath::Vcopy(nqedge, &(deriv_i[0]),1, &(dyde2[offsetregphys]),1);                 
                     Vmath::Zero(nqedge, deriv_i,1);
                     Vmath::Zero(nq2D, derivelmt,1);
                     Vmath::Vcopy(nq2D, &(deriv[0][1][0]),1, &(derivelmt[0]),1);                
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge, derivelmt, deriv_i);
                     Vmath::Vcopy(nqedge, &(deriv_i[0]),1, &(dyde1[offsetregphys]),1);                 
                     Vmath::Zero(nqedge, deriv_i,1);
                     Vmath::Zero(nq2D, derivelmt,1);                
                     Vmath::Vcopy(nq2D, &(deriv[1][0][0]),1, &(derivelmt[0]),1);                
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge, derivelmt, deriv_i);
                     Vmath::Vcopy(nqedge, &(deriv_i[0]),1, &(dxde2[offsetregphys]),1);                         
                  
                     Array<OneD, NekDouble> gmati (nqedge); 
                     Vmath::Vcopy(nq2D, &(gmat[0][0]),1, &(gmat0_2D[0]),1);
                     Vmath::Vcopy(nq2D, &(gmat[3][0]),1, &(gmat3_2D[0]),1); 
                     Vmath::Vcopy(nq2D, &(gmat[1][0]),1, &(gmat1_2D[0]),1);
                     Vmath::Vcopy(nq2D, &(gmat[2][0]),1, &(gmat2_2D[0]),1);                 
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge,gmat0_2D, gmati);
                     Vmath::Vcopy(nqedge, &(gmati[0]),1, &(gmat0[offsetregphys]),1);
                     Vmath::Zero(nqedge, gmati,1);
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge,gmat1_2D, gmati);
                     Vmath::Vcopy(nqedge, &(gmati[0]),1, &(gmat1[offsetregphys]),1);
                     Vmath::Zero(nqedge, gmati,1);
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge,gmat2_2D, gmati);
                     Vmath::Vcopy(nqedge, &(gmati[0]),1, &(gmat2[offsetregphys]),1);
                     Vmath::Zero(nqedge, gmati,1);
                     pressure->GetPlane(0)->GetExp(elmtidreg)->GetEdgePhysVals(localidedge,gmat3_2D, gmati);
                     Vmath::Vcopy(nqedge, &(gmati[0]),1, &(gmat3[offsetregphys]),1);

                }       
                //bwd transform:
		Ilayer->GetPlane(0)->BwdTrans(Recoeffsreg, Rephysreg);
		Ilayer->GetPlane(1)->BwdTrans(Imcoeffsreg, Imphysreg);
		//streak phys values:
		Ilayer->GetPlane(0)->BwdTrans(stcoeffsreg, stphysreg);
		Ilayer->GetPlane(0)->BwdTrans_IterPerExp(stgradxcoeffsreg, stphysgradxreg);
		Ilayer->GetPlane(0)->BwdTrans_IterPerExp(stgradycoeffsreg, stphysgradyreg);
  
            }  
      



            



             //TANGENTS WRONG WITH 3DHOMOGENEOUS 1D...	     
             //calculate the curvature:     

             
	     //attempt to smooth the tangent components:
             Array<OneD, NekDouble> tcoeffs (Nregcoeffs,0.0);
	     //outfieldx->FwdTrans(tx, tcoeffs);
	     //outfieldx->BwdTrans(tcoeffs, tx);
             //Vmath::Zero(Nregcoeffs, tcoeffs,1);
	     //outfieldx->FwdTrans(ty, tcoeffs);
    	     //outfieldx->BwdTrans(tcoeffs, ty);
             //Vmath::Zero(Nregcoeffs, tcoeffs,1);
             Array<OneD, NekDouble> dtx(nq1D,0.0);
             Array<OneD, NekDouble> dty(nq1D,0.0);    
             Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, tx, dtx);
             Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, ty, dty);
	     //attempt to smooth the tangent derivative:
	     outfieldx->FwdTrans(dtx, tcoeffs);
	     outfieldx->BwdTrans(tcoeffs, dtx);  
             Vmath::Zero(Nregcoeffs, tcoeffs,1);
	     outfieldx->FwdTrans(dty, tcoeffs);
	     outfieldx->BwdTrans(tcoeffs, dty);                           
             Array<OneD, NekDouble> f_z(nq1D,0.0);
             Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eX, x1d, f_z);
             Array<OneD, NekDouble> fcoeffs(Nregcoeffs,0.0);
	     //outfieldx->FwdTrans(f_z, fcoeffs);
	     //outfieldx->BwdTrans(fcoeffs, f_z); 
             Vmath::Zero(Nregcoeffs, fcoeffs,1);
             Array<OneD, NekDouble> f_zz(nq1D,0.0);

             Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eX, f_z, f_zz);
	     //outfieldx->FwdTrans(f_zz, fcoeffs);
	     //outfieldx->BwdTrans(fcoeffs, f_zz); 
             Array<OneD, NekDouble> delta(nq1D);
             Array<OneD, NekDouble> curv_unsigned(nq1D,0.0);
             Array<OneD, NekDouble> curv(nq1D,0.0);
        

//cout<<"x     y      tx      ty      dtx        dty       curv"<<endl;     
             for(int t=0; t<nq1D; t++)
	     {

	            curv_unsigned[t]=sqrt(dtx[t]*dtx[t] +dty[t]*dty[t]);     
                    delta[t] = 1+f_z[t]*f_z[t];
                    //Hall's definition...
                    curv[t] = f_zz[t]/sqrt(delta[t]*delta[t]*delta[t]);
//cout<<setw(13)<<x0d[t]<<"     "<<setw(13)<<x1d[t]<<"      "<<setw(13)<<tx[t]<<setw(13)<<"     "<<ty[t]<<"     "
//<<nx[t]<<"    "<<ny[t]<<"      "<<setw(13)<<dtx[t]<<"      "<<dty[t]<<"     "<<curv[t]<<endl;
             }		     
	

	    
           //attempt to smooth the curvature
           Array<OneD, NekDouble> curv_coeffs (Nregcoeffs);
	   outfieldx->FwdTrans(curv, curv_coeffs);
	   outfieldx->BwdTrans(curv_coeffs, curv); 


           //normalise the pressure  norm*sqrt[ (\int Psquare)/Area]=1 =>

           NekDouble norm;
           NekDouble tmp = pressure->GetPlane(0)->L2();
           norm = tmp*tmp;
           tmp = pressure->GetPlane(1)->L2();
           norm += tmp*tmp;
           Array<OneD, NekDouble> I (2*np, 1.0);             
           NekDouble Area = pressure->GetPlane(0)->PhysIntegral(I);           
           norm = sqrt(Area/norm);


           //norm*pressure
	   Vmath::Smul(nq1D,norm,Rephysreg,1,Rephysreg,1);                       
	   Vmath::Smul(nq1D,norm,Imphysreg,1,Imphysreg,1);      
                    	    
	    //print out the fields
            

//cout<<"derivative of the pressure"<<endl;
            Array<OneD, NekDouble> dP_re = Array<OneD, NekDouble>(nq1D,0.0);
            Array<OneD, NekDouble> dP_im = Array<OneD, NekDouble>(nq1D,0.0);
            Array<OneD, NekDouble> dP_square2 = Array<OneD, NekDouble>(nq1D,0.0);

            //attempt to smooth the normal
            //Array<OneD, NekDouble> ncoeffs(Nregcoeffs,0.0);
	    //outfieldx->FwdTrans(nx, ncoeffs);
	    //outfieldx->BwdTrans(ncoeffs, nx);
            //Vmath::Zero(Nregcoeffs, ncoeffs,1);
	    //outfieldx->FwdTrans(ny, ncoeffs);
    	    //outfieldx->BwdTrans(ncoeffs, ny);
            //Vmath::Zero(Nregcoeffs, ncoeffs,1);        

	    Array<OneD, NekDouble> dUreg (nq1D);
	    //Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eN, stphysreg,dUreg);
	    for(int w=0; w<nq1D; w++)
	    {
	         dUreg[w] = nx[w]*stphysgradxreg[w] +ny[w]*stphysgradyreg[w]; 	    
	    }

 

            Array<OneD, NekDouble> dUcoeffs(Nregcoeffs,0.0);
            //outfieldx->FwdTrans(dUreg, dUcoeffs);
            //outfieldx->BwdTrans(dUcoeffs, dUreg);

                   

	    Array<OneD, NekDouble> jaccoeffs(Nregcoeffs);
            Array<OneD, NekDouble> Jacreg (nq1D);   
            Array<OneD, NekDouble> dPrecoeffs (Nregcoeffs);
            //attempt to smooth the jac and the dPre
	    outfieldx->FwdTrans(Jac, jaccoeffs);
	    outfieldx->BwdTrans(jaccoeffs, Jacreg);
            //outfieldx->FwdTrans(dPre, dPrecoeffs);
            //outfieldx->BwdTrans(dPrecoeffs, dPre);
	          
            
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS,Rephysreg,dP_re);          
            Ilayer->GetPlane(1)->PhysDeriv(MultiRegions::eS,Imphysreg,dP_im);   
            //attempt to smooth the derivative:
            Array<OneD, NekDouble> dP_re_coeffs (Nregcoeffs,0.0); 
            Array<OneD, NekDouble> dP_im_coeffs (Nregcoeffs,0.0);   
            Array<OneD, NekDouble> dP_imtest (nq1D);
	    outfieldx->FwdTrans(dP_re, dP_re_coeffs);
	    outfieldx->FwdTrans(dP_im, dP_im_coeffs);
	    outfieldx->BwdTrans(dP_re_coeffs, dP_re);
	    outfieldx->BwdTrans(dP_im_coeffs, dP_im);	             
            
            
            
            Array<OneD, NekDouble> dP_square = Array<OneD, NekDouble>(nq1D, 0.0);         
cout<<"x"<<"  P_re"<<"  dP_re"<<"   streak"<<"   dstreak"<<"   pjump"<<endl;
	    for(int s=0; s<nq1D; s++)
            {
                dP_square[s] = dP_re[s]*dP_re[s] +dP_im[s]*dP_im[s];
	    } 



 
	    //attempt to smooth the pressure     
	    Array<OneD, NekDouble> dPsquare_coeffs (Nregcoeffs,0.0);
	    outfieldx->FwdTrans(dP_square, dPsquare_coeffs);
	    outfieldx->BwdTrans(dPsquare_coeffs, dP_square);	    


	    Array<OneD, NekDouble> mu53  (nq1D,0.0);
	    Array<OneD, NekDouble> d2v  (nq1D,0.0);
	    Array<OneD, NekDouble> prod  (nq1D,0.0);	    	    
		    
	    double pow = 1.0/3.0;
	    double base;
            for(int y=0; y<nq1D; y++)
            {
                base =dUreg[y];
                mu53[y] = std::pow ((base*base),pow);
                mu53[y] = 1/(base*mu53[y]);
                //prod[y] = mu53[y]*dP_square[y];                                
            }
            //attempting to smooth mu53
            Array<OneD, NekDouble> mucoeffs(Nregcoeffs,0.0);
            outfieldx->FwdTrans(mu53, mucoeffs);
            outfieldx->BwdTrans(mucoeffs, mu53);
            Vmath::Vmul(nq1D, mu53,1, dP_square,1, prod,1);
            
	    //attempting to smooth field...	    
	    Array<OneD, NekDouble> prod_coeffs (Nregcoeffs,0.0);
	    outfieldx->FwdTrans(prod, prod_coeffs);
	    outfieldx->BwdTrans(prod_coeffs, prod);
            Vmath::Zero(Nregcoeffs,prod_coeffs,1);	           

	    Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, prod ,d2v);
	    //attempting to smooth the vjump
            Vmath::Zero( Nregcoeffs, prod_coeffs,1);
	    outfieldx->FwdTrans(d2v, prod_coeffs);
	    outfieldx->BwdTrans(prod_coeffs, d2v);
            //test prod deriv... add background..
            //Array<OneD, NekDouble> prod_b(nq1D);
            //Array<OneD, NekDouble> dersprod_b(nq1D);
            //Vmath::Sadd(nq1D,0.001, prod,1,prod_b,1);
	    //Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, prod_b ,dersprod_b);


	    Array<OneD, NekDouble> pjump(nq1D);
	    Array<OneD, NekDouble> vjump(nq1D);


            //test a sin curve....
/*
            Array<OneD, NekDouble> fun1D(nq1D);
            Array<OneD, NekDouble> derfun1D(nq1D);
            NekDouble pi =3.14159265;
            for(int e=0; e<nq1D; e++)
            {
                   fun1D[e]= std::sin((x0d[e]));
            }  
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, fun1D, derfun1D);        

            //add noise to a function
            Array<OneD, NekDouble> fun1D_b(nq1D);
            Array<OneD, NekDouble> dersfun1D_b(nq1D);
            int cnt=0;
            int cnt1=1;

            for(int v=0; v<nq1D/2; v++)
            {

                 fun1D_b[cnt] = fun1D[cnt] +0.000000;
//cout<<cnt<<"x="<<x0d[cnt]<<"  fun1D_b="<<fun1D_b[cnt]<<"   fun1D="<<fun1D[cnt]<<"  diff="<<
//fun1D[cnt]-fun1D_b[cnt]<<endl;
                 //ASSERTL0((fun1D_b[cnt]-fun1D[cnt])==0.0001, "problem");
                 cnt = cnt +2;
                 fun1D_b[cnt1] = fun1D[cnt1] -0.000000;
//cout<<cnt1<<"x="<<x0d[cnt1]<<"  fun1D_b="<<fun1D_b[cnt1]<<"  fun1D="<<fun1D[cnt1]<<"  diff="<<
//fun1D[cnt1]-fun1D_b[cnt1]<<endl;
                 //ASSERTL0((fun1D_b[cnt1]-fun1D[cnt1])==-0.0001, "problem");
                 cnt1 = cnt1 +2;
                 if( (v==((nq1D/2) -1)) && (nq1D%2==1) 
                   )
                 {
                    fun1D_b[cnt] = fun1D[cnt] +0.000000;
//cout<<cnt<<"x="<<x0d[cnt]<<"  fun1D_b="<<fun1D_b[cnt]<<"  fun1D="<<fun1D[cnt]<<"  diff="<<
//fun1D[cnt]-fun1D_b[cnt]<<endl;
                 }
            }

	    Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, fun1D_b ,dersfun1D_b); 
*/
/*

	    int nqedge =nq1D/Icompreg->size();
	    Array<OneD, NekDouble> phys_edge (nqedge);	    
	    Array<OneD, NekDouble> der_edge (nqedge);
            Array<OneD, NekDouble> der_fun1Db_std (nq1D);
            //extract stdderiv....

            for(int w=0; w<Icompreg->size() ; w++)
            {
            	 Vmath::Vcopy(nqedge,&(fun1D_b[w*nqedge]),1,&(phys_edge[0]),1);    
                 //check if the metrix is the same for streak and Ilayer obj
                 StdRegions::StdExpansion1DSharedPtr edgestdexp;
                 edgestdexp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D>  (
                     	     Ilayer->GetPlane(0)->GetExp(w)  );
                 edgestdexp->StdPhysDeriv(phys_edge, 
                     	      der_edge);                 
                 Vmath::Vcopy(nqedge, &(der_edge[0]),1, &(der_fun1Db_std[w*nqedge]),1);           
                 Vmath::Zero(nqedge, phys_edge,1);
                 Vmath::Zero(nqedge, der_edge,1);
/*
                 Vmath::Vcopy(nqedge, &(Rephysreg[w*nqedge]),1, &(phys_edge[0]),1);
                 edgestdexp->StdPhysDeriv(phys_edge, der_edge);
                 Vmath::Vcopy(nqedge, &(der_edge[0]),1, &(der_Pre_std[w*nqedge]),1);
                 Vmath::Zero(nqedge, phys_edge,1);
                 Vmath::Zero(nqedge, der_edge,1);
                 Vmath::Vcopy(nqedge, &(fun1D[w*nqedge]),1, &(phys_edge[0]),1);
                 edgestdexp->StdPhysDeriv(phys_edge, 
                     	      der_edge); 
                 Vmath::Vcopy(nqedge, &(der_edge[0]),1, &(derfun1D[w*nqedge]),1);                                   
                 Vmath::Zero(nqedge, phys_edge,1);
                 Vmath::Zero(nqedge, der_edge,1);

            }

*/


          

       
	   //final jump conditions
           //n0=2*pi*(2/3)^(2/3)*(-2/3)! = 2*pi*(2/3)^(2/3)*gamma(1/3)
           NekDouble gamma13 = 2.6789385347077476336556929409746776441286893779573011009;
           NekDouble n0 = 12.845424015;
           //use the session to get the values of rho,alpha...
           NekDouble rho ;
           //rho = 0.08;
           rho = session->GetParameter("RHO");
           NekDouble alpha;
	   alpha = 1;
           NekDouble alpha53=1;
           alpha53 = std::pow ((alpha*alpha),pow);
           alpha53 = 1/(alpha*alpha53);
           alpha53=1;
	   for(int c=0; c<nq1D; c++)
	   {
	       pjump[c] = -n0*alpha53*curv[c]*prod[c] ;  	   
	   }
	   Vmath::Smul(nq1D, n0,d2v,1, vjump,1);      
	   Vmath::Smul(nq1D, alpha53,vjump,1, vjump,1);  
	   

	   
           Array<OneD, NekDouble> txcoeffs (Nregcoeffs);
           Array<OneD, NekDouble> tycoeffs (Nregcoeffs);
	   //attempt to smooth the tangent components:
	   outfieldx->FwdTrans(tx, txcoeffs);
	   outfieldx->FwdTrans(ty, tycoeffs);
	   outfieldy->BwdTrans(txcoeffs, tx);              
	   outfieldy->BwdTrans(tycoeffs, ty);
	     
	     
//end
//PAY ATTENTION to the sign (vjump -/+ pjump)*t!!!!!!!!!!
            for(int j=0; j<nq1D; j++)
            {
//start   
         	              	    
            	(outfieldx->UpdatePhys())[j] = 
            	  rho*rho*(vjump[j]*tx[j]-pjump[j]);
            	(outfieldy->UpdatePhys())[j] =
            	   rho*rho*(vjump[j]*ty[j]-pjump[j]);		   
//end

//decomment
/*
            	(outfieldx->UpdatePhys())[j] = 
            	  20*sin(2*x0d[j]);
            	(outfieldy->UpdatePhys())[j] =
            	   20*2*cos(2*x0d[j])/3.14;
*/
            }
            //FINAL REFINEMENT:::
//start
      
            Array<OneD, NekDouble> finalcoeffs (Nregcoeffs);
            outfieldx->FwdTrans(outfieldx->GetPhys(), finalcoeffs);
            outfieldx->BwdTrans(finalcoeffs, outfieldx->UpdatePhys());      
            Vmath::Zero(Nregcoeffs,finalcoeffs,1);            
            outfieldy->FwdTrans(outfieldy->GetPhys(), finalcoeffs);
            outfieldy->BwdTrans(finalcoeffs, outfieldy->UpdatePhys());  



/*

	    //calculate J,K
            Array<OneD, NekDouble> lambda(nq1D);            
            Array<OneD, NekDouble> a(nq1D);
            Array<OneD, NekDouble> dPre_dz(nq1D);
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eX, Rephysreg, dPre_dz);
            Array<OneD, NekDouble> dPim_dz(nq1D);
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eX, Imphysreg, dPim_dz);
            Array<OneD, NekDouble> dPre_dz2(nq1D);
            Vmath::Vmul(nq1D, dPre_dz,1, dPre_dz,1, dPre_dz2,1);
            Array<OneD, NekDouble> dPim_dz2(nq1D);
            Vmath::Vmul(nq1D, dPim_dz,1, dPim_dz,1, dPim_dz2,1);
            Array<OneD, NekDouble> dP_dz2(nq1D);
            Vmath::Vadd(nq1D, dPre_dz2,1,dPim_dz2,1,dP_dz2,1);
            Array<OneD, NekDouble> ddP_dz2z(nq1D);
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eX, dP_dz2, ddP_dz2z);


            Array<OneD, NekDouble> delta5(nq1D);
            Array<OneD, NekDouble> a53(nq1D);
            pow = 1./3.;
            for(int e=0; e<nq1D; e++) 
            {
                //delta[e] = 1+ f_z[e]*f_z[e];
                delta5[e] = delta[e]*delta[e]*delta[e]*delta[e]*delta[e];
                lambda[e] = dUreg[e]/sqrt(delta[e]);
                a[e] = lambda[e]*alpha/delta[e];
                a53[e] = std::pow( (a[e]*a[e]), pow);
//cout<<"a^2/3="<<a53[e]<<endl;
                a53[e] = a[e]*a53[e];
            }

            Array<OneD, NekDouble> delta_z(nq1D);
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eX, delta, delta_z);
            Array<OneD, NekDouble> a_z(nq1D);
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eX, a, a_z);

            Array<OneD, NekDouble> J(nq1D);
            Array<OneD, NekDouble> K(nq1D);
            Array<OneD, NekDouble> f1(nq1D);
            Array<OneD, NekDouble> f2(nq1D);
            //NekDouble a53,delta5;
            for(int f=0; f<nq1D; f++)
            {
                 //a53 = std::pow(a[f],5./3.);
                 //delta5 = std::pow( delta[f], 5);
                 J[f] =( (n0*rho*rho)/(a53[f]*delta5[f])   )*
                        (  
                          ( (-7*delta_z[f])/(2*delta[f]) - (5*a_z[f] )/(3*a[f])  )*dP_dz2[f] +
                          ddP_dz2z[f]
                        );
                 K[f] = -(n0*rho*rho)/(a53[f]*delta5[f])*f_zz[f]*dP_dz2[f];
                 f1[f] = lambda[f]*(K[f]-delta[f]*f_z[f]*J[f]);
                 f2[f] = -lambda[f]*(f_z[f]*K[f] +delta[f]*J[f]);
            }


            
            
            //test dermu
            Array<OneD, NekDouble> dermu53(nq1D);
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, mu53, dermu53);

            //test d(dP_square)ds
            Array<OneD, NekDouble> ddP_square_ds(nq1D);
            Ilayer->GetPlane(0)->PhysDeriv(MultiRegions::eS, dP_square, ddP_square_ds);
*/




            NekDouble jactest;

	    for(int g=0; g<nq1D; g++)
	    {
                //pjump[g] = n0*curv[g]*mu53[g]*dP_square[g] ;       
                //pjump[g] = n0*mu53[g]*dP_square[g] ;
                //vjump[g] = n0*d2v[g];
  
//NBBB dPre is the stdDERIV!!! 
                //jactest = (gmat0[g]*tx[g] +gmat2[g]*ty[g])-(1/Jac[g]);;
/*
cout<<x0d[g]<<"       "<<
//stphysreg[g]<<"      "<<lambda[g]<<"      "<<
//delta[g]<<"      "<<
//f_z[g]<<"     "<<a53[g]<<"    "<<
//delta5[g]<<"     "<<delta_z[g]<<"     "<<a_z[g]<<"     "<<
//"      "<<ddP_dz2z[g]<<"        "<<dP_dz2[g]<<"        "<<
//dP_square[g]<<"        "<<dP_dz2[g]/delta[g]<<"       "<<
//ddP_square_ds[g]<<"         "<<ddP_dz2z[g]/(delta[g]*sqrt(delta[g]))<<"         "<<
//ddP_dz2z[g]/delta[g]<<"       "<<
//a[g]<<"        "<<a53[g]<<"       "<<
vjump[g]*tx[g]*sqrt(delta[g])<<"         "<<
J[g]<<"      "<<K[g]<<"       "<<pjump[g]<<"      "<<
f1[g]<<"      "<<f2[g]<<"      "<<
f_zz[g]/(sqrt(delta5[g]))<<"       "<<
//f_zz[g]/sqrt(delta[g]*delta[g]*delta[g])<<"      "
curv[g]<<"        "<<
outfieldx->GetPhys()[g]<<"      "
<<outfieldy->GetPhys()[g]
<<endl;
*/


cout<<setw(14)<<x0d[g]<<"      "<<setw(13)<<x1d[g]<<
"      "<<Rephysreg[g]<<"   "<<dP_re[g]<<"     "
<<stphysreg[g]<<"      "<<mu53[g]<<"       "<<curv[g]
<<setw(13)<<"      "<<dP_square[g]<<"      "<<setw(13)<<prod[g]
<<"       "<<
//(gmat0[g]*tx[g] +gmat2[g]*ty[g])<<"      "<<1/Jac[g]
//<<"       "<<dermu53[g]
//<<"        "<<derfun1D[g]
//tx[g]<<"     "<<ty[g]
//f1[g]<<"        "<<f2[g]
//prod_b[g]<<"      "<<dersprod_b[g]
//fun1D[g]<<"      "<<derfun1D[g]<<"       "<<
//der_fun1Db_std[g]<<"       "<<
//fun1D_b[g]<<"       "<<dersfun1D_b[g]
pjump[g]<<setw(13)<<"        "<<d2v[g]<<"       "
<<outfieldx->GetPhys()[g]<<"      "
<<outfieldy->GetPhys()[g]<<endl;
//cout<<"jactest="<<jactest<<endl;
                ASSERTL0(abs(jactest)<0.5, "jac 1D problem..");


	    }
         
// gamma(1/3)= 2.6789385347077476337            
            // need the tangents related to the expList1D outfieldx[region]





            for(int a=0; a<Elmtid.num_elements(); a++)
            {
cout<<"elmt id="<<Elmtid[a]<<"  edge id="<<Edgeid[a]<<endl;
	    } 
//end
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
    
