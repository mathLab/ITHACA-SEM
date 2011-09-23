#include <cstdio>
#include <cstdlib>

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
    void Manipulate(Array<OneD, int> Iregions, int coordim, SpatialDomains::MeshGraphSharedPtr &mesh,   
    	    SpatialDomains::BoundaryConditions &bcs,
    	    Array<OneD,MultiRegions::ExpListSharedPtr> &infields,  
    	    Array<OneD,MultiRegions::ExpListSharedPtr> &outfieldx,
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &outfieldy);
    void WriteBcs(string variable, int region, string fieldfile, SpatialDomains::MeshGraphSharedPtr &mesh, 
    	    MultiRegions::ExpListSharedPtr &outregionfield);
    void WriteFld(string outfile, SpatialDomains::MeshGraphSharedPtr &mesh, Array<OneD, MultiRegions::ExpListSharedPtr> &fields );
    
    int i,j;
    if(argc != 3)
    {
        fprintf(stderr,"Usage: ./FldCalcBCs  meshfile fieldfile  \n");
        exit(1);
    }
 
    //----------------------------------------------
    string meshfile(argv[argc-2]);
  
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
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------
    
    //Import pressure file.
/*    
    string pfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> pdef;
    vector<vector<NekDouble> > pdata;
    graphShPt->Import(pfile,pdef,pdata);
*/    
    //----------------------------------------------


    // Define Expansion    
    int nfields; 
    nfields = fielddef[0]->m_fields.size(); 
    Array<OneD, MultiRegions::ExpListSharedPtr> fields;  
cout<<nfields<<endl;    
/*   
    if(vSession->GetSolverInfo("SOLVERTYPE")=="CoupledLinearisedNS")
    {
cout<<"Coupled solver nfields does not include the pressure"<<endl;    	    
            //pressure is subtracted
            //nfields = fields.num_elements();
            nfields = fielddef[0]->m_fields.size()-1;            
    }
    else
    {
    	    nfields = fielddef[0]->m_fields.size();         	    
    }
*/


    //fields.num_elements() is empty!!!
    //fielddef[0]->m_fields.size() can counts also pressure be careful!!!

    //subtract pressure for coupledsolver
    //the field has to be 3D homogeneous 1D to separate Re Im coeffs
    SetFields(graphShPt,boundaryConditions,vSession,fields,nfields);
    //----------------------------------------------      
    //In the first I need the boundary conditions
/*    
    MultiRegions::ContField2DSharedPtr firstfield;
    firstfield= MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(vSession, graphShPt,vSession->GetVariable(0),true);
    fields[0]=firstfield;
    //Set 2D ExpList for all the other variables:
    for(i = 1 ; i < nfields; i++)
    {                        	
             fields[i] = MemoryManager<MultiRegions::ContField2DHomogeneous1D>
                 ::AllocateSharedPtr(*firstfield,graphShPt,vSession->GetVariable(i),true);
    }    
*/  
    
 
    // Copy data from file:fill fields with the fielddata
    for(j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j]);
        }  
     
        if(j==3)
        {
	     Array<OneD, NekDouble> Recoeffs = fields[3]->GetPlane(0)->GetCoeffs();
	     Array<OneD, NekDouble> Imcoeffs = fields[3]->GetPlane(1)->GetCoeffs();
	     for(int t=0; t<Imcoeffs.num_elements(); t++)
	     {
cout<<"ncoeff="<<t<<"  A Re="<<Recoeffs[t]<<"  Im="<<Imcoeffs[t]<<endl;
	     }       	        	
        }        
        fields[j]->BwdTrans(fields[j]->GetCoeffs(),fields[j]->UpdatePhys());
        
       
    }
    //----------------------------------------------
                            
    //--------------------------------------------------------------    



    // determine the I regions
    //hypothesis: the number of I regions is the same for all the variables
    //hypothesis: all the I regions have the same nq points
    int nIregions, lastIregion,nq1D; 
    const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConditions  = fields[0]->GetBndConditions();      
    Array<OneD, int> Iregions =Array<OneD, int>(bndConditions.num_elements(),-1);    
    //2 internal layers:
    Array<OneD, int> Ilayers =Array<OneD, int>(2,-1);     
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



    //set output fields (dim=2 as the number of layers)
    int nlayers=2;
    Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldx= fields[0]->GetBndCondExpansions();   
    Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldy= fields[1]->GetBndCondExpansions();
    //Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldx (nlayers);
    //Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldy (nlayers);
    

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
    
    //TESTS'''''''
cout<<"nfields:"<<nfields<<endl;    
    //Array<OneD, MultiRegions::ExpListSharedPtr> bndpressure = fields[3]->GetBndCondExpansions();

    //-----------------------------------------------------------------      

    //set 1D output field:
    
    //manipulate data
   
    //for 2 variables(u,v) only:
    int coordim = graphShPt->GetMeshDimension();       
//    for(int s=0; s< Iregions.num_elements(); s++)
//    {
    	    
  	    
//           if(Iregions[s]!=-1)
//           { 
//cout<<"call manipulate region="<<s<<endl;         	   
       	        Manipulate(Ilayers,coordim, graphShPt, bcs, fields, bndfieldx,bndfieldy);       	       
//      	   }

//    }
    //--------------------------------------------------------------------------------------

/*
cout<<"streak"<<endl; 
    //add y to obtain the streak w=w'+y and write fld file           
    Array<OneD, MultiRegions::ExpListSharedPtr> streak;  
    int nvar=3;
    
   
    SetFields(graphShPt,boundaryConditions,vSession,streak,nvar);
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
      
    //write fld file for the streak
    string   file = fieldfile.substr(0,fieldfile.find_last_of("."));
    file= file+"_streak.fld";
    
 
    WriteFld(file,graphShPt,streak);
*/    	
    //---------------------------------------------------------------
  

    //write bcs files: one for each I region and each variable
    for(int s=0; s<nbnd; s++)
    {    	 
          if(Iregions[s]!=-1)
      	  {
      	        string var="u";
      	        WriteBcs(var,Iregions[s], fieldfile,graphShPt,bndfieldx[Iregions[s]]);
      	        var="v";
      	        WriteBcs(var,Iregions[s], fieldfile,graphShPt,bndfieldy[Iregions[s]]);      	      	      
      	  }
     }
    
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
	    Exp= Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);    
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

                        for(i = 0 ; i < Exp.num_elements(); i++)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous2D>
                                ::AllocateSharedPtr(session,BkeyY,BkeyZ,LhomY,LhomZ,useFFT,mesh,session->GetVariable(i));
                        }
                    }
                    else
                    {
                        for(i = 0 ; i < Exp.num_elements(); i++)
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
 cout<<"Znummodes="<<BkeyZ.GetNumModes()<<endl;
                        for(i = 0 ; i < Exp.num_elements(); i++)
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
                        for(i = 1 ; i < Exp.num_elements(); i++)
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
                            for(i = 1 ; i < Exp.num_elements(); i++)
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

       void  Manipulate(Array<OneD, int> Iregions,int coordim, SpatialDomains::MeshGraphSharedPtr &mesh,
       	       SpatialDomains::BoundaryConditions &bcs,
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &infields, 
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &outfieldx,
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &outfieldy)
        {
            //2 I regions are expected:
            ASSERTL0(Iregions.num_elements()==2, "something wrong with the number of I layers");
            int Iup =Iregions[0];
            int Idown=Iregions[1];
          
            //take the pressure from infields
            MultiRegions::ExpListSharedPtr pressure =infields[3];
            Array<OneD, MultiRegions::ExpListSharedPtr> Ipressure= infields[3]->GetBndCondExpansions();
            MultiRegions::ExpListSharedPtr pressureup = Ipressure[Iup];
            MultiRegions::ExpListSharedPtr pressuredown = Ipressure[Idown];
            
//cout<<"num elmt="<<infields[3]->GetNumElmts()<<endl;


            


        	
            int nq1D= pressureup->GetTotPoints();
            Array<OneD,NekDouble> x0(nq1D);
            Array<OneD,NekDouble> x1(nq1D);  
            Array<OneD,NekDouble> x2(nq1D); 
            pressureup->GetCoords(x0,x1,x2);               
cout<<"region="<<Iregions[0]<<" nq1D="<<nq1D<<endl;


    

/*
	     Array<OneD, NekDouble> Recoeffs = infields[3]->GetPlane(0)->GetCoeffs();
	     Array<OneD, NekDouble> Imcoeffs = infields[3]->GetPlane(1)->GetCoeffs();
	     for(int t=0; t<Imcoeffs.num_elements(); t++)
	     {
cout<<"ncoeff="<<t<<"  A Re="<<Recoeffs[t]<<"  Im="<<Imcoeffs[t]<<endl;
	     }       	        	
     
*/
	    








            Array<OneD, NekDouble> dtest=Array<OneD, NekDouble>(nq1D,0.0);

            //outfieldx[Iregions[0]]->PhysDeriv(MultiRegions::eS,outfieldx[Iregions[0]]->GetPhys(),dtest);            
            Array<OneD, NekDouble> tmp=Array<OneD, NekDouble>(nq1D, 0.0); 
            //Vmath::Vcopy(nq1D,dtest,1,(outfieldx[region]->UpdatePhys()),1);
            //Vmath::Vcopy(nq1D,tmp,1,(outfieldy[region]->UpdatePhys()),1); 


            const SpatialDomains::BoundaryRegionCollection &bregions
                = bcs.GetBoundaryRegions();
            SpatialDomains::Composite Icompup;
            SpatialDomains::Composite Icompdown;

            SpatialDomains::SegGeomSharedPtr segmentGeomup;
            SpatialDomains::SegGeomSharedPtr segmentGeomdown;
            SpatialDomains::ElementEdgeVectorSharedPtr elementup;
            SpatialDomains::ElementEdgeVectorSharedPtr elementdown;
            StdRegions::EdgeOrientation orientup;
            StdRegions::EdgeOrientation orientdown;            
            SpatialDomains::BoundaryRegion::iterator upIt, downIt;
            ASSERTL0(bregions[Iup]->size() ==
            	    bregions[Idown]->size(), "Size of the 2 internal layers should be the same");
            Array<OneD, unsigned int> bmapup, bmapdown;
            //fields to store the Re Im parts in physical space:
       	    Array<OneD,NekDouble> Rephysup (nq1D);
  	    Array<OneD,NekDouble> Imphysup (nq1D);            
           

            for(upIt = bregions[Iup]->begin(),
       	        downIt= bregions[Idown]->begin();
       	        upIt !=bregions[Iup]->end();
       	        ++upIt,++downIt)
       	    {
                Icompup = upIt->second;
                Icompdown=downIt->second;





		//assuming that ncoeffs is equal for all the edges (k=0)
		int nummodes= pressureup->GetExp(0)->GetNcoeffs() ;
		Array<OneD, int> offsetup(Icompup->size());
		offsetup[0]=0;


            	//array coeffs:  ATTENTION DOUBLE ENTRIES
cout<<"test"<<endl;            	
            	//Array<OneD, NekDouble> Recoeffsup =pressureup->GetPlane(0)->GetCoeffs();
            	//Array<OneD, NekDouble> Imcoeffsup =pressureup->GetPlane(1)->GetCoeffs();
            	Array<OneD, NekDouble> Recoeffsup (Icompup->size()*nummodes);
            	Array<OneD, NekDouble> Imcoeffsup (Icompup->size()*nummodes);


//cout<<"n edges uplayer="<<Icompup->size()<<endl;
		//define nq points for each edge
            	int nqedge =nq1D/Icompup->size();

cout<<"n edge /layer up ="<<Icompup->size()<<"  coeffs in each layer(Re)="<<Recoeffsup.num_elements()
				<<"  Ncoeffs(Im)="<<Imcoeffsup.num_elements()<<endl;            	
            	
		Array<OneD,NekDouble> x0edge(nqedge);
		Array<OneD,NekDouble> x1edge(nqedge);    



        	
                for(int k=0; k< Icompup->size(); k++)//loop over segments of each layer
                {
                    
                     if(!(segmentGeomup  = boost::dynamic_pointer_cast<
                                     SpatialDomains::SegGeom>((*Icompup)[k]))||
                             !(segmentGeomdown = boost::dynamic_pointer_cast<
                             	     SpatialDomains::SegGeom>((*Icompdown)[k]))
                     )
		     { ASSERTL0(false, "dynamic cast to a SegGeom failed");  }
		     
		     int EIDup = segmentGeomup->GetEid();
		     int EIDdown = segmentGeomdown->GetEid();
cout<<"E ids    up:"<<EIDup<<"  down:"<<EIDdown<<endl;		     
                     //set the offset(k*nummodes) of the 1D explist
                     offsetup[k+1] = offsetup[k]+nummodes;
cout<<"k="<<k<<"  offset="<<offsetup[k]<<"  nummodes="<<pressureup->GetExp(k)->GetNcoeffs()<<endl;                     

                     elementup = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(mesh)
                                    ->GetElementsFromEdge(segmentGeomup);
                     elementdown = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(mesh)
                                    ->GetElementsFromEdge(segmentGeomdown);
cout<<"elmt size="<<elementup->size()<<" down="<<elementdown->size()<<endl;                                    
	     	     int elmtidup= ((*elementup)[0]->m_Element)->GetGlobalID();   
                     int dimelmtup = ((*elementup)[0]->m_Element)->GetNumEdges();
                     
                     //get edges orientation:
                     orientup = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>
                     	           ((*elementup)[0]->m_Element))->GetEorient((*elementup)[0]
                     	           	   ->m_EdgeIndx);
                     orientdown =(boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>
                     	           ((*elementdown)[0]->m_Element))->GetEorient((*elementdown)[0]
                     	           	   ->m_EdgeIndx);
                     
                     
		     //	setup map between global and local ids
		     map<int, int> EdgeGIDup;
		     int id,cnt,f;
		     for(cnt=f=0; f<dimelmtup; f++)
		     {
			   id = (boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>
                     	           ((*elementdown)[0]->m_Element))->GetEid(f);
			   EdgeGIDup[id]=cnt++;		     	     
cout<<"local id="<<f<<"  Eid="<<id<<endl;			   
		     }
cout<<" GIDup="<<elmtidup<<"  num edges in this element="<<dimelmtup<<endl;
	             int localidup = EdgeGIDup.find(EIDup)->second;
	             localidup = (*elementup)[0]->m_EdgeIndx;
cout<<"check map local id of Eid="<<EIDup<<"  is="<<localidup<<endl;	             
                     int elmtiddown =((*elementdown)[0]->m_Element)->GetGlobalID();
cout<<" GIDdown="<<elmtiddown<<endl;
		     //set bmap 
		     Array<OneD, unsigned int > bmapedgeup;
		     Array<OneD, int > signup;
		     pressure->GetExp(elmtidup)->GetEdgeToElementMap(EdgeGIDup.find(EIDup)->second, orientup,bmapedgeup,signup);
		     		     
	//	     pressureup->GetExp(EIDup)->SetUpPhysTangents(locexp,localidup);


		     //Extracting the coeffs...
		     Array<OneD, NekDouble> Recoeffs = pressure->GetPlane(0)->GetCoeffs();
		     Array<OneD, NekDouble> Imcoeffs = pressure->GetPlane(1)->GetCoeffs();
cout<<" tot coeffs Re="<<Recoeffs.num_elements()<<"   Im="<<Imcoeffs.num_elements()<<endl;
        	     for(int j=0; j<Imcoeffs.num_elements(); j++)
        	     {
cout<<"j="<<j<<"  Re="<<Recoeffs[j]<<"  Im="<<Imcoeffs[j]<<endl;
		     }
                     //Array<OneD, NekDouble> Recoeffsedgeup =pressureup->GetPlane(0)->GetExp(k)->GetCoeffs();
                     //Array<OneD, NekDouble> Imcoeffsedgeup =pressureup->GetPlane(1)->GetExp(k)->GetCoeffs();
                     //TEST:::
/*                     
		     for(int g=0; g<Imcoeffs.num_elements(); g++)
		     {
cout<<"diff Im - Re="<<Imcoeffs[g]-Recoeffs[g]<<endl;
		     }
                     
*/                     
                     //Array<OneD, NekDouble> physup  =pressure->GetExp(elmtidup)->UpdatePhys();

                     //OFFSET EDGE REFERING TO ELEMENT EXPLIST MISSING		     
		     int Reoffset= pressure->GetPlane(0)->GetCoeff_Offset(elmtidup);
		     int Imoffset= pressure->GetPlane(1)->GetCoeff_Offset(elmtidup);
cout<<" Re offset="<<Reoffset<<"   Im offset="<<Imoffset<<endl;		     
		     //int Reoffsetedgeup = (pressure->GetPlane(0)->GetExp(elmtidup)->GetCoeffs()).num_elements();
//cout<<"num coeffs in this element="<<Reoffsetedgeup<<endl;		     
		     //offset = nummodes*k;
                     int ReNumoffset = 	pressureup->GetPlane(0)->GetNcoeffs();
                     int ImNumoffset =  pressureup->GetPlane(1)->GetNcoeffs();
//cout<<" OFFSET TEst Re="<<ReNumoffset<<"  Im="<<ImNumoffset<<endl;                     
		     int Reoffsetedgeup = pressureup->GetPlane(0)->GetCoeff_Offset(k);
		     
cout<<"k="<<k<<"  Offsetedge="<<Reoffsetedgeup<<endl;		     
		     
		     //hypothesis: 2 planes (0 Re, 1 Im)
		     //int nplanecoeffs = pressureup->GetNcoeffs();
		     int nplanecoeffs = pressure->GetExp(elmtidup)->GetNcoeffs();
		     //ATTENTION!!!! P HAS NUMMODES-2 > bmapedgeup-2????????
		     Array<OneD, NekDouble> Recoeffsedgeup(bmapedgeup.num_elements());	
		     Array<OneD, NekDouble> Imcoeffsedgeup(bmapedgeup.num_elements());
cout<<" Recoeffsedgeup num elements="<<Recoeffsedgeup.num_elements()<<"  n bmapedge="<<bmapedgeup.num_elements()<<endl;	

		     for(int d=0; d<bmapedgeup.num_elements(); d++)
		     {	 
cout<<"map value="<<bmapedgeup[d]<<endl;   	     	     	     
		     	   Recoeffsedgeup[d]=Recoeffs[Reoffset+bmapedgeup[d]];		     	   
		     	   Imcoeffsedgeup[d]=Imcoeffs[Imoffset+bmapedgeup[d]];		     	   		     	   
		     	   Recoeffsup[offsetup[k]+d] = Recoeffsedgeup[d];			     	   
		     	   Imcoeffsup[offsetup[k]+d] = Imcoeffsedgeup[d];
cout<<"Re coeff="<<Recoeffsup[offsetup[k]+d]<<"  Im coeff="<<Imcoeffsup[offsetup[k]+d]<<endl; 		     	   
		     }

		     
		     //pressureup->GetExp(k)->BwdTrans(Recoeffsedgeup,Rephysedgeup);		     
		     //pressure->GetExp(EIDup)->BwdTrans(Imcoeffsedgeup,Imphysedgeup);
		     
                     //
                     
	     
		     //Array<OneD,NekDouble> x2(nq1D); 
		     //(*elementup)[0]->GetCoords(x0edge,x1edge);		     
		     for(int w=0; w<nqedge;w++)
		     {
		   	// if(  	     
		     }

		     
		     
/*		     
		     //store the physical values:
		     for(int e=0; e<nq1D; e++)
		     {
		     }
*/		     
		     
		     
		     
	     
		     
		     
		     
		     
		     
/*	     
cout<<"extract tangent for edge="<<k<<endl;		
                     //TANGENTS WRONG WITH 3DHOMOGENEOUS 1D...
		     Array<OneD, Array<OneD, NekDouble> > tangents;
		     tangents = Array<OneD, Array<OneD, NekDouble> >(coordim);		     
		     for(int k=0; k<coordim; ++k)
		     {
              	          tangents[k]= Array<OneD, NekDouble>(nqedge); 
              	     }   		     
              	     //k is the number of the edge according to the composite list
		     LocalRegions::SegExpSharedPtr  bndSegExp = 
		     				boost::dynamic_pointer_cast<LocalRegions::SegExp>(Ipressure[Iup]->GetExp(k)); 		     	     

	             tangents = (bndSegExp)->GetMetricInfo()->GetEdgeTangent();
		     //calculate the curvature:
		     Array<OneD, NekDouble> curv=Array<OneD, NekDouble>(nqedge,0.0);
		     Array<OneD, NekDouble> dtx (nqedge);
		     Array<OneD, NekDouble> dty (nqedge);
		     outfieldx[Iregions[0]]->PhysDeriv(MultiRegions::eS,tangents[0],dtx);
		     outfieldx[Iregions[1]]->PhysDeriv(MultiRegions::eS,tangents[1],dty);
                     for(int t=0; t<nqedge; t++)
		     {
			  curv[t]=(dtx[t]**2 +dty[t]**2)**0.5;
	             }
           for(int w=0; w<nqedge; w++)
           {       
cout<<"tangent tx="<<tangents[0][w]<<" ty="<<tangents[1][w]<<endl;	
           }		     
*/		     

                }
cout<<"comp closed"<<endl;         
                //bwd transform:
                pressureup->BwdTrans(Recoeffsup,Rephysup);
                pressureup->BwdTrans(Imcoeffsup,Imphysup);
            }  
cout<<"dim Rephysup="<<Rephysup.num_elements()<<"  nq1D="<<nq1D<<endl;            


            //fill output fields:
            for(int g=0; g<nq1D; g++)
            {
            	              	    
            	(outfieldx[Iup]->UpdatePhys())[g] = 
            	   Rephysup[g];
            	(outfieldy[Iup]->UpdatePhys())[g] =
            	   Imphysup[g];		   
		
            }
cout<<"x"<<"   y"<<"    Re p"<<"    Im p"<<endl;                      	    
	    //print out the fields
	    for(int g=0; g<nq1D; g++)
	    {

cout<<(x0[g]*6.28)/1.6<<"    "<<x1[g]<<"    "<<Rephysup[g]<<"    "<<Imphysup[g]<<endl;
	    }
            
            
            
            for(int j=0; j<nq1D; j++)
            {
/*            	              	    
            	(outfieldx[region]->UpdatePhys())[j] = 
            	   20*sin(2*x0[j]);
            	(outfieldy[region]->UpdatePhys())[j] =
            	   20*2*cos(2*x0[j])/3.14159265;		   
*/		
            }
          

             
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
	
	void WriteBcs(string variable, int region, string fieldfile,SpatialDomains::MeshGraphSharedPtr &mesh,
		MultiRegions::ExpListSharedPtr &outregionfield)
	{			
		string   outfile = fieldfile.substr(0,fieldfile.find_last_of("."));
		outfile +="_"+variable+"_";
    		char ibnd[16]="";
    		sprintf(ibnd,"%d",region);
    		outfile +=ibnd;		
    		string   endfile(".bc");    		
    		outfile += endfile;
		std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                = outregionfield->GetFieldDefinitions();
                std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
                // copy Data into FieldData and set variable                
            	FieldDef[0]->m_fields.push_back(variable);
            	//fields->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs[j]);
            	outregionfield->AppendFieldData(FieldDef[0], FieldData[0]);
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
                    fields[i]->FwdTrans_IterPerExp(fields[i]->GetPhys(),fields[i]->UpdateCoeffs());
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
    
