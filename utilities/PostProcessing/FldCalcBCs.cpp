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
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>


using namespace Nektar;

int main(int argc, char *argv[])
{
    void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
		SpatialDomains::BoundaryConditionsSharedPtr &boundaryConditions,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables,
		LibUtilities::CommSharedPtr comm);   	    
    void Manipulate(int region, int coordim, Array<OneD,MultiRegions::ExpListSharedPtr> &infields,  
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
    //----------------------------------------------
 
    // Define Expansion   
    Array<OneD, MultiRegions::ExpListSharedPtr> fields;   
    int nfields; 
/*   
    if(boundaryConditions->GetSolverInfo("SOLVERTYPE")=="CoupledLinearisedNS")
    {
            //pressure is subtracted
            //nfields = fields.num_elements();
            nfields = fielddef[0]->m_fields.size();            
    }
    else
    {
    	    //nfields = fielddef[0]->m_fields.size();   
            nfields = fielddef[0]->m_fields.size()-1;       	    
    }
*/
    //fields.num_elements() is empty!!!
    //fielddef[0]->m_fields.size() can counts also pressure be careful!!!
    nfields= fielddef[0]->m_fields.size();
cout<<"nfields"<<nfields<<endl;     
    SetFields(graphShPt,boundaryConditions,fields,nfields,vComm);
    //----------------------------------------------   
     
    // Copy data from file:fill fields with the fielddata
    for(j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j]);
        }
        fields[j]->BwdTrans(fields[j]->GetCoeffs(),fields[j]->UpdatePhys());
    }
    //----------------------------------------------    
    
    // determine the I regions
    //hypothesis: the number of I regions is the same for all the variables
    //hypothesis: all the I regions have the same nq points
    int nIregions, lastIregion,nq1D; 
    const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConditions  = fields[0]->GetBndConditions();      
    Array<OneD, int> Iregions =Array<OneD, int>(bndConditions.num_elements(),-1);    
  
    nIregions=0;
    int nbnd= bndConditions.num_elements();
  
    for(int r=0; r<nbnd; r++)
    {
    	  if(bndConditions[r]->GetUserDefined().GetEquation()=="CalcBC")
    	  {
    	  	  lastIregion=r;
    	  	  Iregions[r]=r;
    	  	  nIregions++;
    	  }    	  
    } 
    ASSERTL0(nIregions>0,"there is any boundary region with the tag USERDEFINEDTYPE=""CalcBC"" specified");
    
    //set output fields
    Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldx= fields[0]->GetBndCondExpansions();   
    Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldy= fields[1]->GetBndCondExpansions(); 
    
    //--------------------------------------------------------
    
    //manipulate data 
    //for 2 variables(u,v) only:
    int coordim = graphShPt->GetMeshDimension();       
    for(int s=0; s< Iregions.num_elements(); s++)
    {
  	    
           if(Iregions[s]!=-1)
           {             			
       	        Manipulate(Iregions[s],coordim, fields, bndfieldx,bndfieldy);       	       
      	   }

    }
    //--------------------------------------------------------------------------------------
  
/*
    //add y to obtain the streak u=u'+y and write fld file
    int nq=fields[0]->GetTotPoints();
    Array<OneD, NekDouble> x(nq,0.0);
    Array<OneD, NekDouble> y(nq,0.0);
    fields[0]->GetCoords(x,y); 
    for(int i=0; i<nq; i++)
    {
    	(fields[1]->UpdatePhys())[i]=(fields[1]->GetPhys())[i]+y[i];   	
    }
      
      //write fld file for the streak
    string   file = fieldfile.substr(0,fieldfile.find_last_of("."));
    file= file+"_streak.fld";
    WriteFld(file,graphShPt,fields);
    	
    //---------------------------------------------------------------
*/    

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
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables,
		LibUtilities::CommSharedPtr comm)
	{
			
		// Setting parameteres for homogenous problems
        	MultiRegions::GlobalSysSolnType solnType;
		NekDouble LhomX;           ///< physical length in X direction (if homogeneous) 
		NekDouble LhomY;           ///< physical length in Y direction (if homogeneous)
		NekDouble LhomZ;           ///< physical length in Z direction (if homogeneous)
		
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
		
		if(boundaryConditions->SolverInfoExists("HOMOGENEOUS"))
		{
			std::string HomoStr = boundaryConditions->GetSolverInfo("HOMOGENEOUS");
			//m_spacedim          = 3;
			
			if((HomoStr == "HOMOGENEOUS1D")||(HomoStr == "Homogeneous1D")||
			   (HomoStr == "1D")||(HomoStr == "Homo1D"))
			{
				HomogeneousType = eHomogeneous1D;
				npointsZ        = boundaryConditions->GetParameter("HomModesZ");
				LhomZ           = boundaryConditions->GetParameter("LZ");
				HomoDirec       = 1;
			}
			
			if((HomoStr == "HOMOGENEOUS2D")||(HomoStr == "Homogeneous2D")||
			   (HomoStr == "2D")||(HomoStr == "Homo2D"))
			{
				HomogeneousType = eHomogeneous2D;
				npointsY        = boundaryConditions->GetParameter("HomModesY");
				LhomY           = boundaryConditions->GetParameter("LY");
				npointsZ        = boundaryConditions->GetParameter("HomModesZ");
				LhomZ           = boundaryConditions->GetParameter("LZ");
				HomoDirec       = 2;
			}
			
			if((HomoStr == "HOMOGENEOUS3D")||(HomoStr == "Homogeneous3D")||
			   (HomoStr == "3D")||(HomoStr == "Homo3D"))
			{
				HomogeneousType = eHomogeneous3D;
				npointsX        = boundaryConditions->GetParameter("HomModesX");
				LhomX           = boundaryConditions->GetParameter("LX");
				npointsY        = boundaryConditions->GetParameter("HomModesY");
				LhomY           = boundaryConditions->GetParameter("LY");
				npointsZ        = boundaryConditions->GetParameter("HomModesZ");
				LhomZ           = boundaryConditions->GetParameter("LZ");
				HomoDirec       = 3;
			}
			
			if(boundaryConditions->SolverInfoExists("USEFFT"))
			{
				useFFT = true;
			}
		}		

		
	    int i;		
	    int expdim   = mesh->GetMeshDimension();
	    Exp= Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);	     
            switch(expdim)
            {
                case 1:
                {
					if(HomogeneousType == eHomogeneous2D)
					{
						SpatialDomains::MeshGraph1DSharedPtr mesh1D;
						
						if(!(mesh1D = boost::dynamic_pointer_cast<
							 SpatialDomains::MeshGraph1D>(mesh)))
						{
							ASSERTL0(false,"Dynamics cast failed");
						}
						
						const LibUtilities::PointsKey PkeyY(npointsY,LibUtilities::eFourierEvenlySpaced);
						const LibUtilities::BasisKey  BkeyY(LibUtilities::eFourier,npointsY,PkeyY);
						const LibUtilities::PointsKey PkeyZ(npointsZ,LibUtilities::eFourierEvenlySpaced);
						const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,npointsZ,PkeyZ);
						
						for(i = 0 ; i < Exp.num_elements(); i++)
						{
							Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous2D>
							::AllocateSharedPtr(comm,BkeyY,BkeyZ,LhomY,LhomZ,useFFT,*mesh1D,*boundaryConditions,i);
						}
					}
					else 
					{
						SpatialDomains::MeshGraph1DSharedPtr mesh1D;
						
						if( !(mesh1D = boost::dynamic_pointer_cast<
							  SpatialDomains::MeshGraph1D>(mesh)) )
						{
							ASSERTL0(false,"Dynamics cast failed");
						}
						
						for(i = 0 ; i < Exp.num_elements(); i++)
						{
							Exp[i] = MemoryManager<MultiRegions::ContField1D>
							::AllocateSharedPtr(comm,*mesh1D,
                                                *boundaryConditions,i);
						}
						
					}
					
                    break;
                }
                case 2:
                {                	
                    if(HomogeneousType == eHomogeneous1D)
                    {
                        SpatialDomains::MeshGraph2DSharedPtr mesh2D;
			
                        if(!(mesh2D = boost::dynamic_pointer_cast<
                             SpatialDomains::MeshGraph2D>(mesh)))
                        {
                            ASSERTL0(false,"Dynamics cast failed");
                        }
			
                        const LibUtilities::PointsKey PkeyZ(npointsZ,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,npointsZ,PkeyZ);
                        
                        for(i = 0 ; i < Exp.num_elements(); i++)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
                                ::AllocateSharedPtr(comm,BkeyZ,LhomZ,useFFT,*mesh2D,*boundaryConditions,i);
                        }
                    }
                    else
                    {                	    
                        SpatialDomains::MeshGraph2DSharedPtr mesh2D;
			
                        if(!(mesh2D = boost::dynamic_pointer_cast<
                             SpatialDomains::MeshGraph2D>(mesh)))
                        {
                            ASSERTL0(false,"Dynamics cast failed");
                        }
/*
			i=0;
			MultiRegions::ContField2DSharedPtr firstfield =
			     MemoryManager<MultiRegions::ContField2D>
			     ::AllocateSharedPtr(comm,*mesh2D,*boundaryConditions,i);
			Exp[0] = firstfield;
*/						
                        for(i = 0 ; i < Exp.num_elements(); i++)
                        {                          	
                            Exp[i] = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(comm,
                                                    *mesh2D,*boundaryConditions,i);
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
						SpatialDomains::MeshGraph3DSharedPtr mesh3D;

						if(!(mesh3D = boost::dynamic_pointer_cast<
										SpatialDomains::MeshGraph3D>(mesh)))
						{
							ASSERTL0(false,"Dynamics cast failed");
						}
/*						
						i=0;
						MultiRegions::ContField3DSharedPtr firstfield =
						MemoryManager<MultiRegions::ContField3D>
						::AllocateSharedPtr(comm,*mesh3D,*boundaryConditions,i);
						Exp[0] = firstfield;
*/						
						for(i = 0 ; i < Exp.num_elements(); i++)
						{
							Exp[i] = MemoryManager<MultiRegions::ContField3D>
											::AllocateSharedPtr(comm,
												*mesh3D,*boundaryConditions,i);
						}
					}
                    break;
                }
            default:
                ASSERTL0(false,"Expansion dimension not recognised");
                break;
            }   
           
        }

       void  Manipulate(int region,int coordim, Array<OneD,MultiRegions::ExpListSharedPtr> &infields, 
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &outfieldx,
       	       Array<OneD,MultiRegions::ExpListSharedPtr> &outfieldy)
        {        	
            int nq1D= outfieldx[region]->GetTotPoints();
            Array<OneD, NekDouble> dtest=Array<OneD, NekDouble>(nq1D,0.0);
            outfieldx[region]->PhysDeriv(MultiRegions::eS,outfieldx[region]->GetPhys(),dtest);            
            Array<OneD, NekDouble> tmp=Array<OneD, NekDouble>(nq1D, 0.0); 
            //Vmath::Vcopy(nq1D,dtest,1,(outfieldx[region]->UpdatePhys()),1);
            //Vmath::Vcopy(nq1D,tmp,1,(outfieldy[region]->UpdatePhys()),1); 
            
            
            
            Array<OneD,NekDouble> x0(nq1D);
            Array<OneD,NekDouble> x1(nq1D);  
            outfieldx[region]->GetCoords(x0,x1);            
            //u,v equation:
cout<<"region="<<region<<" nq1D="<<nq1D<<endl;            
            for(int j=0; j<nq1D; j++)
            {
/*            	              	    
            	(outfieldx[region]->UpdatePhys())[j] = 
            	   20*sin(2*x0[j]);
            	(outfieldy[region]->UpdatePhys())[j] =
            	   20*2*cos(2*x0[j])/3.14159265;		   
*/		
            }

            
            Array<OneD, Array<OneD, NekDouble> > tangents;
            tangents = Array<OneD, Array<OneD, NekDouble> >(coordim);

            for(int k=0; k<coordim; ++k)
            {
              	    tangents[k]= Array<OneD, NekDouble>(nq1D); 
            }   
            
            // need the tangents related to the expList1D outfieldx[region]

           LocalRegions::SegExpSharedPtr  bndSegExp =  boost::dynamic_pointer_cast<LocalRegions::SegExp>(outfieldx[region]->GetExp(0)); 
           tangents = (bndSegExp)->GetMetricInfo()->GetEdgeTangent();
cout<<"tangent x="<<tangents[0][0]<<" y="<<tangents[1][0]<<endl;
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
            else
            {
            	    ASSERTL0(false, " a 2D fld solution of a SteadyAdvectionDiffusion equation  is expected ");
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
    