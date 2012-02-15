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


    void SetFields(SpatialDomains::MeshGraphSharedPtr &graphShPt,
    	        vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables,
                bool homogeneous);
    void Readflddef(string fieldfile, 
                    Array<OneD, std::string> &variables, bool &homogeneous);

    void GenerateField(MultiRegions::ExpListSharedPtr field0,
    	               Array<OneD, NekDouble> x1,
    	               Array<OneD, NekDouble> y1,
    	               MultiRegions::ExpListSharedPtr field1);
    bool Checkbndmeshes(      Array<OneD, NekDouble> x0,
    	                      Array<OneD, NekDouble> y0,       	       
    	                      Array<OneD, NekDouble> x1,
    	                      Array<OneD, NekDouble> y1); 
    void Writefield(LibUtilities::SessionReaderSharedPtr vSession,
                    Array<OneD, std::string> &variables,
    	            string fieldfile, SpatialDomains::MeshGraphSharedPtr &graph,  	    
    	            Array<OneD, MultiRegions::ExpListSharedPtr> &outfield);    

    if(argc != 5)
    {
        fprintf(stderr,"Usage: ./FieldToField  meshfile0 fieldfile0  meshfile1  fieldfile1\n");
        exit(1);
    }


    //----------------------------------------------
    string meshfile0(argv[argc-4]); 
    string fieldfile0(argv[argc-3]);
    //argc=nfiles0 =2 
    //in this way only the mesh0 is taken to create vSession
    int nfiles0 = 2;
    //nfiles0=5;
    std::vector<std::string> filenames0;
    filenames0.push_back(meshfile0);
    filenames0.push_back(fieldfile0);
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames0);
            //= LibUtilities::SessionReader::CreateInstance(2, argv);
    // Read in mesh from input file0
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile0);
    //----------------------------------------------      
    
    // Import fieldfile0.

    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile0,fielddef,fielddata);
    //----------------------------------------------    

    //read info from fldfile
    Array<OneD, std::string> variables ;
cout<<fieldfile0<<endl;
    bool homo=false;
    Readflddef(fieldfile0, variables, homo);



    // Define Expansion    
    int nfields; 
    //nfields = vSession->GetVariables().size();       
    nfields = variables.num_elements();
    nfields=1;
    Array<OneD, MultiRegions::ExpListSharedPtr> fields; 
    fields= Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);  
    SetFields(graphShPt,fielddef, vSession,fields,nfields,homo);
    int nq = fields[0]->GetTotPoints();
    //-----------------------------------------------
cout<<"nfields="<<nfields<<endl;   
    // Copy data from file:fill fields with the fielddata
    for(int j = 0; j < nfields; ++j)
    {  	    
        for(int i = 0; i < fielddata.size(); ++i)
        {
cout<<"fielddef[i]->m_fields[j]="<<fielddef[i]->m_fields[j]<<endl;
            fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j]);
        }             
        fields[j]->BwdTrans_IterPerExp(fields[j]->GetCoeffs(),fields[j]->UpdatePhys());      
    }
    //----------------------------------------------    

    // store mesh0 quadrature points    
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> y0(nq);  
    fields[0]->GetCoords(x0,y0);    
  
/* test   
    for(int u=0; u<nq; u++)
    {
cout<<"x0="<<x0[u]<<"      y0="<<y0[u]<<endl;    	    
    }
    cout<<"finished"<<meshfile1<<endl; 
*/    

    //----------------------------------------------    
  
    //Read in the mesh1
    string meshfile1(argv[argc-2]);    
    //Name of the output fieldfile
    string fieldfile1(argv[argc-1]);           
    //-----------------------------------------------    
    
    //define the output field:
    Array<OneD, MultiRegions::ExpListSharedPtr> outfield; 
    outfield = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);   
    //set output fields over the graphShPt1 graph
    SpatialDomains::MeshGraphSharedPtr graphShPt1;     
    //remark: homo cases malloc() error or  segmentation fault..
    //remember: there is static cnt in stfields function to define 
    // the homo quantities only for the mesh0 case
/*
    if(session->DefinesSolverInfo("HOMOGENEOUS"))
    if(vSession->GetSolverInfo("HOMOGENEOUS")=="1D")
    {
    	    cout<<"homo case"<<endl;
    	    std::vector<std::string> filenames1;
    	    filenames1.push_back(meshfile1);
    	    //filenames1.push_back(fieldfile1);  
    	    argv[1] = argv[3];    	  
    	    const LibUtilities::CommSharedPtr comm = vSession->GetComm();    	    
    	    LibUtilities::SessionReaderSharedPtr vSession1
               // = LibUtilities::SessionReader::CreateInstance( 2 , &argv[3]);    
            = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames1, vSession->GetComm());    
            graphShPt1 = SpatialDomains::MeshGraph::Read(meshfile1);  	    
            SetFields(graphShPt1,fielddef, vSession1, outfield,nfields);
    }
    else
    {
    	    cout<<"2D case"<<endl;
*/    	    
            graphShPt1 = SpatialDomains::MeshGraph::Read(meshfile1);  	     	    
            SetFields(graphShPt1,fielddef, vSession, outfield,nfields,homo);    	    
//    }
    //------------------------------------------------ 


    // store the new points:
    int nq1 = outfield[0]->GetTotPoints();
    Array<OneD, NekDouble> x1(nq1);
    Array<OneD, NekDouble> y1(nq1);
    outfield[0]->GetCoords(x1,y1);
/*    
    for(int u=0; u<nq; u++)
    {
cout<<"x1="<<x1[u]<<"      y1="<<y1[u]<<endl;    	    
    }
*/    
    //------------------------------------------------
    
    //check 2Dmeshes compatibilities
    // same max min values x,y...
    bool check;    
    check = Checkbndmeshes(x0, y0, x1, y1);    
    ASSERTL0(check, "meshes not compatible (different borders)");   
    //-----------------------------------------------
    
    //generate the new fields
    for(int t=0; t< nfields; t++)
    {
    	 GenerateField(fields[t], x1, y1, outfield[t]);
    }
    //------------------------------------------------
/*
    //smooth the field
    if( graphShPt->GetMeshDimension()==2)
    {
         MultiRegions::ContField2DSharedPtr contfield;
           contfield = MemoryManager<MultiRegions::ContField2D>
           ::AllocateSharedPtr(vSession,graphShPt1,
                  vSession->GetVariable(0),true);
         int ncoeffs = outfield[0]->GetNcoeffs();
         Array<OneD, NekDouble> coeffs(ncoeffs);
         for(int s=0; s<nfields; s++)
         {
             contfield->FwdTrans(outfield[s]->GetPhys(), coeffs);
             contfield->BwdTrans(coeffs, outfield[s]->UpdatePhys());
         }
    }
*/
    //------------------------------------------------
    
    //write fieldfile
    Writefield(vSession, variables, fieldfile1, graphShPt1,outfield);        
    //------------------------------------------------
}


	// Define Expansion       		
	void SetFields(SpatialDomains::MeshGraphSharedPtr &graphShPt,
		vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables,
                bool homogeneous)
	{
                //session reader stuff has to be evaluated only from the
		// first session which refers to mesh0
                static int cnt=0;		
		// Setting parameteres for homogenous problems
        	MultiRegions::GlobalSysSolnType solnType;
		NekDouble static LhomX;           ///< physical length in X direction (if homogeneous) 
		NekDouble static LhomY;           ///< physical length in Y direction (if homogeneous)
		NekDouble static LhomZ;           ///< physical length in Z direction (if homogeneous)
		
		bool static DeclareCoeffPhysArrays = true;		
		int static npointsX;              ///< number of points in X direction (if homogeneous)
		int static npointsY;              ///< number of points in Y direction (if homogeneous)
                int static npointsZ;              ///< number of points in Z direction (if homogeneous)	
		int static HomoDirec       = 0;
		bool static useFFT = false;	        
		///Parameter for homogeneous expansions		
		enum HomogeneousType
		{
			eHomogeneous1D,
			eHomogeneous2D,
			eHomogeneous3D,
			eNotHomogeneous
		};
	
		enum HomogeneousType HomogeneousType = eNotHomogeneous;
cout<<"cnt="<<cnt<<" Homodir="<<HomoDirec<<endl;		
                if(cnt==0)
                { 	
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
		}		
		cnt++;
		int i;		
		int expdim   = graphShPt->GetMeshDimension();
		//Exp= Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);    
		// I can always have 3 variables in a 2D mesh (oech vel component i a function which can depend on 1-3 var)
		// Continuous Galerkin projection

        	switch(expdim)
        	{
        	     case 1:
        	     {
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

                               Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::
                                   AllocateSharedPtr(session,Bkey,ly,useFFT,graphShPt);
                               Exp[0] = Exp2DH1;

                               for(i = 1; i < nvariables; ++i)
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
                                const LibUtilities::PointsKey PkeyY(nylines+1,LibUtilities::ePolyEvenlySpaced);
                                const LibUtilities::BasisKey  BkeyY(fielddef[0]->m_basis[1],nylines,PkeyY);
				
                                const LibUtilities::PointsKey PkeyZ(nzlines+1,LibUtilities::ePolyEvenlySpaced);
                                const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);
                
                                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
                                NekDouble lz = fielddef[0]->m_homogeneousLengths[1];
				
                                Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::
                                      AllocateSharedPtr(session,BkeyY,BkeyZ,ly,lz,useFFT,graphShPt);
                                Exp[0] = Exp3DH2;
				
                                for(i = 1; i < nvariables; ++i)
                                {
                                      Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(*Exp3DH2);
                                }
                          }
                          else
                          {
                          	MultiRegions::ExpList1DSharedPtr Exp1D;
                          	Exp1D = MemoryManager<MultiRegions::ExpList1D>::
                          	         AllocateSharedPtr(session,graphShPt);
                          	         Exp[0] = Exp1D;
                          	for(i = 1; i < nvariables; ++i)
                          	{
                          	       Exp[i] = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(*Exp1D);
                          	}
                          }
                     }
                     break;
             	     case 2:
             	     {
cout<<"setfields"<<endl;             	     	     
                          ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

                          //if(fielddef[0]->m_numHomogeneousDir == 1)
                          if(homogeneous==true)
                          {
                          	MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                          	// Define Homogeneous expansion
                          	int nplanes = fielddef[0]->m_numModes[2];

                          	// choose points to be at evenly spaced points at
                          	// nplanes + 1 points
                          	const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                          	const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                          	NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                          	Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                          	     AllocateSharedPtr(session,Bkey,lz,useFFT,graphShPt,fielddef[0]->m_fields[0]);
                          	Exp[0] = Exp3DH1;
                          	for(i = 1; i < nvariables; ++i)
                          	{                         		
                          		Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                          		    AllocateSharedPtr(*Exp3DH1);
cout<<"set field="<<i<<endl;                           		    
                          	}
                          }
                          else
                          {
                                 MultiRegions::ExpList2DSharedPtr Exp2D;
                                 Exp2D = MemoryManager<MultiRegions::ExpList2D>::
                                   AllocateSharedPtr(session,graphShPt,true,fielddef[0]->m_fields[0]);
                                 Exp[0] =  Exp2D;

                                 for(i = 1; i < nvariables; ++i)
                                 {
                                       Exp[i] = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*Exp2D);
                                 }
                          }
                     }
                     break;
             	     case 3:
             	     {
                          MultiRegions::ExpList3DSharedPtr Exp3D;
                          Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(session,graphShPt);
                          Exp[0] =  Exp3D;

                          for(i = 1; i < nvariables; ++i)
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
        }

	void Readflddef(string fieldfile, Array<OneD, std::string> &variables,
                        bool &homogeneous)
        {
	    TiXmlDocument doc(fieldfile);
	    bool loadOkay = doc.LoadFile(); 
            TiXmlHandle docHandle(&doc);
            TiXmlElement* master = NULL;    // Master tag within which all data is contained.

            master = doc.FirstChildElement("NEKTAR");
            ASSERTL0(master, "Unable to find NEKTAR tag in file.");
            TiXmlElement* element = master->FirstChildElement("ELEMENTS");
            ASSERTL0(element, "Unable to find ELEMENTS tag within nektar tag.");
            //determine if the field is homogeneous
            const char *shape = element->Attribute("SHAPE");
            string homostr = string(shape);
            int pos_dash = homostr.find_first_of("-");
            if( pos_dash < homostr.length())
            {
                 homostr = homostr.substr(homostr.find_first_of("-"),homostr.length());

                 if(homostr=="HomogenousExp1D")
                 {
                      homogeneous = true;
                 }
            }           
            else
            {
                 homogeneous = false;
            }


            const char *vars = element->Attribute("FIELDS");   
            //convert char into string object 
            string varstr = string(vars);       
cout<<"char="<<vars<<endl;           
            if(varstr=="u" || varstr=="v" || varstr=="w" || varstr=="p")
            {
                 variables = Array<OneD, std::string>(1);
                 variables[0] = varstr;
cout<<variables[0]<<endl;
            }
            else if(varstr=="u,v" || varstr=="v,w" || varstr=="u,w")
            {
                 variables = Array<OneD, std::string>(2);
		 string   v0 = varstr.substr(0,varstr.find_first_of(","));
                 variables[0] = v0;
		 string   v1 = varstr.substr(varstr.find_first_of(","), varstr.length());
                 variables[1] = v1;
            }
            else if(varstr=="u,v,p")
            {
                 variables = Array<OneD, std::string>(3);
                 variables[0] = "u";
                 variables[1] = "v";
                 variables[2] = "p";
            }
        }

        void GenerateField(MultiRegions::ExpListSharedPtr field0,
    	                   Array<OneD, NekDouble> x1,
    	                   Array<OneD, NekDouble> y1,
    	                   MultiRegions::ExpListSharedPtr field1)
	{
             Array<OneD, NekDouble> coords(2);
	     int nq1 = field1->GetTotPoints();
             int elmtid, offset;
             for(int r=0; r< nq1; r++)
             {
                   coords[0] = x1[r];
                   coords[1] = y1[r];
                  
                   elmtid = field0->GetExpIndex(coords, 0.00001);
                   offset = field0->GetPhys_Offset(elmtid);
                   field1->UpdatePhys()[r] = field0->GetExp(elmtid)->
                           PhysEvaluate(coords, field0->GetPhys() +offset);    
                   if( isnan(field1->UpdatePhys()[r]) )
                   {            
cout<<"x="<<x1[r]<<"   y="<<y1[r]<<"    offset="<<offset<<"  elmtid="<<elmtid<<endl;                  
cout<<"new val="<<field1->UpdatePhys()[r]<<endl;
                       //ASSERTL0( abs(field1->UpdatePhys()[r])<10000000000, "interp failed");
                   }

             }        
	}	

       bool Checkbndmeshes(  Array<OneD, NekDouble> x0,
    	                      Array<OneD, NekDouble> y0,       	       
    	                      Array<OneD, NekDouble> x1,
    	                      Array<OneD, NekDouble> y1)
       {
       	       NekDouble x0min,x0max,y0min, y0max;
       	       NekDouble x1min,x1max,y1min, y1max;       	       
       	       NekDouble tol = 0.0000001;
       	       x0min = Vmath::Vmin(x0.num_elements(),x0,1);
       	       x0max = Vmath::Vmax(x0.num_elements(),x0,1);
       	       y0min = Vmath::Vmin(y0.num_elements(),y0,1);
       	       y0max = Vmath::Vmax(y0.num_elements(),y0,1);
       	       
       	       x1min = Vmath::Vmin(x1.num_elements(),x1,1);
       	       x1max = Vmath::Vmax(x1.num_elements(),x1,1);
       	       y1min = Vmath::Vmin(y1.num_elements(),y1,1);
       	       y1max = Vmath::Vmax(y1.num_elements(),y1,1);       	       

//cout<<"abs(x0min-x1min )="<<abs(x0min-x1min )<<endl;
//cout<<"abs(x0max-x1max)="<<abs(x0max-x1max)<<endl;
//cout<<"abs(y0min-y1min)="<<abs(y0min-y1min)<<endl;
//cout<<"abs(y0max-y1max)="<<abs(y0max-y1max)<<endl;

               if(  abs(x0min-x1min )< tol
                    && abs(x0max-x1max)< tol
               	    && abs(y0min-y1min)< tol
                    && abs(y0max-y1max)< tol
                 )
               {
               	   return true;
               }
               else
               { 
               	   return false;
               }
               	      

       }       	       


	void Writefield(LibUtilities::SessionReaderSharedPtr vSession,
                    Array<OneD, std::string> &variables,
    	            string fieldfile, SpatialDomains::MeshGraphSharedPtr &graph,  	    
    	            Array<OneD, MultiRegions::ExpListSharedPtr> &outfield)
    	{
    		string var;
    		std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
    			= outfield[0]->GetFieldDefinitions();  			
    		std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());    		
    		Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(outfield.num_elements());   	
		
                for(int j=0; j< fieldcoeffs.num_elements(); ++j)
		{  
		     outfield[j]->FwdTrans(outfield[j]->GetPhys(),outfield[j]->UpdateCoeffs()); 		     
		     fieldcoeffs[j] = outfield[j]->UpdateCoeffs();			
		     for(int i=0; i< FieldDef.size(); i++)
		     {		     	     
		     	   //var = vSession->GetVariable(j);		     	   		    
                           var =  variables[j];	   
		     	   FieldDef[i]->m_fields.push_back(var);   
		     	   outfield[0]->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs[j]);  
		     }
		}
		graph->Write(fieldfile,FieldDef,FieldData);		
	}    		
