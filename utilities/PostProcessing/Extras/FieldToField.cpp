#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>


#include <MultiRegions/ExpList.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>

#include <tinyxml/tinyxml.h>

#include <boost/math/special_functions/fpclassify.hpp>

using namespace Nektar;

//@todo : read the field info (i.e. if Homogeneous) from the
// fldfile instead of the sessionfile

int main(int argc, char *argv[])
{


    void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
    	        vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables,
                Array<OneD, std::string> variables, bool homogeneous);
    void Readflddef(string fieldfile, 
                    Array<OneD, std::string> &variables, bool &homogeneous);

    void GenerateField(MultiRegions::ExpListSharedPtr field0,
    	               Array<OneD, NekDouble> x1,
    	               Array<OneD, NekDouble> y1,
    	               Array<OneD, NekDouble> xvert,
    	               Array<OneD, NekDouble> yvert,
    	               MultiRegions::ExpListSharedPtr field1);
    void GenerateFieldHomo(MultiRegions::ExpListSharedPtr field0,
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
    //nfiles0=5;
    std::vector<std::string> filenames0;
    filenames0.push_back(meshfile0);
    filenames0.push_back(meshfile0);
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
    bool homo=true;
    Readflddef(fieldfile0, variables, homo);

    // Define Expansion    
    int nfields; 
    nfields = variables.num_elements();
    Array<OneD, MultiRegions::ExpListSharedPtr> fields; 
    fields= Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);  
    SetFields(graphShPt,fielddef, vSession,fields,nfields,variables,homo);
    int nq;//pointsper plane
    //-----------------------------------------------
    // Copy data from file:fill fields with the fielddata
    if(fielddef[0]->m_numHomogeneousDir == 1)
    {
        nq = fields[0]->GetPlane(1)->GetTotPoints();
        //THE IM PHYS VALUES ARE WRONG USING bwdTrans !!!
        for(int j = 0; j < nfields; ++j)
        {  
            for(int i = 0; i < fielddata.size(); ++i)
            {
                fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j], fields[j]->UpdateCoeffs());
            } 

            //bwd plane 0
            fields[j]->GetPlane(0)->BwdTrans_IterPerExp(fields[j]->GetPlane(0)->GetCoeffs(),  
                          fields[j]->GetPlane(0)->UpdatePhys() );
            
            //bwd plane 1
            fields[j]->GetPlane(1)->BwdTrans_IterPerExp(fields[j]->GetPlane(1)->GetCoeffs(), 
                          fields[j]->GetPlane(1)->UpdatePhys() );
            

        }    
    }
    else
    {
        nq = fields[0]->GetTotPoints();
        for(int j = 0; j < nfields; ++j)
        {  	    
            for(int i = 0; i < fielddata.size(); ++i)
            {
                fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j], fields[j]->UpdateCoeffs());
            }             
            fields[j]->BwdTrans_IterPerExp(fields[j]->GetCoeffs(),fields[j]->UpdatePhys());      
        }
    }


    // store mesh0 quadrature points    
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> y0(nq);  
    Array<OneD, NekDouble> z0(nq);

    if(
        //vSession->DefinesSolverInfo("HOMOGENEOUS")
        fielddef[0]->m_numHomogeneousDir == 1

      )
    {
       fields[0]->GetPlane(1)->GetCoords(x0,y0,z0);    
    }
    else
    {
       fields[0]->GetCoords(x0,y0);  
    }
  
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
  	    
    graphShPt1 = SpatialDomains::MeshGraph::Read(meshfile1);  


    //define second session1..	
    std::vector<std::string> filenames1;       
    filenames1.push_back(meshfile1);   
    filenames1.push_back(meshfile1);     
//cout<<"mesh1 read"<<endl;
    argc=2;
    argv[1]=argv[3];
    argv[2]=argv[4];
//cout<<argv[1]<<"  a="<<argv[2]<<endl;
    LibUtilities::SessionReaderSharedPtr vSession1
        = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames1, vSession->GetComm());
    //--------------------------------------------------

    //set explist 1

    SetFields(graphShPt1,fielddef, vSession1, outfield,nfields,variables,homo);    	    
    //------------------------------------------------ 
    // store the new points:
    int nq1;
    if(
        //vSession->DefinesSolverInfo("HOMOGENEOUS")
        fielddef[0]->m_numHomogeneousDir == 1

      )
    {
        nq1 = outfield[0]->GetPlane(1)->GetTotPoints();
    }
    else
    {
        nq1 = outfield[0]->GetTotPoints();
    }
    Array<OneD, NekDouble> x1(nq1);
    Array<OneD, NekDouble> y1(nq1);
    Array<OneD, NekDouble> z1(nq1);
    
    if(
        //vSession->DefinesSolverInfo("HOMOGENEOUS")
        fielddef[0]->m_numHomogeneousDir == 1

      )
    {
       outfield[0]->GetPlane(1)->GetCoords(x1,y1,z1);    
    }
    else
    {
       outfield[0]->GetCoords(x1,y1);
    }
    
    //------------------------------------------------
    
    //check 2Dmeshes compatibilities
    // same max min values x,y...
    bool check;    
    check = Checkbndmeshes(x0, y0, x1, y1);    
    ASSERTL0(check, "meshes not compatible (different borders)");   
    //-----------------------------------------------
    
    //generate the new fields
    int nVertTot = graphShPt->GetNvertices();
    //arrays to check with the new values
    Array<OneD, NekDouble> xvert(nVertTot);
    Array<OneD, NekDouble> yvert(nVertTot);
    for(int d=0; d< nVertTot; d++)
    {
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(d);
         NekDouble x,y,z;
         vertex->GetCoords(x,y,z); 
         xvert[d] =x;
         yvert[d] =y;         
    }



    if(
        //vSession->DefinesSolverInfo("HOMOGENEOUS")
        fielddef[0]->m_numHomogeneousDir == 1
      )
    {
        for(int t=0; t< nfields; t++)
        {        
    	    GenerateFieldHomo(fields[t], x1, y1, outfield[t]);
        }
    }
    else
    {
        for(int t=0; t< nfields; t++)
        {        
    	    GenerateField(fields[t], x1, y1,xvert,yvert, outfield[t]);
        }
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
        void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
    	        vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables,
                Array<OneD, std::string> variables, bool homogeneous)
	{
                //session reader stuff has to be evaluated only from the
		// first session which refers to mesh0
                static int cnt=0;		
		// Setting parameteres for homogenous problems
		NekDouble static LhomY;           ///< physical length in Y direction (if homogeneous)
		NekDouble static LhomZ;           ///< physical length in Z direction (if homogeneous)
		
		bool static DeclareCoeffPhysArrays = true;		
		int static npointsY;              ///< number of points in Y direction (if homogeneous)
                int static npointsZ;              ///< number of points in Z direction (if homogeneous)	
		bool static useFFT = false;
		bool static deal = false;
	
                if(cnt==0)
                { 	
                if(
                     //vSession->DefinesSolverInfo("HOMOGENEOUS")
                     fielddef[0]->m_numHomogeneousDir == 1
                  )
		{	
                        //only homo1D is working
                    /*
           		HomogeneousType = eHomogeneous1D;
                        npointsZ = fielddef[0]->m_numModes[2];
                        LhomZ = fielddef[0]->m_homogeneousLengths[0];
                    */
/*	
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
*/
	
		}
		}		
		cnt++;
		int i;		
		int expdim   = mesh->GetMeshDimension();
                //Exp= Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);  
		//Exp= Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);    
		// I can always have 3 variables in a 2D mesh (oech vel component i a function which can depend on 1-3 var)
		// Continuous Galerkin projection

        	switch(expdim)
        	{
                case 1:
                {
                    if(fielddef[0]->m_numHomogeneousDir == 1)
                    {
                        const LibUtilities::PointsKey PkeyY(npointsY,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyY(LibUtilities::eFourier,npointsY,PkeyY);

                        for(i = 0 ; i < nvariables; i++)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
                                ::AllocateSharedPtr(session,BkeyY,LhomY,useFFT,deal,mesh,variables[i]);
                        }
                    }
                    else
                    {
                        for(i = 0 ; i < nvariables; i++)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(session,mesh,variables[i]);
                        }
                    }

                    break;
                }
            case 2:
                {   
                    if(fielddef[0]->m_numHomogeneousDir == 1)
                    {
                        cout<<"homogeneous case"<<endl;
                        const LibUtilities::PointsKey PkeyZ(npointsZ,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,npointsZ,PkeyZ);

                        //patch to avoid the lack of p in coupledsolver sessionfile
                        if(variables[0]=="p")
                        {
                             variables[0]="u";
                        }

                        MultiRegions::ContField3DHomogeneous1DSharedPtr firstfield = 
                             MemoryManager<MultiRegions::ContField3DHomogeneous1D>
                          ::AllocateSharedPtr(session,BkeyZ,LhomZ,useFFT,deal,mesh,variables[0]); 
                        Exp[0] = firstfield;
                        for(i = 1 ; i < nvariables; i++)
                        {                        	
                            Exp[i]=MemoryManager<MultiRegions::
                                   ContField3DHomogeneous1D>::AllocateSharedPtr
                                     (*firstfield);                                   
                        }
                    }
                    else
                    {    
                        i = 0;
                        MultiRegions::ContField2DSharedPtr firstfield;
                        firstfield = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(session,mesh,variables[i],DeclareCoeffPhysArrays);

                        Exp[0] = firstfield;
                        for(i = 1 ; i < nvariables; i++)
                        {                        	
                            Exp[i] = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(*firstfield,mesh,variables[i],DeclareCoeffPhysArrays);
                        }
                    }

                    break;
                }
                case 3:
                    {
                        if(fielddef[0]->m_numHomogeneousDir == 1)
                        {
                            ASSERTL0(false,"3D fully periodic problems not implemented yet");
                        }
                        else
                        {
                            i = 0;
                            MultiRegions::ContField3DSharedPtr firstfield =
                                MemoryManager<MultiRegions::ContField3D>
                                ::AllocateSharedPtr(session,mesh,variables[i]);

                            Exp[0] = firstfield;
                            for(i = 1 ; i < nvariables; i++)
                            {
                                Exp[i] = MemoryManager<MultiRegions::ContField3D>
                                    ::AllocateSharedPtr(*firstfield,mesh,variables[i]);
                            }
                        }
                        break;
                    }
            default:
                ASSERTL0(false,"Expansion dimension not recognised");
                break;
            }   
  
        }

	void Readflddef(string fieldfile, Array<OneD, std::string> &variables,
                        bool &homogeneous)
        {
	    TiXmlDocument doc(fieldfile);
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
            cout<<"Fields to Interpolate: "<<vars<<endl;           
            if(varstr=="u" || varstr=="v" || varstr=="w" || varstr=="p")
            {
                 variables = Array<OneD, std::string>(1);
                 variables[0] = varstr;
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
            else if(varstr=="u,v,w,p")
            {
                 variables = Array<OneD, std::string>(4);
                 variables[0] = "u";
                 variables[1] = "v";
                 variables[2] = "w";
                 variables[3] = "p";
            }

        }

        void GenerateField(MultiRegions::ExpListSharedPtr field0,
    	                   Array<OneD, NekDouble> x1,
    	                   Array<OneD, NekDouble> y1,
    	                   Array<OneD, NekDouble> xvert,
    	                   Array<OneD, NekDouble> yvert,
    	                   MultiRegions::ExpListSharedPtr field1)
	{
             Array<OneD, NekDouble> coords(2);
	     int nq1 = field1->GetTotPoints();
             int elmtid, offset;
             for(int r=0; r< nq1; r++)
             {
                   coords[0] = x1[r];
                   coords[1] = y1[r];
                   elmtid = field0->GetExpIndex(coords, 0.001);
/*
                   if(elmtid <0)
                   {
                       NekDouble max;
          
                       //try another guess
                       Array<OneD, NekDouble> tmp(xvert.num_elements());
                       Array<OneD, NekDouble> dinstance(xvert.num_elements());
                       Array<OneD, NekDouble> tmpcoords(2);
                       int ind0=-1;
                       int ind1=-1;
                       int ind2=-1;
                       int ind3=-1;
                       //find neigh verts based on the dinstance:
                       Vmath::Vcopy(xvert.num_elements(), xvert,1,tmp,1);
                       Vmath::Sadd(xvert.num_elements(), -x1[r],tmp,1,tmp,1);
                       //tmp=(x-x0)^2
                       Vmath::Vmul(xvert.num_elements(), tmp,1,tmp,1,tmp,1);

                       Vmath::Vcopy(xvert.num_elements(), yvert,1,dinstance,1);
                       Vmath::Sadd(xvert.num_elements(), -y1[r],dinstance,1,dinstance,1);
                       //dinstance=(y-y0)^2
                       Vmath::Vmul(xvert.num_elements(), dinstance,1,dinstance,1,dinstance,1);
                       //dinstance=(x-x0)^2 +(y-y0)^2
                       Vmath::Vadd(xvert.num_elements(), tmp,1,dinstance,1,dinstance,1);
                       
                       max = Vmath::Vmax(xvert.num_elements(), dinstance,1);


                       ind0 = Vmath::Imin(xvert.num_elements(), dinstance,1); 
                       dinstance[ind0] = max+1000;
                       ind1 = Vmath::Imin(xvert.num_elements(), dinstance,1); 
                       dinstance[ind1] = max+1000;
                       ind2 = Vmath::Imin(xvert.num_elements(), dinstance,1); 
                       dinstance[ind2] = max+1000;
                       ind3 = Vmath::Imin(xvert.num_elements(), dinstance,1); 
                       dinstance[ind3] = max+1000;                       
                        
                     

                       elmtid = field0->GetExpIndex(tmpcoords, 0.001);

                   }
*/

                   if(elmtid<0)
                   {
cout<<"x="<<coords[0]<<"   y="<<coords[1]<<endl;
                        ASSERTL0(elmtid>=0, "elmtid not found"); 
                        //cout<<"warning: elmtid not found"<<endl;
                        //field1->UpdatePhys()[r]= field1->UpdatePhys()[r-1];
                        
                   }
                   else
                   {
                  
                   offset = field0->GetPhys_Offset(elmtid);

                   field1->UpdatePhys()[r] = field0->GetExp(elmtid)->
                           PhysEvaluate(coords, field0->GetPhys() +offset);    
                   }
                
//cout<<r<<"new val="<<field1->UpdatePhys()[r]<<endl;
                   if( (boost::math::isnan)(field1->UpdatePhys()[r]) )
                   {            
cout<<"x="<<x1[r]<<"   y="<<y1[r]<<"    offset="<<offset<<"  elmtid="<<elmtid<<endl;                  
cout<<"new val="<<field1->UpdatePhys()[r]<<endl;
                       ASSERTL0(false, "new value is not a number");
                       //ASSERTL0( abs(field1->UpdatePhys()[r])<10000000000, "interp failed");
                   }

             }        
	}	

        void GenerateFieldHomo(MultiRegions::ExpListSharedPtr field0,
    	               Array<OneD, NekDouble> x1,
    	               Array<OneD, NekDouble> y1,
    	               MultiRegions::ExpListSharedPtr field1)
        {
             Array<OneD, NekDouble> coords(2);
	     int nq1 = field1->GetPlane(1)->GetTotPoints();
             ASSERTL0(nq1 == field1->GetPlane(1)->GetTotPoints(), "problem");
             int elmtid, offset;

             //plane 0
             for(int r=0; r< nq1; r++)
             {
                   coords[0] = x1[r];
                   coords[1] = y1[r];
                  
                   elmtid = field0->GetPlane(0)->GetExpIndex(coords, 0.001);
                   offset = field0->GetPlane(0)->GetPhys_Offset(elmtid);
                   field1->GetPlane(0)->UpdatePhys()[r] = field0->GetPlane(0)->GetExp(elmtid)->
                           PhysEvaluate(coords, field0->GetPlane(0)->GetPhys() +offset);    
                   if( (boost::math::isnan)(field1->GetPlane(0)->UpdatePhys()[r]) )
                   {            
cout<<"x="<<x1[r]<<"   y="<<y1[r]<<"    offset="<<offset<<"  elmtid="<<elmtid<<endl;                  
cout<<"new val="<<field1->GetPlane(0)->UpdatePhys()[r]<<endl;
                       //ASSERTL0( abs(field1->UpdatePhys()[r])<10000000000, "interp failed");
                   }

              } 


             //plane1
             for(int r=0; r< nq1; r++)
             {
                   coords[0] = x1[r];
                   coords[1] = y1[r];
                  
                   elmtid = field0->GetPlane(1)->GetExpIndex(coords, 0.001);
                   offset = field0->GetPlane(1)->GetPhys_Offset(elmtid);
                   field1->GetPlane(1)->UpdatePhys()[r] = field0->GetPlane(1)->GetExp(elmtid)->
                           PhysEvaluate(coords, field0->GetPlane(1)->GetPhys() +offset);    
                   if( (boost::math::isnan)(field1->GetPlane(1)->UpdatePhys()[r]) )
                   {            
cout<<"x="<<x1[r]<<"   y="<<y1[r]<<"    offset="<<offset<<"  elmtid="<<elmtid<<endl;                  
cout<<"new val="<<field1->GetPlane(1)->UpdatePhys()[r]<<endl;
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
               NekDouble tol1 = 0.00001;
       	       x0min = Vmath::Vmin(x0.num_elements(),x0,1);
       	       x0max = Vmath::Vmax(x0.num_elements(),x0,1);
       	       y0min = Vmath::Vmin(y0.num_elements(),y0,1);
       	       y0max = Vmath::Vmax(y0.num_elements(),y0,1);
       	       
       	       x1min = Vmath::Vmin(x1.num_elements(),x1,1);
       	       x1max = Vmath::Vmax(x1.num_elements(),x1,1);
       	       y1min = Vmath::Vmin(y1.num_elements(),y1,1);
       	       y1max = Vmath::Vmax(y1.num_elements(),y1,1);       	       


int Ixmin1 = Vmath::Imin(x1.num_elements(),x1,1);
cout<<std::setprecision(8)<<"x0max="<<x0max<<endl;
cout<<std::setprecision(8)<<"x1max="<<x1max<<endl;
cout<<std::setprecision(8)<<"x0min="<<x0min<<endl;
cout<<std::setprecision(8)<<"x1min="<<x1min<<"   Ixmin1="<<Ixmin1<<"   y_xmin="<<y1[Ixmin1]<<endl;
cout<<std::setprecision(8)<<"y0max="<<y0max<<endl;
cout<<std::setprecision(8)<<"y1max="<<y1max<<endl;
cout<<"abs(x0min-x1min )="<<abs(x0min-x1min )<<endl;
cout<<"abs(x0max-x1max)="<<abs(x0max-x1max)<<endl;
cout<<"abs(y0min-y1min)="<<abs(y0min-y1min)<<endl;
cout<<"abs(y0max-y1max)="<<abs(y0max-y1max)<<endl;

               if(  abs(x0min-x1min )< tol1
                    && abs(x0max-x1max)< tol1
               	    && abs(y0min-y1min)< tol1
                    && abs(y0max-y1max)< tol1
                 )
               {
                   if(abs(x0min-x1min )> tol
                    || abs(x0max-x1max)> tol
               	    || abs(y0min-y1min)> tol
                    || abs(y0max-y1max)> tol  )
                   {
                        cout<<"Warning: mesh boundary points differ more than 10^-7"<<endl;
                   }

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
                        if(
                              //vSession->DefinesSolverInfo("HOMOGENEOUS")
                              FieldDef[0]->m_numHomogeneousDir == 1

                          )
                     {
                         //plane 0
                         outfield[j]->GetPlane(0)->FwdTrans_IterPerExp(outfield[j]->GetPlane(0)->GetPhys(),outfield[j]->GetPlane(0)->UpdateCoeffs());

                         //plane 1
                         outfield[j]->GetPlane(1)->FwdTrans_IterPerExp(outfield[j]->GetPlane(1)->GetPhys(),outfield[j]->GetPlane(1)->UpdateCoeffs());
                         
                     }
                     else
                     {  
		         outfield[j]->FwdTrans_IterPerExp(outfield[j]->GetPhys(),outfield[j]->UpdateCoeffs());
                     }
 		     
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
