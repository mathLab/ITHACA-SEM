#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
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
//#include </PreProcessing/MeshConvert/Convert.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
		SpatialDomains::BoundaryConditionsSharedPtr &boundaryConditions,
		LibUtilities::SessionReaderSharedPtr &session,
		Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,int nvariables);
    void OrderVertices(int nedges,SpatialDomains::MeshGraphSharedPtr graphShPt,
    	    MultiRegions::ExpListSharedPtr & bndfield, 
        	Array<OneD, int>& Vids, int v1,int v2 , NekDouble x_connect,
        	int & lastedge,
        	Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y);
    NekDouble yMove(NekDouble y, NekDouble yold_up, NekDouble Deltaold);
    void Replacevertices(string filename, Array<OneD, NekDouble> newx, 
    	               Array<OneD,  NekDouble> newy);
    
    
    
    int i,j;
    if(argc != 4)
    {
        fprintf(stderr,"Usage: ./FldCalcBCs  meshfile fieldfile  changefile\n");
        exit(1);
    }
//ATTEnTION !!! with argc=2 you impose that vSession refers to is argv[1]=meshfile!!!!! 
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(2, argv);
    //----------------------------------------------
   
    // Read in mesh from input file
    string meshfile(argv[argc-3]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);
    //---------------------------------------------- 

    // Also read and store the boundary conditions
    SpatialDomains::BoundaryConditionsSharedPtr boundaryConditions;        
    boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                                        ::AllocateSharedPtr(vSession,graphShPt);
    //----------------------------------------------

    // Define Expansion   
    Array<OneD, MultiRegions::ExpListSharedPtr> fields;   
    int nfields;  

    //the mesh file should have 2 component: set output fields
    //fields has to be of the SAME dimension of the mesh (that's why there is
    //the changefile as an input)
    nfields=2;          
    SetFields(graphShPt,boundaryConditions,vSession,fields,nfields);
    //---------------------------------------------------------------

    
    // store name of the file to change
    string changefile(argv[argc-1]);
    //----------------------------------------------
    
    // Import field file.
    string fieldfile(argv[argc-2]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------      
   
//cout<<"dim st="<<fielddef[0]->m_fields.size()<<endl;
    //fill a vector with the streak sol
    Array<OneD, MultiRegions::ExpListSharedPtr> streak;      
    SetFields(graphShPt, boundaryConditions, vSession, streak, 1);        
    for(int i = 0; i < fielddata.size(); ++i)
    {
         streak[0]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[0]);
    }    
    streak[0]->BwdTrans(streak[0]->GetCoeffs(), streak[0]->UpdatePhys());    
    //------------------------------------------------    
           

    //----------------------------------------------   
/*     
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
*/  
    // determine the I regions (2 regions expected)
    //hypothesis: the number of I regions is the same for all the variables
    //hypothesis: all the I regions have the same nq points
    int nIregions, lastIregion; 
    const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConditions  = fields[0]->GetBndConditions();    
    Array<OneD, int> Iregions =Array<OneD, int>(bndConditions.num_elements(),-1);    
   
    nIregions=0;
    int nbnd= bndConditions.num_elements();
    for(int r=0; r<nbnd; r++)
    {
    	  if(bndConditions[r]->GetUserDefined().GetEquation()=="I")
    	  {
    	  	  lastIregion=r;
    	  	  Iregions[r]=r;
    	  	  nIregions++;
    	  }    	  
    } 
    ASSERTL0(nIregions>0,"there is any boundary region with the tag USERDEFINEDTYPE=""I"" specified");
   
    //set expansion along a layers
    Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldx= fields[0]->GetBndCondExpansions();        
    //--------------------------------------------------------
     
    //determine the points in the lower and upper curve...
    int nq1D= bndfieldx[lastIregion]->GetTotPoints();
    int nedges = bndfieldx[lastIregion]->GetExpSize();
    int nvertl = nedges +1 ;
    Array<OneD, int> Vids_low(nvertl,-10);    
    Array<OneD, NekDouble> xold_low(nvertl);
    Array<OneD, NekDouble> yold_low(nvertl);    
    Array<OneD, NekDouble> zi(nvertl);

    //order the ids on the lower curve lastIregion starting from the id on x=0
    NekDouble x_connect;
    //first point for x_connect=0
    int lastedge=-1;
    int v1,v2;
    v1=0;
    v2=1;
    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion], 
       	Vids_low, v1, v2 , 0 ,lastedge, xold_low,yold_low);    
    SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_low[v2]);    
    NekDouble xt,yt,zt;
    //update x_connect    
    vertex->GetCoords(x_connect,yt,zt);
      
    i=2;
    while(i<nvertl)
    { 
         v1=i;
         OrderVertices(nedges, graphShPt, bndfieldx[lastIregion], 
        	Vids_low, v1, v2 , x_connect, lastedge, xold_low, yold_low );          
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_low[v1]);     
         //update x_connect  (lastedge is updated on the OrderVertices function) 
         vertex->GetCoords(x_connect,yt,zt);
         i++;           
    }
    //------------------------------------------------------------------------
    
    //order in the same way the id of the upper curve lastIregion-1 starting from x=0:
    Array<OneD, int> Vids_up(nvertl,-10);   
    Array<OneD,NekDouble> xold_up(nvertl);
    Array<OneD,NekDouble> yold_up(nvertl);     
    //first point for x_connect=0
    lastedge=-1;

    v1=0;
    v2=1;
    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-1 ], 
        	Vids_up, v1, v2 , 0 ,lastedge, xold_up, yold_up);    
    SpatialDomains::VertexComponentSharedPtr vertexU = graphShPt->GetVertex(Vids_up[v2]);    

    //update x_connect    
    vertexU->GetCoords(x_connect,yt,zt);
      
    i=2;
    while(i<nvertl)
    { 
         v1=i;
         OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-1], 
         	Vids_up, v1, v2 , x_connect, lastedge, xold_up, yold_up );          
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_up[v1]);  
   
         //update x_connect  (lastedge is updated on the OrderVertices function) 
         vertex->GetCoords(x_connect,yt,zt);
         i++;           
     }   
    //-----------------------------------------------------------------------------------
     
    //find the points where u=0 and calculate deltaold for each point
    int nq= fields[0]->GetTotPoints(); 
    Array<OneD, NekDouble> x(nq);
    Array<OneD,NekDouble> y(nq);    
    fields[0]->GetCoords(x,y);         
    Array<OneD, NekDouble> x_c(nvertl);
    Array<OneD,NekDouble> y_c(nvertl);        
    Array<OneD, NekDouble> Deltaold(nvertl,-200);
    for(int i=0; i<nvertl; i++)
    {
      //NekDouble u_min=100;      	    
      ASSERTL0(yold_up[i]>yold_low[i],"problem with the curve points");
      for(int j=0; j<nq; j++)
      {   
      	 //u is along y in this case..    NB:streak=u+y???
         if(xold_up[i]== x[j]  &&  (streak[0]->GetPhys()[j]+y[j]) < 0.001  && y[j]>=yold_low[i] && y[j]<=yold_up[i])
         {
              //if(   (streak[0]->GetPhys()[j]+y[j])< abs(u_min))
              //{
      	      y_c[i] = y[j];
      	      x_c[i] = x[j];
      	      Deltaold[i]= yold_up[i] - y_c[i];    
      	      //u_min= (streak[0]->GetPhys()[j]+y[j]);
      	      //}
      	      
         }         	          	 
      }
      if(Deltaold[i]==-200)
      {
             cout<<i<<"x="<<x[i]<<endl;      	      
             ASSERTL0(false, "error: cannot find a point where u=0 for this value of x=");

      }
cout<<"yold_up="<<yold_up[i]<<"  y_c="<<y_c[i]<<endl;   
cout<<"Deltaold="<<Deltaold[i]<<endl;         
    }
    //------------------------------------------------------------
    
    //calculate the new cordinates of the vertices
    int nVertTot = graphShPt->GetNvertices();
    Array<OneD, NekDouble> xnew(nVertTot);
    Array<OneD, NekDouble> ynew(nVertTot,-20);
    for( int i=0; i<nVertTot; i++)
    {
      bool mvpoint=false; 	    
      SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(i);
      NekDouble x,y,z;
      vertex->GetCoords(x,y,z); 
      //x coord doesn't change
      xnew[i]=x;
      //find the corresponding yold_l and deltaold
      for(int j=0; j<nvertl; j++)
      {	      
         if(xold_up[j]==x)
         {        	 
            ynew[i]=yMove(y,yold_up[j],Deltaold[j]);
            mvpoint=true;
         }
      }
      if(mvpoint==false)   
      {         
            //determine the closer xold_up
            NekDouble diff=1000;
            int qp_closer;
            for(int k=0; k<nvertl; k++)
            {
               if(abs(x-xold_up[k]) < diff)
               {
                   diff = abs(x-xold_up[k]);
                   qp_closer=k;
               }                                  
            }       
            ynew[i]=yMove(y,yold_up[qp_closer],Deltaold[qp_closer]);             	             
       }

    }
    for(int i=0; i<nVertTot; i++)
    {
cout<<"id="<<i<<"  x="<<xnew[i]<<"  y="<<ynew[i]<<endl;
    }
    
    //replace the vertices with the new ones
    Replacevertices(changefile, xnew , ynew);
    //Replacevertices(meshfile, xnew , ynew);
          	       
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
  
        void OrderVertices(int nedges, SpatialDomains::MeshGraphSharedPtr graphShPt,
        	MultiRegions::ExpListSharedPtr & bndfield, 
        	Array<OneD, int>& Vids, int v1,int v2, NekDouble x_connect, int & lastedge, 
        	Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y)
        {
            int nvertl = nedges+1;
            int edge;
            Array<OneD, int> Vids_temp(nvertl,-10);         	
            for(int j=0; j<nedges; j++)
            {
   	        LocalRegions::SegExpSharedPtr  bndSegExplow = 
   	             boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndfield->GetExp(j)) ;   	
   	        edge = (bndSegExplow->GetGeom1D())->GetEid();
   	   
   	        for(int k=0; k<2; k++)
   	        {
   	            Vids_temp[j+k]=(bndSegExplow->GetGeom1D())->GetVid(k);   
   	            SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_temp[j+k]);
   	            NekDouble x1,y1,z1;
   	            vertex->GetCoords(x1,y1,z1);     	   	   
   	            if(x1==x_connect && edge!=lastedge)
   	            {  
   	            	 //first 2 points   
			 if(x_connect==0)
			 {			 	 
			     Vids[v1]=Vids_temp[j+k]; 
			     x[v1]=x1;
			     y[v1]=y1;
			     
			     if(k==0 )
			     {
			     	   Vids_temp[j+1]=(bndSegExplow->GetGeom1D())->GetVid(1);     	            	 	 
			     	   Vids[v2]=Vids_temp[j+1]; 
			     	   SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids[v2]);
			     	   NekDouble x2,y2,z2;
			     	   vertex->GetCoords(x2,y2,z2); 
				   x[v2]=x2;
				   y[v2]=y2;				   			     	   
			     }
			     if(k==1)
			     {
			     	   Vids_temp[j+0]=(bndSegExplow->GetGeom1D())->GetVid(0);   	       	         	 
			     	   Vids[v2]=Vids_temp[j+0];
			     	   SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids[v2]);
			     	   NekDouble x2,y2,z2;
			     	   vertex->GetCoords(x2,y2,z2); 
				   x[v2]=x2;
				   y[v2]=y2;			     	   
			     } 	       	  

			 }
			 else //other vertices
			 {
			     if(k==0 )
			     {
			     	   Vids_temp[j+1]=(bndSegExplow->GetGeom1D())->GetVid(1);     	            	 	 
			     	   Vids[v1]=Vids_temp[j+1]; 
			     	   SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids[v1]);
			     	   NekDouble x1,y1,z1;
			     	   vertex->GetCoords(x1,y1,z1); 
				   x[v1]=x1;
				   y[v1]=y1;			     	   
			     }
			     if(k==1)
			     {
			     	   Vids_temp[j+0]=(bndSegExplow->GetGeom1D())->GetVid(0);   	       	         	 
			     	   Vids[v1]=Vids_temp[j+0];
			     	   SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids[v1]);
			     	   NekDouble x1,y1,z1;
			     	   vertex->GetCoords(x1,y1,z1); 
				   x[v1]=x1;
				   y[v1]=y1;			     	   
			     } 
			 }
			 break;			 
 	       	    }
 	       	} 
 	       	if(Vids[v1]!=-10)
 	       	{	
 	       	   lastedge = edge;
 	       	   break;
 	       	}
 	    }        	     	
        	
	}        	
        
        NekDouble yMove(NekDouble y, NekDouble yold_up, NekDouble Deltaold)
        {
            NekDouble ynew;        	
            //algorithm for yold_l>0	 
            if(yold_up>=0)
            {
            	    
                if(y>=0 && y<=1)
                {	
            	    ynew= y - Deltaold*(1-y)/(1-yold_up) ;
            	}
            	else if(y>=-1 && y<0)
                {
                    ynew= y - Deltaold*(1+y)/(1-yold_up);
                }
                else if(y>1 && y<-1)
                {                	
                    ASSERTL0(false, "y range wrong!!!");
                }               
            }
            //algorithm for yold_l<0
            else if(yold_up<0)
            {
                if(y>=0 && y<=1)
                {	
            	    ynew= y - Deltaold*(1-y)/(1+yold_up) ;
            	}
            	else if(y>=-1 && y<0)
                {
                    ynew= y - Deltaold*(1+y)/(1+yold_up);                     
                }
                else if(y>1 && y<-1)
                {                	
                    ASSERTL0(false, "y range wrong!!!");
                }              
            }
            else
            {
            	    ASSERTL0(false,"the upper curve corresponds to the critical layer");            	 
            }  
            return ynew;
	}        	
        
        
        
        void Replacevertices(string filename, Array<OneD, NekDouble> newx, Array<OneD, NekDouble> newy)
	{     
cout<<"OOOK"<<endl;		
	    //load existing file
            string newfile;
	    TiXmlDocument doc(filename); 
	    bool loadOkay = doc.LoadFile();
	    
            // Save a new XML file.	          
            newfile = filename.substr(0, filename.find_last_of("."))+"_moved.xml";
 
            doc.SaveFile( newfile );   
            
            //write the new vertices
	    TiXmlDocument docnew(newfile);  
	    bool loadOkaynew = docnew.LoadFile();

	    std::string errstr = "Unable to load file: ";
	    errstr += newfile;
	    ASSERTL0(loadOkaynew, errstr.c_str());

	    TiXmlHandle docHandlenew(&docnew);    
	    TiXmlNode* nodenew = NULL;
	    TiXmlElement* meshnew = NULL;
	    TiXmlElement* masternew = NULL;    // Master tag within which all data is contained.

      
	    masternew = docnew.FirstChildElement("NEKTAR");
	    ASSERTL0(masternew, "Unable to find NEKTAR tag in file.");

	    // Find the Mesh tag and same the dim and space attributes
	    meshnew = masternew->FirstChildElement("GEOMETRY");

	    ASSERTL0(meshnew, "Unable to find GEOMETRY tag in file.");
	    TiXmlAttribute *attrnew = meshnew->FirstAttribute();
	    // Now read the vertices
	    TiXmlElement* elementnew = meshnew->FirstChildElement("VERTEX");
	    ASSERTL0(elementnew, "Unable to find mesh VERTEX tag in file.");
	    TiXmlElement *vertexnew = elementnew->FirstChildElement("V");
 
      	    
	    int indx;
	    int err;
	    int nextVertexNumber = -1;	
	    
	    while (vertexnew)
	    {
	    	   nextVertexNumber++;
	    	   //delete the old one
	    	   TiXmlAttribute *vertexAttr = vertexnew->FirstAttribute();	    	 
	    	   std::string attrName(vertexAttr->Name());	
	    	   ASSERTL0(attrName == "ID", (std::string("Unknown attribute name: ") + attrName).c_str());

	    	   err = vertexAttr->QueryIntValue(&indx);
	    	   ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");	    	 
	    	   ASSERTL0(indx == nextVertexNumber, "Vertex IDs must begin with zero and be sequential.");

	    	   std::string vertexBodyStr;   
	    	   // Now read body of vertex
       	           TiXmlNode *vertexBody = vertexnew->FirstChild();
       	           // Accumulate all non-comment body data.
       	           if (vertexBody->Type() == TiXmlNode::TEXT)
       	           {
       	                vertexBodyStr += vertexBody->ToText()->Value();
       	                vertexBodyStr += " ";
       	           }       
       	           ASSERTL0(!vertexBodyStr.empty(), "Vertex definitions must contain vertex data.");	
       	           vertexnew->RemoveChild(vertexBody);
                   //write the new one   
cout<<"writing.. v:"<<nextVertexNumber<<endl;	    	    	
	    	   stringstream s;
	    	   s << std::scientific << std::setprecision(3) <<  newx[nextVertexNumber] << "   "
	    	   << newy[nextVertexNumber] << "   " << 0.0;
	    	   vertexnew->LinkEndChild(new TiXmlText(s.str()));      
	    	   //TiXmlNode *newvertexBody = vertexnew->FirstChild();
	    	   //string newvertexbodystr= newvertexBody->SetValue(s.str());                     
	    	   //vertexnew->ReplaceChild(vertexBody,new TiXmlText(newvertexbodystr));
	    	 
	    	   vertexnew = vertexnew->NextSiblingElement("V");  
       	   }	
       	    //meshnew->LinkEndChild(elementnew);

   	    
       	    docnew.SaveFile( newfile ); 
       	    
       	    cout<<"new file:  "<<newfile<<endl;
	   
	}
