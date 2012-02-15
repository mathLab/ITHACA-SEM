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
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <boost/lexical_cast.hpp>

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
    void Computestreakpositions(int nvertl, MultiRegions::ExpListSharedPtr &streak,
    	        Array<OneD, NekDouble> &x,  Array<OneD, NekDouble> &y,
    	        Array<OneD, NekDouble> &xold_up, Array<OneD, NekDouble> &yold_up,
    	        Array<OneD, NekDouble> &xold_low, Array<OneD, NekDouble> &yold_low,
    	        Array<OneD, NekDouble> &xold_c, Array<OneD, NekDouble> &yold_c,      	        
    	        Array<OneD, NekDouble> &xc,  Array<OneD, NekDouble> &yc); 
    void GenerateAddPointsNewtonIt( NekDouble &xi, NekDouble &yi,NekDouble &x0, NekDouble &y0,
    	        MultiRegions::ExpListSharedPtr &function, Array<OneD, NekDouble> &derfunction);
    void GenerateCurve(int npoints, int npused, Array<OneD, NekDouble> &x_c, 
    	    Array<OneD, NekDouble> &y_c, Array<OneD, NekDouble> &curve);
    void GenerateAddPoints(int region, SpatialDomains::MeshGraphSharedPtr &mesh,int np, int npused,    	      
    	         Array<OneD, NekDouble> &curve, MultiRegions::ExpListSharedPtr & bndfield, 
    	         Array<OneD, NekDouble>& outx, Array<OneD, NekDouble>& outy,
    	         Array<OneD, int>&Eids);
    NekDouble yMove(NekDouble y, NekDouble yold_up, NekDouble Deltaold);   
    void Replacevertices(string filename, Array<OneD, NekDouble> newx, 
    	               Array<OneD,  NekDouble> newy,
    	               Array<OneD, NekDouble> x_lay, Array<OneD, NekDouble> y_lay,
    	               Array<OneD, NekDouble> & Pcurvx, 
    	               Array<OneD, NekDouble> & Pcurvy,
    	               Array<OneD, int>&Eids, int Npoints);
    
    
    
    int i,j;
    if(argc != 4)
    {
        fprintf(stderr,"Usage: ./MoveMesh  meshfile fieldfile  changefile\n");
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
    nfields=3;          
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
/*    
    Array<OneD, MultiRegions::ExpListSharedPtr> streak;      
    SetFields(graphShPt, boundaryConditions, vSession, streak, 1); 
    for(int i = 0; i < fielddata.size(); ++i)
    {
         streak[0]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[0]);
    }
    streak[0]->BwdTrans(streak[0]->GetCoeffs(), streak[0]->UpdatePhys());       
*/    
    MultiRegions::ExpListSharedPtr streak; 
    streak = MemoryManager<MultiRegions::ContField2D>
          ::AllocateSharedPtr(vSession, graphShPt, "w",true);    

    for(int i=0; i<fielddata.size(); i++)
    {
    	    streak->ExtractDataToCoeffs(fielddef[i], fielddata[i], fielddef[i]->m_fields[0]);
    }
    streak->BwdTrans(streak->GetCoeffs(), streak->UpdatePhys());
    
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
    // determine the I regions (3 region expected)
    //hypothesys: the layes have the same number of points

    int nIregions, lastIregion; 
    const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConditions  = fields[0]->GetBndConditions();    
    Array<OneD, int> Iregions =Array<OneD, int>(bndConditions.num_elements(),-1);    
   
    nIregions=0;
    int nbnd= bndConditions.num_elements();
    for(int r=0; r<nbnd; r++)
    {
    	  if(bndConditions[r]->GetUserDefined()==SpatialDomains::eCalcBC)
    	  {
    	  	  lastIregion=r;
    	  	  Iregions[r]=r;
    	  	  nIregions++;
    	  }    	  
    } 
    ASSERTL0(nIregions>0,"there is any boundary region with the tag USERDEFINEDTYPE=""CalcBC"" specified");
cout<<"nIregions="<<nIregions<<endl;   
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
    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-1], 
       	Vids_low, v1, v2 , 0 ,lastedge, xold_low,yold_low);    
    SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_low[v2]);    
    NekDouble xt,yt,zt;
    //update x_connect    
    vertex->GetCoords(x_connect,yt,zt);
      
    i=2;
    while(i<nvertl)
    {     	    
         v1=i;
         OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-1], 
        	Vids_low, v1, v2 , x_connect, lastedge, xold_low, yold_low );          
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_low[v1]);     
//cout<<"Vids low="<<Vids_low[v1]<<endl;         
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
    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-2 ], 
        	Vids_up, v1, v2 , 0 ,lastedge, xold_up, yold_up);    
    SpatialDomains::VertexComponentSharedPtr vertexU = graphShPt->GetVertex(Vids_up[v2]);    
//cout<<"VIdup="<<Vids_up[v2]<<endl;
    //update x_connect    
    vertexU->GetCoords(x_connect,yt,zt);
      
    i=2;
    while(i<nvertl)
    { 
         v1=i;
         OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-2], 
         	Vids_up, v1, v2 , x_connect, lastedge, xold_up, yold_up );          
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_up[v1]);  
//cout<<"VIdup="<<Vids_up[v1]<<endl;   
         //update x_connect  (lastedge is updated on the OrderVertices function) 
         vertex->GetCoords(x_connect,yt,zt);
         i++;           
     }   
    
    //-----------------------------------------------------------------------------------
    
    
    //order in the same way the id of the layer curve lastIregion starting from x=0:
    Array<OneD, int> Vids_c(nvertl,-10);   
    Array<OneD,NekDouble> xold_c(nvertl);
    Array<OneD,NekDouble> yold_c(nvertl);     
    //first point for x_connect=0
    lastedge=-1;

    v1=0;
    v2=1;
    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion], 
        	Vids_c, v1, v2 , 0 ,lastedge, xold_c, yold_c);    
    SpatialDomains::VertexComponentSharedPtr vertexc = graphShPt->GetVertex(Vids_c[v2]);    

    //update x_connect    
    vertexc->GetCoords(x_connect,yt,zt);
      
    i=2;
    while(i<nvertl)
    { 
         v1=i;
         OrderVertices(nedges, graphShPt, bndfieldx[lastIregion], 
         	Vids_c, v1, v2 , x_connect, lastedge, xold_c, yold_c );          
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_c[v1]);  
   cout<<"Vids cl="<<Vids_low[v1]<<endl;  
         //update x_connect  (lastedge is updated on the OrderVertices function) 
         vertex->GetCoords(x_connect,yt,zt);
         i++;           
     }   
     
     //calculate the distances between the layer and the upper/lower curve
     
     Array<OneD, NekDouble> Deltaup (nvertl, -200);
     Array<OneD, NekDouble> Deltalow (nvertl, -200);    

     for(int r=0; r<nvertl; r++)
     {
     	     //Always positive!!!
//cout<<"i="<<r<<"  yup="<<yold_up[r]<<"   yc="<<yold_c[r]<<"   ylow="<<yold_low[r]<<endl;     	     
     	     Deltaup[r] = yold_up[r] - yold_c[r];
             Deltalow[r] = yold_c[r] - yold_low[r];           
             ASSERTL0(Deltaup[r]>0, "distance between upper and layer curve is not positive");
             ASSERTL0(Deltalow[r]>0, "distance between lower and layer curve is not positive");
     }           
     //------------------------------------------------------------------------
         
    //-----------------------------------------------------------------------------------    
     
    //find the points where u=0 and determine the sign of the shift and the delta
    int nq= fields[0]->GetTotPoints(); 
    Array<OneD, NekDouble> x(nq);
    Array<OneD,NekDouble> y(nq);    
    fields[0]->GetCoords(x,y);         
    Array<OneD, NekDouble> x_c(nvertl);
    Array<OneD,NekDouble> y_c(nvertl,-200);       
    Array<OneD, NekDouble> tmp_w (nvertl, 200);
    Array<OneD, int> Sign (nvertl,1);   
    Array<OneD, NekDouble> Delta_c(nvertl,-200);
   
    //calculate the dU_dy 
    Array<OneD, NekDouble> dU(nq);
    streak->PhysDeriv(MultiRegions::eY, streak->GetPhys(), dU);
    Computestreakpositions(nvertl, streak, x,y,xold_up, yold_up,
    	                   xold_low, yold_low, xold_c, yold_c, x_c, y_c);    
    // if the curve is low the old layer point, it has to shift down    
    for(int q=0; q<nvertl; q++)
    {
         if(y_c[q] < yold_c[q])
         {
             Sign[q] = -1;
         }
         //calculate delta
         Delta_c[q] = abs(yold_c[q]-y_c[q]);
    }
    //------------------------------------------------------------------
    

    
    int ncurves=1;
    //array storing ALL the curvecoeffs 
    //hypothesis: totcoeffs< 2*nvertl; ncurves< nedges;
    Array<OneD, NekDouble> totcurvecoeffs (2*nvertl,0.0);
    Array<OneD, int> ncurvepoints (nedges);
    int ntotcurvecoeffs=0;
    //additional points arrays
    Array<OneD, NekDouble> Cpointsx (nedges);
    Array<OneD, NekDouble> Cpointsy (nedges, 0.0);
    Array<OneD, int> Eids (nedges);    

/*
    //calculate the approximated curve coeffs
    //since the number of degree of freedom is too large 
    // we split the layer in as many curves as many changing of the rate 
    // \Delta( (\Delta y /\Delta x) ) (second derivative)
    //for each change of rate, calculate one curve and the additional point for
    // the curved tag of the session file    
    
    //change of curve second derivative sign quantities
    NekDouble der,tmp;
    tmp = (y_c[1]-y_c[0])/(x_c[1]-x_c[0]) ;
    der = (y_c[2]-y_c[1])/(x_c[2]-x_c[1]) ;
    int sign2der=1;    
    int Npoints=0;
    int Npused=0;
    int cnt=0;
    if(der < tmp)
    {
      sign2der =-1;
    }
//cout<<"sign="<<sign2der<<"  tmp="<<tmp<<"   der="<<der<<endl;    
    for(int r=2; r< nedges; r++)
    {
        tmp = (y_c[r+1]-y_c[r])/(x_c[r+1]-x_c[r]);

        if(sign2der==1)
        { 	
           if(tmp < der)
           {
              cnt++;
              //Array<OneD, NekDouble> curvecoeffs (nvertl);               
//cout<<"change 1"<<endl;   
              Npoints = r-(ncurves-1)-Npoints;
	      ncurvepoints[ncurves-1] = Npoints;              
              ncurves++;
              sign2der = -1*sign2der;

//cout<<"Npoints="<<Npoints<<"  r="<<r<<endl;   
              if(Npused==0)
              {
                 ntotcurvecoeffs += Npoints;              	      
                 GenerateCurve(Npoints, Npused, x_c, y_c, totcurvecoeffs);
                 GenerateAddPoints(lastIregion, graphShPt, Npoints, Npused, totcurvecoeffs,  
                 	     bndfieldx[lastIregion],    Cpointsx, Cpointsy, Eids);
                 
              }
              else
              {
                 ntotcurvecoeffs += Npoints+1;              	      
              	 GenerateCurve(Npoints+1, Npused-1, x_c, y_c, totcurvecoeffs);   
                 GenerateAddPoints(lastIregion, graphShPt, Npoints+1, Npused-1, totcurvecoeffs,  
                 	     bndfieldx[lastIregion],    Cpointsx, Cpointsy, Eids);              	 
              }
              Npused +=Npoints;              
           }
        }
        else
        {
           if(tmp > der)
           {
              cnt++;
//cout<<"change -1"<<endl;           	   
              Npoints = r-(ncurves-1)-Npoints;
	      ncurvepoints[ncurves-1] = Npoints;              
              ncurves++;
              sign2der = -1*sign2der;
//cout<<"npoints="<<Npoints<<"  r="<<r<<endl;              
              if(Npused==0)
              {
                 ntotcurvecoeffs += Npoints;              	      
                 GenerateCurve(Npoints, Npused, x_c, y_c, totcurvecoeffs);
                 GenerateAddPoints(lastIregion, graphShPt, Npoints, Npused, totcurvecoeffs,  
                 	     bndfieldx[lastIregion],    Cpointsx, Cpointsy, Eids);
              }
              else
              {
                 ntotcurvecoeffs += Npoints+1;              	      
              	 GenerateCurve(Npoints+1, Npused-1, x_c, y_c, totcurvecoeffs);     
                 GenerateAddPoints(lastIregion, graphShPt, Npoints+1, Npused-1, totcurvecoeffs,  
                 	     bndfieldx[lastIregion],    Cpointsx, Cpointsy, Eids);              	 
              }               
              Npused +=Npoints;       
           }
        }
//cout<<"sign="<<sign2der<<"  tmp="<<tmp<<"   der="<<der<<endl;    
        
        der =tmp;        
    }
    
    // generate last curve adding a point to fill the gap
    //between the previus calculated curve    
    int Nplast = nvertl - Npused;
    if(Npused==0)
    {
         ntotcurvecoeffs +=  Nplast;   
         ncurvepoints[ncurves-1] = Nplast; 
         GenerateCurve(Npoints, Npused, x_c, y_c, totcurvecoeffs);
         GenerateAddPoints(lastIregion, graphShPt, Npoints, Npused, totcurvecoeffs,  
                  bndfieldx[lastIregion],    Cpointsx, Cpointsy, Eids);
    }
    else
    {
         ntotcurvecoeffs +=  Nplast+1;    
//cout<<"last curve   npoints="<<Nplast<<"   Npused="<<Npused<<endl; 
//cout<<"ntotcurvecoeffs="<<ntotcurvecoeffs<<endl;
      	 GenerateCurve(Nplast+1, Npused-1, x_c, y_c, totcurvecoeffs); 
         GenerateAddPoints(lastIregion, graphShPt, Nplast+1, Npused-1, totcurvecoeffs,  
                  bndfieldx[lastIregion],    Cpointsx, Cpointsy, Eids);       	 
    }  
//cout<<" final curves coeffs values"<<endl;
//cout<<"ntotcurvecoeffs="<<ntotcurvecoeffs<<endl;
    ASSERTL0( ntotcurvecoeffs <= 2*nvertl, "numer of curves coeffs exceeds 2*nvertices");
    for(int u=0; u <2*nvertl; u++)
    {
//cout<<"coeffs = "<<totcurvecoeffs[u]<<endl;
    }

*/    
    
    //generate additional points using the Newton iteration
    //determine the xposition for every edge (in the middle even if it 
    // is not necessary
    //PARAMETER which determines the number of points @todo put as an input
    int npedge=5;
    //additional points arrays
/*
    //NB Double array LEAK!!!!
    Array<OneD, Array<OneD, NekDouble> >Addpointsx;
    Array<OneD, Array<OneD, NekDouble> >Addpointsy;
    Addpointsx = Array<OneD, Array<OneD, NekDouble> > (nedges);
    Addpointsy = Array<OneD, Array<OneD, NekDouble> > (nedges);
    for(int q=0; q< npedge-2; q++)
    {
        // npedge-2 num points per edge
        Addpointsx[q] = Array<OneD, NekDouble> (npedge-2, 0.0);
        Addpointsy[q] = Array<OneD, NekDouble> (npedge-2, 0.0);       
    }
    
*/    
    Array<OneD, NekDouble> Addpointsx (nedges*(npedge-2), 0.0);
    Array<OneD, NekDouble> Addpointsy (nedges*(npedge-2), 0.0);    
    //Array<OneD, NekDouble> Ycoords (nedges*npedge-2, 0.0);


    Array<OneD, NekDouble> derstreak (streak->GetTotPoints());
    streak->PhysDeriv(MultiRegions::eY, streak->GetPhys(), derstreak);
    
    LocalRegions::SegExpSharedPtr bndSegExp;
            
    int Eid,id1,id2;
    NekDouble x1,y1,z1;
    NekDouble x2,y2,z2;
    SpatialDomains::VertexComponentSharedPtr vertex1;
    SpatialDomains::VertexComponentSharedPtr vertex2;    
    for(int r=0; r<nedges; r++)
    {
    	    
         bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndfieldx[lastIregion]->GetExp(r));   
         Eid = (bndSegExp->GetGeom1D())->GetEid();
         id1 = (bndSegExp->GetGeom1D())->GetVid(0);
         id2 = (bndSegExp->GetGeom1D())->GetVid(1);  
         vertex1 = graphShPt->GetVertex(id1);      
	 vertex2 = graphShPt->GetVertex(id2);     
	 vertex1->GetCoords(x1,y1,z1);
	 vertex2->GetCoords(x2,y2,z2);	
//cout<<"edge="<<r<<"  x1="<<x1<<"  x2="<<x2<<endl;
//cout<<"edge="<<r<<"  y1="<<y1<<"  y2="<<y2<<endl;
	 if(x2>x1)
	 {
	     Cpointsx[r] = x1 +(x2-x1)/2;             
//cout<<"edge="<<r<<"  x1="<<x1<<"  x2="<<x2<<"   Cx="<<Cpointsx[r]<<endl;
//cout<<"edge="<<r<<"  y1="<<y1<<"  y2="<<y2<<endl;
	     if( Cpointsx[r]>x2 || Cpointsx[r]< x1)
	     {
	     	  Cpointsx[r] = -Cpointsx[r];		  
	     }    
             for(int w=0; w< npedge-2; w++)
             {    
                 Addpointsx[r*(npedge-2) +w] = x1 +((x2-x1)/(npedge - 1))*(w+1);   
                 if( Addpointsx[r*(npedge-2) +w] > x2 || Addpointsx[r*(npedge-2) +w] < x1)
	         {
	              Addpointsx[r*(npedge-2) +w] = -Addpointsx[r*(npedge-2) +w];
	         }                   	 
	         Addpointsy[r*(npedge-2) +w] = y1 + ((y2-y1)/(x2-x1))*(Addpointsx[r*(npedge-2) +w]-x1);	         
	         GenerateAddPointsNewtonIt( Addpointsx[r*(npedge-2) +w], Addpointsy[r*(npedge-2) +w],
	              	      Addpointsx[r*(npedge-2) +w],  Addpointsy[r*(npedge-2) +w], streak, derstreak); 
             }

	 }
	 else if(x1>x2)
	 {	         	 
	     Cpointsx[r] = x2+ (x1-x2)/2;
	     if( Cpointsx[r] > x1 || Cpointsx[r] < x2)
	     {
	          Cpointsx[r] = -Cpointsx[r];
	     }
             for(int w=0; w< npedge-2; w++)
             { 
                 Addpointsx[r*(npedge-2) +w] = x2 +((x1-x2)/(npedge - 1))*(w+1);
                 if( Addpointsx[r*(npedge-2) +w] > x1 || Addpointsx[r*(npedge-2) +w] < x2)
	         {
	              Addpointsx[r*(npedge-2) +w] = -Addpointsx[r*(npedge-2) +w];	              
	         }
	         Addpointsy[r*(npedge-2) +w] = y2 + ((y1-y2)/(x1-x2))*(Addpointsx[r*(npedge-2) +w]-x2);	         	         
	         GenerateAddPointsNewtonIt( Addpointsx[r*(npedge-2) +w], Addpointsy[r*(npedge-2) +w], 
	               Addpointsx[r*(npedge-2) +w], Addpointsy[r*(npedge-2) +w], streak, derstreak); 
             }	      
	 }
	 else
	 {
	      ASSERTL0(false, "point not generated"); 	 
	 }    	    
//cout<<"calculate cpoints coords"<<endl;	 
         Cpointsy[r] = y1 + (y2-y1)/2;
         GenerateAddPointsNewtonIt( Cpointsx[r], Cpointsy[r],Cpointsx[r], Cpointsy[r],
    	       streak, derstreak); 
         NekDouble diff = Cpointsy[r]-Addpointsy[r*(npedge-2)];
//cout<<"diff="<<diff<<endl;         
	 Eids[r] = Eid;

    }      

    
    //generate the closest curve to the critical layer positions taking 
    // npedge-2 points per edge
    Array<OneD, NekDouble> Crlay_pointsx(nedges*(npedge-2)+nedges+1);
    Array<OneD, NekDouble> Crlay_pointsy(nedges*(npedge-2)+nedges+1);
    int cnt=0;
    int cnt1=0;
    int cnt2=0;
    //HYPOTHESIS: Cpoints ALREADY generated
cout<<"Cr x--y"<<endl;    
    for(int u=0; u<Crlay_pointsx.num_elements(); u++)
    {
    	  if( u== cnt)// u pari  
          {
               Crlay_pointsx[u] = x_c[cnt1];
               Crlay_pointsy[u] = y_c[cnt1];
               cnt =cnt +npedge-1;
               cnt1++;
          }
          else//u dispari
          {
    	       //Crlay_pointsx[u] = Cpointsx[cnt1];
    	       //Crlay_pointsy[u] = Cpointsy[cnt1];
    	       Crlay_pointsx[u] = Addpointsx[cnt2];
    	       Crlay_pointsy[u] = Addpointsy[cnt2];    	       
    	       cnt2++;
          }
//cout<<u<<"      "<<Crlay_pointsx[u]<<"       "<<Crlay_pointsy[u]<<endl;          
    }
//cout<<"num px="<<Crlay_pointsx.num_elements()<<endl;   
/* 
    for(int r=0; r< Crlay_pointsx.num_elements(); r++)
    {
cout<<r<<"     "<<Crlay_pointsx[r]<<"                "<<Crlay_pointsy[r]
<<"       "<<sqrt((Crlay_pointsx[r+1]- Crlay_pointsx[r])
*(Crlay_pointsx[r+1]- Crlay_pointsx[r]) +
(Crlay_pointsy[r+1]- Crlay_pointsy[r])
*(Crlay_pointsy[r+1]- Crlay_pointsy[r])) <<endl;
    }
    
*/   
    
    
    Array<OneD, NekDouble> curve_coeffs (nedges*(npedge-2)+nedges+1);
    //generate the curve
    //NB:: for a high number of points (>20)  it does not work for all 
    // points
//cout<<"CALL generate curve"<<nedges*(npedge-2)+nedges+1<<endl;    
    //GenerateCurve(nedges*(npedge-2)+nedges+1, 0, Crlay_pointsx, Crlay_pointsy, curve_coeffs);    
    //calculate the coords of the additional points as y=poly(x)
    //HYPOTHESES: 
    //x coords ALREADY CALCULATED
    //npedge is the number of points per edge
    //Eid[s] ALREADY GENERATED
/*    
    int polorder; 
    double xl; 
    //put zero Addpointsy!!!!!
    Vmath::Zero(nedges*(npedge-2), Addpointsy,1);    
    for(int e=0; e< nedges*(npedge-2); e++)              
    {
    	  	 xl = Addpointsx[e];
    	  	 polorder = curve_coeffs.num_elements()-1; 	  

                 for(int g= 0; g< curve_coeffs.num_elements(); g++)
	         {	         	 
	             Addpointsy[e] += curve_coeffs[g]*(std::pow( xl , polorder));
//cout<<g<<"  coeff="<<curve_coeffs[g]<<"   x^exp="<<(std::pow( xl , polorder))<<"  exp="<<polorder<<endl;
                     ASSERTL0(polorder >=0, " polynomial with one negative exponent");
                     polorder--;
	         }
cout<<xl<<"     "<<Addpointsy[e]<<"      "<<Crlay_pointsx[e]<<"       "
<<Crlay_pointsy[e]<<"     "
<<curve_coeffs[e]<<endl;	         

    }
cout<<"LAST coeff==y(x==0)="<<curve_coeffs[nedges*(npedge-2)+nedges] <<endl;
*/
    //------------------------------------------------------------    
  
    //------------------------------------------------------------
    
    //calculate the new cordinates of the vertices    
    int nVertTot = graphShPt->GetNvertices();
    Array<OneD, NekDouble> xnew(nVertTot);
    Array<OneD, NekDouble> ynew(nVertTot,-20);
    for(int i=0; i<nVertTot; i++)
    {
         bool mvpoint =false;
       	 NekDouble ratio;  
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(i);
         NekDouble x,y,z;
         vertex->GetCoords(x,y,z); 
         //x coord doesn't change
         xnew[i]=x;
         //find the corresponding yold_l and deltaold
         for(int j=0; j<nvertl; j++)
         {	      
             if((xold_up[j]==x)&&(yold_up[j]==y))
             {        	 
                 //ynew[i]=y_c[j] + Sign[j]*Deltaup[j];
                 ratio = (1-y)*(1+y)/( (1-yold_c[j])*(1+yold_c[j]) );
                 ynew[i] = y + Sign[j]*Delta_c[j]*ratio;
                 //ynew[i] = y_c[j] +(y-yold_c[j])*ratio;
                 mvpoint=true;
             }
             if((xold_low[j]==x)&&(yold_low[j]==y))
             {
             	 //ynew[i]= y_c[j] + Sign[j]*Deltalow[j];
                 ratio = (1-y)*(1+y)/( (1-yold_c[j])*(1+yold_c[j]) );
                 ynew[i] = y + Sign[j]*Delta_c[j]*ratio;
                 //ynew[i] = y_c[j] +(y-yold_c[j])*ratio;             	 
             	 mvpoint=true;
             }
             if((xold_c[j]==x)&&(yold_c[j]==y))
             {
             	ynew[i] = y_c[j];
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
             if( y>yold_up[qp_closer] )
             {	        
                //ynew[i] = y + Sign[qp_closer]*Deltaup[qp_closer];
                 ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                 ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;
//cout<<"upper zone y="<<y<<"  ratio="<<ratio<<endl;                    
                 //ynew[i] = y_c[qp_closer] +(y-yold_c[qp_closer])*ratio;                 
             }
             else if(y<yold_low[qp_closer])
             {
             	//ynew[i] = y + Sign[qp_closer]*Deltalow[qp_closer];
                 ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                 ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;              
                 //ynew[i] = y_c[qp_closer] +(y-yold_c[qp_closer])*ratio;                  
             }
             else if ( y>yold_c[qp_closer] && y < yold_up[qp_closer])
             {
                //ynew[i] = y_c[qp_closer] +Sign[qp_closer]*abs(yold_up[qp_closer]-y);
                 ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                 ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;
                 //ynew[i] = y_c[qp_closer] +(y-yold_c[qp_closer])*ratio;                  
             }
             else if (y<yold_c[qp_closer] && y > yold_low[qp_closer])
             {
             	//ynew[i] = y_c[qp_closer] +Sign[qp_closer]*abs(yold_low[qp_closer]-y);
                 ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                 ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;
                 //ynew[i] = y_c[qp_closer] +(y-yold_c[qp_closer])*ratio;                  
             }
             else
             {
             	cout<<"point x="<<xnew[i]<<"  y="<<y<<"  closer x="<<xold_up[qp_closer]<<endl;   
             	ASSERTL0(false, "impossible to shift the point");
             }
             if(ynew[i]>1 || ynew[i]<-1)
             {
             	cout<<"point x="<<xnew[i]<<"  y="<<y<<"  closer x="<<xold_up[qp_closer]<<endl;              	     
             	ASSERTL0(false, "shifting out of range");
             }
         }         
//cout<<"x="<<x<<"  y="<<y<<"  ynew="<<ynew[i]<<endl;          
    }
  
    //------------------------------------------------------------------
    
    //replace the vertices with the new ones
    //Replacevertices(changefile, xnew , ynew, x_c, y_c, Cpointsx, Cpointsy, Eids, npedge);
    Replacevertices(changefile, xnew , ynew, x_c, y_c, Addpointsx, Addpointsy, Eids, npedge);
          	       
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
/*                	
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
*/                    	    
                        for(i = 0 ; i < Exp.num_elements(); i++)
                        {
                            Exp[i] = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(session,mesh,session->GetVariable(i));
                        }
                    //}

                    break;
                }
            case 2:
                {
/*                	
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
*/                    
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
 //                   }

                    break;
                }
                case 3:     	
                    {
/*                    	    
                        if(HomogeneousType == eHomogeneous3D)
                        {
                            ASSERTL0(false,"3D fully periodic problems not implemented yet");
                        }
                        else
                        {
*/                            
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
                        //}
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
cout<<" edge="<<edge<<endl;   	   
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

	void Computestreakpositions(int nvertl, MultiRegions::ExpListSharedPtr &streak,
    	        Array<OneD, NekDouble> &x,  Array<OneD, NekDouble> &y,
    	        Array<OneD, NekDouble> &xold_up, Array<OneD, NekDouble> &yold_up,
    	        Array<OneD, NekDouble> &xold_low, Array<OneD, NekDouble> &yold_low,
    	        Array<OneD, NekDouble> &xold_c, Array<OneD, NekDouble> &yold_c,    	        
    	        Array<OneD, NekDouble> &xc,  Array<OneD, NekDouble> &yc)   	        
	{
    
	     int nq = streak->GetTotPoints();	
             Array<OneD, NekDouble> coord(2);
             //Array<OneD, NekDouble> stvalues(nvertl,-10);
             Array<OneD, NekDouble> derstreak(nq);
             streak->PhysDeriv(MultiRegions::eY, streak->GetPhys(), derstreak);
             int elmtid, offset;    
             NekDouble U,dU;
             NekDouble F=1000;
             for(int q=0; q<nvertl; q++)
             {

                  NekDouble streaktmp =200;
                  NekDouble streaktmppos=200;
                  NekDouble streaktmpneg=-200;
                  int ipos,ineg, it =0;                  
                  NekDouble weightpos, weightneg;
/*         
           //algorithm to determine the closest streak positions using Newton iteration
           for(int j=0; j<nq; j++)
           {
                if(x[j]==xold_up[q]   && y[j]<= yold_up[q] && y[j]>= yold_low[q])
                {
              	   if(abs(streak->GetPhys()[j])< abs(streaktmp))
              	   {
              	      streaktmp = streak->GetPhys()[j];
              	      it =j;
              	   }
//cout<<" x="<<x[j]<<"  y="<<y[j]<<"   streak="<<streak->GetPhys()[j]<<endl;
                   if(streak->GetPhys()[j]< streaktmppos && streak->GetPhys()[j]>0)
                   {	                 	      
                      streaktmppos = streak->GetPhys()[j];
                      ipos =j;
                   }
                   if(streak->GetPhys()[j]>streaktmpneg && streak->GetPhys()[j]<0)
                   {
              	      streaktmpneg = streak->GetPhys()[j];
              	      ineg =j;
                   }
                }
            }
            x_c[q] = x[it]; 
           //the sign
            
           y_c[q]  = y[it] - (streaktmp)/dU[it];
           ASSERTL0( y_c[q]> y[ineg] && y_c[q]<y[ipos], " wrong position");
cout<<" streak x="<<x_c[q]<<"   y="<<y_c[q]<<" y_pos="<<y[ipos]<<"   y_neg="<<y[ineg]<<endl;           
*/         
cout<<q<<"   xup="<<xold_up[q]<<"   yup="<<yold_up[q]<<"    ydown="<<yold_low[q]<<endl;                 
                    //algorithm to determine the closest points to the streak positions
                    //using the weighted mean         
                    for(int j=0; j<nq; j++)
                    {
                    	 if(x[j]==xold_up[q] && y[j]<= yold_up[q] && y[j]>= yold_low[q])
                    	 {
                             if(streak->GetPhys()[j]< streaktmppos && streak->GetPhys()[j]>0)
                             {	                 	      
                                 streaktmppos = streak->GetPhys()[j];
                                 ipos =j;
                             }
                             if(streak->GetPhys()[j]>streaktmpneg && streak->GetPhys()[j]<0)
                             {
              	                 streaktmpneg = streak->GetPhys()[j];
              	                 ineg =j;
              	             }
              
//cout<<" x="<<x[j]<<"  y="<<y[j]<<"   streak="<<streak->GetPhys()[j]<<endl;
			 }
		    }
cout<<"ipos="<<ipos<<"  ineg="<<ineg<<endl;
//cout<<"closer streak points ypos="<<y[ipos]<<"   yneg="<<y[ineg]<<endl;
	          //determine the streak y position as the result of the weighted mean
	            xc[q]= x[ipos];
	            weightpos = 1/(streaktmppos*streaktmppos);
	            weightneg = 1/(streaktmpneg*streaktmpneg);	            
	            yc[q]= ( (y[ipos]*weightpos) + (y[ineg]*weightneg) )/(weightpos+weightneg);
cout<<" streak x="<<xc[q]<<"   y="<<yc[q]<<endl;
     

         //algorithm to determine the points with u of order 0.08
/*         
                for(int j=0; j<nq; j++)
                {
                       if(x[j]==xold_up[q] && y[j]<= yold_up[q] && y[j]>= yold_low[q])
                       {
                            if(  streak->GetPhys()[j] >0.08  
                  	     && abs(streak->GetPhys()[j] - 0.08) < abs(streaktmppos -0.08)   
                            )
                           {
//cout<<"q="<<q<<" streakpos ="<<streak->GetPhys()[j]<<endl;                  	  
                                 streaktmppos = streak->GetPhys()[j];
                                 ipos =j;                   
                           }
                  
                           if(  streak->GetPhys()[j] < 0.08  && streak->GetPhys()[j]>0 
                  	          && abs(streak->GetPhys()[j] - 0.08) < abs(streaktmpneg -0.08)   
                            )
                           {
//cout<<"q="<<q<<" streakneg ="<<streak->GetPhys()[j]<<endl;                  	  
                                streaktmpneg = streak->GetPhys()[j];
                                ineg =j;                   
                           }               
                 
                        }            
                  }
       
                  //determine streak position as the weghted mean         

                  if(streaktmppos == 200)
                  {         	 
                      y_c[q] = y[ineg];
                      x_c[q] = x[ineg];             
                  }
                  else if(streaktmpneg == -200)
                  {         	 
                       y_c[q] = y[ipos];
                       x_c[q] = x[ipos];             
                  }
                  else
                  {
                       x_c[q] = x[ipos];
	               weightpos = 1/abs( (streaktmppos-0.08)*(streaktmppos-0.08)  );
	               weightneg = 1/abs( (streaktmpneg-0.08)*(streaktmpneg-0.08)  );
	               y_c[q]= ( (y[ipos]*weightpos) + (y[ineg]*weightneg) )/(weightpos+weightneg);
	          }
	 
cout<<" streak x="<<x_c[q]<<"   y="<<y_c[q]<<" streak_p="<<streaktmppos<<"   streak_n="<<streaktmpneg<<endl;          
                  if(q>0 && q< nvertl-1)
                  {       	 
                      ASSERTL0(y_c[q+1] < y_c[q], " the critical layer is oscillating");         
                  }
*/         
              }
              for(int e=0; e<nvertl; e++)
              {
              	   coord[0] =xc[e];
              	   coord[1] =yc[e];
                	   
              	   elmtid = streak->GetExpIndex(coord,0.00001);
           	   offset = streak->GetPhys_Offset(elmtid);
           	   F = 1000;
		   while( abs(F)> 0.00000001)
		   {
		   	U = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);
		   	dU  = streak->GetExp(elmtid)->PhysEvaluate(coord, derstreak + offset);
		   	coord[1] = coord[1] - U/dU;   
		   	F = U;   
		   	ASSERTL0( coord[0]==xc[e], " x coordinate must remain the same");
              	        //stvalues[e] = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() +offset );
//cout<<"elmtid="<<elmtid<<"  x="<<coord[0]<<"   y="<<coord[1]<<"    stvalue="<<U<<endl;
	           }
                   yc[e] = coord[1];
	           //Utilities::Zerofunction(coord[0], coord[1], xtest, ytest, streak, derstreak);
//cout<<"result x="<<xc[e]<<"  y="<<yc[e]<<"   streak="<<U<<endl;	           
              }
           
		
	}
	void GenerateAddPointsNewtonIt( NekDouble &xi, NekDouble &yi,NekDouble &x0, NekDouble &y0,
    	        MultiRegions::ExpListSharedPtr &function, Array<OneD, NekDouble> &derfunction)
	{
		int elmtid,offset;
		NekDouble F,U,dU;
		Array<OneD, NekDouble> coords(2);	
		coords[0] = xi;
		coords[1] = yi;	
//cout<<"generate newton it xi="<<xi<<"  yi="<<yi<<endl;			
		elmtid = function->GetExpIndex(coords, 0.00001);
                //@to do if GetType(elmtid)==triangular WRONG!!!
//cout<<"gen newton xi="<<xi<<"  yi="<<yi<<"  elmtid="<<elmtid<<endl;			
		offset = function->GetPhys_Offset(elmtid);
		F =1000;
		while( abs(F)> 0.00000001)
	        {
		     U = function->GetExp(elmtid)->PhysEvaluate(coords, function->GetPhys() + offset);
		     dU  = function->GetExp(elmtid)->PhysEvaluate(coords, derfunction + offset);
		     coords[1] = coords[1] - U/dU;   
		     F = U;   
		     ASSERTL0( coords[0]==xi, " x coordinate must remain the same");	             	
	        }
	        x0 = xi;
	        y0 = coords[1];
//cout<<"NewtonIt result  x="<<x0<<"  y="<<coords[1]<<"   U="<<U<<endl;	        
	}

        void GenerateCurve(int npoints, int npused, Array<OneD, NekDouble> &x_c, 
    	    Array<OneD, NekDouble> &y_c, Array<OneD, NekDouble> &curve)
    	{
    	    int N= npoints;
    	    int totpoints = npoints + npused;
//cout<<"totppoints="<<totpoints<<endl;    	    
    	    Array<OneD, NekDouble> A (N*N,1.0);
    	    Array<OneD, NekDouble> b (N);
            int row=0;
    	    //fill column by column
    	    for(int e=0; e<N; e++)
    	    {
    	    	row=0;
    	    	for(int w=npused; w < totpoints; w++)
    	    	{         
                     A[N*e+row] = std::pow( x_c[w], N-1-e);
                     row++;
                }	    
            }
            row=0;
            for(int r= npused; r< totpoints; r++)
            {
    	        b[row] =   y_c[r];
//cout<<"b="<<b[row]<<"  y_c="<<y_c[r]<<endl;     	        
    	        row++;       
	    }
	    
//cout<<"A elements="<<A.num_elements()<<endl;
            Array<OneD, int> ipivot (N);  
            int info =0;
            //Lapack::Dgesv( N, 1, A.get(), N, ipivot.get(),  b.get(), N, info);     
            Lapack::Dgetrf( N, N, A.get(), N, ipivot.get(), info);
    
            if( info < 0 )
            {
               std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
                    "th parameter had an illegal parameter for dgetrf";
               ASSERTL0(false, message.c_str());
            }
            else if( info > 0 )
            {
               std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
                    boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
               ASSERTL0(false, message.c_str());
            }
   
            // N means no transponse (direct matrix)
            int ncolumns_b =1;
            Lapack::Dgetrs( 'N', N, ncolumns_b , A.get() , N, ipivot.get(), b.get(), N, info);
            if( info < 0 )
            {
                std::string message = "ERROR: The " + boost::lexical_cast<std::string>(-info) + 
                    "th parameter had an illegal parameter for dgetrf";
                ASSERTL0(false, message.c_str());
            }
            else if( info > 0 )
            {
                std::string message = "ERROR: Element u_" + boost::lexical_cast<std::string>(info) +
                    boost::lexical_cast<std::string>(info) + " is 0 from dgetrf";
                ASSERTL0(false, message.c_str());
            }
//cout<<"coeffs a,b,c,d..."<<endl;
	    //fill the coeffs array
	    if(npused==0)
	    {
	        Vmath::Vcopy(totpoints, &(b[0]), 1, &(curve[npused]), 1);	
//cout<<b[1]<<"   ipiv="<<ipivot[1]<<endl;    	        
	    }
	    else
	    {
	        row=0;
    	        for(int a= npused+1; a< totpoints+1; a++)
    	        {
    	    	   curve[a] = b[row];
    	    	   row++;
//cout<<curve[a]<<"   ipiv="<<ipivot[row-1]<<endl;    	    
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
	
	
        void GenerateAddPoints(int region, SpatialDomains::MeshGraphSharedPtr &mesh,int np, int npused,
    	         Array<OneD, NekDouble> &curve, MultiRegions::ExpListSharedPtr & bndfield, 
    	         Array<OneD, NekDouble>& outx, Array<OneD, NekDouble>& outy, Array<OneD, int>&Eids)
    	{
            LocalRegions::SegExpSharedPtr bndSegExp;
            
            int Eid,id1,id2;
            NekDouble x1,y1,z1;
            NekDouble x2,y2,z2;
            SpatialDomains::VertexComponentSharedPtr vertex1;
            SpatialDomains::VertexComponentSharedPtr vertex2;

            int firstedge=0;
	    int firstcoeff=0;      	    
            if(npused!=0)
            {
              firstedge  = npused;
              firstcoeff = firstedge+1;

            }
            int lastedge  = firstedge + (np-1);
            int lastcoeff  = firstcoeff + np;
//cout<<"firstedge="<<firstedge<<"  lastedge="<<lastedge<<"  firstcoeff="<<firstcoeff
//                              <<"  lastcoeff="<<lastcoeff<<endl;
            int polorder;	    
            for(int s= firstedge; s< lastedge; s++)
            {
            	 bndSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndfield->GetExp(s));   
            	 Eid = (bndSegExp->GetGeom1D())->GetEid();
            	 id1 = (bndSegExp->GetGeom1D())->GetVid(0);
            	 id2 = (bndSegExp->GetGeom1D())->GetVid(1);  
                 vertex1 = mesh->GetVertex(id1);      
		 vertex2 = mesh->GetVertex(id2);     
		 vertex1->GetCoords(x1,y1,z1);
		 vertex2->GetCoords(x2,y2,z2);	
//cout<<"x1="<<x1<<"  x2="<<x2<<endl;		 
		 if(x2>x1)
		 {
		     outx[s] = x1 +(x2-x1)/2;
		     if( outx[s]>x2 || outx[s]< x1)
		     {
		     	  outx[s] = -outx[s];
		     }
		 }
		 else if(x1>x2)
	         {	         	 
	             outx[s] = x2+ (x1-x2)/2;
	             if( outx[s] > x1 || outx[s] <x2)
	             {
	             	  outx[s] = -outx[s];
	             }
	         }
	         else
	         {
	             ASSERTL0(false, "point not generated"); 	 
	         }
                 polorder = np-1;
	         for(int g= firstcoeff; g< lastcoeff; g++)
	         {	         	 
	             outy[s] += curve[g]*(std::pow(outx[s], polorder));
//cout<<"coeff*x^i="<<outy[s]<<"  coeff="<<curve[g]<<"  exp="<<polorder<<endl;
                     ASSERTL0(polorder >=0, " polynomial with one negative exponent");
                     polorder--;
	         }
	         Eids[s] = Eid;
	         
//cout<<"eid="<<Eids[s]<<"  xadd="<<outx[s]<<"   yadd="<<outy[s]<<endl;            	 
            }
            
        }
        
        void Replacevertices(string filename, Array<OneD, NekDouble> newx, 
    	               Array<OneD,  NekDouble> newy,
    	               Array<OneD, NekDouble> x_lay, Array<OneD, NekDouble> y_lay,
    	               Array<OneD, NekDouble> & Pcurvx, 
    	               Array<OneD, NekDouble> & Pcurvy,
    	               Array<OneD, int>&Eids, int Npoints)	
	{     
	    //load existing file
            string newfile;
	    TiXmlDocument doc(filename); 
	    bool loadOkay = doc.LoadFile();
            //load xscale parameter (if exists)
	    TiXmlElement* master = doc.FirstChildElement("NEKTAR");
	    TiXmlElement* mesh = master->FirstChildElement("GEOMETRY");
	    TiXmlElement* element = mesh->FirstChildElement("VERTEX");
            NekDouble xscale = 1.0;
            LibUtilities::AnalyticExpressionEvaluator expEvaluator;
            const char *xscal = element->Attribute("XSCALE");
            if(xscal)
            {
                 std::string xscalstr = xscal;
                 int expr_id  = expEvaluator.DefineFunction("",xscalstr);
                 xscale = expEvaluator.Evaluate0(expr_id);
            }
            

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
            //set xscale 1!!
            if(xscale!=1.0)
            {
                 elementnew->SetAttribute("XSCALE",1.0);
            }
	    TiXmlElement *vertexnew = elementnew->FirstChildElement("V");

      	    
	    int indx;
	    int err, numPts;
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
       	           //remove the old coordinates
       	           vertexnew->RemoveChild(vertexBody);
                   //write the new one   
//cout<<"writing.. v:"<<nextVertexNumber<<endl;	    	    	
	    	   stringstream s;
		   //we need at least 5 digits (setprecision 5) to get the streak position with
		   // precision 10^-10
	    	   s << std::scientific << std::setprecision(8) <<  newx[nextVertexNumber] << "   "
	    	   << newy[nextVertexNumber] << "   " << 0.0;
	    	   vertexnew->LinkEndChild(new TiXmlText(s.str()));      
	    	   //TiXmlNode *newvertexBody = vertexnew->FirstChild();
	    	   //string newvertexbodystr= newvertexBody->SetValue(s.str());                     
	    	   //vertexnew->ReplaceChild(vertexBody,new TiXmlText(newvertexbodystr));
	    	 
	    	   vertexnew = vertexnew->NextSiblingElement("V");  
       	   }	
       	    //meshnew->LinkEndChild(elementnew);
       	    
       	    //read the curved tag
            TiXmlElement* curvednew = meshnew->FirstChildElement("CURVED");
	    ASSERTL0(curvednew, "Unable to find mesh CURVED tag in file.");            
   	    TiXmlElement *edgenew = curvednew->FirstChildElement("E");
   	    int cnt =-1;
            //ID is different from index...
   	    std::string charindex;
   	    int eid;
            int index;
   	    int indexeid, v1,v2;
   	    NekDouble x1,x2,y1,y2;
   	    while (edgenew)
   	    {
   	    	   indexeid =-1;
   	    	   cnt++;
                   //get the index...
   	    	   TiXmlAttribute *edgeAttr = edgenew->FirstAttribute();                   
   	    	   std::string attrName(edgeAttr->Name());
   	    	   charindex = edgeAttr->Value();
   	    	   std::istringstream iss(charindex);
   	    	   iss >> std::dec >> index; 
		   //get the eid
                   edgenew->QueryIntAttribute("EDGEID", &eid);
//cout<<"eid="<<eid<<" neid="<<Eids.num_elements()<<endl;
	           //find the corresponding index curve point
	           for(int u=0; u<Eids.num_elements(); u++)
	           {
//cout<<"Eids="<<Eids[u]<<endl;	           	   
	               if(Eids[u]==eid)
	               {
	                  indexeid = u;
	               }	               
	           }
	           if(indexeid==-1)
   	    	   {
   	    	      ASSERTL0(false, "edge to update not found");	   
   	    	   }


                   std::string edgeBodyStr;
		   //read the body of the edge
		   TiXmlNode *edgeBody = edgenew->FirstChild();
		   if(edgeBody->Type() == TiXmlNode::TEXT)
		   {
		      edgeBodyStr += edgeBody->ToText()->Value();
		      edgeBodyStr += " ";
		   }
   	    	   ASSERTL0(!edgeBodyStr.empty(), "Edge definitions must contain edge data");
   	    	   //remove the old coordinates
   	    	   edgenew->RemoveChild(edgeBody);
   	    	   //write the new points coordinates
   	    	   if(x_lay[cnt]< x_lay[cnt+1])
   	    	   {
   	    	      //x1 = x_lay[cnt];
   	    	      //x2 = x_lay[cnt+1];
   	    	      //y1 = y_lay[cnt];
   	    	      //y2 = y_lay[cnt+1];
   	    	      v1 = cnt;
   	    	      v2 = cnt+1;
   	    	   }
   	    	   else
   	    	   {
   	    	      //x1 = x_lay[cnt+1];
   	    	      //x2 = x_lay[cnt];
   	    	      //y1 = y_lay[cnt+1];
   	    	      //y2 = y_lay[cnt];
   	    	      v1 = cnt+1;
   	    	      v2 = cnt;
                      cout<<"Warning: the edges can have the vertices in the reverse order";   	    	      
   	    	   }
		   //we need at least 5 digits (setprecision 5) to get the streak position with
		   // precision 10^-10
		   
                   //Determine the number of points
                   err = edgenew->QueryIntAttribute("NUMPOINTS", &numPts);
                   ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute NUMPOINTS.");
                   edgenew->SetAttribute("NUMPOINTS", Npoints);
		   stringstream st;
		   st << std::scientific << std::setprecision(8) << x_lay[v1] << "   "
	    	   << y_lay[v1] << "   " << 0.000<<"   ";

		   for(int a=0; a< Npoints-2; a++)
		   {		       	   
		      st << std::scientific << std::setprecision(8) <<
		      "    "<<Pcurvx[indexeid*(Npoints-2) +a]<<"    "<<Pcurvy[indexeid*(Npoints-2) +a]
		      <<"    "<<0.000<<"   "; 
	      
		   }
		   st << std::scientific << std::setprecision(8) <<
		   "    "<<x_lay[v2]<<"   "<< y_lay[v2] <<"   "<< 0.000;
	    	   edgenew->LinkEndChild(new TiXmlText(st.str()));    		   

//cout<<st.str()<<endl;		
/*
	    	   stringstream s;
	    	   s << std::scientific << std::setprecision(5) <<  x_lay[v1] << "   "
	    	   << y_lay[v1] << "   " << 0.000<<"        "<< Pcurvx[indexeid]<< "   "<<
	    	   Pcurvy[indexeid]<<"   "<<0.000 <<"       "<<x_lay[v2]<<"   "<< y_lay[v2] 
	    	   <<"   "<< 0.000;
	    	   edgenew->LinkEndChild(new TiXmlText(s.str()));    
*/  
/*	    	   
                   if(x_lay[cnt] < x_lay[cnt+1])
                   {
                      cout<<"Warning: the edges can have the vertices in the reverse order";
                   }
*/                   
   	    	   edgenew = edgenew->NextSiblingElement("E");
   	    	    
   	    }
   	    
   	    
   	    
   	    
   	    
       	    docnew.SaveFile( newfile ); 
       	    
       	    cout<<"new file:  "<<newfile<<endl;
	   
	}    	
	
