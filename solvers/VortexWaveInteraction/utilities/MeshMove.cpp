///////////////////////////////////////////////////////////////////////////////
//
// File MoveMeshToCriticalLayer.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: MoveMesh to critical layer
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
//#include <sstream>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
//#include <MultiRegions/ExpList0D.h>
//#include <MultiRegions/ExpList2D.h>
//#include <MultiRegions/ExpList3D.h>
//#include <MultiRegions/ExpList2DHomogeneous1D.h>
//#include <MultiRegions/ExpList3DHomogeneous1D.h>
//#include <MultiRegions/ExpList1DHomogeneous2D.h>
//#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
//#include <MultiRegions/ContField3D.h>
//#include <MultiRegions/ContField3DHomogeneous1D.h>
//#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
//#include </PreProcessing/MeshConvert/Convert.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
//#include <LibUtilities/Foundations/GaussPoints.h>
#include <boost/lexical_cast.hpp>
#include <tinyxml/tinyxml.h>


void OrderVertices(int nedges,SpatialDomains::MeshGraphSharedPtr graphShPt,
    	    MultiRegions::ExpListSharedPtr & bndfield, 
        	Array<OneD, int>& Vids, int v1,int v2 , NekDouble x_connect,
        	int & lastedge, 
        	Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y);
void Computestreakpositions(int nvertl, MultiRegions::ExpListSharedPtr streak,
    	        Array<OneD, NekDouble> xold_up, Array<OneD, NekDouble> yold_up,
    	        Array<OneD, NekDouble> xold_low, Array<OneD, NekDouble> yold_low,
    	        Array<OneD, NekDouble> xold_c, Array<OneD, NekDouble> yold_c,      	        
    	        Array<OneD, NekDouble> &xc,  Array<OneD, NekDouble> &yc, NekDouble cr,
                bool verts); 
void GenerateAddPointsNewtonIt( NekDouble xi, NekDouble yi,NekDouble &xout, NekDouble &yout,
    	        MultiRegions::ExpListSharedPtr function, Array<OneD, NekDouble> derfunction,
                NekDouble cr);

void GenerateMapEidsv1v2(MultiRegions::ExpListSharedPtr field,
                 Array<OneD, int> &V1, Array<OneD, int> &V2);

void MappingEVids(Array<OneD, NekDouble> xoldup, Array<OneD, NekDouble> yoldup,
                 Array<OneD, NekDouble> xolddown, Array<OneD, NekDouble> yolddown,
                 Array<OneD, NekDouble> xcold, Array<OneD, NekDouble> ycold,
                 Array<OneD, int> Vids_c,
                 SpatialDomains::MeshGraphSharedPtr mesh,
                 MultiRegions::ExpListSharedPtr streak,
                 Array<OneD, int> V1, Array<OneD, int> V2,
                 int & nlays,  Array<OneD, Array<OneD, int> >& Eids_lay, 
                 Array<OneD, Array<OneD, int> >& Vids_lay);
bool checkcommonvert(Array<OneD, int> Vids_laybefore, Array<OneD, int> Vids_c, int Vid);

void  Cutrepetitions(int nedges,Array<OneD, NekDouble> inarray,
                 Array<OneD, NekDouble>& outarray);

int DetermineclosePointxindex(NekDouble x,Array<OneD, NekDouble> xArray);

void GenerateNeighbourArrays(int index, int neighpoints,Array<OneD, NekDouble> xArray,
                 Array<OneD, NekDouble> yArray,Array<OneD, NekDouble>& Neighbour_x,
                 Array<OneD, NekDouble>& Neighbour_y);

NekDouble LagrangeInterpolant(NekDouble x, int npts, 
                 Array<OneD,NekDouble>  xpts, Array<OneD, NekDouble> funcvals);

void EvaluateTangent(int npoints, Array<OneD, NekDouble> xcQedge, 
                 Array<OneD, NekDouble> coeffsinterp,
                 Array<OneD, NekDouble> & txQedge, Array<OneD, NekDouble> & tyQedge);
void PolyInterp( Array<OneD, NekDouble> xpol, Array<OneD, NekDouble> ypol,
                 Array<OneD, NekDouble> & coeffsinterp,
                 Array<OneD, NekDouble> & xcout, Array<OneD, NekDouble> & ycout,  
                 int edge, int npedge);                
void PolyFit(int polyorder,int npoints,
                 Array<OneD, NekDouble>  xin, Array<OneD, NekDouble>  fin,
                 Array<OneD, NekDouble> & coeffsinterp,
                 Array<OneD, NekDouble> & xout, Array<OneD, NekDouble> & fout,  
                 int npout);
void  Orderfunctionx(Array<OneD, NekDouble> inarray_x,
                 Array<OneD, NekDouble> inarray_y, Array<OneD, NekDouble>& outarray_x,
                 Array<OneD, NekDouble>& outarray_y);

void MoveLayersvertically(int nlays, int nvertl, int cntlow, int cntup,
    	         Array<OneD, Array<OneD, int > > lay_Vids,  Array<OneD, NekDouble> xc,
    	         Array<OneD, NekDouble> yc, Array<OneD, int> Down, Array<OneD, int> Up,
    	         Array<OneD, NekDouble >& xnew, Array<OneD, NekDouble>& ynew,
    	         Array<OneD, Array<OneD, NekDouble > >& layers_x,
    	         Array<OneD, Array<OneD, NekDouble > >& layers_y);

void MoveLayerNfixedxpos(int nvertl, int npedge, Array<OneD, NekDouble> xcPhys,
	         Array<OneD, NekDouble> tmpx_lay, Array<OneD, NekDouble> tmpy_lay,
	         Array<OneD, int> Vids,
	         Array<OneD, NekDouble> &xlay, Array<OneD, NekDouble> &ylay,
	         Array<OneD, NekDouble> &xnew, Array<OneD, NekDouble> &ynew);

void MoveLayerNnormpos(int nvertl, int npedge, Array<OneD, NekDouble> xcPhys,
	         Array<OneD, NekDouble> tmpx_lay, Array<OneD, NekDouble> tmpy_lay,
	         Array<OneD, int> Vids,
	         Array<OneD, NekDouble> &xlay, Array<OneD, NekDouble> &ylay,
	         Array<OneD, NekDouble> &xnew, Array<OneD, NekDouble> &ynew);

void MoveOutsidePointsfixedxpos(int npedge, SpatialDomains::MeshGraphSharedPtr mesh,
	         Array<OneD, NekDouble> xcold,Array<OneD, NekDouble> ycold,
      	         Array<OneD, NekDouble> xolddown,Array<OneD, NekDouble> yolddown,
	         Array<OneD, NekDouble> xoldup,Array<OneD, NekDouble> yoldup,     
	         Array<OneD, NekDouble> ylaydown,Array<OneD, NekDouble> ylayup, 	         
                 Array<OneD, NekDouble>& xnew,Array<OneD, NekDouble>& ynew);

void MoveOutsidePointsNnormpos(int npedge, SpatialDomains::MeshGraphSharedPtr mesh,
	         Array<OneD, NekDouble> xcold,Array<OneD, NekDouble> ycold,
      	         Array<OneD, NekDouble> xolddown,Array<OneD, NekDouble> yolddown,
	         Array<OneD, NekDouble> xoldup,Array<OneD, NekDouble> yoldup,     
	         Array<OneD, NekDouble> xlaydown,Array<OneD, NekDouble> ylaydown,
	         Array<OneD, NekDouble> xlayup,Array<OneD, NekDouble> ylayup, 
	         Array<OneD, NekDouble> nxPhys,Array<OneD, NekDouble> nyPhys,	         
                 Array<OneD, NekDouble>& xnew,Array<OneD, NekDouble>& ynew)  ;

void CheckSingularQuads( MultiRegions::ExpListSharedPtr Exp, 
                 Array<OneD, int> V1, Array<OneD, int> V2,	         
	         Array<OneD, NekDouble>& xnew,Array<OneD, NekDouble>& ynew);

void Replacevertices(string filename, Array<OneD, NekDouble> newx, 
    	               Array<OneD,  NekDouble> newy,
    	               Array<OneD, NekDouble> xcPhys, Array<OneD, NekDouble> ycPhys,
    	               Array<OneD, int>Eids, int Npoints, string s_alp,
                       Array<OneD, Array<OneD, NekDouble> > x_lay,
                       Array<OneD, Array<OneD, NekDouble> > y_lay,
                       Array<OneD, Array<OneD, int > >lay_eids, bool curv_lay)	;

int main(int argc, char *argv[])
{
    NekDouble cr;
    //set cr =0
    cr=0;
    //change argc from 6 to 5 allow the loading of cr to be optional
    if(argc > 6 || argc < 5)
    {
        fprintf(stderr,
            "Usage: ./MoveMesh  meshfile fieldfile  changefile   alpha  cr(optional)\n");
        exit(1);
    }
    
    //ATTEnTION !!! with argc=2 you impose that vSession refers to is argv[1]=meshfile!!!!! 
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(2, argv);
    //----------------------------------------------
        
    if( argc == 6 &&
        vSession->DefinesSolverInfo("INTERFACE")
        && vSession->GetSolverInfo("INTERFACE")=="phase" )
    {
        cr = boost::lexical_cast<NekDouble>(argv[argc-1]);
        argc=5;
    }

    // Read in mesh from input file
    string meshfile(argv[argc-4]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile);
    //---------------------------------------------- 

    // Also read and store the boundary conditions
    SpatialDomains::BoundaryConditionsSharedPtr boundaryConditions;        
    boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                                        ::AllocateSharedPtr(vSession,graphShPt);
    SpatialDomains::BoundaryConditions bcs(vSession, graphShPt); 
    //----------------------------------------------

  
    //the mesh file should have 2 component: set output fields
    //fields has to be of the SAME dimension of the mesh (that's why there is
    //the changefile as an input) 
    //a contfield2D is needed to extract boundary conditions!!!

    // store name of the file to change
    string changefile(argv[argc-2]);
    //----------------------------------------------
      
    //store the value of alpha
    string charalp (argv[argc-1]);
    //NekDouble alpha = boost::lexical_cast<NekDouble>(charalp);
    cout<<"read alpha="<<charalp<<endl;

    //---------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-3]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    //---------------------------------------------- 

    MultiRegions::ExpListSharedPtr streak; 
    streak = MemoryManager<MultiRegions::ContField2D>
          ::AllocateSharedPtr(vSession, graphShPt, "w",true);    

    for(int i=0; i<fielddata.size(); i++)
    {
        streak->ExtractDataToCoeffs(fielddef[i], fielddata[i], fielddef[i]->m_fields[0], streak->UpdateCoeffs());
    }
    streak->BwdTrans_IterPerExp(streak->GetCoeffs(), streak->UpdatePhys());
 
    //------------------------------------------------   
    // determine the I regions (3 region expected)
    // hypothesys: the layes have the same number of points

    int nIregions, lastIregion; 
    const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConditions  = streak->GetBndConditions();    
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
    Array<OneD, MultiRegions::ExpListSharedPtr> bndfieldx= streak->GetBndCondExpansions();        
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
    NekDouble x0,y0,z0,xt,yt=0,zt=0;
    int lastedge=-1;
    int v1,v2;
    //first point for x_connect=0(or-1.6 for the full mesh (-pi,pi)  )
    x_connect=0;
    SpatialDomains::VertexComponentSharedPtr vertex0 =
       graphShPt->GetVertex
       (
       ( (boost::dynamic_pointer_cast<LocalRegions
           ::SegExp>(bndfieldx[lastIregion]->GetExp(0))
         )->GetGeom1D()
       )
       ->GetVid(0)
       );
    vertex0->GetCoords(x0,y0,z0);
    if( x0 != 0.0)
    {
cout<<"WARNING x0="<<x0<<endl;
       x_connect=x0;       
    }
    v1=0;
    v2=1;
    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-1], 
       	Vids_low, v1, v2 , x_connect ,lastedge, xold_low,yold_low);    
    ASSERTL0(Vids_low[v2]!=-10, "Vids_low[v2] is wrong");
    SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_low[v2]);    

    //update x_connect    
cout<<"x_conn="<<x_connect<<"   yt="<<yt<<"  zt="<<zt<<" vid="<<Vids_low[v2]<<endl;
    vertex->GetCoords(x_connect,yt,zt);
      
    int i=2;
    while(i<nvertl)
    {     	    
         v1=i;
         OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-1], 
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
    //first point for x_connect=0 (or-1.6)
    x_connect=0;
    vertex0 =
       graphShPt->GetVertex
       (
       ( (boost::dynamic_pointer_cast<LocalRegions
           ::SegExp>(bndfieldx[lastIregion]->GetExp(0))
         )->GetGeom1D()
       )
       ->GetVid(0)
       );
    vertex0->GetCoords(x0,y0,z0);
    if( x0 != 0.0)
    {
cout<<"WARNING x0="<<x0<<endl;
       x_connect=x0;       
    }
    lastedge=-1;

    v1=0;
    v2=1;
    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion-2 ], 
        	Vids_up, v1, v2 , x_connect ,lastedge, xold_up, yold_up);    
    SpatialDomains::VertexComponentSharedPtr vertexU = graphShPt->GetVertex(Vids_up[v2]);    
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
    //first point for x_connect=0(or-1.6)
    x_connect=0;
    vertex0 =
        graphShPt->GetVertex(((boost::dynamic_pointer_cast<LocalRegions
                               ::SegExp>(bndfieldx[lastIregion]->GetExp(0)))->GetGeom1D()
                              )->GetVid(0));
    vertex0->GetCoords(x0,y0,z0);
    if( x0 != 0.0)
    {
        cout<<"WARNING x0="<<x0<<endl;
       x_connect=x0;       
    }
    lastedge=-1;

    v1=0;
    v2=1;

    OrderVertices(nedges, graphShPt, bndfieldx[lastIregion], 
        	Vids_c, v1, v2 , x_connect ,lastedge, xold_c, yold_c);    
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
//cout<<"Vids cl="<<Vids_low[v1]<<endl;  
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


    //fieds to force continuity:
    const SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
    MultiRegions::ContField1DSharedPtr  Cont_y;
    MultiRegions::ExpList1DSharedPtr yexp;

    yexp = MemoryManager<MultiRegions::ExpList1D>
    		::AllocateSharedPtr(*(bregions[lastIregion]), graphShPt, true);
    Cont_y = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(vSession, *yexp);  
    //--------------------------------------
    /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/ 
    //generate additional points using the Newton iteration                
    //determine the xposition for every edge (in the middle even if it 
    // is not necessary
    //PARAMETER which determines the number of points per edge @todo put as an input      
    int npedge;
    if(vSession->DefinesParameter("npedge"))
    {
          npedge = (int)vSession->GetParameter("npedge"); 
    }
    else
    {
          npedge = 5;//default value
    }
    /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/ 
    //find the points where u=0 and determine the sign of the shift and the delta
    int nq= streak->GetTotPoints(); 
    Array<OneD, NekDouble> x(nq);
    Array<OneD,NekDouble> y(nq);    
    streak->GetCoords(x,y);         
    Array<OneD, NekDouble> x_c(nvertl);
    Array<OneD,NekDouble> y_c(nvertl,-200);       
    Array<OneD, NekDouble> tmp_w (nvertl, 200);
    Array<OneD, int> Sign (nvertl,1);   
    Array<OneD, NekDouble> Delta_c(nvertl,-200);
   

    Computestreakpositions(nvertl, streak, xold_up, yold_up,
    	                   xold_low, yold_low, xold_c, yold_c, x_c, y_c,cr,true);    
    // if the curve is low the old layer point, it has to shift down  
    NekDouble shift;  
    for(int q=0; q<nvertl; q++)
    {
         if(y_c[q] < yold_c[q])
         {
             Sign[q] = -1;
         }
         //calculate delta
         Delta_c[q] = abs(yold_c[q]-y_c[q]);
         //check the shifting of the layer:
         shift+= Delta_c[q];         
         cout<<x_c[q]<<"    "<<y_c[q]<<endl;
    }
    //cout<<"shift="<<shift<<endl;
    if(shift<0.001)
    {
         cout<<"Warning: the critical layer is stationary"<<endl;
    }
    //------------------------------------------------------------------


    //additional points arrays
    Array<OneD, NekDouble> Cpointsx (nedges);
    Array<OneD, NekDouble> Cpointsy (nedges, 0.0);
    Array<OneD, int> Eids (nedges);    
    Array<OneD, NekDouble> Addpointsx (nedges*(npedge-2), 0.0);
    Array<OneD, NekDouble> Addpointsy (nedges*(npedge-2), 0.0); 
    //calculate the dU_dy 
    Array<OneD, NekDouble> dU(streak->GetTotPoints());
    streak->PhysDeriv(MultiRegions::eY, streak->GetPhys(), dU);

   
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
         cout<<"edge="<<r<<"  x1="<<x1<<"  y1="<<y1<<"   x2="<<x2<<"  y2="<<y2<<endl;
	 if(x2>x1)
	 {
	     Cpointsx[r] = x1 +(x2-x1)/2;             
             //cout<<"edge="<<r<<"  x1="<<x1<<"  x2="<<x2<<"   Cx="<<Cpointsx[r]<<endl;
             //cout<<"edge="<<r<<"  x1="<<x1<<"  y1="<<y1<<"   x2="<<x2<<"  y2="<<y2<<endl;
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
                 //initial guess along the line defined by the NEW verts y_c                	 
	         Addpointsy[r*(npedge-2) +w] = y_c[r] + ((y_c[r+1]-y_c[r])/(x_c[r+1]-x_c[r]))*(Addpointsx[r*(npedge-2) +w]-x1);
	         //Addpointsy[r*(npedge-2) +w] = y1 + ((y2-y1)/(x2-x1))*(Addpointsx[r*(npedge-2) +w]-x1);	         
	         GenerateAddPointsNewtonIt( Addpointsx[r*(npedge-2) +w], Addpointsy[r*(npedge-2) +w],
	              	      Addpointsx[r*(npedge-2) +w],  Addpointsy[r*(npedge-2) +w], streak, dU,cr); 

               // Lay_x->UpdatePhys()[r*npedge +1 +w]= Addpointsx[r*(npedge-2) +w];
               // Lay_y->UpdatePhys()[r*npedge +1 +w]= Addpointsy[r*(npedge-2) +w];                
             }
             //Lay_x->UpdatePhys()[r*npedge +0] =x1;
             //Lay_y->UpdatePhys()[r*npedge +0] =y1;
             //Lay_x->UpdatePhys()[r*npedge +npedge-1] =x2;
             //Lay_y->UpdatePhys()[r*npedge +npedge-1] =y2;

	 }
	 else if(x1>x2)
	 {	         	 
	     Cpointsx[r] = x2+ (x1-x2)/2;
//cout<<"edge="<<r<<"  y1="<<y1<<"  y2="<<y2<<endl;
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

                 //initial guess along the line defined by the NEW verts y_c                	 
	         Addpointsy[r*(npedge-2) +w] = y_c[r+1] + ((y_c[r]-y_c[r+1])/(x_c[r]-x_c[r+1]))*(Addpointsx[r*(npedge-2) +w]-x2);

	         //Addpointsy[r*(npedge-2) +w] = y2 + ((y1-y2)/(x1-x2))*(Addpointsx[r*(npedge-2) +w]-x2);	         	         
	         GenerateAddPointsNewtonIt( Addpointsx[r*(npedge-2) +w], Addpointsy[r*(npedge-2) +w], 
	               Addpointsx[r*(npedge-2) +w], Addpointsy[r*(npedge-2) +w], streak, dU,cr); 
               // Lay_x->UpdatePhys()[r*npedge +1]= Addpointsx[r*(npedge-2) +w];
               // Lay_y->UpdatePhys()[r*npedge +1]= Addpointsy[r*(npedge-2) +w];   
             }
             //Lay_x->UpdatePhys()[r*npedge +0] =x2;
             //Lay_y->UpdatePhys()[r*npedge +0] =y2;
             //Lay_x->UpdatePhys()[r*npedge +npedge-1] =x1;
             //Lay_y->UpdatePhys()[r*npedge +npedge-1] =y1;	      
	 }
	 else
	 {
	      ASSERTL0(false, "point not generated"); 	 
	 }    	    
//cout<<"calculate cpoints coords"<<endl;	 
         //Cpointsy[r] = y1 + (y2-y1)/2;
//cout<<"central point:"<<endl;
         //GenerateAddPointsNewtonIt( Cpointsx[r], Cpointsy[r],Cpointsx[r], Cpointsy[r],
    	 //      streak, dU,cr); 
         //NekDouble diff = Cpointsy[r]-Addpointsy[r*(npedge-2)];
//cout<<"diff="<<diff<<endl;         
	 Eids[r] = Eid;

    }      
    //-------------------------------------------------------------


    //fill the xPhys,yPhys array( may necessary after)    
    Array<OneD, NekDouble> xcPhys (nedges*npedge, 0.0);
    Array<OneD, NekDouble> ycPhys (nedges*npedge, 0.0);

    for(int a=0; a<nedges; a++)
    {
    	 //v1
         xcPhys[a*npedge+0] = x_c[a];       
         ycPhys[a*npedge+0] = y_c[a];    
         //v2
         xcPhys[a*npedge+npedge-1] = x_c[a+1];       
         ycPhys[a*npedge+npedge-1] = y_c[a+1];           
    	    
         for(int b=0; b<npedge-2; b++)
         {
               xcPhys[a*npedge +b+1] = Addpointsx[a*(npedge-2)+b];
               ycPhys[a*npedge +b+1] = Addpointsy[a*(npedge-2)+b];	
         }
    }

cout<<"xc,yc before tanevaluate"<<endl;
for(int v=0; v< xcPhys.num_elements(); v++)
{
cout<<xcPhys[v]<<"     "<<ycPhys[v]<<endl;
}

    //-------------------------------------------------

    //V1[eid],V2[eid] vertices associate with the edge Id=eid
    Array<OneD, int> V1;
    Array<OneD, int> V2;   

    GenerateMapEidsv1v2(streak,V1,V2);
    Array<OneD, Array<OneD, int> > lay_Eids;
    Array<OneD, Array<OneD, int> > lay_Vids;
    int nlays=0;
    MappingEVids(xold_up, yold_up, xold_low, yold_low, xold_c, yold_c, Vids_c,
                 graphShPt,streak, V1, V2, nlays,  lay_Eids, lay_Vids);
    
    
    
cout<<"nlays="<<nlays<<endl;
    Array<OneD, Array<OneD, NekDouble> > layers_y(nlays);
    Array<OneD, Array<OneD, NekDouble> > layers_x(nlays);
    //initialise layers_y,lay_eids
    for(int g=0; g<nlays; g++)
    {
        layers_y[g]= Array<OneD, NekDouble> ( (nvertl-1)*npedge );
    }
      
  
    
    /////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    //££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //@todo set Delta0 from session file
    NekDouble Delta0;
    if(vSession->DefinesParameter("Delta"))
    {
          Delta0 = vSession->GetParameter("Delta"); 
    }
    else
    {
          Delta0 = 0.1;//default value
    }

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££
    ////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////  
    //save the coords of the old vertices
    int nVertTot = graphShPt->GetNvertices();
    //arrays to check with the new values
    Array<OneD, NekDouble> xold(nVertTot);
    Array<OneD, NekDouble> yold(nVertTot);
    //calculate the new cordinates of the vertices   
    Array<OneD, NekDouble> xnew(nVertTot);
    Array<OneD, NekDouble> ynew(nVertTot,-20);    
    Array<OneD, int> Up(nvertl);//Vids lay Up
    Array<OneD, int> Down(nvertl);//Vids lay Down   
    int cntup=0;
    int cntlow=0;
    //Vids needed only if a layers has to be moved
    NekDouble bleft=-10;
    NekDouble tright = 10;
    NekDouble bright = -10;
    NekDouble tleft  = 10;
    int cnt=0;
    int bottomleft, topright,bottomright,topleft;
    for(int i=0; i<nVertTot; i++)
    {
         bool mvpoint =false;
       	 NekDouble ratio;  
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(i);
         NekDouble x,y,z;
         vertex->GetCoords(x,y,z); 
         
         xold[i] = x;
         yold[i] = y;
         
         
         //x coord doesn't change
         xnew[i]=x;
//cout<<"x="<<x<<"  y="<<y<<endl;
         //bottom, left (x=0, y<ydown)        
         if(x==0 && y< yold_low[0] 
            && y> bleft)
         {
              bleft = y;
              bottomleft = i;
         }
         //top, right
         if(x== xold_c[nvertl-1] && y> yold_up[nvertl-1]
            && y<tright)
         {
              tright = y;
              topright = i;
         }
         //bottom, right]
         if(x==xold_c[nvertl-1] && y<yold_low[nvertl-1]
            && y> bright)
         {
              bright = y;
              bottomright = i;
         } 
         //top, left
         if(x== 0 && y> yold_up[0]
            && y<tleft)
         {
              tleft = y;
              topleft = i;
         }
         //find the corresponding yold_l and deltaold
         for(int j=0; j<nvertl; j++)
         {	      
             if((xold_up[j]==x)&&(yold_up[j]==y))
             {        	 
                 //ratio = (1-y)*(1+y)/( (1-yold_c[j])*(1+yold_c[j]) );
                 //ynew[i] = y + Sign[j]*Delta_c[j]*ratio;
                 ynew[i] = y_c[j] +Delta0;
                 mvpoint=true;
                 Up[j] = i;
             }
             if((xold_low[j]==x)&&(yold_low[j]==y))
             {
                 //ratio = (1-y)*(1+y)/( (1-yold_c[j])*(1+yold_c[j]) );
                 //ynew[i] = y + Sign[j]*Delta_c[j]*ratio;          	 
                 ynew[i] = y_c[j] -Delta0;
             	 mvpoint=true;
                 Down[j] = i;
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
             if( y>yold_up[qp_closer] && y< 1)
             {	        

                 //ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                 //ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;

                 //ratio = (1-y)*(1+y)/(  (1-yold_up[qp_closer])*(1+yold_up[qp_closer]) );
                 //distance prop to layerup
                 ynew[i] = y_c[qp_closer] +(y-yold_c[qp_closer])*
                   (1-y_c[qp_closer])/(1-yold_c[qp_closer]);
//cout<<"upper zone y="<<y<<"  ratio="<<ratio<<endl;         

           
            
             }
             else if(y<yold_low[qp_closer] && y> -1)
             {

                 //ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                 //ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;              

                 //ratio = (1-y)*(1+y)/(  (1-yold_low[qp_closer])*(1+yold_low[qp_closer]) );
                 //distance prop to layerlow
                 ynew[i] = y_c[qp_closer] + (y-yold_c[qp_closer] )*
                         (-1-y_c[qp_closer])/(-1-yold_c[qp_closer]);
               
             }

             else if ( y>yold_c[qp_closer] && y < yold_up[qp_closer])
             {
                 if(x==0){ cntup++; }
                 //ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                 //ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;
                
             }
             else if (y<yold_c[qp_closer] && y > yold_low[qp_closer])
             {
                if(x==0){ cntlow++; }
                // ratio = (1-y)*(1+y)/( (1-yold_c[qp_closer])*(1+yold_c[qp_closer]) );
                // ynew[i] = y + Sign[qp_closer]*Delta_c[qp_closer]*ratio;
           
             }
             else if( y==1 || y==-1)//bcs don't move
             {
                ynew[i] =y;
             	//cout<<"point x="<<xnew[i]<<"  y="<<y<<"  closer x="<<xold_up[qp_closer]<<endl;   
             }

             //internal layers are not moved yet so...
             if( (ynew[i]>1 || ynew[i]<-1)
                 && ( y>yold_up[qp_closer]  || y<yold_low[qp_closer]) )
             {
             	cout<<"point x="<<xnew[i]<<"  y="<<y<<"  closer x="<<xold_up[qp_closer]<<" ynew="<<ynew[i]<<endl;              	     
             	ASSERTL0(false, "shifting out of range");
             }

         }         
             
    }




    int nqedge = streak->GetExp(0)->GetNumPoints(0);    
    int nquad_lay = (nvertl-1)*nqedge;
    Array<OneD, NekDouble> coeffstmp(Cont_y->GetNcoeffs(),0.0);
    // curve the edges around the NEW critical layer (bool to turn on/off)
    bool curv_lay=true;
    bool move_norm=true;
    int np_lay = (nvertl-1)*npedge;//nedges*npedge (Eq. Points!!!)

    Array<OneD, NekDouble> xcQ(nqedge*nedges,0.0);
    Array<OneD, NekDouble> ycQ(nqedge*nedges,0.0);
    Array<OneD, NekDouble> zcQ(nqedge*nedges,0.0);  
    Array<OneD, NekDouble> nxPhys(npedge*nedges,0.0);
    Array<OneD, NekDouble> nyPhys(npedge*nedges,0.0);  
    Array<OneD, NekDouble> nxQ(nqedge*nedges,0.0);
    Array<OneD, NekDouble> nyQ(nqedge*nedges,0.0);   
    
    if( move_norm==true )
    {
       //np_lay = (nvertl-1)*nqedge;//nedges*nqedge   (Q points!!!)
       //extract crit lay normals (through tangents):
       Array<OneD, NekDouble> xcPhysMOD(xcPhys.num_elements());
       Array<OneD, NekDouble> ycPhysMOD(ycPhys.num_elements());
       //copy(temporary the xcPhys,ycPhys into xcPhysMOD, ycPhysMOD)
       Vmath::Vcopy(xcPhys.num_elements(),xcPhys,1,xcPhysMOD,1);
       Vmath::Vcopy(xcPhys.num_elements(),ycPhys,1,ycPhysMOD,1);

       Array<OneD, StdRegions::StdExpansion1DSharedPtr> Edge_newcoords(2);   

cout<<"nquad per edge="<<nqedge<<endl;
       for(int l=0; l<2; l++)
       {
           Edge_newcoords[l] = boost::dynamic_pointer_cast<StdRegions::StdExpansion1D>
                (bndfieldx[lastIregion]->GetExp(0));
       }   
       Array<OneD, NekDouble> xnull(nqedge);
       Array<OneD, NekDouble> ynull(nqedge);
       Array<OneD, NekDouble> xcedgeQ(nqedge,0.0);
       Array<OneD, NekDouble> ycedgeQ(nqedge,0.0);
       Array<OneD, NekDouble> txedgeQ(nqedge,0.0);
       Array<OneD, NekDouble> tyedgeQ(nqedge,0.0);
       Array<OneD, NekDouble> normsQ(nqedge,0.0); 

       Array<OneD, NekDouble> txQ(nqedge*nedges,0.0);
       Array<OneD, NekDouble> tyQ(nqedge*nedges,0.0); 
       Array<OneD, NekDouble> tx_tyedgeQ(nqedge,0.0);    
       Array<OneD, NekDouble> nxedgeQ(nqedge,0.0);
       Array<OneD, NekDouble> nyedgeQ(nqedge,0.0);
       Array<OneD, const NekDouble> ntempQ(nqedge) ; 

       Array<OneD, NekDouble> nxedgePhys(npedge,0.0);
       Array<OneD, NekDouble> nyedgePhys(npedge,0.0);


       Array<OneD, NekDouble> CoordsPhys(2);


       bndfieldx[lastIregion]->GetCoords(xcQ, ycQ, zcQ);
       //determine the NEW crit lay quad points values(lagrangeInterpolant):
       //interp(from xcPhys, ycPhys).
       Array<OneD, NekDouble>x_tmp(np_lay-(nedges-1),0.0);
       Array<OneD, NekDouble>y_tmp(np_lay-(nedges-1),0.0);
       Array<OneD, NekDouble>closex(4,0.0);
       Array<OneD, NekDouble>closey(4,0.0);
       Cutrepetitions(nedges, xcPhysMOD,x_tmp);
       Cutrepetitions(nedges, ycPhysMOD,y_tmp);
       for(int l=0; l< xcQ.num_elements(); l++)
       {
	    Cutrepetitions(nedges, xcPhysMOD,x_tmp);
	    Cutrepetitions(nedges, ycPhysMOD,y_tmp);
            int indclose = DetermineclosePointxindex( xcQ[l], x_tmp);
            
            //generate neighbour arrays
            GenerateNeighbourArrays(indclose, 4,x_tmp,y_tmp,closex,closey);

            ycQ[l]=  LagrangeInterpolant(
                          xcQ[l],4,closex,closey  );



       }
       
       
       
       //force continuity

       Vmath::Vcopy(nquad_lay, ycQ,1,  Cont_y->UpdatePhys(),1);
       Array<OneD, NekDouble> coeffsy(Cont_y->GetNcoeffs(),0.0);

       Cont_y->FwdTrans_IterPerExp(Cont_y->GetPhys(), coeffsy);
       Cont_y->BwdTrans_IterPerExp( coeffsy, Cont_y->UpdatePhys()); 
       Vmath::Vcopy(nquad_lay, Cont_y->GetPhys(),1, ycQ,1);    

cout<<"xcQ, ycQ"<<endl;
for(int s=0; s<xcQ.num_elements(); s++)
{
cout<<xcQ[s]<<"     "<<ycQ[s]<<endl;
}         
//ASSERTL0(false, "dsdfs");        
       bool evaluatetan=false;
       NekDouble incratio;
       Array<OneD, int> edgeinterp(2);
       int wrong=0;
       for(int k=0; k<nedges; k++)
       {
            Vmath::Vcopy(nqedge, &xcQ[k*nqedge],1,&xcedgeQ[0],1);
            Vmath::Vcopy(nqedge, &ycQ[k*nqedge],1,&ycedgeQ[0],1);
            //calc the NEW tangent values
            Edge_newcoords[0]->StdPhysDeriv(xcedgeQ,txedgeQ);
            Edge_newcoords[1]->StdPhysDeriv(ycedgeQ,tyedgeQ);
            //norms=tx*tx
            Vmath::Vmul(nqedge,txedgeQ,1,txedgeQ,1,normsQ,1);
            //norms=tx*tx+ty*ty
            Vmath::Vvtvp(nqedge,tyedgeQ,1,tyedgeQ,1,normsQ,1,normsQ,1);
            
            Vmath::Vsqrt(nqedge, normsQ,1,normsQ,1);
            Vmath::Sdiv(nqedge,1.0,normsQ,1,normsQ,1);
            
            Vmath::Vmul(nqedge,txedgeQ,1,normsQ,1,txedgeQ,1);
            Vmath::Vmul(nqedge,tyedgeQ,1,normsQ,1,tyedgeQ,1);            
	        
            //try evaluate tangent if the incremental ratio is high





            evaluatetan=false;
            for(int u=0; u<nqedge-1; u++)
            {
               incratio   = (ycedgeQ[u+1]- ycedgeQ[u])/(xcedgeQ[u+1]- xcedgeQ[u]);                
cout<<"incratio="<<incratio<<endl;
               if(abs(incratio)> 4.0 && evaluatetan==false )
               {
cout<<"wrong="<<wrong<<endl;
                    evaluatetan=true;
                    ASSERTL0(wrong<2, "number edges to change is too high!!");
                    edgeinterp[wrong]=k;
                    wrong++;
               }
            }  
            
            if(evaluatetan)
            {
cout<<"tan bef"<<endl;
for(int e=0; e< nqedge; e++)
{
cout<<xcedgeQ[e]<<"     "<<ycedgeQ[e]<<"      "<<txedgeQ[e]<<endl;
}

                 //OR: fit
                 int polyorder =3;
                 Array<OneD, NekDouble>  coeffsinterp(polyorder+1);      
                 Array<OneD, NekDouble> yPedges(npedge,0.0);
                 Array<OneD, NekDouble> xPedges(npedge,0.0);
                 Vmath::Vcopy(npedge, &xcPhysMOD[k*npedge+0],1,&xPedges[0],1);
                 Vmath::Vcopy(npedge, &ycPhysMOD[k*npedge+0],1,&yPedges[0],1);

                 PolyFit(polyorder,nqedge, xcedgeQ,ycedgeQ, coeffsinterp, xPedges,yPedges, 
                  npedge);
                 //update values
                 Vmath::Vcopy(npedge, &xPedges[0],1, &xcPhysMOD[k*npedge+0],1);
                 Vmath::Vcopy(npedge, &yPedges[0],1, &ycPhysMOD[k*npedge+0],1);

                 EvaluateTangent(nqedge,xcedgeQ, coeffsinterp,txedgeQ, tyedgeQ);


            }
            //copy also the tx,ty
            Vmath::Vcopy(nqedge, &(txedgeQ[0]), 1, &(txQ[nqedge*k]),1);
            Vmath::Vcopy(nqedge, &(tyedgeQ[0]), 1, &(tyQ[nqedge*k]),1);



        }            
        
       Array<OneD, NekDouble> fz(nedges*nqedge,0.0);
       bndfieldx[lastIregion]->PhysDeriv(MultiRegions::eX,ycQ,fz);
for(int w=0; w< fz.num_elements(); w++)
{
       txQ[w] = cos(atan(fz[w]));
       tyQ[w] = sin(atan(fz[w]));
cout<<xcQ[w]<<"    "<<ycQ[w]<<"       "<<fz[w]<<endl;
}
//ASSERTL0(false, "bobo");
       //force continuity tx,ty
       //tx
       Vmath::Vcopy(nquad_lay, txQ,1, Cont_y->UpdatePhys(),1);
       Cont_y->FwdTrans_IterPerExp(Cont_y->GetPhys(), coeffsy);
       Cont_y->BwdTrans_IterPerExp( coeffsy, Cont_y->UpdatePhys()); 
       Vmath::Vcopy(nquad_lay, Cont_y->GetPhys(),1, txQ,1);  
       //ty
       Vmath::Vcopy(nquad_lay, tyQ,1, Cont_y->UpdatePhys(),1);
       Cont_y->FwdTrans_IterPerExp(Cont_y->GetPhys(), coeffsy);
       Cont_y->BwdTrans_IterPerExp( coeffsy, Cont_y->UpdatePhys()); 
       Vmath::Vcopy(nquad_lay, Cont_y->GetPhys(),1, tyQ,1); 

        //check if the tan points have the same curvature otherwise interp           
        NekDouble inc,incbefore;

        //build-up the fit for the tan using the edge with 
        //the same derivative sign (before or after)
        
        int edgebef;
        int edgeaft;        
        
        
        for(int q=0; q<2; q++)
        {
             edgebef = edgeinterp[q]-1;
             edgeaft = edgeinterp[q]+1;
             incbefore = (txQ[edgebef*nqedge+nqedge-1]-txQ[edgebef*nqedge])/
                        (xcQ[edgebef*nqedge+nqedge-1]-xcQ[edgebef*nqedge]);
             inc =  (txQ[edgeinterp[q]*nqedge+nqedge-1]-txQ[edgeinterp[q]*nqedge])/
                        (xcQ[edgeinterp[q]*nqedge+nqedge-1]-xcQ[edgeinterp[q]*nqedge]);
             int npoints = 2*nqedge;

             Array<OneD, NekDouble> yQedges(npoints,0.0);
             Array<OneD, NekDouble> xQedges(npoints,0.0);
             Array<OneD, NekDouble> tyQedges(npoints,0.0);
             Array<OneD, NekDouble> txQedges(npoints,0.0);
cout<<"inc="<<inc<<"    incbef="<<incbefore<<endl;
             if(    (inc/incbefore)>0.           )             
             {
cout<<"before!!"<<edgebef<<endl;
                 int polyorder =2;
                 Array<OneD, NekDouble>  coeffsinterp(polyorder+1); 
                 Vmath::Vcopy(npoints, &xcQ[edgebef*nqedge+0],1,&xQedges[0],1);
                 Vmath::Vcopy(npoints, &ycQ[edgebef*nqedge+0],1,&yQedges[0],1);
                 Vmath::Vcopy(npoints, &txQ[edgebef*nqedge+0],1,&txQedges[0],1);
                 Vmath::Vcopy(npoints, &tyQ[edgebef*nqedge+0],1,&tyQedges[0],1);
                      
                 PolyFit(polyorder, npoints,
                   xQedges,txQedges, 
                   coeffsinterp, xQedges,txQedges, npoints);

                 //copy back the values:
                 Vmath::Vcopy(npoints,&txQedges[0],1, &txQ[edgebef*nqedge+0],1);   

                 PolyFit(polyorder, npoints,
                   xQedges,tyQedges, 
                   coeffsinterp, xQedges,tyQedges, npoints);

                 //copy back the values:
                 Vmath::Vcopy(npoints,&tyQedges[0],1, &tyQ[edgebef*nqedge+0],1);   

             }
             else
             {
cout<<"after!!"<<endl;
                 int polyorder =2;
                 Array<OneD, NekDouble>  coeffsinterp(polyorder+1); 
                 Vmath::Vcopy(npoints, &xcQ[edgeinterp[q]*nqedge+0],1,&xQedges[0],1);
                 Vmath::Vcopy(npoints, &ycQ[edgeinterp[q]*nqedge+0],1,&yQedges[0],1);
                 Vmath::Vcopy(npoints, &txQ[edgeinterp[q]*nqedge+0],1,&txQedges[0],1);
                 Vmath::Vcopy(npoints, &tyQ[edgeinterp[q]*nqedge+0],1,&tyQedges[0],1);

                      
                 PolyFit(polyorder, npoints,
                   xQedges,txQedges, 
                   coeffsinterp, xQedges,txQedges, npoints);

                 //copy back the values:
                 Vmath::Vcopy(npoints,&txQedges[0],1, &txQ[edgeinterp[q]*nqedge+0],1);   

                 PolyFit(polyorder, npoints,
                   xQedges,tyQedges, 
                   coeffsinterp, xQedges,tyQedges, npoints);

                 //copy back the values:
                 Vmath::Vcopy(npoints,&tyQedges[0],1, &tyQ[edgeinterp[q]*nqedge+0],1); 

            
             }


        }
        //force continuity of the tangent
        //tyQ
        Vmath::Vcopy(nquad_lay, tyQ,1,  Cont_y->UpdatePhys(),1);        
        Cont_y->FwdTrans_IterPerExp(Cont_y->GetPhys(), coeffstmp);
        Cont_y->BwdTrans_IterPerExp( coeffstmp, Cont_y->UpdatePhys()); 
        Vmath::Vcopy(nquad_lay, Cont_y->GetPhys(),1, tyQ,1);    
        //txQ
        Vmath::Vcopy(nquad_lay, txQ,1,  Cont_y->UpdatePhys(),1);        
        Cont_y->FwdTrans_IterPerExp(Cont_y->GetPhys(), coeffstmp);
        Cont_y->BwdTrans_IterPerExp( coeffstmp, Cont_y->UpdatePhys()); 
        Vmath::Vcopy(nquad_lay, Cont_y->GetPhys(),1, txQ,1);         
        
        for(int k=0; k<nedges; k++)
        {
            //determine the normal from eqs(APART FROM SIGN):
            //tx^2 +ty^2= 1 = nx^2 + ny^2;  
            //t\cdot n=0= tx*nx +ty*ny
            //result: nx = ( 1+(tx/ty)^2 )^(-1/2)


            Vmath::Vcopy(nqedge, &(txQ[nqedge*k]),1, &(txedgeQ[0]), 1);
            Vmath::Vcopy(nqedge, &(tyQ[nqedge*k]),1, &(tyedgeQ[0]), 1);

            Vmath::Vdiv(nqedge, txedgeQ,1,tyedgeQ,1,tx_tyedgeQ,1);
            Vmath::Vmul(nqedge, tx_tyedgeQ,1,tx_tyedgeQ,1,tx_tyedgeQ,1);
            Vmath::Sadd(nqedge,1.0,tx_tyedgeQ,1,nxedgeQ,1);
            Vmath::Vsqrt(nqedge,nxedgeQ,1,nxedgeQ,1);
            Vmath::Sdiv(nqedge,1.0,nxedgeQ,1,nxedgeQ,1);
            //normal DOWNWARDS!!! mult by -1
            Vmath::Smul(nqedge, -1.0,nxedgeQ,1,nxedgeQ,1);
            Vmath::Vcopy(nqedge, &(nxedgeQ[0]),1, &(nxQ[nqedge*k]),1);
            //ny = (1-nx ^2)^(1/2)
            Vmath::Vmul(nqedge, nxedgeQ,1,nxedgeQ,1,nyedgeQ,1);
            Vmath::Smul(nqedge, -1.0,nyedgeQ,1,nyedgeQ,1);
            Vmath::Sadd(nqedge,1.0,nyedgeQ,1,nyedgeQ,1);
            Vmath::Vsqrt(nqedge,nyedgeQ,1,nyedgeQ,1);
            //normal DOWNWARDS!!! mult by -1
            Vmath::Smul(nqedge, -1.0,nyedgeQ,1,nyedgeQ,1);
            Vmath::Vcopy(nqedge, &(nyedgeQ[0]), 1, &(nyQ[nqedge*k]),1);




cout<<"edge:"<<k<<endl;
cout<<"tan/normal"<<endl;
for(int r=0; r<nqedge; r++)
{
cout<<xcQ[k*nqedge+r]<<"     "<<txedgeQ[r]<<"      "<<tyedgeQ[r]<<"    "
<<nxedgeQ[r]<<"      "<<nyedgeQ[r]<<endl;
}


       }  
       
       //force continuity:
       //REMEMBER: the Fwd/Bwd operation get wrong with the border values!!!
       Vmath::Vcopy(nquad_lay, nyQ,1,  Cont_y->UpdatePhys(),1);
        
       Cont_y->FwdTrans_IterPerExp(Cont_y->GetPhys(), coeffstmp);
       Cont_y->BwdTrans_IterPerExp( coeffstmp, Cont_y->UpdatePhys()); 
       Vmath::Vcopy(nquad_lay, Cont_y->GetPhys(),1, nyQ,1);        
       
       Vmath::Zero(nquad_lay,Cont_y->UpdatePhys(),1);
       Vmath::Zero(Cont_y->GetNcoeffs(),Cont_y->UpdateCoeffs(),1);
       Vmath::Vcopy(nquad_lay, nxQ,1,  Cont_y->UpdatePhys(),1);      
       Cont_y->FwdTrans_IterPerExp(Cont_y->GetPhys(), coeffstmp);
       Cont_y->BwdTrans_IterPerExp( coeffstmp, Cont_y->UpdatePhys()); 
       Vmath::Vcopy(nquad_lay, Cont_y->GetPhys(),1, nxQ,1);   

       //force the normal at interface point to be equal
       for(int k=0; k<nedges; k++)
       {
             if(k>0)
             {
                  //nyPhys[f*npedge +0] = 
                  //          (nyPhys[(f-1)*npedge+npedge-1]+nyPhys[f*npedge+0])/2.;
                  nyQ[(k-1)*nqedge+nqedge-1]=
                            nyQ[k*nqedge+0];
                  //nx= (1-ny^2)^{1/2}
                  //nxPhys[f*npedge+0]= 
                  //           sqrt(1- nyPhys[f*npedge+0]*nyPhys[f*npedge+0]);
                  nxQ[(k-1)*nqedge+nqedge-1]=
                            nxQ[k*nqedge+0];
              }
       }

       Array<OneD, NekDouble>x_tmpQ(nquad_lay-(nedges-1));
       //Array<OneD, NekDouble>tmpnxQ(nquad_lay-(nedges-1));
       Array<OneD, NekDouble>tmpnyQ(nquad_lay-(nedges-1));    
       
cout<<"nx,yQbefore"<<endl;
for(int u=0; u<xcQ.num_elements(); u++)
{
cout<<xcQ[u]<<"      "<<nyQ[u]<<"     "<<txQ[u]<<endl;
}

       Cutrepetitions(nedges, xcQ,x_tmpQ);
       //Cutrepetitions(nedges, nxQ, tmpnxQ);
       Cutrepetitions(nedges, nyQ, tmpnyQ);       
cout<<"nx,yQ"<<endl;
for(int u=0; u<x_tmpQ.num_elements(); u++)
{
cout<<x_tmpQ[u]<<"      "<<tmpnyQ[u]<<endl;
}

       
       
       //interpolate the values into phys points(curved points)
       int segid,offset;       
       for(int k=0; k<nedges; k++)
       {
//cout<<"edge:"<<k<<endl;       	       
            for(int a=0; a<npedge; a++)
            {
                if(a==0)//verts pos no interp necessary
                {
                     nxPhys[k*npedge +a]= nxQ[k*nqedge +0];
                     nyPhys[k*npedge +a]= nyQ[k*nqedge +0];

                }
                else if(a== npedge-1)//verts pos no interp necessary
                {
                     nxPhys[k*npedge +a]= nxQ[k*nqedge +nqedge-1];
                     nyPhys[k*npedge +a]= nyQ[k*nqedge +nqedge-1];
//cout<<":last"<<nyQ[k*nqedge+a]<<endl;                     
                     
                }

                else
                {
                     //use lagrange interpolant to get the 
                     //normal at phys(equispaced points)

                     //order normal functions(cut out repetitions)
                     //QUAD POINTS


                     int index; 
                     //determine closest index:
             
                     index= 
                         DetermineclosePointxindex( Addpointsx[k*(npedge-2) +a-1], x_tmpQ);

                     Array<OneD, NekDouble> Pxinterp(4);
                     Array<OneD, NekDouble> Pyinterp(4);

                     //generate neighbour arrays (y):                 
                     GenerateNeighbourArrays(index, 4,x_tmpQ,tmpnyQ,Pxinterp,Pyinterp);
                     //interp the new normal components(y)                     
                     nyPhys[k*npedge +a]=
                      LagrangeInterpolant(Addpointsx[k*(npedge-2) +a-1],4,Pxinterp,Pyinterp  );
/*                      
                     //generate neighbour arrays (x):
                     GenerateNeighbourArrays(index,4,x_tmpQ,tmpnxQ,Pxinterp,Pyinterp);
                     //interp the new normal components(x)                     
                     nxPhys[k*npedge +a]=                     
                      LagrangeInterpolant(Addpointsx[k*(npedge-2) +a],4,Pxinterp,Pyinterp  );
*/
                     //nx=-(1-ny*ny){1/2} the normal is DOWNWARDS!!!  
                     nxPhys[k*npedge +a]= -sqrt(abs(1- nyPhys[k*npedge +a]*nyPhys[k*npedge +a]));
/*
                     //(put the middle points as quad points)
                     //GaussLobattoLegendre points in the middle:
                     nxPhys[k*npedge +a] = nxedgeQ[a];
                     nyPhys[k*npedge +a] = nyedgeQ[a];
                     ASSERTL0(npedge< nqedge," quad points too low");
*/                     
                }


                //force the normal at interface point to be equal
                if(k>0)
                {
                     //nyPhys[f*npedge +0] = 
                     //          (nyPhys[(f-1)*npedge+npedge-1]+nyPhys[f*npedge+0])/2.;
                     nyPhys[(k-1)*npedge+npedge-1]=
                            nyPhys[k*npedge+0];
                     //nx= (1-ny^2)^{1/2}
                     //nxPhys[f*npedge+0]= 
                     //           sqrt(1- nyPhys[f*npedge+0]*nyPhys[f*npedge+0]);
                     nxPhys[(k-1)*npedge+npedge-1]=
                            nxPhys[k*npedge+0];
                }
                 
     

            }


      


        }
cout<<"xcPhys,,"<<endl;    
for(int s=0; s<np_lay; s++)
{	

cout<<xcPhysMOD[s]<<"     "<<ycPhysMOD[s]<<"     "<<nxPhys[s]<<"     "<<nyPhys[s]<<endl;

}          

     //determine the new coords of the vertices and the curve points 
     //for each edge
       
     //int np_lay = (nvertl-1)*npedge;//nedges*npedge
       
     //NB delta=ynew-y_c DEPENDS ON the coord trnaf ynew= y+Delta_c*ratio!!!
     Array<OneD, NekDouble> delta(nlays);
     Array<OneD, NekDouble>tmpy_lay(np_lay);
     Array<OneD, NekDouble>tmpx_lay(np_lay);
     for(int m=0; m<nlays; m++)
     {
         //delta[m] = (ynew[lay_Vids[m][0]] - y_c[0])/1.0;
         
         //depends on Delta0 
         if(m< cntlow+1)
         {
              delta[m]  =  -(cntlow+1-m)*Delta0/(cntlow+1);
         }
         else
         {
              delta[m]  = (  m-(cntlow)  )*Delta0/(cntlow+1); 
         }
                            
                            
         layers_x[m]= Array<OneD, NekDouble>(np_lay);
//cout<<"delta="<<delta[m]<<" cntlow="<<cntlow<<endl;

         for(int h=0; h< nvertl; h++)
         {
             //shift using the dinstance delta from the crit layer AT x=0
             //for each layer

  
//cout<<m<<"Vid:"<<lay_Vids[m][h]<<"  mod from y="<<ynew[lay_Vids[m][h] ]<<"  to y="<<y_c[h]    +delta[m]<<endl;     
             if(move_norm==false)
             {
                 ynew[lay_Vids[m][h] ]= y_c[h] +delta[m];      
                 xnew[lay_Vids[m][h] ]= x_c[h];  
             }
             else
             {
                 if(h==0 || h==nvertl-1 )//nx=0,ny=1 at the borders
                 { 
                     ynew[lay_Vids[m][h] ]= y_c[h] +delta[m];   
                     xnew[lay_Vids[m][h] ]= x_c[h];  
                 }
                 else
                 {
                     ynew[lay_Vids[m][h] ]= y_c[h] +delta[m]*abs(nyPhys[h*npedge+0]);   
                     xnew[lay_Vids[m][h] ]= x_c[h] +delta[m]*abs(nxPhys[h*npedge+0]); 
                 }
             }
cout<<"Vid x="<<xnew[lay_Vids[m][h] ]<<"   y="<<ynew[lay_Vids[m][h] ]<<endl;
             if(h< nedges 
                 //&& curv_lay==true
               )
             {
cout<<"edge=="<<h<<endl;             	     
                 if(h>0)//check normal consistency
                 {
                 ASSERTL0( nyPhys[h*npedge+0]==nyPhys[(h-1)*npedge+npedge-1]," normaly wrong");
                 ASSERTL0( nxPhys[h*npedge+0]==nxPhys[(h-1)*npedge+npedge-1]," normalx wrong");

                 }
                 if(move_norm==false)
                 {
                     //v1
                     layers_y[m][h*npedge +0] = y_c[h] +delta[m];
                     layers_x[m][h*npedge +0] = xnew[lay_Vids[m][h] ];
                     //v2             
                     layers_y[m][h*npedge +npedge-1] = y_c[h+1] +delta[m]; 
                     layers_x[m][h*npedge +npedge-1] = xnew[lay_Vids[m][h+1] ];
                     //middle points (shift crit lay points by delta):
                     for(int d=0; d< npedge-2; d++)
                     {
                         layers_y[m][h*npedge +d+1]=  ycPhysMOD[h*npedge +d+1] +delta[m];
                                  //Addpointsy[h*(npedge-2) +d] +delta[m];    
                         layers_x[m][h*npedge +d+1]=  xcPhysMOD[h*npedge +d+1];
                                  //Addpointsx[h*(npedge-2) +d];
                     }
                 }
                 else
                 {
                     if(h==0) //nx=0,ny=1 at the borders
                     {
                         //v1
                         tmpy_lay[h*npedge +0] = y_c[h] +delta[m];
                         tmpx_lay[h*npedge +0] = xnew[lay_Vids[m][h] ];
                         //v2             
                         tmpy_lay[h*npedge +npedge-1] =
                                        y_c[h+1] +delta[m]*abs(nyPhys[h*npedge +npedge-1]);
                         tmpx_lay[h*npedge +npedge-1] = 
                                        x_c[h+1] +delta[m]*abs(nxPhys[h*npedge +npedge-1]);
                     }
                     else if(h==nedges-1)//nx=0,ny=1 at the borders
                     {
                         //v1
                         tmpy_lay[h*npedge +0] = 
                                        y_c[h] +delta[m]*abs(nyPhys[h*npedge +0]);
                         tmpx_lay[h*npedge +0] = 
                                        x_c[h] +delta[m]*abs(nxPhys[h*npedge +0]);
                         //v2             
                         tmpy_lay[h*npedge +npedge-1] = y_c[h+1] +delta[m];
                         tmpx_lay[h*npedge +npedge-1] = xnew[lay_Vids[m][h+1] ];
                     }
                     else
                     {
                         //v1
                         tmpy_lay[h*npedge +0] = 
                                        y_c[h] +delta[m]*abs(nyPhys[h*npedge +0]);
                         tmpx_lay[h*npedge +0] = 
                                        x_c[h] +delta[m]*abs(nxPhys[h*npedge +0]);
                         //v2             
                         tmpy_lay[h*npedge +npedge-1] = 
                                        y_c[h+1] +delta[m]*abs(nyPhys[h*npedge +npedge-1]);
                         tmpx_lay[h*npedge +npedge-1] =
                                        x_c[h+1] +delta[m]*abs(nxPhys[h*npedge +npedge-1]); 
                     }

                     //middle points 
                     for(int d=0; d< npedge-2; d++)
                     {

                         tmpy_lay[h*npedge +d+1] = ycPhysMOD[h*npedge +d+1] +
                                  delta[m]*abs(nyPhys[h*npedge +d+1]);  
                                //Addpointsy[h*(npedge-2) +d] +
                                //  delta[m]*abs(nyPhys[h*npedge +d+1]);  
                         tmpx_lay[h*npedge +d+1]=  xcPhysMOD[h*npedge +d+1] +
                                  delta[m]*abs(nxPhys[h*npedge +d+1]);  
                                //Addpointsx[h*(npedge-2) +d] +
                                //  delta[m]*abs(nxPhys[h*npedge +d+1]);  

                         //NB ycQ,xcQ refers to nqedge NOT npedge!!!
                         //tmpy_lay[h*npedge +d+1] =  ycQ[h*nqedge +d+1] +
                           //       delta[m]*abs(nyPhys[h*npedge +d+1]);  
                         //tmpx_lay[h*npedge +d+1]=  xcQ[h*nqedge +d+1] +
                          //        delta[m]*abs(nxPhys[h*npedge +d+1]);  
//cout<<"xmoved="<<tmpx_lay[h*npedge +d+1]<<"  xold="<<xcQ[h*nqedge +d+1]<<endl;
                         //ASSERTL0(tmpx_lay[h*npedge +d+1]>0," middle point with x<0")

                     }

                 }
                 
                 
                 
             }//close edges
         }//close verts h                 

for(int s=0; s<np_lay; s++)
{	

cout<<tmpx_lay[s]<<"     "<<tmpy_lay[s]<<endl;

}
         Orderfunctionx(tmpx_lay,tmpy_lay, tmpx_lay,  tmpy_lay);      
cout<<"fisrt interp"<<endl;
for(int s=0; s<np_lay; s++)
{	

cout<<tmpx_lay[s]<<"     "<<tmpy_lay[s]<<endl;

}

         
         
         



         
//ASSERTL0(false, "dasd");



                 //ASSERTL0(tmpx_lay[h*npedge +0]>=0," vert 0 x<0");
                 //ASSERTL0(tmpx_lay[h*npedge +npedge-1]>0," vert 1 x<0");


   



         //check if the x coord is 'outofbound' and calculate the 
         //number of outofbound points
         
         //determine which boudn has been overcome:
         NekDouble boundleft = xcPhysMOD[0];
         NekDouble boundright = xcPhysMOD[np_lay-1];
         bool outboundleft= false;
         bool outboundright=false;
         if(tmpx_lay[1]< boundleft )
         {  
              outboundleft = true;
         }
         if(tmpx_lay[np_lay-2] > boundright )
         {
              outboundright = true;
         }
         
         
         
         int outvert=0;
         int outmiddle=0;
         int outcount=0;         

         Array<OneD, int> vertout(nvertl);                    
         for(int r=0; r< nedges; r++)
         {
              //check point outofboundleft
              if(tmpx_lay[r*npedge + npedge-1]< boundleft && outboundleft==true )//assume the neg coords start from 0
              {
                    vertout[outvert]=r;                     	     
                    outvert++;
              
                    if(r<nedges-1 )
                    {
                         //check if after the last negvert there are neg points  
                    	 if( tmpx_lay[(r+1)*npedge + npedge-1]> boundleft )
                    	 {
                           	   
                              for(int s=0; s<npedge-2; s++)
                              {
                                   if(tmpx_lay[(r+1)*npedge + s+1]<  boundleft)
                                   {
                                        outmiddle++;	     
                                   }
                              }
                                	
                         }                          
                    }
              }
              //check point outofboundright
              if(tmpx_lay[r*npedge + 0]> boundright && outboundright==true )//assume the neg coords start from 0
              {
                   vertout[outvert]=r;     
                   outvert++;
              
                   if( r> 0)
                   {
              	        //check if after the last outvert there are out points  
                        if( tmpx_lay[(r-1)*npedge + 0]< boundright )
                        {
                           	   
                             for(int s=0; s<npedge-2; s++)
                             {
                                  if(tmpx_lay[(r-1)*npedge + s+1]>  boundright)
                                  {                               	  
                                       outmiddle++;	     
                                  }
                             }
                                	
                        }                          
                   }   
              }
         }

         //calc number of point to replace
         outcount = outvert*npedge+1+ outmiddle;
         //determine from(until) which index the replacement will start 
         int replacepointsfromindex=0;
         for(int c=0; c<nedges; c++)
         {
              //assume at least 1 middle point per edge         	 
              if(xcPhysMOD[c*npedge+npedge-1] <= tmpx_lay[c*(npedge-(npedge-2)) +2] && outboundright==true)
              {
              	    replacepointsfromindex =   c*(npedge-(npedge-2))+2;
                    break;
              }

              
              //assume at least 1 middle point per edge         	 
              if(xcPhysMOD[(nedges-1 -c)*npedge+0] >= tmpx_lay[np_lay-1 -(c*(npedge-(npedge-2)) +2)] 
              	      && outboundleft==true)
              {
              	    replacepointsfromindex =   np_lay-1 -(c*(npedge-(npedge-2)) +2);
                    break;
              }
            
              
         }         
         
         
         

//cout<<"out="<<outcount<<endl;                 
//cout<<"replacefrom="<<replacepointsfromindex<<endl;                 
               
                 
	 //if xcoord is neg find the first positive xcoord 


	 if(outcount>1)
	 {
              //determine x new coords:
              //distribute the point all over the layer
              int pstart,shift;
              NekDouble increment;

	      if( outboundright==true)
	      {
	      	   pstart = replacepointsfromindex;
	      	   shift = np_lay-outcount;	      	      
	      	   increment =  (xcPhysMOD[np_lay-outcount]-xcPhysMOD[pstart])/(outcount+1);
                   outcount = outcount-1;  	      	   
	      	   ASSERTL0(tmpx_lay[np_lay-outcount]>xcPhysMOD[(nedges-1)*npedge+0], "no middle points in the last edge");
	      }
	      else
	      {
                   shift=1;
                   pstart= outcount-1;
                   increment = (xcPhysMOD[replacepointsfromindex]-xcPhysMOD[pstart])/(outcount+1);  
                   ASSERTL0(tmpx_lay[pstart]<xcPhysMOD[0*npedge +npedge-1], "no middle points in the first edge");                 
 	      }
	      
	      //interp to points between  posindex and posindex-1
	      Array<OneD, NekDouble> replace_x(outcount);
	      Array<OneD, NekDouble> replace_y(outcount);
	      //order normal functions(cut out repetitions)
	      Array<OneD, NekDouble>x_tmp(np_lay-(nedges-1));
	      Array<OneD, NekDouble>y_tmp(np_lay-(nedges-1));
	      Array<OneD, NekDouble>tmpny(np_lay-(nedges-1));
	      Cutrepetitions(nedges, xcPhysMOD,x_tmp);
	      Cutrepetitions(nedges, ycPhysMOD,y_tmp);
              Cutrepetitions(nedges, nyPhys, tmpny);   
	      //init neigh arrays
      	      Array<OneD, NekDouble>closex(4);
	      Array<OneD, NekDouble>closey(4);                     
              Array<OneD, NekDouble>closeny(4); 
                     NekDouble xctmp,ycinterp,nxinterp,nyinterp;



                     for(int v=0; v<outcount;v++)
                     {
                          xctmp = xcPhysMOD[pstart]+(v+1)*increment;


  
                          //determine closest point index:
                          int index = 
                          DetermineclosePointxindex( xctmp, x_tmp);
//cout<<"  vert="<<index<<endl;


                          //generate neighbour arrays (ny)
                          GenerateNeighbourArrays(index, 4,x_tmp,tmpny,closex,closeny);

                          
                          //interp:
                          nyinterp =
                                 LagrangeInterpolant(
                          xctmp,4,closex,closeny  );
                          
                          //calc nxinterp
                          nxinterp = sqrt(abs(1-nyinterp*nyinterp));
                          
                          //generata neighbour arrays (yc)
                          GenerateNeighbourArrays(index, 4,x_tmp,y_tmp,closex,closey);  
                          //interp:
                          ycinterp  = LagrangeInterpolant(xctmp,4, closex,closey);
                          //calc new coord
                          replace_x[v] = xctmp +delta[m]*abs(nxinterp);
                          replace_y[v] = ycinterp +delta[m]*abs(nyinterp);
                          tmpx_lay[ v+shift ] = replace_x[v];
                          tmpy_lay[ v+shift ] = replace_y[v];                          

//cout<<"xinterp="<<replace_x[v]<<"     yinterp="<<replace_y[v]<<endl;                          

                         
                          

                     }
                 }//end outcount if
                 
                 
                 
/*   
for(int s=0; s<np_lay; s++)
{	

cout<<tmpx_lay[s]<<"     "<<tmpy_lay[s]<<endl;

}                 
if(m== 0)
{	
//ASSERTL0(false, "ssa"); 
}
*/







         int closepoints = 4;

         Array<OneD, NekDouble> Pyinterp(closepoints);
         Array<OneD, NekDouble> Pxinterp(closepoints);

         //check if any edge has less than npedge points
         int edgecount=0;
         int pointscount=0;
         for(int q=0; q<np_lay; q++)
         {
              for(int e=0; e<nedges; e++)
              {
                  if(tmpx_lay[q]<= x_c[e+1]  && tmpx_lay[q]>= x_c[e])
                  {
                      pointscount++;                      
                  }
                  if(q == e*npedge +npedge-1 && pointscount!=npedge )
                  {
//cout<<"edge with few points :"<<e<<endl;                     
                     pointscount=0;
                  }
                  else if(q == e*npedge +npedge-1)
                  {
                     pointscount=0;
                  }
              }
         }
         //----------------------------------------------------------
/*
cout<<"notordered"<<endl;
for(int g=0; g<tmpx_lay.num_elements(); g++)
{
cout<<tmpx_lay[g]<<"    "<<tmpy_lay[g]<<endl;
}         
*/


//cout<<nedges<<"nedges"<<npedge<<" np_lay="<<np_lay<<endl;

         
         
         
         //calc lay coords
         //MoveLayerNfixedxpos(nvertl, npedge, xcPhysMOD, tmpx_lay, tmpy_lay,
         //	 lay_Vids[m], layers_x[m], layers_y[m],xnew,ynew);
         MoveLayerNnormpos(nvertl, npedge, xcPhysMOD, tmpx_lay, tmpy_lay,
         	 lay_Vids[m], layers_x[m], layers_y[m],xnew,ynew);         
       
/*  
         //generate x,y arrays without lastedgepoint
         //(needed to interp correctly)
         Array<OneD, NekDouble>tmpx(np_lay-(nedges-1));
         Array<OneD, NekDouble>tmpy(np_lay-(nedges-1));

         Cutrepetitions(nedges, tmpx_lay, tmpx);
         Cutrepetitions(nedges, tmpy_lay, tmpy);         


         //order points in x:
         int index;
         Array<OneD, NekDouble> copyarray_x(tmpx.num_elements());
         Array<OneD, NekDouble> copyarray_y(tmpx.num_elements());         
         Orderfunctionx(tmpx, tmpy, tmpx, tmpy);        	
        	




         //determine the neighbour points (-3;+3)
         for(int g=0; g< nvertl; g++)
         {
             //verts
//cout<<"determine value for vert x="<<x_c[g]<<endl;
             //determine closest index:
             
             index= 
                 DetermineclosePointxindex( x_c[g], tmpx);
             //generate neighbour arrays:
             GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
             //write vert coords
             ynew[lay_Vids[m][g] ]= LagrangeInterpolant(x_c[g],closepoints,Pxinterp,Pyinterp  );
             xnew[lay_Vids[m][g] ]= x_c[g]; 
             

             if(g<nedges)
             {



             //v1
             layers_y[m][g*npedge +0] = ynew[lay_Vids[m][g] ];
             layers_x[m][g*npedge +0] = xnew[lay_Vids[m][g] ];
             //v2             

             //determine closest index:             
             index= 
                 DetermineclosePointxindex( x_c[g+1], tmpx);
             //generate neighbour arrays:
             GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
             layers_y[m][g*npedge +npedge-1] = 
                       LagrangeInterpolant(x_c[g+1],closepoints,Pxinterp,Pyinterp  );
             layers_x[m][g*npedge +npedge-1] = x_c[g+1];



             //middle points
             for(int r=0; r< npedge-2; r++)
             {

                 //determine closest point index:
                 index = 
                 DetermineclosePointxindex( xcPhysMOD[g*npedge +r+1], tmpx);
//cout<<"  vert+"<<index<<endl;

                 ASSERTL0( index<= tmpy.num_elements()-1, " index wrong");
                 //generate neighbour arrays Pyinterp,Pxinterp
                 GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
                 
                 layers_y[m][g*npedge +r+1]=
                                 LagrangeInterpolant(
                          xcPhysMOD[g*npedge +r+1],closepoints,Pxinterp,Pyinterp  );
//cout<<"x value="<<xcPhysMOD[g*npedge +r+1]<<endl;
                 layers_x[m][g*npedge +r+1]=  xcPhysMOD[g*npedge +r+1];

             }


             }//if edge closed g
         }// nvertl closed
*/         
         //check if there are points out of range:
//cout<<Vmath::Vmax(np_lay,layers_y[m],1)<<endl;
         if(curv_lay==true)
         {
         //ASSERTL0(Vmath::Vmax(np_lay,layers_y[m],1)< Vmath::Vmax(nVertTot,yold,1),"point>ymax");
         //ASSERTL0(Vmath::Vmin(np_lay,layers_y[m],1)> Vmath::Vmin(nVertTot,yold,1),"point<ymin");
         }
   

         //force Polycontinuity of the layer(3 order poly every 2 edges):
         int npoints = npedge;
         Array<OneD, NekDouble>  xPedges(npoints); 
         Array<OneD, NekDouble>  yPedges(npoints); 
         for(int f=0; f<nedges; f++)
         {
             int polyorder=2;

             Array<OneD, NekDouble>  coeffsinterp(polyorder+1); 

             Vmath::Vcopy(npoints, &layers_x[m][(f)*npedge+0],1,&xPedges[0],1);
             Vmath::Vcopy(npoints, &layers_y[m][(f)*npedge+0],1,&yPedges[0],1);


                      
             PolyFit(polyorder, npoints,
                   xPedges,yPedges, 
                   coeffsinterp, xPedges,yPedges, npoints);

             //copy back the values:
             Vmath::Vcopy(npoints,&yPedges[0],1, &layers_y[m][(f)*npedge+0],1);
             
             //verts still the same:
             layers_y[m][f*npedge+0]= ynew[lay_Vids[m][f]];
             layers_y[m][f*npedge+npedge-1]= ynew[lay_Vids[m][f+1]];
             

         }

cout<<" xlay    ylay lay:"<<m<<endl;
for(int l=0; l<np_lay; l++)
{
//cout<<tmpx_lay[l]<<"    "<<tmpy_lay[l]<<endl;
cout<<std::setprecision(8)<<layers_x[m][l]<<"    "<<layers_y[m][l]<<endl;
}
/*
cout<<"nverts"<<endl;
for(int l=0; l<nvertl; l++)
{
cout<<std::setprecision(8)<<xnew[lay_Vids[m][l] ]<<"    "<<ynew[lay_Vids[m][l] ]<<endl;
}
*/

//ASSERTL0(false, "as");



     //if the layers coords are too steep  use two edges verts to get an 
     //poly interp third order
/*
     bool polyinterp=false;
     for(int b=0; b< nedges; b++)
     {
           for(int u=0; u<npedge; u++)
           {
               if(
       abs(layers_y[m][b*npedge+u+1]- layers_y[m][b*npedge +u])/(layers_x[m][b*npedge+u+1]-layers_x[m][b*npedge+u]))>4.0 )
               {
                    polyinterp=true;
                    break;
               }                               
cout<<"incratio="<<incratio<<endl;   
           }
           //
           //Evaluatelay_tan

     }
*/



cout<<"lay="<<m<<endl;
         ASSERTL0(Vmath::Vmin(nVertTot, yold,1)< Vmath::Vmin(np_lay,layers_y[m],1),
             "  different layer ymin val");
         ASSERTL0(Vmath::Vmax(nVertTot, yold,1)> Vmath::Vmax(np_lay,layers_y[m],1),
             "  different layer ymax val");
         ASSERTL0(Vmath::Vmin(nVertTot, xold,1)== Vmath::Vmin(np_lay,layers_x[m],1),
             "  different layer xmin val");
         ASSERTL0(Vmath::Vmax(nVertTot, xold,1)== Vmath::Vmax(np_lay,layers_x[m],1),
             "  different layer xmax val");


     }//close layers!!! m index

     //MoveOutsidePointsfixedxpos(npedge, graphShPt,xold_c, yold_c, xold_low, yold_low,
     //	         xold_up, yold_up, layers_y[0], layers_y[nlays-1], xnew, ynew);
     
     //lastIregion -1 = laydown     
     //lastIregion -2 = layup
     MoveOutsidePointsNnormpos(npedge, graphShPt, xold_c,yold_c, xold_low,yold_low,xold_up,yold_up, 
     layers_x[0], layers_y[0], layers_x[nlays-1], layers_y[nlays-1],nxPhys, nyPhys,xnew, ynew);     
     
    
    
/*
     //update vertices coords outside layers region
     NekDouble ratio;
     for(int n=0; n<nVertTot; n++)
     {
          NekDouble ratio;  
          SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(n);
          NekDouble x,y,z;
          vertex->GetCoords(x,y,z); 
          int qp_closer;
          NekDouble diff;
          int qp_closerup, qp_closerdown;
          NekDouble diffup, diffdown;
          //determine the closer xold_up
          NekDouble tmp=1000;
          diffup =1000;
          diffdown = 1000;
          for(int k=0; k<nvertl; k++)
          {     
              if(abs(x-xold_c[k]) < tmp)
              {
                  tmp = abs(x-xold_c[k]);
                  qp_closer=k;
              }                   
          }
          
          //find nplay_closer
          int nplay_closer;
          if(qp_closer==0)
          {
              nplay_closer=0;//first vert
          }
          else
          {
              nplay_closer= (qp_closer-1)*npedge +npedge-1;
          }


          if(  y>yold_up[qp_closer] && y<1 )//nlays-1 is layerup
          {	        

//              ratio = (1-layers_y[nlays-1][qp_closer])*(1-y_c[qp_closer])/
//                    (  (1-yold_up[n])*(1-yold_c[qp_closer]) );
              ratio = (1-layers_y[nlays-1][nplay_closer])/
                    (  (1-yold_up[qp_closer]) );
              //distance prop to layerup
              ynew[n] = layers_y[nlays-1][nplay_closer] 
                      + (y-yold_up[qp_closer])*ratio;  
              xnew[n] = x;
            
          }
          else if(   y< yold_low[qp_closer]   && y>-1  )//0 is layerdown
          {

              ratio = (1+layers_y[0][nplay_closer])/
                    (  (1+yold_low[qp_closer]) );
              //distance prop to layerlow
              ynew[n] = layers_y[0][nplay_closer] 
                      + (y-yold_low[qp_closer])*ratio;  
              xnew[n] = x;          
          }

     }
*/     

/*
     //update verts coords of critlay(in case EvaluateLayTnaget has been called)
     //it's not necessary !!!
     for(int e=0; e<nvertl; e++)
     {
          ynew[CritLay[e]] = y_c[e];
     }
*/


     }//move_norm bool
     else//move vertically
     {
          MoveLayersvertically(nlays, nvertl, cntlow, cntup,
    	         lay_Vids, x_c, y_c, Down, Up, xnew, ynew, layers_x, layers_y);     	          	     
   	     
     }
     

     //check singular quads:
     CheckSingularQuads(streak, V1, V2,xnew, ynew);     
     //check borders of the new mesh verts:
//cout<<std::setprecision(8)<<"yoldmax="<<Vmath::Vmax(nVertTot, yold,1)<<endl;
//cout<<std::setprecision(8)<<"ynewmax="<<Vmath::Vmax(nVertTot,ynew,1)<<endl;

cout<<std::setprecision(8)<<"xmin="<<Vmath::Vmin(nVertTot, xnew,1)<<endl;
     ASSERTL0(Vmath::Vmin(nVertTot, xold,1)== Vmath::Vmin(nVertTot,xnew,1),
             "  different xmin val");
     ASSERTL0(Vmath::Vmin(nVertTot, yold,1)== Vmath::Vmin(nVertTot,ynew,1),
             "  different ymin val");
     ASSERTL0(Vmath::Vmax(nVertTot, xold,1)== Vmath::Vmax(nVertTot,xnew,1),
             "  different xmax val");
     ASSERTL0(Vmath::Vmax(nVertTot, yold,1)== Vmath::Vmax(nVertTot,ynew,1),
             "  different ymax val");


    
    //replace the vertices with the new ones

    Replacevertices(changefile, xnew , ynew, xcPhys, ycPhys, Eids, npedge, charalp, layers_x,layers_y, lay_Eids, curv_lay);
    

}


void OrderVertices (int nedges, SpatialDomains::MeshGraphSharedPtr graphShPt,
        	MultiRegions::ExpListSharedPtr & bndfield, 
        	Array<OneD, int>& Vids, int v1,int v2, NekDouble x_connect, int & lastedge, 
        	Array<OneD, NekDouble>& x, Array<OneD, NekDouble>& y)
{
      int nvertl = nedges+1;
      int edge;
      Array<OneD, int> Vids_temp(nvertl,-10);      
      NekDouble x0layer =  0.0;
      for(int j=0; j<nedges; j++)
      {
          LocalRegions::SegExpSharedPtr  bndSegExplow = 
          boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndfield->GetExp(j)) ;   	
          edge = (bndSegExplow->GetGeom1D())->GetEid();
//cout<<" edge="<<edge<<endl;   	   
          for(int k=0; k<2; k++)
          {
              Vids_temp[j+k]=(bndSegExplow->GetGeom1D())->GetVid(k);   
              SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(Vids_temp[j+k]);
              NekDouble x1,y1,z1;
              vertex->GetCoords(x1,y1,z1);     	   	   
              if(x1==x_connect && edge!=lastedge)
              {  
              	 //first 2 points   
        	 if(x_connect==x0layer)
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


void Computestreakpositions(int npoints, MultiRegions::ExpListSharedPtr streak,
    	        Array<OneD, NekDouble> xold_up, Array<OneD, NekDouble> yold_up,
    	        Array<OneD, NekDouble> xold_low, Array<OneD, NekDouble> yold_low,
    	        Array<OneD, NekDouble> xold_c, Array<OneD, NekDouble> yold_c,    	        
    	        Array<OneD, NekDouble> &xc,  Array<OneD, NekDouble> &yc, NekDouble cr,
                bool verts)   	        
{
cout<<"Computestreakpositions"<<endl;    
     int nq = streak->GetTotPoints();	
     Array<OneD, NekDouble> coord(2);
     //Array<OneD, NekDouble> stvalues(nvertl,-10);
     Array<OneD, NekDouble> derstreak(nq);
     streak->PhysDeriv(MultiRegions::eY, streak->GetPhys(), derstreak);
     int elmtid, offset;    
     NekDouble U,dU;
     NekDouble F=1000;
           
     if(verts==true)//only for verts makes sense to init the coord values..
     {

           //start guess
           //yc= (yup+ydown)/2
           Vmath::Vadd(xc.num_elements(), yold_up,1,yold_low,1, yc,1);
           Vmath::Smul(xc.num_elements(), 0.5,yc,1,yc,1);
           Vmath::Vcopy(xc.num_elements(),xold_c,1,xc,1);
     }
     else//case of xQ,yQ
     {
           Vmath::Vcopy(xc.num_elements(), xold_c,1,xc,1);
           Vmath::Vcopy(xc.num_elements(), yold_c,1,yc,1);
     }
            
     int its;
     int attempt;
     NekDouble tol = 1e-3;
     NekDouble ytmp;
     for(int e=0; e<npoints; e++)
     {         
     	   coord[0] =xc[e];
       	   coord[1] =yc[e];
            	   
       	   elmtid = streak->GetExpIndex(coord,0.00001);
       	   offset = streak->GetPhys_Offset(elmtid);
       	   F = 1000;
       	   its=0;
           attempt=0;     
           ytmp = coord[1];
//cout<<"start guess:  x="<<xc[e]<<"    y="<<yc[e]<<endl;
           while( abs(F)> 0.000000001)
	   {

      	        elmtid = streak->GetExpIndex(coord,0.00001);

              	        //stvalues[e] = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() +offset );
//cout<<"elmtid="<<elmtid<<"  x="<<coord[0]<<"   y="<<coord[1]<<"    stvalue="<<U<<"   dU="<<dU<<endl;

                if( (abs(coord[1])>1 || elmtid==-1) 
                          && attempt==0  && verts==true
                )
                {                          
                     //try the old crit lay position:
                     coord[1] = yold_c[e];
                     attempt++;
                }
                else if( (abs(coord[1])>1 || elmtid==-1)  )
                {
                     coord[1] = ytmp +0.01;
                     elmtid = streak->GetExpIndex(coord,0.001);
       	             offset = streak->GetPhys_Offset(elmtid);
           	     NekDouble Utmp = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);
                     NekDouble dUtmp = streak->GetExp(elmtid)->PhysEvaluate(coord, derstreak + offset);
   	   	     coord[1] = coord[1] - (Utmp-cr)/dUtmp;
                     if( (abs(Utmp-cr)>abs(F))||(abs(coord[1])>1)  )
                     {
                          coord[1] = ytmp -0.01;
                     }
                             
                     attempt++;
                }
                else
                {
                     ASSERTL0(abs(coord[1])<= 1, " y value out of bound +/-1");
                }
      	        offset = streak->GetPhys_Offset(elmtid);
         	U = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);
        	dU  = streak->GetExp(elmtid)->PhysEvaluate(coord, derstreak + offset);
	   	coord[1] = coord[1] - (U-cr)/dU;   
	   	F = U-cr;   
	   	ASSERTL0( coord[0]==xc[e], " x coordinate must remain the same");                

                its++;
                if(its>200 && abs(F)<0.00001)
                {
                      	cout<<"warning streak position obtained with precision:"<<F<<endl;
                      	break;
                }               
                else if(its>1000 && abs(F)< 0.0001)
                {
                      	cout<<"warning streak position obtained with precision:"<<F<<endl;
                      	break;
                }
                else if(its>1000)
                {
                    ASSERTL0(false, "no convergence after 1000 iterations");
                }
          }
          yc[e] = coord[1] - (U-cr)/dU;
          ASSERTL0( U<= cr + tol, "streak wrong+");
          ASSERTL0( U>= cr -tol, "streak wrong-");
          //Utilities::Zerofunction(coord[0], coord[1], xtest, ytest, streak, derstreak);
cout<<"result streakvert x="<<xc[e]<<"  y="<<yc[e]<<"   streak="<<U<<endl;	           
          }
          
}


void GenerateAddPointsNewtonIt( NekDouble xi, NekDouble yi,NekDouble &xout, NekDouble &yout,
    	        MultiRegions::ExpListSharedPtr function, Array<OneD, NekDouble> derfunction,
                NekDouble cr)
{
	int elmtid,offset;
	NekDouble F,U,dU;
	Array<OneD, NekDouble> coords(2);	

        coords[0] = xi;
        coords[1] = yi;
	F =1000;
        int attempt=0;
        int its=0;
        NekDouble ytmp;
        ytmp = coords[1];
        while( abs(F)> 0.00000001)
        {

//cout<<"generate newton it xi="<<xi<<"  yi="<<yi<<endl;			
	     elmtid = function->GetExpIndex(coords, 0.01);
                //@to do if GetType(elmtid)==triangular WRONG!!!
cout<<"gen newton xi="<<xi<<"  yi="<<coords[1]<<"  elmtid="<<elmtid<<"  F="<<F<<endl;			
  
             if(  (abs(coords[1])>1 || elmtid==-1)     )
             {

                  coords[1] = ytmp +0.01;
                  elmtid = function->GetExpIndex(coords,0.01);
      	          offset = function->GetPhys_Offset(elmtid);
                  NekDouble Utmp = function->GetExp(elmtid)->PhysEvaluate(coords, function->GetPhys() + offset);
                  NekDouble dUtmp = function->GetExp(elmtid)->PhysEvaluate(coords, derfunction + offset);
	   	  coords[1] = coords[1] - (Utmp-cr)/dUtmp;
cout<<"attempt:"<<coords[1]<<endl;
                  if( (abs(Utmp-cr)>abs(F))||(abs(coords[1])>1.01)  )
                  {
                        coords[1] = ytmp -0.01;
                  }
                             
                  attempt++;
             }
             else if( abs(coords[1])<1.01 &&attempt==0)
             {
                  coords[1] =1.0;
                  attempt++;
             }
             else
             {
                  ASSERTL0(abs(coords[1])<= 1.00, " y value out of bound +/-1");
             }
	     offset = function->GetPhys_Offset(elmtid);
	     U = function->GetExp(elmtid)->PhysEvaluate(coords, function->GetPhys() + offset);
	     dU  = function->GetExp(elmtid)->PhysEvaluate(coords, derfunction + offset);
	     coords[1] = coords[1] - (U-cr)/dU;   
cout<<cr<<"U-cr="<<U-cr<<"  tmp result y:"<<coords[1]<<"  dU="<<dU<<endl;
	     F = U-cr;             
	     
             its++;
             if(its>200 && abs(F)<0.00001)
             {
                   cout<<"warning streak position obtained with precision:"<<F<<endl;
                   break;
             }               
             else if(its>1000 && abs(F)< 0.0001)
             {
                   cout<<"warning streak position obtained with precision:"<<F<<endl;
                   break;
             }
             else if(its>1000)
             {
                  ASSERTL0(false, "no convergence after 1000 iterations");
             }


             ASSERTL0( coords[0]==xi, " x coordinate must remain the same");	             	
        }
        xout = xi;
        yout = coords[1] - (U-cr)/dU;
cout<<"NewtonIt result  x="<<xout<<"  y="<<coords[1]<<"   U="<<U<<endl;	        
}

void GenerateMapEidsv1v2(MultiRegions::ExpListSharedPtr field,
                 Array<OneD, int> &V1, Array<OneD, int> &V2)
{

      const boost::shared_ptr<StdRegions::StdExpansionVector> exp2D = field->GetExp();
      int nel        = exp2D->size();
      LocalRegions::QuadExpSharedPtr locQuadExp;
      LocalRegions::TriExpSharedPtr  locTriExp;
      SpatialDomains::Geometry1DSharedPtr SegGeom;
      int id;
      int cnt=0;
      Array<OneD, int> V1tmp(4*nel, 10000);
      Array<OneD, int> V2tmp(4*nel, 10000);
      for(int i=0; i<nel; i++)
      { 
           if((locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>((*exp2D)[i])))
           {
                for(int j = 0; j < locQuadExp->GetNedges(); ++j)
                {
                    SegGeom = (locQuadExp->GetGeom2D())->GetEdge(j);
                    id = SegGeom->GetEid();
                    if( V1tmp[id] == 10000)
                    {
                         V1tmp[id]= SegGeom->GetVertex(0)->GetVid();
                         V2tmp[id]= SegGeom->GetVertex(1)->GetVid();
                         cnt++;
                    }
                }         
           }
           //in the future the tri edges may be not necessary (if the nedges is known)

           else if((locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>((*exp2D)[i])))
           {
                for(int j = 0; j < locTriExp->GetNedges(); ++j)
                {
                     SegGeom = (locTriExp->GetGeom2D())->GetEdge(j);
                     id = SegGeom->GetEid();

                     if( V1tmp[id] == 10000)
                     {
                          V1tmp[id]= SegGeom->GetVertex(0)->GetVid();
                          V2tmp[id]= SegGeom->GetVertex(1)->GetVid();
                          cnt++;
                     }
                 }

           }
           
      }  

      V1 = Array<OneD, int>(cnt);
      V2 = Array<OneD, int>(cnt);
      //populate the output vectors V1,V2
      for(int g=0; g<cnt; g++)
      {
           V1[g] = V1tmp[g];
           ASSERTL0(V1tmp[g]!=10000, "V1 wrong");
           V2[g] = V2tmp[g];
           ASSERTL0(V2tmp[g]!=10000, "V2 wrong");
      }




}





void MappingEVids(Array<OneD, NekDouble> xoldup, Array<OneD, NekDouble> yoldup,
                 Array<OneD, NekDouble> xolddown, Array<OneD, NekDouble> yolddown,
                 Array<OneD, NekDouble> xcold, Array<OneD, NekDouble> ycold,
                 Array<OneD, int> Vids_c,
                 SpatialDomains::MeshGraphSharedPtr mesh,
                 MultiRegions::ExpListSharedPtr streak,
                 Array<OneD, int> V1, Array<OneD, int> V2,
                 int & nlays,  Array<OneD, Array<OneD, int> >& Eids_lay, 
                 Array<OneD, Array<OneD, int> >& Vids_lay)
{

      int nlay_Eids = xcold.num_elements()-1;
      int nlay_Vids = xcold.num_elements();
     
      int nVertsTot = mesh->GetNvertices();
cout<<"nverttot="<<nVertsTot<<endl;
      //find the first vert of each layer 
      //int qp_closerup,qp_closerdown;
      //NekDouble diffup,diffdown;
cout<<"init nlays="<<nlays<<endl;
      //tmp first vert array (of course <100!!)
      Array<OneD, int> tmpVids0(100);
      Array<OneD, NekDouble> tmpx0(100);
      Array<OneD, NekDouble> tmpy0(100);
      Array<OneD, NekDouble> x2D(nVertsTot);  
      Array<OneD, NekDouble> y2D(nVertsTot);  
cout<<"yoldup="<<yoldup[0]<<endl;
cout<<"yolddown="<<yolddown[0]<<endl;

      for(int r=0; r< nVertsTot; r++)
      {

           SpatialDomains::VertexComponentSharedPtr vertex = mesh->GetVertex(r);
           NekDouble x,y,z;
           vertex->GetCoords(x,y,z); 
           x2D[r] = x;
           y2D[r] = y;
           if( xcold[0]==x  )
           {
                //check if the vert is inside layer zone
                if(
                   y<= yoldup[0] && y>= yolddown[0] 
                   && y!= ycold[0]
                  )
                {
//cout<<"x="<<x<<"  y="<<y<<endl;
                    tmpVids0[nlays] =r;  
                    tmpx0[nlays] = x;
                    tmpy0[nlays] = y;                 
                    nlays++;
                }
               
           }
      }
cout<<"nlays="<<nlays<<endl;

      //reorder first verts from xlower to xhigher
      Array<OneD, NekDouble> tmpx(nlays);
      Array<OneD, NekDouble> tmpf(nlays);
      Array<OneD, int> tmpV(nlays);
      //local copy to prevent overwriting
      Vmath::Vcopy(nlays, tmpx0,1,tmpx,1);
      Vmath::Vcopy(nlays, tmpy0,1,tmpf,1);     
      Vmath::Vcopy(nlays, tmpVids0,1,tmpV,1); 
      NekDouble max = Vmath::Vmax(tmpf.num_elements(), tmpf,1);
      int index=0;
      for(int w=0; w< nlays; w++)
      {
           index = Vmath::Imin(tmpf.num_elements(), tmpf,1);
           tmpx0[w]= tmpx[index];
           tmpy0[w]= tmpf[index];
           tmpVids0[w] = tmpV[index];
           tmpf[index] = max+1000;
                           
      }





      //initialise layers arrays
      Eids_lay = Array<OneD, Array<OneD,int> >(nlays);
      Vids_lay = Array<OneD, Array<OneD,int> >(nlays);
      for(int m=0; m<nlays; m++)
      {
           Eids_lay[m] = Array<OneD, int> (nlay_Eids);
           Vids_lay[m] = Array<OneD, int> (nlay_Vids);
      }
      
      //determine the others verts and edge for each layer
      NekDouble normbef, normtmp,xbef,ybef,xtmp,ytmp,normnext,xnext,ynext,diff;
      NekDouble Ubef, Utmp, Unext,diffU;
      Array<OneD, NekDouble> coord(2);
      int elmtid,offset;     
      int nTotEdges = V1.num_elements();
      Array<OneD, int> edgestmp(6);
      for(int m=0; m<nlays; m++)
      {
          for(int g=0; g<nlay_Eids; g++)
          {
             if(g==0)
             {
             for(int h=0; h< nTotEdges; h++)
             {

                 if( tmpVids0[m]== V1[h] )
                 {
                      SpatialDomains::VertexComponentSharedPtr vertex = mesh->GetVertex(V2[h]);
                      NekDouble x,y,z;
                      vertex->GetCoords(x,y,z); 
                      if(x!=tmpx0[m])
                      {
                            Eids_lay[m][0] = h;
                            Vids_lay[m][0] = V1[h];
                            Vids_lay[m][1] = V2[h];
                            SpatialDomains::VertexComponentSharedPtr vertex1 
                                         = mesh->GetVertex(V1[h]);
                            NekDouble x1,y1,z1;
                            vertex1->GetCoords(x1,y1,z1);
                            normbef= sqrt( (y-y1)*(y-y1)+(x-x1)*(x-x1)  );
                            ybef = (y-y1);
                            xbef = (x-x1);
                            coord[0] = x;
                            coord[1] = y;
                            elmtid = streak->GetExpIndex(coord,0.00001);
                            offset = streak->GetPhys_Offset(elmtid);
                            Ubef = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);

                      }
 
                 }
                 if( tmpVids0[m]== V2[h] )
                 {                    
                      SpatialDomains::VertexComponentSharedPtr vertex = mesh->GetVertex(V1[h]);
                      NekDouble x,y,z;
                      vertex->GetCoords(x,y,z); 
                      if(x!=tmpx0[m])
                      {
                            Eids_lay[m][0] = h;
                            Vids_lay[m][0] = V2[h];
                            Vids_lay[m][1] = V1[h];                         
                            SpatialDomains::VertexComponentSharedPtr vertex2 
                                         = mesh->GetVertex(V2[h]);
                            NekDouble x2=0.0,y2=0.0,z2=0.0;
                            normbef= sqrt( (y-y2)*(y-y2)+(x-x2)*(x-x2)  );
                            ybef = (y-y2);
                            xbef = (x-x2);
                            coord[0] = x;
                            coord[1] = y;
                            elmtid = streak->GetExpIndex(coord,0.00001);
                            offset = streak->GetPhys_Offset(elmtid);
                            Ubef = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);                        

                      }
                 }

             }

cout<<"Eid="<<Eids_lay[m][0]<<"   Vids_lay0="<<Vids_lay[m][0]<<"   Vidslay1="<<Vids_lay[m][1]<<endl;

             }
             else
             {
             //three verts(edges) candidates
             int cnt =0;
             for(int h=0; h< nTotEdges; h++)
             {
//cout<<Vids_lay[m][g]<<"    V1="<<V1[h]<<"     V2[h]="<<V2[h]<<endl;
                  if( (Vids_lay[m][g]==V1[h] || Vids_lay[m][g]==V2[h]) && h!= Eids_lay[m][g-1])
                  {
cout<<"edgetmp="<<h<<endl;
                       ASSERTL0(cnt<=6, "wrong number of candidates"); 
                       edgestmp[cnt]= h;
                       cnt++;
                  }
             }

             diff =1000;
             diffU =1000;
             Array<OneD, NekDouble > diffarray(cnt);             
             Array<OneD, NekDouble > diffUarray(cnt);
cout<<"normbef="<<normbef<<endl;
cout<<"Ubefcc="<<Ubef<<endl;
             //choose the right candidate
             for(int e=0; e< cnt; e++)
             {
                 SpatialDomains::VertexComponentSharedPtr vertex1 = mesh->GetVertex(V1[edgestmp[e]]);
                 NekDouble x1,y1,z1;
                 vertex1->GetCoords(x1,y1,z1); 
                 SpatialDomains::VertexComponentSharedPtr vertex2 = mesh->GetVertex(V2[edgestmp[e]]);
                 NekDouble x2,y2,z2;
                 vertex2->GetCoords(x2,y2,z2); 
                 
                 normtmp= sqrt( (y2-y1)*(y2-y1)+(x2-x1)*(x2-x1)  );  

cout<<"edgetmp1="<<edgestmp[e]<<endl;
cout<<"V1 x1="<<x1<<"  y1="<<y1<<endl;
cout<<"V2 x2="<<x2<<"  y2="<<y2<<endl;
                 if( Vids_lay[m][g]==V1[edgestmp[e]] )
                 {


                      ytmp = (y2-y1);
                      xtmp = (x2-x1);
                      coord[0] = x2;
                      coord[1] = y2;
                      elmtid = streak->GetExpIndex(coord,0.00001);
                      offset = streak->GetPhys_Offset(elmtid);
                      Utmp = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);  
                      diffarray[e] = abs((xtmp*xbef+ytmp*ybef)/(normtmp*normbef)-1);
                      diffUarray[e] = abs(Ubef-Utmp);
cout<<"   normtmp="<<normtmp<<endl;
cout<<"   Utmpcc="<<Utmp<<endl;
cout<<xtmp<<"  ytmp="<<ytmp<<"    diff="<<abs(((xtmp*xbef+ytmp*ybef)/(normtmp*normbef))-1)<<endl;
                      if( 
                      	  abs(  (xtmp*xbef+ytmp*ybef)/(normtmp*normbef)-1)<diff  
                      	  && y2<= yoldup[g+1]  &&   y2>= yolddown[g+1]
                      	  && y1<= yoldup[g]  &&   y1>= yolddown[g]
                      	)
                      {

                          Eids_lay[m][g] = edgestmp[e];
                          Vids_lay[m][g+1] = V2[edgestmp[e]];
                          diff = abs((xtmp*xbef+ytmp*ybef)/(normtmp*normbef)-1);
                          normnext =normtmp;
                          ynext = ytmp;
                          xnext = xtmp; 
                          Unext = Utmp;
                          diffU = abs(Ubef-Utmp);
                      }
                 }
                 else if( Vids_lay[m][g]==V2[edgestmp[e]] )
                 {


                      ytmp = (y1-y2);
                      xtmp = (x1-x2);
                      coord[0] = x1;
                      coord[1] = y1;
                      elmtid = streak->GetExpIndex(coord,0.00001);
                      offset = streak->GetPhys_Offset(elmtid);
                      Utmp = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);  
                      diffarray[e] = abs((xtmp*xbef+ytmp*ybef)/(normtmp*normbef)-1);                      
                      diffUarray[e] = abs(Ubef-Utmp);                      
cout<<"   normtmp="<<normtmp<<endl;
cout<<"   Utmpcc="<<Utmp<<endl;
cout<<xtmp<<"  ytmp="<<ytmp<<"    diff="<<abs(((xtmp*xbef+ytmp*ybef)/(normtmp*normbef))-1)<<endl;
                      if( 
                      	  abs((xtmp*xbef+ytmp*ybef)/(normtmp*normbef)-1)<diff 
                      	  && y2<= yoldup[g]  &&   y2>= yolddown[g]
                      	  && y1<= yoldup[g+1]  &&   y1>= yolddown[g+1]                      	      
                      	)
                      {
                          Eids_lay[m][g] = edgestmp[e];
                          Vids_lay[m][g+1] = V1[edgestmp[e]];
                          diff = abs((xtmp*xbef+ytmp*ybef)/(normtmp*normbef)-1);
                          normnext =normtmp;
                          ynext = ytmp;
                          xnext = xtmp; 
                          Unext = Utmp;
                          diffU = abs(Ubef-Utmp);                                             
                      }
                      	     
                 }                 
                 else
                 {
                      ASSERTL0(false, "eid not found");
                 }
                 


             }
cout<<"Eid before check="<<Eids_lay[m][g]<<endl;
for(int q=0; q<cnt; q++)
{
cout<<q<<"   diff"<<diffarray[q]<<endl;	
}
             //check if the eid has a vert in common with another layer
             bool check =false;
             if(m>0 && m< nlays)
             {	     
                  check = checkcommonvert(Vids_lay[m-1],Vids_c,Vids_lay[m][g+1]);                  
             }
             if(check == true)
             {
cout<<"COMMON VERT"<<endl;             	     
             	  int eid = Vmath::Imin(cnt, diffarray,1);
             	  diffarray[eid]=1000;
             	  eid = Vmath::Imin(cnt,diffarray,1);
             	  
             	  
                  SpatialDomains::VertexComponentSharedPtr vertex1 = mesh->GetVertex(V1[edgestmp[eid]]);
                  NekDouble x1,y1,z1;
                  vertex1->GetCoords(x1,y1,z1); 
                  SpatialDomains::VertexComponentSharedPtr vertex2 = mesh->GetVertex(V2[edgestmp[eid]]);
                  NekDouble x2,y2,z2;
                  vertex2->GetCoords(x2,y2,z2);                    
                   
                  normtmp= sqrt( (y2-y1)*(y2-y1)+(x2-x1)*(x2-x1)  );  
                  
                  Eids_lay[m][g] = edgestmp[eid]; 
                  if(Vids_lay[m][g] == V1[edgestmp[eid]])
                  {
                       coord[0] = x2;
                       coord[1] = y2;
                       elmtid = streak->GetExpIndex(coord,0.00001);
                       offset = streak->GetPhys_Offset(elmtid);
                       Utmp = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);                    	  
                       Vids_lay[m][g+1] = V2[edgestmp[eid]];
                       normnext =normtmp;
                       ynext = (y2-y1);
                       xnext = (x2-x1);; 
                       Unext = Utmp;
                              
                              
                  }
                  if(Vids_lay[m][g] == V2[edgestmp[eid]])
                  {
                       coord[0] = x1;
                       coord[1] = y1;
                       elmtid = streak->GetExpIndex(coord,0.00001);
                       offset = streak->GetPhys_Offset(elmtid);
                       Utmp = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);                    	  
                       Vids_lay[m][g+1] = V1[edgestmp[eid]];
                       normnext =normtmp;
                       ynext = (y1-y2);
                       xnext = (x1-x2);; 
                       Unext = Utmp;
                              
                              
                  }                   	   
             	     
             }          
             
cout<<m<<"edge aft:"<<Eids_lay[m][g]<<"   Vid="<<Vids_lay[m][g+1]<<endl;
             normbef = normnext;
             ybef = ynext;
             xbef = xnext;
             Ubef = Unext;
             
cout<<"endelse"<<normtmp<<endl;
             } //end else

         }//end g it
            

          



      } //end m it     

for(int w=0; w< nlays; w++)
{
for(int f=0; f< nlay_Eids; f++)
{
cout<<"check="<<w<<"   Eid:"<<Eids_lay[w][f]<<endl;
}
}

}


bool checkcommonvert(Array<OneD, int> Vids_laybefore, Array<OneD, int> Vids_c, int Vid)
{
     bool check=false;
     for(int u=0; u< Vids_laybefore.num_elements(); u++)
     {
           if( Vids_laybefore[u]==Vid || Vids_c[u]==Vid)
           {
                check =true;
           }   
cout<<Vid<<"    Vert test="<<Vids_laybefore[u]<<endl;           
     }
     return check;
	
}


void  Cutrepetitions(int nedges,Array<OneD, NekDouble> inarray,
                 Array<OneD, NekDouble>& outarray)
{
        	
      //determine npedge:
      int np_lay = inarray.num_elements();
           
      int npedge = np_lay/nedges;
      ASSERTL0(inarray.num_elements()%nedges==0," something on number npedge");
      //cut out the repetitions:

      int cnt=0;
      for(int w=0; w< np_lay; w++)
      {

           if(w< np_lay-1)
           {
                if(inarray[w] ==inarray[w+1])
                {
                     continue;
                }
           }                     
           outarray[cnt]= inarray[w];
           cnt++;
      }


      ASSERTL0( cnt== np_lay-(nedges-1), "wrong cut");
        	
}

int DetermineclosePointxindex(NekDouble x,Array<OneD, NekDouble> xArray)
{
      int npts = xArray.num_elements();
      Array<OneD, NekDouble> xcopy(npts);
      Vmath::Vcopy(npts,xArray,1,xcopy,1);
      //subtract xpoint and abs
            
      Vmath::Sadd(npts, -x, xcopy,1, xcopy,1);
      Vmath::Vabs(npts, xcopy,1,xcopy,1);

      int index = Vmath::Imin(npts, xcopy,1);
      if(xArray[index]> x)//assume always x[index]< x(HYPHOTHESIS)
      {
            index = index-1;
      }
      return index;
}

void GenerateNeighbourArrays(int index,int neighpoints, Array<OneD, NekDouble> xArray,
                 Array<OneD, NekDouble> yArray,Array<OneD, NekDouble>& Neighbour_x,
                 Array<OneD, NekDouble>& Neighbour_y)
{
      ASSERTL0( neighpoints%2==0,"number of neighbour points should be even");
      int leftpoints = (neighpoints/2)-1;//-1 because xArray[index]< x
      int rightpoints = neighpoints/2;
      int diff ;
      int start;
//cout<<"index="<<index<<"  left="<<leftpoints<<"  right="<<rightpoints<<endl;
      if(index-leftpoints<0)
      {
//cout<"case0"<<endl;
            diff = index-leftpoints; 
            start= 0;
            Vmath::Vcopy(neighpoints, &yArray[0],1,&Neighbour_y[0],1);
            Vmath::Vcopy(neighpoints, &xArray[0],1,&Neighbour_x[0],1);
      }
      else if( (yArray.num_elements()-1)-index < rightpoints)
      {
//cout"case1 closest="<<xArray[index]<<endl;
            int rpoints = (yArray.num_elements()-1)-index;//
            diff = rightpoints-rpoints;
//cout<"start index="<<index-leftpoints-diff<<endl;
            start = index-leftpoints-diff;
            Vmath::Vcopy(neighpoints, &yArray[start],1,&Neighbour_y[0],1);
            Vmath::Vcopy(neighpoints, &xArray[start],1,&Neighbour_x[0],1);
      }
      else
      {
//cout<<"caseaa"<<endl;
            start = index-leftpoints;
            Vmath::Vcopy(neighpoints, &yArray[start],1,&Neighbour_y[0],1);
            Vmath::Vcopy(neighpoints, &xArray[start],1,&Neighbour_x[0],1);
      }
/*
for(int t= start; t<start+neighpoints; t++)
{
cout<<"Px="<<xArray[t]<<"    "<<yArray[t]<<endl;
}
*/
      //check if there is any repetition
      for(int f=1; f< neighpoints; f++)
      {
            ASSERTL0(Neighbour_x[f]!=Neighbour_x[f-1]," repetition on NeighbourArrays");
      }
                      
}

NekDouble LagrangeInterpolant(NekDouble x, int npts, 
                 Array<OneD,NekDouble>  xpts, Array<OneD, NekDouble> funcvals)
{
      NekDouble sum = 0.0;
      NekDouble LagrangePoly;
//cout<<"lagrange"<<endl;
      for(int pt=0;pt<npts;++pt)
      {
           NekDouble h=1.0;
            
           for(int j=0;j<pt; ++j)
           { 
                h = h * (x - xpts[j])/(xpts[pt]-xpts[j]);
           }

           for(int k=pt+1;k<npts;++k)
           {
                h = h * (x - xpts[k])/(xpts[pt]-xpts[k]);
           }  
           LagrangePoly=h;              

           sum += funcvals[pt]*LagrangePoly;
      }
//cout<<"result :"<<sum<<endl;
      return sum;
}


void EvaluateTangent(int npoints, Array<OneD, NekDouble> xcQedge, 
                 Array<OneD, NekDouble> coeffsinterp,
                 Array<OneD, NekDouble> & txQedge, Array<OneD, NekDouble> & tyQedge)
{
      Array<OneD, NekDouble>  yprime(npoints,0.0);
      int np_pol= coeffsinterp.num_elements();
cout<<"evaluatetan with "<<np_pol<<endl;

      //calc derivative poly (cut last entry on b array that is the const)
      int derorder;
      Array<OneD, NekDouble> yinterp(npoints,0.0);
      int polorder;
      for(int q=0; q< npoints; q++)
      {
            polorder=np_pol-1;
            derorder=np_pol-2;
            yprime[q] =0;
            for(int d=0; d< np_pol-1; d++)
            {
                 yprime[q] += (derorder +1)*coeffsinterp[d]*std::pow(xcQedge[q],derorder);
                 derorder--;                        
            }
            //test
            for(int a=0; a< np_pol; a++)
            {
                 yinterp[q] += coeffsinterp[a]*std::pow(xcQedge[q],polorder);
//cout<<"coeff*x^n="<<b[a]*std::pow(xcQedge[q],polorder)<<" sum="<<yinterp[q]<<endl;
                 polorder--;                        
            }

      }    
              
      //transf yprime into tx,ty:
      for(int n=0; n< npoints; n++)
      {
            //ABS???!!
            txQedge[n]=0;
            txQedge[n] = cos((atan((yprime[n]))));
            tyQedge[n] = sin((atan((yprime[n]))));
cout<<xcQedge[n]<<"      "<<yinterp[n]<<"      "<<yprime[n]<<"      "<<txQedge[n]<<"       "<<tyQedge[n]<<endl;
      }
}

void PolyInterp( Array<OneD, NekDouble> xpol, Array<OneD, NekDouble> ypol,
                 Array<OneD, NekDouble> & coeffsinterp,
                 Array<OneD, NekDouble> & xcout, Array<OneD, NekDouble> & ycout,  
                 int edge, int npedge)
{
      int np_pol = xpol.num_elements();
      int N = np_pol;    
      Array<OneD, NekDouble> A (N*N,1.0);
      Array<OneD, NekDouble> b (N);
      int row=0;
      //fill column by column
      for(int e=0; e<N; e++)
      {
           row=0;
           for(int w=0; w < N; w++)
           {         
                 A[N*e+row] = std::pow( xpol[w], N-1-e);
                 row++;
           }	    
       }
       row=0;
       for(int r= 0; r< np_pol; r++)
       {
    	   b[row] =   ypol[r];
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
/*
for(int h=0; h<np_pol; h++)
{
cout<<"coeff:"<<b[h]<<endl;
} 
*/     



        //ovewrite the ycPhysMOD:
        int polorder;
        for(int c=0; c< npedge; c++)
        {
              polorder=np_pol-1;
              //put ycPhysMOD to zero
              ycout[edge*(npedge)+c+1]=0;
              for(int d=0; d< np_pol; d++)
              {
                   ycout[edge*(npedge)+c+1] += b[d]
                            *std::pow(xcout[edge*(npedge)+c+1],polorder);
                   polorder--;                        
              } 
        }

        //write coeffs
        Vmath::Vcopy(np_pol, b,1,coeffsinterp,1);

}


void PolyFit(int polyorder,int npoints,
                 Array<OneD, NekDouble> xin, Array<OneD, NekDouble>  fin,
                 Array<OneD, NekDouble> & coeffsinterp,
                 Array<OneD, NekDouble> & xout, Array<OneD, NekDouble> & fout,  
                 int npout)
{

        int N = polyorder+1;
        Array<OneD, NekDouble> A (N*N,0.0);
        Array<OneD, NekDouble> b (N,0.0);
cout<<npoints<<endl;
for(int u=0; u<npoints; u++)
{
cout<<"c="<<xin[u]<<"     "<<
fin[u]<<endl;
}
        //fill column by column
        //e counts cols
        for(int e=0; e<N; e++)
        {
    	          
             for(int row=0; row<N; row++)
             {
       	    	   for(int w=0; w < npoints; w++)
    	           {         
                        A[N*e+row] += std::pow( xin[w], e+row);
                   }	    
             }
        }
              
        for(int row= 0; row< N; row++)
        {
             for(int h=0; h< npoints; h++)
             {
    	         b[row] +=   fin[h]*std::pow(xin[h],row);
//cout<<b="<<b[row]<<"  y_c="<<y_c[r]<<endl;     	        
             }
    	                 
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

       //the lower coeff the lower is the grad
       //reverse:
       Array<OneD, NekDouble> tmp(N);
       Vmath::Vcopy(N,b,1,tmp,1);
       int cnt =N;
       for(int j=0; j<N; j++)
       {
              b[j]= tmp[cnt-1];
              cnt--;
       }

for(int h=0; h<N; h++)
{
cout<<"coeff:"<<b[h]<<endl;
} 

       //ovewrite the ycPhysMOD:
       int polorder;
       for(int c=0; c< npout; c++)
       {
             polorder=polyorder;
             //put ycPhysMOD to zero
             fout[c]=0;
             for(int d=0; d< N; d++)
             {
                        
                  fout[c] += b[d]
                                *std::pow(xout[c],polorder);
                   polorder--;                        

             }

       }

       //write coeffs
       Vmath::Vcopy(N, b,1,coeffsinterp,1);
             
}


void  Orderfunctionx(Array<OneD, NekDouble> inarray_x,
                 Array<OneD, NekDouble> inarray_y, Array<OneD, NekDouble>& outarray_x,
                 Array<OneD, NekDouble>& outarray_y)
{
       Array<OneD, NekDouble>tmpx(inarray_x.num_elements());
       Array<OneD, NekDouble>tmpy(inarray_x.num_elements());
       //local copy to prevent overwriting
       Vmath::Vcopy(inarray_x.num_elements() , inarray_x,1,tmpx,1);
       Vmath::Vcopy(inarray_x.num_elements() , inarray_y,1,tmpy,1);

       //order function with respect to x
       int index;
       NekDouble max = Vmath::Vmax(tmpx.num_elements(), tmpx,1);
       for(int w=0; w<tmpx.num_elements(); w++)
       {
            index = Vmath::Imin(tmpx.num_elements(), tmpx,1);
            outarray_x[w]= tmpx[index];
            outarray_y[w]= tmpy[index];
            if(w< tmpx.num_elements()-1)//case of repetitions
            {
                 if(tmpx[index] == tmpx[index+1])
                 {
                      outarray_x[w+1]= tmpx[index+1];
                      outarray_y[w+1]= tmpy[index+1];
                      tmpx[index+1] = max+1000;
                      w++;
                 }
            }
/*
                 if(w>0)//case of repetitions
                 {
                     if(inarray_x[index] == tmpx[index-1])
                     {
                          outarray_x[w+1]= tmpx[index-1];
                          outarray_y[w+1]= tmpy[index-1];
                          tmpx[index-1] = max+1000;
                          w++;
                     }
                 }
*/
            tmpx[index] = max+1000;
                             
      }
}          

void MoveLayersvertically(int nlays, int nvertl, int cntlow, int cntup,
    	         Array<OneD, Array<OneD, int > > lay_Vids,  Array<OneD, NekDouble> xc,
    	         Array<OneD, NekDouble> yc, Array<OneD, int> Down, Array<OneD, int> Up,
    	         Array<OneD, NekDouble >& xnew,Array<OneD, NekDouble>& ynew,
    	         Array<OneD, Array<OneD, NekDouble > >& layers_x,
    	         Array<OneD, Array<OneD, NekDouble > >& layers_y)    
{
       int np_lay = layers_y[0].num_elements();
       // 0<h<nlays-1 to fill only the 'internal' layers(no up,low);
       for(int h=1; h<nlays-1; h++)
       {
             layers_x[h]= Array<OneD, NekDouble>(np_lay);             	     
             for(int s=0; s<nvertl; s++)
             {
                   //check if ynew is still empty
                   ASSERTL0(ynew[ lay_Vids[h][s] ]==-20, "ynew layers not empty");
                   if(h<cntlow+1)//layers under the crit lay
                   {
                         //y= ylow+delta
                         ynew[ lay_Vids[h][s] ]  = ynew[Down[s]]+  h*abs(ynew[Down[s]] - yc[s])/(cntlow+1);
                         //put the layer vertical
                         xnew[lay_Vids[h][s]  ] = xc[s];
//cout<<"ynew="<<ynew[ lay_Vids[h][s] ]<<" ydown="<<ynew[Down[s]]<<
//" delta="<<abs(ynew[Down[s]] - y_c[s])/(cntlow+1)<<endl;
                  	 //until now layers_y=yold
                	 layers_y[h][s] = ynew[ lay_Vids[h][s] ];
                	 layers_x[h][s] = xnew[ lay_Vids[h][s] ];                	    
                   }
                   else
                   {
                         //y = yc+delta
                         ynew[ lay_Vids[h][s] ]  = yc[s] + (h-cntlow)*abs(ynew[Up[s]] - yc[s])/(cntup+1);
                         //put the layer vertical
                         xnew[lay_Vids[h][s]  ] = xc[s];
                         //until now layers_y=yold
                         layers_y[h][s] = ynew[ lay_Vids[h][s] ];
                	 layers_x[h][s] = xnew[ lay_Vids[h][s] ];                         	    
                   } 
             }
       }        	
        	
}


void MoveLayerNfixedxpos(int nvertl, int npedge, Array<OneD, NekDouble> xcPhys,
	         Array<OneD, NekDouble> tmpx_lay, Array<OneD, NekDouble> tmpy_lay,
	         Array<OneD, int> Vids,
	         Array<OneD, NekDouble> &xlay, Array<OneD, NekDouble> &ylay,
	         Array<OneD, NekDouble> &xnew, Array<OneD, NekDouble> &ynew)
{
       int np_lay  = xcPhys.num_elements();
       int nedges = nvertl-1;
       Array<OneD, NekDouble>tmpx(np_lay-(nedges-1));
       Array<OneD, NekDouble>tmpy(np_lay-(nedges-1));	
       Cutrepetitions(nedges, tmpx_lay, tmpx);
       Cutrepetitions(nedges, tmpy_lay, tmpy);   
       //order points in x:
       int index;
       int closepoints = 4;
       Array<OneD, NekDouble>Pxinterp(closepoints);
       Array<OneD, NekDouble>Pyinterp(closepoints);	       
       Orderfunctionx(tmpx, tmpy, tmpx, tmpy);         
       //determine the neighbour points (-3;+3)
       for(int g=0; g< nedges; g++)
       {
           //write vert coords
             //v1
           index= 
               DetermineclosePointxindex( xcPhys[g*npedge+0], tmpx);
           //generate neighbour arrays:
           GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
           ynew[Vids[g] ]= LagrangeInterpolant(xcPhys[g*npedge+0],closepoints,Pxinterp,Pyinterp  );
           xnew[Vids[g] ]= xcPhys[g*npedge+0]; 
           ylay[g*npedge +0] = ynew[ Vids[g] ];
           xlay[g*npedge +0] = xnew[ Vids[g] ];
           
           //v2             
           //determine closest index:             
           index= 
                 DetermineclosePointxindex( xcPhys[g*npedge +npedge-1], tmpx);
           //generate neighbour arrays:
           GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
           ynew[Vids[g+1] ]=  LagrangeInterpolant(xcPhys[g*npedge +npedge-1],closepoints,Pxinterp,Pyinterp  );
           xnew[Vids[g+1] ]=  xcPhys[g*npedge +npedge-1];
           ylay[g*npedge +npedge-1] = ynew[Vids[g+1] ];                       
           xlay[g*npedge +npedge-1] = xnew[Vids[g+1] ];



           //middle points
           for(int r=0; r< npedge-2; r++)
           {

               //determine closest point index:
               index = 
               DetermineclosePointxindex( xcPhys[g*npedge +r+1], tmpx);
//cout<<"  vert+"<<index<<endl;

               ASSERTL0( index<= tmpy.num_elements()-1, " index wrong");
               //generate neighbour arrays Pyinterp,Pxinterp
               GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
                 
               ylay[g*npedge +r+1]=
                                 LagrangeInterpolant(
                                 xcPhys[g*npedge +r+1],closepoints,Pxinterp,Pyinterp  );
//cout<<"x value="<<xcPhysMOD[g*npedge +r+1]<<endl;
               xlay[g*npedge +r+1]=  xcPhys[g*npedge +r+1];




/*
for(int t=0; t<6; t++)
{
cout<<"Px="<<Pxinterp[t]<<"    "<<Pyinterp[t]<<endl;
}
*/
            }
        
       }//edge closed g
	
}

void MoveLayerNnormpos(int nvertl, int npedge, Array<OneD, NekDouble> xcPhys,
	         Array<OneD, NekDouble> tmpx_lay, Array<OneD, NekDouble> tmpy_lay,
	         Array<OneD, int> Vids,
	         Array<OneD, NekDouble> &xlay, Array<OneD, NekDouble> &ylay,
	         Array<OneD, NekDouble> &xnew, Array<OneD, NekDouble> &ynew)
{
       int np_lay  = xcPhys.num_elements();
       int nedges = nvertl-1;
       NekDouble x0,x1, xtmp;
       Array<OneD, NekDouble>tmpx(np_lay-(nedges-1));
       Array<OneD, NekDouble>tmpy(np_lay-(nedges-1));	
       Cutrepetitions(nedges, tmpx_lay, tmpx);
       Cutrepetitions(nedges, tmpy_lay, tmpy);   
       //order points in x:
       int index;
       int closepoints = 4;
       Array<OneD, NekDouble>Pxinterp(closepoints);
       Array<OneD, NekDouble>Pyinterp(closepoints);	       
       Orderfunctionx(tmpx, tmpy, tmpx, tmpy);         
       //determine the neighbour points (-3;+3)
       for(int g=0; g< nedges; g++)
       {
           //write vert coords
           //v1
           ynew[Vids[g] ]= tmpy_lay[g*npedge+0]; 
           xnew[Vids[g] ]= tmpx_lay[g*npedge+0]; 


           ylay[g*npedge +0] = ynew[ Vids[g] ];
           xlay[g*npedge +0] = xnew[ Vids[g] ];
           
           //v2                       
           ynew[Vids[g+1] ]=  tmpy_lay[g*npedge+npedge-1]; 
           xnew[Vids[g+1] ]=  tmpx_lay[g*npedge+npedge-1]; 
           ylay[g*npedge +npedge-1] = ynew[Vids[g+1] ];                       
           xlay[g*npedge +npedge-1] = xnew[Vids[g+1] ];



           //middle points           
           for(int r=0; r< npedge-2; r++)
           {
               x0 = xlay[g*npedge +0];
               x1 = xlay[g*npedge +npedge-1];
               xtmp = x0 + r*(x1-x0)/(npedge-1);
               //determine closest point index:
               index = 
               DetermineclosePointxindex( xtmp, tmpx);
//cout<<"  vert+"<<index<<endl;

               ASSERTL0( index<= tmpy.num_elements()-1, " index wrong");
               //generate neighbour arrays Pyinterp,Pxinterp
               GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
                 
               ylay[g*npedge +r+1]=
                                 LagrangeInterpolant(
                                 xtmp,closepoints,Pxinterp,Pyinterp  );
//cout<<"x value="<<xcPhysMOD[g*npedge +r+1]<<endl;
               xlay[g*npedge +r+1]=  xtmp;


            }
        
       }//edge closed g	
}

void MoveOutsidePointsfixedxpos(int npedge, SpatialDomains::MeshGraphSharedPtr mesh,
	         Array<OneD, NekDouble> xcold,Array<OneD, NekDouble> ycold,
      	         Array<OneD, NekDouble> xolddown,Array<OneD, NekDouble> yolddown,
	         Array<OneD, NekDouble> xoldup,Array<OneD, NekDouble> yoldup,     
	         Array<OneD, NekDouble> ylaydown,Array<OneD, NekDouble> ylayup, 	         
                 Array<OneD, NekDouble>& xnew,Array<OneD, NekDouble>& ynew)
{
     //update vertices coords outside layers region
     NekDouble ratio;
     int nvertl = ycold.num_elements();
     int nVertTot =  mesh->GetNvertices();
     for(int n=0; n<nVertTot; n++)
     {
          NekDouble ratio;  
          SpatialDomains::VertexComponentSharedPtr vertex = mesh->GetVertex(n);
          NekDouble x,y,z;
          vertex->GetCoords(x,y,z); 
          int qp_closer,diff;
          //determine the closer xold_up
          NekDouble tmp=1000;

          for(int k=0; k<nvertl; k++)
          {     
              if(abs(x-xcold[k]) < tmp)
              {
                  tmp = abs(x-xcold[k]);
                  qp_closer=k;
              }         
     
          } 
          //find nplay_closer
          int nplay_closer;
          if(qp_closer==0)
          {
              nplay_closer=0;//first vert
          }
          else
          {
              nplay_closer= (qp_closer-1)*npedge +npedge-1;
          }


          if(  y>yoldup[qp_closer] && y<1 )//nlays-1 is layerup
          {	        

//              ratio = (1-layers_y[nlays-1][qp_closer])*(1-y_c[qp_closer])/
//                    (  (1-yold_up[n])*(1-yold_c[qp_closer]) );
              ratio = (1-ylayup[nplay_closer])/
                    (  (1-yoldup[qp_closer]) );
              //distance prop to layerup
              ynew[n] = ylayup[nplay_closer] 
                      + (y-yoldup[qp_closer])*ratio;  
              xnew[n] = x;
            
          }
          else if(   y< yolddown[qp_closer]   && y>-1  )//0 is layerdown
          {

              ratio = (1+ylaydown[nplay_closer])/
                    (  (1+yolddown[qp_closer]) );
              //distance prop to layerlow
              ynew[n] = ylaydown[nplay_closer] 
                      + (y-yolddown[qp_closer])*ratio;  
              xnew[n] = x;          
          }

     }     	
}

void MoveOutsidePointsNnormpos(int npedge, SpatialDomains::MeshGraphSharedPtr mesh,
	         Array<OneD, NekDouble> xcold,Array<OneD, NekDouble> ycold,
      	         Array<OneD, NekDouble> xolddown,Array<OneD, NekDouble> yolddown,
	         Array<OneD, NekDouble> xoldup,Array<OneD, NekDouble> yoldup,     
	         Array<OneD, NekDouble> xlaydown,Array<OneD, NekDouble> ylaydown,
	         Array<OneD, NekDouble> xlayup,Array<OneD, NekDouble> ylayup, 
	         Array<OneD, NekDouble> nxPhys,Array<OneD, NekDouble> nyPhys,	         
                 Array<OneD, NekDouble>& xnew,Array<OneD, NekDouble>& ynew) 
{	
/*	
     int nq1D =bndfieldup->GetTotPoints();
     Array<OneD, NekDouble> xlayoldup(nq1D);
     Array<OneD, NekDouble> xlayolddown(nq1D);	
     Array<OneD, NekDouble> ylayoldup(nq1D);
     Array<OneD, NekDouble> ylayolddown(nq1D);
     Array<OneD, NekDouble> zlayoldup(nq1D);
     Array<OneD, NekDouble> zlayolddown(nq1D);
     bndfielddown->GetCoords( xlayolddown,  ylayolddown,zlayolddown);     
     bndfieldup->GetCoords( xlayoldup,  ylayoldup,zlayoldup);
     
     NekDouble xmax = Vmath::Vmax(nq1D, xlayoldup,1);
     NekDouble xmin = Vmath::Vmin(nq1D, xlayoldup,1);     
*/     
     //determine the new verts up/down pos:
     int nvertl = xoldup.num_elements();
     int nedges = nvertl-1;
     Array<OneD, NekDouble> xnew_down(nvertl);
     Array<OneD, NekDouble> ynew_down(nvertl);
     Array<OneD, NekDouble> xnew_up(nvertl);
     Array<OneD, NekDouble> ynew_up(nvertl);
     Array<OneD, NekDouble> nxvert(nvertl);
     Array<OneD, NekDouble> nyvert(nvertl);
     Array<OneD, NekDouble> norm(nvertl);
     Array<OneD, NekDouble> tmp(nvertl);     
     for(int a=0; a< nedges;a++)
     {
     	  if(a==0)
     	  {
     	  //v1
          xnew_down[a] = xlaydown[a*npedge+0];
          ynew_down[a] = ylaydown[a*npedge+0];
          xnew_up[a] = xlayup[a*npedge+0];
          ynew_up[a] = ylayup[a*npedge+0];
          nxvert[a] = nxPhys[a*npedge+0];
          nyvert[a] = nyPhys[a*npedge+0];          
          //v2
          xnew_down[a+1] = xlaydown[a*npedge+npedge-1];
          ynew_down[a+1] = ylaydown[a*npedge+npedge-1];
          xnew_up[a+1] = xlayup[a*npedge+npedge-1];
          ynew_up[a+1] = ylayup[a*npedge+npedge-1];
          nxvert[a+1] = nxPhys[a*npedge+npedge-1];
          nyvert[a+1] = nyPhys[a*npedge+npedge-1];          
          }
          else
          {
          //v2
          xnew_down[a+1] = xlaydown[a*npedge+npedge-1];
          ynew_down[a+1] = ylaydown[a*npedge+npedge-1];
          xnew_up[a+1] = xlayup[a*npedge+npedge-1];
          ynew_up[a+1] = ylayup[a*npedge+npedge-1];   
          nxvert[a+1] = nxPhys[a*npedge+npedge-1];
          nyvert[a+1] = nyPhys[a*npedge+npedge-1];           
          }
          
     }
     
     NekDouble xmax = Vmath::Vmax(nvertl, xoldup,1);
     NekDouble xmin = Vmath::Vmin(nvertl, xoldup,1);       
     
     //update vertices coords outside layers region
     NekDouble ratio, ratiox,ratioy;

     int nVertTot =  mesh->GetNvertices();
     for(int n=0; n<nVertTot; n++)
     {
          NekDouble ratio;  
          SpatialDomains::VertexComponentSharedPtr vertex = mesh->GetVertex(n);
          NekDouble x,y,z;
          vertex->GetCoords(x,y,z); 
          int qp_closeroldup, qp_closerolddown;
          NekDouble diffup, diffdown;
          //determine the closer xold_up,down
          diffdown =1000;
          diffup = 1000;



          
          for(int k=0; k<nvertl; k++)
          {     
              if(abs(x-xolddown[k]) < diffdown)
              {
                  diffdown = abs(x-xolddown[k]);
                  qp_closerolddown=k;
              }
              if(abs(x-xoldup[k]) < diffup)
              {
                  diffup = abs(x-xoldup[k]);
                  qp_closeroldup=k;
              }                 
     
          }
          
          //find nplay_closer
          diffdown =1000;
          diffup = 1000;
          
          int qp_closerup, qp_closerdown;
          
          for(int f=0; f< nvertl; f++)
          {
              if(abs(x-xnew_down[f]) < diffdown)
              {
                  diffdown = abs(x-xnew_down[f]);
                  qp_closerdown=f;
              }
              if(abs(x-xnew_up[f]) < diffup)
              {
                  diffup = abs(x-xnew_up[f]);
                  qp_closerup=f;
              }   
          }
          
          int qp_closernormoldup;
          Vmath::Sadd(nvertl, -x,xoldup,1, tmp,1);
          Vmath::Vmul(nvertl, tmp,1,tmp,1,tmp,1);
          Vmath::Sadd(nvertl,-y,yoldup,1,norm,1);
          Vmath::Vmul(nvertl, norm,1,norm,1,norm,1);
          Vmath::Vadd(nvertl, norm,1,tmp,1,norm,1);
          qp_closernormoldup = Vmath::Imin(nvertl, norm,1);
          
          Vmath::Zero(nvertl, norm,1);
          Vmath::Zero(nvertl, tmp,1);          
          
          int qp_closernormolddown;          
          Vmath::Sadd(nvertl, -x,xolddown,1, tmp,1);
          Vmath::Vmul(nvertl, tmp,1,tmp,1,tmp,1);
          Vmath::Sadd(nvertl,-y,yolddown,1,norm,1);
          Vmath::Vmul(nvertl, norm,1,norm,1,norm,1);
          Vmath::Vadd(nvertl, norm,1,tmp,1,norm,1);     
          qp_closernormolddown = Vmath::Imin(nvertl, norm,1);
          
          Vmath::Zero(nvertl, norm,1);
          Vmath::Zero(nvertl, tmp,1); 

          int qp_closernormup; 
          Vmath::Sadd(nvertl, -x,xnew_up,1, tmp,1);
          Vmath::Vmul(nvertl, tmp,1,tmp,1,tmp,1);
          Vmath::Sadd(nvertl,-y,ynew_up,1,norm,1);
          Vmath::Vmul(nvertl, norm,1,norm,1,norm,1);
          Vmath::Vadd(nvertl, norm,1,tmp,1,norm,1);
          qp_closernormup = Vmath::Imin(nvertl, norm,1);  
          
          Vmath::Zero(nvertl, norm,1);
          Vmath::Zero(nvertl, tmp,1); 

          int qp_closernormdown; 
          Vmath::Sadd(nvertl, -x,xnew_down,1, tmp,1);
          Vmath::Vmul(nvertl, tmp,1,tmp,1,tmp,1);
          Vmath::Sadd(nvertl,-y,ynew_down,1,norm,1);
          Vmath::Vmul(nvertl, norm,1,norm,1,norm,1);
          Vmath::Vadd(nvertl, norm,1,tmp,1,norm,1);
          qp_closernormdown = Vmath::Imin(nvertl, norm,1);           
          
          


          if(  y>yoldup[qp_closeroldup] && y<1 )
          {	        

//              ratio = (1-layers_y[nlays-1][qp_closer])*(1-y_c[qp_closer])/
//                    (  (1-yold_up[n])*(1-yold_c[qp_closer]) );
              ratio = (1-ynew_up[qp_closerup])/
                    (  (1-yoldup[qp_closeroldup]) );
                    
              ratioy = (1-ynew_up[qp_closernormup])/
                    (  (1-yoldup[qp_closernormoldup]) ); 
              //distance prop to layerup
              ynew[n] = ynew_up[qp_closerup] 
                      + (y-yoldup[qp_closeroldup])*ratio;  
              //ynew[n] = y +abs(nyvert[qp_closernormup])*(ynew_up[qp_closeroldup]-yoldup[qp_closeroldup])*ratioy;
                                            
              //ynew[n] = y + 0.3*(ynew_up[qp_closerup]-yoldup[qp_closerup]);
              //xnew[n] = x + abs(nxvert[qp_closeroldup])*(xnew_up[qp_closeroldup]-xoldup[qp_closeroldup]);
              
              if(x> (xmax-xmin)/2. && x< xmax)
              {
          	  ratiox = (xmax-xnew_up[qp_closernormup])/
          	           (xmax-xoldup[qp_closernormup])  ;
          	  if( (xmax-xoldup[qp_closernormup])==0)
          	  {
          	        ratiox = 1.0;	  
          	  }          	           
          	        
          	  //xnew[n] = xnew_up[qp_closerup]
          	  //        + (x-xoldup[qp_closerup])*ratiox;
              xnew[n] = x + abs(nxvert[qp_closernormup])*(xnew_up[qp_closeroldup]-xoldup[qp_closeroldup])*ratiox;
              ASSERTL0(x>xmin," x value <xmin up second half");
              ASSERTL0(x<xmax," x value >xmax up second  half");               
              }
              else if( x> xmin && x<= (xmax-xmin)/2.)
              {
//cout<<"up  close normold="<<qp_closernormoldup<<"   closenorm="<<qp_closernormup<<endl;              	      
          	  ratiox = (xnew_up[qp_closernormup]-xmin)/
          	        (  (xoldup[qp_closernormup]-xmin)  );
          	  if( (xoldup[qp_closernormup]-xmin)==0)
          	  {
          	        ratiox = 1.0;	  
          	  }
          	  //xnew[n] = xnew_up[qp_closerup]
          	  //        + (x-xoldup[qp_closeroldup])*ratiox;
              xnew[n] = x + abs(nxvert[qp_closernormup])*(xnew_up[qp_closeroldup]-xoldup[qp_closeroldup])*ratiox;  
//cout<<"up xold="<<x<<"  xnew="<<xnew[n]<<endl;              
              ASSERTL0(x>xmin," x value <xmin up first half");
              ASSERTL0(x<xmax," x value >xmax up first half");               
              }   
          	      
           
          }
          else if(   y< yolddown[qp_closerolddown]   && y>-1  )
          {

              ratio = (1+ynew_down[qp_closerdown])/
                    (  (1+yolddown[qp_closerolddown]) );
                    
              ratioy = (1-ynew_down[qp_closernormdown])/
                    (  (1-yolddown[qp_closernormolddown]) );      
               
              //distance prop to layerlow
              ynew[n] = ynew_down[qp_closerdown] 
                      + (y-yolddown[qp_closerolddown])*ratio;
              //ynew[n] = y +abs(nyvert[qp_closernormdown])*
              // (ynew_down[qp_closerolddown]-yolddown[qp_closerolddown])*ratioy;              
              //ynew[n] = y + 0.3*(ynew_down[qp_closerdown]-yolddown[qp_closerdown]);
              //xnew[n] = x + abs(nxvert[qp_closerolddown])*(xnew_down[qp_closerolddown]-xolddown[qp_closerolddown]);
/*              
if(n==74)
{
cout<<qp_closerolddown<<"    nplaydown="<<qp_closerdown<<endl;	
cout<<"xolddown="<<xolddown[qp_closerolddown]<<"   xnewdown="<<xnew_down[qp_closerdown]<<endl;
cout<<"xold+"<<x<<"   xnew+"<<xnew[n]<<endl;
} 
*/
              
              if(x> (xmax-xmin)/2.  && x <xmax)
              {
          	  ratiox = (xmax-xnew_down[qp_closernormdown])/
          	        (  (xmax-xolddown[qp_closernormdown])  );
          	  if( (xmax-xolddown[qp_closernormdown])==0)
          	  {
          	        ratiox = 1.0;	  
          	  }            	        
          	  //xnew[n] = xnew_down[qp_closerdown]
          	  //        + (x-xolddown[qp_closerolddown])*ratiox;
              xnew[n] = x + 
              abs(nxvert[qp_closernormdown])*(xnew_down[qp_closerolddown]-xolddown[qp_closerolddown])*ratiox;  
              ASSERTL0(x>xmin," x value <xmin down second half");
              ASSERTL0(x<xmax," x value >xmax down second half");               
              }
              else if( x>xmin  && x<= (xmax-xmin)/2.)
              {
          	  ratiox = (xnew_down[qp_closernormdown]-xmin)/
          	        (  (xolddown[qp_closernormdown]-xmin)  );
          	  if( (xolddown[qp_closernormdown]-xmin)==0)
          	  {
          	        ratiox = 1.0;	  
          	  }          	        
          	  //xnew[n] = xnew_down[qp_closerdown]
          	  //        + (x-xolddown[qp_closerolddown])*ratiox;
              xnew[n] = x + 
              abs(nxvert[qp_closernormdown])*(xnew_down[qp_closerolddown]-xolddown[qp_closerolddown])*ratiox;   
              ASSERTL0(x>xmin," x value <xmin down first half");
              ASSERTL0(x<xmax," x value >xmax down first half");  
              }
              
          }
          
cout<<"xold"<<x<<"   xnew="<<xnew[n]<<endl;       
     ASSERTL0(xnew[n] >= xmin, "newx < xmin");
     ASSERTL0(xnew[n]<= xmax, "newx > xmax");

     }// verts closed     	
}

void CheckSingularQuads( MultiRegions::ExpListSharedPtr Exp, 
                 Array<OneD, int> V1, Array<OneD, int> V2,	         
	         Array<OneD, NekDouble>& xnew,Array<OneD, NekDouble>& ynew)
{
      const boost::shared_ptr<StdRegions::StdExpansionVector> exp2D = Exp->GetExp();
      int nel        = exp2D->size();
      LocalRegions::QuadExpSharedPtr locQuadExp;
      LocalRegions::TriExpSharedPtr  locTriExp;
      SpatialDomains::Geometry1DSharedPtr SegGeom;
      int idbef, idnext;
      NekDouble xV1, yV1, xV2,yV2;
      NekDouble slopebef,slopenext,slopenew;
      Array<OneD, int> locEids(4);
      for(int i=0; i<nel; i++)
      { 
           if((locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>((*exp2D)[i])))
           {
                SegGeom = (locQuadExp->GetGeom2D())->GetEdge(0);
                idbef = SegGeom->GetEid();
                if(xnew[  V1[idbef] ]<= xnew[  V2[idbef] ])
                {
                xV1 = xnew[  V1[idbef] ];
                yV1 = ynew[  V1[idbef] ];
                xV2 = xnew[  V2[idbef] ];
                yV2 = ynew[  V2[idbef] ];
                slopebef = (yV2 -yV1)/(xV2 -xV1);
                }
                else
                {
                xV1 =  xnew[  V2[idbef] ];
                yV1 =  ynew[  V2[idbef] ];
                xV2 = xnew[  V1[idbef] ];
                yV2 = ynew[  V1[idbef] ];
                slopebef = (yV2 -yV1)/(xV2 -xV1);                	
                }
//cout<<"00 V1 x="<<xnew[  V1[idbef] ]<<"   y="<<ynew[  V1[idbef] ]<<endl;
//cout<<"00 V2 x="<<xnew[  V2[idbef] ]<<"   y="<<ynew[  V2[idbef] ]<<endl;                
                for(int j = 1; j < locQuadExp->GetNedges(); ++j)
                {
                    SegGeom = (locQuadExp->GetGeom2D())->GetEdge(j);
                    idnext = SegGeom->GetEid();
//cout<<"id="<<idnext<<" locid="<<j<<endl;
//cout<<" V1 x="<<xnew[  V1[idnext] ]<<"   y="<<ynew[  V1[idnext] ]<<endl;
//cout<<" V2 x="<<xnew[  V2[idnext] ]<<"   y="<<ynew[  V2[idnext] ]<<endl;
                    if(xV1 == xnew[  V1[idnext] ] && yV1 == ynew[  V1[idnext] ]  )
                    {
                    xV1 = xnew[  V1[idnext] ];
                    yV1 = ynew[  V1[idnext] ];
                    xV2 = xnew[  V2[idnext] ];
                    yV2 = ynew[  V2[idnext] ];
                    slopenext = (yV2 -yV1)/(xV2 -xV1);
if(i==23)
{
cout<<"case1 x0="<<xV1<<"   x1="<<xV2<<endl;
cout<<idnext<<"  11slope bef ="<<slopebef<<"  slopenext="<<slopenext<<endl;	
}
                          //compare with slope before
                          if( slopebef/slopenext>0.84 && slopebef/slopenext <1.18)
                          {
                               xnew[ V1[idnext]  ] =  xnew[ V1[idnext]  ] -0.01;
                               slopenew = (yV2-yV1)/(xV2- xnew[ V1[idnext]  ]);

                               if( abs(slopebef-slopenew) < abs(slopebef-slopenext) )
                               {
                                      xnew[ V1[idnext]  ] =  xnew[ V1[idnext]  ] +0.02;
                                      slopenew = (yV2-yV1)/(xV2- xnew[ V1[idnext]  ]);
                               }
                               slopenext = slopenew;
cout<<"slopenew="<<slopenew<<endl;                                       
cout<<"moved x="<<xnew[ V1[idnext]  ]<<endl;	    
                          }
                    }
                    else if(xV2 == xnew[  V2[idnext] ] && yV2 == ynew[  V2[idnext] ]  )
                    {
                    xV1 = xnew[  V2[idnext] ];
                    yV1 = ynew[  V2[idnext] ];
                    xV2 = xnew[  V1[idnext] ];
                    yV2 = ynew[  V1[idnext] ];
                    slopenext = (yV2 -yV1)/(xV2 -xV1);	 
if(i==23)
{
cout<<"case2 x0="<<xV1<<"   x1="<<xV2<<endl;
cout<<idnext<<"  22slope bef ="<<slopebef<<"  slopenext="<<slopenext<<endl;	
}                    
                          //compare with slope before
                          if( slopebef/slopenext>0.84 && slopebef/slopenext <1.18)
                          {
                               xnew[ V2[idnext]  ] = xnew[ V2[idnext]  ] -0.01;   
                               slopenew = (yV2-yV1)/(xV2- xnew[ V2[idnext]  ]);

                               if( abs(slopebef-slopenew) < abs(slopebef-slopenext) )
                               {
                                      xnew[ V2[idnext]  ] =  xnew[ V2[idnext]  ] +0.02;
                                      slopenew = (yV2-yV1)/(xV2- xnew[ V2[idnext]  ]);

                               }
                               slopenext = slopenew;
cout<<"slopenew="<<slopenew<<endl;
cout<<"moved x="<<xnew[ V2[idnext]  ]<<endl;	    
                          }                    
                    }
                    else if(xV1 == xnew[ V2[idnext] ] && yV1 == ynew[  V2[idnext] ]  )
                    {
                    xV1 = xnew[  V2[idnext] ];
                    yV1 = ynew[  V2[idnext] ];
                    xV2 = xnew[  V1[idnext] ];
                    yV2 = ynew[  V1[idnext] ];
                    slopenext = (yV2 -yV1)/(xV2 -xV1);	 
if(i==23)
{
cout<<"case3 x0="<<xV1<<"   x1="<<xV2<<endl;
cout<<idnext<<"  22slope bef ="<<slopebef<<"  slopenext="<<slopenext<<endl;	
}                    
                          //compare with slope before
                          if( slopebef/slopenext>0.84 && slopebef/slopenext <1.18)
                          {
                               xnew[ V2[idnext]  ] = xnew[ V2[idnext]  ] -0.01;   
                               slopenew = (yV2-yV1)/(xV2- xnew[ V2[idnext]  ]);

                               if( abs(slopebef-slopenew) < abs(slopebef-slopenext) )
                               {
                                      xnew[ V2[idnext]  ] =  xnew[ V2[idnext]  ] +0.02;
                                      slopenew = (yV2-yV1)/(xV2- xnew[ V2[idnext]  ]);
                               }
                               slopenext = slopenew;	
cout<<"slopenew="<<slopenew<<endl;
cout<<"moved x="<<xnew[ V2[idnext]  ]<<endl;
                          }      

                    }
                    else if(xV2 == xnew[ V1[idnext] ] && yV2 == ynew[  V1[idnext] ]  )
                    {
                    xV1 = xnew[  V1[idnext] ];
                    yV1 = ynew[  V1[idnext] ];
                    xV2 = xnew[  V2[idnext] ];
                    yV2 = ynew[  V2[idnext] ];
                    slopenext = (yV2 -yV1)/(xV2 -xV1);	 
if(i==23)
{
cout<<"case4 x0="<<xV1<<"   x1="<<xV2<<endl;
cout<<idnext<<"  22slope bef ="<<slopebef<<"  slopenext="<<slopenext<<endl;	
}                    
                          //compare with slope before
                          if( slopebef/slopenext>0.84 && slopebef/slopenext <1.18)
                          {
                               xnew[ V1[idnext]  ] = xnew[ V1[idnext]  ] -0.01;  
                               slopenew = (yV2-yV1)/(xV2- xnew[ V1[idnext]  ]);

                               if( abs(slopebef-slopenew) < abs(slopebef-slopenext) )
                               {
                                      xnew[ V1[idnext]  ] =  xnew[ V1[idnext]  ] +0.02;
                                      slopenew = (yV2-yV1)/(xV2- xnew[ V1[idnext]  ]);
                               } 
                               slopenext = slopenew;
cout<<"slopenew="<<slopenew<<endl;
cout<<"moved x="<<xnew[ V1[idnext]  ]<<endl;
                          }      

                    }
                    else
                    {
                          ASSERTL0(false, "edge not connected");	    
                    }                    
                    slopebef = slopenext;

                    
                    
                }         
           }  
      }
}

void Replacevertices(string filename, Array<OneD, NekDouble> newx, 
    	               Array<OneD,  NekDouble> newy,
    	               Array<OneD, NekDouble> xcPhys, Array<OneD, NekDouble> ycPhys,
    	               Array<OneD, int>Eids, int Npoints, string s_alp,
                       Array<OneD, Array<OneD, NekDouble> > x_lay,
                       Array<OneD, Array<OneD, NekDouble> > y_lay,
                       Array<OneD, Array<OneD, int > >lay_eids, bool curv_lay)	
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
            xscale = expEvaluator.Evaluate(expr_id);
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
       TiXmlElement* masternew = NULL;    
       TiXmlElement* condnew = NULL; 
       TiXmlElement* Parsnew = NULL; 
       TiXmlElement* parnew = NULL; 

       // Master tag within which all data is contained.

      
       masternew = docnew.FirstChildElement("NEKTAR");
       ASSERTL0(masternew, "Unable to find NEKTAR tag in file.");

       //set the alpha value
       string alphastring;
       condnew = masternew->FirstChildElement("CONDITIONS");
       Parsnew = condnew->FirstChildElement("PARAMETERS");
cout<<"alpha="<<s_alp<<endl;
       parnew = Parsnew->FirstChildElement("P");
       while(parnew)
       {
            TiXmlNode *node = parnew->FirstChild();
            if (node)
            {
                  // Format is "paramName = value"
                  std::string line = node->ToText()->Value();
                  std::string lhs;
                  std::string rhs;  
                  /// Pull out lhs and rhs and eliminate any spaces.
                  int beg = line.find_first_not_of(" ");
                  int end = line.find_first_of("=");
                  // Check for no parameter name
                  if (beg == end) throw 1;
                  // Check for no parameter value
                  if (end != line.find_last_of("=")) throw 1;
                  // Check for no equals sign
                  if (end == std::string::npos) throw 1;
                  lhs = line.substr(line.find_first_not_of(" "), end-beg);
                  lhs = lhs.substr(0, lhs.find_last_not_of(" ")+1);

                      //rhs = line.substr(line.find_last_of("=")+1);
                      //rhs = rhs.substr(rhs.find_first_not_of(" "));
                      //rhs = rhs.substr(0, rhs.find_last_not_of(" ")+1);  

                  boost::to_upper(lhs); 
                  if(lhs == "ALPHA")
                  {
                      alphastring = "Alpha   =  "+ s_alp;
                      parnew->RemoveChild(node);
                      parnew->LinkEndChild(new TiXmlText(alphastring) );                           
                  }     
             }
                 
	     parnew = parnew->NextSiblingElement("P");  
      }


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
      int neids_lay = lay_eids[0].num_elements();
      //if edgenew belongs to the crit lay replace it, else delete it.
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
//cout<<"Eids="<<Eids[u]<<"  eid="<<eid<<endl;	           	   
	         if(Eids[u]==eid)
	         {
	             indexeid = u;
	         }	               
	    }
	    if(indexeid==-1)
            {
                  curvednew->RemoveChild(edgenew);
   	          //ASSERTL0(false, "edge to update not found");	   
   	    }
            else
            {

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
		  //we need at least 5 digits (setprecision 5) to get the streak position with
		  // precision 10^-10
		   
                  //Determine the number of points
                  err = edgenew->QueryIntAttribute("NUMPOINTS", &numPts);
                  ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute NUMPOINTS.");

                  stringstream st;
                  edgenew->SetAttribute("NUMPOINTS", Npoints);
                  for(int u=0; u< Npoints; u++)
                  {
                        st << std::scientific << 
                        std::setprecision(8) <<xcPhys[cnt*Npoints+u]
                             << "   " << ycPhys[cnt*Npoints+u] << "   " << 0.000<<"   ";
                  }

                  edgenew->LinkEndChild(new TiXmlText(st.str())); 

                        
/*
		        st << std::scientific << std::setprecision(8) << x_crit[v1] << "   "
	    	        << y_crit[v1] << "   " << 0.000<<"   ";
		        for(int a=0; a< Npoints-2; a++)
                        {		       	   
		             st << std::scientific << std::setprecision(8) <<
		             "    "<<Pcurvx[indexeid*(Npoints-2) +a]<<"    "<<Pcurvy[indexeid*(Npoints-2) +a]
		             <<"    "<<0.000<<"   "; 
	      
		        }
		        st << std::scientific << std::setprecision(8) <<
		        "    "<<x_crit[v2]<<"   "<< y_crit[v2] <<"   "<< 0.000;	
                             edgenew->LinkEndChild(new TiXmlText(st.str()));
  
*/	   
 

            }
                
   	    edgenew = edgenew->NextSiblingElement("E");
   	    	    
      }

      //write also the others layers curve points
      if(curv_lay == true)
      {
cout<<"write other curved edges"<<endl;
            TiXmlElement * curved = meshnew->FirstChildElement("CURVED");
            int idcnt = 300;
            int nlays = lay_eids.num_elements();

            //TiXmlComment * comment = new TiXmlComment();
            //comment->SetValue("   new edges  ");
            //curved->LinkEndChild(comment);
            for (int g=0; g< nlays; ++g)  
            {
                  for(int p=0; p< neids_lay; p++)
                  {
                        stringstream st;
                        TiXmlElement * e = new TiXmlElement( "E" );
                        e->SetAttribute("ID",        idcnt++);
                        e->SetAttribute("EDGEID",    lay_eids[g][p]);
                        e->SetAttribute("NUMPOINTS", Npoints);
                        e->SetAttribute("TYPE", "PolyEvenlySpaced");
                        for(int c=0; c< Npoints; c++)
                        {
                             st << std::scientific << std::setprecision(8) <<x_lay[g][p*Npoints +c]
                             << "   " << y_lay[g][p*Npoints +c] << "   " << 0.000<<"   ";

                        }

                        TiXmlText * t0 = new TiXmlText(st.str());
                        e->LinkEndChild(t0);
                        curved->LinkEndChild(e);
                  }
 
             }  	    
          }	    
   	    
   	    
   	    
       	  docnew.SaveFile( newfile ); 
       	    
       	  cout<<"new file:  "<<newfile<<endl;
}       	  
