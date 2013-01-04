#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <boost/lexical_cast.hpp>
#include <tinyxml/tinyxml.h>


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
    	        Array<OneD, NekDouble> x,  Array<OneD, NekDouble> y,
    	        Array<OneD, NekDouble> &xold_up, Array<OneD, NekDouble> &yold_up,
    	        Array<OneD, NekDouble> &xold_low, Array<OneD, NekDouble> &yold_low,
    	        Array<OneD, NekDouble> &xold_c, Array<OneD, NekDouble> &yold_c,      	        
    	        Array<OneD, NekDouble> &xc,  Array<OneD, NekDouble> &yc, NekDouble cr,
                bool verts); 
    void GenerateAddPointsNewtonIt( NekDouble &xi, NekDouble &yi,NekDouble &x0, NekDouble &y0,
    	        MultiRegions::ExpListSharedPtr &function, Array<OneD, NekDouble> &derfunction,
                NekDouble cr);
    void GenerateCurveSp(Array<OneD, NekDouble> &x_c, Array<OneD, NekDouble> &y_c,
                Array<OneD, NekDouble> &x_lay, Array<OneD, NekDouble> &y_lay);
    void GenerateCurve(int npoints, int npused, Array<OneD, NekDouble> &x_c, 
    	    Array<OneD, NekDouble> &y_c, Array<OneD, NekDouble> &curve);
    void GenerateAddPoints(int region, SpatialDomains::MeshGraphSharedPtr &mesh,int np, int npused,    	      
    	         Array<OneD, NekDouble> &curve, MultiRegions::ExpListSharedPtr & bndfield, 
    	         Array<OneD, NekDouble>& outx, Array<OneD, NekDouble>& outy,
    	         Array<OneD, int>&Eids);
    void MapEdgeVertices(MultiRegions::ExpListSharedPtr field, Array<OneD, int> &V1,
                 Array<OneD, int> &V2); 
    void Findlay_eids(Array<OneD, Array<OneD, NekDouble> >& layers_y,
                 Array<OneD, int> V1, Array<OneD, int> V2, 
                 Array<OneD, Array<OneD, int > >& lay_Vids,
                 Array<OneD, Array<OneD, int > >& lay_eids);
    void MoveLayersvertically(int nlays, int nvertl, int cntlow, int cntup,
    	         Array<OneD, Array<OneD, int > > lay_Vids,  Array<OneD, NekDouble> xc,
    	         Array<OneD, NekDouble> yc, Array<OneD, int> Down, Array<OneD, int> Up,
    	         Array<OneD, NekDouble >& xnew, Array<OneD, NekDouble>& ynew,
    	         Array<OneD, Array<OneD, NekDouble > >& layers_x,
    	         Array<OneD, Array<OneD, NekDouble > >& layers_y);    	         
    void  Cutrepetitions(int nedges,Array<OneD, NekDouble> inarray,
                 Array<OneD, NekDouble>& outarray);
    void  Orderfunctionx(Array<OneD, NekDouble> inarray_x,
                 Array<OneD, NekDouble> inarray_y, Array<OneD, NekDouble>& outarray_x,
                 Array<OneD, NekDouble>& outarray_y);
    int DetermineclosePointxindex(NekDouble x,Array<OneD, NekDouble> xArray);
    void GenerateNeighbourArrays(int index, int neighpoints,Array<OneD, NekDouble> xArray,
                 Array<OneD, NekDouble> yArray,Array<OneD, NekDouble>& Neighbour_x,
                 Array<OneD, NekDouble>& Neighbour_y);
    NekDouble LagrangeInterpolant(NekDouble x, int npts, 
                 Array<OneD,NekDouble>  xpts, Array<OneD, NekDouble> funcvals);
    void ChangeLayerspos(Array<OneD, NekDouble> & ynew,
                 Array<OneD, NekDouble>  yc,
                 Array<OneD, NekDouble> Addpointsy,
                 Array<OneD, Array<OneD, NekDouble > >& layers_y,
                 Array<OneD, Array<OneD, int> >lay_Vids,
                 NekDouble delt, NekDouble delt_opp,string movelay, int npedge);

    void Replacevertices(string filename, Array<OneD, NekDouble> newx, 
    	               Array<OneD,  NekDouble> newy,
    	               Array<OneD, NekDouble> x_crit, Array<OneD, NekDouble> y_crit,
    	               Array<OneD, NekDouble> & Pcurvx, 
    	               Array<OneD, NekDouble> & Pcurvy,
    	               Array<OneD, int>&Eids, int Npoints, string s_alp,
                       Array<OneD, Array<OneD, NekDouble> > x_lay,
                       Array<OneD, Array<OneD, NekDouble> > y_lay,
                       Array<OneD, Array<OneD, int > >lay_eids, bool curv_lay);
    
    
    
    int i,j;

//ATTEnTION !!! with argc=2 you impose that vSession refers to is argv[1]=meshfile!!!!! 
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(2, argv);
    //----------------------------------------------
    NekDouble cr;
    //set cr =0
    cr=0;
    //change argc from 6 to 5 allow the loading of cr to be optional
    if(argc > 6 || argc <5 )
    {
        fprintf(stderr,
            "Usage: ./MoveMesh  meshfile fieldfile  changefile   alpha  cr(optional)\n");
        exit(1);
    }
    else if( argc == 6 &&
               vSession->DefinesSolverInfo("INTERFACE")
               && vSession->GetSolverInfo("INTERFACE")=="phase" )
    {
        cr = boost::lexical_cast<double>(argv[argc-1]);
        argc=5;
    }

    // Read in mesh from input file
    string meshfile(argv[argc-4]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //---------------------------------------------- 

    // Also read and store the boundary conditions
    SpatialDomains::BoundaryConditionsSharedPtr boundaryConditions;        
    boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                                        ::AllocateSharedPtr(vSession,graphShPt);
    SpatialDomains::BoundaryConditions bcs(vSession, graphShPt); 
    //----------------------------------------------

    // Define Expansion   
    Array<OneD, MultiRegions::ExpListSharedPtr> fields;   
    int nfields;  
cout<<"o2k"<<endl;   
    //the mesh file should have 2 component: set output fields
    //fields has to be of the SAME dimension of the mesh (that's why there is
    //the changefile as an input) 
    //a contfield2D is needed to extract boundary conditions!!!
    nfields=3;        
    SetFields(graphShPt,boundaryConditions,vSession,fields,nfields);
    //---------------------------------------------------------------
cout<<"ok"<<endl;
    // store name of the file to change
    string changefile(argv[argc-2]);
    //----------------------------------------------
      
    //store the value of alpha
    string charalp (argv[argc-1]);
    //NekDouble alpha = boost::lexical_cast<double>(charalp);
cout<<"read alpha="<<charalp<<endl;
    // Import field file.
    string fieldfile(argv[argc-3]);
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
        streak->ExtractDataToCoeffs(fielddef[i], fielddata[i], fielddef[i]->m_fields[0], streak->UpdateCoeffs());
    }
    streak->BwdTrans_IterPerExp(streak->GetCoeffs(), streak->UpdatePhys());
 
    //------------------------------------------------  
/*  
static int cnt=0;
int nquad = streak->GetTotPoints();
Array<OneD, NekDouble> xs(nquad);
Array<OneD, NekDouble> ys(nquad);
streak->GetCoords(xs,ys);
if(  
//abs(streak->GetPhys()[j])< 0.01
x[j]==1.6
)
{
cnt++;
cout<<"cnt="<<cnt<<"   x="<<xs[j]<<"  y="<<ys[j]<<"  U="<<streak->GetPhys()[j]<<endl;
}           
*/
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
    NekDouble x0,y0,z0,xt,yt,zt;
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
    //first point for x_connect=0(or-1.6)
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
         
    //-----------------------------------------------------------------------------------   

 
    //fieds to force continuity:
/*
    MultiRegions::ContField1DSharedPtr  Lay_x;    		
    MultiRegions::ContField1DSharedPtr  Lay_y;
    //initialie fields
    const SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
    MultiRegions::ExpList1DSharedPtr xtmp;    		
    MultiRegions::ExpList1DSharedPtr ytmp;
    xtmp = MemoryManager<MultiRegions::ExpList1D>
                ::AllocateSharedPtr(*(bregions[lastIregion]), graphShPt, true);    
    ytmp = MemoryManager<MultiRegions::ExpList1D>
    		::AllocateSharedPtr(*(bregions[lastIregion]), graphShPt, true);

    Lay_x = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(vSession, *xtmp);

    Lay_y = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(vSession, *ytmp);  
*/
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
    fields[0]->GetCoords(x,y);         
    Array<OneD, NekDouble> x_c(nvertl);
    Array<OneD,NekDouble> y_c(nvertl,-200);       
    Array<OneD, NekDouble> tmp_w (nvertl, 200);
    Array<OneD, int> Sign (nvertl,1);   
    Array<OneD, NekDouble> Delta_c(nvertl,-200);
   

    Computestreakpositions(nvertl, streak, x,y,xold_up, yold_up,
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
         //fill layer arrays:
/*
         if(q==0  || q == nvertl-1)
         {
             Lay_x->UpdatePhys()[q]=x_c[q];
             Lay_y->UpdatePhys()[q]=y_c[q];
         }
         else if(q<nvertl && q!=nvertl-1)
         {
             Lay_x->UpdatePhys()[ q*(npedge-1) ]=x_c[q];
             Lay_y->UpdatePhys()[ q*(npedge-1) ]=y_c[q];
             Lay_x->UpdatePhys()[ q*(npedge-1) +1 ]=x_c[q];
             Lay_y->UpdatePhys()[ q*(npedge-1) +1 ]=y_c[q];
         }
*/
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

/*
    int ncurves=1;

    int ntotcurvecoeffs=0;
    //array storing ALL the curvecoeffs 
    //hypothesis: totcoeffs< 2*nvertl; ncurves< nedges;
    Array<OneD, NekDouble> totcurvecoeffs (2*nvertl,0.0);
    Array<OneD, int> ncurvepoints (nedges);

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
    //int npedge=5;
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
   
    //Array<OneD, NekDouble> Ycoords (nedges*npedge-2, 0.0);

 
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

    //force the continuity of the layer

/*
cout<<"force"<<endl;
    int Nregcoeffs = Lay_x->GetNcoeffs();
    Array<OneD, NekDouble> coeffs (Nregcoeffs);
    //Lay_x->FwdTrans(Lay_x->GetPhys(), coeffs);
    //Lay_x->BwdTrans(coeffs, Lay_x->UpdatePhys());      
    //Array<OneD, NekDouble> dery(npedge*nedges);
    //Lay_y->PhysDeriv(MultiRegions::eX, Lay_y->GetPhys(), dery);


    Vmath::Zero(Nregcoeffs,coeffs,1);            

    //Lay_y->FwdTrans(Lay_y->GetPhys(), coeffs);
    //Lay_y->BwdTrans(coeffs, Lay_y->UpdatePhys());  
    for(int r=0; r<nedges; r++)
    {
      
          for(int w=0; w< npedge-2; w++)
          {         
             //   Addpointsx[r*(npedge-2) +w] = Lay_x->GetPhys()[r*npedge +1 +w];
             //   Addpointsy[r*(npedge-2) +w] = Lay_y->GetPhys()[r*npedge +1 +w];
          }
    }    
*/
    Array<OneD, NekDouble> curve(nvertl);
cout<<"gen curve"<<endl;
    //GenerateCurve(nvertl,0, x_c, y_c, curve);


    //-------------------------------------------------



/*  
    //generate the closest curve to the critical layer positions taking 
    // npedge-2 points per edge
    Array<OneD, NekDouble> Crlay_pointsx(nedges*(npedge-2)+nedges+1);
    Array<OneD, NekDouble> Crlay_pointsy(nedges*(npedge-2)+nedges+1);
    int cnt=0;
    int cnt1=0;
    int cnt2=0;
    //HYPOTHESIS: Cpoints ALREADY generated
//cout<<"Cr x--y"<<endl;    
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
*/
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

    // determine the number of layers
    int nlays=0;

  
    int cnt=0;
    //find coords and Vids of the layers vertices
    int nVertTot = graphShPt->GetNvertices();
    //arrays to check with the new values
    Array<OneD, NekDouble> xold(nVertTot);
    Array<OneD, NekDouble> yold(nVertTot);
    cnt=0;
    for(int n=0; n<nVertTot; n++)
    {
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(n);
         NekDouble x,y,z;
         vertex->GetCoords(x,y,z); 
         xold[n] =x;
         yold[n] =y;
         //tollerance 5e-2 to get the points
         if(x<= xold_c[0] +5e-2 && x >= xold_c[0] -5e-2
         && y<= yold_up[0] && y>=yold_low[0]
         && y!= yold_c[0])
         {     
             cnt++;
//cout<<"vert index="<<m<<"   vertID="<<n<<endl;                  
         }
    }

    nlays =cnt;
cout<<"nlays="<<nlays<<endl;
    Array<OneD, Array<OneD, NekDouble> > layers_y(nlays);
    Array<OneD, Array<OneD, NekDouble> > layers_x(nlays);
    Array<OneD, Array<OneD, int > >lay_Vids(nlays);    
    Array<OneD, Array<OneD, int > >lay_eids(nlays);  
    //initialise layers_y,lay_eids
    for(int g=0; g<nlays; g++)
    {
        layers_y[g]= Array<OneD, NekDouble> ( (nvertl-1)*npedge );
        lay_Vids[g]= Array<OneD, int> (nvertl);
        lay_eids[g]= Array<OneD, int> (nvertl-1);
    }
  
    //------------------------------------------------------------
    
    //calculate the new cordinates of the vertices    
    Array<OneD, NekDouble> xnew(nVertTot);
    Array<OneD, NekDouble> ynew(nVertTot,-20);
    Array<OneD, int> Up(nvertl);//Vids lay Up
    Array<OneD, int> Down(nvertl);//Vids lay Down
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
    int cntup=0;
    int cntlow=0;
    Array<OneD, int> Acntlay(nvertl,0);
    //Vids needed only if a layers has to be moved
    NekDouble bleft=-10;
    NekDouble tright = 10;
    NekDouble bright = -10;
    NekDouble tleft  = 10;
    cnt=0;
    int bottomleft, topright,bottomright,topleft;
    for(int i=0; i<nVertTot; i++)
    {
         bool mvpoint =false;
       	 NekDouble ratio;  
         SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(i);
         NekDouble x,y,z;
         vertex->GetCoords(x,y,z); 
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
//cout<<i<<"        "<<xnew[i]<<"     "<<ynew[i]<<endl;     
         //count number of layers nlays(again) initialise vectors
         //(crit lay not incl but uplay and down lay included)
         //nb tollerance 10-5 to get the points
         for(int m=0; m<nvertl; m++)
         {
 
             if( x>= x_c[m] -5e-2 && x<= x_c[m] +5e-2 
              && y<= yold_up[m] && y>=yold_low[m]
              && y!= yold_c[m])
             {     
//cout<<"LAY  x"<<x<<"  y="<<y<<endl;
//cout<<"yup="<<yold_up[4]<<"  ylow="<<yold_low[4]<<endl;

                 //need to fill layers_y to order it(through Findlay_eids)
                 layers_y[Acntlay[m]][m]= y;
                 lay_Vids[Acntlay[m]][m]= i;
                 
                 Acntlay[m]++;
             }

         }
             
    }
cout<<"cntlow="<<cntlow<<endl;
    ASSERTL0(Acntlay[4]==nlays, "something wrong with the number of layers");
    //order layers coords:
    //V1[eid],V2[eid] vertices associate with the edge Id=eid
    Array<OneD, int> V1;
    Array<OneD, int> V2;    
    MapEdgeVertices(streak,V1,V2);
    //find eids which belong to each layer(essentially order the 
    // lay_Vids[n][m] on m layers)
    Findlay_eids(layers_y, V1, V2, lay_Vids, lay_eids);
    //determine the new coords of the vertices and the curve points 
    //for each edge
 

    //------------------------------------------------------------------

    
    
    
    
/*    
//possible interface to alglib!!!!!!!!
    Array<OneD, NekDouble> xalg(5, 1.0); //= "[-1.0,-0.5,0.0,+0.5,+1.0]";
    Array<OneD, NekDouble> yalg(5, 1.0); //= "[+1.0,0.25,0.0,0.25,+1.0]";
    Array<OneD, NekDouble> calg(1,3.0);    
    double t = 0.25;
    double v;
    double diffstep = 0.001;
    lsfitstate state;
    alglib::Lsfitcreatef(xalg,yalg,calg,diffstep,state);
    
    //alglib::spline1dinterpolant s;
    
*/    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    // curve the edges around the NEW critical layer (bool to turn on/off)
    bool curv_lay=true;
    bool move_norm=true;
    int np_lay = (nvertl-1)*npedge;//nedges*npedge (Eq. Points!!!)
    int nqedge = streak->GetExp(0)->GetNumPoints(0);
    Array<OneD, NekDouble> xcQ(nqedge*nedges,0.0);
    Array<OneD, NekDouble> ycQ(nqedge*nedges,0.0);
    Array<OneD, NekDouble> zcQ(nqedge*nedges,0.0);  
    Array<OneD, NekDouble> nxPhys(npedge*nedges,0.0);
    Array<OneD, NekDouble> nyPhys(npedge*nedges,0.0);  
    Array<OneD, NekDouble> nxQ(nqedge*nedges,0.0);
    Array<OneD, NekDouble> nyQ(nqedge*nedges,0.0);       
    if( move_norm==true)
    {
       //np_lay = (nvertl-1)*nqedge;//nedges*nqedge   (Q points!!!)
       //extract crit lay normals (through tangents):

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
/*
       for(int w=0; w< nqedge*nedges; w++)
       {
            
       }
*/
       //determine the NEW crit lay quad points values(Newtonit):
cout<<"xcQ, ycQ"<<endl;
       Computestreakpositions(nqedge*nedges, streak, xcQ, ycQ,xnull,ynull,xnull,xnull,
                       xnull,ynull, xcQ,ycQ, cr,false);
/*
for(int s=0; s<xcQ.num_elements(); s++)
{
cout<<xcQ[s]<<"     "<<ycQ[s]<<endl;
}         
ASSERTL0(false, "dsdfs"); 
*/

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
	        
            
            //determine the normal from eqs(APART FROM SIGN):
            //tx^2 +ty^2= 1 = nx^2 + ny^2;  
            //t\cdot n=0= tx*nx +ty*ny
            //result: nx = ( 1+(tx/ty)^2 )^(-1/2)
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
/*
cout<<"edge:"<<k<<endl;
cout<<"tan/normal"<<endl;
for(int r=0; r<nqedge; r++)
{
cout<<xcQ[k*nqedge+r]<<"     "<<txedgeQ[r]<<"      "<<tyedgeQ[r]<<"    "
<<nxedgeQ[r]<<"      "<<nyedgeQ[r]<<endl;
}
*/
             //force the normal at interface point to be equal
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

       int nquad_lay = (nvertl-1)*nqedge;
       Array<OneD, NekDouble>x_tmpQ(nquad_lay-(nedges-1));
       //Array<OneD, NekDouble>tmpnxQ(nquad_lay-(nedges-1));
       Array<OneD, NekDouble>tmpnyQ(nquad_lay-(nedges-1));

       Cutrepetitions(nedges, xcQ,x_tmpQ);
       //Cutrepetitions(nedges, nxQ, tmpnxQ);
       Cutrepetitions(nedges, nyQ, tmpnyQ);       
/*
for(int u=0; u<x_tmpQ.num_elements(); u++)
{
cout<<x_tmpQ[u]<<"      "<<tmpnyQ[u]<<endl;
}
*/       
       
       
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
/*       
for(int s=0; s<np_lay; s++)
{	

cout<<xcPhys[s]<<"     "<<nxPhys[s]<<"     "<<nyPhys[s]<<endl;

}        
cout<<"xcPhys,,"<<endl;
*/


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
//cout<<"Vid x="<<xnew[lay_Vids[m][h] ]<<"   y="<<ynew[lay_Vids[m][h] ]<<endl;
             if(h< nedges && curv_lay==true)
             {
//cout<<"edge=="<<h<<endl;             	     
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
                         layers_y[m][h*npedge +d+1]=  Addpointsy[h*(npedge-2) +d] +delta[m];    
                         layers_x[m][h*npedge +d+1]=  Addpointsx[h*(npedge-2) +d];
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

                         tmpy_lay[h*npedge +d+1] =  Addpointsy[h*(npedge-2) +d] +
                                  delta[m]*abs(nyPhys[h*npedge +d+1]);  
                         tmpx_lay[h*npedge +d+1]=  Addpointsx[h*(npedge-2) +d] +
                                  delta[m]*abs(nxPhys[h*npedge +d+1]);  

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
/*
for(int s=0; s<np_lay; s++)
{	

cout<<tmpx_lay[s]<<"     "<<tmpy_lay[s]<<endl;

}
//ASSERTL0(false, "dasd");
*/


                 //ASSERTL0(tmpx_lay[h*npedge +0]>=0," vert 0 x<0");
                 //ASSERTL0(tmpx_lay[h*npedge +npedge-1]>0," vert 1 x<0");


        //Orderfunctionx(tmpx_lay,tmpy_lay, Array<OneD, NekDouble>& outarray_x,
        //         Array<OneD, NekDouble>& outarray_y);         



         //check if the x coord is 'outofbound' and calculate the 
         //number of outofbound points
         
         //determine which boudn has been overcome:
         NekDouble boundleft = xcPhys[0];
         NekDouble boundright = xcPhys[np_lay-1];
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
              if(xcPhys[c*npedge+npedge-1] <= tmpx_lay[c*(npedge-(npedge-2)) +2] && outboundright==true)
              {
              	    replacepointsfromindex =   c*(npedge-(npedge-2))+2;
                    break;
              }

              
              //assume at least 1 middle point per edge         	 
              if(xcPhys[(nedges-1 -c)*npedge+0] >= tmpx_lay[np_lay-1 -(c*(npedge-(npedge-2)) +2)] 
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
	      	   increment =  (xcPhys[np_lay-outcount]-xcPhys[pstart])/(outcount+1);
                   outcount = outcount-1;  	      	   
	      	   ASSERTL0(tmpx_lay[np_lay-outcount]>xcPhys[(nedges-1)*npedge+0], "no middle points in the last edge");
	      }
	      else
	      {
                   shift=1;
                   pstart= outcount-1;
                   increment = (xcPhys[replacepointsfromindex]-xcPhys[pstart])/(outcount+1);  
                   ASSERTL0(tmpx_lay[pstart]<xcPhys[0*npedge +npedge-1], "no middle points in the first edge");                 
 	      }
	      
	      //interp to points between  posindex and posindex-1
	      Array<OneD, NekDouble> replace_x(outcount);
	      Array<OneD, NekDouble> replace_y(outcount);
	      //order normal functions(cut out repetitions)
	      Array<OneD, NekDouble>x_tmp(np_lay-(nedges-1));
	      Array<OneD, NekDouble>y_tmp(np_lay-(nedges-1));
	      Array<OneD, NekDouble>tmpny(np_lay-(nedges-1));
	      Cutrepetitions(nedges, xcPhys,x_tmp);
	      Cutrepetitions(nedges, ycPhys,y_tmp);
              Cutrepetitions(nedges, nyPhys, tmpny);   
	      //init neigh arrays
      	      Array<OneD, NekDouble>closex(4);
	      Array<OneD, NekDouble>closey(4);                     
              Array<OneD, NekDouble>closeny(4); 
                     NekDouble xctmp,ycinterp,nxinterp,nyinterp;



                     for(int v=0; v<outcount;v++)
                     {
                          xctmp = xcPhys[pstart]+(v+1)*increment;


  
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
                 }
                 
                 
                 
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


cout<<nedges<<"nedges"<<npedge<<" np_lay="<<np_lay<<endl;
         //generate x,y arrays without lastedgepoint
         //(needed to interp correctly)
         Array<OneD, NekDouble>tmpx(np_lay-(nedges-1));
         Array<OneD, NekDouble>tmpy(np_lay-(nedges-1));
/*         
         for(int b=0; b<nedges; b++)
         {
               Vmath::Vcopy( npedge-1, &tmpy_lay[b*(npedge)],1,&tmpy[b*(npedge-1)],1);
               Vmath::Vcopy( npedge-1, &tmpx_lay[b*(npedge)],1,&tmpx[b*(npedge-1)],1);
               if(b== nedges-1)
               {
                    tmpy[b*(npedge-1)+npedge-1] = tmpy_lay[b*npedge+npedge-1];
                    tmpx[b*(npedge-1)+npedge-1] = tmpx_lay[b*npedge+npedge-1];
               }
         }

         //check if there is any repetition
         for(int f=1; f< tmpx.num_elements(); f++)
         {
               string error= "repetition on arrays interp:"+boost::lexical_cast<string>(m);
               ASSERTL0(tmpx[f]!=tmpx[f-1],error);
         }
*/

         Cutrepetitions(nedges, tmpx_lay, tmpx);
         Cutrepetitions(nedges, tmpy_lay, tmpy);         


         //order points in x:
         int index;
         Array<OneD, NekDouble> copyarray_x(tmpx.num_elements());
         Array<OneD, NekDouble> copyarray_y(tmpx.num_elements());         
/*         
         //NB: it may be possible that some edges have 
         //more/less points than npedge
         //@todo: write an orderfunction

         Vmath::Vcopy(tmpx.num_elements(), tmpx,1, copyarray_x,1);
         Vmath::Vcopy(tmpx.num_elements(), tmpy,1, copyarray_y,1);
         NekDouble max = Vmath::Vmax(tmpx.num_elements(), tmpx,1);
         for(int w=0; w<tmpx.num_elements(); w++)
         {
            index = Vmath::Imin(tmpx.num_elements(), copyarray_x,1);
            tmpx[w]= copyarray_x[index];
            tmpy[w]= copyarray_y[index];
            copyarray_x[index] = max+1000;            
         }
*/


         Orderfunctionx(tmpx, tmpy, tmpx, tmpy);        	
/*        	
cout<<"ordered"<<endl;
for(int g=0; g<tmpx.num_elements(); g++)
{
cout<<tmpx[g]<<"    "<<tmpy[g]<<endl;
}
*/



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
//cout<<"determine value for x="<<Addpointsx[g*(npedge-2) +r]<<endl;
                 //determine closest point index:
                 index = 
                 DetermineclosePointxindex( Addpointsx[g*(npedge-2) +r], tmpx);
//cout<<"  vert+"<<index<<endl;

                 ASSERTL0( index<= tmpy.num_elements()-1, " index wrong");
                 //generate neighbour arrays Pyinterp,Pxinterp
                 GenerateNeighbourArrays(index, closepoints,tmpx,tmpy,Pxinterp,Pyinterp);
                 
                 layers_y[m][g*npedge +r+1]=
                                 LagrangeInterpolant(
                          Addpointsx[g*(npedge-2) +r],closepoints,Pxinterp,Pyinterp  );
//cout<<"x value="<<Addpointsx[g*(npedge-2) +r]<<endl;
                 layers_x[m][g*npedge +r+1]=  Addpointsx[g*(npedge-2) +r];




/*
for(int t=0; t<6; t++)
{
cout<<"Px="<<Pxinterp[t]<<"    "<<Pyinterp[t]<<endl;
}
*/
             }


             }//if edge closed g
         }
         //check if there are points out of range:
         ASSERTL0(Vmath::Vmax(np_lay,layers_y[m],1)< Vmath::Vmax(nVertTot,yold,1),"point>ymax");
         ASSERTL0(Vmath::Vmin(np_lay,layers_y[m],1)> Vmath::Vmin(nVertTot,yold,1),"point<ymin");

   
/*
cout<<" xlay    ylay"<<endl;
for(int l=0; l<np_lay; l++)
{
//cout<<tmpx_lay[l]<<"    "<<tmpy_lay[l]<<endl;
cout<<std::setprecision(8)<<layers_x[m][l]<<"    "<<layers_y[m][l]<<endl;
}

cout<<"nverts"<<endl;
for(int l=0; l<nvertl; l++)
{
cout<<std::setprecision(8)<<xnew[lay_Vids[m][l] ]<<"    "<<ynew[lay_Vids[m][l] ]<<endl;
}
*/

//ASSERTL0(false, "as");

     }//close layers!!! m index



     //update vertices coords outside layers region
     NekDouble ratio;
     for(int n=0; n<nVertTot; n++)
     {
          NekDouble ratio;  
          SpatialDomains::VertexComponentSharedPtr vertex = graphShPt->GetVertex(n);
          NekDouble x,y,z;
          vertex->GetCoords(x,y,z); 
          int qp_closer,diff;
          //determine the closer xold_up
          NekDouble tmp=1000;

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
     
     }//move_norm bool
     else//move vertically
     {
          MoveLayersvertically(nlays, nvertl, cntlow, cntup,
    	         lay_Vids, x_c, y_c, Down, Up, xnew, ynew, layers_x, layers_y);     	          	     
   	     
     }
     

     //check borders of the new mesh verts:
//cout<<std::setprecision(8)<<"yoldmax="<<Vmath::Vmax(nVertTot, yold,1)<<endl;
//cout<<std::setprecision(8)<<"ynewmax="<<Vmath::Vmax(nVertTot,ynew,1)<<endl;
     ASSERTL0(Vmath::Vmin(nVertTot, xold,1)== Vmath::Vmin(nVertTot,xnew,1),
             "  different xmin val");
     ASSERTL0(Vmath::Vmin(nVertTot, yold,1)== Vmath::Vmin(nVertTot,ynew,1),
             "  different ymin val");
     ASSERTL0(Vmath::Vmax(nVertTot, xold,1)== Vmath::Vmax(nVertTot,xnew,1),
             "  different xmax val");
     ASSERTL0(Vmath::Vmax(nVertTot, yold,1)== Vmath::Vmax(nVertTot,ynew,1),
             "  different ymax val");



     //check x,y coords
/*
     NekDouble tmpdiff;
     for(int f=0; f< nlays; f++)
     {
//cout<<"delta lay="<<f<<"  is:"<<delta[f]<<endl;
         for(int u=0; u< nvertl; u++)
         {          
               if(u==0) 
               {
//cout<<"vert x="<<xnew[lay_Vids[f][u] ]<<"  curv x="<<layers_x[f][0]<<endl;
                  ASSERTL0(xnew[lay_Vids[f][u] ]==layers_x[f][0],"wro");
//cout<<"vert y="<<ynew[lay_Vids[f][u] ]<<"  curv y="<<layers_y[f][0]<<endl;
                  ASSERTL0(ynew[lay_Vids[f][u] ]==layers_y[f][0],"wroy");
               }
               else if(u== nvertl-1)
               {
//cout<<"vert x="<<xnew[lay_Vids[f][u] ]<<"  curv x="<<layers_x[f][(nedges-1)*npedge +npedge-1]<<endl;
                    ASSERTL0(
                       xnew[lay_Vids[f][u] ]==layers_x[f][(nedges-1)*npedge +npedge-1],"wroLAST");


//cout<<"vert y="<<ynew[lay_Vids[f][u] ]<<"  curv y="<<layers_y[f][(nedges-1)*npedge +npedge-1]<<endl;
                    ASSERTL0(ynew[lay_Vids[f][u] ]==layers_y[f][(nedges-1)*npedge +npedge-1],"wroyLAST");
               }
               else
               {
//cout<<"vert x="<<xnew[lay_Vids[f][u] ]<<"  curv x1="<<layers_x[f][u*npedge +0]<<" curv x2="
//<<layers_x[f][(u-1)*npedge +npedge-1]<<endl;
                    ASSERTL0(xnew[lay_Vids[f][u] ]== layers_x[f][u*npedge +0],"wrongg");
                    ASSERTL0(layers_x[f][u*npedge +0]== layers_x[f][(u-1)*npedge +npedge-1],"wrong1");

//cout<<"vert y="<<ynew[lay_Vids[f][u] ]<<"  curv y1="<<layers_y[f][u*npedge +0]<<" curv y2="
//<<layers_y[f][(u-1)*npedge +npedge-1]<<endl;
                    ASSERTL0(ynew[lay_Vids[f][u] ]== layers_y[f][u*npedge +0],"wronggy");
                    ASSERTL0(layers_y[f][u*npedge +0]== layers_y[f][(u-1)*npedge +npedge-1],"wrong1y");
               }
          }
      }
*/


    //------------------------------------------------------------

    //check if the layers need to be shifted
/*
    string movelay;
    int Vid_di;
    NekDouble delt, delt_opp;
    NekDouble tmpleft=10;
    NekDouble tmpright = -10;

    //check top left (x=0, y>yup)
    //last layer nlays-1
    if(
         (  abs(ynew[topleft]-ynew[ lay_Vids[nlays-1][0] ])<0.01)||
         (  ynew[topleft]< ynew[lay_Vids[nlays-1][0] ])
      )
    {

        movelay = "layup";
        //calculate the new delta based on the space at the border:
        delt = abs(y_c[0]-ynew[topleft])/(nlays/2);
        delt_opp = abs(y_c[0]-ynew[bottomleft])/(nlays/2);
cout<<"move layerup to down:"<<endl;
    }
    //check bottom right (x=pi/2, y<ydown)
    //first lay: 0
    if(
        (  abs(ynew[bottomright]-ynew[ lay_Vids[0][(nvertl-1)] ])<0.01 ) ||
        (  ynew[bottomright]> ynew[ lay_Vids[0][(nvertl-1)] ]  )
      )
      {

         movelay= "laydown";
         //calculate the new delta based on the space at the border:
         delt = abs(y_c[nvertl-1]-ynew[bottomright])/(nlays/2);
         delt_opp = abs(y_c[nvertl-1]-ynew[topright])/(nlays/2);
cout<<"move layerdown to up:"<<endl;
            
      }
        
    
    //ChangeLayerspos(ynew, y_c, Addpointsy, layers_y, lay_Vids, delt, delt_opp, movelay, npedge);
*/
    //--------------------------------------------------------------


    
    //replace the vertices with the new ones

    Replacevertices(changefile, xnew , ynew, x_c, y_c, Addpointsx, Addpointsy, Eids, npedge, charalp, layers_x,layers_y, lay_eids, curv_lay);
          	       
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

	void Computestreakpositions(int npoints, MultiRegions::ExpListSharedPtr &streak,
    	        Array<OneD, NekDouble> x,  Array<OneD, NekDouble> y,
    	        Array<OneD, NekDouble> &xold_up, Array<OneD, NekDouble> &yold_up,
    	        Array<OneD, NekDouble> &xold_low, Array<OneD, NekDouble> &yold_low,
    	        Array<OneD, NekDouble> &xold_c, Array<OneD, NekDouble> &yold_c,    	        
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
/* 
//we can't rely on the layer down/up position to calc the streak position !!!(criti lay 
//EXTERNAL to the layer zone)              	     
             for(int q=0; q<npoints; q++)
             {
                  NekDouble streaktmp =200;
                  NekDouble streaktmppos=200;
                  NekDouble streaktmpneg=-200;
                  int ipos,ineg, it =0;                  
                  NekDouble weightpos, weightneg;
  
cout<<"nq="<<nq<<" xold_down="<<xold_low[q]<<"  xold_up="<<xold_up[q]<<"  xold_c="<<xold_c[q]<<endl;    

                    ASSERTL0(xold_low[q]==xold_up[q], "layer region non valid");
cout<<q<<"   xup="<<xold_up[q]<<"   yup="<<yold_up[q]<<"    ydown="<<yold_low[q]<<endl;                 
                    //algorithm to determine the closest points to the streak positions
                    //using the weighted mean         
                    for(int j=0; j<nq; j++)
                    {
                    	 if(x[j]==xold_c[q] && y[j]<= yold_up[q] && y[j]>= yold_low[q])
                    	 {
                             if(streak->GetPhys()[j]< streaktmppos && streak->GetPhys()[j]>cr)
                             {	                 	      
                                 streaktmppos = streak->GetPhys()[j];
                                 ipos =j;
                             }
                             static int cnt=0;

                             if( abs(streak->GetPhys()[j]) < -streaktmpneg && streak->GetPhys()[j]<cr)
                             {
//cout<<"test streakneg="<<streak->GetPhys()[j]<<endl;
              	                 streaktmpneg = streak->GetPhys()[j];
              	                 ineg =j;
              	             }
              
//cout<<" x="<<x[j]<<"  y="<<y[j]<<"   streak="<<streak->GetPhys()[j]<<endl;
			 }
		    }
cout<<"ipos="<<ipos<<"  ineg="<<ineg<<endl;
//cout<<"closer streak points ypos="<<y[ipos]<<"   yneg="<<y[ineg]<<endl;
	          //determine the streak y position as the result of the weighted mean
                    if(streaktmppos< 200 && streaktmpneg>-200)
                    {
   	                 xc[q]= x[ipos];
	                 weightpos = 1/(streaktmppos*streaktmppos);
	                 weightneg = 1/(streaktmpneg*streaktmpneg);	            
	                 yc[q]= ( (y[ipos]*weightpos) + (y[ineg]*weightneg) )/(weightpos+weightneg);
                    }
                    else if(streaktmppos< 200)
                    {
   	                 xc[q]= x[ipos];
   	                 yc[q]= y[ipos];
                    }
                    else if(streaktmpneg> -200)
                    {
   	                 xc[q]= x[ineg];
   	                 yc[q]= y[ineg];
                    }
                    else
                    {
                         ASSERTL0(false, "streak not found");
                    }

cout<<" streak x="<<xc[q]<<"   y="<<yc[q]<<endl;
              }
*/
                   //start guess
                   //yc= (yup+ydown)/2
                   Vmath::Vadd(xc.num_elements(), yold_up,1,yold_low,1, yc,1);
                   Vmath::Smul(xc.num_elements(), 0.5,yc,1,yc,1);
                   Vmath::Vcopy(xc.num_elements(),xold_c,1,xc,1);
              }
              else//case of xQ,yQ
              {
                   Vmath::Vcopy(xc.num_elements(), x,1,xc,1);
                   Vmath::Vcopy(xc.num_elements(), y,1,yc,1);
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
//cout<<"start guess:  x="<<xc[e]<<"    y="<<yc[e]<<endl;
		   while( abs(F)> 0.000000001)
		   {
                	ytmp = coord[1];
              	        elmtid = streak->GetExpIndex(coord,0.00001);
           	        offset = streak->GetPhys_Offset(elmtid);
		   	U = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() + offset);
		   	dU  = streak->GetExp(elmtid)->PhysEvaluate(coord, derstreak + offset);
		   	coord[1] = coord[1] - (U-cr)/dU;   
		   	F = U-cr;   
		   	ASSERTL0( coord[0]==xc[e], " x coordinate must remain the same");
              	        //stvalues[e] = streak->GetExp(elmtid)->PhysEvaluate(coord, streak->GetPhys() +offset );
//cout<<"elmtid="<<elmtid<<"  x="<<coord[0]<<"   y="<<coord[1]<<"    stvalue="<<U<<"   dU="<<dU<<endl;
                        if(abs(coord[1])>1 
                           && attempt==0  && verts==true
                          )
                        {
                             
                             //try the old crit lay position:
                             coord[1] = yold_c[e];
                             attempt++;
                        }
                        else if(abs(coord[1])>1 )
                        {
                             coord[1] = ytmp +0.01;
                             elmtid = streak->GetExpIndex(coord,0.00001);
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

                        its++;
                        if(its>1000 && F< 0.0001)
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


	void GenerateAddPointsNewtonIt( NekDouble &xi, NekDouble &yi,NekDouble &x0, NekDouble &y0,
    	        MultiRegions::ExpListSharedPtr &function, Array<OneD, NekDouble> &derfunction,
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
		while( abs(F)> 0.00000001)
	        {
	             ytmp = coords[1];
//cout<<"generate newton it xi="<<xi<<"  yi="<<yi<<endl;			
		     elmtid = function->GetExpIndex(coords, 0.00001);
                //@to do if GetType(elmtid)==triangular WRONG!!!
//cout<<"gen newton xi="<<xi<<"  yi="<<coords[1]<<"  elmtid="<<elmtid<<"  F="<<F<<endl;			
		     offset = function->GetPhys_Offset(elmtid);

		     U = function->GetExp(elmtid)->PhysEvaluate(coords, function->GetPhys() + offset);
		     dU  = function->GetExp(elmtid)->PhysEvaluate(coords, derfunction + offset);
		     coords[1] = coords[1] - (U-cr)/dU;   
//cout<<cr<<"U-cr="<<U-cr<<"  tmp result y:"<<coords[1]<<"  dU="<<dU<<endl;
		     F = U-cr;   

                     if(  abs(coords[1])>1 
                           //&& attempt==0 
                       )
                     {

                          coords[1] = ytmp +0.01;
                          elmtid = function->GetExpIndex(coords,0.00001);
             	          offset = function->GetPhys_Offset(elmtid);
		          NekDouble Utmp = function->GetExp(elmtid)->PhysEvaluate(coords, function->GetPhys() + offset);
                          NekDouble dUtmp = function->GetExp(elmtid)->PhysEvaluate(coords, derfunction + offset);
   		   	  coords[1] = coords[1] - (Utmp-cr)/dUtmp;
cout<<"attempt:"<<coords[1]<<endl;
                          if( (abs(Utmp-cr)>abs(F))||(abs(coords[1])>1)  )
                          {
                               coords[1] = ytmp -0.01;
                          }
                             
                          attempt++;
                     }
                     else
                     {
                          ASSERTL0(abs(coords[1])<= 1, " y value out of bound +/-1");
                     }

                     its++;
                     if(its>1000 && F< 0.0001)
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
	        x0 = xi;
	        y0 = coords[1] - (U-cr)/dU;
cout<<"NewtonIt result  x="<<x0<<"  y="<<coords[1]<<"   U="<<U<<endl;	        
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


        void GenerateCurveSp(Array<OneD, NekDouble> &x_c, Array<OneD, NekDouble> &y_c,
                Array<OneD, NekDouble> &x_lay, Array<OneD, NekDouble> &y_lay)
        {
               Array<OneD, NekDouble> k(x_c.num_elements());
               Array<OneD, NekDouble> b(x_c.num_elements());
               NekDouble dx,dx1,dy,dy1;
               //determine the values of k
               //remark fix the derivative=0 at the border of the domain: k_0=k_n=0
               for(int s=0; s< k.num_elements(); s++)
               {
                   dx = x_c[s] - x_c[s-1];
                   dy = y_c[s] - y_c[s-1];
                   dx1 = x_c[s+1] - x_c[s];
                   dy1 = y_c[s+1] - y_c[s];
                   if(s==0 || s==k.num_elements()-1)
                   {
                        k[s]=0;
                        
                   }
                   else
                   {
                       b[s] = 3*(   dy/(dx*dx) +dy1/(dx1*dx1)  );
                   }
               }
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

        void MapEdgeVertices(MultiRegions::ExpListSharedPtr field, Array<OneD, int> &V1,
                 Array<OneD, int> &V2)
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

        void Findlay_eids(Array<OneD, Array<OneD, NekDouble> >& layers_y,
                 Array<OneD, int> V1, Array<OneD, int> V2, 
                 Array<OneD, Array<OneD, int > >& lay_Vids,
                 Array<OneD, Array<OneD, int > >& lay_eids)
        {
             int neids = V1.num_elements();
             int nlays = layers_y.num_elements();
             int nvertl = lay_Vids[0].num_elements();
             Array<OneD, NekDouble > tmpy(nlays,-1000);
             Array<OneD, int > tmpVids(nlays,-1000);
             Array<OneD, int > tmpIndex(nlays,-1000);

             //sort layers_y low to high
             for(int a=0; a<nvertl; a++)
             {
                  tmpy = Array<OneD, NekDouble> (nlays,-1000);
                  //copy into array:(not possible use Vmath with matrices)
                  for(int t=0; t<nlays; t++)
                  {
                       tmpy[t] = layers_y[t][a];
                       tmpVids[t] = lay_Vids[t][a];                       
                  }
                  

                  for(int b=0; b<nlays; b++)
                  {
                       tmpIndex[b] = Vmath::Imin(nlays,tmpy,1);
                       layers_y[b][a] = tmpy[tmpIndex[b]];
                       lay_Vids[b][a] = tmpVids[tmpIndex[b]];  
                       //tmpy[b] =  Vmath::Vmin(nlays,&(layers_y[0][a]),1);
                       ASSERTL0(layers_y[b][a]==Vmath::Vmin(nlays,tmpy,1),
                       "sort failed ");
                       tmpy[tmpIndex[b]] =  1000;
                     
                  }

             }

             for(int a=0; a< nvertl-1; a++)
             {
                  for(int b=0; b<nlays; b++)
                  {
                      for(int p=0; p< neids; p++)
                      {
                          if( 
                              (lay_Vids[b][a] == V1[p] && lay_Vids[b][a+1]==V2[p])||
                              (lay_Vids[b][a] == V2[p] && lay_Vids[b][a+1]==V1[p])
                            )
                            {
                                lay_eids[b][a] = p;  
                                break;
                            }
                      }
                  } 
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
        void  Cutrepetitions(int nedges,Array<OneD, NekDouble> inarray,
                 Array<OneD, NekDouble>& outarray)
        {
        	
             //determine npedge:
             int np_lay = inarray.num_elements();
           
             int npedge = np_lay/nedges;
             ASSERTL0(inarray.num_elements()%nedges==0," something on number npedge");
             //cut out the repetitions:

           cout<<nedges<<"nedges"<<npedge<<" nplay="<<np_lay<<endl;
             //generate x,y arrays without lastedgepoint
             //(needed to interp correctly)

             for(int b=0; b<nedges; b++)
             {
                  Vmath::Vcopy( npedge-1, &inarray[b*(npedge)],1,&outarray[b*(npedge-1)],1);
                  if(b== nedges-1)
                  {
                    outarray[b*(npedge-1)+npedge-1] = inarray[b*npedge+npedge-1];

                  }
             }        	
        	
        	
        	
        }

        void  Orderfunctionx(Array<OneD, NekDouble> inarray_x,
                 Array<OneD, NekDouble> inarray_y, Array<OneD, NekDouble>& outarray_x,
                 Array<OneD, NekDouble>& outarray_y)
        {

             Array<OneD, NekDouble>tmpx(inarray_x.num_elements());
             //local copy to prevent overwriting
             Vmath::Vcopy(inarray_x.num_elements() , inarray_x,1,tmpx,1);

             //order function with respect to x
             int index;

             NekDouble max = Vmath::Vmax(tmpx.num_elements(), tmpx,1);
             for(int w=0; w<tmpx.num_elements(); w++)
             {
                 index = Vmath::Imin(tmpx.num_elements(), tmpx,1);
                 outarray_x[w]= tmpx[index];
                 outarray_y[w]= inarray_y[index];
                 tmpx[index] = max+1000;            
             }

             
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
//cout<<"case0"<<endl;
                  diff = index-leftpoints; 
                  start= 0;
                  Vmath::Vcopy(neighpoints, &yArray[0],1,&Neighbour_y[0],1);
                  Vmath::Vcopy(neighpoints, &xArray[0],1,&Neighbour_x[0],1);
             }
             else if( (yArray.num_elements()-1)-index < rightpoints)
             {
//cout<<"case1 closest="<<xArray[index]<<endl;
                  int rpoints = (yArray.num_elements()-1)-index;//
                  diff = rightpoints-rpoints;
//cout<<"start index="<<index-leftpoints-diff<<endl;
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
        
        void ChangeLayerspos(Array<OneD, NekDouble> & ynew,
                 Array<OneD, NekDouble> yc,
                 Array<OneD, NekDouble> Addpointsy,
                 Array<OneD, Array<OneD, NekDouble > >& layers_y,
                 Array<OneD, Array<OneD, int> >lay_Vids,
                 NekDouble delt, NekDouble delt_opp,string movelay, int npedge)
        {
             int nvertl = yc.num_elements();
             int nedges = nvertl-1;
             int nlays = layers_y.num_elements();
             int np_lay = layers_y[0].num_elements();
             Array<OneD, Array<OneD, NekDouble> > tmpy(layers_y.num_elements());
             Array<OneD, Array<OneD, NekDouble> > tmpynew(layers_y.num_elements());
             int cntup=1;
             NekDouble delta;
             //cp layers_y  into tmpy and move the first one     
cout<<"movelay"<<movelay<<endl;        
             if(movelay == "laydown")
             {             
                 delta = abs(layers_y[nlays-1][(nvertl-2)*npedge+npedge-1]-yc[nvertl-1]);

                 for(int g=0; g<  nlays;  g++)
                 {
                      tmpy[g] = Array<OneD, NekDouble> (np_lay);
                      tmpynew[g] = Array<OneD, NekDouble>(nvertl);
                      for(int h=0; h< nvertl; h++)
                      {

                          //move layers_up
                          if(g==0 )
                          {  
                               tmpynew[g][h]= yc[h] + delta + delta/(nlays/2);

                               if(h<nedges)
                               {

                                  //v1
                                  tmpy[g][h*npedge +0] = yc[h] + delta+ delta/(nlays/2);
                                  //v2
                                  tmpy[g][h*npedge +npedge-1] = yc[h+1] + delta+delta/(nlays/2);
                                  //middle points (shift crit lay points by delta):
                                  for(int d=0; d< npedge-2; d++)
                                  {
                                        tmpy[g][h*npedge +d+1]=  
                                        Addpointsy[h*(npedge-2) +d] +delta + delta/(nlays/2);    

                                  }
                              }                             
                           }
                           else
                           {
                               tmpynew[g][h]= ynew[ lay_Vids[g][h] ];
                               if(h<nedges)
                               {
                                  //v1
                                  tmpy[g][h*npedge +0] = layers_y[g][h*npedge +0];
                                  //v2
                                  tmpy[g][h*npedge +npedge-1] = layers_y[g][h*npedge +npedge-1];
                                  //middle points (shift crit lay points by delta):
                                  for(int d=0; d< npedge-2; d++)
                                  {
                                       tmpy[g][h*npedge +d+1]=  layers_y[g][h*npedge +d+1];
                                  }                                
                               } 
                           }                              

                      }
                 }

                 //copy back in a different order
                 int cnt=-1;
                 int vcnt=-1;
                 for(int g=0; g<  nlays;  g++)
                 {
cout<<"move y coords"<<g<<endl;
                      for(int h=0; h<np_lay; h++)
                      {
                          //first becomes last
                          if(g==0)
                          {
                               layers_y[nlays-1][h]= tmpy[0][h];
                          }
                          else
                          { 
                               layers_y[cnt][h] = tmpy[g][h];
                          }                          
                      }
                      cnt++;

                      for(int v=0; v<nvertl; v++)
                      {
                          if(g==0)
                          {
                               ynew[lay_Vids[nlays-1][v]] = tmpynew[0][v];
                          }
                          else
                          {
                               ynew[lay_Vids[vcnt][v] ] = tmpynew[g][v];          
                          }
                      }
                      vcnt++;


                      
                  }
             }
             if( movelay =="layup")
             {
                  ASSERTL0(false, "not implemented layup");
             }
                


             
	}      
    
        void Replacevertices(string filename, Array<OneD, NekDouble> newx, 
    	               Array<OneD,  NekDouble> newy,
    	               Array<OneD, NekDouble> x_crit, Array<OneD, NekDouble> y_crit,
    	               Array<OneD, NekDouble> & Pcurvx, 
    	               Array<OneD, NekDouble> & Pcurvy,
    	               Array<OneD, int>&Eids, int Npoints, string s_alp,
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
   	    	        if(x_crit[cnt]< x_crit[cnt+1])
   	    	        {
   	    	            //x1 = x_crit[cnt];
   	    	            //x2 = x_crit[cnt+1];
   	    	            //y1 = y_crit[cnt];
   	    	            //y2 = y_crit[cnt+1];
   	    	            v1 = cnt;
   	    	            v2 = cnt+1;
   	    	        }
   	    	        else
   	    	        {
   	    	            //x1 = x_crit[cnt+1];
   	    	            //x2 = x_crit[cnt];
   	    	            //y1 = y_crit[cnt+1];
   	    	            //y2 = y_crit[cnt];
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
		        st << std::scientific << std::setprecision(8) << x_crit[v1] << "   "
	    	        << y_crit[v1] << "   " << 0.000<<"   ";
//cout<<x_crit[v1]<<"       "<<y_crit[v1]<<endl;
		        for(int a=0; a< Npoints-2; a++)
                        {		       	   
		             st << std::scientific << std::setprecision(8) <<
		             "    "<<Pcurvx[indexeid*(Npoints-2) +a]<<"    "<<Pcurvy[indexeid*(Npoints-2) +a]
		             <<"    "<<0.000<<"   "; 
//cout<<Pcurvx[indexeid*(Npoints-2) +a]<<"        "<<Pcurvy[indexeid*(Npoints-2)+a]<<endl;	      
		        }
		        st << std::scientific << std::setprecision(8) <<
		        "    "<<x_crit[v2]<<"   "<< y_crit[v2] <<"   "<< 0.000;
//cout<<x_crit[v2]<<"       "<<y_crit[v2]<<endl;
                        edgenew->LinkEndChild(new TiXmlText(st.str()));    		   

//cout<<st.str()<<endl;		
  
                        //if(x_crit[cnt] < x_crit[cnt+1])
                        //{
                        //   cout<<"Warning: the edges can have the vertices in the reverse order";
                        //}
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
               int neids_lay = lay_eids[0].num_elements();
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
	
