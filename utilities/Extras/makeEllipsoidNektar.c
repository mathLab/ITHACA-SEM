#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#define PI 3.14159265358979323846

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Dr. Laurent Risser (Imperial College London)
//l.risser@imperial.ac.uk
//25.08.2009
//
//Command line to compile:  gcc -o makeEllipsoidNektar makeEllipsoidNektar.c -O -lm
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// ---------------------------------------------------------------------------------------
//                               structures definition
// ---------------------------------------------------------------------------------------

typedef struct
{
double X;
double Y;     // coordinates of a vertex
double Z;
} Vertex;

typedef struct
{
int V1;   //id of the first vertex
int V2;   //id of the second vertex
Vertex * SubiVertex;
int NbSubdiv;
} Edge;

typedef struct
{
int E1;    //id of the first edge
int E2;    //id of the second edge
int E3;    //id of the third edge
} Element;

typedef struct
{
int NbVertexes;
Vertex * Vertexes;
int NbEdges;
Edge * Edges;
int NbElements;
Element * Elements;
} Mesh;




// ---------------------------------------------------------------------------------------
//                                      functions
// ---------------------------------------------------------------------------------------


//create a simple cubic mesh 
Mesh CreateBasicCubicMesh(){
	Mesh LocMesh;
	int i;
	
	
	//allocations...
	LocMesh.NbVertexes=8;
	LocMesh.Vertexes=(Vertex *)malloc(LocMesh.NbVertexes*sizeof(Vertex));
	
	LocMesh.NbEdges=18;
	LocMesh.Edges=(Edge *)malloc(LocMesh.NbEdges*sizeof(Edge));
        
	LocMesh.NbElements=12;
	LocMesh.Elements=(Element *)malloc(LocMesh.NbElements*sizeof(Element));
	
	//fill the coordinates of each point
	LocMesh.Vertexes[0].X=0; LocMesh.Vertexes[0].Y=0; LocMesh.Vertexes[0].Z=0;
	LocMesh.Vertexes[1].X=1; LocMesh.Vertexes[1].Y=0; LocMesh.Vertexes[1].Z=0;
	LocMesh.Vertexes[2].X=1; LocMesh.Vertexes[2].Y=1; LocMesh.Vertexes[2].Z=0;
	LocMesh.Vertexes[3].X=0; LocMesh.Vertexes[3].Y=1; LocMesh.Vertexes[3].Z=0;
	LocMesh.Vertexes[4].X=0; LocMesh.Vertexes[4].Y=0; LocMesh.Vertexes[4].Z=1;
	LocMesh.Vertexes[5].X=1; LocMesh.Vertexes[5].Y=0; LocMesh.Vertexes[5].Z=1;
	LocMesh.Vertexes[6].X=1; LocMesh.Vertexes[6].Y=1; LocMesh.Vertexes[6].Z=1;
	LocMesh.Vertexes[7].X=0; LocMesh.Vertexes[7].Y=1; LocMesh.Vertexes[7].Z=1;
	
	//fill each edge
	LocMesh.Edges[0].V1=0;  LocMesh.Edges[0].V2=1;
	LocMesh.Edges[1].V1=1;  LocMesh.Edges[1].V2=2;
	LocMesh.Edges[2].V1=2;  LocMesh.Edges[2].V2=3;
	LocMesh.Edges[3].V1=3;  LocMesh.Edges[3].V2=0;
	LocMesh.Edges[4].V1=4;  LocMesh.Edges[4].V2=5;
	LocMesh.Edges[5].V1=5;  LocMesh.Edges[5].V2=6;
	LocMesh.Edges[6].V1=6;  LocMesh.Edges[6].V2=7;
	LocMesh.Edges[7].V1=7;  LocMesh.Edges[7].V2=4;
	LocMesh.Edges[8].V1=0;  LocMesh.Edges[8].V2=4;
	LocMesh.Edges[9].V1=1;  LocMesh.Edges[9].V2=5;
	LocMesh.Edges[10].V1=2;  LocMesh.Edges[10].V2=6;
	LocMesh.Edges[11].V1=3;  LocMesh.Edges[11].V2=7;
	LocMesh.Edges[12].V1=1;  LocMesh.Edges[12].V2=4;
	LocMesh.Edges[13].V1=1;  LocMesh.Edges[13].V2=3;
	LocMesh.Edges[14].V1=1;  LocMesh.Edges[14].V2=6;
	LocMesh.Edges[15].V1=7;  LocMesh.Edges[15].V2=0;
	LocMesh.Edges[16].V1=7;  LocMesh.Edges[16].V2=2;
	LocMesh.Edges[17].V1=7;  LocMesh.Edges[17].V2=5;
        
        for (i=0;i<LocMesh.NbEdges;i++) LocMesh.Edges[i].NbSubdiv=1;
	
	//fill each Element
	LocMesh.Elements[0].E1=0;  LocMesh.Elements[0].E2=8;  LocMesh.Elements[0].E3=12;
	LocMesh.Elements[1].E1=9;  LocMesh.Elements[1].E2=4;  LocMesh.Elements[1].E3=12;
	LocMesh.Elements[2].E1=0;  LocMesh.Elements[2].E2=3;  LocMesh.Elements[2].E3=13;
	LocMesh.Elements[3].E1=1;  LocMesh.Elements[3].E2=2;  LocMesh.Elements[3].E3=13;
	LocMesh.Elements[4].E1=1;  LocMesh.Elements[4].E2=10;  LocMesh.Elements[4].E3=14;
	LocMesh.Elements[5].E1=9;  LocMesh.Elements[5].E2=5;  LocMesh.Elements[5].E3=14;
	LocMesh.Elements[6].E1=3;  LocMesh.Elements[6].E2=11;  LocMesh.Elements[6].E3=15;
	LocMesh.Elements[7].E1=8;  LocMesh.Elements[7].E2=7;  LocMesh.Elements[7].E3=15;
	LocMesh.Elements[8].E1=2;  LocMesh.Elements[8].E2=11;  LocMesh.Elements[8].E3=16;
	LocMesh.Elements[9].E1=10;  LocMesh.Elements[9].E2=6;  LocMesh.Elements[9].E3=16;
	LocMesh.Elements[10].E1=4;  LocMesh.Elements[10].E2=7;  LocMesh.Elements[10].E3=17;
	LocMesh.Elements[11].E1=5;  LocMesh.Elements[11].E2=6;  LocMesh.Elements[11].E3=17;
	
        return LocMesh;
}


//subdivide into 2 edges the edge 'IdEdge' in the mesh 'LocMesh'
void CutEdge(Mesh * LocMesh,int IdEdge){
  int i;
  int IdNgbhVortex1;
  int IdNgbhVortex2;
  int IdNgbhVortex3;
  int IdNgbhVortex4;
  int IdNewVortex;
  int IdSubdiviedEdge;
  int IdNewEdge1;
  int IdNewEdge2;
  int IdNewEdge3;
  int IdNgbhEdge1;
  int IdNgbhEdge2;
  int IdNgbhEdge3;
  int IdNgbhEdge4;
  int IdSubdiviedElement1;
  int IdSubdiviedElement2;
  int IdNewElement2;
  int IdNewElement1;
  
  //reallocations in LocMesh
  LocMesh->NbVertexes++;
  LocMesh->NbEdges=LocMesh->NbEdges+3;
  LocMesh->NbElements=LocMesh->NbElements+2;
  
  LocMesh->Vertexes=(Vertex *)realloc(LocMesh->Vertexes,LocMesh->NbVertexes*sizeof(Vertex));
  LocMesh->Edges=(Edge *)realloc(LocMesh->Edges,LocMesh->NbEdges*sizeof(Edge));
  LocMesh->Elements=(Element *)realloc(LocMesh->Elements,LocMesh->NbElements*sizeof(Element));
  
  //identifiers to take into account
  IdNgbhVortex1=LocMesh->Edges[IdEdge].V1;
  IdNgbhVortex2=LocMesh->Edges[IdEdge].V2;
  IdNewVortex=LocMesh->NbVertexes-1;
  
  IdSubdiviedEdge=IdEdge;
  IdNewEdge1=LocMesh->NbEdges-3;
  IdNewEdge2=LocMesh->NbEdges-2;
  IdNewEdge3=LocMesh->NbEdges-1;
  
  IdSubdiviedElement1=-1;
  IdSubdiviedElement2=-1;
  for (i=0;i<LocMesh->NbElements-2;i++) if ((LocMesh->Elements[i].E1==IdEdge)||(LocMesh->Elements[i].E2==IdEdge)||(LocMesh->Elements[i].E3==IdEdge)){
    if (IdSubdiviedElement1==-1){
      IdSubdiviedElement1=i;
      IdNewElement1=LocMesh->NbElements-2;
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex1))){
        IdNgbhEdge1=LocMesh->Elements[i].E1;
        if (LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1) IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E1].V2;
        else IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E1].V1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex1))){
        IdNgbhEdge1=LocMesh->Elements[i].E2;
        if (LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1) IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E2].V2;
        else IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E2].V1;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex1))){
        IdNgbhEdge1=LocMesh->Elements[i].E3;
        if (LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1) IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E3].V2;
        else IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E3].V1;
      }
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex2))){
        IdNgbhEdge3=LocMesh->Elements[i].E1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex2))){
        IdNgbhEdge3=LocMesh->Elements[i].E2;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex2))){
        IdNgbhEdge3=LocMesh->Elements[i].E3;
      }
    }
    else{
      IdSubdiviedElement2=i;
      IdNewElement2=LocMesh->NbElements-1;
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex1))){
        IdNgbhEdge2=LocMesh->Elements[i].E1;
        if (LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1) IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E1].V2;
        else IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E1].V1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex1))){
        IdNgbhEdge2=LocMesh->Elements[i].E2;
        if (LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1) IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E2].V2;
        else IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E2].V1;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex1))){
        IdNgbhEdge2=LocMesh->Elements[i].E3;
        if (LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1) IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E3].V2;
        else IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E3].V1;
      }
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex2))){
        IdNgbhEdge4=LocMesh->Elements[i].E1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex2))){
        IdNgbhEdge4=LocMesh->Elements[i].E2;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex2))){
        IdNgbhEdge4=LocMesh->Elements[i].E3;
      }
    }
  }
  IdNewElement1=LocMesh->NbElements-2;
  IdNewElement2=LocMesh->NbElements-1;
  
  //Coordinates of the new vertex
  LocMesh->Vertexes[IdNewVortex].X=(LocMesh->Vertexes[IdNgbhVortex1].X+LocMesh->Vertexes[IdNgbhVortex2].X)/2.;
  LocMesh->Vertexes[IdNewVortex].Y=(LocMesh->Vertexes[IdNgbhVortex1].Y+LocMesh->Vertexes[IdNgbhVortex2].Y)/2.;
  LocMesh->Vertexes[IdNewVortex].Z=(LocMesh->Vertexes[IdNgbhVortex1].Z+LocMesh->Vertexes[IdNgbhVortex2].Z)/2.;
  
  //update the vertex identifiers of the edge ends
  LocMesh->Edges[IdSubdiviedEdge].V1=IdNgbhVortex1;
  LocMesh->Edges[IdSubdiviedEdge].V2=IdNewVortex;
  
  LocMesh->Edges[IdNewEdge1].V1=IdNewVortex;
  LocMesh->Edges[IdNewEdge1].V2=IdNgbhVortex2;
  LocMesh->Edges[IdNewEdge1].NbSubdiv=1;
  
  LocMesh->Edges[IdNewEdge2].V1=IdNewVortex;
  LocMesh->Edges[IdNewEdge2].V2=IdNgbhVortex4;
  LocMesh->Edges[IdNewEdge2].NbSubdiv=1;

  LocMesh->Edges[IdNewEdge3].V1=IdNewVortex;
  LocMesh->Edges[IdNewEdge3].V2=IdNgbhVortex3;
  LocMesh->Edges[IdNewEdge3].NbSubdiv=1;


  //update the edges identifiers of the elements
  LocMesh->Elements[IdSubdiviedElement1].E1=IdSubdiviedEdge;
  LocMesh->Elements[IdSubdiviedElement1].E2=IdNewEdge3;
  LocMesh->Elements[IdSubdiviedElement1].E3=IdNgbhEdge1;
  
  LocMesh->Elements[IdSubdiviedElement2].E1=IdNgbhEdge2;
  LocMesh->Elements[IdSubdiviedElement2].E2=IdNewEdge2;
  LocMesh->Elements[IdSubdiviedElement2].E3=IdSubdiviedEdge;
  
  LocMesh->Elements[IdNewElement1].E1=IdNewEdge1;
  LocMesh->Elements[IdNewElement1].E2=IdNgbhEdge3;
  LocMesh->Elements[IdNewElement1].E3=IdNewEdge3;
  
  LocMesh->Elements[IdNewElement2].E1=IdNgbhEdge4;
  LocMesh->Elements[IdNewElement2].E2=IdNewEdge1;
  LocMesh->Elements[IdNewElement2].E3=IdNewEdge2;
}

//Project de coordinates of a mesh on a ellispoid of radius 'Rx,Ry,Rz'
void ProjectMeshEllipsoid(Mesh * LocMesh, double Rx, double Ry, double Rz){
  int i;
  double Xmax,Xmin,Ymax,Ymin,Zmax,Zmin,Xmean,Ymean,Zmean;
  double LocNorm;
  
  //translate the mesh centre to the origin
  Xmax=LocMesh->Vertexes[0].X;  Xmin=LocMesh->Vertexes[0].X;
  Ymax=LocMesh->Vertexes[0].Y;  Ymin=LocMesh->Vertexes[0].Y;
  Zmax=LocMesh->Vertexes[0].Z;  Zmin=LocMesh->Vertexes[0].Z;
  
  for (i=1;i<LocMesh->NbVertexes;i++){
    if (LocMesh->Vertexes[i].X>Xmax) Xmax=LocMesh->Vertexes[i].X;
    if (LocMesh->Vertexes[i].X<Xmin) Xmin=LocMesh->Vertexes[i].X;
    if (LocMesh->Vertexes[i].Y>Ymax) Ymax=LocMesh->Vertexes[i].Y;
    if (LocMesh->Vertexes[i].Y<Ymin) Ymin=LocMesh->Vertexes[i].Y;
    if (LocMesh->Vertexes[i].Z>Zmax) Zmax=LocMesh->Vertexes[i].Z;
    if (LocMesh->Vertexes[i].Z<Zmin) Zmin=LocMesh->Vertexes[i].Z;
  }
  
  Xmean=(Xmax+Xmin)/2;
  Ymean=(Ymax+Ymin)/2;
  Zmean=(Zmax+Zmin)/2;
  
  for (i=0;i<LocMesh->NbVertexes;i++){
    LocMesh->Vertexes[i].X-=Xmean;
    LocMesh->Vertexes[i].Y-=Ymean;
    LocMesh->Vertexes[i].Z-=Zmean;
  }
  
  //project each vertex on a ellipsoid of radius 'Rx,Ry,Rz'
  for (i=0;i<LocMesh->NbVertexes;i++){
    LocNorm=sqrt((pow(LocMesh->Vertexes[i].X,2.0)/pow(Rx,2.0))+
                 (pow(LocMesh->Vertexes[i].Y,2.0)/pow(Ry,2.0))+
                 (pow(LocMesh->Vertexes[i].Z,2.0)/pow(Rz,2.0)));
    LocMesh->Vertexes[i].X=LocMesh->Vertexes[i].X/LocNorm;
    LocMesh->Vertexes[i].Y=LocMesh->Vertexes[i].Y/LocNorm;
    LocMesh->Vertexes[i].Z=LocMesh->Vertexes[i].Z/LocNorm;
  }
}


//order the edges of all element so that this order is clockwise if we see it from outside the shape
void MakeClockwiseOrder(Mesh * LocMesh){
  int i,j,jswap;
  int Ve1,Ed1,Ve2,Ed2,Ve3,Ed3;  //Order: Vertex 1 -> Edge1 -> Vertex 2 -> Edge 2 -> Vertex 3 -> Edge 3 -> Vertex 1 ...
  int temp;
  double x1,x2,x3,y1,y2,y3,z1,z2,z3; //coordinates
  double u1,u2,u3,v1,v2,v3;  //vectors
  double vp1,vp2,vp3;  //cross product
  double tempDbl;
  
  //test all elements of the mesh
  for (i=0;i<LocMesh->NbElements;i++){
    //order the vertexes and edges
    Ed1=LocMesh->Elements[i].E1;
    Ed2=LocMesh->Elements[i].E2;
    Ed3=LocMesh->Elements[i].E3;
    
    if ((LocMesh->Edges[Ed1].V2!=LocMesh->Edges[Ed2].V1)&&(LocMesh->Edges[Ed1].V2!=LocMesh->Edges[Ed2].V2)){ //swap edges
      temp=Ed2; Ed2=Ed3; Ed3=temp;
    }
    
    if (LocMesh->Edges[Ed1].V2!=LocMesh->Edges[Ed2].V1){ //swap vertexes
      temp=LocMesh->Edges[Ed2].V1; LocMesh->Edges[Ed2].V1=LocMesh->Edges[Ed2].V2; LocMesh->Edges[Ed2].V2=temp;
      if (LocMesh->Edges[Ed2].NbSubdiv>1)
        for (j=0;j<=LocMesh->Edges[Ed2].NbSubdiv/2;j++){
          jswap=LocMesh->Edges[Ed2].NbSubdiv-j;
          tempDbl=LocMesh->Edges[Ed2].SubiVertex[j].X; LocMesh->Edges[Ed2].SubiVertex[j].X=LocMesh->Edges[Ed2].SubiVertex[jswap].X;  LocMesh->Edges[Ed2].SubiVertex[jswap].X=tempDbl;
          tempDbl=LocMesh->Edges[Ed2].SubiVertex[j].Y; LocMesh->Edges[Ed2].SubiVertex[j].Y=LocMesh->Edges[Ed2].SubiVertex[jswap].Y;  LocMesh->Edges[Ed2].SubiVertex[jswap].Y=tempDbl;
          tempDbl=LocMesh->Edges[Ed2].SubiVertex[j].Z; LocMesh->Edges[Ed2].SubiVertex[j].Z=LocMesh->Edges[Ed2].SubiVertex[jswap].Z;  LocMesh->Edges[Ed2].SubiVertex[jswap].Z=tempDbl;
        }
    }
    
    if (LocMesh->Edges[Ed2].V2!=LocMesh->Edges[Ed3].V1){ //swap vertexes
      temp=LocMesh->Edges[Ed3].V1; LocMesh->Edges[Ed3].V1=LocMesh->Edges[Ed3].V2; LocMesh->Edges[Ed3].V2=temp;
      if (LocMesh->Edges[Ed3].NbSubdiv>1)
        for (j=0;j<=LocMesh->Edges[Ed3].NbSubdiv/2;j++){
          jswap=LocMesh->Edges[Ed3].NbSubdiv-j;
          tempDbl=LocMesh->Edges[Ed3].SubiVertex[j].X; LocMesh->Edges[Ed3].SubiVertex[j].X=LocMesh->Edges[Ed3].SubiVertex[jswap].X;  LocMesh->Edges[Ed3].SubiVertex[jswap].X=tempDbl;
          tempDbl=LocMesh->Edges[Ed3].SubiVertex[j].Y; LocMesh->Edges[Ed3].SubiVertex[j].Y=LocMesh->Edges[Ed3].SubiVertex[jswap].Y;  LocMesh->Edges[Ed3].SubiVertex[jswap].Y=tempDbl;
          tempDbl=LocMesh->Edges[Ed3].SubiVertex[j].Z; LocMesh->Edges[Ed3].SubiVertex[j].Z=LocMesh->Edges[Ed3].SubiVertex[jswap].Z;  LocMesh->Edges[Ed3].SubiVertex[jswap].Z=tempDbl;
        }
    }
    
    if (LocMesh->Edges[Ed3].V2!=LocMesh->Edges[Ed1].V1){
      printf("This is not good!!!\n");
    }
    
    //invert the order if necessary
    x1=LocMesh->Vertexes[LocMesh->Edges[Ed1].V1].X;  y1=LocMesh->Vertexes[LocMesh->Edges[Ed1].V1].Y;  z1=LocMesh->Vertexes[LocMesh->Edges[Ed1].V1].Z;
    x2=LocMesh->Vertexes[LocMesh->Edges[Ed1].V2].X;  y2=LocMesh->Vertexes[LocMesh->Edges[Ed1].V2].Y;  z2=LocMesh->Vertexes[LocMesh->Edges[Ed1].V2].Z;
    x3=LocMesh->Vertexes[LocMesh->Edges[Ed3].V1].X;  y3=LocMesh->Vertexes[LocMesh->Edges[Ed3].V1].Y;  z3=LocMesh->Vertexes[LocMesh->Edges[Ed3].V1].Z;
    
    u1=x2-x1; u2=y2-y1; u3=z2-z1;
    v1=x3-x1; v2=y3-y1; v3=z3-z1;
    
    vp1=u2*v3-u3*v2;
    vp2=u3*v1-u1*v3;
    vp3=u1*v2-u2*v1;
    
    //printf("%3.2lf %3.2lf %3.2lf | %3.2lf %3.2lf %3.2lf | %3.2lf %3.2lf %3.2lf || %3.2lf %3.2lf %3.2lf || %3.2lf %3.2lf %3.2lf | %lf\n",x1,y1,z1,x2,y2,z2,x3,y3,z3,u1,u2,u3,v1,v2,v3,x1*vp1+y1*vp2+z1*vp3);
    //printf("%3.2lf %3.2lf %3.2lf | %3.2lf %3.2lf %3.2lf | %3.2lf\n",x1,y1,z1,vp1,vp2,vp3,x1*vp1+y1*vp2+z1*vp3);
    
    if (x1*vp1+y1*vp2+z1*vp3<0){//we consider that the origin is within the spheric shape so the vector (x1,y1,z1) points out of the shape
      //swap edges
      temp=Ed3; Ed3=Ed2; Ed2=temp;
      //swap vertexes
      temp=LocMesh->Edges[Ed1].V1; LocMesh->Edges[Ed1].V1=LocMesh->Edges[Ed1].V2; LocMesh->Edges[Ed1].V2=temp;
      if (LocMesh->Edges[Ed1].NbSubdiv>1)
        for (j=0;j<=LocMesh->Edges[Ed1].NbSubdiv/2;j++){
          jswap=LocMesh->Edges[Ed1].NbSubdiv-j;
          tempDbl=LocMesh->Edges[Ed1].SubiVertex[j].X; LocMesh->Edges[Ed1].SubiVertex[j].X=LocMesh->Edges[Ed1].SubiVertex[jswap].X;  LocMesh->Edges[Ed1].SubiVertex[jswap].X=tempDbl;
          tempDbl=LocMesh->Edges[Ed1].SubiVertex[j].Y; LocMesh->Edges[Ed1].SubiVertex[j].Y=LocMesh->Edges[Ed1].SubiVertex[jswap].Y;  LocMesh->Edges[Ed1].SubiVertex[jswap].Y=tempDbl;
          tempDbl=LocMesh->Edges[Ed1].SubiVertex[j].Z; LocMesh->Edges[Ed1].SubiVertex[j].Z=LocMesh->Edges[Ed1].SubiVertex[jswap].Z;  LocMesh->Edges[Ed1].SubiVertex[jswap].Z=tempDbl;
        }
    }
    //new order of the elements
    LocMesh->Elements[i].E1=Ed1;
    LocMesh->Elements[i].E2=Ed2;
    LocMesh->Elements[i].E3=Ed3;
  }

}

//subdivide into 2 edges all edges in the mesh 'LocMesh'
void RefineMesh(Mesh * LocMesh){
  int i,j,k;
  int * order;
  double SqLengthI,SqLengthJ;
  int OrigNbEdges;
  
  order=(int*)malloc(LocMesh->NbEdges*sizeof(int));
  
  for (i=0;i<LocMesh->NbEdges;i++) order[i]=i;
  
  for (i=0;i<LocMesh->NbEdges-1;i++)
    for (j=i+1;j<LocMesh->NbEdges;j++){
    SqLengthI=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].X,2.0);
    SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Y,2.0);
    SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Z,2.0);
    SqLengthJ=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].X,2.0);
    SqLengthJ+=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].Y,2.0);
    SqLengthJ+=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].Z,2.0);
    
    if (SqLengthI<SqLengthJ){
      k=order[i];
      order[i]=order[j];
      order[j]=k;
    }
    }
  
    OrigNbEdges=LocMesh->NbEdges;
    for (i=0;i<OrigNbEdges;i++)  CutEdge(LocMesh,order[i]);
}


//Project de coordinates of a mesh on an ellipsoid of radius 'Rx,Ry,Rz' and refine the mesh until 
//the longest edge has a length smaller than 'MaxEdgeLength'.
//If NbSubdiv>1, each edge is curved into 'NbSubdiv' parts in the end
void ProjectMeshEllipsoidAndRefine(Mesh * LocMesh, double Rx,double Ry,double Rz,double MaxEdgeLength,int NbSubdiv){
  int i,j,k;
  int * order;
  double SqLengthI,SqLengthJ;
  int OrigNbEdges;
  int OK;
  double LocNorm;
  double dj,dNbSubdiv;
  int jm1;
  double EdgeV1_X,EdgeV1_Y,EdgeV1_Z,EdgeV2_X,EdgeV2_Y,EdgeV2_Z;
  
  //initialize
  order=(int*)malloc(LocMesh->NbEdges*sizeof(int));
  for (i=0;i<LocMesh->NbEdges;i++) order[i]=i;
  OK=1;
  
  //mesh projection on an ellipsoid
  ProjectMeshEllipsoid(LocMesh,Rx,Ry,Rz);
  
  //Big loop to subdivide the mesh until each edge has a length smaller than 'MaxEdgeLength'
  while (OK==1){
    OK=0;
    
    order=(int*)realloc(order,LocMesh->NbEdges*sizeof(int));
    for (i=0;i<LocMesh->NbEdges;i++) order[i]=i;

    //order the edges as a function of their length
    for (i=0;i<LocMesh->NbEdges-1;i++) for (j=i+1;j<LocMesh->NbEdges;j++){
      SqLengthI=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].X,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Y,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Z,2.0);
      SqLengthJ=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].X,2.0);
      SqLengthJ+=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].Y,2.0);
      SqLengthJ+=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].Z,2.0);
      
      if (SqLengthI<SqLengthJ){ k=order[i]; order[i]=order[j]; order[j]=k;}
    }
    
    //subdivide the edges if their length is larger than 'MaxEdgeLength'
    OrigNbEdges=LocMesh->NbEdges;
    for (i=0;i<OrigNbEdges;i++){
      SqLengthI=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].X,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Y,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Z,2.0);
      
      if (sqrt(SqLengthI)>MaxEdgeLength){
        //subdvision
        OK=1;
        CutEdge(LocMesh,order[i]);    //rem: 'LocMesh' is changed here. In particular LocMesh->NbEdges is larger.
      }
    }
    
    //projection of the new edges on an ellipsoid
    ProjectMeshEllipsoid(LocMesh,Rx,Ry,Rz);
  }
  
  
  
  //subdivide the edges end project them on the ellipsoid if required
  if (NbSubdiv>1){
    for (i=0;i<LocMesh->NbEdges;i++){
      //init
      LocMesh->Edges[i].NbSubdiv=NbSubdiv;
      LocMesh->Edges[i].SubiVertex=(Vertex *)malloc((LocMesh->Edges[i].NbSubdiv+1)*sizeof(Vertex));
      
      //coordinates of the subdivised edge
      EdgeV1_X=LocMesh->Vertexes[LocMesh->Edges[i].V1].X;
      EdgeV1_Y=LocMesh->Vertexes[LocMesh->Edges[i].V1].Y;
      EdgeV1_Z=LocMesh->Vertexes[LocMesh->Edges[i].V1].Z;
      EdgeV2_X=LocMesh->Vertexes[LocMesh->Edges[i].V2].X;
      EdgeV2_Y=LocMesh->Vertexes[LocMesh->Edges[i].V2].Y;
      EdgeV2_Z=LocMesh->Vertexes[LocMesh->Edges[i].V2].Z;
      dNbSubdiv=(double)LocMesh->Edges[i].NbSubdiv;
      
      for (j=0;j<LocMesh->Edges[i].NbSubdiv+1;j++){
        dj=(double)j;
        
        //interpolation
        LocMesh->Edges[i].SubiVertex[j].X=EdgeV1_X*(dNbSubdiv-dj)/dNbSubdiv+EdgeV2_X*dj/dNbSubdiv;
        LocMesh->Edges[i].SubiVertex[j].Y=EdgeV1_Y*(dNbSubdiv-dj)/dNbSubdiv+EdgeV2_Y*dj/dNbSubdiv;
        LocMesh->Edges[i].SubiVertex[j].Z=EdgeV1_Z*(dNbSubdiv-dj)/dNbSubdiv+EdgeV2_Z*dj/dNbSubdiv;
        
        //projection on an ellipsoid
        LocNorm=sqrt((pow(LocMesh->Edges[i].SubiVertex[j].X,2.0)/pow(Rx,2.0))+
                     (pow(LocMesh->Edges[i].SubiVertex[j].Y,2.0)/pow(Ry,2.0))+
                     (pow(LocMesh->Edges[i].SubiVertex[j].Z,2.0)/pow(Rz,2.0)));
        LocMesh->Edges[i].SubiVertex[j].X=LocMesh->Edges[i].SubiVertex[j].X/LocNorm;
        LocMesh->Edges[i].SubiVertex[j].Y=LocMesh->Edges[i].SubiVertex[j].Y/LocNorm;
        LocMesh->Edges[i].SubiVertex[j].Z=LocMesh->Edges[i].SubiVertex[j].Z/LocNorm;
      }
    }
  }
}


//Extract a submesh of 'LocMesh' in the bounding box [Xmin,Xmax],[Ymin,Ymax],[Zmin,Zmax]  (NEW)
void ExtractSubmesh(Mesh * LocMesh,double Xmin,double Xmax,double Ymin,double Ymax,double Zmin,double Zmax){
  int i,j,NbCurvedEdge;
  int * VertexesToKeep;
  int * EdgesToKeep;
  int * ElementsToKeep;
  int NewValue;
  int NewI;
  double epsilon;
  Mesh TempMesh;

  //1) allocations and initializations
  epsilon=0.001;
  
  VertexesToKeep=(int*)malloc(LocMesh->NbVertexes*sizeof(int));
  for (i=0;i<LocMesh->NbVertexes;i++) VertexesToKeep[i]=0;

  EdgesToKeep=(int*)malloc(LocMesh->NbEdges*sizeof(int));
  for (i=0;i<LocMesh->NbEdges;i++) EdgesToKeep[i]=0;

  ElementsToKeep=(int*)malloc(LocMesh->NbElements*sizeof(int));
  for (i=0;i<LocMesh->NbElements;i++) ElementsToKeep[i]=0;

  //2) find the vertexes to keep
  for (i=0;i<LocMesh->NbVertexes;i++)
    if ((LocMesh->Vertexes[i].X>Xmin-epsilon)&&(LocMesh->Vertexes[i].X<Xmax+epsilon)&&
        (LocMesh->Vertexes[i].Y>Ymin-epsilon)&&(LocMesh->Vertexes[i].Y<Ymax+epsilon)&&
        (LocMesh->Vertexes[i].Z>Zmin-epsilon)&&(LocMesh->Vertexes[i].Z<Zmax+epsilon)) VertexesToKeep[i]=1;

  //3) find the edges to keep
  for (i=0;i<LocMesh->NbEdges;i++) if ((VertexesToKeep[LocMesh->Edges[i].V1]==1)&&(VertexesToKeep[LocMesh->Edges[i].V2]==1))
    EdgesToKeep[i]=1;

  //4) find the elements to keep
  for (i=0;i<LocMesh->NbElements;i++) if ((EdgesToKeep[LocMesh->Elements[i].E1]==1)&&(EdgesToKeep[LocMesh->Elements[i].E2]==1)&&(EdgesToKeep[LocMesh->Elements[i].E3]==1))
    ElementsToKeep[i]=1;

  
  //5) extraction of the submesh
  
  //5.1) new IDs   (ID+1 actually)
  NewValue=1;
  for (i=0;i<LocMesh->NbVertexes;i++) if (VertexesToKeep[i]==1){
    VertexesToKeep[i]=NewValue;
    NewValue++;
  }

  NewValue=1;
  for (i=0;i<LocMesh->NbEdges;i++) if (EdgesToKeep[i]==1){
    EdgesToKeep[i]=NewValue;
    NewValue++;
  }

  NewValue=1;
  for (i=0;i<LocMesh->NbElements;i++) if (ElementsToKeep[i]==1){
    ElementsToKeep[i]=NewValue;
    NewValue++;
  }
  
  //5.2) copy the mesh in a temporary mesh
  TempMesh.NbVertexes=LocMesh->NbVertexes;
  TempMesh.Vertexes=(Vertex *)malloc(TempMesh.NbVertexes*sizeof(Vertex));
  for (i=0;i<TempMesh.NbVertexes;i++){
    TempMesh.Vertexes[i].X=LocMesh->Vertexes[i].X;
    TempMesh.Vertexes[i].Y=LocMesh->Vertexes[i].Y;
    TempMesh.Vertexes[i].Z=LocMesh->Vertexes[i].Z;
  }
 	
  TempMesh.NbEdges=LocMesh->NbEdges;
  TempMesh.Edges=(Edge *)malloc(TempMesh.NbEdges*sizeof(Edge));
  for (i=0;i<TempMesh.NbEdges;i++){
    TempMesh.Edges[i].V1=LocMesh->Edges[i].V1;
    TempMesh.Edges[i].V2=LocMesh->Edges[i].V2;
  }
  for (i=0;i<LocMesh->NbEdges;i++) if (LocMesh->Edges[i].NbSubdiv>1){
    TempMesh.Edges[i].NbSubdiv=LocMesh->Edges[i].NbSubdiv;
    TempMesh.Edges[i].SubiVertex=(Vertex *)malloc((TempMesh.Edges[i].NbSubdiv+1)*sizeof(Vertex));
    for (j=0;j<TempMesh.Edges[i].NbSubdiv+1;j++){
      TempMesh.Edges[i].SubiVertex[j].X=LocMesh->Edges[i].SubiVertex[j].X;
      TempMesh.Edges[i].SubiVertex[j].Y=LocMesh->Edges[i].SubiVertex[j].Y;
      TempMesh.Edges[i].SubiVertex[j].Z=LocMesh->Edges[i].SubiVertex[j].Z;
    }
  }
  
  TempMesh.NbElements=LocMesh->NbElements;
  TempMesh.Elements=(Element *)malloc(TempMesh.NbElements*sizeof(Element));
  for (i=0;i<TempMesh.NbElements;i++){
    TempMesh.Elements[i].E1=LocMesh->Elements[i].E1;
    TempMesh.Elements[i].E2=LocMesh->Elements[i].E2;
    TempMesh.Elements[i].E3=LocMesh->Elements[i].E3;
  }

  //5.3) change the mesh
  NewI=0;
  for (i=0;i<TempMesh.NbVertexes;i++) if (VertexesToKeep[i]>0){
    LocMesh->Vertexes[NewI].X=TempMesh.Vertexes[i].X;
    LocMesh->Vertexes[NewI].Y=TempMesh.Vertexes[i].Y;
    LocMesh->Vertexes[NewI].Z=TempMesh.Vertexes[i].Z;
    //printf("Vertex %d - %f %f %f\n",NewI,LocMesh->Vertexes[NewI].X,LocMesh->Vertexes[NewI].Y,LocMesh->Vertexes[NewI].Z);
    NewI++;
  }
  LocMesh->NbVertexes=NewI;

  NewI=0;
  for (i=0;i<TempMesh.NbEdges;i++) if (EdgesToKeep[i]>0){
    LocMesh->Edges[NewI].V1=VertexesToKeep[TempMesh.Edges[i].V1]-1;
    LocMesh->Edges[NewI].V2=VertexesToKeep[TempMesh.Edges[i].V2]-1;
    if (TempMesh.Edges[i].NbSubdiv>1) for (j=0;j<TempMesh.Edges[i].NbSubdiv+1;j++){
      LocMesh->Edges[NewI].SubiVertex[j].X=TempMesh.Edges[i].SubiVertex[j].X;
      LocMesh->Edges[NewI].SubiVertex[j].Y=TempMesh.Edges[i].SubiVertex[j].Y;
      LocMesh->Edges[NewI].SubiVertex[j].Z=TempMesh.Edges[i].SubiVertex[j].Z;
    } //we suppose that all edges have the same number of subdivisions
    //printf("Edge %d - %d %d\n",NewI,LocMesh->Edges[NewI].V1,LocMesh->Edges[NewI].V2);
    NewI++;
  }
  LocMesh->NbEdges=NewI;

  NewI=0;
  for (i=0;i<TempMesh.NbElements;i++) if (ElementsToKeep[i]>0){
    LocMesh->Elements[NewI].E1=EdgesToKeep[TempMesh.Elements[i].E1]-1;
    LocMesh->Elements[NewI].E2=EdgesToKeep[TempMesh.Elements[i].E2]-1;
    LocMesh->Elements[NewI].E3=EdgesToKeep[TempMesh.Elements[i].E3]-1;
    //printf("Element %d - %d %d %d \n",NewI,LocMesh->Elements[NewI].E1,LocMesh->Elements[NewI].E2,LocMesh->Elements[NewI].E3);
    NewI++;
  }
  LocMesh->NbElements=NewI;
  }



//write the geometry of mesh contained in a Mesh structure in a XML-Nektar file  (NEW STUFFS THERE)
void WriteMesh(Mesh * LocMesh,char FileName[256]){
  FILE * XmlMeshGeomFile;
  int i,j,NbCurvedEdge;
  int * BoundaryEdgesList;
  int FirstBoundary;
  
  //to have a proper order of the edges within the elements
  MakeClockwiseOrder(LocMesh);
  
  //open file
  XmlMeshGeomFile = fopen(FileName,"w");
  
  //write header
  fprintf(XmlMeshGeomFile,"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
  fprintf(XmlMeshGeomFile,"\n");
  fprintf(XmlMeshGeomFile,"<NEKTAR>\n");
  fprintf(XmlMeshGeomFile,"<!-- Embed a 2-dimensional object in a 2-dimensional space -->\n");
  fprintf(XmlMeshGeomFile,"<!-- DIM <= SPACE -->\n");
  fprintf(XmlMeshGeomFile,"<!-- This provides a method of optimizing code for a 1-D curve embedded in 3-space. -->\n");
  fprintf(XmlMeshGeomFile,"<GEOMETRY DIM=\"2\" SPACE=\"3\">\n");
  fprintf(XmlMeshGeomFile,"\n");
  
  //write vertex
  fprintf(XmlMeshGeomFile,"  <VERTEX>\n");
  fprintf(XmlMeshGeomFile,"    <!-- Always must have four values per entry. -->\n");
  for (i=0;i<LocMesh->NbVertexes;i++)
    fprintf(XmlMeshGeomFile,"    <V ID=\"%d\">  %f   %f   %f  </V>\n",i,LocMesh->Vertexes[i].X,LocMesh->Vertexes[i].Y,LocMesh->Vertexes[i].Z);
  fprintf(XmlMeshGeomFile,"  </VERTEX>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  //write edges
  fprintf(XmlMeshGeomFile,"  <EDGE>\n");
  fprintf(XmlMeshGeomFile,"    <!--Edges are vertex pairs -->\n");
  for (i=0;i<LocMesh->NbEdges;i++)
    fprintf(XmlMeshGeomFile,"    <E ID=\"%d\">    %d  %d   </E>\n",i,LocMesh->Edges[i].V1,LocMesh->Edges[i].V2);
  fprintf(XmlMeshGeomFile,"  </EDGE>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  //write Element
  fprintf(XmlMeshGeomFile,"  <ELEMENT>\n");
  for (i=0;i<LocMesh->NbElements;i++)
    fprintf(XmlMeshGeomFile,"    <T ID=\"%d\">    %d     %d     %d </T>\n",i,LocMesh->Elements[i].E1,LocMesh->Elements[i].E2,LocMesh->Elements[i].E3);
  fprintf(XmlMeshGeomFile,"  </ELEMENT>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  
  
  //write curved elements
  fprintf(XmlMeshGeomFile,"  <CURVED>\n");
  NbCurvedEdge=0;
  for (i=0;i<LocMesh->NbEdges;i++) if (LocMesh->Edges[i].NbSubdiv>1){
    fprintf(XmlMeshGeomFile,"    <E ID=\"%d\" EDGEID=\"%d\" TYPE=\"PolyEvenlySpaced\" NUMPOINTS=\"%d\">\n",NbCurvedEdge,i,LocMesh->Edges[i].NbSubdiv+1);
    for (j=0;j<LocMesh->Edges[i].NbSubdiv+1;j++){
      fprintf(XmlMeshGeomFile,"    %lf %lf %lf\n",LocMesh->Edges[i].SubiVertex[j].X,LocMesh->Edges[i].SubiVertex[j].Y,LocMesh->Edges[i].SubiVertex[j].Z);
    }
    fprintf(XmlMeshGeomFile,"   </E>\n\n");
    NbCurvedEdge++;
  }
  fprintf(XmlMeshGeomFile,"   </CURVED>\n");
  
  
  //find the boundaries   (points at 0 in BoundaryEdgesList)
  BoundaryEdgesList=(int *)malloc(LocMesh->NbEdges*sizeof(int));
  
  for (i=0;i<LocMesh->NbEdges;i++) BoundaryEdgesList[i]=-1;
  
  for (i=0;i<LocMesh->NbElements;i++){
    if (BoundaryEdgesList[LocMesh->Elements[i].E1]==1) printf("Problem between the elements and the edges\n");
    if (BoundaryEdgesList[LocMesh->Elements[i].E2]==1) printf("Problem between the elements and the edges\n");
    if (BoundaryEdgesList[LocMesh->Elements[i].E3]==1) printf("Problem between the elements and the edges\n");
    if (BoundaryEdgesList[LocMesh->Elements[i].E1]==0) BoundaryEdgesList[LocMesh->Elements[i].E1]=1;
    if (BoundaryEdgesList[LocMesh->Elements[i].E2]==0) BoundaryEdgesList[LocMesh->Elements[i].E2]=1;
    if (BoundaryEdgesList[LocMesh->Elements[i].E3]==0) BoundaryEdgesList[LocMesh->Elements[i].E3]=1;
    if (BoundaryEdgesList[LocMesh->Elements[i].E1]==-1) BoundaryEdgesList[LocMesh->Elements[i].E1]=0;
    if (BoundaryEdgesList[LocMesh->Elements[i].E2]==-1) BoundaryEdgesList[LocMesh->Elements[i].E2]=0;
    if (BoundaryEdgesList[LocMesh->Elements[i].E3]==-1) BoundaryEdgesList[LocMesh->Elements[i].E3]=0;
  }

  //write the end (could be more evolved)
  fprintf(XmlMeshGeomFile,"<!-- V - vertex, E - edge, F - face, L - element -->\n");
  fprintf(XmlMeshGeomFile,"  <COMPOSITE>\n");
  fprintf(XmlMeshGeomFile,"    <C ID=\"0\"> T[0-%d]\n",LocMesh->NbElements-1);
  fprintf(XmlMeshGeomFile,"    </C>\n");
  fprintf(XmlMeshGeomFile,"    <C ID=\"1\"> E[");
  FirstBoundary=1;
  for (i=0;i<LocMesh->NbEdges;i++) if (BoundaryEdgesList[i]==0){
    if (FirstBoundary==1) FirstBoundary=0;
    else fprintf(XmlMeshGeomFile,",");
    fprintf(XmlMeshGeomFile,"%d",i);
  }
  fprintf(XmlMeshGeomFile,"]\n");
  fprintf(XmlMeshGeomFile,"    </C>\n");
  fprintf(XmlMeshGeomFile,"  </COMPOSITE>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  fprintf(XmlMeshGeomFile,"  <DOMAIN> C[0] </DOMAIN>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  //EOF
  fprintf(XmlMeshGeomFile,"</GEOMETRY>\n");
  
  //close file
  fclose(XmlMeshGeomFile);
}





// ---------------------------------------------------------------------------------------
//                                       main function
// ---------------------------------------------------------------------------------------


void usage()
{
  printf("Usage: makeEllipsoidNektar <options>\n");
  printf("Where <options> are one or more of the following:\n");
  printf("\t<-NameFile n>     Name of the output file containg the ellipsoid (default=\"toto.xml\")\n");
  printf("\t<-MaxLength n>    Maximum length allowed for an edge (default=0.5)\n");
  printf("\t<-NbSubdiv n>     Number of subdivisions of each edge to make it curved (default=1)\n");
  printf("\t<-RadiusX n>      Radius of the ellipsoid in the X direction  (default=1.)\n");
  printf("\t<-RadiusY n>      Radius of the ellipsoid in the Y direction  (default=1.)\n");
  printf("\t<-RadiusZ n>      Radius of the ellipsoid in the Z direction  (default=1.)\n");
  printf("\t<-MinX n>         Bounding box - min X                        (default=-1000.)\n");
  printf("\t<-MaxX n>         Bounding box - max X                        (default=1000.)\n");
  printf("\t<-MinY n>         Bounding box - min Y                        (default=-1000.)\n");
  printf("\t<-MaxY n>         Bounding box - max Y                        (default=1000.)\n");
  printf("\t<-MinZ n>         Bounding box - min Z                        (default=-1000.)\n");
  printf("\t<-MaxZ n>         Bounding box - max Z                        (default=1000.)\n");
  exit(1);
}

int main(int argc, char **argv)
{
  int ok,i;
  char *output_name = "toto.xml";
  double RadiusX=1.;
  double RadiusY=1.;
  double RadiusZ=1.;
  double MinX=-1000;
  double MaxX=1000;
  double MinY=-1000;
  double MaxY=1000;
  double MinZ=-1000;
  double MaxZ=1000;
  double Prec=0.5;
  int NbSubdiv=1;
  Mesh LocMesh;


  //1) Parse parameters
  while (argc > 1) {
    ok = 0;
    if ((ok == 0) && (strcmp(argv[1], "help") == 0)) {
     usage();
    }
    if ((ok == 0) && (strcmp(argv[1], "-NameFile") == 0)) {
      argc--;
      argv++;
      output_name  = argv[1];
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-RadiusX") == 0)) {
      argc--;
      argv++;
      RadiusX = atof(argv[1]);
      if (RadiusX<0) RadiusX=1.;
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-RadiusY") == 0)) {
      argc--;
      argv++;
      RadiusY = atof(argv[1]);
      if (RadiusY<0) RadiusY=1.;
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-RadiusZ") == 0)) {
      argc--;
      argv++;
      RadiusZ = atof(argv[1]);
      if (RadiusZ<0) RadiusZ=1.;
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-MinX") == 0)) {
      argc--;
      argv++;
      MinX = atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-MaxX") == 0)) {
      argc--;
      argv++;
      MaxX = atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-MinY") == 0)) {
      argc--;
      argv++;
      MinY = atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-MaxY") == 0)) {
      argc--;
      argv++;
      MaxY = atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-MinZ") == 0)) {
      argc--;
      argv++;
      MinZ = atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-MaxZ") == 0)) {
      argc--;
      argv++;
      MaxZ = atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-MaxLength") == 0)) {
      argc--;
      argv++;
      Prec = atof(argv[1]);
      if ((Prec<(RadiusX+RadiusY+RadiusZ)*0.001/3.)) Prec=(RadiusX+RadiusY+RadiusZ)*0.001/3;
      argc--;
      argv++;
      ok = 1;
    }
    if ((ok == 0) && (strcmp(argv[1], "-NbSubdiv") == 0)) {
      argc--;
      argv++;
      NbSubdiv = atoi(argv[1]);
      if ((NbSubdiv<1)) NbSubdiv=1;
      argc--;
      argv++;
      ok = 1;
    }
    if (ok == 0) {
      printf("Unknown option: \n");
      usage();
    }
  }
  
  printf("Output file: %s / Radius: %lf / Radius Y: %lf / Radius Z: %lf / Precision: %lf / Subdivisions: %d\n",output_name,RadiusX,RadiusY,RadiusZ,Prec,NbSubdiv);
  
  
  //2) run functions
  LocMesh=CreateBasicCubicMesh();
  ProjectMeshEllipsoidAndRefine(&LocMesh,RadiusX,RadiusY,RadiusZ,Prec,NbSubdiv);
  ExtractSubmesh(&LocMesh,MinX,MaxX,MinY,MaxY,MinZ,MaxZ);
  WriteMesh(&LocMesh,output_name);
  
  return 0;
}
