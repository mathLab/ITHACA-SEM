#include <cstdlib>
#include <math.h>

#include <StdRegions/StdExpansion2D.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>

#include <LibUtilities/Foundations/Foundations.hpp>

using namespace Nektar;

NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2);
NekDouble Quad_sol(NekDouble x, NekDouble y, int order1, int order2,
                   LibUtilities::BasisType btype1, LibUtilities::BasisType btype2);

// This routine projects a polynomial or trigonmetric functions which
// has energy in all mdoes of the expansions and reports and error

int main(int argc, char *argv[])
{
    int           i;
    
    int           order1,order2, nq1,nq2;
    LibUtilities::PointsType    Qtype1,Qtype2;
    LibUtilities::BasisType     btype1 =   LibUtilities::eOrtho_A;
    LibUtilities::BasisType     btype2 =   LibUtilities::eOrtho_B;
    LibUtilities::PointsType    NodalType = LibUtilities::eNodalTriElec;
    LibUtilities::ShapeType     regionshape;
    StdRegions::StdExpansion2D *E;
    Array<OneD, NekDouble> sol;
    Array<OneD, NekDouble> coords(8);
    StdRegions::Orientation edgeDir = StdRegions::eForwards;
    
    
    if((argc != 16)&&(argc != 14))
    {
        //       arg[0]    arg[1]   arg[2]  arg[3] arg[4] arg[5] arg[6] arg[7] arg[8] arg[9] arg[10] arg[11] arg[12] arg[13]
        fprintf(stderr,"Usage: Project2D RegionShape Type1 Type2 order1 order2  nq1    nq2     x1,    y1,      x2,     y2,     x3,    y3 \n");
        
        fprintf(stderr,"Example : ./LocProject2D-g 2 4 5 3 3 4 4 .1 -.5 .6 .1 .3 .2 \n" );
        
        fprintf(stderr,"Where RegionShape is an integer value which "
                "dictates the region shape:\n");
        fprintf(stderr,"\t Triangle      = 3\n");
        fprintf(stderr,"\t Quadrilateral = 4\n");
        
        fprintf(stderr,"Where type is an integer value which "
                "dictates the basis as:\n");
        
        fprintf(stderr,"\t Ortho_A    = 1\n");
        fprintf(stderr,"\t Ortho_B    = 2\n");
        fprintf(stderr,"\t Modified_A = 4\n");
        fprintf(stderr,"\t Modified_B = 5\n");
        fprintf(stderr,"\t Fourier    = 7\n");
        fprintf(stderr,"\t Lagrange   = 8\n");
        fprintf(stderr,"\t Legendre   = 9\n");
        fprintf(stderr,"\t Chebyshev  = 10\n");
        fprintf(stderr,"\t Nodal tri (Electro) = 11\n");
        fprintf(stderr,"\t Nodal tri (Fekete)  = 12\n");
        
        fprintf(stderr,"Note type = 3,6 are for three-dimensional basis\n");
        
        fprintf(stderr,"The last series of values are the coordinates\n");
        exit(1);
    }
    
    regionshape = (LibUtilities::ShapeType) atoi(argv[1]);
    
    // Check to see if 2D region
    if((regionshape != LibUtilities::eTriangle)&&(regionshape != LibUtilities::eQuadrilateral))
    {
        NEKERROR(ErrorUtil::efatal,"This shape is not a 2D region");
    }
  
    int btype1_val = atoi(argv[2]);
    int btype2_val = atoi(argv[3]);
    
    if(( btype1_val <= 10)&&( btype2_val <= 10))
    {
        btype1 =   (LibUtilities::BasisType) btype1_val;
        btype2 =   (LibUtilities::BasisType) btype2_val;
    }
    else if(( btype1_val >=11)&&(btype2_val <= 12))
    {
        btype1 =   LibUtilities::eOrtho_A;
        btype2 =   LibUtilities::eOrtho_B;
        
        if(btype1_val == 11)
        {
            NodalType = LibUtilities::eNodalTriElec;
        }
        else
        {
            NodalType = LibUtilities::eNodalTriFekete;
        }
    }
    
    
    // Check to see that correct Expansions are used
    switch(regionshape)
    {
    case LibUtilities::eTriangle:
        if((btype1 == LibUtilities::eOrtho_B)||(btype1 == LibUtilities::eModified_B))
        {
            NEKERROR(ErrorUtil::efatal,
                     "Basis 1 cannot be of type Ortho_B or Modified_B");
        }
        
        break;
    case LibUtilities::eQuadrilateral:
        if((btype1 == LibUtilities::eOrtho_B)||(btype1 == LibUtilities::eOrtho_C)||
           (btype1 == LibUtilities::eModified_B)||(btype1 == LibUtilities::eModified_C))
        {
            NEKERROR(ErrorUtil::efatal,
                     "Basis 1 is for 2 or 3D expansions");
        }
        
        if((btype2 == LibUtilities::eOrtho_B)||(btype2 == LibUtilities::eOrtho_C)||
           (btype2 == LibUtilities::eModified_B)||(btype2 == LibUtilities::eModified_C))
        {
            NEKERROR(ErrorUtil::efatal,
                     "Basis 2 is for 2 or 3D expansions");
        }
        break;
    default:
        NEKERROR(ErrorUtil::efatal, "Not a valid 2D expansion.");
    }
    
    
    order1 =   atoi(argv[4]);
    order2 =   atoi(argv[5]);
    nq1    =   atoi(argv[6]);
    nq2    =   atoi(argv[7]);
    
    sol = Array<OneD, NekDouble>(nq1*nq2);
    
    if(btype1 != LibUtilities::eFourier)
    {
        Qtype1 = LibUtilities::eGaussLobattoLegendre;
    }
    else
    {
        Qtype1 = LibUtilities::eFourierEvenlySpaced;
    }
    
    if(btype2 != LibUtilities::eFourier)
    {
        if (regionshape == LibUtilities::eTriangle) {
            Qtype2 = LibUtilities::eGaussRadauMAlpha1Beta0;
        }
        else
        {
            Qtype2 = LibUtilities::eGaussLobattoLegendre;
      }
    }
    else
    {
        Qtype2 = LibUtilities::eFourierEvenlySpaced;
    }
    
    //-----------------------------------------------
    // Define a 2D expansion based on basis definition
    
    switch(regionshape)
    {
    case LibUtilities::eTriangle:
        {
            coords[0]    =   atof(argv[8]);
            coords[1]    =   atof(argv[9]);
            coords[2]    =   atof(argv[10]);
            coords[3]    =   atof(argv[11]);
            coords[4]    =   atof(argv[12]);
            coords[5]    =   atof(argv[13]);
            
            // Set up coordinates
            SpatialDomains::VertexComponentSharedPtr verts[3];
            const int zero = 0;
            const int one=1;
            const int two=2;
            const double dZero = 0.0;
            verts[0] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,zero,coords[0],coords[1],dZero);
            verts[1] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,one,coords[2],coords[3],dZero);
            verts[2] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,two,coords[4],coords[5],dZero);
            
            // Set up Edges
            SpatialDomains::SegGeomSharedPtr edges[3];
            edges[0] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero,verts[0],verts[1]);
            edges[1] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(one,verts[1],verts[2]);
            edges[2] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(two,verts[2],verts[0]);
            
            StdRegions::Orientation eorient[3];
            eorient[0] = edgeDir;
            eorient[1] = edgeDir;
            eorient[2] = edgeDir;
            
            SpatialDomains::TriGeomSharedPtr geom = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(zero,verts,edges,eorient);
            geom->SetOwnData();
            
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);

            if(btype1_val >= 10)
            {
                E = new LocalRegions::NodalTriExp(Bkey1,Bkey2,NodalType,geom);
            }
            else
            {
                E = new LocalRegions::TriExp(Bkey1,Bkey2,geom);
            }
            
            Array<OneD,NekDouble> x = Array<OneD,NekDouble>(nq1*nq2);
            Array<OneD,NekDouble> y = Array<OneD,NekDouble>(nq1*nq2);
            E->GetCoords(x,y);
            
            //----------------------------------------------
            // Define solution to be projected
            for(i = 0; i < nq1*nq2; ++i)
            {
                sol[i]  = Tri_sol(x[i],y[i],order1,order2);
            }
            //----------------------------------------------
            
        }
        break;
    case  LibUtilities::eQuadrilateral:
        {
            // Gather coordinates
            coords[0]    =   atof(argv[8]);
            coords[1]    =   atof(argv[9]);
            coords[2]    =   atof(argv[10]);
            coords[3]    =   atof(argv[11]);
            coords[4]    =   atof(argv[12]);
            coords[5]    =   atof(argv[13]);
            coords[6]    =   atof(argv[14]);
            coords[7]    =   atof(argv[15]);
            
            // Set up coordinates
            const int zero=0;
            const int one=1;
            const int two=2;
            const int three=3;
            const double dZero=0.0;
            SpatialDomains::VertexComponentSharedPtr verts[4];
            verts[0] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,zero,coords[0],coords[1],dZero);
            verts[1] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,one,coords[2],coords[3],dZero);
            verts[2] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,two,coords[4],coords[5],dZero);
            verts[3] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,three,coords[6],coords[7],dZero);
            
            // Set up Edges
            SpatialDomains::SegGeomSharedPtr edges[4];
            edges[0] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero,verts[0],verts[1]);
            edges[1] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(one,verts[1],verts[2]);
            edges[2] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(two,verts[2],verts[3]);
            edges[3] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(three,verts[3],verts[0]);
            
            StdRegions::Orientation eorient[4];
            eorient[0] = edgeDir;
            eorient[1] = edgeDir;
            eorient[2] = edgeDir;
            eorient[3] = edgeDir;
            
            SpatialDomains::QuadGeomSharedPtr geom = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(zero,verts,edges,eorient);
            geom->SetOwnData();
            
            const LibUtilities::PointsKey Pkey1(nq1,Qtype1);
            const LibUtilities::PointsKey Pkey2(nq2,Qtype2);
            const LibUtilities::BasisKey  Bkey1(btype1,order1,Pkey1);
            const LibUtilities::BasisKey  Bkey2(btype2,order2,Pkey2);
            
            E = new LocalRegions::QuadExp(Bkey1,Bkey2,geom);
            
            //----------------------------------------------
            // Define solution to be projected
            Array<OneD, NekDouble> x = Array<OneD, NekDouble>(nq1*nq2);
            Array<OneD, NekDouble> y = Array<OneD, NekDouble>(nq1*nq2);
            E->GetCoords(x,y);
            
            for(i = 0; i < nq1*nq2; ++i)
            {
                sol[i]  = Quad_sol(x[i],y[i],order1,order2,btype1,btype2);
            }
            //---------------------------------------------
        }
        
        break;
    default:
        NEKERROR(ErrorUtil::efatal, "Not a valid 2D expansion.");
    }

    //---------------------------------------------
    // Project onto Expansion
    E->FwdTrans(sol,E->UpdateCoeffs());
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    E->BwdTrans(E->GetCoeffs(),E->UpdatePhys());
    //-------------------------------------------
    
    //--------------------------------------------
    // Write solution
    ofstream outfile("ProjectFile2D.dat");
    E->WriteToFile(outfile,eTecplot);
    outfile.close();
    //-------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error
    cout << "L infinity error: " << E->Linf(sol) << endl;
    cout << "L 2 error:        " << E->L2  (sol) << endl;
    //--------------------------------------------
    
    //-------------------------------------------
    // Evaulate solution at x = y =0  and print error
    Array<OneD, NekDouble> x(2);
    x[0] = (coords[0] + coords[2])*0.5;
    x[1] = (coords[1] + coords[5])*0.5;
    
    if(regionshape == LibUtilities::eTriangle)
    {
        sol[0] = Tri_sol(x[0],x[1],order1,order2);
    }
    else
    {
        sol[0] = Quad_sol(x[0],x[1],order1,order2,btype1,btype2);
    }
    
    NekDouble nsol = E->PhysEvaluate(x);
    cout << "error at x = (" <<x[0] <<","<<x[1] <<"): " << nsol - sol[0] << endl;
    
    return 0;
}


NekDouble Tri_sol(NekDouble x, NekDouble y, int order1, int order2){
    int    l,k;
    NekDouble sol = 0.0;

    for(k = 0; k < order1; ++k)
    {
        for(l = 0; l < order2-k; ++l)
        {
            sol += pow(x,k)*pow(y,l);
        }
    }

    return sol;
}

NekDouble Quad_sol(NekDouble x, NekDouble y, int order1, int order2,
                   LibUtilities::BasisType btype1,
                   LibUtilities::BasisType btype2)
{
    int k,l;
    NekDouble sol = 0.0;

    if(btype1 != LibUtilities::eFourier)
    {
        if(btype2 != LibUtilities::eFourier)
        {
            for(k = 0; k < order1; ++k)
            {
                for(l = 0; l < order2; ++l)
                {
                    sol += pow(x,k)*pow(y,l);
                }
            }
        }
        else
        {
            for(k = 0; k < order1; ++k)
            {
                for(l = 0; l < order2/2; ++l)
                {
                    sol += pow(x,k)*sin(M_PI*l*y) + pow(x,k)*cos(M_PI*l*y);
                }
            }
        }
    }
    else
    {
        if(btype2 != LibUtilities::eFourier)
        {
            for(k = 0; k < order1/2; ++k)
            {
                for(l = 0; l < order2; ++l)
                {
                    sol += sin(M_PI*k*x)*pow(y,l) + cos(M_PI*k*x)*pow(y,l);
                }
            }
        }
        else
        {
            for(k = 0; k < order1/2; ++k)
            {
                for(l = 0; l < order2/2; ++l)
                {
                    sol += sin(M_PI*k*x)*sin(M_PI*l*y)
                        + sin(M_PI*k*x)*cos(M_PI*l*y)
                        + cos(M_PI*k*x)*sin(M_PI*l*y)
                        + cos(M_PI*k*x)*cos(M_PI*l*y);
                }
            }
        }
    }

    return sol;
}
