#include <cstdio>
#include <cstdlib>

#include <LocalRegions/NodalTriExp.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 4 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "--- EXERCISE 2 ----" <<endl;
    cout << "-------------------" <<endl;

    // We will not solve this exercise without using the Nektar++ library

    cout << "METHOD 1: Nektar++" << endl;
    {
        int i,j;
        // Input: specify the details of the expansion
        int nModesDir1 = 9;
        int nModesDir2 = 9;
        int nQuadPointsDir1 = 10;
        int nQuadPointsDir2 = 10;

        // Specify the geometry of the triangular element 
        // step 1: specify the coordinates of the vertices
        Array<OneD, NekDouble> coords(6);
        coords[0]    =   1.0;
        coords[1]    =   0.0;
        coords[2]    =   2.0;
        coords[3]    =   1.0;
        coords[4]    =   1.0;
        coords[5]    =   1.0;   
    
        // step 2: set up the vertices
        SpatialDomains::VertexComponentSharedPtr verts[3];
        const int zero  = 0;
        const int one   = 1;
        const int two   = 2;
        const double dZero = 0.0;
        verts[0] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,zero,coords[0],coords[1],dZero);
        verts[1] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,one,coords[2],coords[3],dZero);
        verts[2] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,two,coords[4],coords[5],dZero);
    
        // step 3: set up the edges
        SpatialDomains::VertexComponentSharedPtr verticescouple1[2] = {verts[0],verts[1]};
        SpatialDomains::VertexComponentSharedPtr verticescouple2[2] = {verts[1],verts[2]};
        SpatialDomains::VertexComponentSharedPtr verticescouple3[2] = {verts[0],verts[2]};
        SpatialDomains::SegGeomSharedPtr edges[3];
        edges[0] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero,  two, verticescouple1);
        edges[1] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(one,   two, verticescouple2);
        edges[2] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(two,   two, verticescouple3);
    
        // step 4: set up the edge orientation
        StdRegions::EdgeOrientation eorient[3];    
        eorient[0] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]);
        eorient[1] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]); 
        eorient[2] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[0]); 
    
        // step 5: set up the geometry object
        int indx = 0;
        SpatialDomains::TriGeomSharedPtr triGeom = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(indx,edges,eorient);
        triGeom->SetOwnData();

        LibUtilities::PointsType nodalPointsType = LibUtilities::eNodalTriElec;
        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eOrtho_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eOrtho_B;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModesDir1,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModesDir2,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // local triangular region
        LocalRegions::NodalTriExpSharedPtr triExpansion = 
            MemoryManager<LocalRegions::NodalTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,nodalPointsType,triGeom);
        
        int nTotModes      = triExpansion->GetNcoeffs();
        int nTotQuadPoints = triExpansion->GetTotPoints();  

        // Evaluate the forcing function at the quadrature points
        Array<OneD, NekDouble> x1(nTotQuadPoints);
        Array<OneD, NekDouble> x2(nTotQuadPoints);
        Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
        triExpansion->GetCoords(x1,x2);

        for(i = 0; i < nTotQuadPoints; i++)
        {
            forcingFunction[i] = pow(x1[i],6)*pow(x2[i],6);
        }     
        
        // Do the projection to obtain the coefficients of the expansion
        // The result is stored in the data member m_coeffs of the NodalTriExp object triExpansion.
        triExpansion->FwdTrans(forcingFunction,triExpansion->UpdateCoeffs());

        // Calculate the exact slution at the nodal points
        Array<OneD, NekDouble> nodalPointsXi1(nTotModes);
        Array<OneD, NekDouble> nodalPointsXi2(nTotModes);
        Array<OneD, NekDouble> nodalPointsX1(nTotModes);
        Array<OneD, NekDouble> nodalPointsX2(nTotModes);
        Array<OneD, NekDouble> exactSolutionAtNodes(nTotModes);
        Array<OneD, NekDouble> coordin(2);
        Array<OneD, NekDouble> coordout(2);
        triExpansion->GetNodalPoints(nodalPointsXi1,nodalPointsXi2);

        for(i = 0; i < nTotModes; i++)
        {
            coordin[0] = nodalPointsXi1[i];
            coordin[1] = nodalPointsXi2[i];

            triExpansion->GetCoord(coordin,coordout);

            nodalPointsX1[i] = coordout[0];
            nodalPointsX2[i] = coordout[1];

            exactSolutionAtNodes[i] =  pow(nodalPointsX1[i],6)*pow(nodalPointsX2[i],6);
        }

        // Display the output
        for(i = 0; i < nTotModes; i++)
        {               
            cout << "Nodal point " << i << ": Result = ";
            cout << (triExpansion->GetCoeffs())[i] << endl;
            cout << "          Exact Result = "<< exactSolutionAtNodes[i] << endl; 
        }
        cout << endl;         
    }
}
