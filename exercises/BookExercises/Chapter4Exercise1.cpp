#include <cstdio>
#include <cstdlib>

#include <StdRegions/StdQuadExp.h>
#include <LocalRegions/QuadExp.h>
#include <StdRegions/StdTriExp.h>
#include <LocalRegions/TriExp.h>
#include <LibUtilities/Polylib/Polylib.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 4 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(a) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            
            // Declare variables to hold the quadrature zeros and weights
            double* quadZerosDir1   = new double[nQuadPointsDir1];
            double* quadWeightsDir1 = new double[nQuadPointsDir1];
            double* quadZerosDir2   = new double[nQuadPointsDir2];
            double* quadWeightsDir2 = new double[nQuadPointsDir2];
            
            // Calculate the GLL-quadrature zeros and weights using the corresponding
            // Polylib routine.
            const double alpha = 0.0;
            const double beta  = 0.0;
            Polylib::zwglj(quadZerosDir1,quadWeightsDir1,nQuadPointsDir1,alpha,beta);
            Polylib::zwglj(quadZerosDir2,quadWeightsDir2,nQuadPointsDir2,alpha,beta);
            
            // Apply the Gaussian quadrature technique
            int i,j;
            double result = 0.0;
            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {
                    result += pow(quadZerosDir1[i],6) * pow(quadZerosDir2[j],6) * 
                        quadWeightsDir1[i] * quadWeightsDir2[j];
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << " (exact solution = ";
            cout << 4.0/49.0 << ")" << endl;
            
            // Deallocate the dynamic memory
            delete[] quadZerosDir1;
            delete[] quadWeightsDir1;
            delete[] quadZerosDir2;
            delete[] quadWeightsDir2;
        }
        cout << endl;
    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    { 
        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;
            
            // Declare variables to hold the quadrature zeros and weights
            Array<OneD, NekDouble> quadZerosDir1;
            Array<OneD, NekDouble> quadWeightsDir1;
            Array<OneD, NekDouble> quadZerosDir2;
            Array<OneD, NekDouble> quadWeightsDir2;
            
            // Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
            // Step 1: Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
            
            // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
            // the PointsManager
            quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
            quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();
            
            quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
            quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
            
            
            // Apply the Gaussian quadrature technique
            int i,j;
            NekDouble result = 0.0;
            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {
                    result += pow(quadZerosDir1[i],6) * pow(quadZerosDir2[j],6) * 
                        quadWeightsDir1[i] * quadWeightsDir2[j];
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << " (exact solution = ";
            cout << 4.0/49.0 << ")" << endl;
        }
        cout << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;
            
            // Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

            // Declare a BasisKey which uniquely defines the (discrete) one-dimensional
            // spectral/hp expansion basis.
            // As we will use the basis only to evaluate an integral, the type and order of the basis
            // can be chosen arbitrarily. 
            const LibUtilities::BasisKey basisKeyDir1(LibUtilities::eModified_A,2,quadPointsKeyDir1);
            const LibUtilities::BasisKey basisKeyDir2(LibUtilities::eModified_A,2,quadPointsKeyDir2);
            
            // Using these keys, now define an spectral/hp expansion of corresponding type on a 
            // standard quadrilateral region
            StdRegions::StdQuadExpSharedPtr quadExpansion = 
                MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);
            
            // Calculate the coordinates xi1 and xi2 of all the quadrature points of this (discrete)
            // expansion
            int nTotQuadPoints = quadExpansion->GetTotPoints();
            Array<OneD,NekDouble> xi1(nTotQuadPoints);
            Array<OneD,NekDouble> xi2(nTotQuadPoints);        
            quadExpansion->GetCoords(xi1,xi2);
            
            // Calculate the values of the function to be integrated 
            int i;
            Array<OneD,NekDouble> integrand(nTotQuadPoints);   
            for(i = 0; i < nTotQuadPoints; i++)
            {
                integrand[i] = pow(xi1[i],6)*pow(xi2[i],6);
            }
            
            // Do the numerical integration by calling the corresponding Nektar++ routine of 
            // the StdExpansion class
            NekDouble result = quadExpansion->Integral(integrand);

            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << " (exact solution = ";
            cout << 4.0/49.0 << ")" << endl;
        }
        cout << endl;
    }

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(b) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            
            // Declare variables to hold the quadrature zeros and weights
            double* quadZerosDir1   = new double[nQuadPointsDir1];
            double* quadWeightsDir1 = new double[nQuadPointsDir1];
            double* quadZerosDir2   = new double[nQuadPointsDir2];
            double* quadWeightsDir2 = new double[nQuadPointsDir2];

            // Calculate the quadrature zeros and weights using the corresponding
            // Polylib routine.
            // In direction 1, GLL points will be used.
            // In direction 2, Gauss-Radau-Jacobi points will be used (this type 
            // absorbes the Jacobian of the collapsed coordinate transformation
            // and does not include the singularity at the collapsed vertex)
            Polylib::zwglj (quadZerosDir1,quadWeightsDir1,nQuadPointsDir1,0.0,0.0);
            Polylib::zwgrjm(quadZerosDir2,quadWeightsDir2,nQuadPointsDir2,1.0,0.0);
            
            // Apply the Gaussian quadrature technique
            int i,j;
            double result = 0.0;
            double xi1;
            double xi2;
            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {
                    // calculate the coordinates xi1 and xi2
                    xi1 = 0.5*(1+quadZerosDir1[i])*(1-quadZerosDir2[j]) - 1.0;
                    xi2 = quadZerosDir2[j];

                    result += pow(xi1,6) * pow(xi2,6) * 
                        quadWeightsDir1[i] * 0.5 * quadWeightsDir2[j];
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << " (exact solution = ";
            cout << 2.0/49.0 << ")" << endl;
            
            // Deallocate the dynamic memory
            delete[] quadZerosDir1;
            delete[] quadWeightsDir1;
            delete[] quadZerosDir2;
            delete[] quadWeightsDir2;
        }
        cout << endl;
    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    { 
        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;
            
            // Declare variables to hold the quadrature zeros and weights
            Array<OneD, NekDouble> quadZerosDir1;
            Array<OneD, NekDouble> quadWeightsDir1;
            Array<OneD, NekDouble> quadZerosDir2;
            Array<OneD, NekDouble> quadWeightsDir2;
            
            // Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
            // Step 1: Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
            
            // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
            // the PointsManager
            quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
            quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();
            
            quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
            quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
            
            
            // Apply the Gaussian quadrature technique
            int i,j;
            NekDouble result = 0.0;
            double xi1;
            double xi2;
            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {
                    // calculate the coordinates xi1 and xi2
                    xi1 = 0.5*(1+quadZerosDir1[i])*(1-quadZerosDir2[j]) - 1.0;
                    xi2 = quadZerosDir2[j];

                    result += pow(xi1,6) * pow(xi2,6) * 
                        quadWeightsDir1[i] * 0.5 * quadWeightsDir2[j];
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << " (exact solution = ";
            cout << 2.0/49.0 << ")" << endl;
        }
        cout << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;
            
            // Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

            // Declare a BasisKey which uniquely defines the (discrete) one-dimensional
            // spectral/hp expansion basis.
            // As we will use the basis only to evaluate an integral, the type and order of the basis
            // can be chosen arbitrarily. 
            const LibUtilities::BasisKey basisKeyDir1(LibUtilities::eModified_A,2,quadPointsKeyDir1);
            const LibUtilities::BasisKey basisKeyDir2(LibUtilities::eModified_A,2,quadPointsKeyDir2);
            
            // Using these keys, now define an spectral/hp expansion of corresponding type on a 
            // standard quadrilateral region
            StdRegions::StdTriExpSharedPtr triExpansion = 
                MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);
            
            // Calculate the coordinates xi1 and xi2 of all the quadrature points of this (discrete)
            // expansion
            int nTotQuadPoints = triExpansion->GetTotPoints();
            Array<OneD,NekDouble> xi1(nTotQuadPoints);
            Array<OneD,NekDouble> xi2(nTotQuadPoints);        
            triExpansion->GetCoords(xi1,xi2);
            
            // Calculate the values of the function to be integrated 
            int i;
            Array<OneD,NekDouble> integrand(nTotQuadPoints);   
            for(i = 0; i < nTotQuadPoints; i++)
            {
                integrand[i] = pow(xi1[i],6)*pow(xi2[i],6);
            }
            
            // Do the numerical integration by calling the corresponding Nektar++ routine of 
            // the StdExpansion class
            NekDouble result = triExpansion->Integral(integrand);

            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << " (exact solution = ";
            cout << 2.0/49.0 << ")" << endl;
        }
        cout << endl;
    }

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(c) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Specify the coordinates of the quadrilateral element
        double x1_A = 0.0;
        double x2_A = 0.0;
        double x1_B = 1.0;
        double x2_B = 0.0;
        double x1_C = 2.0;
        double x2_C = 1.0;
        double x1_D = 0.0;
        double x2_D = 1.0;

        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            
            // Declare variables to hold the quadrature zeros and weights
            double* quadZerosDir1   = new double[nQuadPointsDir1];
            double* quadWeightsDir1 = new double[nQuadPointsDir1];
            double* quadZerosDir2   = new double[nQuadPointsDir2];
            double* quadWeightsDir2 = new double[nQuadPointsDir2];

            // Calculate the GLL-quadrature zeros and weights using the corresponding
            // Polylib routine.
            const double alpha = 0.0;
            const double beta  = 0.0;
            Polylib::zwglj(quadZerosDir1,quadWeightsDir1,nQuadPointsDir1,alpha,beta);
            Polylib::zwglj(quadZerosDir2,quadWeightsDir2,nQuadPointsDir2,alpha,beta);

            // Apply the Gaussian quadrature technique
            int i,j;
            double result = 0.0;
            double x1;
            double x2;
            double jacobian;
            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                { 
                    // Calculate the local coordinates of the quadrature zeros using the 
                    // mapping from reference to local element
                    x1 = x1_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x1_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x1_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                        x1_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                    x2 = x2_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x2_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x2_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                        x2_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                    // Analytically evaluate the Jacobian
                    jacobian = (0.25 * (1-quadZerosDir2[j]) * (x1_B-x1_A) + 0.25 * (1+quadZerosDir2[j]) * (x1_C-x1_D)) *
                        (0.25 * (1-quadZerosDir1[i]) * (x2_D-x2_A) + 0.25 * (1+quadZerosDir1[i]) * (x2_C-x2_B)) - 
                        (0.25 * (1-quadZerosDir1[i]) * (x1_D-x1_A) + 0.25 * (1+quadZerosDir1[i]) * (x1_C-x1_B)) *
                        (0.25 * (1-quadZerosDir2[j]) * (x2_B-x2_A) + 0.25 * (1+quadZerosDir2[j]) * (x2_C-x2_D));
                    
                    result += pow(x1,6) * pow(x2,6) * 
                        quadWeightsDir1[i] * quadWeightsDir2[j] * fabs(jacobian);
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << endl;
            
            // Deallocate the dynamic memory
            delete[] quadZerosDir1;
            delete[] quadWeightsDir1;
            delete[] quadZerosDir2;
            delete[] quadWeightsDir2;
        }
        cout << endl;
    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    { 
        // Specify the coordinates of the quadrilateral element
        double x1_A = 0.0;
        double x2_A = 0.0;
        double x1_B = 1.0;
        double x2_B = 0.0;
        double x1_C = 2.0;
        double x2_C = 1.0;
        double x1_D = 0.0;
        double x2_D = 1.0;

        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;
            
            // Declare variables to hold the quadrature zeros and weights
            Array<OneD, NekDouble> quadZerosDir1;
            Array<OneD, NekDouble> quadWeightsDir1;
            Array<OneD, NekDouble> quadZerosDir2;
            Array<OneD, NekDouble> quadWeightsDir2;
            
            // Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
            // Step 1: Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
            
            // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
            // the PointsManager
            quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
            quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();
            
            quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
            quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
            
            
            // Apply the Gaussian quadrature technique
            int i,j;
            double x1;
            double x2;
            double jacobian;
            NekDouble result = 0.0;
            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {   
                    // Calculate the local coordinates of the quadrature zeros using the 
                    // mapping from reference to local element
                    x1 = x1_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x1_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x1_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                        x1_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                    x2 = x2_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x2_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                        x2_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                        x2_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                    // Analytically evaluate the Jacobian
                    jacobian = (0.25 * (1-quadZerosDir2[j]) * (x1_B-x1_A) + 0.25 * (1+quadZerosDir2[j]) * (x1_C-x1_D)) *
                        (0.25 * (1-quadZerosDir1[i]) * (x2_D-x2_A) + 0.25 * (1+quadZerosDir1[i]) * (x2_C-x2_B)) - 
                        (0.25 * (1-quadZerosDir1[i]) * (x1_D-x1_A) + 0.25 * (1+quadZerosDir1[i]) * (x1_C-x1_B)) *
                        (0.25 * (1-quadZerosDir2[j]) * (x2_B-x2_A) + 0.25 * (1+quadZerosDir2[j]) * (x2_C-x2_D));
                    
                    result += pow(x1,6) * pow(x2,6) * 
                        quadWeightsDir1[i] * quadWeightsDir2[j] * fabs(jacobian);
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << endl;
        }
        cout << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        // Specify the geometry of the quadrilateral element 
        // step 1: specify the coordinates of the vertices
        Array<OneD, NekDouble> coords(8);    
        coords[0]    =   0.0;
        coords[1]    =   0.0;
        coords[2]    =   1.0;
        coords[3]    =   0.0;
        coords[4]    =   2.0;
        coords[5]    =   1.0;
        coords[6]    =   0.0;
        coords[7]    =   1.0;
    
        // step 2: set up the vertices
        SpatialDomains::VertexComponentSharedPtr verts[4];
        const int zero  = 0;
        const int one   = 1;
        const int two   = 2;
        const int three = 3;
        const double dZero = 0.0;
        verts[0] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,zero,coords[0],coords[1],dZero);
        verts[1] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,one,coords[2],coords[3],dZero);
        verts[2] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,two,coords[4],coords[5],dZero);
        verts[3] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,three,coords[6],coords[7],dZero);
    
        // step 3: set up the edges
        SpatialDomains::VertexComponentSharedPtr verticescouple1[2] = {verts[0],verts[1]};
        SpatialDomains::VertexComponentSharedPtr verticescouple2[2] = {verts[1],verts[2]};
        SpatialDomains::VertexComponentSharedPtr verticescouple3[2] = {verts[2],verts[3]};
        SpatialDomains::VertexComponentSharedPtr verticescouple4[2] = {verts[0],verts[3]};
        SpatialDomains::SegGeomSharedPtr edges[4];
        edges[0] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero,  two, verticescouple1);
        edges[1] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(one,   two, verticescouple2);
        edges[2] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(two,   two, verticescouple3);
        edges[3] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(three, two, verticescouple4);
    
        // step 4: set up the edge orientation
        StdRegions::EdgeOrientation eorient[4];    
        eorient[0] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]);
        eorient[1] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]); 
        eorient[2] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[3]); 
        eorient[3] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[3], *edges[0]);
    
        // step 5: set up the geometry object
        int indx;
        SpatialDomains::QuadGeomSharedPtr quadGeom = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(indx,edges,eorient);
        quadGeom->SetOwnData();

        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;
            
            // Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

            // Declare a BasisKey which uniquely defines the (discrete) one-dimensional
            // spectral/hp expansion basis.
            // As we will use the basis only to evaluate an integral, the type and order of the basis
            // can be chosen arbitrarily. 
            const LibUtilities::BasisKey basisKeyDir1(LibUtilities::eModified_A,2,quadPointsKeyDir1);
            const LibUtilities::BasisKey basisKeyDir2(LibUtilities::eModified_A,2,quadPointsKeyDir2);
            
            // Using these keys, now define an spectral/hp expansion of corresponding type on a 
            // standard quadrilateral region
            LocalRegions::QuadExpSharedPtr quadExpansion = 
                MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,quadGeom);
            
            // Calculate the coordinates x1 and x2 of all the quadrature points of this (discrete)
            // expansion
            int nTotQuadPoints = quadExpansion->GetTotPoints();
            Array<OneD,NekDouble> x1(nTotQuadPoints);
            Array<OneD,NekDouble> x2(nTotQuadPoints);        
            quadExpansion->GetCoords(x1,x2);
            
            // Calculate the values of the function to be integrated 
            int i;
            Array<OneD,NekDouble> integrand(nTotQuadPoints);   
            for(i = 0; i < nTotQuadPoints; i++)
            {
                integrand[i] = pow(x1[i],6)*pow(x2[i],6);
            }
            
            // Do the numerical integration by calling the corresponding Nektar++ routine of 
            // the StdExpansion class
            NekDouble result = quadExpansion->Integral(integrand);

            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << endl;
        }
        cout << endl;
    }


    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(d) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Specify the coordinates of the triangular element
        double x1_A = 1.0;
        double x2_A = 0.0;
        double x1_B = 2.0;
        double x2_B = 1.0;
        double x1_C = 1.0;
        double x2_C = 1.0;

        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            
            // Declare variables to hold the quadrature zeros and weights
            double* quadZerosDir1   = new double[nQuadPointsDir1];
            double* quadWeightsDir1 = new double[nQuadPointsDir1];
            double* quadZerosDir2   = new double[nQuadPointsDir2];
            double* quadWeightsDir2 = new double[nQuadPointsDir2];


            // Calculate the quadrature zeros and weights using the corresponding
            // Polylib routine.
            // In direction 1, GLL points will be used.
            // In direction 2, Gauss-Radau-Jacobi points will be used (this type 
            // absorbes the Jacobian of the collapsed coordinate transformation
            // and does not include the singularity at the collapsed vertex)
            Polylib::zwglj (quadZerosDir1,quadWeightsDir1,nQuadPointsDir1,0.0,0.0);
            Polylib::zwgrjm(quadZerosDir2,quadWeightsDir2,nQuadPointsDir2,1.0,0.0);

            // Apply the Gaussian quadrature technique
            int i,j;
            double result = 0.0;
            double xi1;
            double xi2;
            double x1;
            double x2;
            double jacobian;

            // Analytically evaluate the Jacobian.
            // For a straight-sided triangle, the Jacobian is constant
            jacobian = 0.25 * (x1_B - x1_A) * (x2_C - x2_A) + 
                0.25 * (x1_C - x1_A) * (x2_B - x2_A);

            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                { 
                    // calculate the coordinates xi1 and xi2
                    xi1 = 0.5*(1+quadZerosDir1[i])*(1-quadZerosDir2[j]) - 1.0;
                    xi2 = quadZerosDir2[j];

                    // Calculate the local coordinates of the quadrature zeros using the 
                    // mapping from reference to local element
                    x1 = x1_A * 0.5 * (-xi2 - xi1) + 
                        x1_B * 0.5 * (1 + xi1) + 
                        x1_C * 0.5 * (1 + xi2); 

                    x2 = x2_A * 0.5 * (-xi2 - xi1) + 
                        x2_B * 0.5 * (1 + xi1) + 
                        x2_C * 0.5 * (1 + xi2); 
                    
                    result += pow(x1,6) * pow(x2,6) * 
                        quadWeightsDir1[i] * 0.5 * quadWeightsDir2[j] * fabs(jacobian);
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << endl;
            
            // Deallocate the dynamic memory
            delete[] quadZerosDir1;
            delete[] quadWeightsDir1;
            delete[] quadZerosDir2;
            delete[] quadWeightsDir2;
        }
        cout << endl;
    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    { 
        // Specify the coordinates of the quadrilateral element
        double x1_A = 1.0;
        double x2_A = 0.0;
        double x1_B = 2.0;
        double x2_B = 1.0;
        double x1_C = 1.0;
        double x2_C = 1.0;

        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;
            
            // Declare variables to hold the quadrature zeros and weights
            Array<OneD, NekDouble> quadZerosDir1;
            Array<OneD, NekDouble> quadWeightsDir1;
            Array<OneD, NekDouble> quadZerosDir2;
            Array<OneD, NekDouble> quadWeightsDir2;
            
            // Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
            // Step 1: Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
            
            // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
            // the PointsManager
            quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
            quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();
            
            quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
            quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
            
            
            // Apply the Gaussian quadrature technique
            int i,j;
            NekDouble result = 0.0;
            double xi1;
            double xi2;
            double x1;
            double x2;
            double jacobian;

            // Analytically evaluate the Jacobian.
            // For a straight-sided triangle, the Jacobian is constant
            jacobian = 0.25 * (x1_B - x1_A) * (x2_C - x2_A) + 
                0.25 * (x1_C - x1_A) * (x2_B - x2_A);

            for(i = 0; i < nQuadPointsDir1; i++)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {   
                    // calculate the coordinates xi1 and xi2
                    xi1 = 0.5*(1+quadZerosDir1[i])*(1-quadZerosDir2[j]) - 1.0;
                    xi2 = quadZerosDir2[j];

                    // Calculate the local coordinates of the quadrature zeros using the 
                    // mapping from reference to local element
                    x1 = x1_A * 0.5 * (-xi2 - xi1) + 
                        x1_B * 0.5 * (1 + xi1) + 
                        x1_C * 0.5 * (1 + xi2); 

                    x2 = x2_A * 0.5 * (-xi2 - xi1) + 
                        x2_B * 0.5 * (1 + xi1) + 
                        x2_C * 0.5 * (1 + xi2);  
                    
                    result += pow(x1,6) * pow(x2,6) * 
                        quadWeightsDir1[i] * 0.5 * quadWeightsDir2[j] * fabs(jacobian);
                }
            }
            
            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << endl;
        }
        cout << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
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
        int indx=0;
        SpatialDomains::TriGeomSharedPtr triGeom = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(indx,edges,eorient);
        triGeom->SetOwnData();

        for(int Q = 4; Q < 7; Q++)
        {
            // Input: specify the number and type of quadrature points in both directions
            int nQuadPointsDir1 = Q;
            int nQuadPointsDir2 = Q;
            LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
            LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;
            
            // Declare a PointsKey which uniquely defines the quadrature points
            const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
            const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

            // Declare a BasisKey which uniquely defines the (discrete) one-dimensional
            // spectral/hp expansion basis.
            // As we will use the basis only to evaluate an integral, the type and order of the basis
            // can be chosen arbitrarily. 
            const LibUtilities::BasisKey basisKeyDir1(LibUtilities::eModified_A,2,quadPointsKeyDir1);
            const LibUtilities::BasisKey basisKeyDir2(LibUtilities::eModified_A,2,quadPointsKeyDir2);
            
            // Using these keys, now define an spectral/hp expansion of corresponding type on a 
            // standard quadrilateral region
            LocalRegions::TriExpSharedPtr triExpansion = 
                MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,triGeom);
            
            // Calculate the coordinates x1 and x2 of all the quadrature points of this (discrete)
            // expansion
            int nTotQuadPoints = triExpansion->GetTotPoints();
            Array<OneD,NekDouble> x1(nTotQuadPoints);
            Array<OneD,NekDouble> x2(nTotQuadPoints);        
            triExpansion->GetCoords(x1,x2);
            
            // Calculate the values of the function to be integrated 
            int i;
            Array<OneD,NekDouble> integrand(nTotQuadPoints);   
            for(i = 0; i < nTotQuadPoints; i++)
            {
                integrand[i] = pow(x1[i],6)*pow(x2[i],6);
            }

            
            // Do the numerical integration by calling the corresponding Nektar++ routine of 
            // the StdExpansion class
            NekDouble result = triExpansion->Integral(integrand);

            // Display the output
            cout << "Q = " << Q << " => Result = "<< result << endl;
        }
        cout << endl;
    }

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(e) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        cout << "Part 1: Quadrilateral element" << endl;
        {
            // Specify the coordinates of the quadrilateral element
            double x1_A = 0.0;
            double x2_A = 0.0;
            double x1_B = 1.0;
            double x2_B = 0.0;
            double x1_C = 2.0;
            double x2_C = 1.0;
            double x1_D = 0.0;
            double x2_D = 1.0;

            // Calculate the exact solution
            double exactResult = (1-cos(1.))*(1-cos(1.)) + 
                0.25*cos(1.)*(5*cos(2.)-cos(4.)-4*cos(1.)) + 
                0.25*sin(1.)*(2+sin(2.)-sin(4.));

            for(int Q = 2; Q < 9; Q++)
            {
                // Input: specify the number of quadrature points in both directions
                int nQuadPointsDir1 = Q;
                int nQuadPointsDir2 = Q;
            
                // Declare variables to hold the quadrature zeros and weights
                double* quadZerosDir1   = new double[nQuadPointsDir1];
                double* quadWeightsDir1 = new double[nQuadPointsDir1];
                double* quadZerosDir2   = new double[nQuadPointsDir2];
                double* quadWeightsDir2 = new double[nQuadPointsDir2];

                // Calculate the GLL-quadrature zeros and weights using the corresponding
                // Polylib routine.
                const double alpha = 0.0;
                const double beta  = 0.0;
                Polylib::zwglj(quadZerosDir1,quadWeightsDir1,nQuadPointsDir1,alpha,beta);
                Polylib::zwglj(quadZerosDir2,quadWeightsDir2,nQuadPointsDir2,alpha,beta);

                // Apply the Gaussian quadrature technique
                int i,j;
                double result = 0.0;
                double x1;
                double x2;
                double jacobian;
                for(i = 0; i < nQuadPointsDir1; i++)
                {
                    for(j = 0; j < nQuadPointsDir2; j++)
                    { 
                        // Calculate the local coordinates of the quadrature zeros using the 
                        // mapping from reference to local element
                        x1 = x1_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x1_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x1_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                            x1_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                        x2 = x2_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x2_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x2_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                            x2_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                        // Analytically evaluate the Jacobian
                        jacobian = (0.25 * (1-quadZerosDir2[j]) * (x1_B-x1_A) + 0.25 * (1+quadZerosDir2[j]) * (x1_C-x1_D)) *
                            (0.25 * (1-quadZerosDir1[i]) * (x2_D-x2_A) + 0.25 * (1+quadZerosDir1[i]) * (x2_C-x2_B)) - 
                            (0.25 * (1-quadZerosDir1[i]) * (x1_D-x1_A) + 0.25 * (1+quadZerosDir1[i]) * (x1_C-x1_B)) *
                            (0.25 * (1-quadZerosDir2[j]) * (x2_B-x2_A) + 0.25 * (1+quadZerosDir2[j]) * (x2_C-x2_D));
                    
                        result += sin(x1) * sin(x2) * 
                            quadWeightsDir1[i] * quadWeightsDir2[j] * fabs(jacobian);
                    }
                }

                // Display the output
                cout << "Q = " << Q << " => error = "<< fabs(result-exactResult) << endl;
            
                // Deallocate the dynamic memory
                delete[] quadZerosDir1;
                delete[] quadWeightsDir1;
                delete[] quadZerosDir2;
                delete[] quadWeightsDir2;
            }
            cout << endl;
        }

        cout << "Part 2: Triangular element" << endl;
        {
            // Specify the coordinates of the triangular element
            double x1_A = 1.0;
            double x2_A = 0.0;
            double x1_B = 2.0;
            double x2_B = 1.0;
            double x1_C = 1.0;
            double x2_C = 1.0;

            // Calculate the exact solution
            double exactResult =  0.25*cos(1.)*(5*cos(2.)-cos(4.)-4*cos(1.)) + 
                0.25*sin(1.)*(2+sin(2.)-sin(4.));

            for(int Q = 2; Q < 9; Q++)
            {
                // Input: specify the number of quadrature points in both directions
                int nQuadPointsDir1 = Q;
                int nQuadPointsDir2 = Q;
            
                // Declare variables to hold the quadrature zeros and weights
                double* quadZerosDir1   = new double[nQuadPointsDir1];
                double* quadWeightsDir1 = new double[nQuadPointsDir1];
                double* quadZerosDir2   = new double[nQuadPointsDir2];
                double* quadWeightsDir2 = new double[nQuadPointsDir2];


                // Calculate the quadrature zeros and weights using the corresponding
                // Polylib routine.
                // In direction 1, GLL points will be used.
                // In direction 2, Gauss-Radau-Jacobi points will be used (this type 
                // absorbes the Jacobian of the collapsed coordinate transformation
                // and does not include the singularity at the collapsed vertex)
                Polylib::zwglj (quadZerosDir1,quadWeightsDir1,nQuadPointsDir1,0.0,0.0);
                Polylib::zwgrjm(quadZerosDir2,quadWeightsDir2,nQuadPointsDir2,1.0,0.0);

                // Apply the Gaussian quadrature technique
                int i,j;
                double result = 0.0;
                double xi1;
                double xi2;
                double x1;
                double x2;
                double jacobian;

                // Analytically evaluate the Jacobian.
                // For a straight-sided triangle, the Jacobian is constant
                jacobian = 0.25 * (x1_B - x1_A) * (x2_C - x2_A) + 
                    0.25 * (x1_C - x1_A) * (x2_B - x2_A);

                for(i = 0; i < nQuadPointsDir1; i++)
                {
                    for(j = 0; j < nQuadPointsDir2; j++)
                    { 
                        // calculate the coordinates xi1 and xi2
                        xi1 = 0.5*(1+quadZerosDir1[i])*(1-quadZerosDir2[j]) - 1.0;
                        xi2 = quadZerosDir2[j];

                        // Calculate the local coordinates of the quadrature zeros using the 
                        // mapping from reference to local element
                        x1 = x1_A * 0.5 * (-xi2 - xi1) + 
                            x1_B * 0.5 * (1 + xi1) + 
                            x1_C * 0.5 * (1 + xi2); 

                        x2 = x2_A * 0.5 * (-xi2 - xi1) + 
                            x2_B * 0.5 * (1 + xi1) + 
                            x2_C * 0.5 * (1 + xi2); 
                    
                        result += sin(x1) * sin(x2) * 
                            quadWeightsDir1[i] * 0.5 * quadWeightsDir2[j] * fabs(jacobian);
                    }
                }
            
                // Display the output
                cout << "Q = " << Q << " => error = "<< fabs(result-exactResult) << endl;
            
                // Deallocate the dynamic memory
                delete[] quadZerosDir1;
                delete[] quadWeightsDir1;
                delete[] quadZerosDir2;
                delete[] quadWeightsDir2;
            }
            cout << endl;
        }
    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    { 
        cout << "Part 1: Quadrilateral element" << endl;
        {
            // Specify the coordinates of the quadrilateral element
            double x1_A = 0.0;
            double x2_A = 0.0;
            double x1_B = 1.0;
            double x2_B = 0.0;
            double x1_C = 2.0;
            double x2_C = 1.0;
            double x1_D = 0.0;
            double x2_D = 1.0;
        
            // Calculate the exact solution
            double exactResult = (1-cos(1.))*(1-cos(1.)) + 
                0.25*cos(1.)*(5*cos(2.)-cos(4.)-4*cos(1.)) + 
                0.25*sin(1.)*(2+sin(2.)-sin(4.));
        
            for(int Q = 2; Q < 9; Q++)
            {
                // Input: specify the number and type of quadrature points in both directions
                int nQuadPointsDir1 = Q;
                int nQuadPointsDir2 = Q;
                LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
                LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;
            
                // Declare variables to hold the quadrature zeros and weights
                Array<OneD, NekDouble> quadZerosDir1;
                Array<OneD, NekDouble> quadWeightsDir1;
                Array<OneD, NekDouble> quadZerosDir2;
                Array<OneD, NekDouble> quadWeightsDir2;
            
                // Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
                // Step 1: Declare a PointsKey which uniquely defines the quadrature points
                const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
                const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
            
                // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
                // the PointsManager
                quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
                quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();
            
                quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
                quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
            
            
                // Apply the Gaussian quadrature technique
                int i,j;
                double x1;
                double x2;
                double jacobian;
                NekDouble result = 0.0;
                for(i = 0; i < nQuadPointsDir1; i++)
                {
                    for(j = 0; j < nQuadPointsDir2; j++)
                    {   
                        // Calculate the local coordinates of the quadrature zeros using the 
                        // mapping from reference to local element
                        x1 = x1_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x1_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x1_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                            x1_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                        x2 = x2_A * 0.25 * (1-quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x2_B * 0.25 * (1+quadZerosDir1[i]) * (1-quadZerosDir2[j]) + 
                            x2_D * 0.25 * (1-quadZerosDir1[i]) * (1+quadZerosDir2[j]) + 
                            x2_C * 0.25 * (1+quadZerosDir1[i]) * (1+quadZerosDir2[j]); 
                    
                        // Analytically evaluate the Jacobian
                        jacobian = (0.25 * (1-quadZerosDir2[j]) * (x1_B-x1_A) + 0.25 * (1+quadZerosDir2[j]) * (x1_C-x1_D)) *
                            (0.25 * (1-quadZerosDir1[i]) * (x2_D-x2_A) + 0.25 * (1+quadZerosDir1[i]) * (x2_C-x2_B)) - 
                            (0.25 * (1-quadZerosDir1[i]) * (x1_D-x1_A) + 0.25 * (1+quadZerosDir1[i]) * (x1_C-x1_B)) *
                            (0.25 * (1-quadZerosDir2[j]) * (x2_B-x2_A) + 0.25 * (1+quadZerosDir2[j]) * (x2_C-x2_D));
                    
                        result += sin(x1) * sin(x2) * 
                            quadWeightsDir1[i] * quadWeightsDir2[j] * fabs(jacobian);
                    }
                }
            
                // Display the output
                cout << "Q = " << Q << " => error = "<< fabs(result-exactResult) << endl;
            }
            cout << endl;
        }

        cout << "Part 2: Triangular element" << endl;
        {
            // Specify the coordinates of the quadrilateral element
            double x1_A = 1.0;
            double x2_A = 0.0;
            double x1_B = 2.0;
            double x2_B = 1.0;
            double x1_C = 1.0;
            double x2_C = 1.0;

            // Calculate the exact solution
            double exactResult =  0.25*cos(1.)*(5*cos(2.)-cos(4.)-4*cos(1.)) + 
                0.25*sin(1.)*(2+sin(2.)-sin(4.));

            for(int Q = 2; Q < 9; Q++)
            {
                // Input: specify the number and type of quadrature points in both directions
                int nQuadPointsDir1 = Q;
                int nQuadPointsDir2 = Q;
                LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
                LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;
            
                // Declare variables to hold the quadrature zeros and weights
                Array<OneD, NekDouble> quadZerosDir1;
                Array<OneD, NekDouble> quadWeightsDir1;
                Array<OneD, NekDouble> quadZerosDir2;
                Array<OneD, NekDouble> quadWeightsDir2;
            
                // Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
                // Step 1: Declare a PointsKey which uniquely defines the quadrature points
                const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
                const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);
            
                // Step 2: Using this key, the quadrature zeros and weights can now be retrieved through
                // the PointsManager
                quadZerosDir1   = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetZ();
                quadWeightsDir1 = LibUtilities::PointsManager()[quadPointsKeyDir1]->GetW();
            
                quadZerosDir2   = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetZ();
                quadWeightsDir2 = LibUtilities::PointsManager()[quadPointsKeyDir2]->GetW();
            
            
                // Apply the Gaussian quadrature technique
                int i,j;
                NekDouble result = 0.0;
                double xi1;
                double xi2;
                double x1;
                double x2;
                double jacobian;

                // Analytically evaluate the Jacobian.
                // For a straight-sided triangle, the Jacobian is constant
                jacobian = 0.25 * (x1_B - x1_A) * (x2_C - x2_A) + 
                    0.25 * (x1_C - x1_A) * (x2_B - x2_A);

                for(i = 0; i < nQuadPointsDir1; i++)
                {
                    for(j = 0; j < nQuadPointsDir2; j++)
                    {   
                        // calculate the coordinates xi1 and xi2
                        xi1 = 0.5*(1+quadZerosDir1[i])*(1-quadZerosDir2[j]) - 1.0;
                        xi2 = quadZerosDir2[j];

                        // Calculate the local coordinates of the quadrature zeros using the 
                        // mapping from reference to local element
                        x1 = x1_A * 0.5 * (-xi2 - xi1) + 
                            x1_B * 0.5 * (1 + xi1) + 
                            x1_C * 0.5 * (1 + xi2); 

                        x2 = x2_A * 0.5 * (-xi2 - xi1) + 
                            x2_B * 0.5 * (1 + xi1) + 
                            x2_C * 0.5 * (1 + xi2);  
                    
                        result += sin(x1) * sin(x2) * 
                            quadWeightsDir1[i] * 0.5 * quadWeightsDir2[j] * fabs(jacobian);
                    }
                }
            
                // Display the output
                cout << "Q = " << Q << " => error = "<< fabs(result-exactResult) << endl;
            }
            cout << endl;
        }
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        cout << "Part 1: Quadrilateral element" << endl;
        {
            // Specify the geometry of the quadrilateral element 
            // step 1: specify the coordinates of the vertices
            Array<OneD, NekDouble> coords(8);    
            coords[0]    =   0.0;
            coords[1]    =   0.0;
            coords[2]    =   1.0;
            coords[3]    =   0.0;
            coords[4]    =   2.0;
            coords[5]    =   1.0;
            coords[6]    =   0.0;
            coords[7]    =   1.0;
    
            // step 2: set up the vertices
            SpatialDomains::VertexComponentSharedPtr verts[4];
            const int zero  = 0;
            const int one   = 1;
            const int two   = 2;
            const int three = 3;
            const double dZero = 0.0;
            verts[0] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,zero,coords[0],coords[1],dZero);
            verts[1] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,one,coords[2],coords[3],dZero);
            verts[2] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,two,coords[4],coords[5],dZero);
            verts[3] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(two,three,coords[6],coords[7],dZero);
    
            // step 3: set up the edges
            SpatialDomains::VertexComponentSharedPtr verticescouple1[2] = {verts[0],verts[1]};
            SpatialDomains::VertexComponentSharedPtr verticescouple2[2] = {verts[1],verts[2]};
            SpatialDomains::VertexComponentSharedPtr verticescouple3[2] = {verts[2],verts[3]};
            SpatialDomains::VertexComponentSharedPtr verticescouple4[2] = {verts[0],verts[3]};
            SpatialDomains::SegGeomSharedPtr edges[4];
            edges[0] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero,  two, verticescouple1);
            edges[1] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(one,   two, verticescouple2);
            edges[2] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(two,   two, verticescouple3);
            edges[3] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(three, two, verticescouple4);
    
            // step 4: set up the edge orientation
            StdRegions::EdgeOrientation eorient[4];    
            eorient[0] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]);
            eorient[1] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]); 
            eorient[2] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[3]); 
            eorient[3] = SpatialDomains::SegGeom::GetEdgeOrientation(*edges[3], *edges[0]);
    
            // step 5: set up the geometry object
            int indx=0;
            SpatialDomains::QuadGeomSharedPtr quadGeom = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(indx,edges,eorient);
            quadGeom->SetOwnData();

            // Calculate the exact solution
            double exactResult = (1-cos(1.))*(1-cos(1.)) + 
                0.25*cos(1.)*(5*cos(2.)-cos(4.)-4*cos(1.)) + 
                0.25*sin(1.)*(2+sin(2.)-sin(4.));

            for(int Q = 2; Q < 9; Q++)
            {
                // Input: specify the number and type of quadrature points in both directions
                int nQuadPointsDir1 = Q;
                int nQuadPointsDir2 = Q;
                LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
                LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;
            
                // Declare a PointsKey which uniquely defines the quadrature points
                const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
                const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

                // Declare a BasisKey which uniquely defines the (discrete) one-dimensional
                // spectral/hp expansion basis.
                // As we will use the basis only to evaluate an integral, the type and order of the basis
                // can be chosen arbitrarily. 
                const LibUtilities::BasisKey basisKeyDir1(LibUtilities::eModified_A,2,quadPointsKeyDir1);
                const LibUtilities::BasisKey basisKeyDir2(LibUtilities::eModified_A,2,quadPointsKeyDir2);
            
                // Using these keys, now define an spectral/hp expansion of corresponding type on a 
                // standard quadrilateral region
                LocalRegions::QuadExpSharedPtr quadExpansion = 
                    MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,quadGeom);
            
                // Calculate the coordinates x1 and x2 of all the quadrature points of this (discrete)
                // expansion
                int nTotQuadPoints = quadExpansion->GetTotPoints();
                Array<OneD,NekDouble> x1(nTotQuadPoints);
                Array<OneD,NekDouble> x2(nTotQuadPoints);        
                quadExpansion->GetCoords(x1,x2);
            
                // Calculate the values of the function to be integrated 
                int i;
                Array<OneD,NekDouble> integrand(nTotQuadPoints);   
                for(i = 0; i < nTotQuadPoints; i++)
                {
                    integrand[i] = sin(x1[i])*sin(x2[i]);
                }
            
                // Do the numerical integration by calling the corresponding Nektar++ routine of 
                // the StdExpansion class
                NekDouble result = quadExpansion->Integral(integrand);

                // Display the output
                cout << "Q = " << Q << " => error = "<< fabs(result-exactResult) << endl;
            }
            cout << endl;
        }     
        cout << "Part 2: Triangular element" << endl;
        {
            // Specify the geometry of the quadrilateral element 
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

            // Calculate the exact solution
            double exactResult =  0.25*cos(1.)*(5*cos(2.)-cos(4.)-4*cos(1.)) + 
                0.25*sin(1.)*(2+sin(2.)-sin(4.));

            for(int Q = 2; Q < 9; Q++)
            {
                // Input: specify the number and type of quadrature points in both directions
                int nQuadPointsDir1 = Q;
                int nQuadPointsDir2 = Q;
                LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
                LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;
            
                // Declare a PointsKey which uniquely defines the quadrature points
                const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
                const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

                // Declare a BasisKey which uniquely defines the (discrete) one-dimensional
                // spectral/hp expansion basis.
                // As we will use the basis only to evaluate an integral, the type and order of the basis
                // can be chosen arbitrarily. 
                const LibUtilities::BasisKey basisKeyDir1(LibUtilities::eModified_A,2,quadPointsKeyDir1);
                const LibUtilities::BasisKey basisKeyDir2(LibUtilities::eModified_A,2,quadPointsKeyDir2);
            
                // Using these keys, now define an spectral/hp expansion of corresponding type on a 
                // standard quadrilateral region
                LocalRegions::TriExpSharedPtr triExpansion = 
                    MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,triGeom);
            
                // Calculate the coordinates x1 and x2 of all the quadrature points of this (discrete)
                // expansion
                int nTotQuadPoints = triExpansion->GetTotPoints();
                Array<OneD,NekDouble> x1(nTotQuadPoints);
                Array<OneD,NekDouble> x2(nTotQuadPoints);        
                triExpansion->GetCoords(x1,x2);
            
                // Calculate the values of the function to be integrated 
                int i;
                Array<OneD,NekDouble> integrand(nTotQuadPoints);   
                for(i = 0; i < nTotQuadPoints; i++)
                {
                    integrand[i] = sin(x1[i])*sin(x2[i]);
                }

            
                // Do the numerical integration by calling the corresponding Nektar++ routine of 
                // the StdExpansion class
                NekDouble result = triExpansion->Integral(integrand);

                // Display the output
                cout << "Q = " << Q << " => error = "<< fabs(result-exactResult) << endl;
            }
            cout << endl;
        }
    }

}

