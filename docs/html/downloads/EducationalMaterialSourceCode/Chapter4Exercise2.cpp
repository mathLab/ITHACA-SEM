#include <cstdio>
#include <cstdlib>

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
    cout << "--- EXERCISE 2 ----" <<endl;
    cout << "-------------------" <<endl;

    // Exercise 2(a), 2(b) and 2(c) will be solved simultaneously
    // It will only be solved for the triangular case with an C0 
    // continuous modal expansion.

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Specify the coordinates of the triangular element
        double x1_A = 1.0;
        double x2_A = 0.0;
        double x1_B = 2.0;
        double x2_B = 1.0;
        double x1_C = 1.0;
        double x2_C = 1.0;
        
        // Input: specify the details of the expansion
        int nModesDir1      =  9;
        int nModesDir2      =  9;
        int nQuadPointsDir1 =  10;
        int nQuadPointsDir2 =  10;
        int nTotModes       =  nModesDir1*(nModesDir1+1)/2 + nModesDir1*(nModesDir2-nModesDir1);
        int nTotQuadPoints  =  nQuadPointsDir1*nQuadPointsDir1;     

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

        // Allocate memory for the matrices which will hold the values of the
        // (one-dimensional) bases at the quadrature points
        double**  base1 = new double*  [nModesDir1];
        double*** base2 = new double** [nModesDir1];

        int i,j,k,p,q;
        for(p = 0; p < nModesDir1; p++)
        {
            base1[p] = new double[nQuadPointsDir1];
        }        
        for(p = 0; p < nModesDir1; p++)
        {
            base2[p] = new double* [nModesDir2];
            
            for(q = 0; q < nModesDir2; q++)
            {
                base2[p][q] = new double [nQuadPointsDir2];
            }
        }

        // Fill the matrix base0 with the Modal C0 basis \psi_i
        double* JacobiPolDir1 = new double [nQuadPointsDir1];
        double* JacobiPolDir2 = new double [nQuadPointsDir2];
        for(i = 0; i < nModesDir1; i++)
        {
            for(j = 0; j < nQuadPointsDir1; j++)
            {
                if(i==0)
                {
                    base1[i][j] = 0.5 * (1.0 - quadZerosDir1[j]);
                }
                else if(i==nModesDir1-1)
                {
                    base1[i][j] = 0.5 * (1.0 + quadZerosDir1[j]);
                }
                else
                {
                    // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
                    Polylib::jacobfd(nQuadPointsDir1, quadZerosDir1, JacobiPolDir1, NULL, i-1, 1.0, 1.0);
                    base1[i][j] = 0.5 * (1.0 - quadZerosDir1[j]) * 0.5 * (1.0 + quadZerosDir1[j]) * JacobiPolDir1[j];
                }
            }
        }

        // Fill the matrix base1 with the Modal C0 basis \psi_{ij}
        // Note that were are storing the values for the entire index
        // range 0 <= i <= P_1 ; 0 <= j <= P_2.
        // Although not entries are required, this facilitates the implementation    
        for(i = 0; i < nModesDir1; i++)
        {
            if( (i==0) || (i==nModesDir1-1))
            {
                for(j = 0; j < nModesDir2; j++)
                {                    
                    if(j==0)
                    {
                        for(k = 0; k < nQuadPointsDir2; k++)
                        {
                            base2[i][j][k] = 0.5 * (1.0 - quadZerosDir2[k]);
                        }
                    }
                    else if(j==nModesDir2-1)
                    {
                        for(k = 0; k < nQuadPointsDir2; k++)
                        {
                            base2[i][j][k] = 0.5 * (1.0 + quadZerosDir2[k]);
                        }
                    }
                    else
                    {                
                        // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
                        Polylib::jacobfd(nQuadPointsDir2, quadZerosDir2, JacobiPolDir2, NULL, j-1, 1.0, 1.0);

                        for(k = 0; k < nQuadPointsDir2; k++)
                        {
                            base2[i][j][k] = 0.5 * (1.0 - quadZerosDir2[k]) * 0.5 * (1.0 + quadZerosDir2[k]) * JacobiPolDir2[k];
                        }
                    }
                }
            } 
            else
            {   
                for(j = 0; j < nModesDir2; j++)
                {
                    if(j == 0)
                    {
                        for(k = 0; k < nQuadPointsDir2; k++)
                        {
                            base2[i][j][k] = pow( 0.5 * (1.0 - quadZerosDir2[k]) , i+1);
                        }
                    }
                    else
                    {
                        // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
                        Polylib::jacobfd(nQuadPointsDir2, quadZerosDir2, JacobiPolDir2, NULL, j-1, 2.0*i + 1.0, 1.0);

                        for(k = 0; k < nQuadPointsDir2; k++)
                        {
                            base2[i][j][k] = pow( 0.5 * (1.0 - quadZerosDir2[k]) , i+1) * 
                                0.5 * (1.0 + quadZerosDir2[k]) * JacobiPolDir2[k];
                        }
                    }
                }            
            }
        }

        // Set up the index map n(p,q)
        // This is done by listing the vertex modes first,
        // then the edge modes and finally the interior modes
        int cnt = 0;
        int** indexmap = new int* [nModesDir1];
        for(i = 0; i < nModesDir1; i++)
        {
            indexmap[i] = new int [nModesDir2];

            // Initialise the indexmap to -1. Doing so, a value of -1
            // for n[p][q] will later indicate that this particular
            // combination of p and q does not correspond to an actual 
            // expansion mode
            for(j = 0; j < nModesDir2; j++)
            {
                indexmap[i][j] = -1;
            }
        }

        indexmap[0][0]                       = 0; // Vertex A
        indexmap[nModesDir1-1][0]            = 1; // Vertex B
        indexmap[0][nModesDir2-1]            = 2; // Vertex C
        indexmap[nModesDir1-1][nModesDir2-1] = 2; // Vertex D (will be collapsed with index C)
        cnt = 3;
        for(p = 1; p < nModesDir1-1; p++)  // Edge AB
        {
            indexmap[p][0] = cnt++;
        }        
        for(q = 1; q < nModesDir2-1; q++)  // Egde AC
        {
            indexmap[0][q] = cnt++;
        }    
        for(q = 1; q < nModesDir2-1; q++)  // Egde BD
        {
            indexmap[nModesDir1-1][q] = cnt++;
        }
        for(p = 1;  p < nModesDir1-1; p++)  // Interior modes
        {
            for(q = 1; q < nModesDir2-1-p; q++)
            {
                indexmap[p][q] = cnt++;
            }
        }
        
        // Analytically evaluate the Jacobian.
        // For a straight-sided triangle, the Jacobian is constant
        double jacobian;
        jacobian = 0.25 * (x1_B - x1_A) * (x2_C - x2_A) + 
            0.25 * (x1_C - x1_A) * (x2_B - x2_A);
        
        // Allocate memory for the mass matrix
        double* M = new double [nTotModes*nTotModes];
        for(i = 0; i < nTotModes*nTotModes; i++)
        {
            // Initialise the mass matrix with zeros
            M[i] = 0.0;
        }  

        // STEP 1: Calculate the mass matrix
        int r,s,m,n;
        double resultDir1;
        double resultDir2;
        for(p = 0; p < nModesDir1; p++)
        { 
            for(q = 0; q < nModesDir2; q++)
            {
                for(r = 0; r < nModesDir1; r++)
                { 
                    for(s = 0; s < nModesDir2; s++)
                    {
                        n = indexmap[p][q];
                        m = indexmap[r][s];

                        // Only if n and m are different from -1, we are dealing with a proper
                        // entrt of the mass matrix
                        if( (n!=-1) && (m!=-1) )
                        {
                            resultDir1 = 0.0;
                            resultDir2 = 0.0;
                            
                            // Integrate in direction 1
                            for(i = 0; i < nQuadPointsDir1; i++)
                            {
                                resultDir1 += base1[p][i] * base1[r][i] * quadWeightsDir1[i];
                            }
                            // Integrate in direction 2
                            for(i = 0; i < nQuadPointsDir2; i++)
                            {
                                resultDir2 += base2[p][q][i] * base2[r][s][i] * 0.5 * quadWeightsDir2[i];
                            }
                                        
                            // We add the result tot the matrix entry as this can account
                            // for the multiple contributions of the collapsed vertex                
                            M[n+m*nTotModes] += fabs(jacobian)*resultDir1*resultDir2;
                        }
                    }
                }
            }
        }

        // STEP 2: calculate the right hand side vector f
        // Allocate memory for the right hand side vector f
        double* f = new double [nTotModes];
        for(i = 0; i < nTotModes; i++)
        {
            // Initialise f with zeros
            f[i] = 0.0;
        }  

        double xi1;
        double xi2;
        double x1;
        double x2;
        // Calculate the rhs
        // Note that we are not using the sum-factorisation technique here
        for(p = 0; p < nModesDir1; p++)
        { 
            for(q = 0; q < nModesDir2; q++)
            {               
                n = indexmap[p][q];
                
                // Only if n is different from -1, we are dealing with a proper
                // entrt of the mass matrix
                if( (n!=-1) )
                {                    
                    // Evaluate the inner product
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

                            f[n] += fabs(jacobian) * pow(x1,6) * pow(x2,6) * 
                                base1[p][i] * quadWeightsDir1[i] *
                                base2[p][q][j] * 0.5 * quadWeightsDir2[j];
                        }
                    }
                }
            }
        }

        // STEP 3: solve the system
        // Factorise the mass matrix
        int*    pivot = new int[nTotModes];
        double* work  = new double[nTotModes];
        int     info  = 0;
                    
        Lapack::Dgetrf(nTotModes, nTotModes, M, nTotModes, pivot, info);
        
        if( info < 0 )
        {
            cerr << "ERROR: The " << -info << "th parameter had an illegal parameter for dgetrf" << endl;
            exit(1);
        }
        else if( info > 0 )
        {
            cerr << "ERROR: Element M[" << info << "] is 0 from dgetrf"<<endl;
            exit(1);
        } 

        // solve the linear system
        Lapack::Dgetrs('N', nTotModes, 1, M, nTotModes, pivot, f, nTotModes, info);
        if( info < 0 )
        {
            cerr << "ERROR: The " << -info << "th parameter had an illegal parameter for dgetrs" << endl;
            exit(1);
        }
        else if( info > 0 )
        {
            cerr << "ERROR: Element M[" << info << "] is 0 from dgetrs"<<endl;
            exit(1);
        }

        // Perform a backward transformation to obtain the solution at the quadrature points
        // Note that we are not using the sum-factorisation technique here 
        double** projectionResult = new double* [nQuadPointsDir1];
        for(i = 0; i < nQuadPointsDir1; i++)
        {
            projectionResult[i] = new double [nQuadPointsDir2];

            // Initialise with zeros
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                projectionResult[i][j] = 0.0;
            }
        }

        for(i = 0; i < nQuadPointsDir1; i++)
        {  
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                for(p = 0; p < nModesDir1; p++)
                { 
                    for(q = 0; q < nModesDir2; q++)
                    {               
                        n = indexmap[p][q];
                        
                        // Only if n is different from -1, we are dealing with a proper
                        // entrt of the mass matrix
                        if( (n!=-1) )
                        {  
                            projectionResult[i][j] += base1[p][i] * base2[p][q][j] * f[n];
                        }
                    }
                }
            }
        }

        // Display the output
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

                cout << "Quadrature point " << i*nQuadPointsDir2+j << ":       Result = ";
                cout << projectionResult[i][j] << endl;
                cout << "                     Exact Result = "<< pow(x1,6)*pow(x2,6) << endl; 
            }
            cout<<endl;
        }
        
        // Deallocate the dynamic memory
        delete[] quadZerosDir1;
        delete[] quadWeightsDir1;
        delete[] quadZerosDir2;
        delete[] quadWeightsDir2;
        delete[] JacobiPolDir1;
        delete[] JacobiPolDir2;
        
        for(i = 0; i < nModesDir1; i++)
        {
            delete[] base1[i];
        }        
        for(i = 0; i < nModesDir1; i++)
        {
            for(j = 0; j < nModesDir2; j++)
            {
                delete[] base2[i][j];
            }

            delete[] base2[i];
        }

        delete[] base1;
        delete[] base2;

        for(i = 0; i < nModesDir1; i++)
        {
            delete[] indexmap[i];
        }

        delete[] indexmap;
        delete[] M;
        delete[] f;

        for(i = 0; i < nQuadPointsDir1; i++)
        {
            delete[] projectionResult[i];
        } 
        delete[] projectionResult;
    }

    // METHOD 2: low-level Nektar++
    // METHOD 3: NEKTAR++
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
        SpatialDomains::TriGeomSharedPtr triGeom = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(zero,edges,eorient);
        triGeom->SetOwnData();

        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eModified_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eModified_B;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModesDir1,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModesDir2,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard triangular region
        StdRegions::StdExpansionSharedPtr triExpansion = 
            MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,triGeom);
        
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

        cout << "METHOD 2: low-level Nektar++" << endl;
        {              
            // Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes,nTotModes);
            
            // Calculate the mass matrix
            for(i = 0; i < nTotModes; i++)
            {
                Array<OneD,NekDouble> base_i(nTotQuadPoints);
                Array<OneD,NekDouble> column_i(nTotModes);
                
                // Fill the Array base_i with the values of the i-th (two-dimensional)
                // expansion mode evaluated at the quadrature points.
                triExpansion->FillMode(i,base_i);

                // Take the inner product of this expansion mode stored in base_i with all
                // the other expansion modes. This gives the i-th column (or row)
                // of the mass matrix
                triExpansion->IProductWRTBase(base_i,column_i);
                
                // Copy the result into the NekMatrix object M. (Note: the class NekMatrix holds
                // the data in a column-major format)
                Vmath::Vcopy(nTotModes, column_i.get(), 1, M.GetRawPtr() + i*nTotModes, 1);
            }

            // Declare and Allocate the right hand side vector 
            Array<OneD, NekDouble> f(nTotModes);

            // Calculate the right hand side vector f
            // Step 1: Multiply the forcing function with the Jacobian
            // The result is stored in the data member m_phys of the TriExp object triExpansion.
            NekDouble jacobian = ((triExpansion->GetMetricInfo())->GetJac())[0];
            Vmath::Smul(nTotQuadPoints, jacobian, forcingFunction, 1, triExpansion->UpdatePhys(), 1);

            // Step 2: Multiply the forcing function with quadrature weights
            // The result is stored in the data member m_phys of the TriExp object triExpansion.
            Array<OneD, const NekDouble> quadWeightsDir1 = (triExpansion->GetBasis(0))->GetW();
            Array<OneD, const NekDouble> quadWeightsDir2 = (triExpansion->GetBasis(1))->GetW();

            for(i = 0; i < nQuadPointsDir2; ++i)
            {
                Vmath::Vmul(nQuadPointsDir1, (triExpansion->GetPhys()).get()+i*nQuadPointsDir1, 1,
                            quadWeightsDir1.get(), 1, (triExpansion->UpdatePhys()).get()+i*nQuadPointsDir1, 1);
            }
            for(i = 0; i < nQuadPointsDir2; ++i)
            {
                Blas::Dscal(nQuadPointsDir1,0.5*quadWeightsDir2[i], (triExpansion->UpdatePhys()).get()+i*nQuadPointsDir1,1);      
            }

            // Step 3: Evaluate B^T*f using the sum-factorisation technique
            // The result is stored in the data member m_coeffs of the TriExp object triExpansion.
            const Array<OneD, const NekDouble> base0 = (triExpansion->GetBasis(0))->GetBdata();
            const Array<OneD, const NekDouble> base1 = (triExpansion->GetBasis(1))->GetBdata();

            Array<OneD, NekDouble>     tmpMat(nQuadPointsDir2*nModesDir1);
            NekMatrix<NekDouble>       in(nQuadPointsDir1,  nQuadPointsDir2, triExpansion->UpdatePhys(),eWrapper);
            NekMatrix<const NekDouble> B0(nQuadPointsDir1,  nModesDir1, base0, eWrapper);
            NekMatrix<NekDouble>       out(nQuadPointsDir2, nModesDir1,tmpMat,eWrapper);
                        
            out = Transpose(in)*B0;

            int mode;
            for(mode=i=0; i < nModesDir1; ++i)
            {
                NekVector<const NekDouble> in2(nQuadPointsDir2,tmpMat + i*nQuadPointsDir2,eWrapper);
                NekVector<NekDouble>       out2(nModesDir2-i,f + mode,eWrapper);
                NekMatrix<const NekDouble> B1(nQuadPointsDir2, nModesDir2-i, base1 + mode*nQuadPointsDir2, eWrapper);
                    
                out2 = Transpose(B1)*in2;
                mode += nModesDir2-i;
            }

            // Do a fix for the collapsed top vertex
            f[1] += Blas::Ddot(nQuadPointsDir2,base1+nQuadPointsDir2,1,
                               tmpMat+nQuadPointsDir2,1);

            // Solve the system
            // The result is stored in the data member m_coeffs of the TriExp object triExpansion.
            NekVector<NekDouble> in3(nTotModes, f, eWrapper);
            NekVector<NekDouble> out3(nTotModes, triExpansion->UpdateCoeffs(), eWrapper);

            M.Invert();
            out3 = M*in3;

            // Perform a backward transformation to obtain the solution at the quadrature points
            // using the sum-factorisation technique.
            // The result is stored in the data member m_phys of the TriExp object triExpansion.
            Array<OneD, NekDouble> tmpMat2(nQuadPointsDir2*nModesDir1);

            for(i = mode = 0; i < nModesDir1; ++i)
            {
                NekVector<const NekDouble> in4(nModesDir2-i,triExpansion->UpdateCoeffs()+mode,eWrapper);
                NekVector<NekDouble>       out4(nQuadPointsDir2,tmpMat2+i*nQuadPointsDir2,eWrapper);
                NekMatrix<const NekDouble> B1_2(nQuadPointsDir2,nModesDir2-i,base1 + mode*nQuadPointsDir2,eWrapper);

                out4 = B1_2*in4;
                mode += nModesDir2-i;
            }

            // Do a fix for the collapsed top vertex
            Blas::Daxpy(nQuadPointsDir2,(triExpansion->UpdateCoeffs())[1],base1.get()+nQuadPointsDir2,1,
                        tmpMat2.get()+nQuadPointsDir2,1);
      

            NekMatrix<const NekDouble> in5(nQuadPointsDir2,nModesDir1,tmpMat2,eWrapper);
            NekMatrix<NekDouble>       out5(nQuadPointsDir1,nQuadPointsDir2,triExpansion->UpdatePhys(),eWrapper);
            NekMatrix<const NekDouble> B0_2(nQuadPointsDir1,nModesDir1,base0,eWrapper);

            out5 = B0_2*Transpose(in5);

            // Display the output
            for(i = 0; i < nTotQuadPoints; i++)
            {               
                cout << "Quadrature point " << i << ":       Result = ";
                cout << (triExpansion->UpdatePhys())[i] << endl;
                cout << "                     Exact Result = "<< pow(x1[i],6)*pow(x2[i],6) << endl; 
            }
            cout << endl;
        }   
        
        cout << "METHOD 3: Nektar++" << endl;
        {  
            // Do the projection to obtain the coefficients of the expansion
            // The result is stored in the data member m_coeffs of the TriExp object triExpansion.
            triExpansion->FwdTrans(forcingFunction,triExpansion->UpdateCoeffs());

            // Perform a backward transformation to obtain the solution at the quadrature points
            // The result is stored in the data member m_phys of the TriExp object triExpansion.
            triExpansion->BwdTrans(triExpansion->GetCoeffs(),triExpansion->UpdatePhys());

            // Display the output
            for(i = 0; i < nTotQuadPoints; i++)
            {               
                cout << "Quadrature point " << i << ":       Result = ";
                cout << (triExpansion->GetPhys())[i] << endl;
                cout << "                     Exact Result = "<< pow(x1[i],6)*pow(x2[i],6) << endl; 
            }
            cout << endl;
        }    
    }
}
