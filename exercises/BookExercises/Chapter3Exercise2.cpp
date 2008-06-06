#include <cstdio>
#include <cstdlib>
#include <sys/time.h>

#include <StdRegions/StdTriExp.h>
#include <LibUtilities/Polylib/Polylib.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 3 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 2(a) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Input: specify the details of the expansion
        int nModesDir1      =  15;
        int nModesDir2      =  15;
        int nQuadPointsDir1 =  16;
        int nQuadPointsDir2 =  16;
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


        // Allocate memory for the mass matrix
        double** M = new double* [nTotModes];
        for(i = 0; i < nTotModes; i++)
        {
            M[i] = new double[nTotModes];

            // Initialise the mass matrix with zeros
            for(j = 0; j < nTotModes; j++)
            {
                M[i][j] = 0.0;
            }
        }  

        // Calculate the mass matrix
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
                            M[n][m] += resultDir1*resultDir2;
                        }
                    }
                }
            }
        }

        // Display the output
        cout << "Mass Matrix structure: " << endl;
        const double tol = 1e-12;
        for(i = 0; i < nTotModes; i++)
        {
            for(j = 0; j < nTotModes; j++)
            {
                if(fabs(M[i][j])<tol)
                {
                    cout<<"-";
                }
                else
                {
                    cout<<"o";
                }   
                
            }
            cout<<endl;
        }
        cout<<endl;
        
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

        for(i = 0; i < nTotModes; i++)
        {
            delete[] M[i];
        }  

        delete[] M;
    }

    // METHOD 2: low-level Nektar++
    // METHOD 3: NEKTAR++
    {
        int i,j;
        // Input: specify the details of the expansion
        int nModesDir1 = 15;
        int nModesDir2 = 15;
        int nQuadPointsDir1 = 16;
        int nQuadPointsDir2 = 16;

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
            MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);

        cout << "METHOD 2: low-level Nektar++" << endl;
        {
            int nTotModes      = triExpansion->GetNcoeffs();
            int nTotQuadPoints = triExpansion->GetTotPoints();    
            
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
            
            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "Mass Matrix structure: " << endl;

            const double tol = 1e-12;            
            for(i = 0; i < nTotModes; i++)
            {
                for(j = 0; j < nTotModes; j++)
                {
                    if(fabs(M(i,j))<tol)
                    {
                        cout<<"-";
                    }
                    else
                    {
                        cout<<"o";
                    }   
                    
                }
                cout<<endl;
            }  
            cout<<endl;
        }   
        cout << "METHOD 3: Nektar++" << endl;
        {
            // Define a MatrixKey which uniquely define the mass matrix   
            StdRegions::StdMatrixKey massMatrixKey(StdRegions::eMass,
                                                   StdRegions::eTriangle,
                                                   *triExpansion);
            
            // Get the mass matrix using the proper Nektar++ routine from the StdExpansion class.
            // This routine asks the StdMatrixManager for the mass matrix. If the matrix has been
            // created before, the manager simply returns the matrix. Otherwise, it calculates, stores
            // and returns the matrix  
            DNekMatSharedPtr M = triExpansion->GetStdMatrix(massMatrixKey);


            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "Mass Matrix structure: " << endl;
            
            const double tol = 1e-12;
            int nTotModes      = triExpansion->GetNcoeffs();            
            for(i = 0; i < nTotModes; i++)
            {
                for(j = 0; j < nTotModes; j++)
                {
                    if(fabs((*M)(i,j))<tol)
                    {
                        cout<<"-";
                    }
                    else
                    {
                        cout<<"o";
                    }   
                    
                }
                cout<<endl;
            } 
            cout<<endl; 
        }    
    }

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 2(b) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Note: this implementation is completely analogous to exercise2(a). The
        // only difference is that the matrices base1 and base2 are now filled with
        // orthogonal basis functions. 

        // Input: specify the details of the expansion
        int nModesDir1      =  15;
        int nModesDir2      =  15;
        int nQuadPointsDir1 =  16;
        int nQuadPointsDir2 =  16;
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

        // Fill the matrix base0 with the ortogonal basis functions \psi_i
        double* JacobiPolDir1 = new double [nQuadPointsDir1];
        double* JacobiPolDir2 = new double [nQuadPointsDir2];
        for(i = 0; i < nModesDir1; i++)
        {
            for(j = 0; j < nQuadPointsDir1; j++)
            {
                // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
                Polylib::jacobfd(nQuadPointsDir1, quadZerosDir1, JacobiPolDir1, NULL, i, 0.0, 0.0);
                base1[i][j] = JacobiPolDir1[j];
            }
        }

        // Fill the matrix base1 with the orthogonal basis functions \psi_{ij}
        // Note that were are storing the values for the entire index
        // range 0 <= i <= P_1 ; 0 <= j <= P_2.
        // Although not entries are required, this facilitates the implementation    
        for(i = 0; i < nModesDir1; i++)
        {
            for(j = 0; j < nModesDir2; j++)
            {                    
                // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
                Polylib::jacobfd(nQuadPointsDir2, quadZerosDir2, JacobiPolDir2, NULL, j, 2.0*i + 1.0, 0.0);
                
                for(k = 0; k < nQuadPointsDir2; k++)
                {
                    base2[i][j][k] = pow( 0.5 * (1.0 - quadZerosDir2[k]) , i) * JacobiPolDir2[k];
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


        // Allocate memory for the mass matrix
        double** M = new double* [nTotModes];
        for(i = 0; i < nTotModes; i++)
        {
            M[i] = new double[nTotModes];

            // Initialise the mass matrix with zeros
            for(j = 0; j < nTotModes; j++)
            {
                M[i][j] = 0.0;
            }
        }  

        // Calculate the mass matrix
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
                            M[n][m] += resultDir1*resultDir2;
                        }
                    }
                }
            }
        }

        // Display the output
        cout << "Mass Matrix structure: " << endl;
        cout << "(Note that the Mass matrix is not completely diagonal." <<endl;
        cout << "This is caused by numerical quadrature errors: As the " <<endl;
        cout << "basisfunction \phi_pq can be of order as high as 2P in" <<endl;
        cout << "the \eta_2 coordinate, a integration order of Q=2P+1 in" <<endl;
        cout << "direction 2 is required for an exact evaluation of all" <<endl;
        cout << "entries in the mass matrix. This can be checked by    " <<endl;
        cout << "setting QuadPointsDir2 = 29.)" <<endl;

        const double tol = 1e-12;
        for(i = 0; i < nTotModes; i++)
        {
            for(j = 0; j < nTotModes; j++)
            {
                if(fabs(M[i][j])<tol)
                {
                    cout<<"-";
                }
                else
                {
                    cout<<"o";
                }   
                
            }
            cout<<endl;
        }
        cout<<endl;
        
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

        for(i = 0; i < nTotModes; i++)
        {
            delete[] M[i];
        }  

        delete[] M;
    }

    // METHOD 2: low-level Nektar++
    // METHOD 3: NEKTAR++
    {        
        // Note: this implementation is completely analogous to exercise1(b). The
        // only difference is that the type of the basis is now eOrtho_A rather than
        // eModified_A.
        int i,j;
        // Input: specify the details of the expansion
        int nModesDir1 = 15;
        int nModesDir2 = 15;
        int nQuadPointsDir1 = 16;
        int nQuadPointsDir2 = 16;

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
        // standard triangular region
        StdRegions::StdExpansionSharedPtr triExpansion = 
            MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);

        cout << "METHOD 2: low-level Nektar++" << endl;
        {
            int nTotModes      = triExpansion->GetNcoeffs();
            int nTotQuadPoints = triExpansion->GetTotPoints();    
            
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
            
            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "(Note that the Mass matrix is not completely diagonal." <<endl;
            cout << "This is caused by numerical quadrature errors: As the " <<endl;
            cout << "basisfunction \phi_pq can be of order as high as 2P in" <<endl;
            cout << "the \eta_2 coordinate, a integration order of Q=2P+1 in" <<endl;
            cout << "direction 2 is required for an exact evaluation of all" <<endl;
            cout << "entries in the mass matrix. This can be checked by    " <<endl;
            cout << "setting QuadPointsDir2 = 29.)" <<endl;
            cout << "Mass Matrix structure: " << endl;

            const double tol = 1e-12;            
            for(i = 0; i < nTotModes; i++)
            {
                for(j = 0; j < nTotModes; j++)
                {
                    if(fabs(M(i,j))<tol)
                    {
                        cout<<"-";
                    }
                    else
                    {
                        cout<<"o";
                    }   
                    
                }
                cout<<endl;
            }  
            cout<<endl;
        }   
        cout << "METHOD 3: Nektar++" << endl;
        {
            // Define a MatrixKey which uniquely define the mass matrix   
            StdRegions::StdMatrixKey massMatrixKey(StdRegions::eMass,
                                                   StdRegions::eTriangle,
                                                   *triExpansion);
            
            // Get the mass matrix using the proper Nektar++ routine from the StdExpansion class.
            // This routine asks the StdMatrixManager for the mass matrix. If the matrix has been
            // created before, the manager simply returns the matrix. Otherwise, it calculates, stores
            // and returns the matrix  
            DNekMatSharedPtr M = triExpansion->GetStdMatrix(massMatrixKey);


            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "Mass Matrix structure: " << endl;
            
            const double tol = 1e-12;
            int nTotModes      = triExpansion->GetNcoeffs();            
            for(i = 0; i < nTotModes; i++)
            {
                for(j = 0; j < nTotModes; j++)
                {
                    if(fabs((*M)(i,j))<tol)
                    {
                        cout<<"-";
                    }
                    else
                    {
                        cout<<"o";
                    }   
                    
                }
                cout<<endl;
            } 
            cout<<endl; 
        }
    }
}
