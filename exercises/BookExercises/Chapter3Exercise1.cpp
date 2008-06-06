#include <cstdio>
#include <cstdlib>
#include <sys/time.h>

#include <StdRegions/StdQuadExp.h>
#include <LibUtilities/Polylib/Polylib.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 3 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(a) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Input: specify the number of quadrature points in both directions
        int nQuadPointsDir1 = 7;
        int nQuadPointsDir2 = 7;

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
        int i;
        double resultDir1 = 0.0;
        double resultDir2 = 0.0;
        double result;
        // Integrate in direction 1
        for(i = 0; i < nQuadPointsDir1; i++)
        {
            resultDir1 += pow(quadZerosDir1[i],4) * quadWeightsDir1[i];
        }
        // Integrate in direction 2
        for(i = 0; i < nQuadPointsDir2; i++)
        {
            resultDir2 += pow(quadZerosDir2[i],5) * quadWeightsDir2[i];
        }

        result = resultDir1*resultDir2;

        // Display the output
        cout << "Result = "<< result << endl << endl;

        // Deallocate the dynamic memory
        delete[] quadZerosDir1;
        delete[] quadWeightsDir1;
        delete[] quadZerosDir2;
        delete[] quadWeightsDir2;
    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    {
        // Input: specify the number and type of quadrature points in both directions
        int nQuadPointsDir1 = 7;
        int nQuadPointsDir2 = 7;
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
        int i;
        NekDouble resultDir1 = 0.0;
        NekDouble resultDir2 = 0.0;
        NekDouble result;
        // Integrate in direction 1
        for(i = 0; i < nQuadPointsDir1; i++)
        {
            resultDir1 += pow(quadZerosDir1[i],4) * quadWeightsDir1[i];
        }
        // Integrate in direction 2
        for(i = 0; i < nQuadPointsDir2; i++)
        {
            resultDir2 += pow(quadZerosDir2[i],5) * quadWeightsDir2[i];
        }

        result = resultDir1*resultDir2;

        // Display the output
        cout << "Result = "<< result << endl << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        // Input: specify the number and type of quadrature points in both directions
        int nQuadPointsDir1 = 7;
        int nQuadPointsDir2 = 7;
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
            integrand[i] = pow(xi1[i],4)*pow(xi2[i],5);
        }

        // Do the numerical integration by calling the corresponding Nektar++ routine of 
        // the StdExpansion class
        NekDouble result = quadExpansion->Integral(integrand);

        // Display the output
        cout << "Result = "<< result << endl << endl;
    }

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(b) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Input: specify the details of the expansion
        int nModesDir1      =  8;
        int nModesDir2      =  8;
        int nQuadPointsDir1 =  9;
        int nQuadPointsDir2 =  9;
        int nTotModes       =  nModesDir1*nModesDir2;
        int nTotQuadPoints  =  nQuadPointsDir1*nQuadPointsDir1;     

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

        // Allocate memory for the matrices which will hold the values of the
        // (one-dimensional) bases at the quadrature points
        double** base1 = new double* [nModesDir1];
        double** base2 = new double* [nModesDir2];

        int i,j;
        for(i = 0; i < nModesDir1; i++)
        {
            base1[i] = new double[nQuadPointsDir1];
        }        
        for(i = 0; i < nModesDir2; i++)
        {
            base2[i] = new double[nQuadPointsDir2];
        }

        // Fill the matrix base0 with the Modal C0 basis
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

        // Fill the matrix base1 with the Modal C0 basis
        for(i = 0; i < nModesDir2; i++)
        {
            if(i==0)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {
                    base2[i][j] = 0.5 * (1.0 - quadZerosDir2[j]);
                }
            }
            
            else if(i==nModesDir2-1)
            {
                for(j = 0; j < nQuadPointsDir2; j++)
                {
                    base2[i][j] = 0.5 * (1.0 + quadZerosDir2[j]);
                }
            }
            else
            {
                // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
                Polylib::jacobfd(nQuadPointsDir2, quadZerosDir2, JacobiPolDir2, NULL, i-1, 1.0, 1.0);
                
                for(j = 0; j < nQuadPointsDir2; j++)
                {
                    base2[i][j] = 0.5 * (1.0 - quadZerosDir2[j]) * 0.5 * (1.0 + quadZerosDir2[j]) * JacobiPolDir2[j];
                }
            }
        }

        // Allocate memory for the mass matrix
        double** M = new double* [nTotModes];
        for(i = 0; i < nTotModes; i++)
        {
            M[i] = new double[nTotModes];
        }  

        // Calculate the mass matrix
        int p,q,r,s,m,n;
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
                            resultDir2 += base2[q][i] * base2[s][i] * quadWeightsDir2[i];
                        }

                        n = p*nModesDir1 + q;
                        m = r*nModesDir2 + s;

                        M[n][m] = resultDir1*resultDir2;
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
        for(i = 0; i < nModesDir2; i++)
        {
            delete[] base2[i];
        }

        delete[] base1;
        delete[] base2;

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
        int nModesDir1 = 8;
        int nModesDir2 = 8;
        int nQuadPointsDir1 = 9;
        int nQuadPointsDir2 = 9;

        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eModified_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eModified_A;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModesDir1,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModesDir2,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard quadrilateral region
        StdRegions::StdExpansionSharedPtr quadExpansion = 
            MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);

        cout << "METHOD 2: low-level Nektar++" << endl;
        {
            int nTotModes      = quadExpansion->GetNcoeffs();
            int nTotQuadPoints = quadExpansion->GetTotPoints();    
            
            // Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes,nTotModes);
            
            // Calculate the mass matrix
            for(i = 0; i < nTotModes; i++)
            {
                Array<OneD,NekDouble> base_i(nTotQuadPoints);
                Array<OneD,NekDouble> column_i(nTotModes);
                
                // Fill the Array base_i with the values of the i-th (two-dimensional)
                // expansion mode evaluated at the quadrature points.
                quadExpansion->FillMode(i,base_i);

                // Take the inner product of this expansion mode stored in base_i with all
                // the other expansion modes. This gives the i-th column (or row)
                // of the mass matrix
                quadExpansion->IProductWRTBase(base_i,column_i);
                
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
                                                   StdRegions::eQuadrilateral,
                                                   *quadExpansion);
            
            // Get the mass matrix using the proper Nektar++ routine from the StdExpansion class.
            // This routine asks the StdMatrixManager for the mass matrix. If the matrix has been
            // created before, the manager simply returns the matrix. Otherwise, it calculates, stores
            // and returns the matrix  
            DNekMatSharedPtr M = quadExpansion->GetStdMatrix(massMatrixKey);


            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "Mass Matrix structure: " << endl;
            
            const double tol = 1e-12;
            int nTotModes      = quadExpansion->GetNcoeffs();            
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
    cout << "-- EXERCISE 1(c) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Note: this implementation is completely analogous to exercise1(b). The
        // only difference is that the matrices base1 and base2 are now filled with
        // orthogonal basis functions.

        // Input: specify the details of the expansion
        int nModesDir1      =  8;
        int nModesDir2      =  8;
        int nQuadPointsDir1 =  9;
        int nQuadPointsDir2 =  9;
        int nTotModes       =  nModesDir1*nModesDir2;
        int nTotQuadPoints  =  nQuadPointsDir1*nQuadPointsDir1;     

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

        // Allocate memory for the matrices which will hold the values of the
        // (one-dimensional) bases at the quadrature points
        double** base1 = new double* [nModesDir1];
        double** base2 = new double* [nModesDir2];

        int i,j;
        for(i = 0; i < nModesDir1; i++)
        {
            base1[i] = new double[nQuadPointsDir1];
        }        
        for(i = 0; i < nModesDir2; i++)
        {
            base2[i] = new double[nQuadPointsDir2];
        }

        // Fill the matrix base0 with the orthogonal expansion
        double* JacobiPolDir1 = new double [nQuadPointsDir1];
        double* JacobiPolDir2 = new double [nQuadPointsDir2];
        for(i = 0; i < nModesDir1; i++)
        {
            for(j = 0; j < nQuadPointsDir1; j++)
            {
                // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
                Polylib::jacobfd(nQuadPointsDir1, quadZerosDir1, JacobiPolDir1, NULL, i, 0.0, 0.0);
                base1[i][j] = sqrt(0.5*(2.0*i+1.0)) * JacobiPolDir1[j];
            }
        }

        // Fill the matrix base1 with the Modal C0 basis
        for(i = 0; i < nModesDir2; i++)
        {
            // Use Polylib to calculate the Jacobi Polynomial at the quadrature points
            Polylib::jacobfd(nQuadPointsDir2, quadZerosDir2, JacobiPolDir2, NULL, i, 0.0, 0.0);
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                base2[i][j] = sqrt(0.5*(2.0*i+1.0)) * JacobiPolDir2[j];
            }
        }

        // Allocate memory for the mass matrix
        double** M = new double* [nTotModes];
        for(i = 0; i < nTotModes; i++)
        {
            M[i] = new double[nTotModes];
        }  

        // Calculate the mass matrix
        int p,q,r,s,m,n;
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
                            resultDir2 += base2[q][i] * base2[s][i] * quadWeightsDir2[i];
                        }

                        n = p*nModesDir1 + q;
                        m = r*nModesDir2 + s;

                        M[n][m] = resultDir1*resultDir2;
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
        for(i = 0; i < nModesDir2; i++)
        {
            delete[] base2[i];
        }

        delete[] base1;
        delete[] base2;

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
        int nModesDir1 = 8;
        int nModesDir2 = 8;
        int nQuadPointsDir1 = 9;
        int nQuadPointsDir2 = 9;

        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eOrtho_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eOrtho_A;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModesDir1,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModesDir2,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard quadrilateral region
        StdRegions::StdExpansionSharedPtr quadExpansion = 
            MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);

        cout << "METHOD 2: low-level Nektar++" << endl;
        {
            int nTotModes      = quadExpansion->GetNcoeffs();
            int nTotQuadPoints = quadExpansion->GetTotPoints();    
            
            // Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes,nTotModes);
            
            // Calculate the mass matrix
            for(i = 0; i < nTotModes; i++)
            {
                Array<OneD,NekDouble> base_i(nTotQuadPoints);
                Array<OneD,NekDouble> column_i(nTotModes);
                
                // Fill the Array base_i with the values of the i-th (two-dimensional)
                // expansion mode evaluated at the quadrature points.
                quadExpansion->FillMode(i,base_i);

                // Take the inner product of this expansion mode stored in base_i with all
                // the other expansion modes. This gives the i-th column (or row)
                // of the mass matrix
                quadExpansion->IProductWRTBase(base_i,column_i);
                
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
                                                   StdRegions::eQuadrilateral,
                                                   *quadExpansion);
            
            // Get the mass matrix using the proper Nektar++ routine from the StdExpansion class.
            // This routine asks the StdMatrixManager for the mass matrix. If the matrix has been
            // created before, the manager simply returns the matrix. Otherwise, it calculates, stores
            // and returns the matrix  
            DNekMatSharedPtr M = quadExpansion->GetStdMatrix(massMatrixKey);


            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "Mass Matrix structure: " << endl;
            
            const double tol = 1e-12;
            int nTotModes      = quadExpansion->GetNcoeffs();            
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
    cout << "-- EXERCISE 1(d) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Note: this implementation is completely analogous to exercise1(b). The
        // only difference is that the matrices base1 and base2 are now filled with
        // nodal basis functions.

        // Input: specify the details of the expansion
        int nModesDir1      =  8;
        int nModesDir2      =  8;
        int nQuadPointsDir1 =  9; //change the number of quadrature points here to toggle the discrete orthogonality
        int nQuadPointsDir2 =  9; //change the number of quadrature points here to toggle the discrete orthogonality
        int nTotModes       =  nModesDir1*nModesDir2;
        int nTotQuadPoints  =  nQuadPointsDir1*nQuadPointsDir1;     

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

        // Allocate memory for the matrices which will hold the values of the
        // (one-dimensional) bases at the quadrature points
        double** base1 = new double* [nModesDir1];
        double** base2 = new double* [nModesDir2];

        int i,j;
        for(i = 0; i < nModesDir1; i++)
        {
            base1[i] = new double[nQuadPointsDir1];
        }        
        for(i = 0; i < nModesDir2; i++)
        {
            base2[i] = new double[nQuadPointsDir2];
        }

        // Fill the matrix base0 with the orthogonal expansion
        double* nodalPointsDir1   = new double[nModesDir1];
        double* nodalWeightsDir1  = new double[nModesDir1];
        double* nodalPointsDir2   = new double[nModesDir2];
        double* nodalWeightsDir2  = new double[nModesDir2];
        Polylib::zwglj(nodalPointsDir1,nodalWeightsDir1,nModesDir1,alpha,beta);
	Polylib::zwglj(nodalPointsDir2,nodalWeightsDir2,nModesDir2,alpha,beta);
        for(i = 0; i < nModesDir1; i++)
        {
            for(j = 0; j < nQuadPointsDir1; j++)
            {
                // Use Polylib to calculate the nodal basis functions at the quadrature points
                base1[i][j] = Polylib::hglj(i,quadZerosDir1[j],nodalPointsDir1,nModesDir1,0.0,0.0);
            }
        }

        // Fill the matrix base1 with the Modal C0 basis
        for(i = 0; i < nModesDir2; i++)
        {
            for(j = 0; j < nQuadPointsDir2; j++)
            {
                // Use Polylib to calculate the nodal basis functions at the quadrature points
                base2[i][j] = Polylib::hglj(i,quadZerosDir2[j],nodalPointsDir2,nModesDir2,0.0,0.0);
            }
        }

        // Allocate memory for the mass matrix
        double** M = new double* [nTotModes];
        for(i = 0; i < nTotModes; i++)
        {
            M[i] = new double[nTotModes];
        }  

        // Calculate the mass matrix
        int p,q,r,s,m,n;
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
                            resultDir2 += base2[q][i] * base2[s][i] * quadWeightsDir2[i];
                        }

                        n = p*nModesDir1 + q;
                        m = r*nModesDir2 + s;

                        M[n][m] = resultDir1*resultDir2;
                    }
                }
            }
        }

        // Display the output
        cout << "Mass Matrix structure: (change the number of quadrature points to toggle the discrete orthogonality)" << endl;
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
        delete[] nodalPointsDir1;
        delete[] nodalWeightsDir1;
        delete[] nodalPointsDir2;
        delete[] nodalWeightsDir2;
        
        for(i = 0; i < nModesDir1; i++)
        {
            delete[] base1[i];
        }        
        for(i = 0; i < nModesDir2; i++)
        {
            delete[] base2[i];
        }

        delete[] base1;
        delete[] base2;

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
        // only difference is that the type of the basis is now eGLL_Lagrange rather than
        // eModified_A.

        int i,j;
        // Input: specify the details of the expansion
        int nModesDir1 = 8;
        int nModesDir2 = 8;
        int nQuadPointsDir1 = 9; //change the number of quadrature points here to toggle the discrete orthogonality
        int nQuadPointsDir2 = 9; //change the number of quadrature points here to toggle the discrete orthogonality

        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussLobattoLegendre;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eGLL_Lagrange;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eGLL_Lagrange;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModesDir1,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModesDir2,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard quadrilateral region
        StdRegions::StdExpansionSharedPtr quadExpansion = 
            MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);

        cout << "METHOD 2: low-level Nektar++" << endl;
        {
            int nTotModes      = quadExpansion->GetNcoeffs();
            int nTotQuadPoints = quadExpansion->GetTotPoints();    
            
            // Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes,nTotModes);
            
            // Calculate the mass matrix
            for(i = 0; i < nTotModes; i++)
            {
                Array<OneD,NekDouble> base_i(nTotQuadPoints);
                Array<OneD,NekDouble> column_i(nTotModes);
                
                // Fill the Array base_i with the values of the i-th (two-dimensional)
                // expansion mode evaluated at the quadrature points.
                quadExpansion->FillMode(i,base_i);

                // Take the inner product of this expansion mode stored in base_i with all
                // the other expansion modes. This gives the i-th column (or row)
                // of the mass matrix
                quadExpansion->IProductWRTBase(base_i,column_i);
                
                // Copy the result into the NekMatrix object M. (Note: the class NekMatrix holds
                // the data in a column-major format)
                Vmath::Vcopy(nTotModes, column_i.get(), 1, M.GetRawPtr() + i*nTotModes, 1);
            }
            
            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "Mass Matrix structure: (change the number of quadrature points to toggle the discrete orthogonality)" << endl;

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
                                                   StdRegions::eQuadrilateral,
                                                   *quadExpansion);
            
            // Get the mass matrix using the proper Nektar++ routine from the StdExpansion class.
            // This routine asks the StdMatrixManager for the mass matrix. If the matrix has been
            // created before, the manager simply returns the matrix. Otherwise, it calculates, stores
            // and returns the matrix  
            DNekMatSharedPtr M = quadExpansion->GetStdMatrix(massMatrixKey);


            // Display the output
            // The matrix structure might be different as in the Polylib version. This is
            // because of a difference in numbering of the expansion modes.
            cout << "Mass Matrix structure: (change the number of quadrature points to toggle the discrete orthogonality)" << endl;
            
            const double tol = 1e-12;
            int nTotModes      = quadExpansion->GetNcoeffs();            
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

