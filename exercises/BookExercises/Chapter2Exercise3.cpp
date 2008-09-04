#include <cstdio>
#include <cstdlib>

#include <StdRegions/StdTriExp.h>
#include <LocalRegions/TriExp.h>
#include <StdRegions/StdSegExp.h>
#include <LocalRegions/SegExp.h>
#include <LibUtilities/Polylib/Polylib.h>


using namespace Nektar;

int main(int argc, char *argv[])
{
	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 3(a) --" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
		// Input: specify the details of the expansion
		int nQuadPoints =  10;
		int nModes = 9;

		// Declare variables to hold the quadrature zeros and weights
        double* quadZeros   = new double[nQuadPoints];
        double* quadWeights = new double[nQuadPoints];

		// Calculate the GLL-quadrature zeros and weights using the corresponding
        // Polylib routine.
        const double alpha = 0.0;
        const double beta  = 0.0;
        Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);
		Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);

		// Allocate memory for the matrix which will hold the values of the
        // (one-dimensional) basis at the quadrature points
        double** base = new double* [nModes];
       
        int i,j;
        for(i = 0; i < nModes; i++)
        {
            base[i] = new double[nQuadPoints];
        }     

		// Fill the matrix base with the Modal C0 basis
		double* JacobiPol = new double [nQuadPoints];
        for(i = 0; i < nModes; i++)
        {
            for(j = 0; j < nQuadPoints; j++)
            {
                if(i==0)
                {
                    base[i][j] = 0.5 * (1.0 - quadZeros[j]);
                }
                else if(i==nModes-1)
                {
                    base[i][j] = 0.5 * (1.0 + quadZeros[j]);
                }
                else
                {
					// Use Polylib to calculate the Jacobi Polynomial at the quadrature points
					Polylib::jacobfd(nQuadPoints, quadZeros, JacobiPol, NULL, i-1, 1.0, 1.0);
					base[i][j] = 0.5 * (1.0 - quadZeros[j]) * 0.5 * (1.0 + quadZeros[j]) * JacobiPol[j];
                }
            }
        }

		// Allocate memory for the mass matrix
        double** M = new double* [nModes];
        for(i = 0; i < nModes; i++)
        {
            M[i] = new double[nModes];
        } 

		// Calculate the mass matrix
        int p,q;
        double result;

		for(p = 0; p < nModes; p++)
        { 
			for(q = 0; q < nModes; q++)
			{
				result = 0.0;
				// Integrate 
                for(i = 0; i < nQuadPoints; i++)
                {
                    result += base[p][i] * base[q][i] * quadWeights[i];
                }
				M[p][q] = result;
			}
		}
		
		// Display the output
        cout << "Mass Matrix structure: " << endl;
		const double tol = 1e-12;
		for(i = 0; i < nModes; i++)
        {
            for(j = 0; j < nModes; j++)
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
        delete[] quadZeros;
        delete[] quadWeights;
        delete[] JacobiPol;
        
        for(i = 0; i < nModes; i++)
        {
            delete[] base[i];
        }        
        
		delete[] base;

        for(i = 0; i < nModes; i++)
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
        int nModes = 9;
        int nQuadPoints = 10;

		// Specify the type of the quadrature points
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
        
		// Specify the type of the basis
        LibUtilities::BasisType basisType = LibUtilities::eModified_A;

		// Declare the PointsKey which uniquely defines the quadrature points
        const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

		// Declare the BasisKey which uniquely defines (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKey(basisType,nModes,quadPointsKey);

		// Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // one dimensional standard region
        StdRegions::StdSegExpSharedPtr segExpansion = 
			MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);

		cout << "METHOD 2: low-level Nektar++" << endl;
        {
			int nTotModes      = segExpansion->GetNcoeffs();
            int nTotQuadPoints = segExpansion->GetTotPoints(); 

			// Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes,nTotModes);
		
			// Calculate the mass matrix
			for(i = 0; i < nTotModes; i++)
			{
				Array<OneD,NekDouble> base_i(nTotQuadPoints);
				Array<OneD,NekDouble> column_i(nTotModes);
	            
				// Fill the Array base_i with the values of the i-th (two-dimensional)
				// expansion mode evaluated at the quadrature points.
				segExpansion->FillMode(i,base_i);

				// Take the inner product of this expansion mode stored in base_i with all
				// the other expansion modes. This gives the i-th column (or row)
				// of the mass matrix
				segExpansion->IProductWRTBase(base_i,column_i);
	            
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
													   StdRegions::eSegment,
													   *segExpansion);
				// Get the mass matrix using the proper Nektar++ routine from the StdExpansion class.
				// This routine asks the StdMatrixManager for the mass matrix. If the matrix has been
				// created before, the manager simply returns the matrix. Otherwise, it calculates, stores
				// and returns the matrix  
				DNekMatSharedPtr M = segExpansion->GetStdMatrix(massMatrixKey);

				// Display the output
				cout << "Mass Matrix structure: " << endl;
	            
				const double tol = 1e-12;
				int nTotModes      = segExpansion->GetNcoeffs();            
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
    cout << "-- EXERCISE 3(b) --" <<endl;
    cout << "-------------------" <<endl;

	 cout << "METHOD 1: PolyLib" << endl;
    {
		// Input: specify the details of the expansion
		int nQuadPoints =  9; // changing the quadrature points to confirm the diagonality of the matrix
		int nModes = 9;

		// Declare variables to hold the quadrature zeros and weights
        double* quadZeros   = new double[nQuadPoints];
        double* quadWeights = new double[nQuadPoints];

		// Calculate the GLL-quadrature zeros and weights using the corresponding
        // Polylib routine.
        const double alpha = 0.0;
        const double beta  = 0.0;
        Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);
		Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);

		// Allocate memory for the matrix which will hold the values of the
        // (one-dimensional) basis at the quadrature points
        double** base = new double* [nModes];
       
        int i,j;
        for(i = 0; i < nModes; i++)
        {
            base[i] = new double[nQuadPoints];
        }     

		// Fill the matrix base0 with the orthogonal expansion
        double* nodalPoints  = new double[nModes];
        double* nodalWeights  = new double[nModes];
		Polylib::zwglj(nodalPoints,nodalWeights,nModes,alpha,beta);
		for(i = 0; i < nModes; i++)
        {
            for(j = 0; j < nQuadPoints; j++)
            {
                // Use Polylib to calculate the nodal basis functions at the quadrature points
                base[i][j] = Polylib::hglj(i,quadZeros[j],nodalPoints,nModes,0.0,0.0);
            }
        }

		// Allocate memory for the mass matrix
        double** M = new double* [nModes];
        for(i = 0; i < nModes; i++)
        {
            M[i] = new double[nModes];
        } 

		// Calculate the mass matrix
        int p,q;
        double result;

		for(p = 0; p < nModes; p++)
        { 
			for(q = 0; q < nModes; q++)
			{
				result = 0.0;
				// Integrate 
                for(i = 0; i < nQuadPoints; i++)
                {
                    result += base[p][i] * base[q][i] * quadWeights[i];
                }
				M[p][q] = result;
			}
		}
		
		// Display the output
        cout << "Mass Matrix structure: " << endl;
		const double tol = 1e-12;
		for(i = 0; i < nModes; i++)
        {
            for(j = 0; j < nModes; j++)
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
        delete[] quadZeros;
        delete[] quadWeights;
        
        for(i = 0; i < nModes; i++)
        {
            delete[] base[i];
        }        
        
		delete[] base;

        for(i = 0; i < nModes; i++)
        {
            delete[] M[i];
        }  

        delete[] M;

	}

	// METHOD 2: low-level Nektar++
    // METHOD 3: NEKTAR++
    {
		// Note: this implementation is completely analogous to exercise3(a). The
        // only difference is that the type of the basis is now eOrtho_A rather than
        // eModified_A.

		int i,j;
        // Input: specify the details of the expansion
        int nModes = 9;
        int nQuadPoints = 9;

		// Specify the type of the quadrature points
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
        
		// Specify the type of the basis
        LibUtilities::BasisType basisType = LibUtilities::eGLL_Lagrange;

		// Declare the PointsKey which uniquely defines the quadrature points
        const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

		// Declare the BasisKey which uniquely defines (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKey(basisType,nModes,quadPointsKey);

		// Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // one dimensional standard region
        StdRegions::StdSegExpSharedPtr segExpansion = 
			MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);

		cout << "METHOD 2: low-level Nektar++" << endl;
        {
			int nTotModes      = segExpansion->GetNcoeffs();
            int nTotQuadPoints = segExpansion->GetTotPoints(); 

			// Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes,nTotModes);
		
			// Calculate the mass matrix
			for(i = 0; i < nTotModes; i++)
			{
				Array<OneD,NekDouble> base_i(nTotQuadPoints);
				Array<OneD,NekDouble> column_i(nTotModes);
	            
				// Fill the Array base_i with the values of the i-th (two-dimensional)
				// expansion mode evaluated at the quadrature points.
				segExpansion->FillMode(i,base_i);

				// Take the inner product of this expansion mode stored in base_i with all
				// the other expansion modes. This gives the i-th column (or row)
				// of the mass matrix
				segExpansion->IProductWRTBase(base_i,column_i);
	            
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
													   StdRegions::eSegment,
													   *segExpansion);
				// Get the mass matrix using the proper Nektar++ routine from the StdExpansion class.
				// This routine asks the StdMatrixManager for the mass matrix. If the matrix has been
				// created before, the manager simply returns the matrix. Otherwise, it calculates, stores
				// and returns the matrix  
				DNekMatSharedPtr M = segExpansion->GetStdMatrix(massMatrixKey);

				// Display the output
				cout << "Mass Matrix structure: " << endl;
	            
				const double tol = 1e-12;
				int nTotModes      = segExpansion->GetNcoeffs();            
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
    cout << "-- EXERCISE 3(c) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: Lapack" << endl << endl;
    {
		// Input: specify the details of the expansion
        int nModes = 9;
        int nQuadPoints = 10;
		
		// Decalare variables to hold the quadrature zeros and weights
		double* quadZeros = new double[nQuadPoints];
		double* quadWeights = new double[nQuadPoints];

		// Calculate the GLL-quadrature zeros and weights using the corresponding
		// Polylib routine
		const double alpha = 0.0;
		const double beta = 0.0;
		Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);

		// Allocate memory for the matrix which will hold the values of the
        // (one-dimensional) basis at the quadrature points
        double** base = new double* [nModes];
       
        int i,j,p,q;
        for(i = 0; i < nModes; i++)
        {
            base[i] = new double[nQuadPoints];
        }     

		// Fill the matrix base0 with the orthogonal expansion
        double* nodalPoints  = new double[nModes];
        double* nodalWeights  = new double[nModes];
		Polylib::zwglj(nodalPoints,nodalWeights,nModes,alpha,beta);
		for(i = 0; i < nModes; i++)
        {
            for(j = 0; j < nQuadPoints; j++)
            {
                // Use Polylib to calculate the nodal basis functions at the quadrature points
                base[i][j] = Polylib::hglj(i,quadZeros[j],nodalPoints,nModes,0.0,0.0);
			}
		}

		// Calculate the right hand side vector fhat
		double* fhat = new double[nModes];
		for(p = 0; p < nModes; p++)
        {	
			fhat[p] = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				fhat[p] += base[p][i] * pow(quadZeros[i],7) * quadWeights[i];
			}
		}

		// Allocate memory for the mass matrix
        double* M = new double [nModes*nModes];
        for(i = 0; i < nModes*nModes; i++)
        {
            // Initialize the mass matrix with zeros
            M[i] = 0.0;
        }  
		
		
		// Calculate the mass matrix
        double result;

		for(p = 0; p < nModes; p++)
        { 
			for(q = 0; q < nModes; q++)
			{
				result = 0.0;
				// Integrate 
                for(i = 0; i < nQuadPoints; i++)
                {
                    result += base[p][i] * base[q][i] * quadWeights[i];
                }
				
				M[q+p*nModes] = result;
			}
		}

		// Calculate the LU factorization of the MASS matrix using the corresponding
		// Lapack routine
		int INFO = 0;
		int* ipiv = new int[nModes];
		Lapack::Dgetrf(nModes,nModes,M,nModes,ipiv,INFO);

		if( INFO < 0 )
        {
            cerr << "ERROR: The " << -INFO << "th parameter had an illegal parameter for dgetrf" << endl;
            exit(1);
        }
        else if( INFO > 0 )
        {
            cerr << "ERROR: Element M[" << INFO << "] is 0 from dgetrf"<<endl;
            exit(1);
        } 

		// Solve the linear system using the correponding Lapack routine
		// The resulting solution vector(uhat) will be stored in the right hand side vector fhat
		int NRHS = 1; // number of columns of the right hand side vector
		Lapack::Dgetrs('N',nModes,NRHS,M,nModes,ipiv,fhat,nModes,INFO);
		
		if( INFO < 0 )
        {
            cerr << "ERROR: The " << -INFO << "th parameter had an illegal parameter for dgetrs" << endl;
            exit(1);
        }
        else if( INFO > 0 )
        {
            cerr << "ERROR: Element M[" << INFO << "] is 0 from dgetrs"<<endl;
            exit(1);
        }

		// Calculate the solution
		double* uDelta = new double[nQuadPoints];
		for(i = 0; i < nQuadPoints; i++)
		{
			uDelta[i] = 0.0;
			// sum
			for(p = 0; p < nModes; p++)
			{
				uDelta[i] += fhat[p]*base[p][i];
			}
		}
		
		// Display the output
        for(i = 0; i < nQuadPoints; i++)
        {               
            cout << "Quadrature point " << i << ":       Result = ";
            cout << uDelta[i] << endl;
            cout << "                     Exact Result = "<< pow(quadZeros[i],7) << endl; 
			cout << endl;
        }
        cout << endl;
		
		// Deallocate the dynamic memory
        delete[] quadZeros;
        delete[] quadWeights;
        delete[] M;
		delete[] fhat;
		delete[] nodalPoints;
		delete[] nodalWeights;
		delete[] uDelta;

		for(i = 0; i < nModes; i++)
		{
			delete[] base[i];
		}
		delete[] base;

	}

	// METHOD 2: low-level Nektar++
    // METHOD 3: NEKTAR++
	{
		int i,j,p;
		// Input: specify the details of the expansion
        int nModes = 9;
        int nQuadPoints = 10;

		// Specify the type of the quadrature points
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		// Specify the type of the basis
        LibUtilities::BasisType basisType = LibUtilities::eGLL_Lagrange;

		// Declare the PointsKey which uniquely defines the quadrature points
        const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

		// Declare the BasisKey which uniquely defines (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKey(basisType,nModes,quadPointsKey);

		// Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // one dimensional standard region
        StdRegions::StdSegExpSharedPtr segExpansion = 
			MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);
		
		int nTotModes      = segExpansion->GetNcoeffs();
        int nTotQuadPoints = segExpansion->GetTotPoints(); 

		// Evaluate the forcing function at the quarature points
		Array<OneD, NekDouble> xi(nTotQuadPoints);
        Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
		segExpansion->GetCoords(xi);

		for(i = 0; i < nTotQuadPoints; i++)
        {
            forcingFunction[i] = pow(xi[i],7);
		}


		cout << "METHOD 2: low-level Nektar++" << endl << endl;
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
				segExpansion->FillMode(i,base_i);

				// Take the inner product of this expansion mode stored in base_i with all
				// the other expansion modes. This gives the i-th column (or row)
				// of the mass matrix
				segExpansion->IProductWRTBase(base_i,column_i);
	            
				// Copy the result into the NekMatrix object M. (Note: the class NekMatrix holds
				// the data in a column-major format)
				Vmath::Vcopy(nTotModes, column_i.get(), 1, M.GetRawPtr() + i*nTotModes, 1);
			}

			// Declare and Allocate the right hand side vector 
            Array<OneD, NekDouble> fhat(nTotModes);

			// Calculate the right hand side vector fhat
			Array<OneD, const NekDouble> quadWeights = (segExpansion->GetBasis(0))->GetW();
			const Array<OneD, const NekDouble> base = (segExpansion->GetBasis(0))->GetBdata();
			
			for(p = 0; p < nModes; p++)
			{	
				fhat[p] = 0.0;
				// Integrate
				for(i = 0; i < nQuadPoints; i++)
				{
					fhat[p] += base[i+p*nQuadPoints] * forcingFunction[i] * quadWeights[i];
				}
			}

			// Solve the system
			// The result will be stored in the array uhat
			Array<OneD, NekDouble> uhat(nTotModes);
			NekVector<NekDouble> in(nTotModes,fhat, eWrapper);
			NekVector<NekDouble> out(nTotModes,uhat, eWrapper);
			
			M.Invert();		
			out = M*in;
			
			// Calculate the solution
			Array<OneD, NekDouble> uDelta(nQuadPoints);
			for(i = 0; i < nQuadPoints; i++)
			{
				uDelta[i] = 0.0;
				// sum
				for(p = 0; p < nModes; p++)
				{
					uDelta[i] += uhat[p]*base[i+p*nQuadPoints];
				}
			}


			// Display the output
			for(i = 0; i < nTotQuadPoints; i++)
			{               
				cout << "Quadrature point " << i << ":       Result = ";
				cout << uDelta[i] << endl;
				cout << "                     Exact Result = "<< pow(xi[i],7) << endl; 
				cout << endl;
			}
			cout << endl;
		}
	
		
		cout << "METHOD 3: Nektar++" << endl << endl;
        { 
			// Do the projection to obtain the coefficients of the expansion
			// The result is stored in the data member m_coeffs of the StdSegExp object segExpansion.
			segExpansion->FwdTrans(forcingFunction,segExpansion->UpdateCoeffs());

			// Perform a backward transformation to obtain the solution at the quadrature points
			// The result is stored in the data member m_phys of the StdSegExp object segExpansion.
			segExpansion->BwdTrans(segExpansion->GetCoeffs(),segExpansion->UpdatePhys());
			
			// Display the output
            for(i = 0; i < nTotQuadPoints; i++)
            {               
                cout << "Quadrature point " << i << ":       Result = ";
                cout << (segExpansion->GetPhys())[i] << endl;
                cout << "                     Exact Result = "<< pow(xi[i],7) << endl; 
				cout << endl;
            }
            cout << endl;

		}
	}

	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 3(d) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: Lapack" << endl << endl;
    {
		// Input: specify the details of the expansion
        int nModes = 9;
        int nQuadPoints = 10;
		
		// Decalare variables to hold the quadrature zeros and weights
		double* quadZeros = new double[nQuadPoints];
		double* quadWeights = new double[nQuadPoints];

		// Calculate the GLL-quadrature zeros and weights using the corresponding
		// Polylib routine
		const double alpha = 0.0;
		const double beta = 0.0;
		Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);

		// Allocate memory for the matrix which will hold the values of the
        // (one-dimensional) basis at the quadrature points
        double** base = new double* [nModes];
       
        int i,j,p,q;
        for(i = 0; i < nModes; i++)
        {
            base[i] = new double[nQuadPoints];
        }     

		// Fill the matrix base0 with the orthogonal expansion
        double* nodalPoints  = new double[nModes];
        double* nodalWeights  = new double[nModes];
		Polylib::zwglj(nodalPoints,nodalWeights,nModes,alpha,beta);
		for(i = 0; i < nModes; i++)
        {
            for(j = 0; j < nQuadPoints; j++)
            {
                // Use Polylib to calculate the nodal basis functions at the quadrature points
                base[i][j] = Polylib::hglj(i,quadZeros[j],nodalPoints,nModes,0.0,0.0);
            }
        }

		// Calculate the right hand side vector fhat (not considering the boundary conditions)
		double* fhat = new double[nModes-1];
		for(p = 0 ; p < nModes-1; p++)
        {	
			fhat[p] = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				fhat[p] += base[p+1][i] * pow(quadZeros[i],7) * quadWeights[i];
			}
		}

		// Calculate the boundary conditions term
		double* bc = new double[nModes-1];
		double* u_d = new double[nQuadPoints];
		for(p = 0; p < nModes-1; p++)
        {	
			bc[p] = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				u_d[i] = (pow(-1.0,7)*base[0][i]+pow(1.0,7)*base[nModes-1][i]);
				bc[p] += quadWeights[i] * base[p+1][i] * u_d[i];
			}
		}

		//updating the right hand side vector fhat (considering the boundary conditions)
		for(p = 0; p < nModes-1; p++)
		{
			fhat[p] -= bc[p];
		}

		// Allocate memory for the mass matrix
        double* M = new double [(nModes-1)*(nModes-1)];
        for(i = 0; i < (nModes-1)*(nModes-1); i++)
        {
            // Initialize the mass matrix with zeros
            M[i] = 0.0;
        }  
		
		// Calculate the mass matrix
        double result;

		for(p = 0; p < nModes-1; p++)
        { 
			for(q = 0; q < nModes-1; q++)
			{
				result = 0.0;
				// Integrate 
                for(i = 0; i < nQuadPoints; i++)
                {
                    result += base[p+1][i] * base[q+1][i] * quadWeights[i];
                }
				
				M[q+p*(nModes-1)] = result;
			}
		}

		// Calculate the LU factorization of the MASS matrix using the corresponding
		// Lapack routine
		int INFO = 0;
		int* ipiv = new int[nModes-1];
		Lapack::Dgetrf(nModes-1,nModes-1,M,nModes-1,ipiv,INFO);

		if( INFO < 0 )
        {
            cerr << "ERROR: The " << -INFO << "th parameter had an illegal parameter for dgetrf" << endl;
            exit(1);
        }
        else if( INFO > 0 )
        {
            cerr << "ERROR: Element M[" << INFO << "] is 0 from dgetrf"<<endl;
            exit(1);
        } 

		// Solve the linear system using the correponding Lapack routine
		// The resulting solution vector(uhat) will be stored in the right hand side vector fhat
		int NRHS = 1; // number of columns of the right hand side vector
		Lapack::Dgetrs('N',nModes-1,NRHS,M,nModes-1,ipiv,fhat,nModes-1,INFO);
		
		if( INFO < 0 )
        {
            cerr << "ERROR: The " << -INFO << "th parameter had an illegal parameter for dgetrs" << endl;
            exit(1);
        }
        else if( INFO > 0 )
        {
            cerr << "ERROR: Element M[" << INFO << "] is 0 from dgetrs"<<endl;
            exit(1);
        }
		
		// Calculate the  solution
		double* uDelta = new double[nQuadPoints];
		for(i = 0; i < nQuadPoints; i++)
		{
			uDelta[i] = 0.0;
			// sum
			for(p = 0; p < nModes-1; p++)
			{
				uDelta[i] += fhat[p]*base[p+1][i];
			}
			uDelta[i] += u_d[i];
		}
		
		// Display the output
        for(i = 0 ; i < nQuadPoints; i++)
        {               
            cout << "Quadrature point " << i << ":       Result = ";
            cout << uDelta[i] << endl;
            cout << "                     Exact Result = "<< pow(quadZeros[i],7) << endl; 
			cout << endl;
        }
        cout << endl;
		
		// Deallocate the dynamic memory
        delete[] quadZeros;
        delete[] quadWeights;
        delete[] M;
		delete[] fhat;
		delete[] nodalPoints;
		delete[] nodalWeights;
		delete[] uDelta;

		for(i = 0; i < nModes; i++)
		{
			delete[] base[i];
		}
		delete[] base;

	}

	// METHOD 2: low-level Nektar++
    // METHOD 3: NEKTAR++
	{
		int i,j,p;
		// Input: specify the details of the expansion
        int nModes = 9;
        int nQuadPoints = 10;

		// Specify the type of the quadrature points
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		// Specify the type of the basis
        LibUtilities::BasisType basisType = LibUtilities::eGLL_Lagrange;

		// Declare the PointsKey which uniquely defines the quadrature points
        const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

		// Declare the BasisKey which uniquely defines (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKey(basisType,nModes,quadPointsKey);

		// Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // one dimensional standard region
        StdRegions::StdSegExpSharedPtr segExpansion = 
			MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);
		
		int nTotModes      = segExpansion->GetNcoeffs();
        int nTotQuadPoints = segExpansion->GetTotPoints(); 

		// Evaluate the forcing function at the quarature points
		Array<OneD, NekDouble> xi(nTotQuadPoints);
        Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
		segExpansion->GetCoords(xi);

		for(i = 0; i < nTotQuadPoints; i++)
        {
            forcingFunction[i] = pow(xi[i],7);
        }


		cout << "METHOD 2: low-level Nektar++" << endl;
        {
			// Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes-1,nTotModes-1);

			// Calculate the mass matrix
			for(i = 0; i < nTotModes-1; i++)
			{
				Array<OneD,NekDouble> base_i(nTotQuadPoints);
				Array<OneD,NekDouble> column_i(nTotModes);
	            
				// Fill the Array base_i with the values of the i+1-th (one-dimensional)
				// expansion mode evaluated at the quadrature points.
				segExpansion->FillMode(i+1,base_i);

				// Take the inner product of this expansion mode stored in base_i with all
				// the other expansion modes. This gives the i-th column (or row)
				// of the mass matrix
				segExpansion->IProductWRTBase(base_i,column_i);
	            
				// Copy the result into the NekMatrix object M. (Note: the class NekMatrix holds
				// the data in a column-major format)
				Vmath::Vcopy(nTotModes-1, column_i.get()+1, 1, M.GetRawPtr() + i*(nTotModes-1), 1);
			}
			
			// Declare and Allocate the right hand side vector 
            Array<OneD, NekDouble> fhat(nTotModes-1);

			// Calculate the right hand side vector fhat (not considering the bounary conditions)
            // Calculate the right hand side vector fhat
			Array<OneD, const NekDouble> quadWeights = (segExpansion->GetBasis(0))->GetW();
			const Array<OneD, const NekDouble> base = (segExpansion->GetBasis(0))->GetBdata();
			
			for(p = 0; p < nModes-1; p++)
			{	
				fhat[p] = 0.0;
				// Integrate
				for(i = 0; i < nQuadPoints; i++)
				{
					fhat[p] += base[i+(p+1)*nQuadPoints] * forcingFunction[i] * quadWeights[i];
				}
			}

			// Declare and Allocate the boundary condition terms
			Array<OneD, NekDouble> bc(nTotModes-1);

			//Calculate the boundary condition vector
			Array<OneD, NekDouble> u_d(nTotQuadPoints);
			for(p = 0; p < nModes-1; p++)
			{	
				bc[p] = 0.0;
				// Integrate
				for(i = 0; i < nQuadPoints; i++)
				{
					u_d[i] = (pow(-1.0,7)*base[i]+pow(1.0,7)*base[i+(nModes-1)*nQuadPoints]);
					bc[p] += quadWeights[i] * base[i+(p+1)*nQuadPoints] * u_d[i];
				}
			}

			// updating the right hand side vector fhat (considering the boundary condition terms)
			for(p = 0; p < nTotModes-1; p++)
			{
				fhat[p] -= bc[p];
			}

			// Solve the system
			// The result will be stored in the array uhat
			Array<OneD, NekDouble> uhat(nTotModes-1);
			NekVector<NekDouble> in(nTotModes-1,fhat,eWrapper);
			NekVector<NekDouble> out(nTotModes-1,uhat,eWrapper);
			
			M.Invert();		
			out = M*in;

			// Calculate the solution
			Array<OneD, NekDouble> uDelta(nQuadPoints);
			for(i = 0; i < nQuadPoints; i++)
			{
				uDelta[i] = 0.0;
				for(p = 0; p < nModes-1; p++)
				{
					uDelta[i] += base[i+(p+1)*nQuadPoints]*uhat[p];
				}
				uDelta[i] += u_d[i];
			}

			// Display the output
			for(i = 0; i < nTotQuadPoints; i++)
			{               
				cout << "Quadrature point " << i << ":       Result = ";
				cout << uDelta[i] << endl;
				cout << "                     Exact Result = "<< pow(xi[i],7) << endl;
				cout << endl;
			}
			cout << endl;
		}
			
		/*cout << "METHOD 3: Nektar++" << endl;
        { 
			// Do the projection to obtain the coefficients of the expansion
			// The result is stored in the data member m_coeffs of the StdSegExp object segExpansion.
			segExpansion->FwdTrans(forcingFunction,segExpansion->UpdateCoeffs());

			// Perform a backward transformation to obtain the solution at the quadrature points
			// The result is stored in the data member m_phys of the StdSegExp object segExpansion.
			segExpansion->BwdTrans(segExpansion->GetCoeffs(),segExpansion->UpdatePhys());
			
			// Display the output
            for(i = 0; i < nTotQuadPoints; i++)
            {               
                cout << "Quadrature point " << i << ":       Result = ";
                cout << (segExpansion->GetPhys())[i] << endl;
                cout << "                     Exact Result = "<< pow(xi[i],7) << endl; 
            }
            cout << endl;

		}*/
	}

	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 3(e) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: Lapack" << endl << endl;
    {
		// Specify the coordinate of the element
		double x1 = 2.0;
		double x2 = 5.0;

		// Input: specify the details of the expansion
        int nModes = 9;
        int nQuadPoints = 10;
		
		// Decalare variables to hold the quadrature zeros and weights
		double* quadZeros = new double[nQuadPoints];
		double* quadWeights = new double[nQuadPoints];

		// Calculate the GLL-quadrature zeros and weights using the corresponding
		// Polylib routine
		const double alpha = 0.0;
		const double beta = 0.0;
		Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);

		// Allocate memory for the matrix which will hold the values of the
        // (one-dimensional) basis at the quadrature points
        double** base = new double* [nModes];
       
        int i,j,p,q;
        for(i = 0; i < nModes; i++)
        {
            base[i] = new double[nQuadPoints];
        }     

		// Fill the matrix base0 with the orthogonal expansion
        double* nodalPoints  = new double[nModes];
        double* nodalWeights  = new double[nModes];
		Polylib::zwglj(nodalPoints,nodalWeights,nModes,alpha,beta);
		for(i = 0; i < nModes; i++)
        {
            for(j = 0; j < nQuadPoints; j++)
            {
                // Use Polylib to calculate the nodal basis functions at the quadrature points
                base[i][j] = Polylib::hglj(i,quadZeros[j],nodalPoints,nModes,0.0,0.0);
            }
        }

		// Analytically evaluate the Jacobian.
			double jacobian;
			jacobian = 0.5 * (-x1) + 0.5 * x2;

		// Calculate the right hand side vector fhat, linear mapping from standard to local element has been applied
		double* fhat = new double[nModes];
		for(p = 0; p < nModes; p++)
        {	
			fhat[p] = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = x1 * (1-quadZeros[i])/2 + x2 * (1+quadZeros[i])/2;

				fhat[p] += base[p][i] * pow(x,7) * quadWeights[i] * jacobian;
			}
		}

		// Allocate memory for the mass matrix
        double* M = new double [nModes*nModes];
        for(i = 0; i < nModes*nModes; i++)
        {
            // Initialize the mass matrix with zeros
            M[i] = 0.0;
        }  
		
		
		// Calculate the mass matrix
        double result;

		for(p = 0; p < nModes; p++)
        { 
			for(q = 0; q < nModes; q++)
			{
				result = 0.0;
				// Integrate 
                for(i = 0; i < nQuadPoints; i++)
                {
                    result += base[p][i] * base[q][i] * quadWeights[i];
                }
				
				M[q+p*nModes] = result * jacobian;
			}
		}
		
		// Calculate the LU factorization of the MASS matrix using the corresponding
		// Lapack routine
		int INFO = 0;
		int* ipiv = new int[nModes];
		Lapack::Dgetrf(nModes,nModes,M,nModes,ipiv,INFO);

		if( INFO < 0 )
        {
            cerr << "ERROR: The " << -INFO << "th parameter had an illegal parameter for dgetrf" << endl;
            exit(1);
        }
        else if( INFO > 0 )
        {
            cerr << "ERROR: Element M[" << INFO << "] is 0 from dgetrf"<<endl;
            exit(1);
        } 

		// Solve the linear system using the correponding Lapack routine
		// The resulting solution vector(uhat) will be stored in the right hand side vector fhat
		int NRHS = 1; // number of columns of the right hand side vector
		Lapack::Dgetrs('N',nModes,NRHS,M,nModes,ipiv,fhat,nModes,INFO);
		
		if( INFO < 0 )
        {
            cerr << "ERROR: The " << -INFO << "th parameter had an illegal parameter for dgetrs" << endl;
            exit(1);
        }
        else if( INFO > 0 )
        {
            cerr << "ERROR: Element M[" << INFO << "] is 0 from dgetrs"<<endl;
            exit(1);
        }

		// Calculate the solution
		double* uDelta = new double[nQuadPoints];
		for(i = 0; i < nQuadPoints; i++)
		{
			uDelta[i] = 0.0;
			// sum
			for(p = 0; p < nModes; p++)
			{
				uDelta[i] += fhat[p]*base[p][i];
			}
		}
		
		// Display the output
        for(i = 0; i < nQuadPoints; i++)
        {               
            cout << "Quadrature point " << i << ":       Result = ";
            cout << uDelta[i] << endl;
			double x = x1 * (1-quadZeros[i])/2 + x2 * (1+quadZeros[i])/2;
            cout << "                     Exact Result = "<< pow(x,7) << endl; 
			cout << endl;
        }
        cout << endl;
		
		// Deallocate the dynamic memory
        delete[] quadZeros;
        delete[] quadWeights;
        delete[] M;
		delete[] fhat;
		delete[] nodalPoints;
		delete[] nodalWeights;
		delete[] uDelta;

		for(i = 0; i < nModes; i++)
		{
			delete[] base[i];
		}
		delete[] base;

	}
	
	// METHOD 2: low-level Nektar++
    // METHOD 3: NEKTAR++
	{
		int i,j,p;
        // Input: specify the details of the expansion
        int nModes = 9;
        int nQuadPoints = 10;
       
        // Specify the geometry of the triangular element 
        // step 1: specify the coordinates of the vertices
        Array<OneD, NekDouble> coords(2);
        coords[0]    =   2.0;
        coords[1]    =   5.0;
            
        // step 2: set up the vertices
        SpatialDomains::VertexComponentSharedPtr verts[2];
        const int zero  = 0;
        const int one   = 1;
        const int two   = 2;
        const double yZero = 0.0;
		const double zZero = 0.0;
        verts[0] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(one,zero,coords[0],yZero,zZero);
        verts[1] = MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr(one,one,coords[1],yZero,zZero);
        
        // step 3: set up the geometry object
        SpatialDomains::Geometry1DSharedPtr segGeom =
			MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(zero,verts[0],verts[1]);
        segGeom->SetOwnData();

        LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
        
        LibUtilities::BasisType basisType = LibUtilities::eModified_A;
        
        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
        
        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKey(basisType,nModes,quadPointsKey);
        
        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard region
        StdRegions::StdExpansionSharedPtr segExpansion = 
            MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(basisKey,segGeom);

		int nTotModes      = segExpansion->GetNcoeffs();
        int nTotQuadPoints = segExpansion->GetTotPoints();  

		// Evaluate the forcing function at the quadrature points
        Array<OneD, NekDouble> x(nTotQuadPoints);
        Array<OneD, NekDouble> forcingFunction(nTotQuadPoints);
        segExpansion->GetCoords(x);

        for(i = 0; i < nTotQuadPoints; i++)
        {
            forcingFunction[i] = pow(x[i],7);
        }

		cout << "METHOD 2: low-level Nektar++" << endl << endl;
        { 
			// Declare and Allocate the mass matrix using an object of the class NekMatrix
            NekMatrix<NekDouble> M(nTotModes,nTotModes);

			// Calculate the mass matrix
            for(i = 0; i < nTotModes; i++)
            {
                Array<OneD,NekDouble> base_i(nTotQuadPoints);
                Array<OneD,NekDouble> column_i(nTotModes);
                
                // Fill the Array base_i with the values of the i-th
                // expansion mode evaluated at the quadrature points.
                segExpansion->FillMode(i,base_i);

                // Take the inner product of this expansion mode stored in base_i with all
                // the other expansion modes. This gives the i-th column (or row)
                // of the mass matrix
                segExpansion->IProductWRTBase(base_i,column_i);

				// Copy the result into the NekMatrix object M. (Note: the class NekMatrix holds
                // the data in a column-major format)
                Vmath::Vcopy(nTotModes, column_i.get(), 1, M.GetRawPtr() + i*nTotModes, 1);
            }

			// Declare and Allocate the right hand side vector 
            Array<OneD, NekDouble> fhat(nTotModes);

			// Calculate the right hand side vector fhat
            // Step 1: Multiply the forcing function with the Jacobian
			// The result is stored in the data member m_phys of the SegExp object segExpansion.
			NekDouble jacobian = ((segExpansion->GetMetricInfo())->GetJac())[0];
            Vmath::Smul(nTotQuadPoints, jacobian, forcingFunction, 1, segExpansion->UpdatePhys(), 1);

			// Step 2: Calculate the inner product of the result with the base expansions
			// Getting the quadrature weights and Base data
			Array<OneD, const NekDouble> quadWeights = (segExpansion->GetBasis(0))->GetW();
			const Array<OneD, const NekDouble> base = (segExpansion->GetBasis(0))->GetBdata();
			
			for(p = 0; p < nModes; p++)
			{	
				fhat[p] = 0.0;
				// Integrate
				for(i = 0; i < nQuadPoints; i++)
				{
					fhat[p] += base[i+(p)*nQuadPoints] * (segExpansion->GetPhys())[i] * quadWeights[i];
				}
			}

			// Solve the system
            // The result is stored in the data member m_coeffs of the SegExp object segExpansion.
            NekVector<NekDouble> in(nTotModes, fhat, eWrapper);
            NekVector<NekDouble> out(nTotModes, segExpansion->UpdateCoeffs(),eWrapper);

            M.Invert();
            out = M*in;
			
			// Calculate the solution
			Array<OneD, NekDouble> uDelta(nQuadPoints);
			for(i = 0; i < nQuadPoints; i++)
			{
				uDelta[i] = 0.0;
				for(p = 0; p < nModes; p++)
				{
					uDelta[i] += base[i+p*nQuadPoints]*(segExpansion->GetCoeffs())[p];
				}
			}

			// Display the output
            for(i = 0; i < nTotQuadPoints; i++)
            {               
                cout << "Quadrature point " << i << ":       Result = ";
                cout << uDelta[i] << endl;
                cout << "                     Exact Result = "<< pow(x[i],7) << endl; 
				cout << endl;
            }
            cout << endl;

		}
		
		cout << "METHOD 3: Nektar++" << endl << endl;
        {  
            // Do the projection to obtain the coefficients of the expansion
            // The result is stored in the data member m_coeffs of the SegExp object segExpansion.
            segExpansion->FwdTrans(forcingFunction,segExpansion->UpdateCoeffs());

            // Perform a backward transformation to obtain the solution at the quadrature points
            // The result is stored in the data member m_phys of the SegExp object segExpansion.
            segExpansion->BwdTrans(segExpansion->GetCoeffs(),segExpansion->UpdatePhys());

            // Display the output
            for(i = 0; i < nTotQuadPoints; i++)
            {               
                cout << "Quadrature point " << i << ":       Result = ";
                cout << (segExpansion->GetPhys())[i]<< endl;
                cout << "                     Exact Result = "<< pow(x[i],7) << endl; 
				cout << endl;
            }
            cout << endl;
        } 

	}

        return 0;
}

