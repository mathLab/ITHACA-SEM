#include <cstdio>
#include <cstdlib>
#include <sys/time.h>

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/Polylib/Polylib.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 3 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 3(a+c)-" <<endl;
    cout << "-------------------" <<endl;

    cout << "METHOD 1: PolyLib" << endl;
    {
        // Input: specify the details of the expansion
        int nModes1D        =  5;
        int nTotModes       =  nModes1D*(nModes1D+1)/2;

        // Declare variables to hold nodal points
        double* nodalPointsXi1  = new double[nTotModes];
        double* nodalPointsXi2  = new double[nTotModes];
        double* nodalPointsEta1 = new double[nTotModes];
        double* nodalPointsEta2 = new double[nTotModes];

        // Calculate the equidistant nodal points
        // Also calcualte the nodal points in eta1-eta2 coordinates
        int i,j;
        const double tol = 1e-12;
        int cnt = 0;
        for(i = 0; i < nModes1D; i++)
        {
            for(j = 0 ; j <= nModes1D-1-i; j++)
            {
                nodalPointsXi1[cnt]   = 2.0*i/(nModes1D-1.0)-1.0;
                nodalPointsXi2[cnt] = 2.0*j/(nModes1D-1.0)-1.0;


                if( fabs(nodalPointsXi2[cnt]-1.0)<tol )
                {
                    nodalPointsEta1[cnt]=0.0;
                }
                else
                {
                    nodalPointsEta1[cnt] = 2*(1+nodalPointsXi1[cnt])/(1-nodalPointsXi2[cnt]);
                }
                nodalPointsEta2[cnt] = nodalPointsXi2[cnt];

                cnt++;
            }
        }

        // Allocate memory for the generalised Vandermonde matrix
        double* V = new double [nTotModes*nTotModes];

        // Calculate the Vandermonde matrix
        int p,q,m;
        double psi_p[1];
        double psi_pq[1];


        cnt = 0;
        for(i = 0; i < nTotModes; i++)
        {            
            for(p = 0; p < nModes1D; p++)
            {
                for(q = 0 ; q <= nModes1D-1-p; q++)
                {                    
                    Polylib::jacobfd(1, &(nodalPointsEta1[i]), psi_p,  NULL, p, 0.0, 0.0);
                    Polylib::jacobfd(1, &(nodalPointsEta2[i]), psi_pq, NULL, q, 2.0*p + 1.0, 0.0);
                    
                    V[cnt] = psi_p[0] * pow((1.0-nodalPointsEta2[i])/2.0 , p) * psi_pq[0];
                    cnt++;
                }
            }
        }

        // Invert the Vandermonde matrix
        int*   pivot = new int[nTotModes];
        double* work = new double[nTotModes];
        int     info = 0;
                    
        Lapack::Dgetrf(nTotModes, nTotModes, V, nTotModes, pivot, info);
        
        if( info < 0 )
        {
            cerr << "ERROR: The " << -info << "th parameter had an illegal parameter for dgetrf" << endl;
            exit(1);
        }
        else if( info > 0 )
        {
            cerr << "ERROR: Element V[" << info << "] is 0 from dgetrf"<<endl;
            exit(1);
        }   
        
        Lapack::Dgetri(nTotModes, V, nTotModes, pivot, work, nTotModes, info);
                    
        if( info < 0 )
        {
            cerr << "ERROR: The " << -info << "th parameter had an illegal parameter for dgetri" << endl;
            exit(1);
        }
        else if( info > 0 )
        {
            cerr << "ERROR: Element V[" << info << "] is 0 from dgetri"<<endl;
            exit(1);
        }   


        // specify the number of points where the Lagrange polynomials
        // should be evaluated
        int nEvalPoints1D    =  9;
        int nTotEvalPoints   =  nEvalPoints1D*(nEvalPoints1D+1)/2;

        // Declare variables to hold evaluation points
        double* evalPointsXi1  = new double[nTotEvalPoints];
        double* evalPointsXi2  = new double[nTotEvalPoints];
        double* evalPointsEta1 = new double[nTotEvalPoints];
        double* evalPointsEta2 = new double[nTotEvalPoints];

        // Calculate the equidistant evaluation points
        // Also calculate the evaluation points in eta1-eta2 coordinates
        cnt = 0;
        for(i = 0; i < nEvalPoints1D; i++)
        {
            for(j = 0 ; j <= nEvalPoints1D-1-i; j++)
            {
                evalPointsXi1[cnt]   = 2.0*i/(nEvalPoints1D-1.0)-1.0;
                evalPointsXi2[cnt] = 2.0*j/(nEvalPoints1D-1.0)-1.0;


                if( fabs(evalPointsXi2[cnt]-1.0)<tol )
                {
                    evalPointsEta1[cnt]=0.0;
                }
                else
                {
                    evalPointsEta1[cnt] = 2*(1+evalPointsXi1[cnt])/(1-evalPointsXi2[cnt]);
                }
                evalPointsEta2[cnt] = evalPointsXi2[cnt];

                cnt++;
            }
        }


        // Declare variables to hold the values of the basis functions
        // at the evaluation points
        double* orthoBasisEvaluated  = new double[nTotModes];
        double* lagBasisEvaluated    = new double[nTotModes*nTotEvalPoints];

        // Evaluate the Lagrange polynomials at the evaluation points
        for(i = 0; i < nTotEvalPoints; i++)
        {         
            // Step 1: Evaluate the Orthogonal Basis at the evaluation points   
            cnt=0;
            for(p = 0; p < nModes1D; p++)
            {
                for(q = 0 ; q <= nModes1D-1-p; q++)
                {                    
                    Polylib::jacobfd(1, &(evalPointsEta1[i]), psi_p,  NULL, p, 0.0, 0.0);
                    Polylib::jacobfd(1, &(evalPointsEta2[i]), psi_pq, NULL, q, 2.0*p + 1.0, 0.0);
                    
                    orthoBasisEvaluated[cnt++] = psi_p[0] * pow((1.0-evalPointsEta2[i])/2.0 , p) * psi_pq[0];
                }
            }

            // Step 2: multiply by the inverse of the Vandermonde matrix to 
            // get the values of the Lagrange polynomial.
            // This matrix-vector multiplication is done by the corresponding
            // BLAS call. 
            Blas::Dgemv('N',nTotModes,nTotModes,1.0,V,nTotModes,orthoBasisEvaluated,1,0.0,lagBasisEvaluated+i,nTotEvalPoints);
        }

        // Exercise 3(c): Calculate the Lebesgue function for these equidistant nodal points
        double* lebesgueFunction  = new double[nTotEvalPoints];
        for(i = 0; i < nTotEvalPoints; i++)
        {
            lebesgueFunction[i] = 0.0;

            for(j = 0; j < nTotModes; j++)
            {
                lebesgueFunction[i] += fabs(lagBasisEvaluated[i+nTotEvalPoints*j]);
            }
        }

        // Display the output
        cout << "The values of the different Lagrange polynomials evaluated at the " << nTotEvalPoints << endl; 
        cout << "equidistant points are now contained in the variable:" << endl;
        cout << "  - lagBasisEvaluated" << endl;
        cout << "The values of the Lebesgue function evaluated at the same points" << endl; 
        cout << "are contained in the variable:" << endl;
        cout << "  - lebesgueFunction" << endl;
        cout << "(If desired, write your own routine to cast these data in such a format" << endl;
        cout<< "that it can be read by your post-processing program of choice.)" << endl << endl;

        // Deallocate the dynamic memory
        delete[] nodalPointsXi1;
        delete[] nodalPointsXi2;
        delete[] nodalPointsEta1;
        delete[] nodalPointsEta2;
        delete[] V;
        delete[] pivot;
        delete[] work;
        delete[] evalPointsXi1;
        delete[] evalPointsXi2;
        delete[] evalPointsEta1;
        delete[] evalPointsEta2;
        delete[] orthoBasisEvaluated;
        delete[] lagBasisEvaluated;
        delete[] lebesgueFunction;

    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    {
        // Input: specify the details of the expansion
        int nModes1D        =  5;
        int nQuadPointsDir1 =  6;
        int nQuadPointsDir2 =  6;
        int nTotModes       =  nModes1D*(nModes1D+1)/2;

        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eOrtho_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eOrtho_B;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModes1D,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModes1D,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard triangular region
        StdRegions::StdTriExpSharedPtr triExpansion = 
            MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);

        // Declare variables to hold nodal points
        Array<OneD, NekDouble> nodalPointsXi1(nTotModes);
        Array<OneD, NekDouble> nodalPointsXi2(nTotModes);

        // Calculate the equidistant nodal points
        // Also calcualte the nodal points in eta1-eta2 coordinates
        int i,j;
        int cnt = 0;
        for(i = 0; i < nModes1D; i++)
        {
            for(j = 0 ; j <= nModes1D-1-i; j++)
            {
                nodalPointsXi1[cnt]   = 2.0*i/(nModes1D-1.0)-1.0;
                nodalPointsXi2[cnt++] = 2.0*j/(nModes1D-1.0)-1.0;
            }
        }

        // Allocate memory for the generalised Vandermonde matrix
        NekMatrix<NekDouble> V(nTotModes,nTotModes);

        // Calculate the Vandermonde matrix
        Array<OneD, NekDouble> coord(2);
        for(i = 0; i < nTotModes; i++)
        {
            // Fill the physical space of triExpansion with orthogonal mode i
            triExpansion->FillMode(i,triExpansion->UpdatePhys());

            // Interpolate mode i to the nodal points j
            for(j = 0; j < nTotModes; j++)
            {
                    coord[0] = nodalPointsXi1[j];
                    coord[1] = nodalPointsXi2[j];
                    V(i,j) = triExpansion->PhysEvaluate(coord);                
            }
        }

        // Invert the Vandermonde matrix
        V.Invert();

        // specify the number of points where the Lagrange polynomials
        // should be evaluated
        int nEvalPoints1D    =  9;
        int nTotEvalPoints   =  nEvalPoints1D*(nEvalPoints1D+1)/2;

        // Declare variables to hold evaluation points
        Array<OneD,NekDouble> evalPointsXi1(nTotEvalPoints);
        Array<OneD,NekDouble> evalPointsXi2(nTotEvalPoints);

        // Calculate the equidistant evaluation points
        cnt = 0;
        for(i = 0; i < nEvalPoints1D; i++)
        {
            for(j = 0 ; j <= nEvalPoints1D-1-i; j++)
            {
                evalPointsXi1[cnt]   = 2.0*i/(nEvalPoints1D-1.0)-1.0;
                evalPointsXi2[cnt++] = 2.0*j/(nEvalPoints1D-1.0)-1.0;
            }
        }

        // Declare variables to hold the values of the basis functions
        // at the evaluation points
        NekVector<NekDouble> orthoBasisEvaluated(nTotModes);
        NekVector<NekDouble> lagBasisEvaluatedtmp(nTotModes);
        NekMatrix<NekDouble> lagBasisEvaluated(nTotModes,nTotEvalPoints);

        // Evaluate the Lagrange polynomials at the evaluation points
        for(i = 0; i < nTotEvalPoints; i++)
        {                   
            coord[0] = evalPointsXi1[i];
            coord[1] = evalPointsXi2[i];

            // Interpolate mode j to the evaluation points i
            for(j = 0; j < nTotModes; j++)
            {
                // Fill the physical space of triExpansion with orthogonal mode i
                triExpansion->FillMode(j,triExpansion->UpdatePhys());

                // Evaluate the ortogonal mode j at the evaluation point j
                orthoBasisEvaluated[j] = triExpansion->PhysEvaluate(coord);                
            }

            // Multiply by the inverse of the Vandermonde matrix to 
            // get the values of the Lagrange polynomial
            lagBasisEvaluatedtmp = V*orthoBasisEvaluated;

            // Store these values in the variable lagBasisEvaluated
            for(j = 0; j < nTotModes; j++)
            {
                lagBasisEvaluated(j,i) = lagBasisEvaluatedtmp[j];
            }
        }

        // Exercise 3(c): Calculate the Lebesgue function for these equidistant nodal points 
        Array<OneD, NekDouble> lebesgueFunction(nTotEvalPoints,0.0);
        for(i = 0; i < nTotEvalPoints; i++)
        {
            for(j = 0; j < nTotModes; j++)
            {
                lebesgueFunction[i] += fabs(lagBasisEvaluated(j,i));
            }
        }

        // Display the output
        cout << "The values of the different Lagrange polynomials evaluated at the " << nTotEvalPoints << endl; 
        cout << "equidistant points are now contained in the variable:" << endl;
        cout << "  - lagBasisEvaluated" << endl;
        cout << "The values of the Lebesgue function evaluated at the same points" << endl; 
        cout << "are contained in the variable:" << endl;
        cout << "  - lebesgueFunction" << endl;
        cout << "(If desired, write your own routine to cast these data in such a format" << endl;
        cout<< "that it can be read by your post-processing program of choice.)" << endl << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        // Input: specify the details of the expansion
        int nModes1D        =  5;
        int nQuadPointsDir1 =  6;
        int nQuadPointsDir2 =  6;

        LibUtilities::PointsType nodalPointsType    = LibUtilities::eNodalTriEvenlySpaced;
        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eOrtho_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eOrtho_B;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModes1D,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModes1D,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard triangular region
        StdRegions::StdNodalTriExpSharedPtr triExpansion = 
            MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,nodalPointsType);

        // Exercise 3(c): Calculate the Lebesgue function for these equidistant nodal points 
        // Specify the number of points where the Lebesgue function 
        // should be evaluated
        int nEvalPoints1D    =  9;
        int nTotEvalPoints   =  nEvalPoints1D*(nEvalPoints1D+1)/2;
        int nTotModes        =  triExpansion->GetNcoeffs();

        // Declare variables to hold evaluation points
        Array<OneD,NekDouble> evalPointsXi1(nTotEvalPoints);
        Array<OneD,NekDouble> evalPointsXi2(nTotEvalPoints);

        // Calculate the equidistant evaluation points
        int i,j;
        int cnt = 0;
        for(i = 0; i < nEvalPoints1D; i++)
        {
            for(j = 0 ; j <= nEvalPoints1D-1-i; j++)
            {
                evalPointsXi1[cnt]   = 2.0*i/(nEvalPoints1D-1.0)-1.0;
                evalPointsXi2[cnt++] = 2.0*j/(nEvalPoints1D-1.0)-1.0;
            }
        }

        Array<OneD, NekDouble> lebesgueFunction(nTotEvalPoints,0.0);
        Array<OneD, NekDouble> coord(2);

        for(i = 0; i < nTotEvalPoints; i++)
        {                   
            coord[0] = evalPointsXi1[i];
            coord[1] = evalPointsXi2[i];
            
            // Interpolate mode j to the evaluation points i
            for(j = 0; j < nTotModes; j++)
            {
                // Fill the physical space of triExpansion with orthogonal mode i
                triExpansion->FillMode(j,triExpansion->UpdatePhys());

                // Evaluate the Lagrange polynomial j at the evaluation point j
                lebesgueFunction[i] += fabs(triExpansion->PhysEvaluate(coord));                
            }
        }

        // Display the output
        // We will use one the Nektar++ output formats to visualise the output.
        // The data of every Lagrange polynomial i are written to the file Chapter3Exercise3_modei.pos
        // which can be opened using the program Gmsh. Given the values of the coefficients of 
        // the expansion, this program then plots the expansion in a high-order fashion.
        for(int i = 0; i < nTotModes; i++)
        {
            stringstream fileName;
            fileName << "Chapter3Exercise3a_mode" << i << ".pos";
            ofstream out((fileName.str()).data());

            Array<OneD,NekDouble> mode(nTotModes,0.0);
            mode[i]=1.0;
            triExpansion->SetCoeffs(mode);
            triExpansion->WriteToFile(out,eGmsh);
        }
        cout << "The different Lagrange Polynomials have been written to the files Chapter3Exercise3a_modei.pos" << endl;
        cout << "Use the program Gmsh to open these files and see a plot of these basis functions." << endl << endl;  
        cout << "The values of the Lebesgue function evaluated at " << nTotEvalPoints << " equidistant points" << endl; 
        cout << "are contained in the variable:" << endl;
        cout << "  - lebesgueFunction" << endl;
        cout << "(If desired, write your own routine to cast these data in such a format" << endl;
        cout<< "that it can be read by your post-processing program of choice.)" << endl;
    }

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 3(b+c)-" <<endl;
    cout << "-------------------" <<endl;

    // The implementation of exercise 3(b) is completely analogue to the implementation of exercise 3(a).
    // The only difference is the use of electrostatic rather than equidistant points for the nodes of the
    // Lagrange polynomials.

    // Before solving the exercises, the Electrostatic points will be calculated.
    int nModesElectrostaticPoints    =  5;
    int nTotModesElectrostaticPoints =  nModesElectrostaticPoints*(nModesElectrostaticPoints+1)/2;
    LibUtilities::PointsType nodalPointsType = LibUtilities::eNodalTriElec;
    LibUtilities::PointsKey electroStaticPointsKey(nModesElectrostaticPoints,nodalPointsType);

    Array<OneD, NekDouble> electroStaticPoints_xi1(nTotModesElectrostaticPoints);
    Array<OneD, NekDouble> electroStaticPoints_xi2(nTotModesElectrostaticPoints);
    LibUtilities::PointsManager()[electroStaticPointsKey]->
        GetPoints(electroStaticPoints_xi1,electroStaticPoints_xi2);    
    
    cout << "METHOD 1: PolyLib" << endl;
    {
        // Input: specify the details of the expansion
        int nModes1D        =  5;
        int nTotModes       =  nModes1D*(nModes1D+1)/2;

        // Declare variables to hold nodal points
        double* nodalPointsXi1  = new double[nTotModes];
        double* nodalPointsXi2  = new double[nTotModes];
        double* nodalPointsEta1 = new double[nTotModes];
        double* nodalPointsEta2 = new double[nTotModes];

        // Calculate the electrostatic nodal points
        // Also calculate the nodal points in eta1-eta2 coordinates
        int i,j;
        const double tol = 1e-12;
        for(i = 0; i < nTotModes; i++)
        {           
            nodalPointsXi1[i] = electroStaticPoints_xi1[i];
            nodalPointsXi2[i] = electroStaticPoints_xi2[i];

            if( fabs(nodalPointsXi2[i]-1.0)<tol )
            {
                nodalPointsEta1[i]=0.0;
            }
            else
            {
                nodalPointsEta1[i] = 2*(1+nodalPointsXi1[i])/(1-nodalPointsXi2[i]);
            }
            nodalPointsEta2[i] = nodalPointsXi2[i];
        }

        // Allocate memory for the generalised Vandermonde matrix
        double* V = new double [nTotModes*nTotModes];

        // Calculate the Vandermonde matrix
        int p,q,m;
        int cnt = 0;
        double psi_p[1];
        double psi_pq[1];

        for(i = 0; i < nTotModes; i++)
        {            
            for(p = 0; p < nModes1D; p++)
            {
                for(q = 0 ; q <= nModes1D-1-p; q++)
                {                    
                    Polylib::jacobfd(1, &(nodalPointsEta1[i]), psi_p,  NULL, p, 0.0, 0.0);
                    Polylib::jacobfd(1, &(nodalPointsEta2[i]), psi_pq, NULL, q, 2.0*p + 1.0, 0.0);
                    
                    V[cnt] = psi_p[0] * pow((1.0-nodalPointsEta2[i])/2.0 , p) * psi_pq[0];
                    cnt++;
                }
            }
        }

        // Invert the Vandermonde matrix
        int*   pivot = new int[nTotModes];
        double* work = new double[nTotModes];
        int     info = 0;
                    
        Lapack::Dgetrf(nTotModes, nTotModes, V, nTotModes, pivot, info);
        
        if( info < 0 )
        {
            cerr << "ERROR: The " << -info << "th parameter had an illegal parameter for dgetrf" << endl;
            exit(1);
        }
        else if( info > 0 )
        {
            cerr << "ERROR: Element V[" << info << "] is 0 from dgetrf"<<endl;
            exit(1);
        }   
        
        Lapack::Dgetri(nTotModes, V, nTotModes, pivot, work, nTotModes, info);
                    
        if( info < 0 )
        {
            cerr << "ERROR: The " << -info << "th parameter had an illegal parameter for dgetri" << endl;
            exit(1);
        }
        else if( info > 0 )
        {
            cerr << "ERROR: Element V[" << info << "] is 0 from dgetri"<<endl;
            exit(1);
        }   


        // specify the number of points where the Lagrange polynomials
        // should be evaluated
        int nEvalPoints1D    =  20;
        int nTotEvalPoints   =  nEvalPoints1D*(nEvalPoints1D+1)/2;

        // Declare variables to hold evaluation points
        double* evalPointsXi1  = new double[nTotEvalPoints];
        double* evalPointsXi2  = new double[nTotEvalPoints];
        double* evalPointsEta1 = new double[nTotEvalPoints];
        double* evalPointsEta2 = new double[nTotEvalPoints];

        // Calculate the equidistant evaluation points
        // Also calculate the evaluation points in eta1-eta2 coordinates
        cnt = 0;
        for(i = 0; i < nEvalPoints1D; i++)
        {
            for(j = 0 ; j <= nEvalPoints1D-1-i; j++)
            {
                evalPointsXi1[cnt] = 2.0*i/(nEvalPoints1D-1.0)-1.0;
                evalPointsXi2[cnt] = 2.0*j/(nEvalPoints1D-1.0)-1.0;


                if( fabs(evalPointsXi2[cnt]-1.0)<tol )
                {
                    evalPointsEta1[cnt]=0.0;
                }
                else
                {
                    evalPointsEta1[cnt] = 2*(1+evalPointsXi1[cnt])/(1-evalPointsXi2[cnt]);
                }
                evalPointsEta2[cnt] = evalPointsXi2[cnt];

                cnt++;
            }
        }


        // Declare variables to hold the values of the basis functions
        // at the evaluation points
        double* orthoBasisEvaluated  = new double[nTotModes];
        double* lagBasisEvaluated    = new double[nTotModes*nTotEvalPoints];

        // Evaluate the Lagrange polynomials at the evaluation points
        for(i = 0; i < nTotEvalPoints; i++)
        {         
            // Step 1: Evaluate the Orthogonal Basis at the evaluation points   
            cnt=0;
            for(p = 0; p < nModes1D; p++)
            {
                for(q = 0 ; q <= nModes1D-1-p; q++)
                {                    
                    Polylib::jacobfd(1, &(evalPointsEta1[i]), psi_p,  NULL, p, 0.0, 0.0);
                    Polylib::jacobfd(1, &(evalPointsEta2[i]), psi_pq, NULL, q, 2.0*p + 1.0, 0.0);
                    
                    orthoBasisEvaluated[cnt++] = psi_p[0] * pow((1.0-evalPointsEta2[i])/2.0 , p) * psi_pq[0];
                }
            }

            // Step 2: multiply by the inverse of the Vandermonde matrix to 
            // get the values of the Lagrange polynomial.
            // This matrix-vector multiplication is done by the corresponding
            // BLAS call. 
            Blas::Dgemv('N',nTotModes,nTotModes,1.0,V,nTotModes,orthoBasisEvaluated,1,0.0,lagBasisEvaluated+i,nTotEvalPoints);
        }

        // Exercise 3(c): Calculate the Lebesgue function for these electrostatic nodal points 
        double* lebesgueFunction  = new double[nTotEvalPoints];
        for(i = 0; i < nTotEvalPoints; i++)
        {
            lebesgueFunction[i] = 0.0;

            for(j = 0; j < nTotModes; j++)
            {
                lebesgueFunction[i] += fabs(lagBasisEvaluated[i+nTotEvalPoints*j]);
            }
        }

        // Display the output
        cout << "The values of the different Lagrange polynomials evaluated at the " << nTotEvalPoints << endl; 
        cout << "equidistant points are now contained in the variable:" << endl;
        cout << "  - lagBasisEvaluated" << endl;
        cout << "The values of the Lebesgue function evaluated at the same points" << endl; 
        cout << "are contained in the variable:" << endl;
        cout << "  - lebesgueFunction" << endl;
        cout << "(If desired, write your own routine to cast these data in such a format" << endl;
        cout<< "that it can be read by your post-processing program of choice.)" << endl << endl;

        // Deallocate the dynamic memory
        delete[] nodalPointsXi1;
        delete[] nodalPointsXi2;
        delete[] nodalPointsEta1;
        delete[] nodalPointsEta2;
        delete[] V;
        delete[] pivot;
        delete[] work;
        delete[] evalPointsXi1;
        delete[] evalPointsXi2;
        delete[] evalPointsEta1;
        delete[] evalPointsEta2;
        delete[] orthoBasisEvaluated;
        delete[] lagBasisEvaluated;
        delete[] lebesgueFunction;

    }

    cout << "METHOD 2: low-level Nektar++" << endl;
    {
        // Input: specify the details of the expansion
        int nModes1D        =  5;
        int nQuadPointsDir1 =  6;
        int nQuadPointsDir2 =  6;
        int nTotModes       =  nModes1D*(nModes1D+1)/2;

        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eOrtho_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eOrtho_B;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModes1D,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModes1D,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard triangular region
        StdRegions::StdTriExpSharedPtr triExpansion = 
            MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2);

        // Declare variables to hold nodal points
        Array<OneD, NekDouble> nodalPointsXi1(nTotModes);
        Array<OneD, NekDouble> nodalPointsXi2(nTotModes);

        // Calculate the equidistant nodal points
        // Also calcualte the nodal points in eta1-eta2 coordinates
        int i,j;
        int cnt = 0;
        for(i = 0; i < nTotModes; i++)
        {
            nodalPointsXi1[i] = electroStaticPoints_xi1[i];
            nodalPointsXi2[i] = electroStaticPoints_xi2[i];
        }

        // Allocate memory for the generalised Vandermonde matrix
        NekMatrix<NekDouble> V(nTotModes,nTotModes);

        // Calculate the Vandermonde matrix
        Array<OneD, NekDouble> coord(2);
        for(i = 0; i < nTotModes; i++)
        {
            // Fill the physical space of triExpansion with orthogonal mode i
            triExpansion->FillMode(i,triExpansion->UpdatePhys());

            // Interpolate mode i to the nodal points j
            for(j = 0; j < nTotModes; j++)
            {
                    coord[0] = nodalPointsXi1[j];
                    coord[1] = nodalPointsXi2[j];
                    V(i,j) = triExpansion->PhysEvaluate(coord);                
            }
        }

        // Invert the Vandermonde matrix
        V.Invert();

        // specify the number of points where the Lagrange polynomials
        // should be evaluated
        int nEvalPoints1D    =  9;
        int nTotEvalPoints   =  nEvalPoints1D*(nEvalPoints1D+1)/2;

        // Declare variables to hold evaluation points
        Array<OneD,NekDouble> evalPointsXi1(nTotEvalPoints);
        Array<OneD,NekDouble> evalPointsXi2(nTotEvalPoints);

        // Calculate the equidistant evaluation points
        cnt = 0;
        for(i = 0; i < nEvalPoints1D; i++)
        {
            for(j = 0 ; j <= nEvalPoints1D-1-i; j++)
            {
                evalPointsXi1[cnt]   = 2.0*i/(nEvalPoints1D-1.0)-1.0;
                evalPointsXi2[cnt++] = 2.0*j/(nEvalPoints1D-1.0)-1.0;
            }
        }

        // Declare variables to hold the values of the basis functions
        // at the evaluation points
        NekVector<NekDouble> orthoBasisEvaluated(nTotModes);
        NekVector<NekDouble> lagBasisEvaluatedtmp(nTotModes);
        NekMatrix<NekDouble> lagBasisEvaluated(nTotModes,nTotEvalPoints);

        // Evaluate the Lagrange polynomials at the evaluation points
        for(i = 0; i < nTotEvalPoints; i++)
        {                   
            coord[0] = evalPointsXi1[i];
            coord[1] = evalPointsXi2[i];

            // Interpolate mode j to the evaluation points i
            for(j = 0; j < nTotModes; j++)
            {
                // Fill the physical space of triExpansion with orthogonal mode i
                triExpansion->FillMode(j,triExpansion->UpdatePhys());

                // Evaluate the ortogonal mode j at the evaluation point j
                orthoBasisEvaluated[j] = triExpansion->PhysEvaluate(coord);                
            }

            // Multiply by the inverse of the Vandermonde matrix to 
            // get the values of the Lagrange polynomial
            lagBasisEvaluatedtmp = V*orthoBasisEvaluated;

            // Store these values in the variable lagBasisEvaluated
            for(j = 0; j < nTotModes; j++)
            {
                lagBasisEvaluated(j,i) = lagBasisEvaluatedtmp[j];
            }
        }

       // Exercise 3(c): Calculate the Lebesgue function for these electrostatic nodal points 
        Array<OneD, NekDouble> lebesgueFunction(nTotEvalPoints,0.0);
        for(i = 0; i < nTotEvalPoints; i++)
        {
            for(j = 0; j < nTotModes; j++)
            {
                lebesgueFunction[i] += fabs(lagBasisEvaluated(j,i));
            }
        }

        // Display the output
        cout << "The values of the different Lagrange polynomials evaluated at the " << nTotEvalPoints << endl; 
        cout << "equidistant points are now contained in the variable:" << endl;
        cout << "  - lagBasisEvaluated" << endl;
        cout << "The values of the Lebesgue function evaluated at the same points" << endl; 
        cout << "are contained in the variable:" << endl;
        cout << "  - lebesgueFunction" << endl;
        cout << "(If desired, write your own routine to cast these data in such a format" << endl;
        cout<< "that it can be read by your post-processing program of choice.)" << endl << endl;
    }

    cout << "METHOD 3: Nektar++" << endl;
    {
        // Input: specify the details of the expansion
        int nModes1D        =  5;
        int nQuadPointsDir1 =  6;
        int nQuadPointsDir2 =  6;

        LibUtilities::PointsType nodalPointsType    = LibUtilities::eNodalTriElec;
        LibUtilities::PointsType quadPointsTypeDir1 = LibUtilities::eGaussLobattoLegendre;
        LibUtilities::PointsType quadPointsTypeDir2 = LibUtilities::eGaussRadauMAlpha1Beta0;

        LibUtilities::BasisType basisTypeDir1 = LibUtilities::eOrtho_A;
        LibUtilities::BasisType basisTypeDir2 = LibUtilities::eOrtho_B;

        // Declare the PointsKey's which uniquely define the quadrature points
        const LibUtilities::PointsKey quadPointsKeyDir1(nQuadPointsDir1, quadPointsTypeDir1);
        const LibUtilities::PointsKey quadPointsKeyDir2(nQuadPointsDir2, quadPointsTypeDir2);

        // Declare the BasisKey's which uniquely define (one-dimensional) expansion bases
        const LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,nModes1D,quadPointsKeyDir1);
        const LibUtilities::BasisKey basisKeyDir2(basisTypeDir2,nModes1D,quadPointsKeyDir2);

        // Using these keys, now define an spectral/hp expansion of corresponding type on a 
        // standard triangular region
        StdRegions::StdNodalTriExpSharedPtr triExpansion = 
            MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(basisKeyDir1,basisKeyDir2,nodalPointsType);

        // Exercise 3(c): Calculate the Lebesgue function for these electrostatic nodal points 
        // Specify the number of points where the Lebesgue function 
        // should be evaluated
        int nEvalPoints1D    =  9;
        int nTotEvalPoints   =  nEvalPoints1D*(nEvalPoints1D+1)/2;
        int nTotModes        =  triExpansion->GetNcoeffs();

        // Declare variables to hold evaluation points
        Array<OneD,NekDouble> evalPointsXi1(nTotEvalPoints);
        Array<OneD,NekDouble> evalPointsXi2(nTotEvalPoints);

        // Calculate the equidistant evaluation points
        int i,j;
        int cnt = 0;
        for(i = 0; i < nEvalPoints1D; i++)
        {
            for(j = 0 ; j <= nEvalPoints1D-1-i; j++)
            {
                evalPointsXi1[cnt]   = 2.0*i/(nEvalPoints1D-1.0)-1.0;
                evalPointsXi2[cnt++] = 2.0*j/(nEvalPoints1D-1.0)-1.0;
            }
        }

        Array<OneD, NekDouble> lebesgueFunction(nTotEvalPoints,0.0);
        Array<OneD, NekDouble> coord(2);

        for(i = 0; i < nTotEvalPoints; i++)
        {                   
            coord[0] = evalPointsXi1[i];
            coord[1] = evalPointsXi2[i];
            
            // Interpolate mode j to the evaluation points i
            for(j = 0; j < nTotModes; j++)
            {
                // Fill the physical space of triExpansion with orthogonal mode i
                triExpansion->FillMode(j,triExpansion->UpdatePhys());

                // Evaluate the Lagrange polynomial j at the evaluation point j
                lebesgueFunction[i] += fabs(triExpansion->PhysEvaluate(coord));                
            }
        }

        // Display the output
        // We will use one the Nektar++ output formats to visualise the output.
        // The data of every Lagrange polynomial i are written to the file Chapter3Exercise3_modei.pos
        // which can be opened using the program Gmsh. Given the values of the coefficients of 
        // the expansion, this program then plots the expansion in a high-order fashion.
        for(int i = 0; i < nTotModes; i++)
        {
            stringstream fileName;
            fileName << "Chapter3Exercise3b_mode" << i << ".pos";
            ofstream out((fileName.str()).data());

            Array<OneD,NekDouble> mode(nTotModes,0.0);
            mode[i]=1.0;
            triExpansion->SetCoeffs(mode);
            triExpansion->WriteToFile(out,eGmsh);
        }
        cout << "The different Lagrange Polynomials have been written to the files Chapter3Exercise3b_modei.pos" << endl;
        cout << "Use the program Gmsh to open these files and see a plot of these basis functions." << endl << endl;  
        cout << "The values of the Lebesgue function evaluated at " << nTotEvalPoints << " equidistant points" << endl; 
        cout << "are contained in the variable:" << endl;
        cout << "  - lebesgueFunction" << endl;
        cout << "(If desired, write your own routine to cast these data in such a format" << endl;
        cout<< "that it can be read by your post-processing program of choice.)" << endl; 
    }
}

