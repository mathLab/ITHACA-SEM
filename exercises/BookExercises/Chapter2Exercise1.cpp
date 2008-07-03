#include <cstdio>
#include <cstdlib>
#include <sys/utime.h>

#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Polylib/Polylib.h>

using namespace Nektar;
#define M_PI 2.0*acos(0.0)

void main(int argc, char *argv[])
{
	cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 2 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(a) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: PolyLib" << endl;
	{
		
		// DEclare variable that hold the number of quadrature points
		int nQuadPoints;

		for(nQuadPoints = 4; nQuadPoints <= 6; nQuadPoints++)
		{
		
			// Decalare variables to hold the quadrature zeros and weights
			double* quadZeros = new double[nQuadPoints];
			double* quadWeights = new double[nQuadPoints];

			// Calculate the GLL-quadrature zeros and weights using the corresponding
			// Polylib routine
			const double alpha = 0.0;
			const double beta = 0.0;
			Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);
			
			// Apply the Gaussian quadrature technique
			int i;
			double result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				result += quadWeights[i] * pow(quadZeros[i],6);
			}
			


			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
			// Deallocate the dynamic memory
			delete[] quadZeros;
			delete[] quadWeights;
			
		}
		cout << endl;
		
	}

	cout << "METHOD2: low-level Nektar++" << endl;
	{
		// Declare variable the holds the number of quadrature points
		int nQuadPoints;
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 4; nQuadPoints <= 6; nQuadPoints++)
		{
			//Decalre variables to hold the quadrature zeros and weights
			Array<OneD, NekDouble> quadZeros;
			Array<OneD, NekDouble> quadWeights;
			
			// Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
			// Step 1: Decalre a PointsKey which uniquele defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			//Step 2:: Using this key, the quadrature zeros and weights can now be retrieved through
			// The PointsManager
			quadZeros = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			quadWeights = LibUtilities::PointsManager()[quadPointsKey]->GetW();

			// Apply the Guassian quadrature technique
			int i;
			NekDouble result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				result += quadWeights[i] * pow(quadZeros[i],6);
			}
			

			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
		}
		cout << endl;

	}

	cout << "METHOD 3: Nektar++" << endl;
	{
		
		// Declare variable the holds the number of quadrature zeros
		int nQuadPoints;
		
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 4; nQuadPoints <= 6; nQuadPoints++)
		{
			// Declare a PointsKey which uniquely defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			// Declare a BasisKey which uniquely defines the (discrete) one-dimensional
			// spectral/hp expansion basis
			// As we will use the basis only to evaluate an integral, the type and order of the basis
			// can be chosen arbitrarily
			const LibUtilities::BasisKey basisKey(LibUtilities::eModified_A,2,quadPointsKey);
			
			// Using these keys, now define an spectral/hp expansion of corresponding type on a
			// standard region
			StdRegions::StdSegExpSharedPtr segExpansion = 
				MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);
			
			// Calculate the coordinates xi of all the quadrature points of this (discrete)
			// expansion
			int nTotQuadPoints = segExpansion->GetTotPoints();
			Array<OneD,NekDouble> xi(nTotQuadPoints);
			segExpansion->GetCoords(xi);

			// Calculate the values of the function to be integrated
			int i;
			Array<OneD,NekDouble> integrand(nTotQuadPoints);
			for(i = 0; i < nTotQuadPoints; i++)
			{
				integrand[i] = pow(xi[i],6);
			}
			
			// Do the numerical integration by calling the ecorresponding Nektar++ routine of
			// the StdExpansion class
			NekDouble result = segExpansion->Integral(integrand);
			
			// Display the outputs
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
		}
		cout << endl;
	}


	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(b) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: PolyLib" << endl;
	{
		// Specify the coordinates of the element
		double x1 = -1.0;
		double x2 = 2.0;
		
		// Declare variable that holds the number of quadrature zeros
		int nQuadPoints;

		for(nQuadPoints = 4; nQuadPoints <= 6; nQuadPoints++)
		{
			// Decalare variables to hold the quadrature zeros and weights
			double* quadZeros = new double[nQuadPoints];
			double* quadWeights = new double[nQuadPoints];

			// Calculate the GLL-quadrature zeros and weights using the corresponding
			// Polylib routine
			const double alpha = 0.0;
			const double beta = 0.0;
			Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);
			
			// Analytically evaluate the Jacobian.
			double jacobian;
			jacobian = 0.5 * (-x1) + 0.5 * x2;

			// Apply the Gaussian quadrature technique
			int i;
			double result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = -1.0 * (1-quadZeros[i])/2 + 2.0 * (1+quadZeros[i])/2;
				
				result += quadWeights[i] * pow(x,6) * jacobian;
			}
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
			// Deallocate the dynamic memory
			delete[] quadZeros;
			delete[] quadWeights;
			
		}
		cout << endl;
	}

	cout << "METHOD 2: low-level Nektar++" << endl;
    {
		// Specify the coordinates of the element
		NekDouble x1 = -1.0;
		NekDouble x2 = 2.0;
		
		// Declare variable the holds the number of quadrature zeros
		int nQuadPoints;
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 4; nQuadPoints <= 6; nQuadPoints++)
		{
			//Decalre variables to hold the quadrature zeros and weights
			Array<OneD, NekDouble> quadZeros;
			Array<OneD, NekDouble> quadWeights;
			
			// Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
			// Step 1: Decalre a PointsKey which uniquele defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			//Step 2:: Using this key, the quadrature zeros and weights can now be retrieved through
			// The PointsManager
			quadZeros = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			quadWeights = LibUtilities::PointsManager()[quadPointsKey]->GetW();

			// Analytically evaluate the Jacobian.
			NekDouble jacobian;
			jacobian = 0.5 * (-x1) + 0.5 * x2;

			// Apply the Guassian quadrature technique
			int i;
			NekDouble result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = -1.0 * (1-quadZeros[i])/2 + 2.0 * (1+quadZeros[i])/2;
				
				result += quadWeights[i] * pow(x,6) * jacobian;
			}
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
		}
		cout << endl;
	}

	cout << "METHOD 3: Nektar++" << endl;
	{
		// Specify the coordinates of the element
		NekDouble x1 = -1.0;
		NekDouble x2 = 2.0;
		
		// Declare variable the holds the number of quadrature zeros
		int nQuadPoints;
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 4; nQuadPoints <= 6; nQuadPoints++)
		{
			// Declare a PointsKey which uniquely defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			// Declare a BasisKey which uniquely defines the (discrete) one-dimensional
			// spectral/hp expansion basis
			// As we will use the basis only to evaluate an integral, the type and order of the basis
			// can be chosen arbitrarily
			const LibUtilities::BasisKey basisKey(LibUtilities::eModified_A,2,quadPointsKey);
			
			// Using these keys, now define an spectral/hp expansion of corresponding type on a
			// standard region
			StdRegions::StdSegExpSharedPtr segExpansion = 
				MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);
			
			// Calculate the coordinates xi of all the quadrature points of this (discrete)
			// expansion
			int nTotQuadPoints = segExpansion->GetTotPoints();
			Array<OneD,NekDouble> xi(nTotQuadPoints);
			segExpansion->GetCoords(xi);

			// Analytically evaluate the Jacobian.
			NekDouble jacobian;
			jacobian = 0.5 * (-x1) + 0.5 * x2;
			
			// Calculate the values of the function to be integrated
			int i;
			Array<OneD,NekDouble> integrand(nTotQuadPoints);
			for(i = 0; i < nTotQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = -1.0 * (1-xi[i])/2 + 2.0 * (1+xi[i])/2;
				
				integrand[i] = pow(x,6);

			}
			
			// Do the numerical integration by calling the ecorresponding Nektar++ routine of
			// the StdExpansion class
			NekDouble result = jacobian * segExpansion->Integral(integrand);
			
			// Display the outputs
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
		}
		cout << endl;
	}

	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 1(c) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: PolyLib" << endl;
	{
		// Specify the coordinate of the element
		double x1 = 0.0;
		double x2 = M_PI/2;

		// Declare variables the hold the number of quadrature zeros
		int nQuadPoints;

		for(nQuadPoints = 2; nQuadPoints <= 8; nQuadPoints++)
		{
			// Decalare variables to hold the quadrature zeros and weights
			double* quadZeros = new double[nQuadPoints];
			double* quadWeights = new double[nQuadPoints];

			// Calculate the GLL-quadrature zeros and weights using the corresponding
			// Polylib routine
			const double alpha = 0.0;
			const double beta = 0.0;
			Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);
			
			// Analytically evaluate the Jacobian.
			double jacobian;
			jacobian = 0.5 * (-x1) + 0.5 * x2;

			// Apply the Gaussian quadrature technique
			int i;
			double result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = x2 * (1+quadZeros[i])/2;

				result += quadWeights[i] * sin(x) * jacobian;
			}
			

			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
			// Deallocate the dynamic memory
			delete[] quadZeros;
			delete[] quadWeights;
			
		}
		cout << endl;
	}

	cout << "METHOD2: low-level Nektar++" << endl;
	{
		// Specify the coordinates of the element
		NekDouble x1 = 0.0;
		NekDouble x2 = M_PI/2;

		// Decalre variable that holds the number of quadrature zeros
		int nQuadPoints;

		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 2; nQuadPoints <= 8; nQuadPoints++)
		{
			//Decalre variables to hold the quadrature zeros and weights
			Array<OneD, NekDouble> quadZeros;
			Array<OneD, NekDouble> quadWeights;
			
			// Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
			// Step 1: Decalre a PointsKey which uniquele defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			//Step 2:: Using this key, the quadrature zeros and weights can now be retrieved through
			// The PointsManager
			quadZeros = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			quadWeights = LibUtilities::PointsManager()[quadPointsKey]->GetW();

			// Analytically evaluate the Jacobian.
			NekDouble jacobian;
			jacobian = 0.5 * (-x1) + 0.5 * x2;
			
			// Apply the Guassian quadrature technique
			int i;
			NekDouble result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = M_PI * (1+quadZeros[i])/2;

				result += quadWeights[i] * sin(x) * jacobian;
			}
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
		}
		cout << endl;
	}

	cout << "METHOD 3: Nektar++" << endl;
	{
		// Specify the coordinate of the element
		NekDouble x1 = 0.0;
		NekDouble x2 = M_PI/2;
		
		// Declare variable the holls the number od quadratue points		
		int nQuadPoints;
		
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 2; nQuadPoints <= 8; nQuadPoints++)
		{
			// Declare a PointsKey which uniquely defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			// Declare a BasisKey which uniquely defines the (discrete) one-dimensional
			// spectral/hp expansion basis
			// As we will use the basis only to evaluate an integral, the type and order of the basis
			// can be chosen arbitrarily
			const LibUtilities::BasisKey basisKey(LibUtilities::eModified_A,2,quadPointsKey);
			
			// Using these keys, now define an spectral/hp expansion of correpoinding type on a
			// standard region
			StdRegions::StdSegExpSharedPtr segExpansion = 
				MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);
			
			// Calculate the coordinates xi of all the quadrature points of this (discrete)
			// expansion
			int nTotQuadPoints = segExpansion->GetTotPoints();
			Array<OneD,NekDouble> xi(nTotQuadPoints);
			segExpansion->GetCoords(xi);

			// Analytically evaluate the Jacobian.
			NekDouble jacobian;
			jacobian = 0.5 * (-x1) + 0.5 * x2;

			// Calculate the values of the function to be integrated
			int i;
			Array<OneD,NekDouble> integrand(nTotQuadPoints);
			for(i = 0; i < nTotQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x =  M_PI * (1+xi[i])/2;

				integrand[i] = sin(x);
			}
			
			// Do the numerical integration by calling the ecorresponding Nektar++ routine of
			// the StdExpansion class
			NekDouble result = jacobian * segExpansion->Integral(integrand);
			
			// Display the outputs
			cout << "Result for Q = " << nQuadPoints << ":" << result << endl;
			
		}
		cout << endl;
	}

}
