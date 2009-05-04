#include <cstdio>
#include <cstdlib>

#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Polylib/Polylib.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
	cout << "=====================================================" <<endl;
    cout << "================ EXERCISES CHAPTER 2 ================" <<endl;
    cout << "=====================================================" <<endl;

    cout << "-------------------" <<endl;
    cout << "-- EXERCISE 2(a) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: PolyLib" << endl;
	{
		
		// Declare variable the holds the 
		int nQuadPoints;

		for(nQuadPoints = 7; nQuadPoints <= 9; nQuadPoints++)
		{
			// Decalare variables to hold the quadrature zeros and weights
			// weights are not used in computing the derivative, we define them to 
			// provide the zwglj routine with enough inputs
			double* quadZeros = new double[nQuadPoints];
			double* quadWeights = new double[nQuadPoints];

			// Calculate the GLL-quadrature zeros and weights using the corresponding
			// Polylib routine
			const double alpha = 0.0;
			const double beta = 0.0;
			Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);
			
			// Decalare variable that holds the values of the derivative matrix
			double* derivative = new double[nQuadPoints*nQuadPoints];
			
			// Calculate the derivative matrix using the quadrature zeros with
			// the corresponding Polylib routine
			Polylib::Dglj(derivative,quadZeros,nQuadPoints,alpha,beta);
			
			// Calculate the derivative of each quadrature point
			int i;
			int j;
			double* result = new double[nQuadPoints];
								
			for(i = 0; i < nQuadPoints; i++)
			{
				result[i] = 0.0;
				for(j = 0; j < nQuadPoints; j++)
				{
					result[i] += derivative[i+nQuadPoints*j]*pow(quadZeros[j],7);
				}
			}
			
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << endl;
			for(i = 0; i < nQuadPoints; i++)
			{
				cout << "Quadrature Zero: " << quadZeros[i] << "  Derivative: " << result[i] << endl;
			}
			
			// Deallocate the dynamic memory
			delete[] quadZeros;
			delete[] quadWeights;
			
		}
		cout << endl;

	}

	cout << "METHOD 2: low-level Nektar++" << endl;
	{
		
		// Declare variable that holds the number of quadrature zeros
		int nQuadPoints;
		
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
		
		for(nQuadPoints = 7; nQuadPoints <= 9; nQuadPoints++)
		{
			//Decalre variables to hold the quadrature zeros
			Array<OneD, NekDouble> quadZeros;
			
			// Calculate the GLL-quadrature zeros. This is done in 2 steps.
			// Step 1: Decalre a PointsKey which uniquele defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			//Step 2:: Using this key, the quadrature zeros can now be retrieved through
			// The PointsManager
			quadZeros = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			
			// Calculate the derivative matrix
			LibUtilities::Points<NekDouble>::MatrixSharedPtrType derivMatrix = 
				LibUtilities::PointsManager()[quadPointsKey]->GetD();
			
			// calculate the derivative of each quadrature point
			Array<OneD, NekDouble> result(nQuadPoints);
			int i;
			int j;
			for(i = 0; i < nQuadPoints; i++)
			{
				result[i] = 0.0;
				for(j = 0; j < nQuadPoints; j++)
				{
					result[i] += derivMatrix->GetValue(i,j)*pow(quadZeros[j],7);
				}
			}
			

			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << endl;
			for(i = 0; i < nQuadPoints; i++)
			{
				cout << "Quadrature Zero: " << quadZeros[i] << "  Derivative: " << result[i] << endl;
			}
			
		}
		cout << endl;

	}
	
	cout << "METHOD 3: Nektar++" << endl;
	{
		
		// Declare variable that holds the number of quadrature zeros
		int nQuadPoints;
		
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 7; nQuadPoints <= 9; nQuadPoints++)
		{
			// Declare a PointsKey which uniquely defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			// Declare a BasisKey which uniquely defines the (discrete) one-dimensional
			// spectral/hp expansion basis
			// As we will use the basis only to evaluate a derivative, the type and order of the basis
			// can be chosen arbitrarily
			const LibUtilities::BasisKey basisKey(LibUtilities::eModified_A,2,quadPointsKey);
			
			// Using these keys, now define an spectral/hp expansion of corresponding type on a
			// one-dimensional standard region
			StdRegions::StdSegExpSharedPtr segExpansion = 
				MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(basisKey);
			
			// Calculate the coordinates xi of all the quadrature points of this (discrete)
			// expansion
			int nTotQuadPoints = segExpansion->GetTotPoints();
			Array<OneD,NekDouble> xi(nTotQuadPoints);
			segExpansion->GetCoords(xi);

			// Calculate the values of the function to be differentiated
			int i;
			Array<OneD,NekDouble> function(nTotQuadPoints);
			for(i = 0; i < nTotQuadPoints; i++)
			{
				function[i] = pow(xi[i],7);
			}
			
			// Do the numerical differentiation by calling the corresponding Nektar++ routine of
			// the StdExpansion class

			Array<OneD, NekDouble> result(nQuadPoints);
			
			segExpansion->PhysDeriv(function,result);
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << endl;
			for(i = 0; i < nQuadPoints; i++)
			{
				cout << "Quadrature Zero: " << xi[i] << "  Derivative: " << result[i] << endl;
			}
			
		}
		cout << endl;
	}

	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 2(b) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: PolyLib" << endl;
	{
		// Specify the coordinates of the element
		double x1 = 2.0;
		double x2 = 10.0;

		// Declare variable that hold the number of quadrature zeros
		int nQuadPoints;

		for(nQuadPoints = 7; nQuadPoints <= 9; nQuadPoints++)
		{
			// Decalare variables to hold the quadrature zeros and weights
			// weights are not used in computing the derivative, we define them to 
			// provide the zwglj routine with enough inputs
			double* quadZeros = new double[nQuadPoints];
			double* quadWeights = new double[nQuadPoints];

			// Calculate the GLL-quadrature zeros and weights using the corresponding
			// Polylib routine
			const double alpha = 0.0;
			const double beta = 0.0;
			Polylib::zwglj(quadZeros,quadWeights,nQuadPoints,alpha,beta);
			
			// Decalare variable that holds the values of the derivative matrix
			double* derivative = new double[nQuadPoints*nQuadPoints];
			
			// Calculate the derivative matrix using the quadrature zeros with
			// the corresponding Polylib routine
			Polylib::Dglj(derivative,quadZeros,nQuadPoints,alpha,beta);
			
			// Analytically evaluate the Inverse of Jacobian.
			double invJacobian;
			invJacobian = 1.0/(0.5 * (-x1) + 0.5 * x2);

			cout << invJacobian << endl;

			// Calculate the derivative of each quadrature point
			int i;
			int j;
			double* result = new double[nQuadPoints];
					
			for(i = 0; i < nQuadPoints; i++)
			{
				result[i] = 0.0;
				for(j = 0; j < nQuadPoints; j++)
				{
					// Calculate the local coordinate of quadrature zeros
					double x = x1 * (1-quadZeros[j])/2 + x2 * (1+quadZeros[j])/2;

					result[i] += derivative[i+nQuadPoints*j]* pow(x,7) * invJacobian;
					
				}
			}
						
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << endl;
			for(i = 0; i < nQuadPoints; i++)
			{
				cout << "Quadrature Zero: " << x1 * (1-quadZeros[i])/2 + x2 * (1+quadZeros[i])/2
						<< "  Derivative: " << result[i] << endl;
			}
			
			// Deallocate the dynamic memory
			delete[] quadZeros;
			delete[] quadWeights;
			
		}
		cout << endl;

	}

	cout << "METHOD 2: low-level Nektar++" << endl;
	{
		// Specify the coordinates of the element
		NekDouble x1 = 2.0;
		NekDouble x2 = 10.0;

		// Declare variable to hold the number of quadrature zeros
		int nQuadPoints;

		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;
		
		for(nQuadPoints = 7; nQuadPoints <= 9; nQuadPoints++)
		{
			//Decalre variables to hold the quadrature zeros
			Array<OneD, NekDouble> quadZeros;
			
			// Calculate the GLL-quadrature zeros. This is done in 2 steps.
			// Step 1: Decalre a PointsKey which uniquele defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			//Step 2:: Using this key, the quadrature zeros can now be retrieved through
			// The PointsManager
			quadZeros = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			
			// Calculate the deriative matrix
			LibUtilities::Points<NekDouble>::MatrixSharedPtrType derivMatrix = 
				LibUtilities::PointsManager()[quadPointsKey]->GetD();
			
			// Analytically evaluate the inverse of the Jacobian.
			NekDouble invJacobian;
			invJacobian = 1.0/(0.5 * (-x1) + 0.5 * x2);
			
			// calculate the derivative of each quadrature point
			Array<OneD, NekDouble> result(nQuadPoints);
			int i;
			int j;
			for(i = 0; i < nQuadPoints; i++)
			{
				result[i] = 0.0;
				for(j = 0; j < nQuadPoints; j++)
				{
					// Calculate the local coordinate of quadrature zeros
					double x = x1 * (1-quadZeros[j])/2 + x2 * (1+quadZeros[j])/2;
					
					result[i] += derivMatrix->GetValue(i,j) * pow(x,7) * invJacobian;
										
				}
			}
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << endl;
			for(i = 0; i < nQuadPoints; i++)
			{
				cout << "Quadrature Zero: " << x1 * (1-quadZeros[i])/2 + x2 * (1+quadZeros[i])/2
											<< "  Derivative: " << result[i] << endl;
			}
			
		}
		cout << endl;

	}

	cout << "METHOD 3: Nektar++" << endl;
	{
		// Specify the coordinates of the element
		NekDouble x1 = 2.0;
		NekDouble x2 = 10.0;

		// Declare variable to hold the number of quadrature zeros
		int nQuadPoints;
		
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 7; nQuadPoints <= 9; nQuadPoints++)
		{
			// Declare a PointsKey which uniquely defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);
			
			// Declare a BasisKey which uniquely defines the (discrete) one-dimensional
			// spectral/hp expansion basis
			// As we will use the basis only to evaluate a derivative, the type and order of the basis
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

			// Analytically evaluate the inverse of the Jacobian.
			NekDouble invJacobian;
			invJacobian = 1.0/(0.5 * (-x1) + 0.5 * x2);

			// Calculate the values of the function to be differentiated
			int i;
			Array<OneD,NekDouble> integrand(nTotQuadPoints);
			for(i = 0; i < nTotQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				double x = x1 * (1-xi[i])/2 + x2 * (1+xi[i])/2;

				integrand[i] = pow(x,7);
			}
			
			// Do the numerical differentiation by calling the corresponding Nektar++ routine of
			// the StdExpansion class

			Array<OneD, NekDouble> result(nQuadPoints);
			
			segExpansion->PhysDeriv(integrand,result);
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ":" << endl;
			for(i = 0; i < nQuadPoints; i++)
			{
				cout << "Quadrature Zero: " << xi[i] << "  Derivative: " << result[i]*invJacobian << endl;
			}
			
		}
		cout << endl;
	}

	cout << "-------------------" <<endl;
    cout << "-- EXERCISE 2(c) --" <<endl;
    cout << "-------------------" <<endl;

	cout << "METHOD 1: PolyLib" << endl;
	{
		// Specify the coordinates of the element
		double x1 = 0.0;
		double x2 = M_PI/2;

		// Input:Declare variable for quadrature points
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

			// Decalare variable that holds the values of the derivative matrix
			double* derivative = new double[nQuadPoints*nQuadPoints];

			// Calculate the derivative matrix using the quadrature zeros with
			// the corresponding Polylib routine
			Polylib::Dglj(derivative,quadZeros,nQuadPoints,alpha,beta);

			// Analytically evaluate the inverse of the Jacobian.
			double invJacobian, jacobian;
			invJacobian = 1.0/(0.5 * x2);
			jacobian = (0.5 * x2);

			// Calculate the derivative of each quadrature point
			// We use the result of this step to calculate the final integral
			int i;
			int j;
			double* deriv_result = new double[nQuadPoints];
							
			for(i = 0; i < nQuadPoints; i++)
			{
				deriv_result[i] = 0.0;
				for(j = 0; j < nQuadPoints; j++)
				{
					// Calculate the local coordinate of quadrature zeros
					double x = x2 * (1+quadZeros[j])/2;

					deriv_result[i] += derivative[i+nQuadPoints*j] * cos(x) * invJacobian;
				}
			}

			// Apply the Gaussian quadrature technique
			double result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				result += quadWeights[i] * deriv_result[i] * jacobian;
			}
			result = -result;

			// Display the output
			cout << "Result for Q = " << nQuadPoints << ": " << result << endl;

			// Deallocate the dynamic memory
			delete[] quadZeros;
			delete[] quadWeights;
			delete[] deriv_result;
			delete[] derivative;
		}
		cout << endl;
	
	}

	cout << "METHOD 2: low-level Nektar++" << endl;
	{
		// Specify the coordinates of the element
		NekDouble x1 = 0.0;
		NekDouble x2 = M_PI/2;

		// Input:Declare variable for quadrature points
		int nQuadPoints;

		// Declare the type of the quadrature points
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 2; nQuadPoints <= 8; nQuadPoints++)
		{
			//Decalre variables to hold the quadrature zeros and weights
			Array<OneD, NekDouble> quadZeros;
			Array<OneD, NekDouble> quadWeights;

			// Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
			// Step 1: Decalre a PointsKey which uniquele defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey1(nQuadPoints, quadPointsType);

			// Calculate the GLL-quadrature zeros and weights. This is done in 2 steps.
			// Step 1: Decalre a PointsKey which uniquele defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

			//Step 2:: Using this key, the quadrature zeros and weights can now be retrieved through
			// The PointsManager
			quadZeros = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			quadWeights = LibUtilities::PointsManager()[quadPointsKey]->GetW();

			// Calculate the derivative matrix
			LibUtilities::Points<NekDouble>::MatrixSharedPtrType derivMatrix = 
				LibUtilities::PointsManager()[quadPointsKey]->GetD();

			// Analytically evaluate the inverse of the Jacobian.
			NekDouble invJacobian, jacobian;
			invJacobian = 1.0/(0.5 * (-x1) + 0.5 * x2);
			jacobian = (0.5 * (-x1) + 0.5 * x2);

			// Calculate the derivative of each quadrature point
			// We use the result of this step to calculate the final integral
			Array<OneD, NekDouble> deriv_result(nQuadPoints);
			int i;
			int j;
			for(i = 0; i < nQuadPoints; i++)
			{
				deriv_result[i] = 0.0;
				for(j = 0; j < nQuadPoints; j++)
				{
					// Calculate the local coordinate of quadrature zeros
					NekDouble x = x1 * (1-quadZeros[j])/2 + x2 * (1+quadZeros[j])/2;
					
					deriv_result[i] += derivMatrix->GetValue(i,j) * cos(x) * invJacobian;
				}
			}

			// Apply the Guassian quadrature technique
			NekDouble result = 0.0;
			// Integrate
			for(i = 0; i < nQuadPoints; i++)
			{
				result += quadWeights[i] * deriv_result[i] * jacobian;
			}

			// Display the output
			cout << "Result for Q = " << nQuadPoints << ": " << -1 * result << endl;
		}
		cout << endl;

	}

	cout << "METHOD 3: Nektar++" << endl;
	{
		// Specify the coordinates of the element
		NekDouble x1 = 0.0;
		NekDouble x2 = M_PI/2;

		// Input:Declare variable for quadrature points
		int nQuadPoints;

		// Declare the type of the quadrature points
		LibUtilities::PointsType quadPointsType = LibUtilities::eGaussLobattoLegendre;

		for(nQuadPoints = 2; nQuadPoints <= 8; nQuadPoints++)
		{
			// Declare a PointsKey which uniquely defines the quadrature points
			const LibUtilities::PointsKey quadPointsKey(nQuadPoints, quadPointsType);

			// Declare a BasisKey which uniquely defines the (discrete) one-dimensional
			// spectral/hp expansion basis
			// As we will use the basis only to evaluate a derivative, the type and order of the basis
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
			
			// Analytically evaluate the inverse of the Jacobian.
			NekDouble invJacobian,jacobian;
			invJacobian = 1.0/(0.5 * (-x1) + 0.5 * x2);
			jacobian = (0.5 * (-x1) + 0.5 * x2);
			
			// Calculate the values of the function to be differentiated
			int i;
			Array<OneD,NekDouble> function(nTotQuadPoints);
			for(i = 0; i < nTotQuadPoints; i++)
			{
				// Calculate the local coordinate of quadrature zeros
				NekDouble x = x1 * (1-xi[i])/2 + x2 * (1+xi[i])/2;

				function[i] = invJacobian * cos(x);
			}

			// Do the numerical differentiation by calling the corresponding Nektar++ routine of
			// the StdExpansion class

			Array<OneD, NekDouble> deriv_result(nQuadPoints);
			
			segExpansion->PhysDeriv(function,deriv_result);

			// Do the numerical integration by calling the ecorresponding Nektar++ routine of
			// the StdExpansion class
			NekDouble result = segExpansion->Integral(deriv_result);
			result = result * jacobian;
			result = -1 * result;
			
			// Display the output
			cout << "Result for Q = " << nQuadPoints << ": " << result << endl;
		}
	}

        return 0;

}
