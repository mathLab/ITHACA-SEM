///////////////////////////////////////////////////////////////////////////////
//
// File: testStdSegExp.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description:  Tests the standard segment expansion routines
//
///////////////////////////////////////////////////////////////////////////////

#include <UnitTests/StdRegions/testStdSegExp.h>

#include <StdRegions/StdSegExp.h>
#include "StdRegions/StdRegions.hpp"
#include "LibUtilities/Foundations/Foundations.hpp"

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/progress.hpp>

#include <sstream>

#define CHECK_CLOSE_ABS_EXP(lhs, rhs, tol, BType, Nmodes, PType, Npoints) \
if(fabs(lhs-rhs)>tol) \
{ \
    std::ostringstream error_message; \
    error_message << "difference between " \
    << BOOST_TEST_STRINGIZE(lhs)<<"{"<<lhs<<"} and " \
    << BOOST_TEST_STRINGIZE(rhs)<<"{"<<rhs<<"} exceeds " \
    << tol << " for the case: " \
    << "BasisType = " << LibUtilities::BasisTypeMap[BType] << "; " \
    << "NumModes = " << Nmodes << "; " \
    << "PointsType = " << LibUtilities::PointsTypeMap[PType] << "; " \
    << "NumPoints = " << Npoints << "." ; \
    BOOST_ERROR(error_message.str());\
}

#define CHECK_CLOSE_ARR_ABS_EXP(lhs, rhs, index, tol, BType, Nmodes, PType, Npoints) \
if(fabs(lhs-rhs)>tol) \
{ \
    std::ostringstream error_message; \
    error_message << "difference between " \
    << BOOST_TEST_STRINGIZE(lhs)<<"{"<<lhs<<"} and " \
    << BOOST_TEST_STRINGIZE(rhs)<<"{"<<rhs<<"} with " \
    << BOOST_TEST_STRINGIZE(index) << " = " << index \
    << ", exceeds " << tol <<" for the case: " \
    << "BasisType = " << LibUtilities::BasisTypeMap[BType] << "; " \
    << "NumModes = " << Nmodes << "; " \
    << "PointsType = " << LibUtilities::PointsTypeMap[PType] << "; " \
    << "NumPoints = " << Npoints << "." ; \
    BOOST_ERROR(error_message.str());\
}


namespace Nektar
{   // IS THIS THE CORRECT NAMESPACE (or should it be StdRegionsUnitTests)
    namespace StdSegExpUnitTests
    {
        using namespace Nektar;

	const int max_nummodes = 15;
	const int max_nq = 15;
	/*	
	const int num_BasisTypes = 5;
	LibUtilities::BasisType TestedBasisTypes[num_BasisTypes] = {LibUtilities::eOrtho_A,
								    LibUtilities::eModified_A,
								    LibUtilities::eGLL_Lagrange,
								    LibUtilities::eLegendre,
								    LibUtilities::eChebyshev};

	const int num_PointsTypes = 11;
	LibUtilities::PointsType TestedPointsTypes[num_PointsTypes] = {LibUtilities::eGaussGaussLegendre,
								       LibUtilities::eGaussRadauMLegendre,
								       LibUtilities::eGaussRadauPLegendre,
								       LibUtilities::eGaussLobattoLegendre,
								       LibUtilities::eGaussGaussChebyshev,
								       LibUtilities::eGaussRadauMChebyshev,
								       LibUtilities::eGaussRadauPChebyshev,
								       LibUtilities::eGaussLobattoChebyshev,
								       LibUtilities::eGaussRadauMAlpha0Beta1,
								       LibUtilities::eGaussRadauMAlpha0Beta2,
								       LibUtilities::ePolyEvenlySpaced};
	*/

	const int num_BasisTypes = 3;
	LibUtilities::BasisType TestedBasisTypes[num_BasisTypes] = {LibUtilities::eOrtho_A,
								    LibUtilities::eModified_A,
								    LibUtilities::eLegendre};

	const int num_PointsTypes = 4;
	LibUtilities::PointsType TestedPointsTypes[num_PointsTypes] = {LibUtilities::eGaussGaussLegendre,
								       LibUtilities::eGaussRadauMLegendre,
								       LibUtilities::eGaussRadauPLegendre,
								       LibUtilities::eGaussLobattoLegendre};
	/*
								       LibUtilities::eGaussGaussChebyshev,
								       LibUtilities::eGaussRadauMChebyshev,
								       LibUtilities::eGaussRadauPChebyshev,
								       LibUtilities::eGaussLobattoChebyshev,
								       LibUtilities::eGaussRadauMAlpha0Beta1,
								       LibUtilities::eGaussRadauMAlpha0Beta2,
								       LibUtilities::ePolyEvenlySpaced};
	*/	

        void testMassMatrix()
	{
	    double exactmatrices[3][36] = {
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,},
		{2.0/3.0, 1.0/3.0, 1.0/6.0, -1.0/15.0, 0.0, 0.0,
		 1.0/3.0, 2.0/3.0, 1.0/6.0, 1.0/15.0, 0.0, 0.0,
		 1.0/6.0, 1.0/6.0, 1.0/15.0, 0.0, -1.0/70.0, 0.0,
		 -1.0/15.0, 1.0/15.0, 0.0, 4.0/105.0, 0.0, -4.0/315.0,
		 0.0, 0.0, -1.0/70.0, 0.0, 2.0/70.0, 0.0,
		 0.0, 0.0, 0.0, -4.0/315.0, 0.0, 16.0/693.0},
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,}
	    };

	    int nummodes = 6;
	    int nq_min[11]={ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    2.0*(nummodes-1.0)+1.0};
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];
		    
		    for(int nq = nq_min[j] ; nq <= max_nq; nq++)
		    {
			const LibUtilities::PointsKey Pkey(nq,Qtype);
			const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 

			DNekMatSharedPtr matrix = E->GenMassMatrix();
			double *result =  &((*matrix).GetPtr())[0];
			double *expected_result = exactmatrices[i];

			double epsilon = 1e-12;
			for(int k = 0; k < 36; k++)
			{ 
			    CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}

        void testLapMatrix()
	{
	    double exactmatrices[3][36] = {
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 3.0, 0.0, sqrt(21.0), 0.0, sqrt(33.0),
		 0.0, 0.0, 15.0, 0.0, sqrt(5.0)*9.0, 0.0,
		 0.0, sqrt(21.0), 0.0, 42.0, 0.0, sqrt(77.0)*6.0,
		 0.0, 0.0, sqrt(5.0)*9.0, 0.0, 90.0, 0.0,
		 0.0, sqrt(33.0), 0.0, sqrt(77.0)*6.0, 0.0, 165.0,},
		{1.0/2.0, -1.0/2.0, 0.0, 0.0, 0.0, 0.0,
		 -1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 1.0/6.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 2.0/5.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 9.0/14.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 0.0, 8.0/9.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 3.0, 0.0, sqrt(21.0), 0.0, sqrt(33.0),
		 0.0, 0.0, 15.0, 0.0, sqrt(5.0)*9.0, 0.0,
		 0.0, sqrt(21.0), 0.0, 42.0, 0.0, sqrt(77.0)*6.0,
		 0.0, 0.0, sqrt(5.0)*9.0, 0.0, 90.0, 0.0,
		 0.0, sqrt(33.0), 0.0, sqrt(77.0)*6.0, 0.0, 165.0,}
	    };

	    int nummodes = 6;
	    int nq_min[11]={ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    2.0*(nummodes-1.0)+1.0};
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];
		    
		    for(int nq = nq_min[j] ; nq <= max_nq; nq++)
		    {
			const LibUtilities::PointsKey Pkey(nq,Qtype);
			const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 

			DNekMatSharedPtr matrix = E->GenLapMatrix();
			double *result =  &((*matrix).GetPtr())[0];
			double *expected_result = exactmatrices[i];

			double epsilon = 1e-11;
			for(int k = 0; k < 36; k++)
			{ 
			    CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}

        void testIntegration()
	{
	    double expected_result = -4.0/3.0;

	    int order_f = 5;
	    int nq_min[11]={ceil((order_f+1.0)/2.0),
			  ceil((order_f+2.0)/2.0),
			  ceil((order_f+2.0)/2.0),
			  ceil((order_f+3.0)/2.0),
			  ceil((order_f+1.0)/2.0),
			  ceil((order_f+2.0)/2.0),
			  ceil((order_f+2.0)/2.0),
			  ceil((order_f+3.0)/2.0),
			  ceil((order_f+2.0)/2.0),
			  ceil((order_f+2.0)/2.0),
			  order_f+1.0};
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];

		    for(int nummodes = 2; nummodes <= max_nummodes; nummodes++)
		    {
			
			for(int nq = nq_min[j] ; nq <= max_nq; nq++)
			{
			    const LibUtilities::PointsKey Pkey(nq,Qtype);
			    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 

			    BstShrDArray z = GetDoubleTmpSpace(nq);
			    BstShrDArray f = GetDoubleTmpSpace(nq);
			    double * tmp[1];
			    tmp[0] = z.get();
			    E->GetCoords(tmp);
			    for(int r = 0; r < nq; ++r)
			    {
				f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
			    }
			    double result = E->Integral(&(f.get())[0]);
			    double epsilon = 1e-12;
			    CHECK_CLOSE_ABS_EXP(result, expected_result, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}

        void testDifferentiation()
	{
	    double expected_result = -4.0;
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];

		    for(int nummodes = 2; nummodes <= max_nummodes; nummodes++)
		    {
			
			for(int nq = 6 ; nq <= max_nq; nq++)
			{

			    const LibUtilities::PointsKey Pkey(nq,Qtype);
			    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey);

			    BstShrDArray z = GetDoubleTmpSpace(nq);
			    BstShrDArray f = GetDoubleTmpSpace(nq);
			    BstShrDArray df = GetDoubleTmpSpace(nq);
			    double * tmp[1];
			    tmp[0] = z.get();
			    E->GetCoords(tmp);
			    for(int r = 0; r < nq; r++)
			    {
				f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
			    }
			    E->PhysDeriv(&(f.get())[0],&(df.get())[0]);
			    double result = E->Integral(&(df.get())[0]);
			    double epsilon = 1e-12;
			    CHECK_CLOSE_ABS_EXP(result, expected_result, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}

        void testIProductWRTBase()
	{
	    double exactmatrices[3][6] = {
		{-sqrt(0.5)*4.0/3.0, -sqrt(1.5)*8.0/7.0, sqrt(2.5)*4.0/15.0, -sqrt(3.5)*4.0/21.0, 0.0, sqrt(5.5)*16.0/231.0},
		{-2.0/21.0, -26.0/21.0, -4.0/15.0, -4.0/21.0, 2.0/35.0, -40.0/693.0},
		{-sqrt(0.5)*4.0/3.0, -sqrt(1.5)*8.0/7.0, sqrt(2.5)*4.0/15.0, -sqrt(3.5)*4.0/21.0, 0.0, sqrt(5.5)*16.0/231.0}
	    };

	    int nummodes = 6;
	    int nq_min[11]={ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    2.0*(nummodes-1.0)+1.0};
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];
		    
		    for(int nq = nq_min[j] ; nq <= max_nq; nq++)
		    {
			const LibUtilities::PointsKey Pkey(nq,Qtype);
			const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey);

			BstShrDArray z = GetDoubleTmpSpace(nq);
			BstShrDArray f = GetDoubleTmpSpace(nq);
			double * tmp[1];
			tmp[0] = z.get();
			E->GetCoords(tmp);
			for(int r = 0; r < nq; r++)
			{
			    f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
			}		
			double result[6];
			E->IProductWRTBase(&(f.get())[0],result);
			double *expected_result = exactmatrices[i];
			
			double epsilon = 1e-12;	
			for(int k = 0; k < 6; k++)
			{ 
			    CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}

        void testFwdTrans()
	{
	    double exactmatrices[3][6] = {
		{-sqrt(2)*2.0/3.0, -sqrt(2.0/3.0)*12.0/7.0, sqrt(0.4)*2.0/3.0, -sqrt(2.0/7.0)*2.0/3.0, 0.0, sqrt(2.0/11.0)*(8.0/21.0)},
		{2.0, -2.0, -4.0, 10.0/7.0, 0.0, -12.0/7.0},
		{-sqrt(2)*2.0/3.0, -sqrt(2.0/3.0)*12.0/7.0, sqrt(0.4)*2.0/3.0, -sqrt(2.0/7.0)*2.0/3.0, 0.0, sqrt(2.0/11.0)*(8.0/21.0)}
	    };

	    int nummodes = 6;
	    int nq_min[11]={ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+1.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+3.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    ceil((2*(nummodes-1.0)+2.0)/2.0),
			    2.0*(nummodes-1.0)+1.0};
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];
		    
		    for(int nq = nq_min[j] ; nq <= max_nq; nq++)
		    {
			const LibUtilities::PointsKey Pkey(nq,Qtype);
			const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey);

			BstShrDArray z = GetDoubleTmpSpace(nq);
			BstShrDArray f = GetDoubleTmpSpace(nq);
			double * tmp[1];
			tmp[0] = z.get();
			E->GetCoords(tmp);
			for(int r = 0; r < nq; r++)
			{
			    f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
			}		
			E->FwdTrans(&(f.get())[0]);
			double *result = &(E->GetCoeffs())[0];
			double *expected_result = exactmatrices[i];

			double epsilon = 1e-12;	
			for(int k = 0; k < 6; k++)
			{ 
			    CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}

        void testBwdTrans()
	{	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];

		    for(int nummodes = 6; nummodes <= max_nummodes; nummodes++)
		    {
			int nq_min[11]={ceil((2*(nummodes-1.0)+1.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+3.0)/2.0),
					ceil((2*(nummodes-1.0)+1.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+3.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					2.0*(nummodes-1.0)+1.0};
			
			for(int nq = nq_min[j] ; nq <= max_nq; nq++)
			{
			    const LibUtilities::PointsKey Pkey(nq,Qtype);
			    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 

			    BstShrDArray z = GetDoubleTmpSpace(nq);
			    BstShrDArray f = GetDoubleTmpSpace(nq);
			    double * tmp[1];
			    tmp[0] = z.get();
			    E->GetCoords(tmp);
			    for(int r = 0; r < nq; ++r)
			    {
				f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
			    }
			    E->FwdTrans(&(f.get())[0]);
			    double *result = new double [nq];
			    E->BwdTrans(result);
			    double *expected_result = f.get();

			    double epsilon = 1e-12;	
			    for(int k = 0; k < nq; k++)
			    { 
				CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
			    }
			}
		    }
		}	
	    }    
	}

        void testEvaluate()
	{
	    double expected_result = -7.0/32.0;
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];

		    for(int nummodes = 2; nummodes <= max_nummodes; nummodes++)
		    {
			
			for(int nq = 6 ; nq <= max_nq; nq++)
			{
			    const LibUtilities::PointsKey Pkey(nq,Qtype);
			    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 

			    BstShrDArray z = GetDoubleTmpSpace(nq);
			    BstShrDArray f = GetDoubleTmpSpace(nq);
			    double * tmp[1];
			    tmp[0] = z.get();
			    E->GetCoords(tmp);
			    for(int r = 0; r < nq; ++r)
			    {
				f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
			    }
			    E->SetPhys(&(f.get())[0]);
			    double x = -0.5;
			    double result = E->Evaluate(&x);

			    double epsilon = 1e-12;
			    CHECK_CLOSE_ABS_EXP(result, expected_result, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}

        void testNorms()
	{
	    double expected_resultL2 = sqrt(1224.0/385.0);
	    double expected_resultL2err = 0.0;
	    double expected_resultLinferr = 0.0;
	    
	    for(int i = 0; i < num_BasisTypes; i++)
	    {
		LibUtilities::BasisType btype = TestedBasisTypes[i];

		for(int j = 0; j < num_PointsTypes; j++)
		{
		    LibUtilities::PointsType Qtype = TestedPointsTypes[j];

		    for(int nummodes = 6; nummodes <= max_nummodes; nummodes++)
		    {
			int nq_min[11]={ceil((2*(nummodes-1.0)+1.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+3.0)/2.0),
					ceil((2*(nummodes-1.0)+1.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+3.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					ceil((2*(nummodes-1.0)+2.0)/2.0),
					2.0*(nummodes-1.0)+1.0};
			
			for(int nq = nq_min[j]; nq <= max_nq; nq++)
			{
			    const LibUtilities::PointsKey Pkey(nq,Qtype);
			    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
			    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 

			    BstShrDArray z = GetDoubleTmpSpace(nq);
			    BstShrDArray f = GetDoubleTmpSpace(nq);
			    double * tmp[1];
			    tmp[0] = z.get();
			    E->GetCoords(tmp);
			    for(int r = 0; r < nq; ++r)
			    {
				f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
			    }
			    E->FwdTrans(&(f.get())[0]);
			    E->BwdTrans(&(E->GetPhys())[0]);
	       
			    double resultL2 = E->L2();
			    double resultL2err = E->L2(&(f.get())[0]);
			    double resultLinferr = E->Linf(&(f.get())[0]);

			    double epsilon = 1e-12;
			    CHECK_CLOSE_ABS_EXP(resultL2, expected_resultL2, epsilon, btype, nummodes, Qtype, nq);
			    CHECK_CLOSE_ABS_EXP(resultL2err, expected_resultL2err, epsilon, btype, nummodes, Qtype, nq);
			    CHECK_CLOSE_ABS_EXP(resultLinferr, expected_resultLinferr, epsilon, btype, nummodes, Qtype, nq);
			}
		    }
		}	
	    }    
	}
    }
}

/**
    $Log: testStdSegExp.cpp,v $
    Revision 1.2  2007/03/10 13:04:17  pvos
    *** empty log message ***

    Revision 1.1  2007/03/08 17:06:41  pvos
    added to repository


 **/
