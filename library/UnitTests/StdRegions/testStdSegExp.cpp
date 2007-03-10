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

namespace Nektar
{   // IS THIS THE CORRECT NAMESPACE (or should it be StdRegionsUnitTests)
    namespace StdSegExpUnitTests
    {
        using namespace Nektar;

        void testMassMatrix()
	{
	    // test the combination of Modified_A BasisType with GLL quad points
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nq = 7;
		int nummodes = 6;
		const LibUtilities::PointsKey Pkey(nq,Qtype);
		const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		
		DNekMatSharedPtr massmatrix = E->GenMassMatrix();
		double *result =  &((*massmatrix).GetPtr())[0];

		double epsilon = 1e-12;
                BOOST_CHECK_SMALL(2.0/3.0 - result[0], epsilon);
                BOOST_CHECK_SMALL(1.0/3.0 - result[1], epsilon);
                BOOST_CHECK_SMALL(1.0/6.0 - result[2], epsilon);
                BOOST_CHECK_SMALL(-1.0/15.0 - result[3], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[4], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[5], epsilon);

                BOOST_CHECK_SMALL(1.0/3.0 - result[6], epsilon);
                BOOST_CHECK_SMALL(2.0/3.0 - result[7], epsilon);
                BOOST_CHECK_SMALL(1.0/6.0 - result[8], epsilon);
                BOOST_CHECK_SMALL(1.0/15.0 - result[9], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[10], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[11], epsilon);

                BOOST_CHECK_SMALL(1.0/6.0 - result[12], epsilon);
                BOOST_CHECK_SMALL(1.0/6.0 - result[13], epsilon);
                BOOST_CHECK_SMALL(1.0/15.0 - result[14], epsilon);
                BOOST_CHECK_SMALL(0.0 -  result[15], epsilon);
                BOOST_CHECK_SMALL(-1.0/70.0 - result[16], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[17], epsilon);

                BOOST_CHECK_SMALL(-1.0/15.0 - result[18], epsilon);
                BOOST_CHECK_SMALL(1.0/15.0 - result[19], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[20], epsilon);
                BOOST_CHECK_SMALL(4.0/105.0 - result[21], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[22], epsilon);
                BOOST_CHECK_SMALL(-4.0/315.0 - result[23], epsilon);

                BOOST_CHECK_SMALL(0.0 - result[24], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[25], epsilon);
                BOOST_CHECK_SMALL(-1.0/70.0 - result[26], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[27], epsilon);
                BOOST_CHECK_SMALL(2.0/70.0 - result[28], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[29], epsilon);

                BOOST_CHECK_SMALL(0.0 - result[30], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[31], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[32], epsilon);
                BOOST_CHECK_SMALL(-4.0/315.0 - result[33], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[34], epsilon);
                BOOST_CHECK_SMALL(16.0/693.0 - result[35], epsilon);		
	    }
	}

	void testLapMatrix()
	{
	    // test the combination of Modified_A BasisType with GLL quad points
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nq = 7;
		int nummodes = 6;
		const LibUtilities::PointsKey Pkey(nq,Qtype);
		const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		
		DNekMatSharedPtr lapmatrix = E->GenLapMatrix();
		double *result =  &((*lapmatrix).GetPtr())[0];

		double epsilon = 1e-12;
                BOOST_CHECK_SMALL(1.0/2.0 - result[0], epsilon);
                BOOST_CHECK_SMALL(-1.0/2.0 - result[1], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[2], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[3], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[4], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[5], epsilon);

                BOOST_CHECK_SMALL(-1.0/2.0 - result[6], epsilon);
                BOOST_CHECK_SMALL(1.0/2.0 - result[7], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[8], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[9], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[10], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[11], epsilon);

                BOOST_CHECK_SMALL(0.0 - result[12], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[13], epsilon);
                BOOST_CHECK_SMALL(1.0/6.0 - result[14], epsilon);
                BOOST_CHECK_SMALL(0.0 -  result[15], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[16], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[17], epsilon);

                BOOST_CHECK_SMALL(0.0 - result[18], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[19], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[20], epsilon);
                BOOST_CHECK_SMALL(2.0/5.0 - result[21], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[22], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[23], epsilon);

                BOOST_CHECK_SMALL(0.0 - result[24], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[25], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[26], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[27], epsilon);
                BOOST_CHECK_SMALL(9.0/14.0 - result[28], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[29], epsilon);

                BOOST_CHECK_SMALL(0.0 - result[30], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[31], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[32], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[33], epsilon);
                BOOST_CHECK_SMALL(0.0 - result[34], epsilon);
                BOOST_CHECK_SMALL(8.0/9.0 - result[35], epsilon);
	    }
	}

	void testIntegration()
	{
	    // test the combination of Modified_A BasisType with GLL quad points (for a range of nq quad points)
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nummodes = 2;
		
		double result[21];
		for(int nq = 4; nq <= 20; nq++)
		{
		    const LibUtilities::PointsKey Pkey(nq,Qtype);
		    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		    BstShrDArray z = GetDoubleTmpSpace(nq);
		    BstShrDArray f = GetDoubleTmpSpace(nq);
		    double * tmp[1];
		    tmp[0] = z.get();
		    E->GetCoords(tmp);
		    for(int i = 0; i < nq; ++i)
		    {
			f[i] = 3*pow(z[i],5)-5*pow(z[i],3)+pow(z[i],2)-1;
		    }
		    result[nq] = E->Integral(&(f.get())[0]);
		}
		
		double expected_result = -4.0/3.0;
		double epsilon = 1e-12;	
		BOOST_CHECK_SMALL(expected_result - result[4], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[5], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[6], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[7], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[8], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[9], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[10], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[11], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[12], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[13], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[14], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[15], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[16], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[17], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[18], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[19], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[20], epsilon);
	    }
	}

	void testDifferentation()
	{
	    // test the combination of Modified_A BasisType with GLL quad points (for a range of nq quad points)
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nummodes = 2;
		
		double result[21];
		for(int nq = 6; nq <= 20; nq++)
		{
		    const LibUtilities::PointsKey Pkey(nq,Qtype);
		    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		    BstShrDArray z = GetDoubleTmpSpace(nq);
		    BstShrDArray f = GetDoubleTmpSpace(nq);
		    BstShrDArray df = GetDoubleTmpSpace(nq);
		    double * tmp[1];
		    tmp[0] = z.get();
		    E->GetCoords(tmp);
		    for(int i = 0; i < nq; ++i)
		    {
			f[i] = 3*pow(z[i],5)-5*pow(z[i],3)+pow(z[i],2)-1;
		    }
		    E->PhysDeriv(&(f.get())[0],&(df.get())[0]);
		    result[nq] = E->Integral(&(df.get())[0]);
		}
		
		double expected_result = -4.0;
		double epsilon = 1e-12;	
		BOOST_CHECK_SMALL(expected_result - result[6], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[7], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[8], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[9], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[10], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[11], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[12], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[13], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[14], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[15], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[16], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[17], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[18], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[19], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[20], epsilon);
	    }
	}

	void testIProductWRTBase()
	{
	    // test the combination of Modified_A BasisType with GLL quad points
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nummodes = 6;
		int nq = 7;
		const LibUtilities::PointsKey Pkey(nq,Qtype);
		const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		BstShrDArray z = GetDoubleTmpSpace(nq);
		BstShrDArray f = GetDoubleTmpSpace(nq);
		double * tmp[1];
		tmp[0] = z.get();
		E->GetCoords(tmp);
		for(int i = 0; i < nq; ++i)
		{
		    f[i] = 3*pow(z[i],5)-5*pow(z[i],3)+pow(z[i],2)-1;
		}

		double result[6];
		E->IProductWRTBase(&(f.get())[0],result);
		
		double epsilon = 1e-12;	
		BOOST_CHECK_SMALL(-2.0/21.0 - result[0], epsilon);
		BOOST_CHECK_SMALL(-26.0/21.0 - result[1], epsilon);
		BOOST_CHECK_SMALL(-4.0/15.0 - result[2], epsilon);
		BOOST_CHECK_SMALL(-4.0/21.0 - result[3], epsilon);
		BOOST_CHECK_SMALL(2.0/35.0 - result[4], epsilon);
		BOOST_CHECK_SMALL(-40.0/693.0 - result[5], epsilon);
	    }
	}

	void testFwdTrans()
	{
	    // test the combination of Modified_A BasisType with GLL quad points
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nummodes = 6;
		int nq = 7;
		const LibUtilities::PointsKey Pkey(nq,Qtype);
		const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		BstShrDArray z = GetDoubleTmpSpace(nq);
		BstShrDArray f = GetDoubleTmpSpace(nq);
		double * tmp[1];
		tmp[0] = z.get();
		E->GetCoords(tmp);
		for(int i = 0; i < nq; ++i)
		{
		    f[i] = 3*pow(z[i],5)-5*pow(z[i],3)+pow(z[i],2)-1;
		}

		E->FwdTrans(&(f.get())[0]);
		double *result = &(E->GetCoeffs())[0];
		
		double epsilon = 1e-12;	
		BOOST_CHECK_SMALL(2.0 - result[0], epsilon);
		BOOST_CHECK_SMALL(-2.0 - result[1], epsilon);
		BOOST_CHECK_SMALL(-4.0 - result[2], epsilon);
		BOOST_CHECK_SMALL(10.0/7.0 - result[3], epsilon);
		BOOST_CHECK_SMALL(0.0 - result[4],epsilon);
		BOOST_CHECK_SMALL(-12.0/7.0 - result[5], epsilon);
	    }
	}

	void testBwdTrans()
	{
	    // test the combination of Modified_A BasisType with GLL quad points
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nummodes = 6;
		int nq = 7;
		const LibUtilities::PointsKey Pkey(nq,Qtype);
		const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		BstShrDArray z = GetDoubleTmpSpace(nq);
		BstShrDArray f = GetDoubleTmpSpace(nq);
		double * tmp[1];
		tmp[0] = z.get();
		E->GetCoords(tmp);
		for(int i = 0; i < nq; ++i)
		{
		    f[i] = 3*pow(z[i],5)-5*pow(z[i],3)+pow(z[i],2)-1;
		}

		double result[7];
		E->FwdTrans(&(f.get())[0]);
		E->BwdTrans(result);

		double epsilon = 1e-12;	
		BOOST_CHECK_SMALL(f[0] - result[0], epsilon);
		BOOST_CHECK_SMALL(f[1] - result[1], epsilon);
		BOOST_CHECK_SMALL(f[2] - result[2], epsilon);
		BOOST_CHECK_SMALL(f[3] - result[3], epsilon);
		BOOST_CHECK_SMALL(f[4] - result[4], epsilon);
		BOOST_CHECK_SMALL(f[5] - result[5], epsilon);
		BOOST_CHECK_SMALL(f[6] - result[6], epsilon);
	    }
	}

	void testEvaluate()
	{
	    // test the combination of Modified_A BasisType with GLL quad points (for a range of nq quad points)
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		int nummodes = 2;
		
		double result[21];
		for(int nq = 6; nq <= 20; nq++)
		{
		    const LibUtilities::PointsKey Pkey(nq,Qtype);
		    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		    BstShrDArray z = GetDoubleTmpSpace(nq);
		    BstShrDArray f = GetDoubleTmpSpace(nq);
		    BstShrDArray df = GetDoubleTmpSpace(nq);
		    double * tmp[1];
		    tmp[0] = z.get();
		    E->GetCoords(tmp);
		    for(int i = 0; i < nq; ++i)
		    {
			f[i] = 3*pow(z[i],5)-5*pow(z[i],3)+pow(z[i],2)-1;
		    }
		    E->SetPhys(&(f.get())[0]);
		    double x = -0.5;
		    result[nq] = E->Evaluate(&x);
		}

		double expected_result = -7.0/32.0;		
		double epsilon = 1e-12;	
		BOOST_CHECK_SMALL(expected_result - result[6], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[7], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[8], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[9], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[10], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[11], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[12], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[13], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[14], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[15], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[16], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[17], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[18], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[19], epsilon);
		BOOST_CHECK_SMALL(expected_result - result[20], epsilon);
	    }
	}

	void testNorms()
	{
	    // test the combination of Modified_A BasisType with GLL quad points (for a range of nq quad points)
	    {
		LibUtilities::PointsType Qtype = LibUtilities::eGaussLobattoLegendre;
		LibUtilities::BasisType btype = LibUtilities::eModified_A;
		
		double resultL2[21];
		double resultL2err[21];
		double resultLinferr[21];
		for(int nummodes = 6; nummodes <= 20; nummodes++)
		{
		    int nq = ((nummodes+7 + (nummodes+7)%2 )/2);
		    const LibUtilities::PointsKey Pkey(nq,Qtype);
		    const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
		    StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); //or StdExpansion(1D)?
		    BstShrDArray z = GetDoubleTmpSpace(nq);
		    BstShrDArray f = GetDoubleTmpSpace(nq);
		    double * tmp[1];
		    tmp[0] = z.get();
		    E->GetCoords(tmp);
		    for(int i = 0; i < nq; ++i)
		    {
			f[i] = 3*pow(z[i],5)-5*pow(z[i],3)+pow(z[i],2)-1;
		    }
		    E->FwdTrans(&(f.get())[0]);
		    E->BwdTrans(&(E->GetPhys())[0]);
	       
		    resultL2[nummodes] = E->L2();
		    resultL2err[nummodes] = E->L2(&(f.get())[0]);
		    resultLinferr[nummodes] = E->Linf(&(f.get())[0]);
		}
		
		double epsilon = 1e-12;	
		double expected_result = sqrt(1224.0/385.0);
		BOOST_CHECK_SMALL(expected_result - resultL2[6], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[7], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[8], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[9], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[10], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[11], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[12], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[13], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[14], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[15], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[16], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[17], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[18], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[19], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2[20], epsilon);

		expected_result = 0.0;
		BOOST_CHECK_SMALL(expected_result - resultL2err[6], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[7], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[8], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[9], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[10], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[11], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[12], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[13], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[14], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[15], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[16], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[17], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[18], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[19], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultL2err[20], epsilon);

		expected_result = 0.0;
		BOOST_CHECK_SMALL(expected_result - resultLinferr[6], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[7], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[8], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[9], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[10], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[11], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[12], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[13], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[14], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[15], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[16], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[17], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[18], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[19], epsilon);
		BOOST_CHECK_SMALL(expected_result - resultLinferr[20], epsilon);
	    }
	}


    }
}
    /*
                BOOST_CHECK_CLOSE(result[0], , epsilon);
                BOOST_CHECK_CLOSE(result[1], , epsilon);
                BOOST_CHECK_CLOSE(result[2], , epsilon);
                BOOST_CHECK_CLOSE(result[3], , epsilon);
                BOOST_CHECK_CLOSE(result[4], , epsilon);
                BOOST_CHECK_CLOSE(result[5], , epsilon);

                BOOST_CHECK_CLOSE(result[6], , epsilon);
                BOOST_CHECK_CLOSE(result[7], , epsilon);
                BOOST_CHECK_CLOSE(result[8], , epsilon);
                BOOST_CHECK_CLOSE(result[9], , epsilon);
                BOOST_CHECK_CLOSE(result[10], , epsilon);
                BOOST_CHECK_CLOSE(result[11], , epsilon);

                BOOST_CHECK_CLOSE(result[12], , epsilon);
                BOOST_CHECK_CLOSE(result[13], , epsilon);
                BOOST_CHECK_CLOSE(result[14], , epsilon);
                BOOST_CHECK_CLOSE(result[15], , epsilon);
                BOOST_CHECK_CLOSE(result[16], , epsilon);
                BOOST_CHECK_CLOSE(result[17], , epsilon);

                BOOST_CHECK_CLOSE(result[18], , epsilon);
                BOOST_CHECK_CLOSE(result[19], , epsilon);
                BOOST_CHECK_CLOSE(result[20], , epsilon);
                BOOST_CHECK_CLOSE(result[21], , epsilon);
                BOOST_CHECK_CLOSE(result[22], , epsilon);
                BOOST_CHECK_CLOSE(result[23], , epsilon);

                BOOST_CHECK_CLOSE(result[24], , epsilon);
                BOOST_CHECK_CLOSE(result[25], , epsilon);
                BOOST_CHECK_CLOSE(result[26], , epsilon);
                BOOST_CHECK_CLOSE(result[27], , epsilon);
                BOOST_CHECK_CLOSE(result[28], , epsilon);
                BOOST_CHECK_CLOSE(result[29], , epsilon);

                BOOST_CHECK_CLOSE(result[30], , epsilon);
                BOOST_CHECK_CLOSE(result[31], , epsilon);
                BOOST_CHECK_CLOSE(result[32], , epsilon);
                BOOST_CHECK_CLOSE(result[33], , epsilon);
                BOOST_CHECK_CLOSE(result[34], , epsilon);
                BOOST_CHECK_CLOSE(result[35], , epsilon);
    */

/**
    $Log: testStdSegExp.cpp,v $
    Revision 1.1  2007/03/08 17:06:41  pvos
    added to repository


 **/
