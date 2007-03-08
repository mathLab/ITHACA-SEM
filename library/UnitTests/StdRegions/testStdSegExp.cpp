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
                BOOST_CHECK_CLOSE(result[0], 2.0/3.0, epsilon);
                BOOST_CHECK_CLOSE(result[1], 1.0/3.0, epsilon);
                BOOST_CHECK_CLOSE(result[2], 1.0/6.0, epsilon);
                BOOST_CHECK_CLOSE(result[3], -1.0/15.0, epsilon);
                BOOST_CHECK_CLOSE(result[4], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[5], 0.0, epsilon);

                BOOST_CHECK_CLOSE(result[6], 1.0/3.0, epsilon);
                BOOST_CHECK_CLOSE(result[7], 2.0/3.0, epsilon);
                BOOST_CHECK_CLOSE(result[8], 1.0/6.0, epsilon);
                BOOST_CHECK_CLOSE(result[9], 1.0/15.0, epsilon);
                BOOST_CHECK_CLOSE(result[10], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[11], 0.0, epsilon);

                BOOST_CHECK_CLOSE(result[12], 1.0/6.0, epsilon);
                BOOST_CHECK_CLOSE(result[13], 1.0/6.0, epsilon);
                BOOST_CHECK_CLOSE(result[14], 1.0/15.0, epsilon);
                BOOST_CHECK_CLOSE(result[15], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[16], -1.0/70.0, epsilon);
                BOOST_CHECK_CLOSE(result[17], 0.0, epsilon);

                BOOST_CHECK_CLOSE(result[18], -1.0/15.0, epsilon);
                BOOST_CHECK_CLOSE(result[19], 1.0/15.0, epsilon);
                BOOST_CHECK_CLOSE(result[20], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[21], 4.0/105.0, epsilon);
                BOOST_CHECK_CLOSE(result[22], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[23], -4.0/315.0, epsilon);

                BOOST_CHECK_CLOSE(result[24], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[25], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[26], -1.0/70.0, epsilon);
                BOOST_CHECK_CLOSE(result[27], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[28], 2.0/70.0, epsilon);
                BOOST_CHECK_CLOSE(result[29], 0.0, epsilon);

                BOOST_CHECK_CLOSE(result[30], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[31], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[32], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[33], -4.0/315.0, epsilon);
                BOOST_CHECK_CLOSE(result[34], 0.0, epsilon);
                BOOST_CHECK_CLOSE(result[35], 16.0/693.0, epsilon);		
	    }
	}

	void testLapMatrix()
	{
	    // test the combination of Modified_A BasisType with GLL quad points
	    {
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
		
		double epsilon = 1e-12;	
		BOOST_CHECK_CLOSE(result[4], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[5], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[6], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[7], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[8], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[9], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[10], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[11], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[12], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[13], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[14], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[15], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[16], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[17], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[18], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[19], -4.0/3.0, epsilon);
		BOOST_CHECK_CLOSE(result[20], -4.0/3.0, epsilon);
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
		
		double epsilon = 1e-12;	
		BOOST_CHECK_CLOSE(result[6], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[7], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[8], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[9], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[10], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[11], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[12], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[13], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[14], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[15], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[16], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[17], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[18], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[19], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[20], -4.0, epsilon);
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
		BOOST_CHECK_CLOSE(result[0], -2.0/21.0, epsilon);
		BOOST_CHECK_CLOSE(result[1], -26.0/21.0, epsilon);
		BOOST_CHECK_CLOSE(result[2], -4.0/15.0, epsilon);
		BOOST_CHECK_CLOSE(result[3], -4.0/21.0, epsilon);
		BOOST_CHECK_CLOSE(result[4], 2.0/35.0, epsilon);
		BOOST_CHECK_CLOSE(result[5], -40.0/693.0, epsilon);
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
		BOOST_CHECK_CLOSE(result[0], 2.0, epsilon);
		BOOST_CHECK_CLOSE(result[1], -2.0, epsilon);
		BOOST_CHECK_CLOSE(result[2], -4.0, epsilon);
		BOOST_CHECK_CLOSE(result[3], 10.0/7.0, epsilon);
		BOOST_CHECK_CLOSE(result[4], 0.0, epsilon);
		BOOST_CHECK_CLOSE(result[5], -12.0/7.0, epsilon);
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
		BOOST_CHECK_CLOSE(result[0], f[0], epsilon);
		BOOST_CHECK_CLOSE(result[1], f[1], epsilon);
		BOOST_CHECK_CLOSE(result[2], f[2], epsilon);
		BOOST_CHECK_CLOSE(result[3], f[3], epsilon);
		BOOST_CHECK_CLOSE(result[4], f[4], epsilon);
		BOOST_CHECK_CLOSE(result[5], f[5], epsilon);
		BOOST_CHECK_CLOSE(result[6], f[6], epsilon);
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
		
		double epsilon = 1e-12;	
		BOOST_CHECK_CLOSE(result[6], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[7], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[8], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[9], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[10], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[11], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[12], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[13], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[14], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[15], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[16], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[17], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[18], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[19], -7.0/32.0, epsilon);
		BOOST_CHECK_CLOSE(result[20], -7.0/32.0, epsilon);
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
		BOOST_CHECK_CLOSE(resultL2[6], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[7], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[8], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[9], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[10], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[11], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[12], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[13], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[14], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[15], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[16], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[17], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[18], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[19], sqrt(1224.0/385.0), epsilon);
		BOOST_CHECK_CLOSE(resultL2[20], sqrt(1224.0/385.0), epsilon);

		BOOST_CHECK_SMALL(resultL2err[6], epsilon);
		BOOST_CHECK_SMALL(resultL2err[7], epsilon);
		BOOST_CHECK_SMALL(resultL2err[8], epsilon);
		BOOST_CHECK_SMALL(resultL2err[9], epsilon);
		BOOST_CHECK_SMALL(resultL2err[10], epsilon);
		BOOST_CHECK_SMALL(resultL2err[11], epsilon);
		BOOST_CHECK_SMALL(resultL2err[12], epsilon);
		BOOST_CHECK_SMALL(resultL2err[13], epsilon);
		BOOST_CHECK_SMALL(resultL2err[14], epsilon);
		BOOST_CHECK_SMALL(resultL2err[15], epsilon);
		BOOST_CHECK_SMALL(resultL2err[16], epsilon);
		BOOST_CHECK_SMALL(resultL2err[17], epsilon);
		BOOST_CHECK_SMALL(resultL2err[18], epsilon);
		BOOST_CHECK_SMALL(resultL2err[19], epsilon);
		BOOST_CHECK_SMALL(resultL2err[20], epsilon);

		BOOST_CHECK_SMALL(resultLinferr[6], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[7], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[8], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[9], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[10], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[11], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[12], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[13], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[14], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[15], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[16], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[17], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[18], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[19], epsilon);
		BOOST_CHECK_SMALL(resultLinferr[20], epsilon);
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

 **/
