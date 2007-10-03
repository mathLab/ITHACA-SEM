/////////////////////////////////////////////////////////////////////////////////
////
//// File: testStdSegExp.cpp
////
//// For more information, please see: http://www.nektar.info
////
//// The MIT License
////
//// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//// Department of Aeronautics, Imperial College London (UK), and Scientific
//// Computing and Imaging Institute, University of Utah (USA).
////
//// License for the specific language governing rights and limitations under
//// Permission is hereby granted, free of charge, to any person obtaining a
//// copy of this software and associated documentation files (the "Software"),
//// to deal in the Software without restriction, including without limitation
//// the rights to use, copy, modify, merge, publish, distribute, sublicense,
//// and/or sell copies of the Software, and to permit persons to whom the
//// Software is furnished to do so, subject to the following conditions:
////
//// The above copyright notice and this permission notice shall be included
//// in all copies or substantial portions of the Software.
////
//// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//// DEALINGS IN THE SOFTWARE.
////
//// Description:  Tests the standard segment expansion routines
////
/////////////////////////////////////////////////////////////////////////////////
//
//#include <UnitTests/StdRegions/testStdSegExp.h>
//#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
//#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
//
//#include <StdRegions/StdSegExp.h>
//#include <StdRegions/StdRegions.hpp>
//#include <LibUtilities/Foundations/Foundations.hpp>
//
//#include <boost/test/auto_unit_test.hpp>
//#include <boost/test/test_case_template.hpp>
//#include <boost/test/floating_point_comparison.hpp>
//#include <boost/test/unit_test.hpp>
//#include <boost/progress.hpp>
//
//#include <sstream>
//
//#define CHECK_CLOSE_ABS_EXP(lhs, rhs, tol, BType, Nmodes, PType, Npoints) \
//if(fabs(lhs-rhs)>tol) \
//{ \
//    std::ostringstream error_message; \
//    error_message << "difference between " \
//    << BOOST_TEST_STRINGIZE(lhs)<<"{"<<lhs<<"} and " \
//    << BOOST_TEST_STRINGIZE(rhs)<<"{"<<rhs<<"} exceeds " \
//    << tol << " for the case: " \
//    << "BasisType = " << LibUtilities::BasisTypeMap[BType] << "; " \
//    << "NumModes = " << Nmodes << "; " \
//    << "PointsType = " << LibUtilities::PointsTypeMap[PType] << "; " \
//    << "NumPoints = " << Npoints << "." ; \
//    BOOST_ERROR(error_message.str());\
//}
//
//#define CHECK_CLOSE_ARR_ABS_EXP(lhs, rhs, index, tol, BType, Nmodes, PType, Npoints) \
//if(fabs(lhs-rhs)>tol) \
//{ \
//    std::ostringstream error_message; \
//    error_message << "difference between " \
//    << BOOST_TEST_STRINGIZE(lhs)<<"{"<<lhs<<"} and " \
//    << BOOST_TEST_STRINGIZE(rhs)<<"{"<<rhs<<"} with " \
//    << BOOST_TEST_STRINGIZE(index) << " = " << index \
//    << ", exceeds " << tol <<" for the case: " \
//    << "BasisType = " << LibUtilities::BasisTypeMap[BType] << "; " \
//    << "NumModes = " << Nmodes << "; " \
//    << "PointsType = " << LibUtilities::PointsTypeMap[PType] << "; " \
//    << "NumPoints = " << Npoints << "." ; \
//    BOOST_ERROR(error_message.str());\
//}
//
//
//namespace Nektar
//{   // IS THIS THE CORRECT NAMESPACE (or should it be StdRegionsUnitTests)
//    namespace StdSegExpUnitTests
//    {
//	    const int max_nummodes = 15;
//	    const int max_nq = 15;
//
//        //const int num_BasisTypes = 5;
//        //LibUtilities::BasisType TestedBasisTypes[num_BasisTypes] = {LibUtilities::eOrtho_A,
//        //    LibUtilities::eModified_A,
//        //    LibUtilities::eGLL_Lagrange,
//        //    LibUtilities::eLegendre,
//        //    LibUtilities::eChebyshev};
//
//        //const int num_PointsTypes = 11;
//        //LibUtilities::PointsType TestedPointsTypes[num_PointsTypes] = {LibUtilities::eGaussGaussLegendre,
//        //    LibUtilities::eGaussRadauMLegendre,
//        //    LibUtilities::eGaussRadauPLegendre,
//        //    LibUtilities::eGaussLobattoLegendre,
//        //    LibUtilities::eGaussGaussChebyshev,
//        //    LibUtilities::eGaussRadauMChebyshev,
//        //    LibUtilities::eGaussRadauPChebyshev,
//        //    LibUtilities::eGaussLobattoChebyshev,
//        //    LibUtilities::eGaussRadauMAlpha0Beta1,
//        //    LibUtilities::eGaussRadauMAlpha0Beta2,
//        //    LibUtilities::ePolyEvenlySpaced};
//
//        const int num_BasisTypes = 3;
//        LibUtilities::BasisType TestedBasisTypes[num_BasisTypes] = {LibUtilities::eOrtho_A,
//            LibUtilities::eModified_A,
//            LibUtilities::eLegendre};
//
//        const int num_PointsTypes = 4;
//        LibUtilities::PointsType TestedPointsTypes[num_PointsTypes] = {LibUtilities::eGaussGaussLegendre,
//            LibUtilities::eGaussRadauMLegendre,
//            LibUtilities::eGaussRadauPLegendre,
//            LibUtilities::eGaussLobattoLegendre};
//
//            //LibUtilities::eGaussGaussChebyshev,
//            //LibUtilities::eGaussRadauMChebyshev,
//            //LibUtilities::eGaussRadauPChebyshev,
//            //LibUtilities::eGaussLobattoChebyshev,
//            //LibUtilities::eGaussRadauMAlpha0Beta1,
//            //LibUtilities::eGaussRadauMAlpha0Beta2,
//            //LibUtilities::ePolyEvenlySpaced};
//
//		
//
//        BOOST_AUTO_TEST_CASE(testMassMatrix)
//        {
//            NekDouble exactmatrices[3][36] = {
//	        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//	         0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
//	         0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
//	         0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
//	         0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
//	         0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
//	        {2.0/3.0, 1.0/3.0, 1.0/6.0, -1.0/15.0, 0.0, 0.0,
//	         1.0/3.0, 2.0/3.0, 1.0/6.0, 1.0/15.0, 0.0, 0.0,
//	         1.0/6.0, 1.0/6.0, 1.0/15.0, 0.0, -1.0/70.0, 0.0,
//	         -1.0/15.0, 1.0/15.0, 0.0, 4.0/105.0, 0.0, -4.0/315.0,
//	         0.0, 0.0, -1.0/70.0, 0.0, 2.0/70.0, 0.0,
//	         0.0, 0.0, 0.0, -4.0/315.0, 0.0, 16.0/693.0},
//	        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//	         0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
//	         0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
//	         0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
//	         0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
//	         0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
//            };
//
//            int nummodes = 6;
//            int nq_min[11]={(int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)(2.0*(nummodes-1.0)+1.0)};
//            
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//        	    
//	            for(int nq = nq_min[j] ; nq <= max_nq; nq++)
//	            {
//		        const LibUtilities::PointsKey Pkey(nq,Qtype);
//		        const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		        StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 
//
//		        DNekMatSharedPtr matrix = E->GenMatrix(StdRegions::eMass);
//		        NekDouble *result =  &((*matrix).GetPtr())[0];
//		        NekDouble *expected_result = exactmatrices[i];
//
//		        NekDouble epsilon = 1e-12;
//		        for(int k = 0; k < 36; k++)
//		        { 
//		            CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
//		        }
//	            }
//	        }	
//            }    
//        }
//
//        BOOST_AUTO_TEST_CASE(testLapMatrix)
//        {
//            NekDouble exactmatrices[3][36] = {
//	        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//	         0.0, 3.0, 0.0, sqrt(21.0), 0.0, sqrt(33.0),
//	         0.0, 0.0, 15.0, 0.0, sqrt(5.0)*9.0, 0.0,
//	         0.0, sqrt(21.0), 0.0, 42.0, 0.0, sqrt(77.0)*6.0,
//	         0.0, 0.0, sqrt(5.0)*9.0, 0.0, 90.0, 0.0,
//	         0.0, sqrt(33.0), 0.0, sqrt(77.0)*6.0, 0.0, 165.0},
//	        {1.0/2.0, -1.0/2.0, 0.0, 0.0, 0.0, 0.0,
//	         -1.0/2.0, 1.0/2.0, 0.0, 0.0, 0.0, 0.0,
//	         0.0, 0.0, 1.0/6.0, 0.0, 0.0, 0.0,
//	         0.0, 0.0, 0.0, 2.0/5.0, 0.0, 0.0,
//	         0.0, 0.0, 0.0, 0.0, 9.0/14.0, 0.0,
//	         0.0, 0.0, 0.0, 0.0, 0.0, 8.0/9.0},
//	        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//	         0.0, 3.0, 0.0, sqrt(21.0), 0.0, sqrt(33.0),
//	         0.0, 0.0, 15.0, 0.0, sqrt(5.0)*9.0, 0.0,
//	         0.0, sqrt(21.0), 0.0, 42.0, 0.0, sqrt(77.0)*6.0,
//	         0.0, 0.0, sqrt(5.0)*9.0, 0.0, 90.0, 0.0,
//	         0.0, sqrt(33.0), 0.0, sqrt(77.0)*6.0, 0.0, 165.0}
//	         };
//
//            int nummodes = 6;
//            int nq_min[11]={(int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)(2.0*(nummodes-1.0)+1.0)};
//            
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//        	    
//	            for(int nq = nq_min[j] ; nq <= max_nq; nq++)
//	            {
//		        const LibUtilities::PointsKey Pkey(nq,Qtype);
//		        const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		        StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 
//
//		        DNekMatSharedPtr matrix = E->GenMatrix(StdRegions::eLaplacian);
//		        NekDouble *result =  &((*matrix).GetPtr())[0];
//		        NekDouble *expected_result = exactmatrices[i];
//
//		        NekDouble epsilon = 1e-122;
//		        for(int k = 0; k < 36; k++)
//		        { 
//		            CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
//		        }
//	            }
//	        }	
//            }  	
//        }
//
//        BOOST_AUTO_TEST_CASE(testIntegration)
//        {
//            NekDouble expected_result = 2.0;//-4.0/3.0;
//
//            int order_f = 5;
//            int nq_min[11]={(int)ceil((order_f+1.0)/2.0),
//		          (int)ceil((order_f+2.0)/2.0),
//		          (int)ceil((order_f+2.0)/2.0),
//		          (int)ceil((order_f+3.0)/2.0),
//		          (int)ceil((order_f+1.0)/2.0),
//		          (int)ceil((order_f+2.0)/2.0),
//		          (int)ceil((order_f+2.0)/2.0),
//		          (int)ceil((order_f+3.0)/2.0),
//		          (int)ceil((order_f+2.0)/2.0),
//		          (int)ceil((order_f+2.0)/2.0),
//		          (int)(order_f+1.0)};
//            
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//
//	            for(int nummodes = 2; nummodes <= max_nummodes; nummodes++)
//	            {
//        		
//		        for(int nq = nq_min[j] ; nq <= max_nq; nq++)
//		        {
//		            const LibUtilities::PointsKey Pkey(nq,Qtype);
//		            const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		            StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 
//
//		            Array<OneD, NekDouble> z = Array<OneD, NekDouble>(nq);
//		            Array<OneD, NekDouble> f = Array<OneD, NekDouble>(nq);
//		            E->GetCoords(z);
//		            for(int r = 0; r < nq; ++r)
//		            {
//			        f[r] = 1.0;//3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
//		            }
//		            NekDouble result = E->Integral(f);
//		            NekDouble epsilon = 1e-12;
//		            CHECK_CLOSE_ABS_EXP(result, expected_result, epsilon, btype, nummodes, Qtype, nq);
//		        }
//	            }
//	        }	
//            }    
//        }
//
//        BOOST_AUTO_TEST_CASE(testDifferentiation)
//	    {
//    // 	    NekDouble expected_result = -4.0;
//    // 	    
//    // 	    for(int i = 0; i < num_BasisTypes; i++)
//    // 	    {
//    // 		    LibUtilities::BasisType btype = TestedBasisTypes[i];
//    // 
//    // 		    for(int j = 0; j < num_PointsTypes; j++)
//    // 		    {
//    // 		        LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//    // 
//    // 			    Array<OneD, NekDouble> z = Array<OneD, NekDouble>(nq);
//    // 			    Array<OneD, NekDouble> f = Array<OneD, NekDouble>(nq);
//    // 			    Array<OneD, NekDouble> df = Array<OneD, NekDouble>(nq);
//    // 			    E->GetCoords(z);
//    // 			    for(int r = 0; r < nq; r++)
//    // 			    {
//    // 				f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
//    // 			    }
//    // 			    E->PhysDeriv(f,df);
//    // 			    NekDouble result = E->Integral(df);
//    // 			    NekDouble epsilon = 1e-12;
//    // 			    CHECK_CLOSE_ABS_EXP(result, expected_result, epsilon, btype, nummodes, Qtype, nq);
//    // 			}
//    // 		    }
//    // 		}	
//    // 	    }    
//	    }
//
//        BOOST_AUTO_TEST_CASE(testIProductWRTBase)
//        {
//            NekDouble exactmatrices[3][6] = {
//	        {-sqrt(0.5)*4.0/3.0, -sqrt(1.5)*8.0/7.0, sqrt(2.5)*4.0/15.0, -sqrt(3.5)*4.0/21.0, 0.0, sqrt(5.5)*16.0/231.0},
//	        {-2.0/21.0, -26.0/21.0, -4.0/15.0, -4.0/21.0, 2.0/35.0, -40.0/693.0},
//	        {-sqrt(0.5)*4.0/3.0, -sqrt(1.5)*8.0/7.0, sqrt(2.5)*4.0/15.0, -sqrt(3.5)*4.0/21.0, 0.0, sqrt(5.5)*16.0/231.0}
//            };
//
//            int nummodes = 6;
//            int nq_min[11]={(int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)(2.0*(nummodes-1.0)+1.0)};
//            
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//        	    
//	            for(int nq = nq_min[j] ; nq <= max_nq; nq++)
//	            {
//		        const LibUtilities::PointsKey Pkey(nq,Qtype);
//		        const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		        StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey);
//
//		        Array<OneD, NekDouble> z = Array<OneD, NekDouble>(nq);
//		        Array<OneD, NekDouble> f = Array<OneD, NekDouble>(nq);
//		        E->GetCoords(z);
//		        for(int r = 0; r < nq; r++)
//		        {
//		            f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
//		        }		
//		        Array<OneD, NekDouble> result = Array<OneD, NekDouble>(6);
//		        E->IProductWRTBase(f,result);
//		        NekDouble *expected_result = exactmatrices[i];
//        		
//		        NekDouble epsilon = 1e-12;	
//		        for(int k = 0; k < 6; k++)
//		        { 
//		            CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
//		        }
//	            }
//	        }	
//            }    
//        }
//
//        BOOST_AUTO_TEST_CASE(testFwdTrans)
//        {
//            NekDouble exactmatrices[3][6] = {
//	        {-sqrt(2.0)*2.0/3.0, -sqrt(2.0/3.0)*12.0/7.0, sqrt(0.4)*2.0/3.0, -sqrt(2.0/7.0)*2.0/3.0, 0.0, sqrt(2.0/11.0)*(8.0/21.0)},
//	        {2.0, -2.0, -4.0, 10.0/7.0, 0.0, -12.0/7.0},
//	        {-sqrt(2.0)*2.0/3.0, -sqrt(2.0/3.0)*12.0/7.0, sqrt(0.4)*2.0/3.0, -sqrt(2.0/7.0)*2.0/3.0, 0.0, sqrt(2.0/11.0)*(8.0/21.0)}
//            };
//
//            int nummodes = 6;
//            int nq_min[11]={(int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//		            (int)(2.0*(nummodes-1.0)+1.0)};
//            
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//        	    
//	            for(int nq = nq_min[j] ; nq <= max_nq; nq++)
//	            {
//		        const LibUtilities::PointsKey Pkey(nq,Qtype);
//		        const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		        StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey);
//
//		        Array<OneD, NekDouble> z = Array<OneD, NekDouble>(nq);
//		        Array<OneD, NekDouble> f = Array<OneD, NekDouble>(nq);
//		        E->GetCoords(z);
//		        for(int r = 0; r < nq; r++)
//		        {
//		            f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
//		        }	
//		        Array<OneD, NekDouble> result = Array<OneD, NekDouble>(nummodes);	
//		        E->FwdTrans(f,result);
//		        NekDouble *expected_result = exactmatrices[i];
//
//		        NekDouble epsilon = 1e-12;	
//		        for(int k = 0; k < 6; k++)
//		        { 
//		            CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
//		        }
//	            }
//	        }	
//            }    
//        }
//
//        BOOST_AUTO_TEST_CASE(testBwdTrans)
//        {	    
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//
//	            for(int nummodes = 6; nummodes <= max_nummodes; nummodes++)
//	            {
//		        int nq_min[11]={(int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)(2.0*(nummodes-1.0)+1.0)};
//        		
//		        for(int nq = nq_min[j] ; nq <= max_nq; nq++)
//		        {
//		            const LibUtilities::PointsKey Pkey(nq,Qtype);
//		            const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		            StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 
//
//		            Array<OneD, NekDouble> z = Array<OneD, NekDouble>(nq);
//		            Array<OneD, NekDouble> f = Array<OneD, NekDouble>(nq);
//		            E->GetCoords(z);
//		            for(int r = 0; r < nq; ++r)
//		            {
//			        f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
//		            }
//		            Array<OneD, NekDouble> tmp = Array<OneD, NekDouble>(nummodes);
//		            Array<OneD, NekDouble> result = Array<OneD, NekDouble>(nq);
//		            E->FwdTrans(f,tmp);
//		            E->BwdTrans(tmp,result);
//		            NekDouble *expected_result = f.get();
//
//		            NekDouble epsilon = 1e-12;	
//		            for(int k = 0; k < nq; k++)
//		            { 
//			        CHECK_CLOSE_ARR_ABS_EXP(result[k], expected_result[k], k, epsilon, btype, nummodes, Qtype, nq);
//		            }
//		        }
//	            }
//	        }	
//            }    
//        }
//
//        BOOST_AUTO_TEST_CASE(testPhysEvaluate)
//        {
//            NekDouble expected_result = -7.0/32.0;
//            
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//
//	            for(int nummodes = 2; nummodes <= max_nummodes; nummodes++)
//	            {
//        		
//		        for(int nq = 6 ; nq <= max_nq; nq++)
//		        {
//		            const LibUtilities::PointsKey Pkey(nq,Qtype);
//		            const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		            StdRegions::StdSegExp *E = new StdRegions::StdSegExp(Bkey); 
//
//		            Array<OneD, NekDouble> z = Array<OneD, NekDouble>(nq);
//		            Array<OneD, NekDouble> f = Array<OneD, NekDouble>(nq);
//		            E->GetCoords(z);
//		            for(int r = 0; r < nq; ++r)
//		            {
//			        f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
//		            }
//		            E->SetPhys(f);
//		            Array<OneD, NekDouble> x = Array<OneD, NekDouble>(1);
//		            x[0] = -0.5;
//		            NekDouble result = E->PhysEvaluate(x);
//
//		            NekDouble epsilon = 1e-12;
//		            CHECK_CLOSE_ABS_EXP(result, expected_result, epsilon, btype, nummodes, Qtype, nq);
//		        }
//	            }
//	        }	
//            }    
//        }
//
//        BOOST_AUTO_TEST_CASE(testNorms)
//        {
//            NekDouble expected_resultL2 = sqrt(1224.0/385.0);
//            NekDouble expected_resultL2err = 0.0;
//            NekDouble expected_resultLinferr = 0.0;
//            
//            for(int i = 0; i < num_BasisTypes; i++)
//            {
//	        LibUtilities::BasisType btype = TestedBasisTypes[i];
//
//	        for(int j = 0; j < num_PointsTypes; j++)
//	        {
//	            LibUtilities::PointsType Qtype = TestedPointsTypes[j];
//
//	            for(int nummodes = 6; nummodes <= max_nummodes; nummodes++)
//	            {
//		        int nq_min[11]={(int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+1.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+3.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)ceil((2*(nummodes-1.0)+2.0)/2.0),
//				        (int)(2.0*(nummodes-1.0)+1.0)};
//        		
//		        for(int nq = nq_min[j]; nq <= max_nq; nq++)
//		        {
//		            const LibUtilities::PointsKey Pkey(nq,Qtype);
//		            const LibUtilities::BasisKey Bkey(btype,nummodes,Pkey);
//		            StdRegions::StdExpansion1D *E = new StdRegions::StdSegExp(Bkey); 
//
//		            Array<OneD, NekDouble> z = Array<OneD, NekDouble>(nq);
//		            Array<OneD, NekDouble> f = Array<OneD, NekDouble>(nq);
//		            E->GetCoords(z);
//		            for(int r = 0; r < nq; ++r)
//		            {
//			        f[r] = 3*pow(z[r],5)-5*pow(z[r],3)+pow(z[r],2)-1;
//		            }
//		            E->SetPhys(f);
//		            E->FwdTrans(*E);
//		            E->BwdTrans(*E);
//               
//		            NekDouble resultL2 = E->L2();
//		            NekDouble resultL2err = E->L2(f);
//		            NekDouble resultLinferr = E->Linf(f);
//
//		            NekDouble epsilon = 1e-12;
//		            CHECK_CLOSE_ABS_EXP(resultL2, expected_resultL2, epsilon, btype, nummodes, Qtype, nq);
//		            CHECK_CLOSE_ABS_EXP(resultL2err, expected_resultL2err, epsilon, btype, nummodes, Qtype, nq);
//		            CHECK_CLOSE_ABS_EXP(resultLinferr, expected_resultLinferr, epsilon, btype, nummodes, Qtype, nq);
//		        }
//	            }
//	        }	
//            }    
//        }
//    }
//}
//
///**
//    $Log: testStdSegExp.cpp,v $
//    Revision 1.14  2007/09/12 03:59:42  bnelson
//    *** empty log message ***
//
//    Revision 1.13  2007/07/16 08:15:09  sherwin
//    Updated for new matrix calls
//
//    Revision 1.12  2007/05/27 22:01:26  bnelson
//    *** empty log message ***
//
//    Revision 1.11  2007/05/27 16:40:04  bnelson
//    *** empty log message ***
//
//    Revision 1.10  2007/05/15 05:21:04  bnelson
//    Updated to use the new Array object.
//
//    Revision 1.9  2007/03/31 00:40:16  bnelson
//    *** empty log message ***
//
//    Revision 1.8  2007/03/29 19:45:46  bnelson
//    *** empty log message ***
//
//    Revision 1.7  2007/03/26 11:17:03  pvos
//    made testStdRegions back working
//
//    Revision 1.6  2007/03/20 03:50:28  bnelson
//    Removed a lot of warning messages.
//
//    Revision 1.5  2007/03/16 14:47:56  pvos
//    *** empty log message ***
//
//    Revision 1.4  2007/03/16 12:10:37  pvos
//    adapted to the new code
//
//    Revision 1.3  2007/03/14 12:01:00  pvos
//    Added testcases
//
//    Revision 1.2  2007/03/10 13:04:17  pvos
//    *** empty log message ***
//
//    Revision 1.1  2007/03/08 17:06:41  pvos
//    added to repository
//
//
// **/
