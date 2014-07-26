///////////////////////////////////////////////////////////////////////////////
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/QuadExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <Collections/collection.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace QuadCollectionTests
    {
        SpatialDomains::SegGeomSharedPtr CreateSegGeom(unsigned int id, 
                                                       SpatialDomains::PointGeomSharedPtr v0,
                                                       SpatialDomains::PointGeomSharedPtr v1)
        {
            SpatialDomains::PointGeomSharedPtr vertices[] = {v0, v1};
            SpatialDomains::SegGeomSharedPtr result(new SpatialDomains::SegGeom(id, 3, vertices));
            return result;
        }
        
        SpatialDomains::QuadGeomSharedPtr CreateQuad(SpatialDomains::PointGeomSharedPtr v0,
                                                     SpatialDomains::PointGeomSharedPtr v1,
                                                     SpatialDomains::PointGeomSharedPtr v2,
                                                     SpatialDomains::PointGeomSharedPtr v3)
        {
            Nektar::SpatialDomains::SegGeomSharedPtr e0 = CreateSegGeom(0, v0, v1);
            Nektar::SpatialDomains::SegGeomSharedPtr e1 = CreateSegGeom(1, v1, v2);
            Nektar::SpatialDomains::SegGeomSharedPtr e2 = CreateSegGeom(2, v2, v3);
            Nektar::SpatialDomains::SegGeomSharedPtr e3 = CreateSegGeom(3, v3, v0);
            
            Nektar::SpatialDomains::SegGeomSharedPtr edges[Nektar::SpatialDomains::QuadGeom::kNedges] =
                {
                    e0, e1, e2, e3
                };
            
            Nektar::StdRegions::Orientation edgeorient[Nektar::SpatialDomains::QuadGeom::kNedges] =
                {
                    Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                    Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                    Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
                    Nektar::SpatialDomains::SegGeom::GetEdgeOrientation(*edges[3], *edges[0])
                };
            
            SpatialDomains::QuadGeomSharedPtr quadGeom(new SpatialDomains::QuadGeom(0,edges,edgeorient));
            return quadGeom;
        }
        
        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(quadGeom);
            
            Collections::Collection c(Exp, GeomVec,Collections::eStdMat);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());
            

            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }
        
        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_StdMat_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(quadGeom);

            Collections::Collection c(Exp, GeomVec,Collections::eStdMat);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());

            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(quadGeom);
            
            Collections::Collection c(Exp, GeomVec,Collections::eIterPerExp);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());
            

            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }
        
        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_IterPerExp_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(quadGeom);

            Collections::Collection c(Exp, GeomVec,Collections::eIterPerExp);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());
            

            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }
        
        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;

            int nelmts = 1;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(quadGeom);
            }
            
            Collections::Collection c(Exp, GeomVec, Collections::eSumFac);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());
            
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_SumFac_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            
            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(quadGeom);
            }
            
            Collections::Collection c(Exp, GeomVec, Collections::eSumFac);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());
            
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_SumFac_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, quadGeom);

            int nelmts = 1;
            
            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(quadGeom);
            }
            
            Collections::Collection c(Exp, GeomVec, Collections::eSumFac);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());
            
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, quadGeom);

            int nelmts = 10;
            
            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(quadGeom);
            }
            
            Collections::Collection c(Exp, GeomVec, Collections::eSumFac);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());
            
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(quadGeom);

            Collections::Collection c(Exp, GeomVec,Collections::eStdMat);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nq);
            Array<OneD, NekDouble> coeffs1(Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            
            Exp->GetCoords(xc, yc, zc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }

            Exp->IProductWRTBase(phys, coeffs1);
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_StdMat_VariableP)
        {

            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));
            
            SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
            
            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            
            Nektar::LocalRegions::QuadExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, quadGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(quadGeom);
            
            Collections::Collection c(Exp, GeomVec,Collections::eStdMat);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nq);
            Array<OneD, NekDouble> coeffs1(Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            
            Exp->GetCoords(xc, yc, zc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }

            Exp->IProductWRTBase(phys, coeffs1);
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }        
    }
}
