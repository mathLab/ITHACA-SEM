///////////////////////////////////////////////////////////////////////////////
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
/// The above copyright notice and this permission notice shall be included
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

#include <LocalRegions/HexExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <Collections/Collection.h>
#include <Collections/CollectionOptimisation.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace HexCollectionTests
    {
        SpatialDomains::SegGeomSharedPtr CreateSegGeom(unsigned int id,
            SpatialDomains::PointGeomSharedPtr v0,
            SpatialDomains::PointGeomSharedPtr v1)
        {
            SpatialDomains::PointGeomSharedPtr vertices[] = {v0, v1};
            SpatialDomains::SegGeomSharedPtr result(new SpatialDomains::SegGeom(id, 3, vertices));
            return result;
        }

        SpatialDomains::HexGeomSharedPtr CreateHex(
            SpatialDomains::PointGeomSharedPtr v0,
            SpatialDomains::PointGeomSharedPtr v1,
            SpatialDomains::PointGeomSharedPtr v2,
            SpatialDomains::PointGeomSharedPtr v3,
            SpatialDomains::PointGeomSharedPtr v4,
            SpatialDomains::PointGeomSharedPtr v5,
            SpatialDomains::PointGeomSharedPtr v6,
            SpatialDomains::PointGeomSharedPtr v7)
        {
            Nektar::SpatialDomains::SegGeomSharedPtr e0 = CreateSegGeom(0, v0, v1);
            Nektar::SpatialDomains::SegGeomSharedPtr e1 = CreateSegGeom(1, v1, v2);
            Nektar::SpatialDomains::SegGeomSharedPtr e2 = CreateSegGeom(2, v2, v3);
            Nektar::SpatialDomains::SegGeomSharedPtr e3 = CreateSegGeom(3, v3, v0);
            Nektar::SpatialDomains::SegGeomSharedPtr e4 = CreateSegGeom(4, v0, v4);
            Nektar::SpatialDomains::SegGeomSharedPtr e5 = CreateSegGeom(5, v1, v5);
            Nektar::SpatialDomains::SegGeomSharedPtr e6 = CreateSegGeom(6, v2, v6);
            Nektar::SpatialDomains::SegGeomSharedPtr e7 = CreateSegGeom(7, v3, v7);
            Nektar::SpatialDomains::SegGeomSharedPtr e8 = CreateSegGeom(8, v4, v5);
            Nektar::SpatialDomains::SegGeomSharedPtr e9 = CreateSegGeom(9, v5, v6);
            Nektar::SpatialDomains::SegGeomSharedPtr e10 = CreateSegGeom(10, v6, v7);
            Nektar::SpatialDomains::SegGeomSharedPtr e11 = CreateSegGeom(11, v4, v7);

            Nektar::SpatialDomains::SegGeomSharedPtr edgesF0[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e0, e1, e2, e3
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF1[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e0, e5, e8, e4
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF2[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e1, e6, e9, e5
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF3[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e2, e6, e10, e7
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF4[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e3, e7, e11, e4
            };
            Nektar::SpatialDomains::SegGeomSharedPtr edgesF5[Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e8, e9, e10, e11
            };

            Nektar::SpatialDomains::QuadGeomSharedPtr face0(new SpatialDomains::QuadGeom(0, edgesF0));
            Nektar::SpatialDomains::QuadGeomSharedPtr face1(new SpatialDomains::QuadGeom(1, edgesF1));
            Nektar::SpatialDomains::QuadGeomSharedPtr face2(new SpatialDomains::QuadGeom(2, edgesF2));
            Nektar::SpatialDomains::QuadGeomSharedPtr face3(new SpatialDomains::QuadGeom(3, edgesF3));
            Nektar::SpatialDomains::QuadGeomSharedPtr face4(new SpatialDomains::QuadGeom(4, edgesF4));
            Nektar::SpatialDomains::QuadGeomSharedPtr face5(new SpatialDomains::QuadGeom(5, edgesF5));

            Nektar::SpatialDomains::QuadGeomSharedPtr qfaces[] = {face0, face1, face2, face3, face4, face5};
            SpatialDomains::HexGeomSharedPtr hexGeom(new SpatialDomains::HexGeom(0,qfaces));
            return hexGeom;
        }

        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());

            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_StdMat_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());


            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());


            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_IterPerExp_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());


            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_IterPerExp_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);


            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());

            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }

            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_NoCollection_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession,
                                                       Collections::eNoCollection);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);


            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1 (Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2 (Exp->GetTotPoints());

            Exp->BwdTrans(coeffs, phys1);

            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            int nelmts = 1;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());

            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_SumFac_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());

            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_SumFac_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 1;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());

            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexBwdTrans_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());

            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_StdMat_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);


            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_NoCollection_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);


            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eNoCollection);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_VariableP_CollAll)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eGLL_Lagrange;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(4, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_VariableP_CollDir02)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eGLL_Lagrange;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(4, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_VariableP_CollDir12)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eGLL_Lagrange;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_StdMat_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);


            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq,
                                     tmp = coeffs1 + i*Exp->GetNcoeffs());
            }

            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_IterPerExp_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);


            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq,
                                     tmp = coeffs1 + i*Exp->GetNcoeffs());
            }

            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }



        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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

            double epsilon = 1.0e-6;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                // clamp values below 1e-16 to zero
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-16)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-16)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

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

            double epsilon = 1.0e-6;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                // clamp values below 1e-16 to zero
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-16)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-16)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq, tmp = coeffs1 + i*Exp->GetNcoeffs());
            }
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-6;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                // clamp values below 1e-16 to zero
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-16)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-16)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(7, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(9, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq, tmp = coeffs1 + i*Exp->GetNcoeffs());
            }
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-4;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                // clamp values below 1e-14 to zero
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_VariableP_MultiElmt_CollDir02)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eGLL_Lagrange;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(4, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq, tmp = coeffs1 + i*Exp->GetNcoeffs());
            }
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-4;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                // clamp values below 1e-14 to zero
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTBase_SumFac_VariableP_MultiElmt_CollDir12)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eGLL_Lagrange;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
                CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq, tmp = coeffs1 + i*Exp->GetNcoeffs());
            }
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-4;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                // clamp values below 1e-14 to zero
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexPhysDeriv_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            Array<OneD, NekDouble> phys(nq),tmp,tmp1,tmp2;
            Array<OneD, NekDouble> diff1(3*nq);
            Array<OneD, NekDouble> diff2(3*nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }

            Exp->PhysDeriv(phys,diff1,tmp = diff1+nq, tmp1 = diff1+2*nq);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2,tmp = diff2+nq, tmp2 = diff2+2*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexPhysDeriv_IterPerExp_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp,tmp1,tmp2;
            Array<OneD, NekDouble> diff1(3*nelmts*nq);
            Array<OneD, NekDouble> diff2(3*nelmts*nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->PhysDeriv(phys, tmp = diff1,
                           tmp1 = diff1+(nelmts)*nq,
                           tmp2 = diff1+(2*nelmts)*nq);
            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                               tmp1 = diff1+(nelmts+i)*nq,
                               tmp2 = diff1+(2*nelmts+i)*nq);
            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2,tmp = diff2 + nelmts*nq,
                            tmp2 = diff2+2*nelmts*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexPhysDeriv_NoCollection_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eNoCollection);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp,tmp1,tmp2;
            Array<OneD, NekDouble> diff1(3*nelmts*nq);
            Array<OneD, NekDouble> diff2(3*nelmts*nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->PhysDeriv(phys, tmp = diff1,
                               tmp1 = diff1+(nelmts)*nq,
                               tmp2 = diff1+(2*nelmts)*nq);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                               tmp1 = diff1+(nelmts+i)*nq,
                               tmp2 = diff1+(2*nelmts+i)*nq);
            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2,tmp = diff2 + nelmts*nq,
                            tmp2 = diff2+2*nelmts*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexPhysDeriv_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 4;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,3,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            Array<OneD, NekDouble> phys(nq),tmp,tmp1,tmp2;
            Array<OneD, NekDouble> diff1(3*nq);
            Array<OneD, NekDouble> diff2(3*nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }

            Exp->PhysDeriv(phys,diff1,tmp = diff1+nq, tmp1 = diff1+2*nq);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2,tmp = diff2+nq, tmp2 = diff2+2*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexPhysDeriv_StdMat_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp,tmp1,tmp2;
            Array<OneD, NekDouble> diff1(3*nelmts*nq);
            Array<OneD, NekDouble> diff2(3*nelmts*nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->PhysDeriv(phys, tmp = diff1,
                           tmp1 = diff1+(nelmts)*nq,
                           tmp2 = diff1+(2*nelmts)*nq);
            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                               tmp1 = diff1+(nelmts+i)*nq,
                               tmp2 = diff1+(2*nelmts+i)*nq);
            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2,tmp = diff2 + nelmts*nq,
                            tmp2 = diff2+2*nelmts*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexPhysDeriv_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            Array<OneD, NekDouble> phys(nq),tmp,tmp1,tmp2;
            Array<OneD, NekDouble> diff1(3*nq);
            Array<OneD, NekDouble> diff2(3*nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }

            Exp->PhysDeriv(phys,diff1,tmp = diff1+nq, tmp1 = diff1+2*nq);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2,tmp = diff2+nq, tmp2 = diff2+2*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexPhysDeriv_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp,tmp1,tmp2;
            Array<OneD, NekDouble> diff1(3*nelmts*nq);
            Array<OneD, NekDouble> diff2(3*nelmts*nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
            }
            Exp->PhysDeriv(phys, tmp = diff1,
                           tmp1 = diff1+(nelmts)*nq,
                           tmp2 = diff1+(2*nelmts)*nq);
            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                               tmp1 = diff1+(nelmts+i)*nq,
                               tmp2 = diff1+(2*nelmts+i)*nq);
            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2,tmp = diff2 + nelmts*nq,
                            tmp2 = diff2+2*nelmts*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTDerivBase_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nq,   0.0);
            Array<OneD, NekDouble> phys2(nq,   0.0);
            Array<OneD, NekDouble> phys3(nq,   0.0);
            Array<OneD, NekDouble> coeffs1(nm, 0.0);
            Array<OneD, NekDouble> coeffs2(nm, 0.0);

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            // Standard routines
            Exp->IProductWRTDerivBase(0, phys1, coeffs1);
            Exp->IProductWRTDerivBase(1, phys2, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);
            Exp->IProductWRTDerivBase(2, phys3, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, phys3, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < nm; ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTDerivBase_IterPerExp_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys2(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys3(nelmts*nq,   0.0);
            Array<OneD, NekDouble> coeffs1(nelmts*nm, 0.0);
            Array<OneD, NekDouble> coeffs2(nelmts*nm, 0.0);
            Array<OneD, NekDouble> tmp;

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
                Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
                Vmath::Vcopy(nq,phys3,1,tmp = phys3+i*nq,1);
            }

            // Standard routines
            for(int i = 0; i < nelmts; ++i)
            {

                Exp->IProductWRTDerivBase(0, phys1 + i*nq, tmp = coeffs1 + i*nm);
                Exp->IProductWRTDerivBase(1, phys2 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
                Exp->IProductWRTDerivBase(2, phys3 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
            }

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, phys3, coeffs2);

            double epsilon = 1.0e-6;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTDerivBase_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nq,   0.0);
            Array<OneD, NekDouble> phys2(nq,   0.0);
            Array<OneD, NekDouble> phys3(nq,   0.0);
            Array<OneD, NekDouble> coeffs1(nm, 0.0);
            Array<OneD, NekDouble> coeffs2(nm, 0.0);

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            // Standard routines
            Exp->IProductWRTDerivBase(0, phys1, coeffs1);
            Exp->IProductWRTDerivBase(1, phys2, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);
            Exp->IProductWRTDerivBase(2, phys3, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, phys3, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestHexIProductWRTDerivBase_StdMat_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys2(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys3(nelmts*nq,   0.0);
            Array<OneD, NekDouble> coeffs1(nelmts*nm, 0.0);
            Array<OneD, NekDouble> coeffs2(nelmts*nm, 0.0);
            Array<OneD, NekDouble> tmp;

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            for (int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
                Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
                Vmath::Vcopy(nq,phys3,1,tmp = phys3+i*nq,1);
            }

            // Standard routines
            for(int i = 0; i < nelmts; ++i)
            {

                Exp->IProductWRTDerivBase(0, phys1 + i*nq, tmp = coeffs1 + i*nm);
                Exp->IProductWRTDerivBase(1, phys2 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
                Exp->IProductWRTDerivBase(2, phys3 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
            }

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, phys3, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTDerivBase_NoCollection_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eNoCollection);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys2(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys3(nelmts*nq,   0.0);
            Array<OneD, NekDouble> coeffs1(nelmts*nm, 0.0);
            Array<OneD, NekDouble> coeffs2(nelmts*nm, 0.0);
            Array<OneD, NekDouble> tmp;
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            for (int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
                Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
                Vmath::Vcopy(nq,phys3,1,tmp = phys3+i*nq,1);
            }

            // Standard routines
            for(int i = 0; i < nelmts; ++i)
            {

                Exp->IProductWRTDerivBase(0, phys1 + i*nq, tmp = coeffs1 + i*nm);
                Exp->IProductWRTDerivBase(1, phys2 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
                Exp->IProductWRTDerivBase(2, phys3 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
            }

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, phys3, coeffs2);

            double epsilon = 1.0e-6;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTDerivBase_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numQuadPoints = 6;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nq,   0.0);
            Array<OneD, NekDouble> phys2(nq,   0.0);
            Array<OneD, NekDouble> phys3(nq,   0.0);
            Array<OneD, NekDouble> coeffs1(nm, 0.0);
            Array<OneD, NekDouble> coeffs2(nm, 0.0);

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            // Standard routines
            Exp->IProductWRTDerivBase(0, phys1, coeffs1);
            Exp->IProductWRTDerivBase(1, phys2, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);
            Exp->IProductWRTDerivBase(2, phys3, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, phys3, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestHexIProductWRTDerivBase_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.5, -1.5, -1.5));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
            SpatialDomains::PointGeomSharedPtr v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
            SpatialDomains::PointGeomSharedPtr v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

            SpatialDomains::HexGeomSharedPtr hexGeom = CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

            Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(5, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(6, quadPointsTypeDir1);
            const Nektar::LibUtilities::PointsKey quadPointsKeyDir3(8, quadPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,6,quadPointsKeyDir2);
            const Nektar::LibUtilities::BasisKey basisKeyDir3(basisTypeDir1,8,quadPointsKeyDir3);

            Nektar::LocalRegions::HexExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3, hexGeom);

            Nektar::StdRegions::StdHexExpSharedPtr stdExp =
                MemoryManager<Nektar::StdRegions::StdHexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, basisKeyDir3);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys2(nelmts*nq,   0.0);
            Array<OneD, NekDouble> phys3(nelmts*nq,   0.0);
            Array<OneD, NekDouble> coeffs1(nelmts*nm, 0.0);
            Array<OneD, NekDouble> coeffs2(nelmts*nm, 0.0);
            Array<OneD, NekDouble> tmp;

            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            for (int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
                Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
                Vmath::Vcopy(nq,phys3,1,tmp = phys3+i*nq,1);
            }

            // Standard routines
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->IProductWRTDerivBase(0, phys1 + i*nq, tmp = coeffs1 + i*nm);
                Exp->IProductWRTDerivBase(1, phys2 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
                Exp->IProductWRTDerivBase(2, phys3 + i*nq, tmp = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1+i*nm ,1,coeffs2+i*nm ,1,tmp = coeffs1 + i*nm,1);
            }

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, phys3, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }
    }
}
