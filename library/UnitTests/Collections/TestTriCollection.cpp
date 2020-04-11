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

#include <LocalRegions/TriExp.h>
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
    namespace TriCollectionTests
    {
        SpatialDomains::SegGeomSharedPtr CreateSegGeom(unsigned int id,
                                                       SpatialDomains::PointGeomSharedPtr v0,
                                                       SpatialDomains::PointGeomSharedPtr v1)
        {
            SpatialDomains::PointGeomSharedPtr vertices[] = {v0, v1};
            SpatialDomains::SegGeomSharedPtr result(new SpatialDomains::SegGeom(id, 2, vertices));
            return result;
        }

        SpatialDomains::TriGeomSharedPtr CreateTri(SpatialDomains::PointGeomSharedPtr   v0,
                                                     SpatialDomains::PointGeomSharedPtr v1,
                                                     SpatialDomains::PointGeomSharedPtr v2)
        {
            Nektar::SpatialDomains::SegGeomSharedPtr e0 = CreateSegGeom(0, v0, v1);
            Nektar::SpatialDomains::SegGeomSharedPtr e1 = CreateSegGeom(1, v1, v2);
            Nektar::SpatialDomains::SegGeomSharedPtr e2 = CreateSegGeom(2, v2, v0);

            Nektar::SpatialDomains::SegGeomSharedPtr edges[Nektar::SpatialDomains::TriGeom::kNedges] =
                {
                    e0, e1, e2
                };

            SpatialDomains::TriGeomSharedPtr triGeom(new SpatialDomains::TriGeom(0,edges));
            return triGeom;
        }

        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

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

        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_StdMat_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

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



        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_StdMat_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
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


        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
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

        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_IterPerExp_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);


            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
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

        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;

            int nelmts = 1;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
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

        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_SumFac_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;

            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
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


        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_SumFac_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);


            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 1;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
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

        BOOST_AUTO_TEST_CASE(TestTriBwdTrans_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
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


        BOOST_AUTO_TEST_CASE(TestTriIProductWRTBase_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nq);
            Array<OneD, NekDouble> coeffs1(Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->IProductWRTBase(phys, coeffs1);
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriIProductWRTBase_StdMat_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nq);
            Array<OneD, NekDouble> coeffs1(Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->IProductWRTBase(phys, coeffs1);
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriIProductWRTBase_StdMat_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }



            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq, tmp = coeffs1 + i*Exp->GetNcoeffs());
            }

            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriIProductWRTBase_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nq);
            Array<OneD, NekDouble> coeffs1(Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->IProductWRTBase(phys, coeffs1);
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestTriIProductWRTBase_SumFac_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nq);
            Array<OneD, NekDouble> coeffs1(Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->IProductWRTBase(phys, coeffs1);
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestTriIProductWRTBase_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> phys(nelmts*nq), tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*Exp->GetNcoeffs());
            Array<OneD, NekDouble> coeffs2(nelmts*Exp->GetNcoeffs());

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }
            Exp->IProductWRTBase(phys, coeffs1);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,&phys[0],1,&phys[i*nq],1);
                Exp->IProductWRTBase(phys +i*nq, tmp = coeffs1 + i*Exp->GetNcoeffs());
            }

            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriPhysDeriv_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nq),tmp,tmp1;
            Array<OneD, NekDouble> diff1(2*nq);
            Array<OneD, NekDouble> diff2(2*nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->PhysDeriv(phys,diff1,tmp = diff1+nq);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2+nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestTriPhysDeriv_IterPerExp_VariableP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nq),tmp,tmp1;
            Array<OneD, NekDouble> diff1(2*nq);
            Array<OneD, NekDouble> diff2(2*nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->PhysDeriv(phys,diff1,tmp = diff1+nq);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2+nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriPhysDeriv_IterPerExp_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp,tmp1;
            Array<OneD, NekDouble> diff1(2*nelmts*nq);
            Array<OneD, NekDouble> diff2(2*nelmts*nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }
            Exp->PhysDeriv(phys, diff1,
                           tmp1 = diff1+(nelmts)*nq);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                               tmp1 = diff1+(nelmts+i)*nq);
            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2+nelmts*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestTriPhysDeriv_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nq),tmp,tmp1;
            Array<OneD, NekDouble> diff1(2*nq);
            Array<OneD, NekDouble> diff2(2*nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->PhysDeriv(phys,diff1,tmp = diff1+nq);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2+nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestTriPhysDeriv_StdMat_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp,tmp1;
            Array<OneD, NekDouble> diff1(2*nelmts*nq);
            Array<OneD, NekDouble> diff2(2*nelmts*nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }
            Exp->PhysDeriv(phys, tmp = diff1,
                           tmp1 = diff1+(nelmts)*nq);

            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                               tmp1 = diff1+(nelmts+i)*nq);

            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2+nelmts*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriPhysDeriv_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nq),tmp,tmp1;
            Array<OneD, NekDouble> diff1(2*nq);
            Array<OneD, NekDouble> diff2(2*nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }

            Exp->PhysDeriv(phys,diff1,tmp = diff1+nq);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2+nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriPhysDeriv_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp,tmp1;
            Array<OneD, NekDouble> diff1(2*nelmts*nq);
            Array<OneD, NekDouble> diff2(2*nelmts*nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }
            Exp->PhysDeriv(phys, tmp = diff1,
                           tmp1 = diff1+(nelmts)*nq);
            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                               tmp1 = diff1+(nelmts+i)*nq);

            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2+nelmts*nq);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.size(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriIProductWRTDerivBase_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nq);
            Array<OneD, NekDouble> phys2(nq);
            Array<OneD, NekDouble> coeffs1(nm);
            Array<OneD, NekDouble> coeffs2(nm);

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i]);
            }

            // Standard routines
            Exp->IProductWRTDerivBase(0, phys1, coeffs1);
            Exp->IProductWRTDerivBase(1, phys2, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriIProductWRTDerivBase_IterPerExp_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq), yc(nq),tmp,tmp1;
            Array<OneD, NekDouble> phys1(nelmts*nq);
            Array<OneD, NekDouble> phys2(nelmts*nq);
            Array<OneD, NekDouble> coeffs1(nelmts*nm);
            Array<OneD, NekDouble> coeffs2(nelmts*nm);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i]);
            }
            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
                Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
            }

            for(int i = 0; i < nelmts; ++i)
            {
                // Standard routines
                Exp->IProductWRTDerivBase(0, phys1 + i*nq,
                                          tmp  = coeffs1 + i*nm);
                Exp->IProductWRTDerivBase(1, phys2 + i*nq,
                                          tmp1 = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1 +i*nm ,1,coeffs2 + i*nm ,1,
                            tmp = coeffs1 + i*nm,1);
            }

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1, phys2, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestTriIProductWRTDerivBase_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nq);
            Array<OneD, NekDouble> phys2(nq);
            Array<OneD, NekDouble> coeffs1(nm);
            Array<OneD, NekDouble> coeffs2(nm);

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i]);
            }

            // Standard routines
            Exp->IProductWRTDerivBase(0, phys1, coeffs1);
            Exp->IProductWRTDerivBase(1, phys2, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1,
                            phys2, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriIProductWRTDerivBase_StdMat_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }


            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);


            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq), yc(nq),tmp,tmp1;
            Array<OneD, NekDouble> phys1(nelmts*nq);
            Array<OneD, NekDouble> phys2(nelmts*nq);
            Array<OneD, NekDouble> coeffs1(nelmts*nm);
            Array<OneD, NekDouble> coeffs2(nelmts*nm);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i]);
            }
            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
                Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
            }

            for(int i = 0; i < nelmts; ++i)
            {
                // Standard routines
                Exp->IProductWRTDerivBase(0, phys1 + i*nq,
                                          tmp  = coeffs1 + i*nm);
                Exp->IProductWRTDerivBase(1, phys2 + i*nq,
                                          tmp1 = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1 +i*nm ,1,coeffs2 + i*nm ,1,
                            tmp = coeffs1 + i*nm,1);
            }

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1, phys2, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestTriIProductWRTDerivBase_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(5, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,4,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            CollExp.push_back(Exp);

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> phys1(nq);
            Array<OneD, NekDouble> phys2(nq);
            Array<OneD, NekDouble> coeffs1(nm);
            Array<OneD, NekDouble> coeffs2(nm);

            Array<OneD, NekDouble> xc(nq), yc(nq);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i]);
            }

            // Standard routines
            Exp->IProductWRTDerivBase(0, phys1, coeffs1);
            Exp->IProductWRTDerivBase(1, phys2, coeffs2);
            Vmath::Vadd(nm,coeffs1,1,coeffs2,1,coeffs1,1);

            c.ApplyOperator(Collections::eIProductWRTDerivBase,
                            phys1, phys2, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.size(); ++i)
            {
                coeffs1[i] = (fabs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
                coeffs2[i] = (fabs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestTriIProductWRTDerivBase_SumFac_VariableP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u, -1.0,  1.0, 0.0));

            SpatialDomains::TriGeomSharedPtr triGeom = CreateTri(v0, v1, v2);

            Nektar::LibUtilities::PointsType triPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir1(5, triPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,triPointsKeyDir1);

            Nektar::LibUtilities::PointsType triPointsTypeDir2 = Nektar::LibUtilities::eGaussRadauMAlpha1Beta0;
            const Nektar::LibUtilities::PointsKey triPointsKeyDir2(7, triPointsTypeDir2);
            Nektar::LibUtilities::BasisType       basisTypeDir2 = Nektar::LibUtilities::eModified_B;
            const Nektar::LibUtilities::BasisKey  basisKeyDir2(basisTypeDir2,6,triPointsKeyDir2);

            Nektar::LocalRegions::TriExpSharedPtr Exp =
                MemoryManager<Nektar::LocalRegions::TriExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir2, triGeom);

            int nelmts = 10;

            std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
            for(int i = 0; i < nelmts; ++i)
            {
                CollExp.push_back(Exp);
            }

            LibUtilities::SessionReaderSharedPtr dummySession;
            Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
            Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(Exp);
            Collections::Collection     c(CollExp, impTypes);

            const int nq = Exp->GetTotPoints();
            const int nm = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq), yc(nq),tmp,tmp1;
            Array<OneD, NekDouble> phys1(nelmts*nq);
            Array<OneD, NekDouble> phys2(nelmts*nq);
            Array<OneD, NekDouble> coeffs1(nelmts*nm);
            Array<OneD, NekDouble> coeffs2(nelmts*nm);

            Exp->GetCoords(xc, yc);

            for (int i = 0; i < nq; ++i)
            {
                phys1[i] = sin(xc[i])*cos(yc[i]);
                phys2[i] = cos(xc[i])*sin(yc[i]);
            }
            for(int i = 1; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
                Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
            }

            for(int i = 0; i < nelmts; ++i)
            {
                // Standard routines
                Exp->IProductWRTDerivBase(0, phys1 + i*nq,
                                          tmp  = coeffs1 + i*nm);
                Exp->IProductWRTDerivBase(1, phys2 + i*nq,
                                          tmp1 = coeffs2 + i*nm);
                Vmath::Vadd(nm,coeffs1 +i*nm ,1,coeffs2 + i*nm ,1,
                            tmp = coeffs1 + i*nm,1);
            }

            c.ApplyOperator(Collections::eIProductWRTDerivBase, phys1, phys2, coeffs2);

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
