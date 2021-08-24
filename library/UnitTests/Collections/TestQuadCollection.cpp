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

#include <LocalRegions/QuadExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <Collections/Collection.h>
#include <Collections/CollectionOptimisation.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

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

        SpatialDomains::QuadGeomSharedPtr quadGeom(new SpatialDomains::QuadGeom(0,edges));
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
        Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
        Array<OneD, NekDouble> phys2(Exp->GetTotPoints());

        Exp->BwdTrans(coeffs, phys1);
        c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

        double epsilon = 1.0e-8;
        for(int i = 0; i < phys1.size(); ++i)
        {
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
        }
    }


    BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_StdMat_VariableP_MultiElmt)
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        int nelmts = 10;
        for(int i = 0; i < nelmts; ++i)
        {
            CollExp.push_back(Exp);
        }

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);


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
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eIterPerExp);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
        Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
        Array<OneD, NekDouble> phys2(Exp->GetTotPoints());

        Exp->BwdTrans(coeffs, phys1);
        c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

        double epsilon = 1.0e-8;
        for(int i = 0; i < phys1.size(); ++i)
        {
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
        Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
        Array<OneD, NekDouble> phys2(Exp->GetTotPoints());

        Exp->BwdTrans(coeffs, phys1);
        c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

        double epsilon = 1.0e-8;
        for(int i = 0; i < phys1.size(); ++i)
        {
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadBwdTrans_MatrixFree_UniformP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 4;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        impTypes[Collections::eBwdTrans] = Collections::eMatrixFree;
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
        Array<OneD, NekDouble> physRef(Exp->GetTotPoints());
        Array<OneD, NekDouble> phys(Exp->GetTotPoints());


        Exp->BwdTrans(coeffs, physRef);
        c.ApplyOperator(Collections::eBwdTrans, coeffs, phys);

        double epsilon = 1.0e-8;
        for(int i = 0; i < physRef.size(); ++i)
        {
            BOOST_CHECK_CLOSE(physRef[i], phys[i], epsilon);
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;

        int nelmts = 1;
        for(int i = 0; i < nelmts; ++i)
        {
            CollExp.push_back(Exp);
        }

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

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
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;

        int nelmts = 10;
        for(int i = 0; i < nelmts; ++i)
        {
            CollExp.push_back(Exp);
        }

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

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
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

        int nelmts = 1;

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        for(int i = 0; i < nelmts; ++i)
        {
            CollExp.push_back(Exp);
        }

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eBwdTrans);

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
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

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
        c.Initialise(Collections::eBwdTrans);

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
            BOOST_CHECK_CLOSE(phys1[i], phys2[i], epsilon);
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);


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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_StdMat_VariableP_MultiElmt)
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

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
        c.Initialise(Collections::eIProductWRTBase);


        const int nq = Exp->GetTotPoints();
        Array<OneD, NekDouble> phys(nelmts*nq),tmp;
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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_IterPerExp_UniformP)
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }


    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_SumFac_UniformP)
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_SumFac_VariableP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 1.0, 0.0));

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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_MatrixFree_UniformP_Undeformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);

        const int nq = Exp->GetTotPoints();
        Array<OneD, NekDouble> phys(nq);
        Array<OneD, NekDouble> coeffsRef(Exp->GetNcoeffs());
        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs());

        Array<OneD, NekDouble> xc(nq), yc(nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys[i] = sin(xc[i])*cos(yc[i]);
        }

        Exp->IProductWRTBase(phys, coeffsRef);
        c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffs.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i],coeffs[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_MatrixFree_UniformP_Deformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 2.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);

        const int nq = Exp->GetTotPoints();
        Array<OneD, NekDouble> phys(nq);
        Array<OneD, NekDouble> coeffsRef(Exp->GetNcoeffs());
        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs());

        Array<OneD, NekDouble> xc(nq), yc(nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys[i] = sin(xc[i])*cos(yc[i]);
        }

        Exp->IProductWRTBase(phys, coeffsRef);
        c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffs.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i],coeffs[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTBase_MatrixFree_UniformP_Deformed_OverInt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 2.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 10;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTBase);

        const int nq = Exp->GetTotPoints();
        Array<OneD, NekDouble> phys(nq);
        Array<OneD, NekDouble> coeffsRef(Exp->GetNcoeffs());
        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs());

        Array<OneD, NekDouble> xc(nq), yc(nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys[i] = sin(xc[i])*cos(yc[i]);
        }

        Exp->IProductWRTBase(phys, coeffsRef);
        c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffs.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i],coeffs[i], epsilon);
        }
    }



    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_IterPerExp_UniformP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::ePhysDeriv);

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
        c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2 + nq);

        double epsilon = 1.0e-8;
        for(int i = 0; i < diff1.size(); ++i)
        {
            BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_IterPerExp_VariableP_MultiElmt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);
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
        c.Initialise(Collections::ePhysDeriv);

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
        Exp->PhysDeriv(phys, diff1, tmp = diff1 + nelmts*nq);
        for(int i = 1; i < nelmts; ++i)
        {
            Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
            Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                           tmp1 = diff1+(nelmts+i)*nq);
        }

        c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2 + nelmts*nq);

        double epsilon = 1.0e-8;
        for(int i = 0; i < diff1.size(); ++i)
        {
            diff1[i] = (std::abs(diff1[i]) < 1e-14)? 0.0: diff1[i];
            diff2[i] = (std::abs(diff2[i]) < 1e-14)? 0.0: diff2[i];
            BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_MatrixFree_UniformP_Undeformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 5;
        unsigned int numModes = 2;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
                basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::ePhysDeriv);

        const int nq = Exp->GetTotPoints();
        Array<OneD, NekDouble> xc(nq), yc(nq);
        Array<OneD, NekDouble> phys(nq), tmp, tmp1;
        Array<OneD, NekDouble> derivRef(2*nq);
        Array<OneD, NekDouble> deriv(2*nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys[i] = sin(xc[i])*cos(yc[i]);
        }

        Exp->PhysDeriv(phys, derivRef, tmp = derivRef+nq);
        c.ApplyOperator(Collections::ePhysDeriv, phys, deriv, tmp = deriv + nq);

        double epsilon = 1.0e-8;
        for (int i = 0; i < derivRef.size(); ++i)
        {
            derivRef[i] = (std::abs(derivRef[i]) < 1e-14)? 0.0: derivRef[i];
            deriv[i] = (std::abs(deriv[i]) < 1e-14)? 0.0: deriv[i];
            BOOST_CHECK_CLOSE(derivRef[i], deriv[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_MatrixFree_UniformP_Deformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -2.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 4;
        unsigned int numModes = 2;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
                basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::ePhysDeriv);

        const int nq = Exp->GetTotPoints();
        Array<OneD, NekDouble> xc(nq), yc(nq);
        Array<OneD, NekDouble> phys(nq), tmp, tmp1;
        Array<OneD, NekDouble> derivRef(2*nq);
        Array<OneD, NekDouble> deriv(2*nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys[i] = sin(xc[i])*cos(yc[i]);
        }

        Exp->PhysDeriv(phys, derivRef, tmp = derivRef + nq);
        c.ApplyOperator(Collections::ePhysDeriv, phys, deriv, tmp = deriv + nq);

        double epsilon = 1.0e-8;
        for (int i = 0; i < derivRef.size(); ++i)
        {
            derivRef[i] = (std::abs(derivRef[i]) < 1e-14)? 0.0: derivRef[i];
            deriv[i] = (std::abs(deriv[i]) < 1e-14)? 0.0: deriv[i];
            BOOST_CHECK_CLOSE(derivRef[i], deriv[i], epsilon);
        }
    }

    
    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_Directional_MatrixFree_UniformP_Undeformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 5;
        unsigned int numModes = 2;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
                basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::ePhysDeriv);

        const int nq = Exp->GetTotPoints();
        Array<OneD, NekDouble> xc(nq), yc(nq);
        Array<OneD, NekDouble> phys(nq), tmp, tmp1;
        Array<OneD, NekDouble> derivRef(2*nq);
        Array<OneD, NekDouble> deriv(2*nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys[i] = sin(xc[i])*cos(yc[i]);
        }

        Exp->PhysDeriv(0, phys, derivRef);
        Exp->PhysDeriv(1, phys, tmp = derivRef+nq);

        c.ApplyOperator(Collections::ePhysDeriv, 0, phys, deriv);
        c.ApplyOperator(Collections::ePhysDeriv, 1, phys, tmp = deriv + nq);

        double epsilon = 1.0e-8;
        for (int i = 0; i < derivRef.size(); ++i)
        {
            derivRef[i] = (std::abs(derivRef[i]) < 1e-14)? 0.0: derivRef[i];
            deriv[i] = (std::abs(deriv[i]) < 1e-14)? 0.0: deriv[i];
            BOOST_CHECK_CLOSE(derivRef[i], deriv[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_StdMat_UniformP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::ePhysDeriv);

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
        c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2 + nq);

        double epsilon = 1.0e-8;
        for(int i = 0; i < diff1.size(); ++i)
        {
            BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_StdMat_VariableP_MultiElmt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

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
        c.Initialise(Collections::ePhysDeriv);

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
        Exp->PhysDeriv(phys, diff1, tmp1 = diff1 + nelmts*nq);
        for(int i = 1; i < nelmts; ++i)
        {
            Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
            Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                           tmp1 = diff1+(nelmts+i)*nq);

        }

        c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2 + nelmts*nq);

        double epsilon = 1.0e-8;
        for(int i = 0; i < diff1.size(); ++i)
        {
            diff1[i] = (std::abs(diff1[i]) < 1e-14)? 0.0: diff1[i];
            diff2[i] = (std::abs(diff2[i]) < 1e-14)? 0.0: diff2[i];
            BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_SumFac_UniformP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::ePhysDeriv);

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
        c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2 + nq);

        double epsilon = 1.0e-8;
        for(int i = 0; i < diff1.size(); ++i)
        {
            BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadPhysDeriv_SumFac_VariableP_MultiElmt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

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
        c.Initialise(Collections::ePhysDeriv);

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
        Exp->PhysDeriv(phys, diff1, tmp1 = diff1 + nelmts*nq);
        for(int i = 1; i < nelmts; ++i)
        {
            Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
            Exp->PhysDeriv(phys, tmp = diff1+i*nq,
                           tmp1 = diff1+(nelmts+i)*nq);

        }

        c.ApplyOperator(Collections::ePhysDeriv, phys, diff2, tmp = diff2 + nelmts*nq);

        double epsilon = 1.0e-8;
        for(int i = 0; i < diff1.size(); ++i)
        {
            diff1[i] = (std::abs(diff1[i]) < 1e-14)? 0.0: diff1[i];
            diff2[i] = (std::abs(diff2[i]) < 1e-14)? 0.0: diff2[i];
            BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_IterPerExp_UniformP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTDerivBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }


    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_IterPerExp_VariableP_MultiElmt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.5, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(7, quadPointsTypeDir1);
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(5, quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,6,quadPointsKeyDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,4,quadPointsKeyDir2);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

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
        c.Initialise(Collections::eIProductWRTDerivBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_MatrixFree_UniformP_Undeformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     cref(CollExp, impTypes);
        cref.Initialise(Collections::eIProductWRTDerivBase);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTDerivBase);

        const int nq = Exp->GetTotPoints();
        const int nm = Exp->GetNcoeffs();
        Array<OneD, NekDouble> phys1(nq);
        Array<OneD, NekDouble> phys2(nq);
        Array<OneD, NekDouble> coeffsRef(nm);
        Array<OneD, NekDouble> coeffs(nm);

        Array<OneD, NekDouble> xc(nq), yc(nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys1[i] = sin(xc[i])*cos(yc[i]);
            phys2[i] = cos(xc[i])*sin(yc[i]);
        }

        // Standard routines
        Exp->IProductWRTDerivBase(0, phys1, coeffsRef);
        Exp->IProductWRTDerivBase(1, phys2, coeffs);
        Vmath::Vadd(nm,coeffsRef,1,coeffs,1,coeffsRef,1);

        cref.ApplyOperator(Collections::eIProductWRTDerivBase,
                        phys1, phys2, coeffs);
        c.ApplyOperator(Collections::eIProductWRTDerivBase,
                        phys1, phys2, coeffs);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffsRef.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i], coeffs[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_MatrixFree_UniformP_Deformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  3.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     cref(CollExp, impTypes);
        cref.Initialise(Collections::eIProductWRTDerivBase);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTDerivBase);

        const int nq = Exp->GetTotPoints();
        const int nm = Exp->GetNcoeffs();
        Array<OneD, NekDouble> phys1(nq);
        Array<OneD, NekDouble> phys2(nq);
        Array<OneD, NekDouble> coeffsRef(nm);
        Array<OneD, NekDouble> coeffs(nm);

        Array<OneD, NekDouble> xc(nq), yc(nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys1[i] = sin(xc[i])*cos(yc[i]);
            phys2[i] = cos(xc[i])*sin(yc[i]);
        }

        // Standard routines
        Exp->IProductWRTDerivBase(0, phys1, coeffsRef);
        Exp->IProductWRTDerivBase(1, phys2, coeffs);
        Vmath::Vadd(nm,coeffsRef,1,coeffs,1,coeffsRef,1);

        cref.ApplyOperator(Collections::eIProductWRTDerivBase,
                        phys1, phys2, coeffs);
        c.ApplyOperator(Collections::eIProductWRTDerivBase,
                        phys1, phys2, coeffs);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffsRef.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i], coeffs[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_MatrixFree_UniformP_Deformed_ThreeD)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u,
            1u,  3.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u,
            2u,  1.0, 1.0, 1.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u,
            3u, -1.0, 1.0, 1.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);
        int nelmts = 10;

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        for(int i = 0; i < nelmts; ++i)
        {
            CollExp.push_back(Exp);
        }

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTDerivBase);

        const int nq = Exp->GetTotPoints();
        const int nm = Exp->GetNcoeffs();
        Array<OneD, NekDouble> phys1(nelmts*nq);
        Array<OneD, NekDouble> phys2(nelmts*nq);
        Array<OneD, NekDouble> phys3(nelmts*nq);
        Array<OneD, NekDouble> coeffsRef(nelmts*nm);
        Array<OneD, NekDouble> coeffs(nelmts*nm);

        Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq), tmp;

        Exp->GetCoords(xc, yc, zc);

        for (int i = 0; i < nq; ++i)
        {
            phys1[i] = sin(xc[i])*cos(yc[i]);
            phys2[i] = cos(xc[i])*sin(yc[i]);
            phys3[i] = cos(xc[i])*sin(zc[i]);
        }

        for(int i = 1; i < nelmts; ++i)
        {
            Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
            Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
            Vmath::Vcopy(nq,phys2,1,tmp = phys3+i*nq,1);
        }

        for(int i = 0; i < nelmts; ++i)
        {
            // Standard routines
            Exp->IProductWRTDerivBase(0, phys1 + i*nq,
                                      tmp = coeffsRef + i*nm);
            Exp->IProductWRTDerivBase(1, phys2 + i*nq,
                                      tmp = coeffs + i*nm);
            Vmath::Vadd(nm,coeffsRef + i*nm, 1, coeffs + i*nm, 1,
                        tmp = coeffsRef + i*nm, 1);
            Exp->IProductWRTDerivBase(2, phys3 + i*nq,
                                      tmp = coeffs + i*nm);
            Vmath::Vadd(nm,coeffsRef + i*nm, 1, coeffs + i*nm, 1,
                        tmp = coeffsRef + i*nm, 1);
        }
        
        c.ApplyOperator(Collections::eIProductWRTDerivBase,
                        phys1, phys2, phys3, coeffs);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffsRef.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i], coeffs[i], epsilon);
        }
    }
    
    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_MatrixFree_UniformP_Deformed_OverInt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  3.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 10;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
            Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     cref(CollExp, impTypes);
        cref.Initialise(Collections::eIProductWRTDerivBase);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTDerivBase);

        const int nq = Exp->GetTotPoints();
        const int nm = Exp->GetNcoeffs();
        Array<OneD, NekDouble> phys1(nq);
        Array<OneD, NekDouble> phys2(nq);
        Array<OneD, NekDouble> coeffsRef(nm);
        Array<OneD, NekDouble> coeffs(nm);

        Array<OneD, NekDouble> xc(nq), yc(nq);

        Exp->GetCoords(xc, yc);

        for (int i = 0; i < nq; ++i)
        {
            phys1[i] = sin(xc[i])*cos(yc[i]);
            phys2[i] = cos(xc[i])*sin(yc[i]);
        }

        // Standard routines
        Exp->IProductWRTDerivBase(0, phys1, coeffsRef);
        Exp->IProductWRTDerivBase(1, phys2, coeffs);
        Vmath::Vadd(nm,coeffsRef,1,coeffs,1,coeffsRef,1);

        cref.ApplyOperator(Collections::eIProductWRTDerivBase,
                        phys1, phys2, coeffs);
        c.ApplyOperator(Collections::eIProductWRTDerivBase,
                        phys1, phys2, coeffs);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffsRef.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i], coeffs[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_StdMat_UniformP)
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

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eStdMat);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTDerivBase);


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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }
    
    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_StdMat_VariableP_MultiElmt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.5, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(7, quadPointsTypeDir1);
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(5, quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,6,quadPointsKeyDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,4,quadPointsKeyDir2);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

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
        c.Initialise(Collections::eIProductWRTDerivBase);


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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

 BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_SumFac_UniformP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0, 2.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir1);

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        CollExp.push_back(Exp);

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession, Collections::eSumFac);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise(Collections::eIProductWRTDerivBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_SumFac_VariableP_MultiElmt)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.5, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(7, quadPointsTypeDir1);
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(5, quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,6,quadPointsKeyDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,4,quadPointsKeyDir2);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1,
            basisKeyDir2);

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
        c.Initialise(Collections::eIProductWRTDerivBase);

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
            coeffs1[i] = (std::abs(coeffs1[i]) < 1e-14)? 0.0: coeffs1[i];
            coeffs2[i] = (std::abs(coeffs2[i]) < 1e-14)? 0.0: coeffs2[i];
            BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadIProductWRTDerivBase_SumFac_VariableP_MultiElmt_threedim)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(3u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(3u, 2u,  1.0,  1.0, 1.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(3u, 3u, -1.0,  1.0, 1.0));
        
        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);
        
        Nektar::LibUtilities::PointsType quadPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(7, quadPointsTypeDir1);
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir2(5, quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,6,quadPointsKeyDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir2(basisTypeDir1,4,quadPointsKeyDir2);
        
        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(basisKeyDir1, basisKeyDir2, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(basisKeyDir1, basisKeyDir2);

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
        c.Initialise(Collections::eIProductWRTDerivBase);
        
        const int nq = Exp->GetTotPoints();
        const int nm = Exp->GetNcoeffs();
        Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq), tmp,tmp1;
        Array<OneD, NekDouble> phys1(nelmts*nq);
        Array<OneD, NekDouble> phys2(nelmts*nq);
        Array<OneD, NekDouble> phys3(nelmts*nq);
        Array<OneD, NekDouble> coeffs1(nelmts*nm);
        Array<OneD, NekDouble> coeffs2(nelmts*nm);
        
        Exp->GetCoords(xc, yc, zc);

        for (int i = 0; i < nq; ++i)
        {
            phys1[i] = sin(xc[i])*cos(yc[i]);
            phys2[i] = cos(xc[i])*sin(yc[i]);
            phys2[i] = cos(xc[i])*sin(zc[i]);
        }
        for(int i = 1; i < nelmts; ++i)
        {
            Vmath::Vcopy(nq,phys1,1,tmp = phys1+i*nq,1);
            Vmath::Vcopy(nq,phys2,1,tmp = phys2+i*nq,1);
            Vmath::Vcopy(nq,phys3,1,tmp = phys3+i*nq,1);
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
            Exp->IProductWRTDerivBase(2, phys3 + i*nq,
                                      tmp1 = coeffs2 + i*nm);
            Vmath::Vadd(nm,coeffs1 +i*nm ,1,coeffs2 + i*nm ,1,
                        tmp = coeffs1 + i*nm,1);
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

    BOOST_AUTO_TEST_CASE(TestQuadHelmholtz_MatrixFree_UniformP)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u,
            0u, -1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u,
            1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u,
            2u,  1.0, 1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u,
            3u, -1.0, 1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        int nelmts = 10;

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        for(int i = 0; i < nelmts; ++i)
        {
            CollExp.push_back(Exp);
        }

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
                                                   Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise (Collections::eHelmholtz);
        StdRegions::ConstFactorMap factors;

        factors[StdRegions::eFactorLambda] = 1.5; 
        
        const int nm = Exp->GetNcoeffs();
        Array<OneD, NekDouble> coeffsIn(nelmts*nm);
        Array<OneD, NekDouble> coeffsRef(nelmts*nm);
        Array<OneD, NekDouble> coeffs(nelmts*nm), tmp;

        for (int i = 0; i < nm; ++i)
        {
            coeffsIn[i] = 1.0; 
        }

        for(int i = 1; i < nelmts; ++i)
        {
            Vmath::Vcopy(nm,coeffsIn,1,tmp = coeffsIn + i*nm,1);
        }

        StdRegions::StdMatrixKey mkey(StdRegions::eHelmholtz,
                                      Exp->DetShapeType(),
                                      *Exp, factors);
        
        for(int i = 0; i < nelmts; ++i)
        {            
            // Standard routines
            Exp->GeneralMatrixOp(coeffsIn + i*nm,
                                 tmp = coeffsRef + i*nm, mkey);
        }
        
        c.ApplyOperator(Collections::eHelmholtz,coeffsIn,coeffs,factors);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffsRef.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i], coeffs[i], epsilon);
        }
    }

    BOOST_AUTO_TEST_CASE(TestQuadHelmholtz_MatrixFree_UniformP_Deformed)
    {
        SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
        SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v2(new SpatialDomains::PointGeom(2u, 2u,  1.0,  1.0, 0.0));
        SpatialDomains::PointGeomSharedPtr v3(new SpatialDomains::PointGeom(2u, 3u, -1.0,  1.0, 0.0));

        SpatialDomains::QuadGeomSharedPtr quadGeom = CreateQuad(v0, v1, v2, v3);

        Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
            Nektar::LibUtilities::eGaussLobattoLegendre;
        Nektar::LibUtilities::BasisType basisTypeDir1 =
            Nektar::LibUtilities::eModified_A;
        unsigned int numQuadPoints = 6;
        unsigned int numModes = 5;
        const Nektar::LibUtilities::PointsKey quadPointsKeyDir1(numQuadPoints,
            quadPointsTypeDir1);
        const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,
            numModes, quadPointsKeyDir1);

        Nektar::LocalRegions::QuadExpSharedPtr Exp =
            MemoryManager<Nektar::LocalRegions::QuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1, quadGeom);

        Nektar::StdRegions::StdQuadExpSharedPtr stdExp =
            MemoryManager<Nektar::StdRegions::StdQuadExp>::AllocateSharedPtr(
            basisKeyDir1, basisKeyDir1);

        int nelmts = 10;

        std::vector<StdRegions::StdExpansionSharedPtr> CollExp;
        for(int i = 0; i < nelmts; ++i)
        {
            CollExp.push_back(Exp);
        }

        LibUtilities::SessionReaderSharedPtr dummySession;
        Collections::CollectionOptimisation colOpt(dummySession,
                                                   Collections::eMatrixFree);
        Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(stdExp);
        Collections::Collection     c(CollExp, impTypes);
        c.Initialise (Collections::eHelmholtz);
        StdRegions::ConstFactorMap factors;

        factors[StdRegions::eFactorLambda] = 0.0; 
        
        const int nm = Exp->GetNcoeffs();
        Array<OneD, NekDouble> coeffsIn(nelmts*nm);
        Array<OneD, NekDouble> coeffsRef(nelmts*nm);
        Array<OneD, NekDouble> coeffs(nelmts*nm), tmp;

        for (int i = 0; i < nm; ++i)
        {
            coeffsIn[i] = 1.0; 
        }

        for(int i = 1; i < nelmts; ++i)
        {
            Vmath::Vcopy(nm,coeffsIn,1,tmp = coeffsIn + i*nm,1);
        }

        StdRegions::StdMatrixKey mkey(StdRegions::eHelmholtz,
                                      Exp->DetShapeType(),
                                      *Exp, factors);
        
        for(int i = 0; i < nelmts; ++i)
        {            
            // Standard routines
            Exp->GeneralMatrixOp(coeffsIn + i*nm,
                                 tmp = coeffsRef + i*nm, mkey);
        }
        
        c.ApplyOperator(Collections::eHelmholtz,coeffsIn,coeffs,factors);

        double epsilon = 1.0e-8;
        for(int i = 0; i < coeffsRef.size(); ++i)
        {
            coeffsRef[i] = (std::abs(coeffsRef[i]) < 1e-14)? 0.0: coeffsRef[i];
            coeffs[i] = (std::abs(coeffs[i]) < 1e-14)? 0.0: coeffs[i];
            BOOST_CHECK_CLOSE(coeffsRef[i], coeffs[i], epsilon);
        }
    }
} // namespace QuadCollectionTests
} // namespace Nektar
