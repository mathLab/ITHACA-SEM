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

#include <LocalRegions/HexExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <LocalRegions/Expansion3D.h>
#include <StdRegions/StdRegions.hpp>
#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
namespace Expansion3DTests
{

    SpatialDomains::SegGeomSharedPtr CreateSegGeom
         (unsigned int id,
          SpatialDomains::PointGeomSharedPtr v0,
          SpatialDomains::PointGeomSharedPtr v1)
    {
        SpatialDomains::PointGeomSharedPtr vertices[] = {v0, v1};
        SpatialDomains::SegGeomSharedPtr result(new SpatialDomains::SegGeom(id, 3, vertices));
        return result;
    }
    
    SpatialDomains::HexGeomSharedPtr CreateHex
           (SpatialDomains::PointGeomSharedPtr v0,
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
        
        Nektar::SpatialDomains::SegGeomSharedPtr edgesF0
            [Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
             e0, e1, e2, e3
            };
        Nektar::SpatialDomains::SegGeomSharedPtr edgesF1
            [Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
             e0, e5, e8, e4
            };
        Nektar::SpatialDomains::SegGeomSharedPtr edgesF2
            [Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e1, e6, e9, e5
            };
        Nektar::SpatialDomains::SegGeomSharedPtr edgesF3
            [Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e2, e6, e10, e7
            };
        Nektar::SpatialDomains::SegGeomSharedPtr edgesF4
            [Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e3, e7, e11, e4
            };
        Nektar::SpatialDomains::SegGeomSharedPtr edgesF5
            [Nektar::SpatialDomains::QuadGeom::kNedges] =
            {
                e8, e9, e10, e11
            };

        Nektar::SpatialDomains::QuadGeomSharedPtr face0
            (new SpatialDomains::QuadGeom(0, edgesF0));
        Nektar::SpatialDomains::QuadGeomSharedPtr face1
            (new SpatialDomains::QuadGeom(1, edgesF1));
        Nektar::SpatialDomains::QuadGeomSharedPtr face2
            (new SpatialDomains::QuadGeom(2, edgesF2));
        Nektar::SpatialDomains::QuadGeomSharedPtr face3
            (new SpatialDomains::QuadGeom(3, edgesF3));
        Nektar::SpatialDomains::QuadGeomSharedPtr face4
            (new SpatialDomains::QuadGeom(4, edgesF4));
        Nektar::SpatialDomains::QuadGeomSharedPtr face5
            (new SpatialDomains::QuadGeom(5, edgesF5));
       
        Nektar::SpatialDomains::QuadGeomSharedPtr qfaces[] =
             {face0, face1, face2, face3, face4, face5};
        SpatialDomains::HexGeomSharedPtr hexGeom
            (new SpatialDomains::HexGeom(0,qfaces));
        return hexGeom;
    }

BOOST_AUTO_TEST_CASE(TestReOrientQuadFacePhysMap)
{
    SpatialDomains::PointGeomSharedPtr
        v0(new SpatialDomains::PointGeom(3u, 0u, -1.0, -1.0, -1.0));
    SpatialDomains::PointGeomSharedPtr
        v1(new SpatialDomains::PointGeom(3u, 1u, 1.0, -1.0, -1.0));
    SpatialDomains::PointGeomSharedPtr
        v2(new SpatialDomains::PointGeom(3u, 2u, 1.0, 1.0, -1.0));
    SpatialDomains::PointGeomSharedPtr
        v3(new SpatialDomains::PointGeom(3u, 3u, -1.0, 1.0, -1.0));
    SpatialDomains::PointGeomSharedPtr
        v4(new SpatialDomains::PointGeom(3u, 4u, -1.0, -1.0, 1.0));
    SpatialDomains::PointGeomSharedPtr
        v5(new SpatialDomains::PointGeom(3u, 5u, 1.0, -1.0, 1.0));
    SpatialDomains::PointGeomSharedPtr
        v6(new SpatialDomains::PointGeom(3u, 6u, 1.0, 1.0, 1.0));
    SpatialDomains::PointGeomSharedPtr
        v7(new SpatialDomains::PointGeom(3u, 7u, -1.0, 1.0, 1.0));

    SpatialDomains::HexGeomSharedPtr hexGeom =
        CreateHex(v0, v1, v2, v3, v4, v5, v6, v7);

    Nektar::LibUtilities::PointsType quadPointsTypeDir1 =
        Nektar::LibUtilities::eGaussLobattoLegendre;
    Nektar::LibUtilities::BasisType basisTypeDir1 =
        Nektar::LibUtilities::eModified_A;
    unsigned int numQuadPoints = 6;
    const Nektar::LibUtilities::PointsKey
        quadPointsKeyDir1(numQuadPoints, quadPointsTypeDir1);
    const Nektar::LibUtilities::BasisKey
        basisKeyDir1(basisTypeDir1,4,quadPointsKeyDir1);

    // Dummy 3D expansion
    Nektar::LocalRegions::HexExpSharedPtr exp3d =
        MemoryManager<Nektar::LocalRegions::HexExp>::AllocateSharedPtr(basisKeyDir1,
                basisKeyDir1, basisKeyDir1, hexGeom);
    
    // Initialiaze 3*3 face id map for quad
    int nq0 = 3;
    int nq1 = 3;
    Array<OneD,int> idmap(nq0*nq1, -1);

    // Test different orientations
    StdRegions::Orientation orient;

    orient = StdRegions::eDir1FwdDir1_Dir2FwdDir2;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 0);
    BOOST_CHECK_EQUAL(idmap[1], 1);
    BOOST_CHECK_EQUAL(idmap[2], 2);
    BOOST_CHECK_EQUAL(idmap[3], 3);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 5);
    BOOST_CHECK_EQUAL(idmap[6], 6);
    BOOST_CHECK_EQUAL(idmap[7], 7);
    BOOST_CHECK_EQUAL(idmap[8], 8);

    orient = StdRegions::eDir1FwdDir1_Dir2BwdDir2;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 6);
    BOOST_CHECK_EQUAL(idmap[1], 7);
    BOOST_CHECK_EQUAL(idmap[2], 8);
    BOOST_CHECK_EQUAL(idmap[3], 3);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 5);
    BOOST_CHECK_EQUAL(idmap[6], 0);
    BOOST_CHECK_EQUAL(idmap[7], 1);
    BOOST_CHECK_EQUAL(idmap[8], 2);

    orient = StdRegions::eDir1BwdDir1_Dir2FwdDir2;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 2);
    BOOST_CHECK_EQUAL(idmap[1], 1);
    BOOST_CHECK_EQUAL(idmap[2], 0);
    BOOST_CHECK_EQUAL(idmap[3], 5);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 3);
    BOOST_CHECK_EQUAL(idmap[6], 8);
    BOOST_CHECK_EQUAL(idmap[7], 7);
    BOOST_CHECK_EQUAL(idmap[8], 6);

    orient = StdRegions::eDir1BwdDir1_Dir2BwdDir2;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 8);
    BOOST_CHECK_EQUAL(idmap[1], 7);
    BOOST_CHECK_EQUAL(idmap[2], 6);
    BOOST_CHECK_EQUAL(idmap[3], 5);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 3);
    BOOST_CHECK_EQUAL(idmap[6], 2);
    BOOST_CHECK_EQUAL(idmap[7], 1);
    BOOST_CHECK_EQUAL(idmap[8], 0);

    orient = StdRegions::eDir1FwdDir2_Dir2FwdDir1;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 0);
    BOOST_CHECK_EQUAL(idmap[1], 3);
    BOOST_CHECK_EQUAL(idmap[2], 6);
    BOOST_CHECK_EQUAL(idmap[3], 1);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 7);
    BOOST_CHECK_EQUAL(idmap[6], 2);
    BOOST_CHECK_EQUAL(idmap[7], 5);
    BOOST_CHECK_EQUAL(idmap[8], 8);

    orient = StdRegions::eDir1FwdDir2_Dir2BwdDir1;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 2);
    BOOST_CHECK_EQUAL(idmap[1], 5);
    BOOST_CHECK_EQUAL(idmap[2], 8);
    BOOST_CHECK_EQUAL(idmap[3], 1);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 7);
    BOOST_CHECK_EQUAL(idmap[6], 0);
    BOOST_CHECK_EQUAL(idmap[7], 3);
    BOOST_CHECK_EQUAL(idmap[8], 6);

    orient = StdRegions::eDir1BwdDir2_Dir2FwdDir1;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 6);
    BOOST_CHECK_EQUAL(idmap[1], 3);
    BOOST_CHECK_EQUAL(idmap[2], 0);
    BOOST_CHECK_EQUAL(idmap[3], 7);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 1);
    BOOST_CHECK_EQUAL(idmap[6], 8);
    BOOST_CHECK_EQUAL(idmap[7], 5);
    BOOST_CHECK_EQUAL(idmap[8], 2);

    orient = StdRegions::eDir1BwdDir2_Dir2BwdDir1;
    exp3d->ReOrientTracePhysMap(orient, idmap, nq0, nq1);
    BOOST_CHECK_EQUAL(idmap[0], 8);
    BOOST_CHECK_EQUAL(idmap[1], 5);
    BOOST_CHECK_EQUAL(idmap[2], 2);
    BOOST_CHECK_EQUAL(idmap[3], 7);
    BOOST_CHECK_EQUAL(idmap[4], 4);
    BOOST_CHECK_EQUAL(idmap[5], 1);
    BOOST_CHECK_EQUAL(idmap[6], 6);
    BOOST_CHECK_EQUAL(idmap[7], 3);
    BOOST_CHECK_EQUAL(idmap[8], 0);
}

} // namespace Expansion3DTests
} // namespace Nektar
