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

#include <LocalRegions/Expansion3D.h>
#include <StdRegions/StdRegions.hpp>
#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
namespace Expansion3DTests
{

BOOST_AUTO_TEST_CASE(TestReOrientQuadFacePhysMap)
{
    // Dummy 3D expansion
    LocalRegions::Expansion3DSharedPtr exp3d;

    // Initialiaze 3*3 face id map for quad
    int nq0 = 3;
    int nq1 = 3;
    int nvert = 4;
    Array<OneD,int> idmap(nq0*nq1, -1);

    // Test different orientations
    StdRegions::Orientation orient;

    orient = StdRegions::eDir1FwdDir1_Dir2FwdDir2;
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
    exp3d->ReOrientFacePhysMap(nvert, orient, nq0, nq1, idmap);
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
