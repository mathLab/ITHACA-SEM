////////////////////////////////////////////////////////////////////////////////
//
//  File: TriangleInterface.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Interface to triangle mesher
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/ExtLibInterface/TriangleInterface.h>

#include <sstream>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void TriangleInterface::Mesh(bool Quiet, bool Quality)
{
    SetUp();

    int numPoints = 0;
    int numSeg = 0;
    for (int i = 0; i < m_boundingloops.size(); i++)
    {
        numSeg += m_boundingloops[i].size();
    }
    numPoints = numSeg + m_stienerpoints.size();

    stringstream ss;
    ss << "3 points required for triangulation, " << numPoints << " provided";

    ASSERTL0(numPoints > 2, ss.str());

    in.numberofpoints          = numPoints;
    in.numberofpointattributes = 0;
    in.pointlist               = new REAL[in.numberofpoints * 2];
    in.pointmarkerlist         = new int [in.numberofpoints];

    int pointc = 0;

    for (int i = 0; i < m_boundingloops.size(); i++)
    {
        for (int j = 0; j < m_boundingloops[i].size(); j++, pointc++)
        {
            nodemap[pointc] = m_boundingloops[i][j];

            Array<OneD, NekDouble> uv =
                m_boundingloops[i][j]->GetCADSurfInfo(sid);
            in.pointlist[pointc * 2 + 0] = uv[0] * m_str;
            in.pointlist[pointc * 2 + 1] = uv[1];
        }
    }

    for (int i = 0; i < m_stienerpoints.size(); i++, pointc++)
    {
        nodemap[pointc] = m_stienerpoints[i];

        Array<OneD, NekDouble> uv    = m_stienerpoints[i]->GetCADSurfInfo(sid);
        in.pointlist[pointc * 2 + 0] = uv[0] * m_str;
        in.pointlist[pointc * 2 + 1] = uv[1];
    }

    in.numberofsegments = numSeg;
    in.segmentlist      = new int[in.numberofsegments * 2];
    pointc = 0;
    for (int i = 0; i < m_boundingloops.size(); i++, pointc++)
    {
        int first = pointc;
        for (int j = 0; j < m_boundingloops[i].size() - 1; j++, pointc++)
        {
            in.segmentlist[pointc * 2 + 0] = pointc;
            in.segmentlist[pointc * 2 + 1] = pointc + 1;
        }
        in.segmentlist[pointc * 2 + 0] = pointc;
        in.segmentlist[pointc * 2 + 1] = first;
    }

    in.numberofregions = 0;
    in.numberofholes   = m_centers.size() - 1;
    in.holelist        = new REAL[in.numberofholes * 2];

    for (int i = 1; i < m_centers.size(); i++)
    {
        in.holelist[(i - 1) * 2 + 0] = m_centers[i][0] * m_str;
        in.holelist[(i - 1) * 2 + 1] = m_centers[i][1];
    }

    if (Quiet && Quality)
    {
        triangulate("pzenqQYY", &in, &out, NULL);
    }
    else if (Quiet && !Quality)
    {
        triangulate("pzenYYQ", &in, &out, NULL);
    }
    else if (!Quiet && Quality)
    {
        triangulate("pzenqYY", &in, &out, NULL);
    }
    else if (!Quiet && !Quality)
    {
        triangulate("pzenYY", &in, &out, NULL);
    }

    // verify the mesh a bit
    if (out.numberofpoints - out.numberofedges + out.numberoftriangles !=
        2 - m_centers.size())
    {
        cout << endl << "epc wrong" << endl;
        cout << out.numberofpoints - out.numberofedges + out.numberoftriangles
             << " " << m_centers.size() << " " << sid << endl;
    }
}

void TriangleInterface::SetUp()
{
    in.pointlist                   = (REAL *)NULL;
    in.pointattributelist          = (REAL *)NULL;
    in.pointmarkerlist             = (int *)NULL;
    in.numberofpoints              = 0;
    in.numberofpointattributes     = 0;
    //
    in.trianglelist                = (int *)NULL;
    in.triangleattributelist       = (REAL *)NULL;
    in.trianglearealist            = (REAL *)NULL;
    in.neighborlist                = (int *)NULL;
    in.numberoftriangles           = 0;
    in.numberofcorners             = 0;
    in.numberoftriangleattributes  = 0;
    //
    in.segmentlist                 = (int *)NULL;
    in.segmentmarkerlist           = (int *)NULL;
    in.numberofsegments            = 0;
    //
    in.holelist                    = (REAL *)NULL;
    in.numberofholes               = 0;
    //
    in.regionlist                  = (REAL *)NULL;
    in.numberofregions             = 0;
    //
    in.edgelist                    = (int *)NULL;
    in.edgemarkerlist              = (int *)NULL;
    in.normlist                    = (REAL *)NULL;
    in.numberofedges               = 0;
    //
    out.pointlist                  = (REAL *)NULL;
    out.pointattributelist         = (REAL *)NULL;
    out.pointmarkerlist            = (int *)NULL;
    out.numberofpoints             = 0;
    out.numberofpointattributes    = 0;
    //
    out.trianglelist               = (int *)NULL;
    out.triangleattributelist      = (REAL *)NULL;
    out.trianglearealist           = (REAL *)NULL;
    out.neighborlist               = (int *)NULL;
    out.numberoftriangles          = 0;
    out.numberofcorners            = 0;
    out.numberoftriangleattributes = 0;
    //
    out.segmentlist                = (int *)NULL;
    out.segmentmarkerlist          = (int *)NULL;
    out.numberofsegments           = 0;
    //
    out.holelist                   = (REAL *)NULL;
    out.numberofholes              = 0;
    //
    out.regionlist                 = (REAL *)NULL;
    out.numberofregions            = 0;
    //
    out.edgelist                   = (int *)NULL;
    out.edgemarkerlist             = (int *)NULL;
    out.normlist                   = (REAL *)NULL;
    out.numberofedges              = 0;
}

void TriangleInterface::Extract(
    std::vector<std::vector<NodeSharedPtr> > &Connec)
{
    Connec.clear();
    for (int i = 0; i < out.numberoftriangles; i++)
    {
        map<int, NodeSharedPtr>::iterator n1, n2, n3;
        n1 = nodemap.find(out.trianglelist[i * 3 + 0]);
        n2 = nodemap.find(out.trianglelist[i * 3 + 1]);
        n3 = nodemap.find(out.trianglelist[i * 3 + 2]);

        ASSERTL0(n1 != nodemap.end() && n2 != nodemap.end() &&
                     n3 != nodemap.end(),
                 "node index error");

        vector<NodeSharedPtr> tri(3);
        tri[0] = n1->second;
        tri[1] = n2->second;
        tri[2] = n3->second;
        Connec.push_back(tri);
    }
}
}
}
