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

void TriangleInterface::Mesh(bool Quality)
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

    dt.in.numberofpoints          = numPoints;
    dt.in.numberofpointattributes = 0;
    dt.in.pointlist               = new double[dt.in.numberofpoints * 2];

    int pointc = 0;

    for (int i = 0; i < m_boundingloops.size(); i++)
    {
        for (int j = 0; j < m_boundingloops[i].size(); j++, pointc++)
        {
            nodemap[pointc] = m_boundingloops[i][j];

            Array<OneD, NekDouble> uv =
                m_boundingloops[i][j]->GetCADSurfInfo(sid);
            dt.in.pointlist[pointc * 2 + 0] = uv[0] * m_str;
            dt.in.pointlist[pointc * 2 + 1] = uv[1];
        }
    }

    for (int i = 0; i < m_stienerpoints.size(); i++, pointc++)
    {
        nodemap[pointc] = m_stienerpoints[i];

        Array<OneD, NekDouble> uv    = m_stienerpoints[i]->GetCADSurfInfo(sid);
        dt.in.pointlist[pointc * 2 + 0] = uv[0] * m_str;
        dt.in.pointlist[pointc * 2 + 1] = uv[1];
    }

    dt.in.numberofsegments = numSeg;
    dt.in.segmentlist      = new int[dt.in.numberofsegments * 2];
    pointc = 0;
    for (int i = 0; i < m_boundingloops.size(); i++, pointc++)
    {
        int first = pointc;
        for (int j = 0; j < m_boundingloops[i].size() - 1; j++, pointc++)
        {
            dt.in.segmentlist[pointc * 2 + 0] = pointc;
            dt.in.segmentlist[pointc * 2 + 1] = pointc + 1;
        }
        dt.in.segmentlist[pointc * 2 + 0] = pointc;
        dt.in.segmentlist[pointc * 2 + 1] = first;
    }

    dt.in.numberofregions = 0;
    dt.in.numberofholes   = m_centers.size() - 1;
    dt.in.holelist        = new double[dt.in.numberofholes * 2];

    for (int i = 1; i < m_centers.size(); i++)
    {
        dt.in.holelist[(i - 1) * 2 + 0] = m_centers[i][0] * m_str;
        dt.in.holelist[(i - 1) * 2 + 1] = m_centers[i][1];
    }

    string cmd;
    if (Quality)
    {
        cmd = "pqzQY";
    }
    else if (!Quality)
    {
        cmd = "pzQY";
    }
    char *cstr = new char[cmd.length() + 1];
    strcpy(cstr, cmd.c_str());

    dt.Run(cstr);
}

void TriangleInterface::SetUp()
{
    dt.in.pointlist                   = (double *)NULL;
    dt.in.pointattributelist          = (double *)NULL;
    dt.in.pointmarkerlist             = (int *)NULL;
    dt.in.numberofpoints              = 0;
    dt.in.numberofpointattributes     = 0;
    //
    dt.in.trianglelist                = (int *)NULL;
    dt.in.triangleattributelist       = (double *)NULL;
    dt.in.trianglearealist            = (double *)NULL;
    dt.in.neighborlist                = (int *)NULL;
    dt.in.numberoftriangles           = 0;
    dt.in.numberofcorners             = 0;
    dt.in.numberoftriangleattributes  = 0;
    //
    dt.in.segmentlist                 = (int *)NULL;
    dt.in.segmentmarkerlist           = (int *)NULL;
    dt.in.numberofsegments            = 0;
    //
    dt.in.holelist                    = (double *)NULL;
    dt.in.numberofholes               = 0;
    //
    dt.in.regionlist                  = (double *)NULL;
    dt.in.numberofregions             = 0;
    //
    dt.in.edgelist                    = (int *)NULL;
    dt.in.edgemarkerlist              = (int *)NULL;
    dt.in.normlist                    = (double *)NULL;
    dt.in.numberofedges               = 0;
    //
    dt.out.pointlist                  = (double *)NULL;
    dt.out.pointattributelist         = (double *)NULL;
    dt.out.pointmarkerlist            = (int *)NULL;
    dt.out.numberofpoints             = 0;
    dt.out.numberofpointattributes    = 0;
    //
    dt.out.trianglelist               = (int *)NULL;
    dt.out.triangleattributelist      = (double *)NULL;
    dt.out.trianglearealist           = (double *)NULL;
    dt.out.neighborlist               = (int *)NULL;
    dt.out.numberoftriangles          = 0;
    dt.out.numberofcorners            = 0;
    dt.out.numberoftriangleattributes = 0;
    //
    dt.out.segmentlist                = (int *)NULL;
    dt.out.segmentmarkerlist          = (int *)NULL;
    dt.out.numberofsegments           = 0;
    //
    dt.out.holelist                   = (double *)NULL;
    dt.out.numberofholes              = 0;
    //
    dt.out.regionlist                 = (double *)NULL;
    dt.out.numberofregions            = 0;
    //
    dt.out.edgelist                   = (int *)NULL;
    dt.out.edgemarkerlist             = (int *)NULL;
    dt.out.normlist                   = (double *)NULL;
    dt.out.numberofedges              = 0;
}

void TriangleInterface::Extract(
    std::vector<std::vector<NodeSharedPtr> > &Connec)
{
    Connec.clear();
    for (int i = 0; i < dt.out.numberoftriangles; i++)
    {
        map<int, NodeSharedPtr>::iterator n1, n2, n3;
        n1 = nodemap.find(dt.out.trianglelist[i * 3 + 0]);
        n2 = nodemap.find(dt.out.trianglelist[i * 3 + 1]);
        n3 = nodemap.find(dt.out.trianglelist[i * 3 + 2]);

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
