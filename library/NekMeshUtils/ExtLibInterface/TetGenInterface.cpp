////////////////////////////////////////////////////////////////////////////////
//
//  File: TetGenInterface.cpp
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
//  Description: tetgen interface methods
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/ExtLibInterface/TetGenInterface.h>

#define TETLIBRARY
#include <tetgen.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

void TetGenInterface::InitialMesh(map<int, NodeSharedPtr>  tgidton,
                                  vector<Array<OneD, int> > tri)
{
    surface.initialize();
    output.initialize();

    // build surface input
    tetgenio::facet *f;
    tetgenio::polygon *p;

    surface.firstnumber    = 0;
    surface.numberofpoints = tgidton.size();
    surface.pointlist      = new REAL[surface.numberofpoints * 3];

    map<int, NodeSharedPtr>::iterator it;
    for (it = tgidton.begin(); it != tgidton.end(); it++)
    {
        Array<OneD, NekDouble> loc = it->second->GetLoc();

        surface.pointlist[it->first * 3 + 0] = loc[0];
        surface.pointlist[it->first * 3 + 1] = loc[1];
        surface.pointlist[it->first * 3 + 2] = loc[2];
    }

    surface.numberoffacets  = tri.size();
    surface.facetlist       = new tetgenio::facet[surface.numberoffacets];
    surface.facetmarkerlist = new int[surface.numberoffacets];

    for (int i = 0; i < tri.size(); i++)
    {
        f                          = &surface.facetlist[i];
        f->numberofpolygons        = 1;
        f->polygonlist             = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes           = 0;
        f->holelist                = NULL;
        p                          = &f->polygonlist[0];
        p->numberofvertices        = 3;
        p->vertexlist              = new int[p->numberofvertices];
        p->vertexlist[0]           = tri[i][0];
        p->vertexlist[1]           = tri[i][1];
        p->vertexlist[2]           = tri[i][2];
        surface.facetmarkerlist[i] = 1;
    }

    surface.numberofholes = m_holes.size();

    if (m_holes.size() > 0)
    {
        surface.holelist = new REAL[m_holes.size() * 3];
        for (int j = 0; j < m_holes.size(); ++j)
        {
            surface.holelist[3 * j + 0] = m_holes[j][0];
            surface.holelist[3 * j + 1] = m_holes[j][1];
            surface.holelist[3 * j + 2] = m_holes[j][2];
        }
    }
    else
    {
        surface.holelist = NULL;
    }

    string cmd = "pYzqQ";
    char *cstr = new char[cmd.length() + 1];
    strcpy(cstr, cmd.c_str());

    tetrahedralize(cstr, &surface, &output);
}

void TetGenInterface::GetNewPoints(int num,
                                   vector<Array<OneD, NekDouble> > &newp)
{
    for (int i = num; i < output.numberofpoints; i++)
    {
        Array<OneD, NekDouble> loc(3);
        loc[0] = output.pointlist[i * 3 + 0];
        loc[1] = output.pointlist[i * 3 + 1];
        loc[2] = output.pointlist[i * 3 + 2];
        newp.push_back(loc);
    }
}

void TetGenInterface::RefineMesh(std::map<int, NekDouble> delta)
{
    input = output;

    input.numberofpointmtrs = 1;

    input.pointmtrlist = new REAL[input.numberofpoints];

    for (int i = 0; i < input.numberofpoints; i++)
    {
        input.pointmtrlist[i] = delta[i];
    }

    string cmd = "pYrmzq1.1/0QO2/7";
    char *cstr = new char[cmd.length() + 1];
    strcpy(cstr, cmd.c_str());

    tetrahedralize(cstr, &input, &output);
}

vector<Array<OneD, int> > TetGenInterface::Extract()
{
    vector<Array<OneD, int> > tets;

    for (int i = 0; i < output.numberoftetrahedra; i++)
    {
        Array<OneD, int> tet(4);
        tet[0] = output.tetrahedronlist[i * 4 + 0];
        tet[1] = output.tetrahedronlist[i * 4 + 1];
        tet[2] = output.tetrahedronlist[i * 4 + 2];
        tet[3] = output.tetrahedronlist[i * 4 + 3];
        tets.push_back(tet);
    }

    return tets;
}

void TetGenInterface::freetet()
{
    surface.deinitialize();
    input.deinitialize();
    output.deinitialize();
}
}
}
