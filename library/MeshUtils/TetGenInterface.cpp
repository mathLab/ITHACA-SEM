////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <MeshUtils/TetGenInterface.h>

using namespace std;
namespace Nektar{
namespace MeshUtils{

void TetGenInterface::Mesh(bool Quiet, bool Quality)
{
    if(meshloaded)
    {
        freetet();
    }
    surface.initialize();
    additional.initialize();
    output.initialize();

    nodemap.clear(); nodemapr.clear();

    //build surface input
    tetgenio::facet *f;
    tetgenio::polygon *p;

    surface.firstnumber = 0;
    surface.numberofpoints = m_nodesinsurface.size();
    surface.pointlist = new REAL[surface.numberofpoints*3];

    int pointc = 0;

    for(int i = 0; i < m_nodesinsurface.size(); i++)
    {
        nodemap[pointc] = m_nodesinsurface[i];
        nodemapr[m_nodesinsurface[i]] = pointc;

        Array<OneD, NekDouble> loc = Nodes[m_nodesinsurface[i]]->GetLoc();

        surface.pointlist[i*3+0] = loc[0];
        surface.pointlist[i*3+1] = loc[1];
        surface.pointlist[i*3+2] = loc[2];

        pointc++;
    }

    //build stiener list in additional
    additional.firstnumber = 0;
    additional.numberofpoints = m_stienerpoints.size();
    additional.pointlist = new REAL[additional.numberofpoints*3];

    for(int i = 0; i < m_stienerpoints.size(); i++)
    {
        nodemap[pointc] = m_stienerpoints[i];
        nodemapr[m_stienerpoints[i]] = pointc;

        Array<OneD, NekDouble> loc = Nodes[m_stienerpoints[i]]->GetLoc();
        additional.pointlist[i*3+0] = loc[0];
        additional.pointlist[i*3+1] = loc[1];
        additional.pointlist[i*3+2] = loc[2];

        pointc++;
    }

    surface.numberoffacets = Tris.size();
    surface.facetlist = new tetgenio::facet[surface.numberoffacets];
    surface.facetmarkerlist = new int[surface.numberoffacets];

    for(int i = 0; i < Tris.size(); i++)
    {
        f = &surface.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];

        Array<OneD, int> n = Tris[i]->GetN();
        p->vertexlist[0] = nodemapr[n[0]];
        p->vertexlist[1] = nodemapr[n[1]];
        p->vertexlist[2] = nodemapr[n[2]];
        surface.facetmarkerlist[i] = 0;

    }

    if(Quiet)
    {
        tetrahedralize("pYizfennQ", &surface, &output, &additional, NULL);
    }
    else
    {
        tetrahedralize("pYizfenn", &surface, &output, &additional, NULL);
    }

}

void TetGenInterface::Extract(int &numtet,
                    Array<OneD, Array<OneD, int> > &tetconnect)
{
    numtet = output.numberoftetrahedra;
    tetconnect = Array<OneD, Array<OneD, int> >(numtet);

    for(int i = 0; i < numtet; i++)
    {
        Array<OneD, int> tet(4);
        tet[0] = nodemap[output.tetrahedronlist[i*4+0]];
        tet[1] = nodemap[output.tetrahedronlist[i*4+1]];
        tet[2] = nodemap[output.tetrahedronlist[i*4+2]];
        tet[3] = nodemap[output.tetrahedronlist[i*4+3]];
        tetconnect[i] = tet;
    }
}


void TetGenInterface::freetet()
{
    if(meshloaded)
    {
        surface.deinitialize();
        additional.deinitialize();
        output.deinitialize();
    }
}

}
}
