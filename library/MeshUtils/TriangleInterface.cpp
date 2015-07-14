////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
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

#include <string>
#include <fstream>

#include <MeshUtils/TriangleInterface.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

    void TriangleInterface::Mesh(bool Quiet, bool Quality)
    {

        if(meshloaded)
        {
            freetri();
        }
        ASSERTL0(meshloaded==false,"Mesh must be cleared before meshing");

        int numPoints = 0;
        int numSeg = 0;
        for(int i = 0; i < m_boundingloops.size(); i++)
        {
            numSeg+=m_boundingloops[i].size();
        }
        numPoints = numSeg + m_stienerpoints.size();

        in.numberofpoints = numPoints;
        in.numberofpointattributes = 0;
        in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
        in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));

        int pointc = 0;

        for(int i = 0; i < m_boundingloops.size(); i++)
        {
            for(int j = 0; j < m_boundingloops[i].size(); j++)
            {
                nodemap[pointc] = m_boundingloops[i][j];
                nodemapr[m_boundingloops[i][j]] = pointc;

                Array<OneD, NekDouble> uv = Nodes[m_boundingloops[i][j]]
                                                ->GetS(sid);

                in.pointlist[pointc*2+0] = uv[0]*m_str;
                in.pointlist[pointc*2+1] = uv[1];
                pointc++;
            }
        }

        for(int i = 0; i < m_stienerpoints.size(); i++)
        {
            nodemap[pointc] = m_stienerpoints[i];
            nodemapr[m_stienerpoints[i]] = pointc;
            Array<OneD, NekDouble> uv = Nodes[m_stienerpoints[i]]
                                                ->GetS(sid);
            in.pointlist[pointc*2+0] = uv[0]*m_str;
            in.pointlist[pointc*2+1] = uv[1];
            pointc++;
        }

        in.numberofsegments = numSeg;
        in.segmentlist = (int *) malloc(in.numberofsegments*2*sizeof(int));
        pointc=0;
        for(int i = 0; i < m_boundingloops.size(); i++)
        {
            for(int j = 0; j < m_boundingloops[i].size()-1; j++)
            {
                in.segmentlist[pointc*2+0] = nodemapr[m_boundingloops[i][j]];
                in.segmentlist[pointc*2+1] = nodemapr[m_boundingloops[i][j+1]];
                pointc++;
            }
            in.segmentlist[pointc*2+0] = nodemapr[m_boundingloops[i].back()];
            in.segmentlist[pointc*2+1] = nodemapr[m_boundingloops[i][0]];
            pointc++;
        }

        in.numberofregions = 0;
        in.numberofholes = m_centers.size()-1;
        in.holelist = (REAL *) malloc(in.numberofholes*2*sizeof(REAL));

        for(int i = 1; i < m_centers.size(); i++)
        {
            in.holelist[(i-1)*2+0] = m_centers[i][0]*m_str;
            in.holelist[(i-1)*2+1] = m_centers[i][1];
        }

        out.pointlist = (REAL *) NULL;
        out.pointattributelist = (REAL *) NULL;
        out.pointmarkerlist = (int *) NULL;
        out.trianglelist = (int *) NULL;
        out.trianglearealist = (REAL *) NULL;
        out.triangleattributelist = (REAL *) NULL;
        out.neighborlist = (int *) NULL;
        out.segmentlist = (int *) NULL;
        out.segmentmarkerlist = (int *) NULL;
        out.edgelist = (int *) NULL;
        out.edgemarkerlist = (int *) NULL;

        if(Quiet && Quality)
        {
            triangulate("pzenqQYY", &in, &out,  NULL);
        }
        else if(Quiet && !Quality)
        {
            triangulate("pzenQYY", &in, &out,  NULL);
        }
        else if(!Quiet && Quality)
        {
            triangulate("pzenqYY", &in, &out,  NULL);
        }
        else if(!Quiet && !Quality)
        {
            triangulate("pzenYY", &in, &out,  NULL);
        }

        for(int i = 0; i < out.numberofpoints; i++)
        {
            out.pointlist[i*2+0] = out.pointlist[i*2+0]/m_str;
        }

    }

    void TriangleInterface::Extract(int &np, int &nt,
                                    Array<OneD, Array<OneD, int> > &Connec)
    {
        Connec = Array<OneD, Array<OneD, int> >(out.numberoftriangles);
        nt = out.numberoftriangles;
        np = out.numberofpoints;
        for(int i = 0; i < out.numberoftriangles; i++)
        {
            Array<OneD, int> tri(3);
            tri[0] = nodemap[out.trianglelist[i*3+0]];
            tri[1] = nodemap[out.trianglelist[i*3+1]];
            tri[2] = nodemap[out.trianglelist[i*3+2]];
            Connec[i] = tri;
        }
    }

    void TriangleInterface::GetEdges(Array<OneD, Array<OneD, int> > &edgelist,
                                     int &num)
    {
        edgelist = Array<OneD, Array<OneD, int> >(out.numberofedges);
        for(int i = 0; i < out.numberofedges; i++)
        {
            Array<OneD, int> el(2);
            el[0] = nodemap[out.edgelist[i*2+0]];
            el[1] = nodemap[out.edgelist[i*2+1]];
            edgelist[i] = el;
        }
        num = out.numberofedges;
    }
    void TriangleInterface::GetNeighbour(
                            Array<OneD, Array<OneD, int> > &neigbourlist)
    {
        neigbourlist = Array<OneD, Array<OneD, int> >(out.numberoftriangles);
        for(int i = 0; i < out.numberoftriangles; i++)
        {
            Array<OneD, int> nl(3);
            nl[0] = out.neighborlist[i*3+0];
            nl[1] = out.neighborlist[i*3+1];
            nl[2] = out.neighborlist[i*3+2];
            neigbourlist[i]=nl;
        }
    }

    void TriangleInterface::freetri()
    {
        if(meshloaded)
        {
            free(in.pointlist);
            free(in.pointmarkerlist);
            free(in.segmentlist);
            free(in.holelist);

            free(out.pointlist);
            free(out.pointattributelist);
            free(out.pointmarkerlist);
            free(out.trianglelist);
            free(out.triangleattributelist);
            free(out.trianglearealist);
            free(out.neighborlist);
            free(out.segmentlist);
            free(out.segmentmarkerlist);
            free(out.edgelist);
            free(out.edgemarkerlist);
        }
        meshloaded = false;
    }

}
}
