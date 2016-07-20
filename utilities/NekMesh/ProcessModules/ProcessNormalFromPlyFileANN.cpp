///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessNormalFromPlyFileANN.cpp
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
//  Description: Method for processing normals from ply file which is
//  part of the ProcessSpherigon method using libANN
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessSpherigon.h"
#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include "ANN/ANN.h"

using namespace std;

namespace Nektar
{
namespace Utilities
{
void ProcessSpherigon::FindNormalFromPlyFile(MeshSharedPtr &plymesh,
                                             map<int,NodeSharedPtr> &surfverts)
{
    int      cnt = 0;
    int      j;
    int      prog,cntmin;
    NodeSet::iterator it;
    map<int, NodeSharedPtr>::iterator vIt;

    int n_neighbs = 5;
    int neighbs_max = 100; 
    int nplypts = plymesh->m_vertexSet.size();
    ANNpointArray dataPts = annAllocPts(nplypts, 3);
    ANNidxArray   nnIdx = new ANNidx [neighbs_max];
    ANNdistArray  dists = new ANNdist[neighbs_max];
    map<int,int>  AnnidtoPlyid; 
        
    //Fill vertex array into libAnn format
    for (j = 0, it = plymesh->m_vertexSet.begin();
         it != plymesh->m_vertexSet.end();
         ++it, ++j)
    {
        dataPts[j][0] = (*it)->m_x;
        dataPts[j][1] = (*it)->m_y;
        dataPts[j][2] = (*it)->m_z;
        AnnidtoPlyid[j] = (*it)->m_id;
    }
    
    //Build tree
    ANNkd_tree* kdTree;
    // build search structure
    kdTree = new ANNkd_tree( dataPts,	  // the data points
                             nplypts,	  // number of points
                             3);	  // dimension of space
    
    //Find neipghbours
    ANNpoint queryPt = annAllocPt(3);
    
    for (cnt = 0, vIt = surfverts.begin(); vIt != surfverts.end();
         ++vIt, ++cnt)
    {
        if(m_mesh->m_verbose)
        {
            prog = LibUtilities::PrintProgressbar(cnt,surfverts.size(),
                                                  "Nearest ply verts",prog);
        }
        
        n_neighbs  = 5; 
        queryPt[0] = vIt->second->m_x; 
        queryPt[1] = vIt->second->m_y;
        queryPt[2] = vIt->second->m_z;
        kdTree->annkSearch(queryPt, n_neighbs, nnIdx, dists, 0); //eps set to zero
        
        ASSERTL1(dists[0] < dists[1],"Assumption that dist values are ordered from smallest to largest is not correct");
        
        cntmin = AnnidtoPlyid[nnIdx[0]];
        
        ASSERTL1(cntmin < plymesh->m_vertexNormals.size(),
                 "cntmin is out of range");
        
        m_mesh->m_vertexNormals[vIt->first] =
            plymesh->m_vertexNormals[cntmin];
    }
}
}
}
