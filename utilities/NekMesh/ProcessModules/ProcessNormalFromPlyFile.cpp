///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessNormalFromPlyFile.cpp
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
//  part of the ProcessSpherigon method
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessSpherigon.h"
#include <LibUtilities/BasicUtils/Progressbar.hpp>

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
    Node     tmp,tmpsav;
    NekDouble mindiff, diff;

    for (vIt = surfverts.begin(); vIt != surfverts.end(); ++vIt, ++cnt)
    {
        if(m_mesh->m_verbose)
        {
            prog = LibUtilities::PrintProgressbar(cnt,surfverts.size(),
                                                  "Nearest ply verts",prog);
        }
        
        mindiff = 1e12;
        
        for (j = 0, it = plymesh->m_vertexSet.begin();
             it != plymesh->m_vertexSet.end();
             ++it, ++j)
        {
            tmp  = *(vIt->second) - *(*it);
            diff = tmp.abs2();
            
            if (diff < mindiff)
            {
                mindiff = diff;
                cntmin  = (*it)->m_id;
                tmpsav  = tmp;
            }
        }
        
        ASSERTL1(cntmin < plymesh->m_vertexNormals.size(),
                 "cntmin is out of range");
        m_mesh->m_vertexNormals[vIt->first] =
            plymesh->m_vertexNormals[cntmin];
        
    }
}
}
}
