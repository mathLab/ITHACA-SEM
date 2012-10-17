////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LibUtilities/Foundations/Graph.cpp,v $
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
//  Description:  
//
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Graph.h>

namespace Nektar
{
    namespace LibUtilities
    {

        int GraphVertexObject::m_nextid = 0;

        GraphVertexObject::~GraphVertexObject()
        {
        }

        Graph::Graph()
        {
        }

        Graph::~Graph()
        {
#if 0       
            std::list<GraphVertexObject *>::const_iterator def; 
            for(def = m_vertset.begin(); def != m_vertset.end(); ++def)
            {
                delete  *def;
            }

            std::list<GraphEdgeObject *>::const_iterator def1; 
            for(def1 = _Edgeset.begin(); def1 != _Edgeset.end(); ++def1)
            {
                delete  *def1;
            }
#endif    
        }

        bool operator  == (const GraphVertexObject &x, const GraphVertexObject &y)
        {
            return (x.m_id == y.m_id);
        }

        bool operator != (const GraphVertexObject &x, const GraphVertexObject &y)
        {
            return (x.m_id != y.m_id);
        }
    } // end of namespace 
} // end of namespace 

//
// $Log: Graph.cpp,v $
// Revision 1.3  2006/08/16 03:26:36  jfrazier
// Optimized the equality and inequality operators.
//
// Revision 1.2  2006/06/01 14:59:58  sherwin
// Modified header file path
//
// Revision 1.1  2006/06/01 14:56:20  kirby
// *** empty log message ***
//
// Revision 1.1  2006/05/04 18:57:42  kirby
// *** empty log message ***
//
// Revision 1.11  2006/03/25 00:52:43  jfrazier
// Minor formatting stuff to correct indenting.
//
// Revision 1.10  2006/03/13 11:17:03  sherwin
//
// First compiing version of Demos in SpatialDomains and LocalRegions. However they do not currently seem to execute properly
//
// Revision 1.9  2006/03/12 07:42:48  sherwin
//
// Updated to meet coding standard. Still has not been compiled
//
//



