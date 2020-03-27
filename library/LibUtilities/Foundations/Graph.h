////////////////////////////////////////////////////////////////////////////////
//
//  File:  Graph.h
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
//  Description:  
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef GRAPH_H
#define GRAPH_H

#include <list>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Foundations/FoundationsFwd.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        class GraphVertexObject
        {
        public:

            GraphVertexObject()
            {
                m_id = m_nextid++;
            }
            //    virtual int ReadFromFile(FILE * fp) = 0;
            //    virtual int WriteToFile(FILE * fp) = 0;

            GraphVertexObject(const GraphVertexID id)
            {
                m_id = id;
            }

            inline int  getid()
            {
                return m_id;
            }

            inline void setid(const GraphVertexID id)
            { 
                m_id = id;
            }

            LIB_UTILITIES_EXPORT virtual ~GraphVertexObject();

            LIB_UTILITIES_EXPORT friend bool operator  == (const GraphVertexObject &x, 
                const GraphVertexObject &y);
            LIB_UTILITIES_EXPORT friend bool operator  != (const GraphVertexObject &x, 
                const GraphVertexObject &y);

        protected:
            GraphVertexID m_id; //vertex object identifier

        private:
            LIB_UTILITIES_EXPORT static GraphVertexID m_nextid;
        }; 

        // ------------------------------------------------------------------------

        class GraphEdgeObject
        {
        public:

        protected:
            int m_gvoid1;
            int m_gvoid2;  //two graph vertex object identifiers which are being connected
        };

        // --------------------------------------------------------------------------

        class Graph
        {

        public:
            LIB_UTILITIES_EXPORT Graph();
            LIB_UTILITIES_EXPORT ~Graph();

        protected:
            int m_curmaxvid;
            std::list<GraphVertexObject*> m_vertset;
            std::list<GraphEdgeObject*>   m_edgeset;

        };

    }  //end of namespace
}  //end of namespace
#endif // GRAPH_H

