////////////////////////////////////////////////////////////////////////////////
//
//  File: Convert.h
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
//  Description: Mesh converter base class and XML writer.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_CONVERT
#define UTILITIES_PREPROCESSING_MESHCONVERT_CONVERT

#include <tinyxml/tinyxml.h>
#include <cstring>
#include <sstream>
#include <vector>
#include <list>

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "MeshElements.h"

namespace Nektar
{
    namespace Utilities
    {
        /// Base class for mesh conversions and XML writer.
        class Convert {
        public:
            /// Read a foreign-format file.
            virtual void ReadFile(const std::string pFilename) = 0;
            /// Process foreign-format file.
            virtual void Process() = 0;
            /// Write native Nektar++ XML formatted file.
            void WriteXmlFile(const std::string pFilename);

        protected:
            /// Dimension of the expansion.
            unsigned int                        m_expDim;
            /// Dimension of the space in which the mesh is defined.
            unsigned int                        m_spaceDim;
            /// List of mesh nodes.
            std::vector<NodeSharedPtr>          m_node;
            /// List of element vertices.
            std::vector<NodeSharedPtr>          m_vertex;
            /// List of element edges.
            std::vector<EdgeSharedPtr>          m_edge;
            /// List of element faces.
            std::vector<FaceSharedPtr>          m_face;
            /// List of elements.
            std::vector<ElementSharedPtr>       m_element;
            /// List of composites.
            std::vector<CompositeSharedPtr>     m_composite;

            /// Returns true if the given node exists in the mesh.
            int FindNodeIndex(const NodeSharedPtr pSrc);
            /// Returns true if the given node is a mesh vertex.
            int FindVertexIndex(const NodeSharedPtr pSrc);
            /// Returns the index of the edge in #m_edge, or -1 if not found.
            int FindEdgeIndex(const EdgeSharedPtr pSrc);
            /// Returns the index of the face in #m_face, or -1 if not found.
            int FindFaceIndex(const FaceSharedPtr pSrc);

        private:
            /// Writes the <NODES> section of the XML file.
            void WriteXmlNodes(TiXmlElement * pRoot);
            /// Writes the <EDGES> section of the XML file.
            void WriteXmlEdges(TiXmlElement * pRoot);
            /// Writes the <FACES> section of the XML file if needed.
            void WriteXmlFaces(TiXmlElement * pRoot);
            /// Writes the <ELEMENTS> section of the XML file.
            void WriteXmlElements(TiXmlElement * pRoot);
            /// Writes the <CURVES> section of the XML file if needed.
            void WriteXmlCurves(TiXmlElement * pRoot);
            /// Writes the <COMPOSITES> section of the XML file.
            void WriteXmlComposites(TiXmlElement * pRoot);
        };
    }
}
#endif
