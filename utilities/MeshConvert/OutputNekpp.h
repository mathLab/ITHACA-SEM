////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputNekpp.h
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_OUTPUTNEKPP
#define UTILITIES_PREPROCESSING_MESHCONVERT_OUTPUTNEKPP

#include <tinyxml.h>
#include "Module.h"

namespace Nektar
{
    namespace Utilities
    {
        /// Converter for Gmsh files.
        class OutputNekpp : public OutputModule
        {
        public:
            /// Creates an instance of this class
            static boost::shared_ptr<Module> create(MeshSharedPtr m) {
                return MemoryManager<OutputNekpp>::AllocateSharedPtr(m);
            }
            static ModuleKey className;
            
            OutputNekpp(MeshSharedPtr m);
            virtual ~OutputNekpp();
            
            /// Write mesh to output file.
            virtual void Process();

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
            /// Writes the <DOMAIN> section of the XML file.
            void WriteXmlDomain(TiXmlElement * pRoot);
            /// Writes the <EXPANSIONS> section of the XML file.
            void WriteXmlExpansions(TiXmlElement * pRoot);
            /// Writes the <CONDITIONS> section of the XML file.
            void WriteXmlConditions(TiXmlElement * pRoot);
        };
    }
}

#endif
