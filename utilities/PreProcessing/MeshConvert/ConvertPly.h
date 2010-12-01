////////////////////////////////////////////////////////////////////////////////
//
//  File: ConvertPly.h
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
//  Description: PLY converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_CONVERTPLY
#define UTILITIES_PREPROCESSING_MESHCONVERT_CONVERTPLY

#include "Convert.h"

namespace Nektar
{
    namespace Utilities
    {
        /// Converter for Gmsh files.
        class ConvertPly : public Convert
        {
        public:
            /// Creates an instance of this class
            static boost::shared_ptr<Convert> create() {
                return MemoryManager<ConvertPly>::AllocateSharedPtr();
            }
            /// Name of class
            static std::string className;

            ConvertPly();
            ConvertPly(const ConvertPly& pSrc);
            virtual ~ConvertPly();

            /// Reads the Gmsh file.
            virtual void ReadFile(const std::string pFilename);
            /// Populate and validate required data structures.
            virtual void Process();

        private:
            /// Extract element vertices
            void ProcessVertices();
            /// Extract element edges
            void ProcessEdges();
            /// Extract element faces
            void ProcessFaces();
            /// Generate element IDs
            void ProcessElements();
            /// Generate composites
            void ProcessComposites();
            /// Returns the number of nodes in a given Gmsh entity.
            int GetNnodes(int GmshEntity);
        };
    }
}

#endif
