////////////////////////////////////////////////////////////////////////////////
//
//  File: InputGmsh.h
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

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_INPUTGMSH
#define UTILITIES_PREPROCESSING_MESHCONVERT_INPUTGMSH

#include "../Module.h"

namespace Nektar
{
    namespace Utilities
    {
        /**
         * Converter for Gmsh files.
         */
        class InputGmsh : public InputModule
        {
        public:
            InputGmsh(MeshSharedPtr m);
            virtual ~InputGmsh();
            virtual void Process();

            /// Creates an instance of this class
            static ModuleSharedPtr create(MeshSharedPtr m) {
                return MemoryManager<InputGmsh>::AllocateSharedPtr(m);
            }
            /// %ModuleKey for class.
            static ModuleKey className;
            static std::map<unsigned int, ElmtConfig> GenElmMap();
            
            /**
             * Element map; takes a msh id to an %ElmtConfig object.
             */
            static std::map<unsigned int, ElmtConfig> elmMap;

        private:
            int GetNnodes(unsigned int InputGmshEntity);
            vector<int> CreateReordering(unsigned int InputGmshEntity);
            vector<int> TriReordering(ElmtConfig conf);
            vector<int> QuadReordering(ElmtConfig conf);
            vector<int> HexReordering(ElmtConfig conf);
            vector<int> PrismReordering(ElmtConfig conf);
            vector<int> TetReordering(ElmtConfig conf);
        };
    }
}

#endif
