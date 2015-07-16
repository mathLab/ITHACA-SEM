////////////////////////////////////////////////////////////////////////////////
//
//  File: InputStarTec.h
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
//  Description: Tecplot (ascii .dat) converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_INPUTTEC
#define UTILITIES_PREPROCESSING_MESHCONVERT_INPUTTEC

#include "../Module.h"
#include "../MeshElements.h"

namespace Nektar
{
    namespace Utilities
    {
        /// Converter for VTK files.
        class InputTec : public InputModule
        {
        public:
            /// Creates an instance of this class
            static ModuleSharedPtr create(MeshSharedPtr m) {
                return MemoryManager<InputTec>::AllocateSharedPtr(m);
            }
            static ModuleKey className;

            InputTec(MeshSharedPtr m);
            virtual ~InputTec();

            /// Populate and validate required data structures.
            virtual void Process();

            void ReadZone(int &nComposite);

        protected:
            
            void GenElement3D(vector<NodeSharedPtr> &Nodes,
                              int i, vector<int> &ElementFaces, 
                              vector<vector<int> >&FaceNodes,
                              int ncomposite,
                              bool DoOrient);


            void GenElement2D(vector<NodeSharedPtr> &Nodes,
                              int i, vector<int> &ElementFaces, 
                              vector<vector<int> >&FaceNodes,
                              int ncomposite); 


            Array<OneD, int> SortEdgeNodes(vector<NodeSharedPtr> &Nodes,
                                           vector<int> &ElementFaces, 
                                           vector<vector<int> >&FaceNodes);

            Array<OneD, int> SortFaceNodes(vector<NodeSharedPtr> &Nodes,
                                           vector<int> &ElementFaces, 
                                           vector<vector<int> >&FaceNodes);

            void ResetNodes(vector<NodeSharedPtr> &Nodes,
                            Array<OneD, vector<int> >&ElementFaces,
                            vector<vector<int> >&FaceNodes);
        };
    }
}

#endif
