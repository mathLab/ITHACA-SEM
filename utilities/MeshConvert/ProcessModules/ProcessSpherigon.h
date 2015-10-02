////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessSpherigon.h
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
//  Description: Refine boundary layer of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_PROCESSJAC
#define UTILITIES_PREPROCESSING_MESHCONVERT_PROCESSJAC


#include "../Module.h"
#include "../InputModules/InputPly.h"

namespace Nektar
{
    namespace Utilities
    {
        class ProcessSpherigon : public ProcessModule
        {
        public:
            /// Creates an instance of this class
            static boost::shared_ptr<Module> create(MeshSharedPtr m) 
            {
                return MemoryManager<ProcessSpherigon>::AllocateSharedPtr(m);
            }
            static ModuleKey className;
            
            ProcessSpherigon(MeshSharedPtr m);
            virtual ~ProcessSpherigon();
            
            /// Write mesh to output file.
            virtual void Process();
            
        protected:
            void   GenerateNormals(vector<ElementSharedPtr> &el,
                                   MeshSharedPtr &mesh);
            double CrossProdMag   (Node &a, Node &b);
            void   UnitCrossProd  (Node &a, Node &b, Node &c);
            double Blend          (double r);
            void   SuperBlend     (vector<double> &r, 
                                   vector<Node>   &Q, 
                                   Node           &P, 
                                   vector<double> &blend);
        };
    }
}

#endif
