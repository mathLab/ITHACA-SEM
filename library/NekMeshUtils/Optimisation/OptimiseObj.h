////////////////////////////////////////////////////////////////////////////////
//
//  File: Curvemesh.h
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
//  Description: object for individual curve meshes.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_OPTIMISATION_OPTIMISEOBJ_H
#define NEKTAR_MESHUTILS_OPTIMISATION_OPTIMISEOBJ_H

namespace Nektar
{
namespace NekMeshUtils
{

class OptiObj
{
    public:

        OptiObj(){};

        virtual ~OptiObj(){};

        virtual NekDouble F(Array<OneD, NekDouble> xitst)
        {
            ASSERTL0(false,"should be implemented in inheriting class");
            return 0.0;
        };

        virtual DNekMat dF(Array<OneD, NekDouble> xitst)
        {
            ASSERTL0(false,"should be implemented in inheriting class");
            return DNekMat(1,1,0.0);
        };

        virtual Array<OneD, NekDouble> Getxi()
        {
            ASSERTL0(false,"should be implemented in inheriting class");
            return Array<OneD,NekDouble>();
        };

        virtual Array<OneD, NekDouble> Getli()
        {
            ASSERTL0(false,"should be implemented in inheriting class");
            return Array<OneD,NekDouble>();
        };

        virtual Array<OneD, NekDouble> Getui()
        {
            ASSERTL0(false,"should be implemented in inheriting class");
            return Array<OneD,NekDouble>();
        };

        virtual void Update(Array<OneD, NekDouble> xinew)
        {
            ASSERTL0(false,"should be implemented in inheriting class");
        };

};
typedef boost::shared_ptr<OptiObj> OptiObjSharedPtr;

}
}
#endif
