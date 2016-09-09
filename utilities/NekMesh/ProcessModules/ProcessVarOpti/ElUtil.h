////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI_ELUTIL
#define UTILITIES_NEKMESH_PROCESSVAROPTI_ELUTIL

#include "../../Module.h"

#include <LibUtilities/BasicUtils/Thread.h>

namespace Nektar
{
namespace Utilities
{

struct DerivUtil;
struct Residual;

typedef boost::shared_ptr<DerivUtil> DerivUtilSharedPtr;
typedef boost::shared_ptr<Residual> ResidualSharedPtr;

class ElUtilJob;

class ElUtil : public boost::enable_shared_from_this<ElUtil>
{
public:
    ElUtil(ElementSharedPtr e, DerivUtilSharedPtr d,
           ResidualSharedPtr, int n);

    ElUtilJob *GetJob();

    int GetId()
    {
        return m_el->GetId();
    }

    std::vector<std::vector<NekDouble *> > nodes;
    std::vector<Array<OneD, NekDouble> > maps;
    NekDouble minJac;
    NekDouble scaledJac;

    void Evaluate();

    ElementSharedPtr GetEl()
    {
        return m_el;
    }

private:

    std::vector<Array<OneD, NekDouble> > MappingIdealToRef();

    ElementSharedPtr m_el;
    int m_dim;
    int m_mode;
    std::map<int,int> m_idmap;

    DerivUtilSharedPtr derivUtil;
    ResidualSharedPtr res;
};
typedef boost::shared_ptr<ElUtil> ElUtilSharedPtr;

class ElUtilJob : public Thread::ThreadJob
{
public:
    ElUtilJob(ElUtil* e) : el(e) {}

    void Run()
    {
        el->Evaluate();
    }
private:
    ElUtil* el;
};

}
}

#endif
