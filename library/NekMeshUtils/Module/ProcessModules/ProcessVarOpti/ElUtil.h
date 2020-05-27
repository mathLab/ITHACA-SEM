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

#include <LibUtilities/BasicUtils/Thread.h>

#include <NekMeshUtils/Module/Module.h>

#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/Interpolator.h>

typedef Nektar::LibUtilities::PtsFieldSharedPtr PtsFieldSharedPtr;

namespace Nektar
{
namespace Utilities
{

using namespace NekMeshUtils;

struct DerivUtil;
struct Residual;

typedef std::shared_ptr<DerivUtil> DerivUtilSharedPtr;
typedef std::shared_ptr<Residual> ResidualSharedPtr;

class ElUtilJob;

class ElUtil : public std::enable_shared_from_this<ElUtil>
{
public:
    ElUtil(ElementSharedPtr e, DerivUtilSharedPtr d,
           ResidualSharedPtr, int n, int o);

    ElUtilJob *GetJob(bool update = false);

    int GetId()
    {
        return m_el->GetId();
    }

    // Leaving these varibles as public for sake of efficiency
    std::vector<std::vector<NekDouble *> > nodes;
    std::vector<std::vector<NekDouble> > maps, mapsStd;

    void Evaluate();
    void InitialMinJac();

    ElementSharedPtr GetEl()
    {
        return m_el;
    }

    int NodeId(int in)
    {
        return m_idmap[in];
    }

    NekDouble GetScaledJac()
    {
        return m_scaledJac;
    }

    NekDouble &GetMinJac()
    {
        return m_minJac;
    }

    void SetScaling(LibUtilities::Interpolator interp)
    {
        m_interp = interp;
        UpdateMapping();
    }

    void UpdateMapping();

private:

    void MappingIdealToRef();

    ElementSharedPtr m_el;
    int m_dim;
    int m_mode;
    int m_order;
    std::map<int,int> m_idmap;

    NekDouble m_scaledJac;
    NekDouble m_minJac;

    PtsFieldSharedPtr m_interpField;
    LibUtilities::Interpolator m_interp;

    DerivUtilSharedPtr m_derivUtil;
    ResidualSharedPtr m_res;

    // Initial maps
    std::vector<std::vector<NekDouble>> m_maps, m_mapsStd;
};
typedef std::shared_ptr<ElUtil> ElUtilSharedPtr;

class ElUtilJob : public Thread::ThreadJob
{
public:
    ElUtilJob(ElUtil* e, bool update) : el(e), m_update(update) {}

    void Run()
    {
        el->Evaluate();

        if (m_update)
        {
            el->UpdateMapping();
        }
    }
private:
    ElUtil* el;
    bool m_update;
};

}
}

#endif
