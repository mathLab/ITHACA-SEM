////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessStreamFunction.cpp
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
//  Description: Computes stream-function for 2D field.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "ProcessStreamFunction.h"

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessStreamFunction::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "streamfunction"),
        ProcessStreamFunction::create,
        "Computes stream-function for 2D field.");

ProcessStreamFunction::ProcessStreamFunction(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_f->m_declareExpansionAsContField = true;
}

ProcessStreamFunction::~ProcessStreamFunction()
{
}

void ProcessStreamFunction::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int nfields  = m_f->m_variables.size();
    // Append field name
    m_f->m_variables.push_back("StreamFunc");

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // Check if dimension is correct
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = expdim + m_f->m_numHomogeneousDir;
    ASSERTL0(spacedim == 2, "Stream function can only be obtained in 2D.");

    // Allocate arrays
    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, NekDouble> vx(npoints);
    Array<OneD, NekDouble> uy(npoints);
    Array<OneD, NekDouble> vortNeg(npoints);

    // Resize expansion
    m_f->m_exp.resize(nfields + 1);
    m_f->m_exp[nfields] = m_f->AppendExpList(m_f->m_numHomogeneousDir);

    // Calculate vorticity: -W_z = Uy - Vx
    m_f->m_exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],
                             m_f->m_exp[0]->GetPhys(),uy);
    m_f->m_exp[1]->PhysDeriv(MultiRegions::DirCartesianMap[0],
                             m_f->m_exp[1]->GetPhys(),vx);

    Vmath::Vsub(npoints, uy, 1, vx, 1, vortNeg, 1);

    // Calculate the Stream Function as the solution of the
    //   Poisson equation: \nabla^2 StreamFunction = -Vorticity
    //   Use same boundary conditions as v
    StdRegions::ConstFactorMap factor;
    factor[StdRegions::eFactorLambda] = 0.0;

    m_f->m_exp[1]->HelmSolve(vortNeg,
                            m_f->m_exp[nfields]->UpdateCoeffs(), factor);
    m_f->m_exp[1]->BwdTrans(m_f->m_exp[nfields]->GetCoeffs(),
                            m_f->m_exp[nfields]->UpdatePhys());
}

}
}
