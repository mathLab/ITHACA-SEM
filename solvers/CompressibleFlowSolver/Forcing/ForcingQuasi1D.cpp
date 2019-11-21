///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingQuasi1D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Forcing for quasi-1D nozzle flow.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/Forcing/ForcingQuasi1D.h>

using namespace std;

namespace Nektar
{
std::string ForcingQuasi1D::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("Quasi1D",
                                    ForcingQuasi1D::create,
                                    "Quasi-1D nozzle Forcing");

ForcingQuasi1D::ForcingQuasi1D(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

void ForcingQuasi1D::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
    m_NumVariable = pNumForcingFields;
    m_varConv     = MemoryManager<VariableConverter>::AllocateSharedPtr(
                    m_session, 1);

    ASSERTL0( pFields[0]->GetGraph()->GetSpaceDimension() == 1,
              "ForcingQuasi1D requires a 1D problem.");

    const TiXmlElement* funcNameElmt = pForce->FirstChildElement("AREAFCN");
    if(!funcNameElmt)
    {
        ASSERTL0(funcNameElmt, "Requires AREAFCN tag "
                 "specifying function name which prescribes nozzle area.");
    }

    string    funcName = funcNameElmt->GetText();
    ASSERTL0(m_session->DefinesFunction(funcName),
             "Function '" + funcName + "' not defined.");

    // Evaluate geometrical term -Ax/A for forcing
    m_geomFactor = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(), 0.0);
    Array<OneD, NekDouble> tmp (pFields[0]->GetTotPoints(), 0.0);

    std::string  sFieldStr   = m_session->GetVariable(0);
    ASSERTL0(m_session->DefinesFunction(funcName, sFieldStr),
             "Variable '" + sFieldStr + "' not defined.");
    GetFunction(pFields, m_session, funcName, true)
        ->Evaluate(sFieldStr, m_geomFactor, 0.0);

    // Check if DADXFCN is defined
    const TiXmlElement* dAFuncNameElmt = pForce->FirstChildElement("DADXFCN");
    if(dAFuncNameElmt)
    {
        funcName = dAFuncNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(funcName),
             "Function '" + funcName + "' not defined.");
        ASSERTL0(m_session->DefinesFunction(funcName, sFieldStr),
             "Variable '" + sFieldStr + "' not defined.");
        GetFunction(pFields, m_session, funcName, true)
            ->Evaluate(sFieldStr, tmp, 0.0);
    }
    else
    {
        // Numerically evaluate dA/dX
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],
                              m_geomFactor, tmp);
    }

    Vmath::Vdiv(pFields[0]->GetTotPoints(), tmp, 1,
                                            m_geomFactor, 1,
                                            m_geomFactor, 1);
    Vmath::Neg(pFields[0]->GetTotPoints(), m_geomFactor, 1);

    // Project m_geomFactor to solution space
    Array<OneD, NekDouble> tmpCoeff (pFields[0]->GetNcoeffs(), 0.0);
    pFields[0]->FwdTrans_IterPerExp(m_geomFactor, tmpCoeff);
    pFields[0]->BwdTrans(tmpCoeff, m_geomFactor);

    m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
    for (int i = 0; i < m_NumVariable; ++i)
    {
        m_Forcing[i] = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(), 0.0);
    }
}

void ForcingQuasi1D::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  pFields,
        const Array<OneD, Array<OneD, NekDouble> >&         inarray,
              Array<OneD, Array<OneD, NekDouble> >&         outarray,
        const NekDouble&                                    time)
{
    boost::ignore_unused(time);

    int nPoints = pFields[0]->GetTotPoints();

    // Get (E+p)
    Array<OneD, NekDouble> tmp (nPoints, 0.0);
    m_varConv->GetPressure(inarray, tmp);
    Vmath::Vadd(nPoints, tmp, 1,
                    inarray[2], 1, tmp, 1);

    // F-rho = -Ax/A *rhou
    Vmath::Vmul(nPoints, m_geomFactor, 1,
                         inarray[1], 1, m_Forcing[0], 1);

    // F-rhou = -Ax/A *rhou*u
    Vmath::Vmul(nPoints, m_geomFactor, 1,
                         inarray[1], 1, m_Forcing[1], 1);
    Vmath::Vmul(nPoints, m_Forcing[1], 1,
                         inarray[1], 1, m_Forcing[1], 1);
    Vmath::Vdiv(nPoints, m_Forcing[1], 1,
                         inarray[0], 1, m_Forcing[1], 1);

    // F-E = -Ax/A *(E+p)*u
    Vmath::Vmul(nPoints, m_geomFactor, 1,
                         inarray[1], 1, m_Forcing[2], 1);
    Vmath::Vmul(nPoints, m_Forcing[2], 1,
                         tmp, 1, m_Forcing[2], 1);
    Vmath::Vdiv(nPoints, m_Forcing[2], 1,
                         inarray[0], 1, m_Forcing[2], 1);

    // Apply forcing
    for (int i = 0; i < m_NumVariable; i++)
    {
        Vmath::Vadd(nPoints, outarray[i], 1,
                    m_Forcing[i], 1, outarray[i], 1);
    }
}

}
