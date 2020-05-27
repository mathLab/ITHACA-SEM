///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingAxiSymmetric.cpp
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
// Description: Forcing for axi-symmetric flow.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/Forcing/ForcingAxiSymmetric.h>

using namespace std;

namespace Nektar
{
std::string ForcingAxiSymmetric::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("AxiSymmetric",
                                ForcingAxiSymmetric::create,
                                "Forcing for axi-symmetric flow (around x=0)");

ForcingAxiSymmetric::ForcingAxiSymmetric(
                const LibUtilities::SessionReaderSharedPtr         &pSession,
                const std::weak_ptr<SolverUtils::EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

void ForcingAxiSymmetric::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
    boost::ignore_unused(pForce);

    int spacedim = pFields[0]->GetGraph()->GetSpaceDimension();
    int nPoints  = pFields[0]->GetTotPoints();

    m_NumVariable = pNumForcingFields;
    m_varConv     = MemoryManager<VariableConverter>::AllocateSharedPtr(
                    m_session, spacedim);

    // Get coordinates
    Array<OneD, Array<OneD, NekDouble> > coords(3);
    for (int i = 0; i < 3; i++)
    {
        coords[i]      = Array<OneD, NekDouble> (nPoints);
    }
    pFields[0]->GetCoords(coords[0], coords[1], coords[2]);

    // Calculate fac = -1/r if r!=0, fac = 0 if r == 0
    m_geomFactor = Array<OneD, NekDouble> (nPoints);
    for (int i = 0; i < nPoints; ++i)
    {
        if (coords[0][i] < NekConstants::kNekZeroTol)
        {
            m_geomFactor[i] = 0;
        }
        else
        {
            m_geomFactor[i] = -1.0/coords[0][i];
        }
    }

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

void ForcingAxiSymmetric::v_Apply(
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
                    inarray[m_NumVariable-1], 1, tmp, 1);

    // F-rho = -1/r *rhou
    Vmath::Vmul(nPoints, m_geomFactor, 1,
                         inarray[1], 1, m_Forcing[0], 1);

    // F-rhou_r = -1/r *rhou_r * u_r and F-rhou_y = -1/r *rhou_y * u_r
    for (int i = 1; i < 3; ++i)
    {
        Vmath::Vmul(nPoints, inarray[1], 1,
                             inarray[i], 1, m_Forcing[i], 1);
        Vmath::Vdiv(nPoints, m_Forcing[i], 1,
                             inarray[0], 1, m_Forcing[i], 1);
        Vmath::Vmul(nPoints, m_Forcing[i], 1,
                             m_geomFactor, 1, m_Forcing[i], 1);
    }

    // F-E = -1/r *(E+p)*u
    Vmath::Vmul(nPoints, inarray[1], 1,
                         tmp, 1, m_Forcing[m_NumVariable-1], 1);
    Vmath::Vdiv(nPoints, m_Forcing[m_NumVariable-1], 1,
                         inarray[0], 1, m_Forcing[m_NumVariable-1], 1);
    Vmath::Vmul(nPoints, m_Forcing[m_NumVariable-1], 1,
                         m_geomFactor, 1, m_Forcing[m_NumVariable-1], 1);

    // Swirl
    if (m_NumVariable == 5)
    {
        // F-rhou_r -= (-1/r) * rho * u_theta * u_theta
        Vmath::Vmul(nPoints, inarray[3], 1,
                             inarray[3], 1, tmp, 1);
        Vmath::Vdiv(nPoints, tmp, 1,
                             inarray[0], 1, tmp, 1);
        Vmath::Vmul(nPoints, tmp, 1,
                             m_geomFactor, 1, tmp, 1);
        Vmath::Vsub(nPoints, m_Forcing[1], 1,
                             tmp, 1, m_Forcing[1], 1);

        // F-rhou_theta = 2 * (-1/r *rhou_theta * u_r)
        Vmath::Vmul(nPoints, inarray[1], 1,
                             inarray[3], 1, m_Forcing[3], 1);
        Vmath::Vdiv(nPoints, m_Forcing[3], 1,
                             inarray[0], 1, m_Forcing[3], 1);
        Vmath::Vmul(nPoints, m_Forcing[3], 1,
                             m_geomFactor, 1, m_Forcing[3], 1);
        Vmath::Smul(nPoints, 2.0,
                             m_Forcing[3], 1, m_Forcing[3], 1);
    }

    // Apply forcing
    for (int i = 0; i < m_NumVariable; i++)
    {
        Vmath::Vadd(nPoints, outarray[i], 1,
                    m_Forcing[i], 1, outarray[i], 1);
    }
}

}
