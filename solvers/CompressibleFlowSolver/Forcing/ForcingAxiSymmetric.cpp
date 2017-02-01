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
// License for the specific language governing rights and limitations under
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

#include <CompressibleFlowSolver/Forcing/ForcingAxiSymmetric.h>

using namespace std;

namespace Nektar
{
std::string ForcingAxiSymmetric::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("AxiSymmetric",
                                ForcingAxiSymmetric::create,
                                "Forcing for axi-symmetric flow (around x=0");

ForcingAxiSymmetric::ForcingAxiSymmetric(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : Forcing(pSession)
{
}

void ForcingAxiSymmetric::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
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
    int nPoints = pFields[0]->GetTotPoints();

    // Get (E+p)
    Array<OneD, NekDouble> tmp (nPoints, 0.0);
    m_varConv->GetPressure(inarray, tmp);
    Vmath::Vadd(nPoints, tmp, 1,
                    inarray[2], 1, tmp, 1);

    // F-rho = -1/r *rhou
    Vmath::Vmul(nPoints, m_geomFactor, 1,
                         inarray[1], 1, m_Forcing[0], 1);

    // F-rhou_i = -1/r *rhou_i * u_r
    for (int i = 1; i < m_NumVariable-1; ++i)
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

    // Apply forcing
    for (int i = 0; i < m_NumVariable; i++)
    {
        Vmath::Vadd(nPoints, outarray[i], 1,
                    m_Forcing[i], 1, outarray[i], 1);
    }
}

}
