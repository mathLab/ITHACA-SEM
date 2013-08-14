///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingProgrammatic.cpp
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
// Description: Programmatic forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingProgrammatic.h>

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingProgrammatic::className = GetForcingFactory().
                                RegisterCreatorFunction("Programamtic",
                                                        ForcingProgrammatic::create,
                                                        "Programmatic Forcing");

    ForcingProgrammatic::ForcingProgrammatic(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : Forcing(pSession)
    {
    }

    Array<OneD, Array<OneD, NekDouble> >& ForcingProgrammatic::UpdateForces()
    {
        return m_Forcing;
    }

    void ForcingProgrammatic::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const TiXmlElement* pForce)
    {
        std::string m_SolverInfo = m_session->GetSolverInfo("SolverType");
        int nvariables = m_session->GetVariables().size();
        int nq         = pFields[0]->GetTotPoints();

        if (m_SolverInfo == "VelocityCorrectionScheme")
        {
            m_NumVariable = nvariables - 1; // e.g. (u v w p) for 3D case
        }
        if (m_SolverInfo == "CoupledLinearisedNS")
        {
            m_NumVariable = nvariables; // e.g. (u v w)  for 3D case
        }

        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (nq, 0.0);
        }
    }

    void ForcingProgrammatic::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
                  Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }
    }

}
}
