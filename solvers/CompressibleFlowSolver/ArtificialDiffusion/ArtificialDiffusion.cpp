///////////////////////////////////////////////////////////////////////////////
//
// File: ArtificialDiffusion.cpp
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
// Description: Abstract base class for compressible solver artificial diffusion
//              used for shock capturing artificial diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#include "ArtificialDiffusion.h"

using namespace std;

namespace Nektar
{
ArtificialDiffusionFactory& GetArtificialDiffusionFactory()
{
    typedef Loki::SingletonHolder<ArtificialDiffusionFactory,
                                  Loki::CreateUsingNew,
                                  Loki::NoDestroy,
                                  Loki::SingleThreaded> Type;
    return Type::Instance();
}

ArtificialDiffusion::ArtificialDiffusion(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const int spacedim)
        : m_session(pSession),
        m_fields(pFields)
{
    m_session->LoadParameter ("FL",            m_FacL,          0.0);
    m_session->LoadParameter ("FH",            m_FacH,          0.0);
    m_session->LoadParameter ("hFactor",       m_hFactor,       1.0);
    m_session->LoadParameter ("C1",            m_C1,            3.0);
    m_session->LoadParameter ("C2",            m_C2,            5.0);
    m_session->LoadParameter ("mu0",           m_mu0,           1.0);
    m_session->LoadParameter ("Skappa",        m_Skappa,        -2.048);
    m_session->LoadParameter ("Kappa",         m_Kappa,         0.0);

    // Create auxiliary object to convert variables
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                m_session, spacedim);

    m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance("LDG", "LDG");
    m_diffusion->SetArtificialDiffusionVector(
                        &ArtificialDiffusion::GetArtificialViscosity, this);
    m_diffusion->InitObject (m_session, m_fields);
}

/**
 *
 */
void ArtificialDiffusion::DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    v_DoArtificialDiffusion(inarray, outarray);
}

/**
 *
 */
void ArtificialDiffusion::v_DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    int i;
    int nvariables = inarray.num_elements();
    int npoints    = m_fields[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

    for (i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff);

    for (i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(npoints,
                    outarray[i], 1,
                    outarrayDiff[i], 1,
                    outarray[i], 1);
    }
}

/**
 *
 */
void ArtificialDiffusion::GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu)
{
    v_GetArtificialViscosity(physfield, mu);
}

}
