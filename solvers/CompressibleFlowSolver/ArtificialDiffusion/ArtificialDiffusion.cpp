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
    static ArtificialDiffusionFactory instance;
    return instance;
}

ArtificialDiffusion::ArtificialDiffusion(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const int spacedim)
        : m_session(pSession),
          m_fields(pFields)
{
    // Create auxiliary object to convert variables
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                m_session, spacedim);

    m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance("LDG", "LDG");
    m_diffusion->SetFluxVector(&ArtificialDiffusion::GetFluxVector, this);
    m_diffusion->InitObject (m_session, m_fields);

    // Get constant scaling
    m_session->LoadParameter("mu0", m_mu0, 1.0);

    // Init h/p scaling
    int nElements = m_fields[0]->GetExpSize();
    m_hOverP = Array<OneD, NekDouble>(nElements, 1.0);
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
    int nvariables = inarray.size();
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

void ArtificialDiffusion::DoArtificialDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    v_DoArtificialDiffusion_coeff(inarray, outarray);
}

void ArtificialDiffusion::v_DoArtificialDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray)
{
    size_t nvariables = inarray.size();
    size_t ncoeffs    = m_fields[0]->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble> > outarrayDiff {nvariables};

    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>{ncoeffs, 0.0};
    }

    m_diffusion->DiffuseCoeffs(nvariables, m_fields, inarray, outarrayDiff);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(ncoeffs,
                    outarray[i], 1,
                    outarrayDiff[i], 1,
                    outarray[i], 1);
    }
}

void ArtificialDiffusion::GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu)
{
    v_GetArtificialViscosity(physfield, mu);
}

/**
 *
 */
void ArtificialDiffusion::SetElmtHP(const Array<OneD, NekDouble> &hOverP)
{
    m_hOverP = hOverP;
}

/**
 * @brief Return the flux vector for the artificial viscosity operator.
 */
void ArtificialDiffusion::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble> > &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&qfield,
          Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&viscousTensor)
{
    unsigned int nDim = qfield.size();
    unsigned int nConvectiveFields = qfield[0].size();
    unsigned int nPts = qfield[0][0].size();

    // Get Artificial viscosity
    Array<OneD, NekDouble> mu{nPts, 0.0};
    GetArtificialViscosity(inarray, mu);

    // Compute viscous tensor
    for (unsigned int j = 0; j < nDim; ++j)
    {
        for (unsigned int i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Vmul(nPts, qfield[j][i], 1, mu , 1, viscousTensor[j][i], 1);
        }
    }
}


}
