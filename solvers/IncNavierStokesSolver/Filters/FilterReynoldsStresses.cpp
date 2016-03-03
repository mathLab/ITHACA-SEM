///////////////////////////////////////////////////////////////////////////////
//
// File FilterReynoldsStresses.cpp
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
// Description: Append Reynolds stresses to the average fields
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Filters/FilterReynoldsStresses.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterReynoldsStresses::className =
    GetFilterFactory().RegisterCreatorFunction("ReynoldsStresses",
                                               FilterReynoldsStresses::create);

/**
 * @class FilterReynoldsStresses
 * 
 * @brief Append Reynolds stresses to the average fields
 * 
 * This class appends the average fields with the Reynolds stresses of the form
 * \f$ \overline{u' v'} \f$. This is achieved by calculating 
 * \f$ C_{n} = \Sigma_{i=1}^{n} (u_i - \bar{u}_n)(v_i - \bar{v}_n)\f$
 * using the recursive relation:
 * 
 * \f[ C_{n} = C_{n-1} + \frac{n}{n-1} (u_n - \bar{u}_n)(v_n - \bar{v}_n) \f]
 * 
 * The FilterAverageFields base class then divides the result by n, leading
 * to the Reynolds stress.
 */
FilterReynoldsStresses::FilterReynoldsStresses(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::map<std::string, std::string> &pParams)
    : FilterSampler(pSession, pParams)
{
}

FilterReynoldsStresses::~FilterReynoldsStresses()
{
}

void FilterReynoldsStresses::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    int dim = pFields.num_elements() - 1;
    int nExtraFields = dim == 2 ? 3 : 6;
    int origFields = pFields.num_elements();

    // Fill name of variables
    for(int n = 0; n < origFields; ++n)
    {
        m_variables.push_back(pFields[n]->GetSession()->GetVariable(n));
    }
    if (dim == 2)
    {
        m_variables.push_back("uu");
        m_variables.push_back("uv");
        m_variables.push_back("vv");
    }
    else if (dim == 3)
    {
        m_variables.push_back("uu");
        m_variables.push_back("uv");
        m_variables.push_back("uw");
        m_variables.push_back("vv");
        m_variables.push_back("vw");
        m_variables.push_back("ww");
    }
    else
    {
        ASSERTL0(false, "Unsupported dimension");
    }

    // Allocate storage
    m_fields.resize(origFields + nExtraFields);
    m_delta.resize(dim);

    for (int n = 0; n < m_fields.size(); ++n)
    {
        m_fields[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }
    for (int n = 0; n < m_delta.size(); ++n)
    {
        m_delta[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }

    // Initialise output arrays
    FilterSampler::v_Initialise(pFields, time);
}

void FilterReynoldsStresses::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    int i, j, n;
    int nq             = pFields[0]->GetTotPoints();
    int dim            = pFields.num_elements() - 1;
    bool waveSpace     = pFields[0]->GetWaveSpace();
    NekDouble nSamples = (NekDouble) m_numSamples;
    NekDouble fac      = nSamples / (nSamples-1);

    Array<OneD, NekDouble> vel(nq);
    Array<OneD, NekDouble> tmp(nq);

    // Update original velocities in phys space and calculate (\bar{u} - u_n)
    for (n = 0; n < dim; ++n)
    {
        if (waveSpace)
        {
            pFields[n]->HomogeneousBwdTrans(pFields[n]->GetPhys(), vel);
        }
        else
        {
            vel = pFields[n]->GetPhys();
        }
        Vmath::Vadd(nq, vel, 1, m_fields[n], 1, m_fields[n], 1);
        Vmath::Svtvm(nq, 1.0/nSamples, m_fields[n], 1,
                     vel, 1, m_delta[n], 1);
    }
    // Update pressure (directly to outFields)
    Vmath::Vadd(m_outFields[dim].num_elements(), pFields[dim]->GetCoeffs(), 1, 
                m_outFields[dim], 1, m_outFields[dim], 1);

    // Ignore Reynolds stress for first sample (its contribution is zero)
    if( m_numSamples == 1)
    {
        return;
    }

    // Calculate C_{n} = C_{n-1} + fac * deltaI * deltaJ
    for (i = 0, n = dim+1; i < dim; ++i)
    {
        for (j = i; j < dim; ++j, ++n)
        {
            Vmath::Vmul(nq, m_delta[i], 1, m_delta[j], 1, tmp, 1);
            Vmath::Smul(nq, fac, tmp, 1, tmp, 1);
            Vmath::Vadd(nq, tmp, 1, m_fields[n], 1, m_fields[n], 1);
        }
    }
}

void FilterReynoldsStresses::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_scale = 1.0/m_numSamples;
    int dim = pFields.num_elements() - 1;

    // Set wavespace to false, as calculations were performed in physical space
    bool waveSpace  = pFields[0]->GetWaveSpace();
    pFields[0]->SetWaveSpace(false);

    // Forward transform and put into m_outFields (except pressure)
    for (int i = 0; i < m_fields.size(); ++i)
    {
        if(i != dim)
        {
            pFields[0]->FwdTrans_IterPerExp(m_fields[i], m_outFields[i]);
        }
    }

    //Restore waveSpace
    pFields[0]->SetWaveSpace(waveSpace);
}

bool FilterReynoldsStresses::v_IsTimeDependent()
{
    return true;
}
}
}
