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
    : FilterAverageFields(pSession, pParams)
{
}

FilterReynoldsStresses::~FilterReynoldsStresses()
{
}

void FilterReynoldsStresses::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    FilterAverageFields::v_Initialise(pFields, time);

    int dim = pFields.num_elements() - 1;
    int nExtraFields = dim == 2 ? 3 : 6;
    int origFields = pFields.num_elements();

    m_fields.resize(nExtraFields);
    m_delta.resize(dim);
    m_avgFields.resize(m_avgFields.size() + nExtraFields);

    for (int n = 0; n < nExtraFields; ++n)
    {
        m_avgFields[n + origFields] =
            Array<OneD, NekDouble>(pFields[0]->GetNcoeffs(), 0.0);
        m_fields[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
    }

    for (int n = 0; n < dim; ++n)
    {
        m_delta[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(), 0.0);
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
}

void FilterReynoldsStresses::v_AddExtraFields(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Ignore first sample (its contribution is zero)
    if( m_numAverages == 1)
    {
        return;
    }

    int i, j, n;
    int nq          = pFields[0]->GetTotPoints();
    int dim         = pFields.num_elements() - 1;
    int extraFields = dim == 2 ? 3 : 6;
    int origFields  = pFields.num_elements();

    Array<OneD, NekDouble> tmpCoeff(pFields[0]->GetNcoeffs());
    Array<OneD, NekDouble> tmpPhys(pFields[0]->GetTotPoints());

    // Constant n/(n-1)
    NekDouble fac = ((NekDouble) m_numAverages) / (m_numAverages-1);

    // Deal with Homogeneous case
    bool      waveSpace = pFields[0]->GetWaveSpace();
    pFields[0]->SetWaveSpace(false);

    // delta_i = (ui_n - \bar{ui}_n)
    for (i = 0; i < dim; ++i)
    {
        // phys values of \bar{ui}
        pFields[0]->BwdTrans(m_avgFields[i], m_delta[i]);
        Vmath::Smul(nq, 1.0/m_numAverages, m_delta[i], 1, m_delta[i], 1);

        // Put new velocity in physical space in homogeneous case
        if (waveSpace)
        {
            pFields[i]->HomogeneousBwdTrans(pFields[i]->GetPhys(), tmpPhys);
        }
        else
        {
            tmpPhys = pFields[i]->GetPhys();
        }
        // Calculate delta
        Vmath::Vsub(nq, tmpPhys, 1, m_delta[i], 1, m_delta[i], 1);
    }

    // Calculate correction: C_{n} - C_{n-1} = fac * deltaI * deltaJ
    for (i = 0, n = 0; i < dim; ++i)
    {
        for (j = 0; j < i; ++j, ++n)
        {
            Vmath::Vmul(nq, m_delta[i], 1, m_delta[j], 1, m_fields[n], 1);
            Vmath::Smul(nq, fac, m_fields[n], 1, m_fields[n], 1);
        }
    }

    // Forward transform and put into m_avgFields
    for (i = 0; i < extraFields; ++i)
    {
        pFields[0]->FwdTrans_IterPerExp(m_fields[i], tmpCoeff);
        Vmath::Vadd(m_avgFields[i + origFields].num_elements(), tmpCoeff, 1,
                    m_avgFields[i + origFields], 1,
                    m_avgFields[i + origFields], 1);
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
