///////////////////////////////////////////////////////////////////////////////
//
// File FilterAverageField.cpp
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
// Description: Average solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterAverageFields.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterAverageFields::className =
        GetFilterFactory().RegisterCreatorFunction(
                "AverageFields", FilterAverageFields::create);

FilterAverageFields::FilterAverageFields(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams)
    : FilterSampler(pSession, pParams)
{
}

FilterAverageFields::~FilterAverageFields()
{
}

void FilterAverageFields::v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    int nfield = pFields.num_elements();
    m_variables.resize(pFields.num_elements());
    // Fill name of variables
    for (int n = 0; n < nfield; ++n)
    {
        m_variables[n] = pFields[n]->GetSession()->GetVariable(n);
    }
    // Now let FilterSampler initialise the output
    FilterSampler::v_Initialise(pFields, time);
}

void FilterAverageFields::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    for(int n = 0; n < pFields.num_elements(); ++n)
    {
        Vmath::Vadd(m_outFields[n].num_elements(),
                    pFields[n]->GetCoeffs(),
                    1,
                    m_outFields[n],
                    1,
                    m_outFields[n],
                    1);
    }
}

void FilterAverageFields::v_PrepareOutput(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_scale = 1.0 / m_numSamples;
}

bool FilterAverageFields::v_IsTimeDependent()
{
    return true;
}
}
}
