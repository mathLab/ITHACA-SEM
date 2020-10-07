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

#include <SolverUtils/Filters/FilterAverageFieldsOutputVariables.h>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterAverageFieldsOutputVariables::className =
        GetFilterFactory().RegisterCreatorFunction(
                "AverageFieldsOutputVariables",
                FilterAverageFieldsOutputVariables::create);

FilterAverageFieldsOutputVariables::FilterAverageFieldsOutputVariables(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams)
    : FilterAverageFields(pSession, pEquation, pParams)
{
}

FilterAverageFieldsOutputVariables::~FilterAverageFieldsOutputVariables()
{
}

void FilterAverageFieldsOutputVariables::v_FillVariablesName(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    int nfield = pFields.size();
    std::vector<Array<OneD, NekDouble> > fieldcoeffs(
                pFields.size());
    std::vector<std::string> variables(pFields.size());
    for (int n = 0; n < nfield; ++n)
    {
        m_variables[n] = pFields[n]->GetSession()->GetVariable(n);
    }

    // Need to create a dummy coeffs vector to get extra variables names...
    std::vector<Array<OneD, NekDouble> > coeffs(nfield);
    for (int n = 0; n < nfield; ++n)
    {
        coeffs[n] = pFields[n]->GetCoeffs();
    }
    // Get extra variables
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");
    equ->ExtraFldOutput(coeffs, m_variables);
}

void FilterAverageFieldsOutputVariables::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    const NekDouble &time)
{
    boost::ignore_unused(time);
    std::vector<std::string> variables(pFields.size());

    ASSERTL0(fieldcoeffs.size()==m_outFields.size(),
        "FilterAverageFieldsOutputVariables:: fieldcoeffs.size()!=m_outFields.size()");
    
    for(int n = 0; n < m_outFields.size(); ++n)
    {
        Vmath::Vadd(m_outFields[n].size(),
                    fieldcoeffs[n],
                    1,
                    m_outFields[n],
                    1,
                    m_outFields[n],
                    1);
    }
}

}
}
