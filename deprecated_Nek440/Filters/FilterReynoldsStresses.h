///////////////////////////////////////////////////////////////////////////////
//
// File FilterReynoldsStresses.h
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

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERREYNOLDSSTRESSES_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERREYNOLDSSTRESSES_H

#include <SolverUtils/Filters/FilterFieldConvert.h>

namespace Nektar
{
namespace SolverUtils
{
class FilterReynoldsStresses : public FilterFieldConvert
{
public:
    friend class MemoryManager<FilterReynoldsStresses>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p =
            MemoryManager<FilterReynoldsStresses>::AllocateSharedPtr(pSession,
                                                                     pParams);
        return p;
    }

    /// Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT
    FilterReynoldsStresses(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const std::map<std::string, std::string> &pParams);
    SOLVER_UTILS_EXPORT ~FilterReynoldsStresses();

protected:
    virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    virtual void v_FillVariablesName(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
    virtual void v_ProcessSample(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    virtual void v_PrepareOutput(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    virtual NekDouble v_GetScale();
    virtual std::string v_GetFileSuffix()
    {
        return "_stress";
    }

    std::vector<Array<OneD, NekDouble> > m_fields;
    std::vector<Array<OneD, NekDouble> > m_delta;
    NekDouble m_alpha;
    bool m_movAvg;
};
}
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERREYNOLDSSTRESES_H */
