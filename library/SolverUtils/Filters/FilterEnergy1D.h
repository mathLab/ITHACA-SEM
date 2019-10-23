///////////////////////////////////////////////////////////////////////////////
//
// File FilterEnergy1D.h
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
// Description: Outputs orthogonal expansion of 1D elements.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERENERGY1D_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERENERGY1D_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{
/**
 * @brief Filter for one-dimensional energy spectrum.
 *
 * This filter transforms the coefficients of each element into an
 * orthogonal basis and outputs these coefficients at each timestep. It
 * can be used to, for example, detect high-order oscillations in the
 * solution field.
 */
class FilterEnergy1D : public Filter
{
public:
    friend class MemoryManager<FilterEnergy1D>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const std::map<std::string, std::string>   &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterEnergy1D>
            ::AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterEnergy1D(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT ~FilterEnergy1D();

protected:
    virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
        const NekDouble                                         &time);
    virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
        const NekDouble                                         &time);
    virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
        const NekDouble                                         &time);
    virtual bool v_IsTimeDependent();

private:
    /// Output file.
    std::ofstream m_out;
    /// Output frequency.
    unsigned int m_outputFrequency;
    /// Current index counter.
    unsigned int m_index;
};
}
}

#endif
