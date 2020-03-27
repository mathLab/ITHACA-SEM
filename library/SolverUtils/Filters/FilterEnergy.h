///////////////////////////////////////////////////////////////////////////////
//
// File FilterEnergy.h
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
// Description: Outputs solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERENERGY_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERENERGY_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{
class FilterEnergy : public Filter
{
public:
    /// Creates an instance of this class
    static SolverUtils::FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams)
    {
        SolverUtils::FilterSharedPtr p = MemoryManager<FilterEnergy>
                            ::AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterEnergy(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const ParamMap                             &pParams);
    SOLVER_UTILS_EXPORT ~FilterEnergy();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pField,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual bool v_IsTimeDependent();

private:
    unsigned int                m_index;
    unsigned int                m_outputFrequency;
    std::ofstream               m_outFile;
    bool                        m_homogeneous;
    NekDouble                   m_homogeneousLength;
    NekDouble                   m_area;
    LibUtilities::CommSharedPtr m_comm;
    Array<OneD, unsigned int>   m_planes;
};
}
}

#endif
