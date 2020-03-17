///////////////////////////////////////////////////////////////////////////////
//
// File FilterBenchmark.h
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
// Description: Outputs time when solution first exceeds a threshold value.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERBENCHMARK_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERBENCHMARK_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{

/// Records activation and repolarisation times
class FilterBenchmark : public SolverUtils::Filter
{
public:
    friend class MemoryManager<FilterBenchmark>;

    /// Creates an instance of this class
    static SolverUtils::FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams)
    {
        SolverUtils::FilterSharedPtr p =
            MemoryManager<FilterBenchmark>::AllocateSharedPtr(pSession,
                                                        pEquation, pParams);
        return p;
    }

    /// Name of the class.
    static std::string className;

    /// Construct the benchmark filter.
    FilterBenchmark(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams);

    /// Destructor for the benchmark filter.
    virtual ~FilterBenchmark();

protected:
    /// Initialises the benchmark filter and allocates storage.
    virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    /// Update recorded times.
    virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    /// Finalises the benchmark filter and write out recorded data.
    virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    /// Identifies that the benchmark filter is time dependent.
    virtual bool v_IsTimeDependent();

private:
    /// Storage for activation and repolarisation times.
    std::vector<Array<OneD, NekDouble> > m_threshold;
    /// Number of activations and repolarisations detected for each point.
    Array<OneD, int>                     m_idx;
    /// Indicates if the previous event was an activation or repolarisation.
    Array<OneD, int>                     m_polarity;
    /// Time at which to start detecting activations and repolarisations.
    NekDouble                            m_startTime;
    /// Value at which tissue is considered active.
    NekDouble                            m_thresholdValue;
    /// Initial time to use in storage array.
    NekDouble                            m_initialValue;
    /// Filename of output files.
    std::string                          m_outputFile;
    /// FieldIO object used for writing output files.
    LibUtilities::FieldIOSharedPtr       m_fld;
};

}

#endif /* FILTERTHRESHOLDMAX_H_ */
