///////////////////////////////////////////////////////////////////////////////
//
// File FilterFieldConvert.h
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
// Description: Base clase for filters performing operations on samples
//              of the field.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERFIELDCONVERT_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERFIELDCONVERT_H

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/Filter.h>
#include <FieldUtils/Module.h>

using namespace Nektar::FieldUtils;

namespace Nektar
{
namespace SolverUtils
{
class FilterFieldConvert : public Filter
{
public:
    friend class MemoryManager<FilterFieldConvert>;

    /// Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const std::map<std::string, std::string> &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterFieldConvert>
                            ::AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    ///Name of the class
    static std::string className;
    
    SOLVER_UTILS_EXPORT FilterFieldConvert(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT virtual ~FilterFieldConvert();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_FillVariablesName(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_ProcessSample(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
              std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_PrepareOutput(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
    {
        boost::ignore_unused(pFields, time);
        // Do nothing by default
    }
    SOLVER_UTILS_EXPORT virtual NekDouble v_GetScale()
    {
        return 1.0;
    }
    SOLVER_UTILS_EXPORT virtual std::string v_GetFileSuffix()
    {
        return "_fc";
    }

    void OutputField(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        int dump = -1);

    SOLVER_UTILS_EXPORT virtual bool v_IsTimeDependent();
    
    void CreateModules(std::vector<std::string> &modcmds);

    void CreateFields(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields);

    void CheckModules(std::vector<ModuleSharedPtr> &modules);

    unsigned int m_numSamples;
    unsigned int m_outputFrequency;
    unsigned int m_sampleFrequency;
    std::string  m_outputFile;
    std::string  m_restartFile;
    unsigned int m_index;
    unsigned int m_outputIndex;

    // Phase sample parameters
    bool         m_phaseSample;
    NekDouble    m_phaseSamplePeriod;
    NekDouble    m_phaseSamplePhase;
    NekDouble    m_phaseTolerance;
    NekDouble    m_dt;

    std::vector<ModuleSharedPtr> m_modules;
    LibUtilities::FieldMetaDataMap m_fieldMetaData;
    std::vector<Array<OneD, NekDouble> > m_outFields;
    std::vector<std::string> m_variables;
    FieldSharedPtr m_f;
};
}
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERFIELDCONVERT_H */
