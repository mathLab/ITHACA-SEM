///////////////////////////////////////////////////////////////////////////////
//
// File FilterSampler.h
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
// Description: Base clase for filters performing operations on samples
//              of the field.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERSAMPLER_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERSAMPLER_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{
class FilterSampler : public Filter
{
public:
    friend class MemoryManager<FilterSampler>;

    SOLVER_UTILS_EXPORT FilterSampler(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const ParamMap &pParams);
    SOLVER_UTILS_EXPORT virtual ~FilterSampler();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_ProcessSample(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) = 0;
    SOLVER_UTILS_EXPORT virtual void v_PrepareOutput(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
    {
        // Do nothing by default
    }
    SOLVER_UTILS_EXPORT virtual std::string v_GetFileSuffix() = 0;

    void OutputField(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        int dump = -1);

    NekDouble m_scale;
    unsigned int m_numSamples;
    unsigned int m_outputFrequency;
    unsigned int m_sampleFrequency;
    unsigned int m_index;
    unsigned int m_outputIndex;
    std::string m_outputFile;
    LibUtilities::FieldIOSharedPtr m_fld;
    LibUtilities::FieldMetaDataMap m_fieldMetaData;
    std::vector<Array<OneD, NekDouble> > m_outFields;
    std::vector<std::string> m_variables;
};
}
}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERSAMPLER_H */
