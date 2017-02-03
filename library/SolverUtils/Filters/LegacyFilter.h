///////////////////////////////////////////////////////////////////////////////
//
// File LegacyFilter.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Base class for legacy filters.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_LEGACYFILTER_H
#define NEKTAR_SOLVERUTILS_LEGACYFILTER_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{

class LegacyFilter : public Filter
{
public:
    SOLVER_UTILS_EXPORT LegacyFilter(
        const LibUtilities::SessionReaderSharedPtr &pSession);

    SOLVER_UTILS_EXPORT virtual ~LegacyFilter();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
        const Array<OneD, const Array<OneD, NekDouble> > &coeffs,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
        const Array<OneD, const Array<OneD, NekDouble> > &coeffs,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
        const Array<OneD, const Array<OneD, NekDouble> > &coeffs,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time);

    SOLVER_UTILS_EXPORT virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) = 0;
    SOLVER_UTILS_EXPORT virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) = 0;
    SOLVER_UTILS_EXPORT virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time) = 0;

private:
    void SolvedVarsOnly(
        const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        Array<OneD, MultiRegions::ExpListSharedPtr> &sFields);
};

inline void LegacyFilter::v_Initialise(
    const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
    const Array<OneD, const Array<OneD, NekDouble> > &coeffs,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    Array<OneD, MultiRegions::ExpListSharedPtr> sFields;
    SolvedVarsOnly(fieldMetaDataMap, pFields, sFields);
    v_Initialise(sFields, time);
}

inline void LegacyFilter::v_Update(
    const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
    const Array<OneD, const Array<OneD, NekDouble> > &coeffs,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    Array<OneD, MultiRegions::ExpListSharedPtr> sFields;
    SolvedVarsOnly(fieldMetaDataMap, pFields, sFields);
    v_Update(sFields, time);
}

inline void LegacyFilter::v_Finalise(
    const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
    const Array<OneD, const Array<OneD, NekDouble> > &coeffs,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    Array<OneD, MultiRegions::ExpListSharedPtr> sFields;
    SolvedVarsOnly(fieldMetaDataMap, pFields, sFields);
    v_Finalise(sFields, time);
}
}
}
#endif /* NEKTAR_SOLVERUTILS_LEGACYFILTER_H */
