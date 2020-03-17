///////////////////////////////////////////////////////////////////////////////
//
// File Filter.h
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
// Description: Base class for filters.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTER_FILTER_H
#define NEKTAR_SOLVERUTILS_FILTER_FILTER_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/EquationSystem.h>

#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
namespace SolverUtils
{
class Filter;
class EquationSystem;

/// A shared pointer to a Driver object
typedef std::shared_ptr<Filter> FilterSharedPtr;

/// Datatype of the NekFactory used to instantiate classes derived from
/// the Driver class.
typedef LibUtilities::NekFactory<
    std::string, Filter,
    const LibUtilities::SessionReaderSharedPtr&,
    const std::weak_ptr<EquationSystem>&,
    const std::map<std::string, std::string>&
    > FilterFactory;
SOLVER_UTILS_EXPORT FilterFactory& GetFilterFactory();

class Filter
{
public:
    typedef std::map<std::string, std::string> ParamMap;
    SOLVER_UTILS_EXPORT Filter(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation);
    SOLVER_UTILS_EXPORT virtual ~Filter();

    SOLVER_UTILS_EXPORT inline void Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);
    SOLVER_UTILS_EXPORT inline void Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);
    SOLVER_UTILS_EXPORT inline void Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);
    SOLVER_UTILS_EXPORT inline bool IsTimeDependent();

protected:
    LibUtilities::SessionReaderSharedPtr  m_session;
    const std::weak_ptr<EquationSystem> m_equ;

    virtual void v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time) = 0;
    virtual void v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time) = 0;
    virtual void v_Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time) = 0;
    virtual bool v_IsTimeDependent() = 0;
};

inline void Filter::Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    v_Initialise(pFields, time);
}

inline void Filter::Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    v_Update(pFields, time);
}

inline void Filter::Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    v_Finalise(pFields, time);
}

inline bool Filter::IsTimeDependent()
{
    return v_IsTimeDependent();
}
}
}
#endif /* NEKTAR_SOLVERUTILS_FILTER_FILTER_H */
