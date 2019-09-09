///////////////////////////////////////////////////////////////////////////////
//
// File SessionFunction.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Session Function
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_SESSIONFUNCTION_H
#define NEKTAR_SOLVERUTILS_SESSIONFUNCTION_H

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/CsvIO.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
namespace FieldUtils
{
class Interpolator;
}

namespace SolverUtils
{

class SessionFunction
{
public:
    /// Representation of a FUNCTION defined in the session xml file.
    SOLVER_UTILS_EXPORT SessionFunction(
        const LibUtilities::SessionReaderSharedPtr &session,
        const MultiRegions::ExpListSharedPtr &field,
        std::string functionName,
        bool toCache = false);

    /// Evaluates a function defined in the xml session file at each quadrature point.
    SOLVER_UTILS_EXPORT void Evaluate(
        Array<OneD, Array<OneD, NekDouble> > &pArray,
        const NekDouble pTime = 0.0,
        const int domain      = 0);

    /// Evaluates a function defined in the xml session file at each quadrature point.
    SOLVER_UTILS_EXPORT void Evaluate(
        std::vector<std::string> pFieldNames,
        Array<OneD, Array<OneD, NekDouble> > &pArray,
        const NekDouble &pTime = 0.0,
        const int domain       = 0);

    /// Evaluates a function defined in the xml session file at each quadrature point.
    SOLVER_UTILS_EXPORT void Evaluate(
        std::vector<std::string> pFieldNames,
        Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &pTime = 0.0,
        const int domain       = 0);

    // Evaluates a function defined in the xml session file at each quadrature point.
    SOLVER_UTILS_EXPORT void Evaluate(std::string pFieldName,
                                      Array<OneD, NekDouble> &pArray,
                                      const NekDouble &pTime = 0.0,
                                      const int domain       = 0);

    /// Provide a description of a function for a given field name.
    SOLVER_UTILS_EXPORT std::string Describe(std::string pFieldName,
                                             const int domain = 0);

    SOLVER_UTILS_EXPORT const LibUtilities::SessionReaderSharedPtr & GetSession()
    {
        return m_session;
    }

    SOLVER_UTILS_EXPORT const MultiRegions::ExpListSharedPtr & GetExpansion()
    {
        return m_field;
    }

private:
    /// The session reader
    LibUtilities::SessionReaderSharedPtr m_session;
    /// The expansion we want to evaluate this function for
    MultiRegions::ExpListSharedPtr m_field;
    // Name of this function
    std::string m_name;
    /// Store resulting arrays (and interpolators)
    bool m_toCache;
    /// Last time the cache for this variable & domain combo was updated
    std::map<std::pair<std::string, int>, NekDouble> m_lastCached;
    /// Interpolator for pts file input for a variable & domain combination
    std::map<std::string, FieldUtils::Interpolator> m_interpolators;
    /// Cached result arrays
    std::map<std::pair<std::string, int>, Array<OneD, NekDouble> > m_arrays;

    // Evaluates a function from expression
    SOLVER_UTILS_EXPORT void EvaluateExp(std::string pFieldName,
                                         Array<OneD, NekDouble> &pArray,
                                         const NekDouble &pTime = 0.0,
                                         const int domain       = 0);

    // Evaluates a function from fld file
    SOLVER_UTILS_EXPORT void EvaluateFld(std::string pFieldName,
                                         Array<OneD, NekDouble> &pArray,
                                         const NekDouble &pTime = 0.0,
                                         const int domain       = 0);

    /// Evaluates a function from pts file
    SOLVER_UTILS_EXPORT void EvaluatePts(std::string pFieldName,
                                         Array<OneD, NekDouble> &pArray,
                                         const NekDouble &pTime = 0.0,
                                         const int domain       = 0);

    SOLVER_UTILS_EXPORT void PrintProgressbar(const int position,
                                              const int goal) const
    {
        LibUtilities::PrintProgressbar(position, goal, "Interpolating");
    }
};

typedef std::shared_ptr<SessionFunction> SessionFunctionSharedPtr;
static SessionFunctionSharedPtr NullSessionFunction;
}
}

#endif
