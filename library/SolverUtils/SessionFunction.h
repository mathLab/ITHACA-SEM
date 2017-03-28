///////////////////////////////////////////////////////////////////////////////
//
// File PtsField.h
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
// Description: Session Function
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_SESSIONFUNCTION_H
#define NEKTAR_SOLVERUTILS_SESSIONFUNCTION_H

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
namespace Solverutils
{

struct loadedFldField
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
    std::vector<std::vector<NekDouble> > fieldData;
};

class SessionFunction : public boost::enable_shared_from_this<SessionFunction>
{
public:
    SOLVER_UTILS_EXPORT SessionFunction(
        LibUtilities::SessionReaderSharedPtr session,
        MultiRegions::ExpListSharedPtr field,
        std::string functionName);

    /// Evaluates a function as specified in the session file.
    SOLVER_UTILS_EXPORT void Evaluate(
        Array<OneD, Array<OneD, NekDouble> > &pArray,
        const NekDouble pTime = 0.0,
        const int domain      = 0);

    /// Populate given fields with the function from session.
    SOLVER_UTILS_EXPORT void Evaluate(
        std::vector<std::string> pFieldNames,
        Array<OneD, Array<OneD, NekDouble> > &pArray,
        const NekDouble &pTime = 0.0,
        const int domain       = 0);

    /// Populate given fields with the function from session.
    SOLVER_UTILS_EXPORT void Evaluate(
        std::vector<std::string> pFieldNames,
        Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &pTime = 0.0,
        const int domain       = 0);

    // Populate an array with a function variable from session.
    SOLVER_UTILS_EXPORT void Evaluate(std::string pFieldName,
                                      Array<OneD, NekDouble> &pArray,
                                      const NekDouble &pTime = 0.0,
                                      const int domain       = 0);

    SOLVER_UTILS_EXPORT std::string Describe(std::string pFieldName,
                                             const int domain = 0);

private:
    /// The session reader
    LibUtilities::SessionReaderSharedPtr m_session;
    /// The expansion we want to evaluate this function for
    MultiRegions::ExpListSharedPtr m_field;
    // Name of this function
    std::string m_functionName;
    // type of this function
    LibUtilities::FunctionType m_type;
    /// interpolator for pts file input
    FieldUtils::Interpolator m_interpolator;

    // Populate an array with a function variable from session.
    SOLVER_UTILS_EXPORT void EvaluateExp(std::string pFieldName,
                                         Array<OneD, NekDouble> &pArray,
                                         const NekDouble &pTime = 0.0,
                                         const int domain       = 0);

    // Populate an array with a function variable from session.
    SOLVER_UTILS_EXPORT void EvaluateFld(std::string pFieldName,
                                         Array<OneD, NekDouble> &pArray,
                                         const NekDouble &pTime = 0.0,
                                         const int domain       = 0);

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
}
}

#endif
