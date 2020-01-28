///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationAlgorithm.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2019 Division of Applied Mathematics, Brown University (USA),
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
// Description: Header file of time integration algorithm base class
//
// The TimeIntegrationAlgorithm class should only be used by the
// TimeIntegrationScheme class.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationTypes.h>

#define LUE LIB_UTILITIES_EXPORT

namespace Nektar
{
namespace LibUtilities
{

class TimeIntegrationAlgorithm
{

public:
    TimeIntegrationAlgorithm(const TimeIntegrationScheme *parent)
        : m_parent(parent)
    {
    }
    virtual ~TimeIntegrationAlgorithm()
    {
    }

    inline TimeIntegrationSchemeType GetIntegrationSchemeType() const
    {
        return m_schemeType;
    }

    /**
     * \brief This function initialises the time integration
     * scheme
     *
     */
    LUE virtual TimeIntegrationSolutionSharedPtr InitializeData(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op) = 0;

    /**
     * \brief Explicit integration of an ODE.
     *
     */
    LUE virtual ConstDoubleArray &TimeIntegrate(
        const NekDouble deltaT, TimeIntegrationSolutionSharedPtr &y,
        const TimeIntegrationSchemeOperators &op) = 0;

    LUE virtual void TimeIntegrate(const NekDouble deltaT, ConstTripleArray &y_old,
                           ConstSingleArray &t_old, TripleArray &y_new,
                           SingleArray &t_new,
                           const TimeIntegrationSchemeOperators &op) = 0;


    /// Parent scheme object
    const TimeIntegrationScheme *m_parent{nullptr};
  
    /// Type of time integration scheme (Explicit, Implicit, IMEX, etc)
    TimeIntegrationSchemeType m_schemeType{eNoTimeIntegrationSchemeType};
  
    int  m_nvars;       ///< The number of variables in integration scheme.
    int  m_npoints;     ///< The size of inner data which is stored for reuse.

}; // end class TimeIntegrationAlgorithm

} // end of namespace LibUtilities
} // end of namespace Nektar
