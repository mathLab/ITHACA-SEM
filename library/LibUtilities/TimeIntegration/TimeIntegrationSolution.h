///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSolution.h
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
// Description: Header file of time integration solution class
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>

#define LUE LIB_UTILITIES_EXPORT

namespace Nektar
{
namespace LibUtilities
{

class TimeIntegrationSolution
{
public:
    typedef Array<OneD, Array<OneD, Array<OneD, NekDouble>>> TripleArray;
    typedef Array<OneD, Array<OneD, NekDouble>> DoubleArray;

    // Constructor for single step methods
    LUE TimeIntegrationSolution(const TimeIntegrationSchemeData *schemeData,
                                const DoubleArray &y, const NekDouble time,
                                const NekDouble timestep);

    // Constructor for multi-step methods
    LUE TimeIntegrationSolution(const TimeIntegrationSchemeData *schemeData,
                                const TripleArray &y,
                                const Array<OneD, NekDouble> &t);

    LUE TimeIntegrationSolution(const TimeIntegrationSchemeData *schemeData,
                                const unsigned int nvar,
                                const unsigned int npoints);

    LUE TimeIntegrationSolution(const TimeIntegrationSchemeData *schemeData);

    inline const TimeIntegrationSchemeData *GetIntegrationSchemeData() const
    {
        return m_schemeData;
    }

    std::string GetName() const;

    inline const TripleArray &GetSolutionVector() const
    {
        return m_solVector;
    }
    inline TripleArray &UpdateSolutionVector()
    {
        return m_solVector;
    }

    inline const DoubleArray &GetSolution() const
    {
        return m_solVector[0];
    }
    inline DoubleArray &UpdateSolution()
    {
        return m_solVector[0];
    }

    inline const Array<OneD, const NekDouble> &GetTimeVector() const
    {
        return m_t;
    }
    inline Array<OneD, NekDouble> &UpdateTimeVector()
    {
        return m_t;
    }

    inline NekDouble GetTime() const
    {
        return m_t[0];
    }
    int GetNsteps()
    {
        return m_schemeData->m_numsteps;
    }

    inline int GetFirstDim() const
    {
        return m_solVector[0].size();
    }
    inline int GetSecondDim() const
    {
        return m_solVector[0][0].size();
    }

    // Return the number of entries in the solution vector that correspond to
    // (multi-step) values.
    inline unsigned int GetNvalues() const
    {
        return m_schemeData->GetNmultiStepValues();
    }

    // Return the number of entries in the solution vector that correspond to
    // (multi-step) derivatives.
    inline unsigned int GetNderivs() const
    {
        return m_schemeData->GetNmultiStepDerivs();
    }

    // Returns an array which indicates to which time-level the entries in the
    // solution vector correspond.
    inline const Array<OneD, const unsigned int> &GetTimeLevelOffset()
    {
        return m_schemeData->GetTimeLevelOffset();
    }

    // returns the entry in the solution vector which corresponds to the
    // (multi-step) value at the time-level with specified offset
    inline DoubleArray &GetValue(const unsigned int timeLevelOffset)
    {
        int nMultiStepVals = m_schemeData->GetNmultiStepValues();
        const Array<OneD, const unsigned int> &offsetvec =
	    m_schemeData->GetTimeLevelOffset();

        for (int i = 0; i < nMultiStepVals; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                return m_solVector[i];
            }
        }
        ASSERTL1(false, "The solution vector of this scheme does not contain a "
                        "value at the requested time-level");
        return m_solVector[0];
    }

    // returns the entry in the solution vector which corresponds to the
    // (multi-step) derivative at the time-level with specified offset
    inline DoubleArray &GetDerivative(const unsigned int timeLevelOffset)
    {
        int nMultiStepVals = m_schemeData->GetNmultiStepValues();
        int size           = m_schemeData->m_numsteps;
        const Array<OneD, const unsigned int> &offsetvec =
	    m_schemeData->GetTimeLevelOffset();

        for (int i = nMultiStepVals; i < size; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                return m_solVector[i];
            }
        }
        ASSERTL1(false, "The solution vector of this scheme does not contain a "
                        "derivative at the requested time-level");
        return m_solVector[0];
    }

    // returns the time associated with the (multi-step) value at the time-level
    // with the given offset
    inline NekDouble GetValueTime(const unsigned int timeLevelOffset)
    {
        int nMultiStepVals = m_schemeData->GetNmultiStepValues();
        const Array<OneD, const unsigned int> &offsetvec =
	     m_schemeData->GetTimeLevelOffset();

        for (int i = 0; i < nMultiStepVals; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                return m_t[i];
            }
        }
        ASSERTL1(false, "The solution vector of this scheme does not contain a "
                        "value at the requested time-level");
        return m_t[0];
    }

    // sets the (multi-step) value and time in the solution
    // vector which corresponds to
    // the value at the time-level with specified offset
    inline void SetValue(const unsigned int timeLevelOffset,
                         const DoubleArray &y, const NekDouble t)
    {
        int nMultiStepVals = m_schemeData->GetNmultiStepValues();
        const Array<OneD, const unsigned int> &offsetvec =
	    m_schemeData->GetTimeLevelOffset();

        for (int i = 0; i < nMultiStepVals; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                m_solVector[i] = y;
                m_t[i]         = t;
                return;
            }
        }
    }

    // sets the (multi-step) derivative and time in the
    // solution vector which corresponds to
    // the derivative at the time-level with specified offset
    inline void SetDerivative(const unsigned int timeLevelOffset,
                              const DoubleArray &y, const NekDouble timestep)
    {
        int nMultiStepVals = m_schemeData->GetNmultiStepValues();
        int size           = m_schemeData->m_numsteps;
        const Array<OneD, const unsigned int> &offsetvec =
	    m_schemeData->GetTimeLevelOffset();

        for (int i = nMultiStepVals; i < size; i++)
        {
            if (timeLevelOffset == offsetvec[i])
            {
                m_solVector[i] = y;
                m_t[i]         = timestep;
                return;
            }
        }
    }

    // sets the soln Vector
    inline void SetSolVector(const int Offset, const DoubleArray &y)
    {
        m_solVector[Offset] = y;
    }

    // Rotate the solution vector
    // (i.e. updating without calculating/inserting new values)
    inline void RotateSolutionVector()
    {
        int nMultiStepVals = m_schemeData->GetNmultiStepValues();
        int size           = m_schemeData->m_numsteps;
        for (int i = (nMultiStepVals - 1); i > 0; i--)
        {
            m_solVector[i] = m_solVector[i - 1];
        }

        for (int i = (size - 1); i > nMultiStepVals; i--)
        {
            m_solVector[i] = m_solVector[i - 1];
        }
    }

private:
    const TimeIntegrationSchemeData *m_schemeData;
    TripleArray m_solVector;
    Array<OneD, NekDouble> m_t;

}; // end class TimeIntegrationSolution

} // end of namespace LibUtilities
} // end of namespace Nektar
