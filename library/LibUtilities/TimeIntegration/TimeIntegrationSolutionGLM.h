///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationSolutionGLM.h
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

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SOLUTION_GLM
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_SOLUTION_GLM

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationAlgorithmGLM.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationTypes.hpp>

namespace Nektar
{
namespace LibUtilities
{

class TimeIntegrationSolutionGLM
{
public:
    // Constructor for single step methods
    LUE TimeIntegrationSolutionGLM(const TimeIntegrationAlgorithmGLM *schemeAlgorithm,
                                   const DoubleArray &y, const NekDouble time,
                                   const NekDouble timestep);

    // Constructor for multi-step methods
    LUE TimeIntegrationSolutionGLM(const TimeIntegrationAlgorithmGLM *schemeAlgorithm,
                                   const TripleArray &y,
                                   const Array<OneD, NekDouble> &t);

    LUE TimeIntegrationSolutionGLM(const TimeIntegrationAlgorithmGLM *schemeAlgorithm,
                                   const unsigned int nvar,
                                   const unsigned int npoints);

    LUE TimeIntegrationSolutionGLM(const TimeIntegrationAlgorithmGLM *schemeAlgorithm);

    inline const TimeIntegrationAlgorithmGLM *GetIntegrationSchemeData() const
    {
        return m_schemeAlgorithm;
    }

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

    // Sets the solution Vector
    inline void SetSolutionVector(const int Offset, const DoubleArray &y)
    {
        m_solVector[Offset] = y;
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
        return m_schemeAlgorithm->m_numsteps;
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
        return m_schemeAlgorithm->GetNmultiStepValues();
    }

    // Return the number of entries in the solution vector that correspond to
    // (multi-step) derivatives.
    inline unsigned int GetNderivs() const
    {
        return m_schemeAlgorithm->GetNmultiStepDerivs();
    }

    // Returns an array which indicates to which time-level the entries in the
    // solution vector correspond.
    inline const Array<OneD, const unsigned int> &GetTimeLevelOffset()
    {
        return m_schemeAlgorithm->GetTimeLevelOffset();
    }

    // Returns the entry in the solution vector which corresponds to the
    // (multi-step) value at the time-level with specified offset
    inline DoubleArray &GetValue(const unsigned int timeLevelOffset)
    {
        int nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        const Array<OneD, const unsigned int> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

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
        int nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        int size           = m_schemeAlgorithm->m_numsteps;
        const Array<OneD, const unsigned int> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

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
        int nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        const Array<OneD, const unsigned int> &offsetvec =
             m_schemeAlgorithm->GetTimeLevelOffset();

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
        int nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        const Array<OneD, const unsigned int> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

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
        int nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        int size           = m_schemeAlgorithm->m_numsteps;
        const Array<OneD, const unsigned int> &offsetvec =
            m_schemeAlgorithm->GetTimeLevelOffset();

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

    // Rotate the solution vector
    // (i.e. updating without calculating/inserting new values)
    inline void RotateSolutionVector()
    {
        int nMultiStepVals = m_schemeAlgorithm->GetNmultiStepValues();
        int size           = m_schemeAlgorithm->m_numsteps;
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
    const TimeIntegrationAlgorithmGLM *m_schemeAlgorithm;

    TripleArray m_solVector;
    Array<OneD, NekDouble> m_t;

}; // end class TimeIntegrationSolutionGLM

} // end of namespace LibUtilities
} // end of namespace Nektar

#endif
