///////////////////////////////////////////////////////////////////////////////
//
// File FilterAeroForcesSPM.h
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
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INCNAVIERSTOKES_FILTERS_FILTERAEROFORCESSPM_H
#define NEKTAR_INCNAVIERSTOKES_FILTERS_FILTERAEROFORCESSPM_H

#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
class FilterAeroForcesSPM : public SolverUtils::Filter
{
public:
    /// Creates an instance of this class
    static SolverUtils::FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr       &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const std::map<std::string, std::string>         &pParams)
    {
        SolverUtils::FilterSharedPtr p = MemoryManager<FilterAeroForcesSPM>::
                            AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    ///Name of the class
    static std::string className;

    FilterAeroForcesSPM(
        const LibUtilities::SessionReaderSharedPtr       &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const std::map<std::string, std::string>         &pParams);

    virtual ~FilterAeroForcesSPM();

    // Calculates the forces and fills the array 'm_Forces' up
    void CalculateForces(
        const Array<OneD, Array<OneD, NekDouble> > &pIntVel,
        const Array<OneD, Array<OneD, NekDouble> > &pUpPrev,
        const MultiRegions::ExpListSharedPtr &pPhi,
        NekDouble time,
        NekDouble dt);

protected:
    unsigned int m_index;
    unsigned int m_outputFrequency;
    std::string m_outputFile;
    std::ofstream m_outputStream;
    // Time when we start calculating the forces
    NekDouble m_startTime;
    /// STL vector containing the names of the different directions
    std::vector<std::string> m_dirNames;
    /// Dimension of the fluid domain
    NekDouble m_spaceDim;
    /// Array storing the last value of the aerodynamic forces
    Array<OneD, NekDouble> m_Forces;

    virtual void v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);

    virtual void v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);

    virtual void v_Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);

    virtual bool v_IsTimeDependent();

private:
};

typedef std::shared_ptr<FilterAeroForcesSPM> FilterAeroForcesSPMSharedPtr;
}

#endif /* NEKTAR_INCNAVIERSTOKES_FILTERS_FILTERAEROFORCESSPM_H */
