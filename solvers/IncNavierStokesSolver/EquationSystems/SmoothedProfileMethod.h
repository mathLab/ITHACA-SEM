///////////////////////////////////////////////////////////////////////////////
//
// File SmoothedProfileMethod.h
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
// Description: Smoothed Profile Method header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SMOOTHEDPROFILEMETHOD_H
#define NEKTAR_SOLVERS_SMOOTHEDPROFILEMETHOD_H

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>
#include <fstream>

namespace Nektar
{
    class SmoothedProfileMethod: public VelocityCorrectionScheme
    {
    public:
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p =
                MemoryManager<SmoothedProfileMethod>::AllocateSharedPtr(
                    pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

        //Constructor
        SmoothedProfileMethod(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph);

        // Destructor
        virtual ~SmoothedProfileMethod();

        virtual void v_InitObject();

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        // Solves the linear part of the velocity correction scheme incluiding
        // the SPM method calculation for 'fs'
        void SolveUnsteadyStokesSystem(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time,
                    const NekDouble a_iixDt)
        {
            v_SolveUnsteadyStokesSystem(inarray, outarray, time, a_iixDt);
        }

    protected:
        /// Correction pressure field for SPM
        MultiRegions::ExpListSharedPtr m_pressureP;
        /// Velocity of the immersed body(ies)
        Array<OneD, Array<OneD, NekDouble> > m_up;
        Array<OneD, Array<OneD, NekDouble> > m_upPrev;
        /// Vector storing the names of the components of \u_p
        std::vector<std::string> m_velName;
        /// Flag signaling if \f$\u_p\f$ depends on time
        bool m_timeDependentUp;
        /// Stiffly-stable scheme \f$\gamma_0\f$ coefficient
        NekDouble m_gamma0;
        /// Forcing function 'f_s'
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fs;
        /// Shape function 'phi' as expansion list
        MultiRegions::ExpListSharedPtr m_phi;
        /// Function that evaluates the values of \Phi
        SolverUtils::SessionFunctionSharedPtr m_phiEvaluator;
        /// Flag that is true when phi depends on time
        bool m_timeDependentPhi;
        /// Flag indicating that phi was defined in a file
        bool m_filePhi;
        /// Position of "AeroForcesSPM" filter in 'm_session->GetFilters()'
        int m_forcesFilter;

        // Interface for 'v_SolveUnsteadyStokesSystem'
        virtual void v_SolveUnsteadyStokesSystem(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    NekDouble time,
                    NekDouble a_iixDt);
        // Sets the parameters and BCs for the Poisson equation
        void SetUpCorrectionPressure(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing);
        // Solves the Poisson equation for the correction pressure
        void SolveCorrectionPressure(
                    const Array<OneD, NekDouble> &Forcing);
        // Explicitly corrects the velocity by using the force 'fs'
        void SolveCorrectedVelocity(
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    Array<OneD, Array<OneD, NekDouble> > &fields,
                    NekDouble dt);
        // Set proper BCs for the corrected pressure 'p_p'
        void SetCorrectionPressureBCs();
        // Calculates the shape function values
        // (only for non-moving boundaries)
        void UpdatePhiUp(NekDouble time);
        // Calculates the virtual force 'fs'
        void UpdateForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    NekDouble dt);
        // Gets time-dependence information from the first elment of 'name'
        bool GetVarTimeDependence(std::string funcName,
                                  std::string attrName);
        // Returns a Tinyxml handle of the requested function element
        TiXmlElement* GetFunctionHdl(std::string functionName);
        // Reads and set the values of Phi from an analyitic func. or a file
        void ReadPhi();

        /// Initialises the expansions for the intermediate pressure, the 'Phi'
        /// function and the IB forcing
        template <typename T>
        void SetUpExpansions(int nvel)
        {
            int iVel = m_velocity[0];
            m_pressureP = MemoryManager<T>::AllocateSharedPtr(
                          *std::dynamic_pointer_cast<T>(m_pressure));
            m_phi = MemoryManager<T>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<T>(m_fields[iVel]));

            m_fs = Array<OneD, MultiRegions::ExpListSharedPtr>(nvel);
            for (int i = 0; i < nvel; ++i)
            {
                m_fs[i] = MemoryManager<T>::AllocateSharedPtr(
                          *std::dynamic_pointer_cast<T>(m_fields[iVel]));
            }

            // Set to wave space if homogeneous case
            if (m_HomogeneousType != EquationSystem::eNotHomogeneous)
            {
                m_pressureP->SetWaveSpace(true);
                m_phi->SetWaveSpace(true);
                for (int i = 0; i < nvel; ++i)
                {
                    m_fs[i]->SetWaveSpace(true);
                }
            }
        }

    private:
    };

    typedef std::shared_ptr<SmoothedProfileMethod>
            SmoothedProfileMethodSharedPtr;

} // end of namespace

#endif // SMOOTHEDPROFILEMETHOD_H
