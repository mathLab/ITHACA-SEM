///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokes.h
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
// Description: Basic Advection Diffusion Reaction Field definition in two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_INCNAVIERSTOKES_H
#define NEKTAR_SOLVERS_INCNAVIERSTOKES_H

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <SolverUtils/UnsteadySystem.h>
#include <IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

namespace Nektar
{     

    enum EquationType
    {
        eNoEquationType,
        eSteadyStokes,
        eSteadyOseen,
        eSteadyLinearisedNS,
        eUnsteadyStokes,
        eUnsteadyLinearisedNS,
        eUnsteadyNavierStokes,
        eSteadyNavierStokes,
        eEquationTypeSize
    };

    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] =
    {
        "NoType",
        "SteadyStokes",
        "SteadyOseen",
        "SteadyLinearisedNS",
        "UnsteadyStokes",
        "UnsteadyLinearisedNS",
        "UnsteadyNavierStokes",
        "SteadyNavierStokes",
    };


    enum AdvectionForm
    {
        eNoAdvectionForm,
        eConvective,
        eNonConservative,
        eLinearised,
        eAdjoint,
        eSkewSymmetric,
        eNoAdvection,
        eAdvectionFormSize
    };

    // Keep this consistent with the enums in EquationType.
    const std::string kAdvectionFormStr[] =
    {
        "NoType",
        "Convective",
        "NonConservative",
        "Linearised",
        "Adjoint",
        "SkewSymmetric"
        "NoAdvection"
    };

    /**
     * \brief This class is the base class for Navier Stokes problems
     *
     */
    class IncNavierStokes: public SolverUtils::UnsteadySystem
    {
    public:
        // Destructor
        virtual ~IncNavierStokes();

        virtual void v_InitObject();


        virtual void v_GetFluxVector(
                const int i,
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numflux);

        AdvectionTermSharedPtr GetAdvObject(void)
        {
            return m_advObject;
        }


        int GetNConvectiveFields(void)
        {
            return m_nConvectiveFields;  
        }

        Array<OneD, int> &GetVelocity(void)
        {
            return  m_velocity; 
        }

        NekDouble GetSubstepTimeStep();
        
        Array<OneD, NekDouble> GetElmtCFLVals(void);
        
        NekDouble GetCFLEstimate(int &elmtid);

        // Mapping of the real convective field on the standard element.
        // This function gives back the convective filed in the standard
        // element to calculate the stability region of the problem in a
        // unique way.
        Array<OneD,NekDouble> GetMaxStdVelocity(
                const Array<OneD, Array<OneD,NekDouble> > inarray);


        // Sub-stepping related methods
        void SubStepAdvection (
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                      Array<OneD, Array<OneD, NekDouble> > &outarray,
                const NekDouble time);

        void SubStepProjection(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                      Array<OneD, Array<OneD, NekDouble> > &outarray,
                const NekDouble time);

        void SubStepExtrapoloteField(
                NekDouble toff,
                Array< OneD, Array<OneD, NekDouble> > &ExtVel);

        void AddAdvectionPenaltyFlux(
                const Array<OneD, const Array<OneD, NekDouble> > &velfield,
                const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                      Array<OneD, Array<OneD, NekDouble> > &outarray);

    protected:
        /// modal energy file
        std::ofstream m_mdlFile;

        LibUtilities::TimeIntegrationSolutionSharedPtr  m_integrationSoln;

        /// bool to identify if using a substepping scheme
        bool m_subSteppingScheme;
        /// bool to identify if advection term smoothing is requested
        bool m_SmoothAdvection;

        LibUtilities::TimeIntegrationWrapperSharedPtr m_subStepIntegrationScheme;
        //LibUtilities::TimeIntegrationSchemeSharedPtr m_subStepIntegrationScheme;
        LibUtilities::TimeIntegrationSchemeOperators m_subStepIntegrationOps;

        Array<OneD, Array<OneD, NekDouble> > m_previousVelFields;

        /// Advection term
        AdvectionTermSharedPtr m_advObject;

        /// Number of fields to be convected;
        int   m_nConvectiveFields;

        /// int which identifies which components of m_fields contains the
        /// velocity (u,v,w);
        Array<OneD, int> m_velocity;

        /// Pointer to field holding pressure field
        MultiRegions::ExpListSharedPtr m_pressure;

        /// Kinematic viscosity
        NekDouble   m_kinvis;
        /// dump energy to file at steps time
        int         m_energysteps;
        /// dump cfl estimate
        int         m_cflsteps;
        /// Check for steady state at step interval
        int         m_steadyStateSteps;
        /// Tolerance to which steady state should be evaluated at
        NekDouble   m_steadyStateTol;

        /// equation type;
        EquationType  m_equationType;

        /// Mapping from BCs to Elmt IDs
        Array<OneD, Array<OneD, int> > m_fieldsBCToElmtID;
        /// Mapping from BCs to Elmt Edge IDs
        Array<OneD, Array<OneD, int> > m_fieldsBCToTraceID;
        /// RHS Factor for Radiation Condition
        Array<OneD, Array<OneD, NekDouble> > m_fieldsRadiationFactor;

        /// Time integration classes
        LibUtilities::TimeIntegrationSchemeOperators m_integrationOps;
        LibUtilities::TimeIntegrationWrapperSharedPtr m_integrationScheme;

        /// Number of time integration steps AND Order of extrapolation for
        /// pressure boundary conditions.
        int m_intSteps;

        /// Constructor.
        IncNavierStokes(const LibUtilities::SessionReaderSharedPtr& pSession);

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }

        void AdvanceInTime(int nsteps);

        void SubStepAdvance   (const int nstep);
        void SubStepSaveFields(const int nstep);

        void EvaluateAdvectionTerms(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                      Array<OneD, Array<OneD, NekDouble> > &outarray,
                      Array<OneD, NekDouble> &wk = NullNekDouble1DArray);

        void WriteModalEnergy(void);

        /// time dependent boundary conditions updating
        void SetBoundaryConditions(NekDouble time);

        /// Set Radiation forcing term
        void SetRadiationBoundaryForcing(int fieldid);

        /// evaluate steady state
        bool CalcSteadyState(void);

        virtual MultiRegions::ExpListSharedPtr v_GetPressure()
        {
            return m_pressure;
        }

        virtual void v_PrintSummary(std::ostream &out)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

        virtual void v_DoInitialise(void)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

        virtual void v_DoSolve(void)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

        virtual void v_TransCoeffToPhys(void)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

        virtual void v_TransPhysToCoeff(void)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

    private:

    };

    typedef boost::shared_ptr<IncNavierStokes> IncNavierStokesSharedPtr;

} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H
