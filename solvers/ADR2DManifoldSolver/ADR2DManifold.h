///////////////////////////////////////////////////////////////////////////////
//
// File ADR2DManifold.h
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
// Description: Basic ADR2DManifold definition
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_SOLVERS_ADVECTIONDIFFUSIONREACTION_ADVECTIONDIFFUSIONREACTION_H
#define NEKTAR_SOLVERS_ADVECTIONDIFFUSIONREACTION_ADVECTIONDIFFUSIONREACTION_H

#include <MultiRegions/DisContField2D.h>
#include <SolverUtils/EquationSystem.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{


    enum EquationType
    {
        eNoEquationType,
        eUnsteadyAdvection,
        eUnsteadyDiffusion,
        eUnsteadyDiffusionReaction,
        eNonlinearMorphogenensis,
        eFHNtesttype1,
        eFHNMONO,
        eAlievPanfilov,
        eBarkley,
        eEquationTypeSize,
    };

    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] =
    {
        "NoType",
        "UnsteadyAdvection",
        "UnsteadyDiffusion",
        "UnsteadyDiffusionReaction",
        "NonlinearMorphogenensis",
        "FHNtesttype1",
        "FHNMONO",
        "AlievPanfilov",
        "Barkley"
    };

    /**
     * \brief This class is the base class for the development of solvers.
     *
     * It is basically a class handling vector valued fields where every field is
     * a DisContField2D class
     */

    class ADR2DManifold: public EquationSystem
    {
    public:

        /**
         * Default constructor.
         *
         */
//        ADR2DManifold();

        /**
         * Constructor.
         * /param
         *
         *
         */
        ADR2DManifold(LibUtilities::SessionReaderSharedPtr& pSession);

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }

        bool  GetExplicitAdvection(void)
        {
            return m_explicitAdvection;
        }

        bool  GetExplicitDiffusion(void)
        {
            return m_explicitDiffusion;
        }

        bool  GetExplicitReaction(void)
        {
            return m_explicitReaction;
        }

        LibUtilities::TimeIntegrationMethod GetTimeIntMethod(void)
        {
            return m_timeIntMethod;
        }

        inline int Getinitialwavetype(void)
        {
            return m_initialwavetype;
        }

        inline int GetUsedirDeriv(void)
        {
            return m_UseDirDeriv;
        }

        // Return true if equation is a steady state problem

        void WeakDGAdvectionDir(const Array<OneD, Array<OneD, NekDouble> >& InField,
                                Array<OneD, Array<OneD, NekDouble> >& OutField,
                                bool NumericalFluxIncludesNormal = true,
                                bool InFieldIsInPhysSpace = false, int nvariables = 0);

        void WeakDGAdvectionDir(const Array<OneD, Array<OneD, NekDouble> >& InField,
                                const Array<OneD, Array<OneD, NekDouble> >& dirForcing,
                                Array<OneD, Array<OneD, NekDouble> >& OutField,
                                bool NumericalFluxIncludesNormal = true,
                                bool InFieldIsInPhysSpace = false, int nvariables = 0);

        void WeakDirectionalDerivative(const Array<OneD, NekDouble> &physfield,
                                       const Array<OneD, Array<OneD, NekDouble> > &dirForcing,
                                       Array<OneD, NekDouble> &outarray);

        void MGBoundaryCondtions(int bcRegion, int cnt, NekDouble time);

        void WallBoundary2D(const int bcRegion, const int cnt,
                            const Array<OneD, const Array<OneD, NekDouble> >&inarray);

        void SetUSERDEFINEDInitialConditions(const int Readfntype, bool SetRestingState = true, NekDouble initialtime=0.0);

        void EvaluateUSERDEFINEDExactSolution(unsigned int field, Array<OneD, NekDouble> &outfield,
                                              const NekDouble time, const int Readfntype);

        NekDouble USERDEFINEDError(int field, const int type, const int Readfntype = 0,
                                     Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);

        void SetUpTangentialVectors();

        void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, Array<OneD, Array<OneD, NekDouble> >&flux);

        void GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield,
                                                 Array<OneD, Array<OneD, NekDouble> > &flux);

        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &numflux);

        void NumericalFluxdir(Array<OneD, Array<OneD, NekDouble> > &physfield,
                              const Array<OneD, Array<OneD, NekDouble> > &dirForcing,
                              Array<OneD, Array<OneD, NekDouble> > &numflux);

        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                           Array<OneD, Array<OneD, NekDouble> > &numfluxY);

        void NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield,
                              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux);

        void NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                              Array<OneD, Array<OneD, NekDouble> >  &qflux);

        NekDouble AdvectionSphere(const NekDouble x0j, const NekDouble x1j,
                                  const NekDouble x2j, const NekDouble time);

        NekDouble Morphogenesis(const int field, const NekDouble x0j, const NekDouble x1j,
                                const NekDouble x2j, const NekDouble time);

        void WeakPenaltyforScalar(const int var,
                                  const Array<OneD, const NekDouble> &physfield,
                                  Array<OneD, NekDouble> &penaltyflux,
                                  NekDouble time=0.0);

        void WeakPenaltyforVector(const int var,
                                  const int dir,
                                  const Array<OneD, const NekDouble> &physfield,
                                  Array<OneD, NekDouble> &penaltyflux,
                                  NekDouble C11,
                                  NekDouble time=0.0);

        void ODErhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);

        void ODEeReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                          Array<OneD, Array<OneD, NekDouble> >&outarray,
                          const NekDouble time);

        void ODEeLinearMGReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                  Array<OneD, Array<OneD, NekDouble> >&outarray,
                                  const NekDouble time);

        void ODEeNonlinearMorphoReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                         Array<OneD, Array<OneD, NekDouble> >&outarray,
                                         const NekDouble time);

        void ODEeReactionFHNtest1(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                  Array<OneD, Array<OneD, NekDouble> >&outarray,
                                  const NekDouble time);

        void ODEeReactionFHNmono(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                 Array<OneD, Array<OneD, NekDouble> >&outarray,
                                 const NekDouble time);

        void ODEeReactionAP(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                 Array<OneD, Array<OneD, NekDouble> >&outarray,
                                 const NekDouble time);

        void ODEeReactionBarkley(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                 Array<OneD, Array<OneD, NekDouble> >&outarray,
                                 const NekDouble time);

        void ODEhelmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                          Array<OneD, Array<OneD, NekDouble> >&outarray,
                          NekDouble time,
                          NekDouble lambda);

        void ODEhelmSolveFHNtest1(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                  Array<OneD, Array<OneD, NekDouble> >&outarray,
                                  NekDouble time,
                                  NekDouble lambda);

        void ODEhelmSolveFHNmono(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                 Array<OneD, Array<OneD, NekDouble> >&outarray,
                                 NekDouble time,
                                 NekDouble lambda);

        void ODEhelmSolveAP(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                 Array<OneD, Array<OneD, NekDouble> >&outarray,
                                 NekDouble time,
                                 NekDouble lambda);

        void GeneralTimeIntegration(int nsteps,
                                   LibUtilities::TimeIntegrationMethod IntMethod,
                                   LibUtilities::TimeIntegrationSchemeOperators ode);

        void SolveHelmholtz(const int indx, const NekDouble kappa, const int m_UseDirDeriv = 1);

        void PlotTangentialVectorMap();

        void Getrestingstate(const NekDouble epsilon, const NekDouble beta,
                             Array<OneD, NekDouble> rstates);

        void Summary(std::ostream &out);

        void AdditionalSessionSummary(std::ostream &out);


    protected:
        int m_initialwavetype;                ///< Type of function for initial condition and exact solutions
        NekDouble m_duration;
        int m_UseDirDeriv;
        int m_Anisotropy;
        NekDouble m_Angularfreq;
        NekDouble m_Angleofaxis;
        NekDouble m_x0c, m_x1c, m_x2c;

    private:
        int          m_infosteps;    ///< dump info to stdout at steps time
        EquationType m_equationType; ///< equation type;

        bool m_explicitAdvection;  ///< Flag to identify explicit Advection
        bool m_explicitDiffusion;  ///< Flag to identify explicit Diffusion
        bool m_explicitReaction;   ///< Flag to identify explicit Reaction

        NekDouble m_a, m_b, m_c, m_d;

        // -- AP variables
        NekDouble mK, mA, mB, mC, mEps, mMu1, mMu2;
        NekDouble mDiam;
        // -- end AP variables

        Array<OneD, NekDouble>   m_epsilon;
        NekDouble   m_beta;

        LibUtilities::TimeIntegrationMethod m_timeIntMethod; /// Time integration method

        Array<OneD, Array<OneD, NekDouble> >  m_velocity;
        Array<OneD, Array<OneD, NekDouble> >  m_velmagnitude;
        Array<OneD, Array<OneD, NekDouble> >  m_vellc;
        Array<OneD, NekDouble> m_Vn;

        Array<OneD, Array<OneD, NekDouble> > m_traceNormals_tbasis; // forward normals in tangential basis

        Array<OneD, Array<OneD,NekDouble> > m_dirForcing; // 2 by dim*nq

        void EvaluateAdvectionVelocity();

        void SetBoundaryConditions(const Array<OneD, const Array<OneD, NekDouble> >&inarray, NekDouble time);

        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                                     Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            GetFluxVector(i,physfield,flux);
        }

        virtual void v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield,
                                         Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            GetFluxVector(i,j,physfield,flux);
        }

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                                     Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
          NumericalFlux(physfield, numflux);
        }

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
          NumericalFlux(physfield, numfluxX, numfluxY);
        }

        virtual void v_NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield,
                                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            NumFluxforScalar(ufield, uflux);
        }

        virtual void v_NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                                        Array<OneD, Array<OneD, NekDouble> >  &qflux)
        {
            NumFluxforVector(ufield, qfield, qflux);
        }

        virtual NekDouble v_AdvectionSphere(const NekDouble x0j, const NekDouble x1j,
                                            const NekDouble x2j, const NekDouble time)
        {
            return AdvectionSphere(x0j,x1j,x2j,time);
        }

        virtual NekDouble v_Morphogenesis(const int field, const NekDouble x0j, const NekDouble x1j,
                                          const NekDouble x2j, const NekDouble time)
        {
            return Morphogenesis(field, x0j, x1j, x2j, time);
        }

    };

    typedef boost::shared_ptr<ADR2DManifold> ADR2DManifoldSharedPtr;

} //end of namespace

#endif //NEKTAR_SOLVERS_ADVECTIONDIFFUSIONREACTION_ADVECTIONDIFFUSIONREACTION_H
