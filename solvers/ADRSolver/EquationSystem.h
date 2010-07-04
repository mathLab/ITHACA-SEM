#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEM_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEM_H

#include <ADRSolver/Factory.hpp>
#include <ADRSolver/SessionReader.h>
#include <Auxiliary/ADRBase.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SpatialDomains/SpatialData.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    class EquationSystem;

    typedef boost::shared_ptr<EquationSystem> EquationSystemSharedPtr;
    typedef Factory< std::string, EquationSystem, SessionReaderSharedPtr&,
        LibUtilities::TimeIntegrationSchemeOperators&> EquationSystemFactory;

    class EquationSystem : public ADRBase
    {
    public:
        virtual ~EquationSystem();

        void doInitialise();
        void printSummary(std::ostream &out);

        void doOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);

        void doImplicitSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD, Array<OneD, NekDouble> >&outarray,
              NekDouble time,
              NekDouble lambda);

        void setPhysForcingFunction();
        void doSolveHelmholtz();

        void GeneralTimeIntegration(int nsteps,
                           LibUtilities::TimeIntegrationMethod IntMethod,
                   LibUtilities::TimeIntegrationSchemeOperators ode);

        bool isSteady();
        bool hasExplicitAdvection();
        bool hasExplicitDiffusion();
        bool hasExplicitReaction();
        LibUtilities::TimeIntegrationMethod getTimeIntMethod();

    protected:
        EquationSystem(SessionReaderSharedPtr& pSession,
                LibUtilities::TimeIntegrationSchemeOperators& pOde);

        void evaluateFunction(Array<OneD, NekDouble>& pArray,
                SpatialDomains::ConstUserDefinedEqnShPtr pEqn);
        void setBoundaryConditions(NekDouble time);

        virtual void v_printSummary(std::ostream &out);

        virtual void v_doOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);
        virtual void v_doImplicitSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD, Array<OneD, NekDouble> >&outarray,
              NekDouble time,
              NekDouble lambda);
        virtual void v_doSolveHelmholtz();
        virtual bool v_isSteady();

        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxY );

        virtual void v_NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield,
                                        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux);

        virtual void v_NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                    Array<OneD, Array<OneD, NekDouble> >  &qflux);

        SessionReaderSharedPtr                  mSession;
        LibUtilities::TimeIntegrationMethod     mTimeIntMethod;
        bool                                    mExplicitAdvection;
        bool                                    mExplicitDiffusion;
        bool                                    mExplicitReaction;
        NekDouble                               mEpsilon;

    private:
        /// dump info to stdout at steps time
        int                                     mInfoSteps;

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
    };

}

#endif
