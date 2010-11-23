#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAIN_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MONODOMAIN_H

#include <ADRSolver/EquationSystems/UnsteadySystem.h>

namespace Nektar
{
    // Forward declaration
    class CellModel;

    /// A shared pointer to an EquationSystem object
    typedef boost::shared_ptr<CellModel> CellModelSharedPtr;
    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the EquationSystem class.
    typedef LibUtilities::NekFactory< std::string, CellModel,
        SessionReaderSharedPtr&, int> CellModelFactory;

    /// A model for cardiac conduction.
    class Monodomain : public UnsteadySystem
    {
    public:
        /// Creates an instance of this class
        static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession)
        {
            return MemoryManager<Monodomain>::AllocateSharedPtr(pSession);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        Monodomain(SessionReaderSharedPtr& pSession);

        /// Desctructor
        virtual ~Monodomain();

    protected:
        /// Solve for the diffusion term.
        void DoImplicitSolve(
                const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                      NekDouble time,
                      NekDouble lambda);

        /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
        void DoOdeRhs(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        /// Sets a custom initial condition.
        virtual void v_SetInitialConditions(NekDouble initialtime,
                                bool dumpInitialConditions);

        /// Prints a summary of the model parameters.
        virtual void v_PrintSummary(std::ostream &out);

    private:
        /// Cell model.
        CellModelSharedPtr m_cell;
    };



    /// Cell model base class.
    class CellModel
    {
    public:
        CellModel(SessionReaderSharedPtr& pSession, const int nq);
        virtual ~CellModel() {}

        virtual void Update(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time) = 0;

        virtual void v_PrintSummary(std::ostream &out) = 0;

    protected:
        /// Spatially varying parameters.
        SpatialDomains::SpatialParametersSharedPtr  m_spatialParameters;
        /// Number of physical points.
        int m_nq;
    };


    /// Aliev Panfilov model.
    class CellModelAlievPanfilov : public CellModel
    {
    public:
        /// Creates an instance of this class
        static CellModelSharedPtr create(SessionReaderSharedPtr& pSession, const int nq)
        {
            return MemoryManager<CellModelAlievPanfilov>::AllocateSharedPtr(pSession, nq);
        }

        /// Name of class
        static std::string className;

        CellModelAlievPanfilov(SessionReaderSharedPtr& pSession, const int nq);
        virtual ~CellModelAlievPanfilov() {}

        virtual void Update(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        virtual void v_PrintSummary(std::ostream &out);

    private:
        /// Trigger parameter a.
        NekDouble m_a;
        /// Scaling parameter k.
        NekDouble m_k;
        /// Restitution parameter \f$\mu_1\f$.
        NekDouble m_mu1;
        /// Restitution parameter \f$\mu_2\f$.
        NekDouble m_mu2;
        /// Restitution parameter \f$\epsilon\f$.
        NekDouble m_eps;

        /// Temporary space for storing \f$u^2\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uu;
        /// Temporary space for storing \f$u^3\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uuu;
        /// Workspace for computing reaction term.
        Array<OneD, NekDouble> m_tmp1;
        /// Workspace for computing reaction term.
        Array<OneD, NekDouble> m_tmp2;
    };


    /// FitzHugh-Nagumo model.
    class CellModelFitzHughNagumo : public CellModel
    {
    public:
        /// Creates an instance of this class
        static CellModelSharedPtr create(SessionReaderSharedPtr& pSession, const int nq)
        {
            return MemoryManager<CellModelFitzHughNagumo>::AllocateSharedPtr(pSession, nq);
        }

        /// Name of class
        static std::string className;

        CellModelFitzHughNagumo(SessionReaderSharedPtr& pSession, const int nq);
        virtual ~CellModelFitzHughNagumo() {}

        virtual void Update(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        virtual void v_PrintSummary(std::ostream &out);

    private:
        NekDouble              m_beta;
        NekDouble              m_epsilon;
        /// Temporary space for storing \f$u^3\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uuu;
    };
}

#endif
