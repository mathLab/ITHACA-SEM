///////////////////////////////////////////////////////////////////////////////
//
// File: Extrapolate.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_EXTRAPOLATE_H
#define NEKTAR_SOLVERS_EXTRAPOLATE_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <SolverUtils/AdvectionSystem.h>


namespace Nektar
{

    enum HBCType
    {
        eNOHBC,
        eHBCNeumann,    // Standard High Order BC
        eOBC,           // High Order outflow BC (Neumann-Dirichlet) from Dong et al JCP 2014
        eConvectiveOBC  // Convective High Order (Robin type) BC from Dong JCP 2015
    };
    
    // Forward declaration
    class Extrapolate;
    typedef std::shared_ptr<Extrapolate> ExtrapolateSharedPtr;
    typedef LibUtilities::NekFactory< std::string, Extrapolate,
        const LibUtilities::SessionReaderSharedPtr& ,
        Array<OneD, MultiRegions::ExpListSharedPtr>& ,
        MultiRegions::ExpListSharedPtr& ,
        const Array<OneD, int>& ,
        const SolverUtils::AdvectionSharedPtr& > ExtrapolateFactory;

    struct HighOrderOutflow;
    typedef std::shared_ptr<HighOrderOutflow> HighOrderOutflowSharedPtr;


    ExtrapolateFactory& GetExtrapolateFactory();

    class Extrapolate
    {
    public:
        Extrapolate(        
            const LibUtilities::SessionReaderSharedPtr  pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
            MultiRegions::ExpListSharedPtr              pPressure,
            const Array<OneD, int>                      pVel,
            const SolverUtils::AdvectionSharedPtr                advObject);
        
        virtual ~Extrapolate();

        void GenerateHOPBCMap(const LibUtilities::SessionReaderSharedPtr& pSsession);
        
        void UpdateRobinPrimCoeff(void);

        inline void SubSteppingTimeIntegration(
            const LibUtilities::TimeIntegrationSchemeSharedPtr &IntegrationScheme);

        inline void SubStepSaveFields(
            const int nstep);
        
        inline void SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const NekDouble Aii_DT,
            NekDouble kinvis);

        inline void SubStepAdvance(
            const LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr &integrationSoln, 
            const int nstep, 
            NekDouble time);

        inline void MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, 
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Advection);

        inline void  EvaluatePressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis);

        
        void AddDuDt(void);

        void AddVelBC(void);

        void ExtrapolatePressureHBCs(void);
        void CopyPressureHBCsToPbndExp(void);

        Array<OneD,NekDouble> GetMaxStdVelocity(
            const Array<OneD, Array<OneD,NekDouble> > inarray);
        

        void CorrectPressureBCs( const Array<OneD, NekDouble>  &pressure);
        
        void IProductNormVelocityOnHBC(const Array<OneD, const Array<OneD, NekDouble> >  &Vel, 
                                       Array<OneD, NekDouble> &IprodVn);
        
        void IProductNormVelocityBCOnHBC(Array<OneD, NekDouble> &IprodVn);
        
        std::string GetSubStepName(void); 

        void ExtrapolateArray( Array<OneD, Array<OneD, NekDouble> > &array);

        void EvaluateBDFArray(Array<OneD, Array<OneD, NekDouble> > &array);

        void AccelerationBDF( Array<OneD, Array<OneD, NekDouble> > &array);

        void ExtrapolateArray(
            Array<OneD, Array<OneD, NekDouble> > &oldarrays,
            Array<OneD, NekDouble>  &newarray,
            Array<OneD, NekDouble>  &outarray);
        
        void AddNormVelOnOBC(const int nbcoeffs, const int nreg,
                             Array<OneD, Array<OneD, NekDouble> > &u);

        void AddPressureToOutflowBCs(NekDouble kinvis);
        
    protected: 
        virtual void v_EvaluatePressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis)=0;

       virtual void v_SubSteppingTimeIntegration(
            const LibUtilities::TimeIntegrationSchemeSharedPtr & IntegrationScheme ) = 0;

        virtual void v_SubStepSaveFields(
            int nstep)=0;

        virtual void v_SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            NekDouble Aii_DT,
            NekDouble kinvis)=0;

        virtual void v_SubStepAdvance(
            const LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr & integrationSoln, 
                  int                                                                     nstep, 
                  NekDouble                                                               time ) = 0;

        virtual void v_MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, 
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Advection)=0;
        
        virtual std::string v_GetSubStepName(void);

        void CalcNeumannPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis)
        {
            v_CalcNeumannPressureBCs( fields, N, kinvis);
        }
                
        virtual void v_CalcNeumannPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis);
            
        virtual void v_CorrectPressureBCs( const Array<OneD, NekDouble>  &pressure);

        virtual void v_AddNormVelOnOBC(const int nbcoeffs, const int nreg,
                                       Array<OneD, Array<OneD, NekDouble> > &u);
        void CalcOutflowBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            NekDouble kinvis);

        void RollOver(Array<OneD, Array<OneD, NekDouble> > &input);

        LibUtilities::SessionReaderSharedPtr m_session;

        LibUtilities::CommSharedPtr m_comm;

        /// Array of type of high order BCs for splitting shemes
        Array<OneD, HBCType> m_hbcType; 
        
        /// Velocity fields 
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;

        /// Pointer to field holding pressure field
        MultiRegions::ExpListSharedPtr m_pressure;

        /// int which identifies which components of m_fields contains the
        /// velocity (u,v,w);
        Array<OneD, int> m_velocity;

        SolverUtils::AdvectionSharedPtr m_advObject;

        Array<OneD, Array<OneD, NekDouble> > m_previousVelFields;

        /// Curl-curl dimensionality
        int m_curl_dim;
        
        /// bounday dimensionality
        int m_bnd_dim;

        /// pressure boundary conditions container
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr> m_PBndConds;

        /// pressure boundary conditions expansion container
        Array<OneD, MultiRegions::ExpListSharedPtr>  m_PBndExp;
        
        /// number of times the high-order pressure BCs have been called
        int m_pressureCalls;

        // Number of degrees of freedom (coefficients) for HOPbc
        int m_numHBCDof;

        // Number of HOPbcs
        int m_HBCnumber;

        /// Maximum points used in pressure BC evaluation
        int m_intSteps;

        NekDouble m_timestep;

        /// Storage for current and previous levels of high order pressure boundary conditions.
        Array<OneD, Array<OneD, NekDouble> > m_pressureHBCs;

        /// Storage for current and previous levels of the inner product of normal velocity
        Array<OneD, Array<OneD, NekDouble> > m_iprodnormvel;

        Array<OneD, Array<OneD, NekDouble> > m_traceNormals;

        // Velocity correction scheme coefficient required for extrapolation.
        static NekDouble StifflyStable_Betaq_Coeffs[3][3];
        static NekDouble StifflyStable_Alpha_Coeffs[3][3];
        static NekDouble StifflyStable_Gamma0_Coeffs[3];

        // data related to high order outflow. 
        HighOrderOutflowSharedPtr m_houtflow; 
        
    private:
        static std::string def;
    };
    
    
    struct HighOrderOutflow
    {
    HighOrderOutflow(const int numHOpts, const int outHBCnumber,const int curldim, const LibUtilities::SessionReaderSharedPtr &pSession):
        m_numOutHBCPts(numHOpts),
            m_outHBCnumber(outHBCnumber)
        {
            m_outflowVel = Array<OneD,
                Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (outHBCnumber);
            
            m_outflowVelBnd = Array<OneD,
                Array<OneD, Array<OneD,
                Array<OneD, NekDouble> > > > (outHBCnumber);
            
            m_UBndExp  = Array<OneD,
                Array<OneD, MultiRegions::ExpListSharedPtr> >(curldim);
            
            pSession->LoadParameter("OutflowBC_theta",  m_obcTheta,  1.0);
            pSession->LoadParameter("OutflowBC_alpha1", m_obcAlpha1, 0.0);
            pSession->LoadParameter("OutflowBC_alpha2", m_obcAlpha2, 0.0);
            
            pSession->LoadParameter("U0_HighOrderBC",    m_U0,1.0);
            pSession->LoadParameter("Delta_HighOrderBC", m_delta,1/20.0);
        }

        virtual ~HighOrderOutflow()
        {};
        
        /// Number of quadrature points for Outflow HOBC
        int m_numOutHBCPts;

        /// Number of Outflow HOBCs
        int m_outHBCnumber;

        /// Parameters for outflow boundary condition
        NekDouble   m_obcTheta;
        NekDouble   m_obcAlpha1;
        NekDouble   m_obcAlpha2;
        NekDouble   m_U0;
        NekDouble   m_delta;
        std::string m_defVelPrimCoeff[3]; 
        
        /// Storage for current and previous velocity fields along the outflow 
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble > > > > m_outflowVel;

        /// Storage for current and previous velocities along the outflow boundary
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble > > > > m_outflowVelBnd;
        
        /// Velocity boundary condition expansions on high order boundaries. 
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >  m_UBndExp;

        /// primitive coefficient for pressure  when using convetive like OBCs
        Array<OneD, NekDouble> m_pressurePrimCoeff; 

        /// primitive coefficient for velocities when using convetive like OBCs
        Array<OneD, Array<OneD, NekDouble> > m_velocityPrimCoeff; 
    };

    /**
     * Evaluate Pressure Boundary Conditions for Standard Extrapolation
     */
    inline void Extrapolate::EvaluatePressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis)
    {
        v_EvaluatePressureBCs(inarray,N,kinvis);
    }

    /**
     *
     */
    inline void Extrapolate::SubSteppingTimeIntegration(
        const LibUtilities::TimeIntegrationSchemeSharedPtr &IntegrationScheme)
    {
        v_SubSteppingTimeIntegration(IntegrationScheme);
    }
    
    /**
     *
     */
    inline void Extrapolate::SubStepSaveFields(
        int nstep)
    {
        v_SubStepSaveFields(nstep);
    }

    /**
     *
     */
    inline void Extrapolate::SubStepSetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        NekDouble Aii_DT,
        NekDouble kinvis)
    {
        v_SubStepSetPressureBCs(inarray,Aii_DT,kinvis);
    }
    
    /**
     *
     */
    inline void Extrapolate::SubStepAdvance(
        const LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr &integrationSoln, 
        int nstep, 
        NekDouble time)
    {
        v_SubStepAdvance(integrationSoln,nstep, time);
    }

    /**
     *
     */
    inline void Extrapolate::MountHOPBCs(
        int HBCdata, 
        NekDouble kinvis, 
        Array<OneD, NekDouble> &Q, 
        Array<OneD, const NekDouble> &Advection)
    {
        v_MountHOPBCs(HBCdata,kinvis,Q,Advection);
    }

    /**
     *
     */
    inline std::string Extrapolate::GetSubStepName(void)
    {
        return v_GetSubStepName();
    }

    /**
     *
     */
    inline void Extrapolate::CorrectPressureBCs(
        const Array<OneD, NekDouble>  &pressure)
    {
        v_CorrectPressureBCs(pressure);
    }

    
    /**
     *
     */
    inline void Extrapolate::AddNormVelOnOBC(const int nbcoeffs, const int nreg,
                                             Array<OneD, Array<OneD, NekDouble> > &u)
    {
        v_AddNormVelOnOBC(nbcoeffs,nreg,u);
    }

}

#endif

