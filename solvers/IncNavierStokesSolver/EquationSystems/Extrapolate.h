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
#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include <SolverUtils/AdvectionSystem.h>

namespace Nektar
{
    struct HBCInfo
    {
        /// Global element ID.
        int m_globalElmtID;
        /// Number of physical points in the element.
        int m_ptsInElmt;
        /// Physical offset of the element in the global list.
        int m_physOffset;
        /// Physical offset of the element in the boundary expansion.
        int m_bndElmtOffset;
        /// Trace ID of the element
        int m_elmtTraceID;
        /// Pressure boundary condition ID.
        int m_bndryElmtID;
        /// Associated element physical offset (\f$ k\f$ and \f$ k_c\f$ are the
        /// real and complex plane).
        int m_assPhysOffset;
        /// Coefficient offset used to locate the acceleration term in the
        /// general m_pressureHBC.
        int m_coeffOffset;
    };

    // Forward declaration
    class Extrapolate;
    typedef boost::shared_ptr<Extrapolate> ExtrapolateSharedPtr;
    typedef LibUtilities::NekFactory< std::string, Extrapolate,
        const LibUtilities::SessionReaderSharedPtr& ,
        Array<OneD, MultiRegions::ExpListSharedPtr>& ,
        MultiRegions::ExpListSharedPtr& ,
        const Array<OneD, int>& ,
        const SolverUtils::AdvectionSharedPtr& > ExtrapolateFactory;

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

        void GenerateHOPBCMap();

        inline void SubSteppingTimeIntegration(
            const int intMethod,
            const LibUtilities::TimeIntegrationWrapperSharedPtr &IntegrationScheme);

        inline void SubStepSaveFields(
            const int nstep);
        
        inline void SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            const NekDouble Aii_DT,
            NekDouble kinvis);

        inline void SubStepAdvance(
            const LibUtilities::TimeIntegrationSolutionSharedPtr &integrationSoln, 
            const int nstep, 
            NekDouble time);

        inline void MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, 
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Advection);

        void EvaluatePressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis);

        Array<OneD,NekDouble> GetMaxStdVelocity(
            const Array<OneD, Array<OneD,NekDouble> > inarray);
        
    protected:
        virtual void v_SubSteppingTimeIntegration(
            int intMethod,        
            const LibUtilities::TimeIntegrationWrapperSharedPtr &IntegrationScheme)=0;

        virtual void v_SubStepSaveFields(
            int nstep)=0;

        virtual void v_SubStepSetPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
            NekDouble Aii_DT,
            NekDouble kinvis)=0;

        virtual void v_SubStepAdvance(
            const LibUtilities::TimeIntegrationSolutionSharedPtr &integrationSoln, 
            int nstep, 
            NekDouble time)=0;

        virtual void v_MountHOPBCs(
            int HBCdata, 
            NekDouble kinvis, 
            Array<OneD, NekDouble> &Q, 
            Array<OneD, const NekDouble> &Advection)=0;
        
        void CalcNeumannPressureBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis);
        
        void CalcOutflowBCs(
            const Array<OneD, const Array<OneD, NekDouble> > &fields,
            const Array<OneD, const Array<OneD, NekDouble> >  &N,
            NekDouble kinvis);

        void RollOver(
            Array<OneD, Array<OneD, NekDouble> > &input);

        void CurlCurl(
            Array<OneD, Array<OneD, const NekDouble> > &Vel,
            Array<OneD, Array<OneD, NekDouble> > &Q,
            const int j);
        
        LibUtilities::SessionReaderSharedPtr m_session;

        LibUtilities::CommSharedPtr m_comm;

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
        
        /// Maximum points used in pressure BC evaluation
        int m_pressureBCsMaxPts;

        /// Maximum points used in Element adjacent to pressure BC evaluation
        int m_pressureBCsElmtMaxPts;

        /// Maximum points used in pressure BC evaluation
        int m_intSteps;

        NekDouble m_timestep;

        /// Flag to determine if single homogeneous mode is used.
        bool m_SingleMode;
        /// Flag to determine if half homogeneous mode is used.
        bool m_HalfMode;
        /// Flag to determine if use multiple homogenenous modes are used.
        bool m_MultipleModes;

        NekDouble m_LhomZ;  ///< physical length in Z direction (if homogeneous)
        
        int m_npointsX;     ///< number of points in X direction (if homogeneous)
        int m_npointsY;     ///< number of points in Y direction (if homogeneous)
        int m_npointsZ;     ///< number of points in Z direction (if homogeneous)

        /// Id of element to which pressure  boundary condition belongs
        Array<OneD, int> m_pressureBCtoElmtID;
        
        /// Id of edge (2D) or face (3D) to which pressure boundary condition belongs
        Array<OneD, int> m_pressureBCtoTraceID;
        
        /// Storage for current and previous levels of high order pressure boundary conditions.
        Array<OneD, Array<OneD, NekDouble> >  m_pressureHBCs;

        /// Storage for current and previous levels of the acceleration term.
        Array<OneD, Array<OneD, NekDouble> >  m_acceleration;
        
        /// data structure to old all the information regarding High order pressure BCs
        Array<OneD, HBCInfo > m_HBCdata;

        /// wave number 2 pi k /Lz
        Array<OneD, NekDouble>  m_wavenumber;
        
        /// minus Square of wavenumber
        Array<OneD, NekDouble>  m_negWavenumberSq;
        
        /// Storage for current and previous velocity fields at the otuflow for high order outflow BCs
        Array<OneD, Array<OneD, Array<OneD, NekDouble > > > m_outflowVel;

        /// Storage for current and previous velocity fields in physical space at the otuflow for high order outflow BCs
        Array<OneD, Array<OneD, Array<OneD, NekDouble > > > m_PhyoutfVel; ///(if homogeneous)

        /// Storage for nonlinear term in physical space at the outflow for high order outflow BCs 
        Array<OneD, NekDouble> m_nonlinearterm_phys; ///(if homogeneous)

        ///    Storage for nonlinear term in wave space at the outflow for high order outflow BCs
        Array<OneD, NekDouble> m_nonlinearterm_coeffs; ///(if homogeneous)

        /// expansion sizes of pressure boundary conditions in each plane 
        /// at the outflow for high order outflow BCs
        Array<OneD, unsigned int> m_expsize_per_plane; ///(if homogeneous)

        /// Storage for Fourier Coeffs of Dirichlet pressure condition from the input file
        Array<OneD, NekDouble> m_PBndCoeffs; ///(if homogeneous)

        /// Storage for Fourier Coeffs of Neumann velocity condition from the input file
        Array<OneD, Array<OneD, NekDouble> > m_UBndCoeffs; ///(if homogeneous)

        int m_totexps_per_plane; ///total number of expansion for each plane (if homogeneous)

    private:
        static std::string def;

        // Velocity correction scheme coefficient required for extrapolation.
        static NekDouble StifflyStable_Betaq_Coeffs[3][3];
        static NekDouble StifflyStable_Alpha_Coeffs[3][3];
        static NekDouble StifflyStable_Gamma0_Coeffs[3];
    };

    /**
     *
     */
    inline void Extrapolate::SubSteppingTimeIntegration(
        int intMethod,
        const LibUtilities::TimeIntegrationWrapperSharedPtr &IntegrationScheme)
    {
        v_SubSteppingTimeIntegration(intMethod, IntegrationScheme);
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
        const LibUtilities::TimeIntegrationSolutionSharedPtr &integrationSoln, 
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
}

#endif

