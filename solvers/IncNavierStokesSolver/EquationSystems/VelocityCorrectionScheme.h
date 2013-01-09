///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionScheme.h
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
// Description: Velocity Correction Scheme header 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEME_H
#define NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEME_H

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

namespace Nektar
{     
    static NekDouble StifflyStable_Betaq_Coeffs[][3] = {{ 1.0,  0.0, 0.0},
                                                        { 2.0, -1.0, 0.0},
                                                        { 3.0, -3.0, 1.0}};

    static NekDouble StifflyStable_Alpha_Coeffs[][3] = {{ 1.0,  0.0, 0.0},
                                                        { 2.0, -0.5, 0.0},
                                                        { 3.0, -1.5, 1.0/3.0}};

    static NekDouble StifflyStable_Gamma0_Coeffs[3]  = {1.0,  1.5, 11.0/6.0};

    /**
     * \brief This class is the base class for the Velocity Correction Scheme
     *
     */
    
    struct HBCInfo
    {
        int m_globalElmtID;  // elements ID in the global ordering
        int m_ptsInElmt;     // number of physical points of the element
        int m_physOffset;    // elmenent physical offset in the global list
        int m_bndElmtOffset; // element offset in the boundary expansion
        int m_elmtTraceID;   // trace ID on the element
        int m_bndryElmtID;   // pressure boundary condition ID
        int m_assPhysOffset; // associated elments physical offset (k and k_c are the real and the complex plane)
        int m_coeffOffset;   // coefficients offset used to locate the acceleration term in the general m_pressureHBC
    };
    
    class VelocityCorrectionScheme: public IncNavierStokes
    {
    public:           
        
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession) {
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<VelocityCorrectionScheme>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
            }
            /// Name of class
            static std::string className;
            
            
            /**
             * Constructor.
             * \param 
             * 
             */
            VelocityCorrectionScheme(const LibUtilities::SessionReaderSharedPtr& pSession);
            
            virtual ~VelocityCorrectionScheme();
            
            virtual void v_InitObject();
            
            void SubStepSetPressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &inarray, const NekDouble Aii_DT);
            
            void EvaluatePressureBCs(const Array<OneD, const Array< OneD,  NekDouble> > &fields, const Array<OneD, const Array< OneD,  NekDouble> > &N, const NekDouble Aii_Dt = NekConstants::kNekUnsetDouble);
            
            
            void SetUpPressureForcing(const Array<OneD, const Array<OneD, NekDouble> > &fields, Array<OneD, Array<OneD, NekDouble> > &Forcing, const NekDouble aii_Dt);
            void SetUpViscousForcing(const Array<OneD, const Array<OneD, NekDouble> > &inarray, Array<OneD, Array<OneD, NekDouble> > &Forcing, const NekDouble aii_Dt);
            
            void SolveUnsteadyStokesSystem(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                           Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                           const NekDouble time,
                                           const NekDouble a_iixDt);
            
            void EvaluateAdvection_SetPressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                  Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                  const NekDouble time);
            
            void AddDuDt  (const Array<OneD, const Array<OneD, NekDouble> >  &N, NekDouble Aii_Dt);
            void AddDuDt2D(const Array<OneD, const Array<OneD, NekDouble> >  &N, NekDouble Aii_Dt);
            void AddDuDt3D(const Array<OneD, const Array<OneD, NekDouble> >  &N, NekDouble Aii_Dt);
            
    protected:
        
    private: 
        int m_pressureCalls;
        int m_pressureBCsMaxPts;                // Maximum points used in pressure BC evaluation 

        bool m_showTimings;                     // Show timings for each step
        bool m_useHomo1DSpecVanVisc;                  // bool to identify if spectral vanishing viscosity is active. 
        Array<OneD, int> m_pressureBCtoElmtID;  // Id of element to which pressure  boundary condition belongs
        Array<OneD, int> m_pressureBCtoTraceID; // Id of edge (2D) or face (3D) to which pressure boundary condition belongs
        
        Array<OneD, Array<OneD, NekDouble> >  m_pressureHBCs; //< Storage for current and previous levels of high order pressure boundary conditions. 
        Array<OneD, Array<OneD, NekDouble> >  m_acceleration;
        
        Array<OneD, HBCInfo > m_HBCdata;  //data structure to old all the information regarding High order pressure BCs
        
        StdRegions::StdExpansionSharedPtr m_elmt; // general standard element used to deaal with HOPBC calculations
        
        Array<OneD, NekDouble>  m_wavenumber;            // wave number 2 pi k /Lz
        Array<OneD, NekDouble>  m_negWavenumberSq;      // minus Square of wavenumber
        
        /**  \brief This function evaluates the normal Neumann pressure
         *  boundary condition for the velocity correction scheme at the
         *  current time level which requires as input the non-linear
         *  terms at the current time level
         *
         *   \f[ \frac{\partial p}{\partial n}^n = 
         *      \left [ {\bf N( u)}^n - kinvis \nabla \times \nabla {\bf u}^n \right]
         *                  \cdot {\bf n} \f] 
         *                                                                           
         * where \f$ {\bf n}\f$ is the unit outward normal along the edge,
         * \f$ {\bf u} \f$ is the velocity field, and \f$ {\bf N(u)}\f$
         * are the non-linear terms in the momentum equation.
         */
        
        void CalcPressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N);
        
        void CalcPressureBCs2D(const Array<OneD, const Array<OneD, NekDouble> > &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N);
        
        void CalcPressureBCs3D(const Array<OneD, const Array<OneD, NekDouble> > &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N);
        
        void FillHOPBCMap(const int HOPBCnumber);
		
        void Rotate(Array<OneD, Array<OneD, NekDouble> > &input);
        
        // Virtual functions 
        virtual void v_PrintSummary(std::ostream &out);
        
        virtual void v_DoSolve(void);
        
        virtual void v_TransCoeffToPhys(void);
        
        virtual void v_TransPhysToCoeff(void);
        
        virtual void v_DoInitialise(void);
        
        virtual Array<OneD, bool> v_GetSystemSingularChecks();
        
        virtual int v_GetForceDimension();
    };
    
    
    typedef boost::shared_ptr<VelocityCorrectionScheme> VelocityCorrectionSchemeSharedPtr;
    
} //end of namespace


#endif //VELOCITY_CORRECTION_SCHEME_H

/**
 * $Log: VelocityCorrectionScheme.h,v $
 * Revision 1.1  2009/09/06 22:31:16  sherwin
 * First working version of Navier-Stokes solver and input files
 *
 **/
