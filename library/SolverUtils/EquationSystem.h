///////////////////////////////////////////////////////////////////////////////
//
// File EquationSystem.h
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
// Description: Base class for individual solvers.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_EQUATIONSYSTEM_H
#define NEKTAR_SOLVERUTILS_EQUATIONSYSTEM_H

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <SpatialDomains/SpatialData.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class EquationSystem;
        
        /// A shared pointer to an EquationSystem object
        typedef boost::shared_ptr<EquationSystem> EquationSystemSharedPtr;
        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef LibUtilities::NekFactory<
        std::string, EquationSystem,
        const LibUtilities::SessionReaderSharedPtr&
        > EquationSystemFactory;
        SOLVER_UTILS_EXPORT EquationSystemFactory& GetEquationSystemFactory();
        
        /// A base class for describing how to solve specific equations.
        class EquationSystem
        {
        public:
            /// Destructor
            SOLVER_UTILS_EXPORT virtual ~EquationSystem();
            
            // Set up trace normals if required
            SOLVER_UTILS_EXPORT void SetUpTraceNormals(void);
            
            /// Initialises the members of this object.
            SOLVER_UTILS_EXPORT inline void InitObject();
            
            /// Perform any initialisation necessary before solving the problem.
            SOLVER_UTILS_EXPORT inline void DoInitialise();
            
            /// Solve the problem.
            SOLVER_UTILS_EXPORT inline void DoSolve();
            
            /// Transform from coefficient to physical space.
            SOLVER_UTILS_EXPORT inline void TransCoeffToPhys();
            
            /// Transform from physical to coefficient space.
            SOLVER_UTILS_EXPORT inline void TransPhysToCoeff();
            
            /// Perform output operations after solve.
            SOLVER_UTILS_EXPORT inline void Output();
            
            /// Linf error computation
            SOLVER_UTILS_EXPORT inline NekDouble LinfError(unsigned int field, const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);
            
            /// Get Session name
            SOLVER_UTILS_EXPORT std::string GetSessionName()
            {
                return m_sessionName;
            }
            
            /// Get Session name
            SOLVER_UTILS_EXPORT LibUtilities::SessionReaderSharedPtr GetSession()
            {
                return m_session;
            }
            
            /// Get pressure field if available
            SOLVER_UTILS_EXPORT MultiRegions::ExpListSharedPtr GetPressure(); 
            
            /// Print a summary of parameters and solver characteristics.
            SOLVER_UTILS_EXPORT inline void PrintSummary(std::ostream &out);
            
            /// Set parameter m_lambda
            SOLVER_UTILS_EXPORT inline void SetLambda(NekDouble lambda);
            
            /// Evaluates a function as specified in the session file.
            SOLVER_UTILS_EXPORT void EvaluateFunction(
                Array<OneD, Array<OneD, NekDouble> >& pArray,
                std::string pFunctionName,
                const NekDouble pTime = 0.0);
            
            /// Populate given fields with the function from session.
            SOLVER_UTILS_EXPORT void EvaluateFunction(
                std::vector<std::string> pFieldNames,
                Array<OneD, Array<OneD, NekDouble> > &pFields,
                const std::string& pName);
            
            /// Populate given fields with the function from session.
            SOLVER_UTILS_EXPORT void EvaluateFunction(
                std::vector<std::string> pFieldNames,
                Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const std::string& pName);
            
            // Populate an array with a function variable from session.
            SOLVER_UTILS_EXPORT void EvaluateFunction(
                std::string pFieldName,
                Array<OneD, NekDouble>& pArray,
                const std::string& pFunctionName,
                const NekDouble& pTime = 0.0);
            
            // Describe a function.
            SOLVER_UTILS_EXPORT std::string DescribeFunction(
                std::string pFieldName,
                const std::string &pFunctionName);
            
            /// Perform initialisation of the base flow.
            SOLVER_UTILS_EXPORT void InitialiseBaseFlow(
                Array<OneD, Array<OneD, NekDouble> > &base);
            
            /// Initialise the data in the dependent fields.
            SOLVER_UTILS_EXPORT inline void SetInitialConditions(
                NekDouble initialtime = 0.0,
                bool dumpInitialConditions = true);
            
            /// Evaluates an exact solution
            SOLVER_UTILS_EXPORT inline void EvaluateExactSolution(
                int                     field,
                Array<OneD, NekDouble> &outfield,
                const NekDouble         time);
            
            /// Compute the L2 error between fields and a given exact
            /// solution.
            SOLVER_UTILS_EXPORT NekDouble L2Error(
                unsigned int                 field,
                const Array<OneD,NekDouble> &exactsoln,
                bool                         Normalised = false);
            
            /// Compute the L2 error of the fields
            SOLVER_UTILS_EXPORT inline NekDouble L2Error(
                unsigned int field, 
                bool         Normalised = false)
            {
                return L2Error(field,NullNekDouble1DArray,Normalised);
            }
            
            /// Compute error (L2 and L_inf) over an larger set of quadrature
            /// points return [L2 Linf]
            SOLVER_UTILS_EXPORT Array<OneD,NekDouble> ErrorExtraPoints(
                unsigned int field);
            
            /// Compute the inner product \f$ (\nabla \phi \cdot F) \f$.
            SOLVER_UTILS_EXPORT void WeakAdvectionGreensDivergenceForm(
                const Array<OneD, Array<OneD, NekDouble> > &F,
                Array<OneD,             NekDouble>   &outarray);
            
            /// Compute the inner product \f$ (\phi, \nabla \cdot F) \f$.
            SOLVER_UTILS_EXPORT void WeakAdvectionDivergenceForm(
                const Array<OneD, Array<OneD, NekDouble> > &F,
                Array<OneD,             NekDouble>   &outarray);
            
            /// Compute the inner product \f$ (\phi, V\cdot \nabla u) \f$.
            SOLVER_UTILS_EXPORT void WeakAdvectionNonConservativeForm(
                const Array<OneD, Array<OneD, NekDouble> > &V,
                const Array<OneD,       const NekDouble>   &u,
                Array<OneD,             NekDouble>   &outarray,
                bool UseContCoeffs = false);
            
            /// Compute the non-conservative advection \f$ (V \cdot \nabla u)
            /// \f$.
            SOLVER_UTILS_EXPORT void AdvectionNonConservativeForm(
                const Array<OneD, Array<OneD, NekDouble> > &V,
                const Array<OneD, const NekDouble> &u,
                Array<OneD,       NekDouble> &outarray,
                Array<OneD,       NekDouble> &wk = NullNekDouble1DArray);
            
            /// Calculate the weak discontinuous Galerkin advection.
            SOLVER_UTILS_EXPORT void WeakDGAdvection(
                const Array<OneD, Array<OneD, NekDouble> >& InField,
                Array<OneD, Array<OneD, NekDouble> >& OutField,
                bool NumericalFluxIncludesNormal = true,
                bool InFieldIsInPhysSpace = false,
                int nvariables = 0);
            
            /// Calculate weak DG Diffusion in the LDG form.
            SOLVER_UTILS_EXPORT void WeakDGDiffusion(
                const Array<OneD, Array<OneD, NekDouble> >& InField,
                Array<OneD, Array<OneD, NekDouble> >& OutField,
                bool NumericalFluxIncludesNormal = true,
                bool InFieldIsInPhysSpace = false);
            
            /// Write checkpoint file of #m_fields.
            SOLVER_UTILS_EXPORT void Checkpoint_Output(const int n);
            
            /// Write checkpoint file of custom data fields.
            SOLVER_UTILS_EXPORT void Checkpoint_Output(
                const int n,
                MultiRegions::ExpListSharedPtr &field,
                Array< OneD, Array<OneD, NekDouble> > &fieldcoeffs,
                Array<OneD, std::string> &variables);
            
            /// Write field data to the given filename.
            SOLVER_UTILS_EXPORT void WriteFld(const std::string &outname);
            
            /// Write input fields to the given filename.
            SOLVER_UTILS_EXPORT void WriteFld(
                const std::string &outname,
                MultiRegions::ExpListSharedPtr &field,
                Array<OneD, Array<OneD, NekDouble> > &fieldcoeffs,
                Array<OneD, std::string> &variables);
            
            /// Input field data from the given file.
            SOLVER_UTILS_EXPORT void ImportFld(
                const std::string &infile,
                Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);
            
            
            /// Output a field.
            /// Input field data into array from the given file.
            SOLVER_UTILS_EXPORT void ImportFld(
                const std::string &infile, 
                std::vector<std::string> &fieldStr, 
                Array<OneD, Array<OneD, NekDouble> > &coeffs);
            
            /// Output a field.
            /// Input field data into ExpList from the given file.
            SOLVER_UTILS_EXPORT void ImportFld(
                const std::string &infile, 
                MultiRegions::ExpListSharedPtr &pField, 
                std::string &pFieldName);
            
            /// Output a field.
            //        void Array_Output(
            //                const int n,
            //                std::string name,
            //                const Array<OneD, const NekDouble>&inarray,
            //                bool IsInPhysicalSpace);
            
            SOLVER_UTILS_EXPORT void WriteTecplotFile(
                const int n, 
                const std::string &name, 
                bool IsInPhysicalSpace);
            
            /// Builds map of which element holds each history point.
            SOLVER_UTILS_EXPORT void ScanForHistoryPoints();
            
            /// Probe each history point and write to file.
            SOLVER_UTILS_EXPORT void WriteHistoryData (std::ostream &out);
            
            /// Write out a full summary.
            SOLVER_UTILS_EXPORT void Summary          (std::ostream &out);
            
            /// Write out a session summary.
            SOLVER_UTILS_EXPORT void SessionSummary   (std::ostream &out);
            
            /// Write out a summary of the time parameters.
            SOLVER_UTILS_EXPORT void TimeParamSummary (std::ostream &out);
            
            SOLVER_UTILS_EXPORT inline Array<
            OneD, MultiRegions::ExpListSharedPtr> &UpdateFields();
            
            SOLVER_UTILS_EXPORT inline Array<
            OneD, MultiRegions::ExpListSharedPtr> &UpdateForces();
            

            /// Get hold of FieldInfoMap so it can be updated
            SOLVER_UTILS_EXPORT inline LibUtilities::FieldMetaDataMap 
                &UpdateFieldMetaDataMap();

            /// Return final time
            SOLVER_UTILS_EXPORT inline NekDouble GetFinalTime();
            
            SOLVER_UTILS_EXPORT inline int GetNcoeffs();
            
            SOLVER_UTILS_EXPORT inline int GetNcoeffs(const int eid);
            
            SOLVER_UTILS_EXPORT inline int GetNumExpModes();
            
            SOLVER_UTILS_EXPORT inline const Array<OneD,int> 
            GetNumExpModesPerExp();
            
            SOLVER_UTILS_EXPORT inline int GetNvariables();
            
            SOLVER_UTILS_EXPORT inline const std::string 
            GetVariable(unsigned int i);
            
            SOLVER_UTILS_EXPORT inline int GetTraceTotPoints();
            
            SOLVER_UTILS_EXPORT inline int GetTraceNpoints();
            
            SOLVER_UTILS_EXPORT inline int GetExpSize();
            
            SOLVER_UTILS_EXPORT inline int GetPhys_Offset(int n);
            
            SOLVER_UTILS_EXPORT inline int GetCoeff_Offset(int n);
            
            SOLVER_UTILS_EXPORT inline int GetTotPoints();
            
            SOLVER_UTILS_EXPORT inline int GetTotPoints(int n);
            
            SOLVER_UTILS_EXPORT inline int GetNpoints();
            
            SOLVER_UTILS_EXPORT inline int GetNumElmVelocity();
            
            SOLVER_UTILS_EXPORT inline int GetSteps();
            
            SOLVER_UTILS_EXPORT inline NekDouble GetTimeStep();
            
            SOLVER_UTILS_EXPORT inline void CopyFromPhysField(const int i,
                                                              Array<OneD, NekDouble> &output);
            
            SOLVER_UTILS_EXPORT inline void CopyToPhysField(const int i,
                                                            Array<OneD, NekDouble> &output);
            
            SOLVER_UTILS_EXPORT inline void SetStepsToOne();
            
            SOLVER_UTILS_EXPORT void ZeroPhysFields();
            
            SOLVER_UTILS_EXPORT void FwdTransFields();
            
            SOLVER_UTILS_EXPORT inline void GetFluxVector(
                const int i,
                Array<OneD, Array<OneD, NekDouble> >&physfield,
                Array<OneD, Array<OneD, NekDouble> >&flux);
            
            SOLVER_UTILS_EXPORT inline void GetFluxVector(
                const int i,
                Array<OneD, Array<OneD, NekDouble> >&physfield,
                Array<OneD, Array<OneD, NekDouble> >&fluxX,
                Array<OneD, Array<OneD, NekDouble> > &fluxY);
            
            SOLVER_UTILS_EXPORT inline void GetFluxVector(
                const int i, 
                const int j,
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &flux);
            
            SOLVER_UTILS_EXPORT inline void NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numflux);
            
            SOLVER_UTILS_EXPORT inline void NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                Array<OneD, Array<OneD, NekDouble> > &numfluxY);
            
            SOLVER_UTILS_EXPORT inline void NumFluxforScalar(
                const Array<OneD, Array<OneD, NekDouble> >         &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux);
            
            SOLVER_UTILS_EXPORT inline void NumFluxforVector(
                const Array<OneD, Array<OneD, NekDouble> >         &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                Array<OneD, Array<OneD, NekDouble> >               &qflux);
            
            SOLVER_UTILS_EXPORT inline void SetModifiedBasis(
                const bool modbasis);
            
            /// Perform a case-insensitive string comparison.
            SOLVER_UTILS_EXPORT int NoCaseStringCompare(
                const string & s1, const string& s2) ;
                
        protected:
            /// Communicator
            LibUtilities::CommSharedPtr                 m_comm;
            /// The session reader
            LibUtilities::SessionReaderSharedPtr        m_session;
            /// Array holding all dependent variables.
            Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
            /// Base fields.
            Array<OneD, MultiRegions::ExpListSharedPtr> m_base;
            /// Array holding force values.
            Array<OneD, MultiRegions::ExpListSharedPtr> m_forces;
            /// Array holding all dependent variables.
            Array<OneD, MultiRegions::ExpListSharedPtr> m_derivedfields;
            /// Pointer to boundary conditions object.
            SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
            /// Pointer to graph defining mesh.
            SpatialDomains::MeshGraphSharedPtr          m_graph;
            SpatialDomains::SpatialParametersSharedPtr  m_spatialParameters;
            /// Filename.
            std::string                                 m_filename;
            /// Name of the session.
            std::string                                 m_sessionName;
            /// Current time of simulation.
            NekDouble                                   m_time;
            /// Finish time of the simulation.
            NekDouble                                   m_fintime;
            /// Time step size
            NekDouble                                   m_timestep;
            /// Lambda constant in real system if one required.
            NekDouble                                   m_lambda;
            /// Time between checkpoints.
            NekDouble                                   m_checktime;
            /// Number of steps to take.
            int                                         m_steps;
            /// Number of steps between checkpoints.
            int                                         m_checksteps;
            /// Spatial dimension (>= expansion dim).
            int                                         m_spacedim;
            /// Expansion dimension.
            int                                         m_expdim;
            /// Flag to determine if single homogeneous mode is used.
            bool                                        m_SingleMode;
            /// Flag to determine if half homogeneous mode is used.
            bool                                        m_HalfMode;
            /// Flag to determine if use multiple homogenenous modes are used.
            bool                                        m_MultipleModes;
            /// Flag to determine if FFT is used for homogeneous transform.
            bool                                        m_useFFT;
            /// Flag to determine if dealiasing is used for homogeneous
            /// simulations.
            bool                                        m_dealiasing;
            /// Flag to determine if dealisising is usde for the
            /// Spectral/hp element discretisation.
            bool                                        m_specHP_dealiasing;
            /// Type of projection; e.g continuous or discontinuous.
            enum MultiRegions::ProjectionType           m_projectionType;
            /// Array holding trace normals for DG simulations in the forwards
            /// direction.
            Array<OneD, Array<OneD, NekDouble> >        m_traceNormals;
            /// 1 x nvariable x nq
            Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_gradtan;
            /// 2 x m_spacedim x nq
            Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_tanbasis;
            /// Flag to indicate if the fields should be checked for
            /// singularity.
            Array<OneD, bool>                           m_checkIfSystemSingular;
            
            /// Map to identify relevant solver info to dump in output fields
            LibUtilities::FieldMetaDataMap            m_fieldMetaDataMap;

            /// Number of Quadrature points used to work out the error
            int  m_NumQuadPointsError;
            
            /// Parameter for homogeneous expansions
            enum HomogeneousType
            {
                eHomogeneous1D,
                eHomogeneous2D,
                eHomogeneous3D,
                eNotHomogeneous
            };
            
            
            
            enum HomogeneousType m_HomogeneousType;
            
            NekDouble m_LhomX;  ///< physical length in X direction (if homogeneous)
            NekDouble m_LhomY;  ///< physical length in Y direction (if homogeneous)
            NekDouble m_LhomZ;  ///< physical length in Z direction (if homogeneous)
            
            int m_npointsX;     ///< number of points in X direction (if homogeneous)
            int m_npointsY;     ///< number of points in Y direction (if homogeneous)
            int m_npointsZ;     ///< number of points in Z direction (if homogeneous)
            
            int m_HomoDirec;    ///< number of homogenous directions
            
            int m_NumMode;      ///< Mode to use in case of single mode analysis
            
            
            /// Initialises EquationSystem class members.
            SOLVER_UTILS_EXPORT EquationSystem( const LibUtilities::SessionReaderSharedPtr& pSession);
            
            // Here for consistency purposes with old version
            int nocase_cmp(const string & s1, const string& s2)
            {
                return NoCaseStringCompare(s1,s2);
            }
            
            SOLVER_UTILS_EXPORT virtual void v_InitObject();
            
            /// Evaluates the boundary conditions at the given time.
            SOLVER_UTILS_EXPORT void SetBoundaryConditions(NekDouble time);
            
            /// Virtual function for initialisation implementation.
            SOLVER_UTILS_EXPORT virtual void v_DoInitialise();
            
            /// Virtual function for solve implementation.
            SOLVER_UTILS_EXPORT virtual void v_DoSolve();
            
            /// Virtual function for the L_inf error computation between fields and a given exact solution.
            SOLVER_UTILS_EXPORT virtual NekDouble v_LinfError(
                unsigned int field,
                const Array<OneD, NekDouble> &exactsoln = NullNekDouble1DArray);
            
            /// Virtual function for the L_2 error computation between fields and a given exact solution.
            SOLVER_UTILS_EXPORT virtual NekDouble v_L2Error(
                unsigned int field, 
                const Array<OneD, NekDouble> &exactsoln = NullNekDouble1DArray, 
                bool Normalised = false);
            
            /// Virtual function for transformation to physical space.
            SOLVER_UTILS_EXPORT virtual void v_TransCoeffToPhys();
            
            /// Virtual function for transformation to coefficient space.
            SOLVER_UTILS_EXPORT virtual void v_TransPhysToCoeff();
            
            /// Virtual function for printing summary information.
            SOLVER_UTILS_EXPORT virtual void v_PrintSummary(std::ostream &out);
            
            SOLVER_UTILS_EXPORT virtual void v_SetInitialConditions(
                NekDouble initialtime = 0.0,
                bool dumpInitialConditions = true);
            
            SOLVER_UTILS_EXPORT virtual void v_EvaluateExactSolution(
                unsigned int field,
                Array<OneD, NekDouble> &outfield,
                const NekDouble time);
            
            //Initialise m_base in order to store the base flow from a file 
            SOLVER_UTILS_EXPORT void SetUpBaseFields(SpatialDomains::MeshGraphSharedPtr &mesh);
            
            // Fill m_base with the values stored in a fld file
            SOLVER_UTILS_EXPORT void ImportFldBase(
                std::string pInfile, 
                SpatialDomains::MeshGraphSharedPtr pGraph);
            
            // Ouptut field information
            SOLVER_UTILS_EXPORT virtual void v_Output(void);
            
            // Get pressure field if available
            SOLVER_UTILS_EXPORT virtual MultiRegions::ExpListSharedPtr v_GetPressure(void); 
            
        private:
            
            SOLVER_UTILS_EXPORT virtual Array<OneD, bool> v_GetSystemSingularChecks();
            SOLVER_UTILS_EXPORT virtual int v_GetForceDimension();
            SOLVER_UTILS_EXPORT virtual void v_GetFluxVector(
                const int i, Array<OneD,
                Array<OneD, NekDouble> >&physfield,
                Array<OneD, Array<OneD, NekDouble> >&flux);
            
            SOLVER_UTILS_EXPORT virtual void v_GetFluxVector(
                const int i, const int j,
                Array<OneD, Array<OneD, NekDouble> >&physfield,
                Array<OneD, Array<OneD, NekDouble> >&flux);
            
            SOLVER_UTILS_EXPORT virtual void v_GetFluxVector(
                const int i, Array<OneD,
                Array<OneD, NekDouble> >&physfield,
                Array<OneD, Array<OneD, NekDouble> >&fluxX,
                Array<OneD, Array<OneD, NekDouble> > &fluxY);
            
            SOLVER_UTILS_EXPORT virtual void v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numflux);
            
            SOLVER_UTILS_EXPORT virtual void v_NumericalFlux(
                Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                Array<OneD, Array<OneD, NekDouble> > &numfluxY);
            
            SOLVER_UTILS_EXPORT virtual void v_NumFluxforScalar(
                const Array<OneD, Array<OneD, NekDouble> >         &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux);
            
            SOLVER_UTILS_EXPORT virtual void v_NumFluxforVector(
                const Array<OneD, Array<OneD, NekDouble> >         &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                Array<OneD, Array<OneD, NekDouble > >              &qflux);
        };
        
        
        /**
         * This is the second part of the two-phase initialisation process.
         * Calls to virtual functions will correctly resolve to the derived class
         * during this phase of the construction.
         */
        inline void EquationSystem::InitObject()
        {
            v_InitObject();
        }
        
        
        /**
         * This allows initialisation of the solver which cannot be completed
         * during object construction (such as setting of initial conditions).
         *
         * Public interface routine to virtual function implementation.
         */
        inline void EquationSystem::DoInitialise()
        {
            v_DoInitialise();
        }
        
        
        /**
         * Performs the transformation from coefficient to physical space.
         *
         * Public interface routine to virtual function implementation.
         */
        inline void EquationSystem::TransCoeffToPhys(void)
        {
            v_TransCoeffToPhys();
        }
        
        /**
         * Performs the transformation from physical to coefficient space.
         *
         * Public interface routine to virtual function implementation.
         */
        inline void EquationSystem::TransPhysToCoeff(void)
        {
            v_TransPhysToCoeff();
        }
        
        
        /**
         * Performs the actual solve.
         *
         * Public interface routine to virtual function implementation.
         */
        inline void EquationSystem::DoSolve(void)
        {
            v_DoSolve();
        }
        
        
        /**
         * Perform output operations after solve.
         */
        inline void EquationSystem::Output(void)
        {
            v_Output();
        }
        
        /**
         * L_inf Error computation
         * Public interface routine to virtual function implementation.
         */
        inline NekDouble EquationSystem::LinfError(unsigned int field, const Array<OneD,NekDouble> &exactsoln)
        {
            return v_LinfError(field, exactsoln);
        }
        
        /**
         * L_2 Error computation
         * Public interface routine to virtual function implementation.
         */
        inline NekDouble EquationSystem::L2Error(unsigned int field, const Array<OneD,NekDouble> &exactsoln, bool Normalised)
        {
            return v_L2Error(field, exactsoln, Normalised);
        }
        
        /**
         * Get Pressure field if available
         */
        inline  MultiRegions::ExpListSharedPtr EquationSystem::GetPressure(void)
        {
            return v_GetPressure();
        }
        
        /**
         * Prints a summary of variables and problem parameters.
         *
         * Public interface routine to virtual function implementation.
         *
         * @param   out             The ostream object to write to.
         */
        inline void EquationSystem::PrintSummary(std::ostream &out)
        {
            if (m_session->GetComm()->GetRank() == 0)
            {
                out << "=======================================================================" << endl;
                out << "\tEquation Type   : " << m_session->GetSolverInfo("EQTYPE") << endl;
                SessionSummary(out);
                
                v_PrintSummary(out);
                
                out << "=======================================================================" << endl;
            }
        }
        
        
        inline void EquationSystem::SetLambda(NekDouble lambda)
        {
            m_lambda = lambda;
        }
        
        inline void EquationSystem::SetInitialConditions(NekDouble initialtime,
                                                         bool dumpInitialConditions)
        {
            v_SetInitialConditions(initialtime,dumpInitialConditions);
        }
        
        /// Evaluates an exact solution
        inline void EquationSystem::EvaluateExactSolution(int field,
                                                          Array<OneD, NekDouble> &outfield,
                                                          const NekDouble time)
        {
            v_EvaluateExactSolution(field, outfield, time);
        }
        
        inline Array<OneD, MultiRegions::ExpListSharedPtr> &EquationSystem::UpdateFields(void)
        {
            return m_fields;
        }
        
        inline Array<OneD, MultiRegions::ExpListSharedPtr> &EquationSystem::UpdateForces(void)
        {
            return m_forces;
        }
        
        /// Return final time
        inline NekDouble EquationSystem::GetFinalTime()
        {
            return m_time;
        }
        
        inline int EquationSystem::GetNcoeffs(void)
        {
            return m_fields[0]->GetNcoeffs();
        }
        
        inline int EquationSystem::GetNcoeffs(const int eid)
        {
            return m_fields[0]->GetNcoeffs(eid);
        }
        
        inline int EquationSystem::GetNumExpModes(void)
        {
            return m_graph->GetExpansions().begin()->second->m_basisKeyVector[0]
            .GetNumModes();
        }
        
        inline const Array<OneD,int> EquationSystem::GetNumExpModesPerExp(void)
        {
            return m_fields[0]->EvalBasisNumModesMaxPerExp();
        }
        
        inline int EquationSystem::GetNvariables(void)
        {
            return m_fields.num_elements();
        }
        
        inline const std::string EquationSystem::GetVariable(unsigned int i)
        {
            return m_boundaryConditions->GetVariable(i);
        }
        
        inline int EquationSystem::GetTraceTotPoints(void)
        {
            return GetTraceNpoints();
        }
        
        
        inline int EquationSystem::GetTraceNpoints(void)
        {
            return m_fields[0]->GetTrace()->GetNpoints();
        }
        
        inline int EquationSystem::GetExpSize(void)
        {
            return m_fields[0]->GetExpSize();
        }
        
        inline int EquationSystem::GetPhys_Offset(int n)
        {
            return m_fields[0]->GetPhys_Offset(n);
        }
        
        inline int EquationSystem::GetCoeff_Offset(int n)
        {
            return m_fields[0]->GetCoeff_Offset(n);
        }
        
        inline int EquationSystem::GetTotPoints(void)
        {
            return m_fields[0]->GetNpoints();
        }
        
        inline int EquationSystem::GetTotPoints(int n)
        {
            return m_fields[0]->GetTotPoints(n);
        }
        
        inline int EquationSystem::GetNpoints(void)
        {
            return m_fields[0]->GetNpoints();
        }
        
        inline int EquationSystem::GetNumElmVelocity(void)
        {
            return (m_fields.num_elements() - 1);
        }
        
        inline int EquationSystem::GetSteps(void)
        {
            return m_steps;
        }
        
        inline NekDouble EquationSystem::GetTimeStep(void)
        {
            return m_timestep;
        }
        
        inline void EquationSystem::SetStepsToOne(void)
        {
            m_steps=1;
        }
        
        inline void EquationSystem::CopyFromPhysField(const int i,
                                                      Array<OneD, NekDouble> &output)
        {
            Vmath::Vcopy(output.num_elements(), m_fields[i]->GetPhys(), 1, output, 1 );
        }
        
        inline void EquationSystem::CopyToPhysField(const int i,
                                                    Array<OneD, NekDouble> &output)
        {
            Vmath::Vcopy(output.num_elements(), output, 1, m_fields[i]->UpdatePhys(), 1 );
        }
        
        inline void EquationSystem::GetFluxVector(const int i,
                                                  Array<OneD, Array<OneD, NekDouble> >&physfield,
                                                  Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            v_GetFluxVector(i,physfield, flux);
        }
        
        inline void EquationSystem::GetFluxVector(const int i,
                                                  Array<OneD, Array<OneD, NekDouble> >&physfield,
                                                  Array<OneD, Array<OneD, NekDouble> >&fluxX,
                                                  Array<OneD, Array<OneD, NekDouble> > &fluxY)
        {
            v_GetFluxVector(i,physfield, fluxX, fluxY);
        }
        
        inline void EquationSystem::GetFluxVector(const int i, const int j,
                                                  Array<OneD, Array<OneD, NekDouble> > &physfield,
                                                  Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            v_GetFluxVector(i,j,physfield,flux);
        }
        
        inline void EquationSystem::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                                                  Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            v_NumericalFlux(physfield, numflux);
        }
        
        inline void EquationSystem::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                                                  Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                                                  Array<OneD, Array<OneD, NekDouble> > &numfluxY)
        {
            v_NumericalFlux(physfield, numfluxX, numfluxY);
        }
        
        inline void EquationSystem::NumFluxforScalar(
            const Array<OneD, Array<OneD, NekDouble> >   &ufield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            v_NumFluxforScalar(ufield, uflux);
        }
        
        inline void EquationSystem::NumFluxforVector(            
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            v_NumFluxforVector(ufield, qfield, qflux);
        }
    }
}

#endif
