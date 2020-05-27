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

#include <SolverUtils/Core/SessionFunction.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Core/Misc.h>

namespace Nektar
{
namespace FieldUtils {
class Interpolator;
}
    namespace SolverUtils
    {
        class EquationSystem;

        /// A shared pointer to an EquationSystem object
        typedef std::shared_ptr<EquationSystem> EquationSystemSharedPtr;
        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef LibUtilities::NekFactory<
            std::string,
            EquationSystem,
            const LibUtilities::SessionReaderSharedPtr&,
            const SpatialDomains::MeshGraphSharedPtr&
        > EquationSystemFactory;
        SOLVER_UTILS_EXPORT EquationSystemFactory& GetEquationSystemFactory();

        /// A base class for describing how to solve specific equations.
        class EquationSystem : public std::enable_shared_from_this<EquationSystem>
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
            SOLVER_UTILS_EXPORT inline NekDouble LinfError(unsigned int field,
                const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);

            /// Get Session name
            SOLVER_UTILS_EXPORT std::string GetSessionName()
            {
                return m_sessionName;
            }

            template<class T>
            std::shared_ptr<T> as()
            {
                return std::dynamic_pointer_cast<T>( shared_from_this() );
            }

            /// Reset Session name
            SOLVER_UTILS_EXPORT void ResetSessionName(std::string newname)
            {
                m_sessionName = newname;
            }

            /// Get Session name
            SOLVER_UTILS_EXPORT LibUtilities::SessionReaderSharedPtr GetSession()
            {
                return m_session;
            }

            /// Get pressure field if available
            SOLVER_UTILS_EXPORT MultiRegions::ExpListSharedPtr GetPressure();

            SOLVER_UTILS_EXPORT inline void ExtraFldOutput(
                std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                std::vector<std::string>             &variables);

            /// Print a summary of parameters and solver characteristics.
            SOLVER_UTILS_EXPORT inline void PrintSummary(std::ostream &out);

            /// Set parameter m_lambda
            SOLVER_UTILS_EXPORT inline void SetLambda(NekDouble lambda);

            /// Get a SessionFunction by name
            SOLVER_UTILS_EXPORT SessionFunctionSharedPtr GetFunction(
                std::string name,
                const MultiRegions::ExpListSharedPtr &field = MultiRegions::NullExpListSharedPtr,
                bool cache = false);

            /// Initialise the data in the dependent fields.
            SOLVER_UTILS_EXPORT inline void SetInitialConditions(
                NekDouble initialtime = 0.0,
                bool dumpInitialConditions = true,
                const int domain = 0);

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

            /// Write checkpoint file of #m_fields.
            SOLVER_UTILS_EXPORT void Checkpoint_Output(const int n);

            /// Write checkpoint file of custom data fields.
            SOLVER_UTILS_EXPORT void Checkpoint_Output(
                const int n,
                MultiRegions::ExpListSharedPtr &field,
                std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                std::vector<std::string> &variables);

            /// Write base flow file of #m_fields.
            SOLVER_UTILS_EXPORT void Checkpoint_BaseFlow(const int n);

            /// Write field data to the given filename.
            SOLVER_UTILS_EXPORT void WriteFld(const std::string &outname);

            /// Write input fields to the given filename.
            SOLVER_UTILS_EXPORT void WriteFld(
                const std::string &outname,
                MultiRegions::ExpListSharedPtr &field,
                std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                std::vector<std::string> &variables);

            /// Input field data from the given file.
            SOLVER_UTILS_EXPORT void ImportFld(
                const std::string &infile,
                Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

            /// Input field data from the given file to multiple domains
            SOLVER_UTILS_EXPORT void ImportFldToMultiDomains(
                const std::string &infile,
                Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const int ndomains);

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

            /// Write out a session summary.
            SOLVER_UTILS_EXPORT void SessionSummary   (SummaryList& vSummary);

            SOLVER_UTILS_EXPORT inline Array<
            OneD, MultiRegions::ExpListSharedPtr> &UpdateFields();


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

            SOLVER_UTILS_EXPORT inline int GetSteps();

            SOLVER_UTILS_EXPORT inline NekDouble GetTimeStep();

            SOLVER_UTILS_EXPORT inline void CopyFromPhysField(const int i,
                    Array<OneD, NekDouble> &output);

            SOLVER_UTILS_EXPORT inline void CopyToPhysField(const int i,
                    Array<OneD, NekDouble> &output);

            SOLVER_UTILS_EXPORT inline void SetSteps(const int steps);

            SOLVER_UTILS_EXPORT void ZeroPhysFields();

            SOLVER_UTILS_EXPORT void FwdTransFields();

            SOLVER_UTILS_EXPORT inline void SetModifiedBasis(
                const bool modbasis);

            SOLVER_UTILS_EXPORT int GetCheckpointNumber()
            {
                return m_nchk;
            }

            SOLVER_UTILS_EXPORT void SetCheckpointNumber(int num)
            {
                m_nchk = num;
            }

            SOLVER_UTILS_EXPORT int GetCheckpointSteps()
            {
                return m_checksteps;
            }

            SOLVER_UTILS_EXPORT void SetCheckpointSteps(int num)
            {
                m_checksteps = num;
            }

            SOLVER_UTILS_EXPORT void SetTime(
                const NekDouble time)
            {
                m_time = time;
            }

            SOLVER_UTILS_EXPORT void SetInitialStep(
                const int step)
            {
                m_initialStep = step;
            }

            /// Evaluates the boundary conditions at the given time.
            SOLVER_UTILS_EXPORT void SetBoundaryConditions(NekDouble time);

            /// Virtual function to identify if operator is negated in DoSolve
            SOLVER_UTILS_EXPORT virtual bool v_NegatedOp();

        protected:
            /// Communicator
            LibUtilities::CommSharedPtr                 m_comm;
            /// The session reader
            LibUtilities::SessionReaderSharedPtr        m_session;
            /// Map of known SessionFunctions
            std::map<std::string, SolverUtils::SessionFunctionSharedPtr> m_sessionFunctions;
            /// Field input/output
            LibUtilities::FieldIOSharedPtr              m_fld;
            /// Array holding all dependent variables.
            Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
            /// Pointer to boundary conditions object.
            SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
            /// Pointer to graph defining mesh.
            SpatialDomains::MeshGraphSharedPtr          m_graph;
            /// Name of the session.
            std::string                                 m_sessionName;
            /// Current time of simulation.
            NekDouble                                   m_time;
            /// Number of the step where the simulation should begin
            int                                         m_initialStep;
            /// Finish time of the simulation.
            NekDouble                                   m_fintime;
            /// Time step size
            NekDouble                                   m_timestep;
            /// Time step size
            NekDouble                                   m_timestepMax = -1.0;
            
            /// Lambda constant in real system if one required.
            NekDouble                                   m_lambda;
            /// Time between checkpoints.
            NekDouble                                   m_checktime;
            /// Number of checkpoints written so far
            int                                         m_nchk;
            /// Number of steps to take.
            int                                         m_steps;
            /// Number of steps between checkpoints.
            int                                         m_checksteps;
            /// Spatial dimension (>= expansion dim).
            int                                         m_spacedim;
            /// Expansion dimension.
            int                                         m_expdim;
            /// Flag to determine if single homogeneous mode is used.
            bool                                        m_singleMode;
            /// Flag to determine if half homogeneous mode is used.
            bool                                        m_halfMode;
            /// Flag to determine if use multiple homogenenous modes are used.
            bool                                        m_multipleModes;
            /// Flag to determine if FFT is used for homogeneous transform.
            bool                                        m_useFFT;
            /**
             * \brief Flag to determine if dealiasing is used for
             * homogeneous simulations.
             */
            bool m_homogen_dealiasing;
            /**
             * \brief Flag to determine if dealisising is usde for the
             * Spectral/hp element discretisation.
             */
            bool                                        m_specHP_dealiasing;
            /// Type of projection; e.g continuous or discontinuous.
            enum MultiRegions::ProjectionType           m_projectionType;
            /// Array holding trace normals for DG simulations in the forwards direction.
            Array<OneD, Array<OneD, NekDouble> >        m_traceNormals;
            /// Flag to indicate if the fields should be checked for singularity.
            Array<OneD, bool>                           m_checkIfSystemSingular;
            /// Map to identify relevant solver info to dump in output fields
            LibUtilities::FieldMetaDataMap              m_fieldMetaDataMap;

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

            /// Initialises EquationSystem class members.
            SOLVER_UTILS_EXPORT EquationSystem(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const SpatialDomains::MeshGraphSharedPtr& pGraph);

            SOLVER_UTILS_EXPORT virtual void v_InitObject();

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

            /// Virtual function for generating summary information.
            SOLVER_UTILS_EXPORT virtual void v_GenerateSummary(SummaryList& l);

            SOLVER_UTILS_EXPORT virtual void v_SetInitialConditions(
                NekDouble initialtime = 0.0,
                bool dumpInitialConditions = true,
                const int domain = 0);

            SOLVER_UTILS_EXPORT virtual void v_EvaluateExactSolution(
                unsigned int field,
                Array<OneD, NekDouble> &outfield,
                const NekDouble time);

            // Ouptut field information
            SOLVER_UTILS_EXPORT virtual void v_Output(void);

            // Get pressure field if available
            SOLVER_UTILS_EXPORT virtual MultiRegions::ExpListSharedPtr v_GetPressure(void);

            SOLVER_UTILS_EXPORT virtual void v_ExtraFldOutput(
                std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                std::vector<std::string>             &variables);

            static std::string equationSystemTypeLookupIds[];

        private:

            SOLVER_UTILS_EXPORT virtual Array<OneD, bool> v_GetSystemSingularChecks();

            SOLVER_UTILS_EXPORT void PrintProgressbar(const int position,
                                                      const int goal) const
            {
                LibUtilities::PrintProgressbar(position, goal, "Interpolating");
            }
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
         * Append the coefficients and name of variables with solver specific
         * extra variables
         *
         * @param fieldcoeffs     Vector with coefficients
         * @param variables       Vector with name of variables
         */
        inline void EquationSystem::ExtraFldOutput(
                std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                std::vector<std::string>             &variables)
        {
            v_ExtraFldOutput(fieldcoeffs, variables);
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
                std::vector<std::pair<std::string, std::string> > vSummary;
                v_GenerateSummary(vSummary);

                out << "=======================================================================" << std::endl;
                for (auto &x : vSummary)
                {
                    out << "\t";
                    out.width(20);
                    out << x.first << ": " << x.second << std::endl;
                }
                out << "=======================================================================" << std::endl;
            }
        }

        inline void EquationSystem::SetLambda(NekDouble lambda)
        {
            m_lambda = lambda;
        }

        inline void EquationSystem::SetInitialConditions(NekDouble initialtime,
                                                         bool dumpInitialConditions,
                                                         const int domain)
        {
            v_SetInitialConditions(initialtime,dumpInitialConditions,domain);
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
            return m_session->GetVariables().size();
        }

        inline const std::string EquationSystem::GetVariable(unsigned int i)
        {
            return m_session->GetVariable(i);
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

        inline int EquationSystem::GetSteps(void)
        {
            return m_steps;
        }

        inline NekDouble EquationSystem::GetTimeStep(void)
        {
            return m_timestep;
        }

        inline void EquationSystem::SetSteps(const int steps)
        {
            m_steps = steps;
        }

        inline void EquationSystem::CopyFromPhysField(const int i,
                                                      Array<OneD, NekDouble> &output)
        {
            Vmath::Vcopy(output.size(), m_fields[i]->GetPhys(), 1, output, 1 );
        }

        inline void EquationSystem::CopyToPhysField(const int i,
                                                    Array<OneD, NekDouble> &output)
        {
            Vmath::Vcopy(output.size(), output, 1, m_fields[i]->UpdatePhys(), 1 );
        }
    }
}

#endif
