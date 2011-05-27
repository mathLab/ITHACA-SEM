///////////////////////////////////////////////////////////////////////////////
//
// File ADRBase.h
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
// Description: Basic class for AdvectionDiffusionReaction class,
// Euler Class, ShallowWaterEquations and BoussinesqEquations
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H
#define NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/HistoryPoints.h>
#include <SpatialDomains/SpatialData.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>

#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>

namespace Nektar
{
    static std::string NekNullString;

    /// Base class for the development of solvers.
    class ADRBase
    {
    public:

        /// Default constructor.
        ADRBase();

        /// Constructor.
        ADRBase(const string &fileStringName,
                bool UseInputFileForProjectionType = false,
                bool UseContinuousField = false);

        /// Initialise dependent variable fields.
        void SetADRBase(SpatialDomains::MeshGraphSharedPtr &graph,
                        int nvariables);

        /// Initialise the data in the dependent fields.
        void SetInitialConditions(NekDouble initialtime = 0.0,
                                  bool dumpInitialConditions = true)
        {
            v_SetInitialConditions(initialtime,dumpInitialConditions);
        }
       /// Decide the kind of forcing functions
        void SetInitialForce(NekDouble initialtime=0.0) ;
	/// populate m_forces
        void CalcForce(Array<OneD, MultiRegions::ExpListSharedPtr> &force);

       ///Initialise the dimendion of forcing functions
        void InitialiseForcingFunctions(
                        Array<OneD, MultiRegions::ExpListSharedPtr> &force);

        /// Populate given fields with the forcing function from session.
        void SetPhysForcingFunctions(
                        Array<OneD, MultiRegions::ExpListSharedPtr> &force);

        /// Populate given field with the exact solution from session.
        void EvaluateExactSolution(int field,
                                   Array<OneD, NekDouble > &exactsoln,
                                   const NekDouble time);

        /// Populate given fields with a user expression from session.
        void EvaluateUserDefinedEqn(
                        Array<OneD, Array<OneD, NekDouble> > &outfield);

        /// Compute the L2 error between fields and a given exact solution.
        NekDouble L2Error(int field,
                          const Array<OneD,NekDouble> &exactsoln, bool Normalised = false);

        /// Compute the L2 error of the fields
        NekDouble L2Error(int field, bool Normalised = false)
        {
            return L2Error(field,NullNekDouble1DArray,Normalised);
        }

        /// Compute the L_inf error between fields and a given exact solution.
        NekDouble LinfError(int field,
                const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);
		
		///Compute error (L2 and L_inf) over an larger set of quadrature points return [L2 Linf]
		Array<OneD,NekDouble> ErrorExtraPoints(int field);

        /// Compute the inner product \f$ (\nabla \phi \cdot F) \f$.
        void WeakAdvectionGreensDivergenceForm(
                const Array<OneD, Array<OneD, NekDouble> > &F,
                Array<OneD, NekDouble> &outarray);

        /// Compute the inner product \f$ (\phi, \nabla \cdot F) \f$.
        void WeakAdvectionDivergenceForm(
                const Array<OneD, Array<OneD, NekDouble> > &F,
                Array<OneD, NekDouble> &outarray);

        /// Compute the inner product \f$ (\phi, V\cdot \nabla u) \f$.
        void WeakAdvectionNonConservativeForm(
                const Array<OneD, Array<OneD, NekDouble> > &V,
                const Array<OneD, const NekDouble> &u,
                Array<OneD, NekDouble> &outarray);

        /// Compute the non-conservative advection \f$ (V \cdot \nabla u) \f$.
        void AdvectionNonConservativeForm(
                const Array<OneD, Array<OneD, NekDouble> > &V,
                const Array<OneD, const NekDouble> &u,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wk = NullNekDouble1DArray);

        /// Calculate the weak discontinuous Galerkin advection.
        void WeakDGAdvection(
                const Array<OneD, Array<OneD, NekDouble> >& InField,
                Array<OneD, Array<OneD, NekDouble> >& OutField,
                bool NumericalFluxIncludesNormal = true,
                bool InFieldIsInPhysSpace = false,
                int nvariables = 0);

        /// Calculate weak DG Diffusion in the LDG form.
        void WeakDGDiffusion(
                const Array<OneD, Array<OneD, NekDouble> >& InField,
                Array<OneD, Array<OneD, NekDouble> >& OutField,
                bool NumericalFluxIncludesNormal = true,
                bool InFieldIsInPhysSpace = false);

        /// Write field to file.
        void Output();
        /// Write field to file.
        void Output(MultiRegions::ExpListSharedPtr &field, Array< OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables);


        /// Write checkpoint file.
        void Checkpoint_Output(const int n);
        /// Write checkpoint file.
        void Checkpoint_Output(const int n, MultiRegions::ExpListSharedPtr &field, Array< OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables);
        
        /// Write field data to the given filename.
        void WriteFld(std::string &outname);

        /// Write input fields to the given filename.
        void WriteFld(std::string &outname, MultiRegions::ExpListSharedPtr &field, Array<OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables);

        /// Input field data from the given file.
        void ImportFld(std::string &infile);
        
        ///Input force data from the given file.
        void ImportFldForce(std::string infile);

        

        /// Output a field.
        void Array_Output(const int n, std::string name,
                          const Array<OneD, const NekDouble>&inarray,
                          bool IsInPhysicalSpace);

	void WriteTecplotFile(const int n, std::string name, bool IsInPhysicalSpace);

        /// Builds map of which element holds each history point.
        void ScanForHistoryPoints();

        /// Probe each history point and write to file.
        void WriteHistoryData (std::ostream &out);

        /// Write out a full summary.
        void Summary          (std::ostream &out);

        /// Write out a session summary.
        void SessionSummary   (std::ostream &out);

        /// Write out a summary of the time parameters.
        void TimeParamSummary (std::ostream &out);

        inline Array<OneD, MultiRegions::ExpListSharedPtr> &UpdateFields(void)
        {
            return m_fields;
        }

	/// Return final time
	inline NekDouble GetFinalTime()
	{
	  return m_time;
	}

        inline int GetNcoeffs(void)
        {
            return m_fields[0]->GetNcoeffs();
        }

        inline int GetNcoeffs(const int eid)
        {
            return m_fields[0]->GetNcoeffs(eid);
        }

        inline int GetNumExpModes(void)
        {
            return m_graph->GetExpansions()[0]->m_basisKeyVector[0]
                                                        .GetNumModes();
        }

	inline const Array<OneD,int> GetNumExpModesPerExp(void)
        {
	  return m_fields[0]->EvalBasisNumModesMaxPerExp();
        }

        inline int GetNvariables(void)
        {
            return m_fields.num_elements();
        }

        inline const std::string &GetVariable(unsigned int i)
        {
            return m_boundaryConditions->GetVariable(i);
        }

        inline int GetTraceTotPoints(void)
        {
            return GetTraceNpoints();
        }

        inline int GetTraceNpoints(void)
        {
            switch(m_expdim)
            {
                case 1:
                    // can't have two &GetTrace in ExpList.h hmm...
                    //return m_fields[0]->GetTrace().num_elements();
                    break;
                case 2:
                case 3:
                    return m_fields[0]->GetTrace()->GetNpoints();
                    break;
                default:
                    ASSERTL0(false,"illegal expansion dimension");
            }
        }

	inline int GetExpSize(void)
	{
	  return m_fields[0]->GetExpSize();
	}

	inline int GetPhys_Offset(int n)
	{
	  return m_fields[0]->GetPhys_Offset(n);
	}

        inline int GetCoeff_Offset(int n)
        {
            return m_fields[0]->GetCoeff_Offset(n);
        }

        inline int GetTotPoints(void)
        {
            return m_fields[0]->GetNpoints();
        }

        inline int GetTotPoints(int n)
        {
            return m_fields[0]->GetTotPoints(n);
        }

        inline int GetNpoints(void)
        {
            return m_fields[0]->GetNpoints();
        }

        inline int GetSteps(void)
        {
            return m_steps;
        }

        inline NekDouble GetParameter(const std::string &parmName)
        {
            return m_boundaryConditions->GetParameter(parmName);
        }

        void ZeroPhysFields(void);

        void FwdTransFields(void);

        /// Type of Galerkin projection.
        enum ProjectionType
        {
            eGalerkin,
            eDiscontinuousGalerkin
        };

        //-----------------------------------------------------------
        // virtual functions wrappers

        void GetFluxVector(const int i,
                            Array<OneD, Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            v_GetFluxVector(i,physfield, flux);
        }

        void GetFluxVector(const int i,
                            Array<OneD, Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&fluxX,
                            Array<OneD, Array<OneD, NekDouble> > &fluxY)
        {
            v_GetFluxVector(i,physfield, fluxX, fluxY);
        }

        virtual void GetFluxVector(const int i, const int j,
                            Array<OneD, Array<OneD, NekDouble> > &physfield,
                            Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            v_GetFluxVector(i,j,physfield,flux);
        }

        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            v_NumericalFlux(physfield, numflux);
        }

        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                           Array<OneD, Array<OneD, NekDouble> > &numfluxY)
        {
            v_NumericalFlux(physfield, numfluxX, numfluxY);
        }

        void NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            v_NumFluxforScalar(ufield, uflux);
        }

        void NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                  Array<OneD, Array<OneD, NekDouble> >  &qflux)
        {
            v_NumFluxforVector(ufield,qfield, qflux);
        }

        NekDouble AdvectionSphere(const NekDouble x0j, const NekDouble x1j,
                            const NekDouble x2j, const NekDouble time)
        {
            return v_AdvectionSphere(x0j, x1j, x2j, time);
        }

        NekDouble Morphogenesis(const int field, const NekDouble x0j,
                            const NekDouble x1j, const NekDouble x2j,
                            const NekDouble time)
        {
            return v_Morphogenesis(field, x0j, x1j, x2j, time);
        }
		
		/// Number of Quadrature points used to work out the error
		int  m_NumQuadPointsError;
		
		bool m_UseContCoeff;
		
		///Parameter for homogeneous expansions
		
		enum HomogeneousType
		{
			eHomogeneous1D,
			eHomogeneous2D,
			eHomogeneous3D,
			eNotHomogeneous
		};
		
		bool m_useFFT;               ///< flag to determine if use or not the FFT for transformations
		
		enum HomogeneousType m_HomogeneousType;
		
		NekDouble m_LhomX;           ///< physical length in X direction (if homogeneous) 
		NekDouble m_LhomY;           ///< physical length in Y direction (if homogeneous)
		NekDouble m_LhomZ;           ///< physical length in Z direction (if homogeneous)
		
		int m_npointsX;              ///< number of points in X direction (if homogeneous)
		int m_npointsY;              ///< number of points in Y direction (if homogeneous)
		int m_npointsZ;              ///< number of points in Z direction (if homogeneous)
		
		int m_HomoDirec;             ///< number of homogenous directions
		
    protected:
        /// Array holding all dependent variables.
        Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
        /// Array holding force values.
        Array<OneD, MultiRegions::ExpListSharedPtr> m_forces;
        /// variable that determine if the force is necessary or not.
        bool m_bforce;
        /// Array holding all dependent variables.
        Array<OneD, MultiRegions::ExpListSharedPtr> m_derivedfields;
        /// dimension force array
        int m_FDim;
        /// Pointer to boundary conditions object.
        SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
        /// Pointer to history data object.
        SpatialDomains::HistorySharedPtr            m_historyPoints;
        /// Pointer to graph defining mesh.
        SpatialDomains::MeshGraphSharedPtr          m_graph;
		
        std::list<std::pair<SpatialDomains::VertexComponentSharedPtr, int> >
                                                    m_historyList;

        SpatialDomains::SpatialParametersSharedPtr  m_spatialParameters;

        std::string m_filename;      ///< Filename
        std::string m_sessionName;   ///< Name of the sessions
        NekDouble m_time;            ///< Continous time
	NekDouble m_fintime;         ///< time to be taken during the simulation
        NekDouble m_timestep;        ///< Time step size
        int m_steps;                 ///< Number of steps to take
        int m_checksteps;            ///< Number of steps between checkpoints
        int m_spacedim;              ///< Spatial dimension (> expansion dim)
        int m_expdim;                ///< Dimension of the expansion

        /// Type of projection, i.e. Galerkin or DG.
        enum ProjectionType m_projectionType;

        /// Array holding the forward normals.
        Array<OneD, Array<OneD, NekDouble> > m_traceNormals;
        /// 1 x nvariable x nq
        Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_gradtan;
        /// 2 x m_spacedim x nq
        Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_tanbasis;

        /// Flag to indicate if the fields should be checked for singularity.
        Array<OneD, bool> checkIfSystemSingular;

        /// Perform a case-insensitive string comparison.
        int NoCaseStringCompare(const string & s1, const string& s2) ;

        // Here for consistency purposes with old version
        int nocase_cmp(const string & s1, const string& s2)
        {
            return NoCaseStringCompare(s1,s2);
        }

        /// Check for and load an integer parameter
        void LoadParameter(std::string name, int &var, int def = 0);

        /// Check for and load a double precision parameter
        void LoadParameter(std::string name, NekDouble &var, NekDouble def= 0.0);

	virtual void v_SetInitialConditions(NekDouble initialtime = 0.0,
					    bool dumpInitialConditions = true);

        virtual void v_EvaluateExactSolution(int field,
					     Array<OneD, NekDouble > &outfield,
					     const NekDouble time = 0.0);

    private:

        virtual void v_GetFluxVector(const int i, Array<OneD,
                            Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            ASSERTL0(false, "v_GetFluxVector: This function is not valid "
                            "for the Base class");
        }

        virtual void v_GetFluxVector(const int i, const int j,
                            Array<OneD, Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&flux)
        {
            ASSERTL0(false, "v_GetqFluxVector: This function is not valid "
                            "for the Base class");
        }

        virtual void v_GetFluxVector(const int i, Array<OneD,
                            Array<OneD, NekDouble> >&physfield,
                            Array<OneD, Array<OneD, NekDouble> >&fluxX,
                            Array<OneD, Array<OneD, NekDouble> > &fluxY)
        {
            ASSERTL0(false, "v_GetFluxVector: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumericalFlux(
                            Array<OneD, Array<OneD, NekDouble> > &physfield,
                            Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            ASSERTL0(false, "v_NumericalFlux: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumericalFlux(
                            Array<OneD, Array<OneD, NekDouble> > &physfield,
                            Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                            Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
            ASSERTL0(false, "v_NumericalFlux: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumFluxforScalar(
                    Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
        {
            ASSERTL0(false, "v_NumFluxforScalar: This function is not valid "
                            "for the Base class");
        }

        virtual void v_NumFluxforVector(
                    Array<OneD, Array<OneD, NekDouble> > &ufield,
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                    Array<OneD, Array<OneD, NekDouble > >  &qflux)
        {
            ASSERTL0(false, "v_NumFluxforVector: This function is not valid "
                            "for the Base class");
        }

        virtual NekDouble v_AdvectionSphere(const NekDouble x0j,
                            const NekDouble x1j, const NekDouble x2j,
                            const NekDouble time)
        {
            ASSERTL0(false, "v_AdvectionSphere: This function is not valid "
                            "for the Base class");
            return 0.0;
        }

        virtual NekDouble v_Morphogenesis(const int field, const NekDouble x0j,
                            const NekDouble x1j, const NekDouble x2j,
                            const NekDouble time)
        {
            ASSERTL0(false, "v_Morphogenesis: This function is not valid "
                            "for the Base class");
            return 0.0;
        }
    };

    /// Pointer to an ADRBase object.
    typedef boost::shared_ptr<ADRBase> ADRBaseSharedPtr;

} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

/**
* $Log: ADRBase.h,v $
* Revision 1.25  2010/02/26 13:52:47  cantwell
* Tested and fixed where necessary Hex/Tet projection and differentiation in
*   StdRegions, and LocalRegions for regular and deformed (where applicable).
* Added SpatialData and SpatialParameters classes for managing spatiall-varying
*   data.
* Added TimingGeneralMatrixOp3D for timing operations on 3D geometries along
*   with some associated input meshes.
* Added 3D std and loc projection demos for tet and hex.
* Added 3D std and loc regression tests for tet and hex.
* Fixed bugs in regression tests in relation to reading OK files.
* Extended Elemental and Global optimisation parameters for 3D expansions.
* Added GNUPlot output format option.
* Updated ADR2DManifoldSolver to use spatially varying data.
* Added Barkley model to ADR2DManifoldSolver.
* Added 3D support to FldToVtk and XmlToVtk.
* Renamed History.{h,cpp} to HistoryPoints.{h,cpp}
*
* Revision 1.24  2010/02/02 13:53:26  cantwell
* Moved reading in of history data to separate SpatialDomains class.
* Updated AlievPanfilov demo to move history specification.
* Replaced FindNektar line and NEKTAR_BIN_DIR def in regressionTests
* CMakeLists.txt as this is required to locate the regression test execs.
*
* Revision 1.23  2010/01/27 15:55:57  cantwell
* Fixed incorrect ordering of history point data.
* Fixed parsing of session name when session filename contains multiple
*   full-stops.
* Removed extra empty composite entries in XML generation from Gmsh.
*
* Revision 1.22  2010/01/27 13:19:13  cantwell
* Added functions to write history/probe data during timestepping.
*
* Revision 1.21  2010/01/26 17:43:08  cantwell
* Updated CMakeLists.txt to build FitzHughNagumoSolver
* Added Aliev-Panfilov model to ADR2DManifoldSolver
*
* Revision 1.20  2010/01/02 04:32:14  sehunchun
* Add additional inline functions
*
* Revision 1.19  2009/12/14 17:59:07  cbiotto
* Adding writing tecplot file
*
* Revision 1.18  2009/12/09 12:36:36  cbiotto
* Update for regression test
*
* Revision 1.17  2009/11/22 19:36:47  bnelson
* Fixed windows compile error.
*
* Revision 1.16  2009/11/15 18:40:57  sehunchun
* Add GetExpSize(), GetPhys_Offset and GetTotPoints(int n)
*
* Revision 1.15  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.14  2009/10/07 16:52:20  cbiotto
* Updating WriteFld function
*
* Revision 1.13  2009/09/07 11:21:58  sherwin
* Updates related to Navier-Stokes solver
*
* Revision 1.12  2009/08/14 09:29:13  cbiotto
* Add WriteVar function
*
* Revision 1.11  2009/07/23 05:23:21  sehunchun
* WeakDiffusion operator is updated
*
* Revision 1.10  2009/07/11 23:39:23  sehunchun
* Move uncommon functions to each solvers
*
* Revision 1.9  2009/07/09 21:29:13  sehunchun
* VS: ----------------------------------------------------------------------
* Add SetUpSurfaceNormal function..
*
* Revision 1.8  2009/07/02 15:57:36  sehunchun
* "ReadBoundaryCondition" options with extenstion to 2D geometry imbedded in 3D
*
* Revision 1.7  2009/04/29 20:45:55  sherwin
* Update for new definition of enum
*
* Revision 1.6  2009/04/27 21:37:14  sherwin
* Updated to dump .fld and .chk file in compressed coefficient format
*
* Revision 1.5  2009/03/10 23:37:14  claes
* Updated the ShallowWaterSolver to work with the general timestepping scheme
*
* Revision 1.4  2009/02/28 22:00:38  sehunchun
*  Explicit Diffusion solver is added
*
* Revision 1.3  2009/02/03 14:33:44  pvos
* Modifications for solvers with time-dependent dirichlet BC's
*
* Revision 1.2  2009/02/02 16:10:16  claes
* Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working
*
* Revision 1.1  2009/01/13 10:59:32  pvos
* added new solvers file
*
* Revision 1.6  2009/01/06 21:11:03  sherwin
* Updates for Virtual ExpList calls
*
* Revision 1.5  2008/11/17 08:10:07  claes
* Removed functions that were no longer used after the solver library was restructured
*
* Revision 1.4  2008/10/31 10:50:10  pvos
* Restructured directory and CMakeFiles
*
* Revision 1.3  2008/10/29 22:51:07  sherwin
* Updates for const correctness and ODEforcing
*
* Revision 1.2  2008/10/19 15:59:20  sherwin
* Added Summary method
*
* Revision 1.1  2008/10/16 15:25:45  sherwin
* Working verion of restructured AdvectionDiffusionReactionSolver
*

**/
