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

#include <ADRSolver/EquationSystem.h>
#include <IncNavierStokesSolver/EquationSystems/AdvectionTerm.h>
#include <IncNavierStokesSolver/EquationSystems/LinearisedAdvection.h>
#include <IncNavierStokesSolver/EquationSystems/NavierStokesAdvection.h>
#include <IncNavierStokesSolver/EquationSystems/AdjointAdvection.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{     

    enum EquationType
    {
        eNoEquationType,
        eSteadyStokes,
        eSteadyOseen,
        eUnsteadyStokes,
        eUnsteadyNavierStokes,
        eEquationTypeSize
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] = 
    {
        "NoType",
        "SteadyStokes",
        "SteadyOseen",
        "UnsteadyStokes",
        "UnsteadyNavierStokes"
    };


    enum AdvectionForm
    {
        eNoAdvectionForm,
        eConvective,
        eNonConservative,
        eLinearised,
		eAdjoint,
		eAdvectionFormSize
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kAdvectionFormStr[] = 
    {
        "NoType",
        "Convective",
        "NonConservative",
        "Linearised",
		"Adjoint"
    };
	
    /**
     * \brief This class is the base class for Navier Stokes problems
     *
     */
    
    class IncNavierStokes: public EquationSystem
    {
    public:           

        // Destructor
        virtual ~IncNavierStokes();

    protected: 

        /// Advection term
        AdvectionTermSharedPtr m_advObject;

        /// Number of fields to be convected; 
        int   m_nConvectiveFields;  

        /// int which identifies which components of m_fields contains the velocity (u,v,w);
        Array<OneD, int> m_velocity; 
 
        /// Pointer to field holding pressure field
        MultiRegions::ExpListSharedPtr m_pressure;  
        
        NekDouble     m_kinvis;        ///< Kinematic viscosity
        int           m_infosteps;     ///< dump info to stdout at steps time
        EquationType  m_equationType;  ///< equation type;
        AdvectionForm m_advectionForm; ///< Form of advection terms. 

        // Time integration classes
        LibUtilities::TimeIntegrationSchemeOperators m_integrationOps;
        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> m_integrationScheme;
        int m_intSteps;  ///< Number of time integration steps AND  Order of extrapolation for pressure boundary conditions.         


        /**
         * Constructor.
         */
        IncNavierStokes(LibUtilities::CommSharedPtr& pComm,
                LibUtilities::SessionReaderSharedPtr& pSession);

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }

        void AdvanceInTime(int nsteps);

        void EvaluateAdvectionTerms(const Array<OneD, 
                                    const Array<OneD, NekDouble> > &inarray, 
                                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                    Array<OneD, NekDouble> &wk = NullNekDouble1DArray);
		
        //time dependent boundary conditions updating
	
        void SetBoundaryConditions(NekDouble time);

        // Virtual functions
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

    private: 
    };
    
    typedef boost::shared_ptr<IncNavierStokes> IncNavierStokesSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H

/**
* $Log: IncNavierStokes.h,v $
* Revision 1.2  2010/01/28 15:17:05  abolis
* Time-Dependent boundary conditions
*
* Revision 1.1  2009/09/06 22:31:15  sherwin
* First working version of Navier-Stokes solver and input files
*
**/
