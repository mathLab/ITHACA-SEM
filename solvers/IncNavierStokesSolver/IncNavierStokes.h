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

#include <MultiRegions/DisContField2D.h>
#include <Auxiliary/ADRBase.h>

namespace Nektar
{     

    enum EquationType
    {
        eNoEquationType,
        eUnsteadyStokes,
        eUnsteadyNavierStokes,
        eEquationTypeSize
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] = 
    {
        "NoType",
        "UnsteadyStokes",
        "UnsteadyNavierStokes"
    };


    enum AdvectionForm
    {
        eNoAdvectionForm,
        eConvective,
        eNonConservative,
        eAdvectionFormSize
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kAdvectionFormStr[] = 
    {
        "NoType",
        "Convective",
        "NonConservative"
    };
	
    /**
     * \brief This class is the base class for Navier Stokes problems
     *
     */
    
    class IncNavierStokes: public ADRBase
    {
    public:           

        /**
         * Default constructor. 
         * 
         */ 
        IncNavierStokes();


        /**
         * Constructor.
         * \param 
         * 
         *
         */
        IncNavierStokes(string &fileStringName);

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }


        void VelocityCorrectionScheme(int nsteps);

        void Summary(std::ostream &out);


        // Wrapper functions in Velocity Correction Scheme
        void AdvanceInTime(int nsteps)
        {
            v_AdvanceInTime(nsteps);
        }

        void EvaluateAdvectionTerms(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                    Array<OneD, NekDouble> &wk = NullNekDouble1DArray);
		
		//time dependent boundary conditions updating
		
		void SetBoundaryConditions(NekDouble time);
		
		
    protected:

        int   m_nConvectiveFields;  /// Number of fields to be convected; 

        Array<OneD, int> m_velocity; ///< int which identifies which components of m_fields contains the velocity (u,v,w);

        MultiRegions::ExpListSharedPtr m_pressure;  ///< Pointer to field holding pressure field

        NekDouble     m_kinvis;                ///< Kinematic viscosity
        int           m_infosteps;             ///< dump info to stdout at steps time
		
    private: 
        EquationType  m_equationType;  ///< equation type;
        AdvectionForm m_advectionForm; ///< Form of advection terms. 


	//void SetBoundaryConditions(NekDouble time); 
				   

        // Virtual functions
        virtual void v_AdvanceInTime(int nsteps)
        {
            ASSERTL0(false,"This method is not defined in this class");
        }

    };
    
    typedef boost::shared_ptr<IncNavierStokes> IncNavierStokesSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H

/**
* $Log: IncNavierStokes.h,v $
* Revision 1.1  2009/09/06 22:31:15  sherwin
* First working version of Navier-Stokes solver and input files
*
**/
