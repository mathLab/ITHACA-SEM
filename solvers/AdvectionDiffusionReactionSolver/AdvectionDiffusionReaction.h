///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionDiffusionReaction.h
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

#ifndef NEKTAR_SOLVERS_ADVECTIONDIFFUSIONREACTION_ADVECTIONDIFFUSIONREACTION_H
#define NEKTAR_SOLVERS_ADVECTIONDIFFUSIONREACTION_ADVECTIONDIFFUSIONREACTION_H

#include <MultiRegions/DisContField2D.h>
#include <Auxiliary/ADRBase.h>

namespace Nektar
{     

    enum EquationType
    {
        eNoEquationType,
        eAdvection,
        eSteadyDiffusion,
        eSteadyDiffusionReaction,
        eLaplace,
        ePoisson,
        eHelmholtz,
        eEquationTypeSize
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] = 
    {
        "NoType",
        "Advection",
        "SteadyDiffusion",
        "SteadyDiffusionReaction",
        "Laplace",
        "Poisson",
        "Helmholtz"
    };

    /**
     * \brief This class is the base class for the development of solvers.
     *
     * It is basically a class handling vector valued fields where every field is
     * a DisContField2D class 
     */
    
    class AdvectionDiffusionReaction: public ADRBase
    {
    public:           

        /**
         * Default constructor. 
         * 
         */ 
        AdvectionDiffusionReaction();


        /**
         * Constructor.
         * /param 
         * 
         *
         */
        AdvectionDiffusionReaction(string &fileStringName);

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }


        void GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> >&physfield, Array<OneD, Array<OneD, NekDouble> >&flux);

        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &numflux);

	void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &numfluxX,
			   Array<OneD, Array<OneD, NekDouble> > &numfluxY);
        
        void ODEforcing(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, Array<OneD, Array<OneD, NekDouble> >&outarray, NekDouble time);

        void ExplicitlyIntegrateAdvection(int nsteps);

        void SolveHelmholtz(NekDouble lambda);

        void Summary(std::ostream &out);


    protected:

    private: 
        int m_infosteps;             ///< dump info to stdout at steps time
        EquationType m_equationType; ///< equation type;

        Array<OneD, Array<OneD, NekDouble> >  m_velocity;

        void EvaluateAdvectionVelocity();

	void SetBoundaryConditions(NekDouble time); 
				   
        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            GetFluxVector(i,physfield,flux);
        }
        
        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            NumericalFlux(physfield, numflux); 
        }
        
	virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numfluxX, Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
	  NumericalFlux(physfield, numfluxX, numfluxY); 
        }
      
    };
    
    typedef boost::shared_ptr<AdvectionDiffusionReaction> AdvectionDiffusionReactionSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_ADVECTIONDIFFUSIONREACTION_ADVECTIONDIFFUSIONREACTION_H

/**
* $Log: AdvectionDiffusionReaction.h,v $
* Revision 1.4  2009/01/06 21:10:34  sherwin
* Updates for virtual calls to IProductWRTBase and introduced reader to handle SOLVERINFO section to specify different solvers
*
* Revision 1.3  2008/11/17 08:20:14  claes
* Temporary fix for CG schemes. 1D CG working (but not for userdefined BC). 1D DG not working
*
* Revision 1.2  2008/11/12 12:12:26  pvos
* Time Integration update
*
* Revision 1.1  2008/10/31 10:50:10  pvos
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
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
