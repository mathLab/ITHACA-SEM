///////////////////////////////////////////////////////////////////////////////
//
// File FHN.h
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
// Description: Basic FHN definition in two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_FHN_FHN_H
#define NEKTAR_SOLVERS_FHN_FHN_H

#include <MultiRegions/DisContField2D.h>
#include <Auxiliary/ADRBase.h>


namespace Nektar
{     

    enum EquationType
    {
        eTestmodel,
         eFHN1961
    };
    
    // Keep this consistent with the enums in EquationType.
    const std::string kEquationTypeStr[] = 
    {
        "Testmodel",
        "FHN1961"
    };

    /**
     * \brief This class is the base class for the development of solvers.
     *
     * It is basically a class handling vector valued fields where every field is
     * a DisContField2D class 
     */
    
    class FHN: public ADRBase
    {
    public:           

        /**
         * Default constructor. 
         * 
         */ 
        FHN();


        /**
         * Constructor.
         * /param 
         * 
         *
         */
        FHN(string &fileStringName);

        EquationType GetEquationType(void)
        {
            return m_equationType;
        }	   

        void ODErhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, 
                          Array<OneD,        Array<OneD, NekDouble> >&outarray, 
                    const NekDouble time);

        void ODETest_Reaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time);

        void ODEFHN_Reaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
                               const NekDouble time);
					
        void ODEhelmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			   Array<OneD, Array<OneD, NekDouble> >&outarray,
			   const NekDouble time, 
                           const NekDouble lambda);

        void ODEFHN_helmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			      Array<OneD, Array<OneD, NekDouble> >&outarray,
			      const NekDouble time, 
                              const NekDouble lambda);

        void GeneralTimeIntegration(int nsteps, 
		                   LibUtilities::TimeIntegrationMethod IntMethod,
				   LibUtilities::TimeIntegrationSchemeOperators ode);

        void Summary(std::ostream &out);

    protected:

    private: 
        int m_infosteps;             ///< dump info to stdout at steps time
        EquationType m_equationType; ///< equation type;
        NekDouble m_epsilon;    /// constant epsilon

        Array<OneD, Array<OneD, NekDouble> >  m_velocity;

        void EvaluateAdvectionVelocity();

	void SetBoundaryConditions(NekDouble time); 
				   
        virtual void v_GetFluxVector(const int i,
                                     Array<OneD, Array<OneD, NekDouble> > &physfield, 
                                     Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            GetFluxVector(i,physfield,flux);
        }
		
	virtual void v_GetFluxVector(const int i, 
                                         const int j, 
                                         Array<OneD, Array<OneD, NekDouble> > &physfield, 
                                         Array<OneD, Array<OneD, NekDouble> > &flux)
        {
            GetFluxVector(i,j,physfield,flux);
        }
        
        virtual void v_NumericalFlux(Array<OneD, 
                                     Array<OneD, NekDouble> > &physfield, 
                                     Array<OneD, Array<OneD, NekDouble> > &numflux)
        {
            NumericalFlux(physfield, numflux); 
        }
        
	virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
                                     Array<OneD, Array<OneD, NekDouble> > &numfluxY )
        {
	    NumericalFlux(physfield, numfluxX, numfluxY); 
        }
		
	virtual void v_NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield, 
				      Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
	{
	    NumFluxforDiff(ufield, uflux);
	}
						   
	virtual void v_NumFluxforDiff(Array<OneD, Array<OneD, NekDouble> > &ufield,
	                              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
				      Array<OneD, Array<OneD, NekDouble> >  &qflux)
	{
	    NumFluxforDiff(ufield, qfield, qflux);
	}
      
    };
    
    typedef boost::shared_ptr<FHN> FHNSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_FHN_FHN_H

/**
* $Log: FHN.h,v $
* Revision 1.1  2009/03/06 16:02:55  sehunchun
* FitzHugh-Nagumo modeling
*

* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Sehun' FHN Solver
*
**/
