///////////////////////////////////////////////////////////////////////////////
//
// File BoussinesqWaterEquations.h
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
// Description: Boussinesq equation routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_AUXILIARY_BOUSSINESQEQUATIONS_H
#define NEKTAR_SOLVERS_AUXILIARY_BOUSSINESQEQUATIONS_H

#include <../solvers/Auxiliary/ShallowWaterEquations.h>

namespace Nektar
{
    /**
     * \brief This class is used for the development of solvers for the depth-integrated Shallow Water Equations (SWE) 
     * \f[ \left[ \begin{array}{l} h \\ hu \\ hv \end{array} \right]_t 
     * + \left[ \begin{array}{l} hu \\ hu^2 + 0.5 g h^2 \\ huv \end{array} \right]_{x_0}
     * + \left[ \begin{array}{l} hv \\ huv \\ hv^2 + 0.5 g h^2 \end{array} \right]_{x_1} = {\bf S} \f]
     *
     * where \f$ h({\bf x,t}) = \eta({\bf x},t) + d({\bf x})\f$ is the total water depth, 
     * \f$ \eta\f$ is the free-surface elevation measured from the still-water level (SWL) and
     * \f$ d \f$ is the SWL water depth.
     */
  
class BoussinesqEquations: public ShallowWaterEquations

    {
    public:
        /**
         * \brief 
         */ 
        BoussinesqEquations();
        
	BoussinesqEquations(const BoussinesqEquations &In);

        /**
         * Constructor.
         *
         * /param
         * /param 
         */
        BoussinesqEquations(SpatialDomains::MeshGraph2D &graph2D,
			    SpatialDomains::BoundaryConditions &bcs,
			    int &nVariables);

      
        /**
         * Computes 
         *
         */
      

	void SetBoundaryConditionsWaveCont(void);
	
	void WallBoundaryWaveCont(int bcRegion);
	
	void NumericalFluxWaveCont(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY);


	void SetBoundaryConditionsSolve(void);

	void WallBoundarySolve(int bcRegion);

	void WaveContSolve(Array<OneD, NekDouble> &fce, NekDouble lambda);

	
	void SetBoundaryConditionsContVariables(void);

	void WallBoundaryContVariables(int bcRegion);
	
	void NumericalFluxConsVariables(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY);

      
	void Madsen92SpatialTerms(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY);

	void Lambda20Spatial(BoussinesqEquations &Uf, Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY,
			     Array<OneD, NekDouble> &Ht);
	
	void FullyNonlinearSpatial(BoussinesqEquations &Uf, Array<OneD, NekDouble> &outX, 
				   Array<OneD, NekDouble> &outY, Array<OneD, NekDouble> &Ht);

	void NumericalFluxGradient(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY, int field_no);

	void NumericalFluxDivergence(Array<OneD, NekDouble> &outX, Array<OneD, NekDouble> &outY, int field_no_1, int field_no_2);

        /**
         * Returns the value of \f$ \alpha_1 \f$
         */
        inline NekDouble GetAlpha1(void)
        {
            return m_alpha1;
        }

        /**
         * Sets the value of \f$ \alpha_1 \f$
         */
        inline void SetAlpha1(NekDouble alpha1)
        {
            m_alpha1 = alpha1;
        }
	/**
         * Returns the value of \f$ \alpha_2 \f$
         */
	inline NekDouble GetAlpha2(void)
        {
            return m_alpha2;
        }

        /**
         * Sets the value of \f$ \alpha_2 \f$
         */
        inline void SetAlpha2(NekDouble alpha2)
        {
            m_alpha2 = alpha2;
        }

    
    protected:
        NekDouble m_alpha1;  ///< enhancement parameter for dispersion
	NekDouble m_alpha2;  ///< enhancement parameter for shoaling
      
    
    private:
               
    };
    
  typedef boost::shared_ptr<BoussinesqEquations>   BoussinesqEquationsSharedPtr;
    
} //end of namespace

#endif // NEKTAR_SOLVERS_AUXILIARY_BOUSSINESQEQUATIONS_H

