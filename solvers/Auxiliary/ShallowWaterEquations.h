///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterEquations.h
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
// Description: shallow water equation routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_AUXILIARY_SHALLOWWATEREQUATIONS_H
#define NEKTAR_SOLVERS_AUXILIARY_SHALLOWWATEREQUATIONS_H

#include <../solvers/Auxiliary/AdvectionDiffusionReaction.h>

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
    
class ShallowWaterEquations: public AdvectionDiffusionReaction
    {
    public:
        /**
         * \brief 
         */ 
        ShallowWaterEquations();
        
	ShallowWaterEquations(ShallowWaterEquations &In);

        /**
         * Constructor.
         *
         * /param
         * /param 
         */
        ShallowWaterEquations(SpatialDomains::MeshGraph2D &graph2D,
                              SpatialDomains::BoundaryConditions &bcs,
			      int &nVariables);

        void ConservativeToPrimitive(void);
	void PrimitiveToConservative(void);

        /**
         * Computes 
         *
         */
        void GetFluxVector(Array<OneD, Array<OneD, NekDouble> >&FX,
                           Array<OneD, Array<OneD, NekDouble> >&FY);
	
	void GetFluxVectorPrimitive(Array<OneD, Array<OneD, NekDouble> >&FX,
				    Array<OneD, Array<OneD, NekDouble> >&FY);
	
	void GetFluxVectorPrimitiveLinear(Array<OneD, Array<OneD, NekDouble> >&FX,
					  Array<OneD, Array<OneD, NekDouble> >&FY);
        /**
         * This function evaluates the numerical flux term \f$ \hat{\bf F} \cdot {\bf n} \f$
         * subject to boundary conditions.
         *
         * \param upX. An array of size where on output the values of  
         * \param upY. An array of size where on output ...
         */
        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &upX, 
                           Array<OneD, Array<OneD, NekDouble> > &upY);
	void NumericalFluxPrimitiveLinear(Array<OneD, Array<OneD, NekDouble> > &upX, 
			   Array<OneD, Array<OneD, NekDouble> > &upY);
	

        /**
         * Solves the Riemann problem by the HLLC approximated solver, see Toro (2002)
         */
        void RiemannSolver(NekDouble hL,NekDouble huL,NekDouble hvL,NekDouble hR,NekDouble huR, 
                           NekDouble hvR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux );
      
        /** 
         * This function applies any user defined boundary conditions in the boundary conditions file.
         * Presently the following user defined BC are implemented:
         *
         * Wall: no-flux slip condition. Implemented through \f$ h^{bc} = h^{in}\,, 
         * \quad {\bf u}^{bc}\cdot {\bf n} = - \left({\bf u}^{in} \cdot {\bf n}\right)\,,
         * \quad {\bf u}^{bc}\cdot {\bf t} =   {\bf u}^{in} \cdot {\bf t} \,.\f$
         */
        void SetBoundaryConditions(void);
        void WallBoundary(int bcRegion);

        /**
         * Returns the value of \f$ g \f$
         */
        inline NekDouble GetGravity(void)
        {
            return m_g;
        }

        /**
         * Sets the value of \f$ g \f$
         */
        inline void SetGravity(NekDouble g)
        {
            m_g = g;
        }

        inline NekDouble GetConstantDepth(void)
        {
	  return m_d;
        }

        /**
         * Sets the value of \f$ g \f$
         */
        inline void SetConstantDepth(NekDouble d)
        {
            m_d = d;
        }

        void SetCoriolis(SpatialDomains::BoundaryConditions &bcs);
      
        inline Array<OneD, NekDouble> &GetCoriolis(void)
	{
            return m_coriolis;
	}

    protected:
        NekDouble m_g;  ///< acceleration of gravity

	NekDouble m_d; ///< constant water depth
      
        Array<OneD, NekDouble> m_coriolis;
      
        Array<OneD, NekDouble> m_depth;

        Array<OneD, NekDouble> m_friction;
    
    private:
               
    };
    
    typedef boost::shared_ptr<ShallowWaterEquations>   ShallowWaterEquationsSharedPtr;
    
} //end of namespace

#endif // NEKTAR_SOLVERS_AUXILIARY_SHALLOWWATEREQUATIONS_H

/**
* $Log: ShallowWaterEquations.h,v $
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
