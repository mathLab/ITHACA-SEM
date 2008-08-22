///////////////////////////////////////////////////////////////////////////////
//
// File EulerEquations.h
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
// Description: two-dimensional Euler DG class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_AUXILIARY_EULEREQUATIONS_H
#define NEKTAR_SOLVERS_AUXILIARY_EULEREQUATIONS_H

#include <../solvers/Auxiliary/AdvectionDiffusionReaction.h>

namespace Nektar
{ 
    /**
     * \brief This class is used for the development of solvers for the compressible Euler Equations 
     * \f[ \left[ \begin{array}{l} \rho \\ \rho u \\  \rho v \\ E \end{array} \right]_t 
     * + \left[ \begin{array}{l} \rho u \\ \rho u^2 + p \\  \rho uv \\ u \left( E + p \right) \end{array} \right]_{x_0}
     * + \left[ \begin{array}{l} \rho v \\ \rho uv \\ \rho v^2 + p \\ v \left( E + p\right) \end{array} \right]_{x_1} = {\bf 0} \f]
     *
     * where \f$ E = \rho \left( \frac{1}{2} (u^2+v^2) + e \right)\f$ is the total energy per unit area and 
     */
      
class EulerEquations: public AdvectionDiffusionReaction
    {
    public:
        /**
         * \brief 
         */ 
        EulerEquations();
      
        EulerEquations(SpatialDomains::MeshGraph2D &graph2D,
                       SpatialDomains::BoundaryConditions &bcs);

      
        void GetFluxVector(Array<OneD, Array<OneD, NekDouble> >&FX,
                           Array<OneD, Array<OneD, NekDouble> >&FY);
      
        /**
         *
         * \param upX 
         * \param upY
         */
        void NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &upX, 
                           Array<OneD, Array<OneD, NekDouble> > &upY);
      
        void RiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
                           NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
                           NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux );  
      
        /**
         * Takes no arguments.
         * This function extracts the values from the traces and
         * sets any user defined boundary types. The values are
         * then copied into the #m_phys array of the #m_bndConstraints
         */
        void SetBoundaryConditions(void);
        void WallBoundary(int bcRegion);
        void IsenTropicVortexBoundary(int bcRegion);

      
        inline NekDouble GetGamma(void)
        {
            return m_gamma;
        }

        inline void SetGamma(NekDouble gamma)
        {
            m_gamma = gamma;
        }

    protected:
      
    private:
        NekDouble m_gamma; 
      
        void GetIsenTropicVortex(NekDouble x, NekDouble y, NekDouble time, NekDouble &rho, 
                                 NekDouble &rhou, NekDouble &rhov, NekDouble &E);
    };
    
    typedef boost::shared_ptr<EulerEquations>   EulerEquationsSharedPtr;
    
} //end of namespace

#endif // NEKTAR_SOLVERS_AUXILIARY_EULEREQUATIONS_H

/**
* $Log: $
**/
