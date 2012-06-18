///////////////////////////////////////////////////////////////////////////////
//
// File: RiemannSolver.cpp
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
// Description: Abstract base class for Riemann solvers with factory.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{
    namespace SolverUtils
    {
        /**
         * Retrieves the singleton instance of the Riemann solver factory.
         */
        RiemannSolverFactory& GetRiemannSolverFactory()
        {
            typedef Loki::SingletonHolder<RiemannSolverFactory,
                                          Loki::CreateUsingNew,
                                          Loki::NoDestroy > Type;
            return Type::Instance();
        }
        
        /**
         * @class RiemannSolver
         * 
         * @brief The RiemannSolver class provides an abstract interface under
         * which solvers for various Riemann problems can be implemented.
         */
        
        RiemannSolver::RiemannSolver() : m_requiresRotation(false)
        {
            
        }
        
        /**
         * @brief Perform the Riemann solve given the forwards and backwards
         * spaces.
         * 
         * This routine calls the virtual function #v_Solve to perform the
         * Riemann solve. If the flag #m_requiresRotation is set, then the
         * velocity field is rotated to the normal direction to perform
         * dimensional splitting, and the resulting fluxes are rotated back to
         * the Cartesian directions before being returned.
         * 
         * @param Fwd   Forwards trace space.
         * @param Bwd   Backwards trace space.
         * @param flux  Resultant flux along trace space.
         */
        void RiemannSolver::Solve(
            const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
            const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
                  Array<OneD,       Array<OneD, NekDouble> > &flux)
        {
            if (m_requiresRotation)
            {
                int nFields = Fwd   .num_elements();
                int nPts    = Fwd[0].num_elements();
                Array<OneD, Array<OneD, NekDouble> > tmp1(nFields);
                Array<OneD, Array<OneD, NekDouble> > tmp2(nFields);
                Array<OneD, Array<OneD, NekDouble> > tmp3(nFields);
                
                for (int i = 0; i < nFields; ++i)
                {
                    tmp1[i] = Array<OneD, NekDouble>(nPts);
                    tmp2[i] = Array<OneD, NekDouble>(nPts);
                    tmp3[i] = Array<OneD, NekDouble>(nPts);
                }
                
                rotateToNormal(Fwd, tmp1);
                rotateToNormal(Bwd, tmp2);
                v_Solve(tmp1, tmp2, tmp3);
                rotateFromNormal(tmp3, flux);
            }
            else
            {
                v_Solve(Fwd, Bwd, flux);
            }
        }

        /**
         * @brief Rotate velocity field to trace normal.
         * 
         * This function performs a rotation of the triad of velocity components
         * provided in inarray so that the first component aligns with the trace
         * normal direction.
         * 
         * In 2D, this is accomplished through the transform:
         * 
         * \f[ (u_x, u_y) = (n_x u_x + n_y u_y, -n_x v_x + n_y v_y) \f]
         * 
         * @todo Implement for 1D and 3D.
         */
        void RiemannSolver::rotateToNormal(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray)
        {
            const Array<OneD, const Array<OneD, NekDouble> > &normals = 
                m_vectors["N"]();
            const Array<OneD, NekDouble> &velLoc = m_scalars["velLoc"]();
            
            switch(normals.num_elements())
            {
                case 2:
                {
                    int vx = (int)velLoc[0];
                    int vy = (int)velLoc[1];
                    int i;
                    
                    for (i = 0; i < normals[0].num_elements(); ++i)
                    {
                        double tmp1 =  inarray[vx][i]*normals[0][i] + 
                            inarray[vy][i]*normals[1][i];
                        double tmp2 = -inarray[vx][i]*normals[1][i] + 
                            inarray[vy][i]*normals[0][i];
                        outarray[vx][i] = tmp1;
                        outarray[vy][i] = tmp2;
                    }
                    
                    for (i = 0; i < inarray.num_elements(); ++i)
                    {
                        if (i == vx || i == vy)
                        {
                            continue;
                        }
                        
                        Vmath::Vcopy(inarray[i].num_elements(), inarray[i], 1, 
                                     outarray[i], 1);
                    }
                    break;
                }
                
                case 1:
                case 3:
                    ASSERTL0(false, "1D and 3D not implemented yet.");
                    break;
                default:
                    ASSERTL0(false, "Invalid space dimension.");
                    break;
            }
        }
        
        /**
         * @brief Rotate velocity field from trace normal.
         * 
         * This function performs a rotation of the triad of velocity components
         * provided in inarray so that the first component aligns with the
         * Cartesian components; it performs the inverse operation of
         * RiemannSolver::rotateToNormal.
         */
        void RiemannSolver::rotateFromNormal(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray)
        {
            const Array<OneD, const Array<OneD, NekDouble> > &normals = 
                m_vectors["N"]();
            const Array<OneD, NekDouble> &velLoc = m_scalars["velLoc"]();
            
            switch(normals.num_elements())
            {
                case 2:
                {
                    int vx = (int)velLoc[0];
                    int vy = (int)velLoc[1];
                    int i;

                    for (i = 0; i < normals[0].num_elements(); ++i)
                    {
                        double tmp1 = inarray[vx][i]*normals[0][i] - 
                            inarray[vy][i]*normals[1][i];
                        double tmp2 = inarray[vx][i]*normals[1][i] + 
                            inarray[vy][i]*normals[0][i];
                        outarray[vx][i] = tmp1;
                        outarray[vy][i] = tmp2;
                    }
                    
                    for (i = 0; i < inarray.num_elements(); ++i)
                    {
                        if (i == vx || i == vy)
                        {
                            continue;
                        }
                        
                        Vmath::Vcopy(inarray[i].num_elements(), inarray[i], 1, 
                                     outarray[i], 1);
                    }
                    break;
                }
                
                case 1:
                case 3:
                    ASSERTL0(false, "1D and 3D not implemented yet.");
                    break;
                default:
                    ASSERTL0(false, "Invalid space dimension.");
                    break;
            }
        }

        /**
         * @brief Determine whether a scalar has been defined in #m_scalars.
         * 
         * @param name  Scalar name.
         */
        bool RiemannSolver::CheckScalars(std::string name)
        {
            std::map<std::string, RSScalarFuncType>::iterator it = 
                m_scalars.find(name);
            
            return it != m_scalars.end();
        }

        /**
         * @brief Determine whether a vector has been defined in #m_vectors.
         * 
         * @param name  Vector name.
         */
        bool RiemannSolver::CheckVectors(std::string name)
        {
            std::map<std::string, RSVecFuncType>::iterator it = 
                m_vectors.find(name);
            
            return it != m_vectors.end();
        }

        /**
         * @brief Determine whether a parameter has been defined in #m_params.
         * 
         * @param name  Parameter name.
         */
        bool RiemannSolver::CheckParams(std::string name)
        {
            std::map<std::string, RSParamFuncType>::iterator it = 
                m_params.find(name);
            
            return it != m_params.end();
        }
    }
}
