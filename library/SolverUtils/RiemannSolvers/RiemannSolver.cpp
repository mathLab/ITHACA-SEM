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

#define EPSILON 0.000001

#define CROSS(dest, v1, v2){                 \
          dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
          dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
          dest[2] = v1[0] * v2[1] - v1[1] * v2[0];}

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2){       \
          dest[0] = v1[0] - v2[0]; \
          dest[1] = v1[1] - v2[1]; \
          dest[2] = v1[2] - v2[2];}

using namespace std;

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
            const Array<OneD, const Array<OneD, NekDouble> > &normals = m_vectors["N"]();
            const Array<OneD, NekDouble> &velLoc = m_scalars["velLoc"]();
            
            switch(normals.num_elements())
            {
                case 1:
                    ASSERTL0(false, "1D not implemented yet.");
                    break;

                case 2:
                {
                    int vx = (int)velLoc[0];
                    int vy = (int)velLoc[1];
                    int i;
                    
                    for (i = 0; i < normals[0].num_elements(); ++i)
                    {
                        NekDouble tmp1 =  inarray[vx][i]*normals[0][i] + 
                            inarray[vy][i]*normals[1][i];
                        NekDouble tmp2 = -inarray[vx][i]*normals[1][i] + 
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
                
                case 3:
                {
                    int vx = (int)velLoc[0];
                    int vy = (int)velLoc[1];
                    int vz = (int)velLoc[2];
                    
                    // Generate matrices if they don't already exist.
                    if (m_rotMatrices.size() == 0)
                    {
                        GenerateRotationMatrices();
                    }
                    
                    // Apply rotation matrices.
                    for (int i = 0; i < normals[0].num_elements(); ++i)
                    {
                        DNekMatSharedPtr m = m_rotMatrices[i];
                        NekDouble tmp1 = inarray[vx][i]*(*m)(0,0) + 
                            inarray[vy][i]*(*m)(0,1) + inarray[vz][i]*(*m)(0,2);
                        NekDouble tmp2 = inarray[vx][i]*(*m)(1,0) + 
                            inarray[vy][i]*(*m)(1,1) + inarray[vz][i]*(*m)(1,2);
                        NekDouble tmp3 = inarray[vx][i]*(*m)(2,0) + 
                            inarray[vy][i]*(*m)(2,1) + inarray[vz][i]*(*m)(2,2);
                        outarray[vx][i] = tmp1;
                        outarray[vy][i] = tmp2;
                        outarray[vz][i] = tmp3;
                    }
                    
                    for (int i = 0; i < inarray.num_elements(); ++i)
                    {
                        if (i == vx || i == vy || i == vz)
                        {
                            continue;
                        }
                        
                        Vmath::Vcopy(inarray[i].num_elements(), inarray[i], 1, 
                                     outarray[i], 1);
                    }
                    break;
                }
                
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
                case 1:
                    ASSERTL0(false, "1D not implemented yet.");
                    break;
                    
                case 2:
                {
                    int vx = (int)velLoc[0];
                    int vy = (int)velLoc[1];
                    int i;

                    for (i = 0; i < normals[0].num_elements(); ++i)
                    {
                        NekDouble tmp1 = inarray[vx][i]*normals[0][i] - 
                            inarray[vy][i]*normals[1][i];
                        NekDouble tmp2 = inarray[vx][i]*normals[1][i] + 
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
                
                case 3:
                {
                    int vx = (int)velLoc[0];
                    int vy = (int)velLoc[1];
                    int vz = (int)velLoc[2];
                    
                    // Generate matrices if they don't already exist.
                    if (m_rotMatrices.size() == 0)
                    {
                        GenerateRotationMatrices();
                    }
                    
                    // Apply rotation matrices.
                    for (int i = 0; i < normals[0].num_elements(); ++i)
                    {
                        DNekMatSharedPtr m = m_invRotMatrices[i];
                        NekDouble tmp1 = inarray[vx][i]*(*m)(0,0) + 
                            inarray[vy][i]*(*m)(0,1) + inarray[vz][i]*(*m)(0,2);
                        NekDouble tmp2 = inarray[vx][i]*(*m)(1,0) + 
                            inarray[vy][i]*(*m)(1,1) + inarray[vz][i]*(*m)(1,2);
                        NekDouble tmp3 = inarray[vx][i]*(*m)(2,0) + 
                            inarray[vy][i]*(*m)(2,1) + inarray[vz][i]*(*m)(2,2);
                        outarray[vx][i] = tmp1;
                        outarray[vy][i] = tmp2;
                        outarray[vz][i] = tmp3;
                    }
                    
                    for (int i = 0; i < inarray.num_elements(); ++i)
                    {
                        if (i == vx || i == vy || i == vz)
                        {
                            continue;
                        }
                        
                        Vmath::Vcopy(inarray[i].num_elements(), inarray[i], 1, 
                                     outarray[i], 1);
                    }
                    break;
                }

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

        /**
         * @brief Generate rotation matrices for 3D expansions.
         */
        void RiemannSolver::GenerateRotationMatrices()
        {
            const Array<OneD, const Array<OneD, NekDouble> > &normals = 
                m_vectors["N"]();

            Array<OneD, NekDouble> xdir(3,0.0);
            Array<OneD, NekDouble> tn  (3);
            xdir[0] = 1.0;
            
            for (int i = 0; i < normals[0].num_elements(); ++i)
            {
                DNekMatSharedPtr mat = MemoryManager<DNekMat>::
                    AllocateSharedPtr(3,3);
                
                // Generate matrix which takes us from (1,0,0) vector to trace
                // normal.
                tn[0] = normals[0][i];
                tn[1] = normals[1][i];
                tn[2] = normals[2][i];
                FromToRotation(tn, xdir, mat);
                
                // Generate inverse rotation matrix.
                DNekMatSharedPtr mat2 = MemoryManager<DNekMat>::
                    AllocateSharedPtr(*mat);
                mat2->Invert();
                
                m_rotMatrices.   push_back(mat);
                m_invRotMatrices.push_back(mat2);
            }
        }

        /**
         * @brief A function for creating a rotation matrix that rotates a
         * vector @a from into another vector @a to.
         * 
         * Authors: Tomas Möller, John Hughes
         *          "Efficiently Building a Matrix to Rotate One Vector to
         *          Another" Journal of Graphics Tools, 4(4):1-4, 1999
         * 
         * @param from  Normalised 3-vector to rotate from.
         * @param to    Normalised 3-vector to rotate to.
         * @param out   Resulting 3x3 rotation matrix (row-major order).
         */
        void RiemannSolver::FromToRotation(
            Array<OneD, const NekDouble> &from,
            Array<OneD, const NekDouble> &to,
            DNekMatSharedPtr              mat)
        {
            NekDouble v[3];
            NekDouble e, h, f;
        
            CROSS(v, from, to);
            e = DOT(from, to);
            f = (e < 0)? -e:e;
            if (f > 1.0 - EPSILON)
            {
                NekDouble u[3], v[3];
                NekDouble x[3];
                NekDouble c1, c2, c3;
                int i, j;
            
                x[0] = (from[0] > 0.0)? from[0] : -from[0];
                x[1] = (from[1] > 0.0)? from[1] : -from[1];
                x[2] = (from[2] > 0.0)? from[2] : -from[2];
            
                if (x[0] < x[1])
                {
                    if (x[0] < x[2])
                    {
                        x[0] = 1.0; x[1] = x[2] = 0.0;
                    }
                    else
                    {
                        x[2] = 1.0; x[0] = x[1] = 0.0;
                    }
                }
                else
                {
                    if (x[1] < x[2])
                    {
                        x[1] = 1.0; x[0] = x[2] = 0.0;
                    }
                    else
                    {
                        x[2] = 1.0; x[0] = x[1] = 0.0;
                    }
                }
            
                u[0] = x[0] - from[0]; 
                u[1] = x[1] - from[1]; 
                u[2] = x[2] - from[2];
                v[0] = x[0] - to  [0];   
                v[1] = x[1] - to  [1];   
                v[2] = x[2] - to  [2];
            
                c1 = 2.0 / DOT(u, u);
                c2 = 2.0 / DOT(v, v);
                c3 = c1 * c2  * DOT(u, v);
            
                for (i = 0; i < 3; i++) {
                    for (j = 0; j < 3; j++) {
                        (*mat)(i,j) =  - c1 * u[i] * u[j]
                            - c2 * v[i] * v[j]
                            + c3 * v[i] * u[j];
                    }
                    (*mat)(i,i) += 1.0;
                }
            }
            else
            {
                NekDouble hvx, hvz, hvxy, hvxz, hvyz;
                h = 1.0/(1.0 + e);
                hvx = h * v[0];
                hvz = h * v[2];
                hvxy = hvx * v[1];
                hvxz = hvx * v[2];
                hvyz = hvz * v[1];
                (*mat)(0,0) = e + hvx * v[0];
                (*mat)(0,1) = hvxy - v[2];
                (*mat)(0,2) = hvxz + v[1];
            
                (*mat)(1,0) = hvxy + v[2];
                (*mat)(1,1) = e + h * v[1] * v[1];
                (*mat)(1,2) = hvyz - v[0];
            
                (*mat)(2,0) = hvxz - v[1];
                (*mat)(2,1) = hvyz + v[0];
                (*mat)(2,2) = e + hvz * v[2];
            }
        }
    }
}
