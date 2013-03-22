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
        
        RiemannSolver::RiemannSolver() : m_requiresRotation(false),
                                         m_rotStorage      (3)
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
                
                if (m_rotStorage[0].num_elements() != nFields ||
                    m_rotStorage[0][0].num_elements() != nPts)
                {
                    for (int i = 0; i < 3; ++i)
                    {
                        m_rotStorage[i] =
                            Array<OneD, Array<OneD, NekDouble> >(nFields);
                        for (int j = 0; j < nFields; ++j)
                        {
                            m_rotStorage[i][j] = Array<OneD, NekDouble>(nPts);
                        }
                    }
                }
                
                rotateToNormal  (Fwd, m_rotStorage[0]);
                rotateToNormal  (Bwd, m_rotStorage[1]);
                v_Solve         (m_rotStorage[0], m_rotStorage[1],
                                 m_rotStorage[2]);
                rotateFromNormal(m_rotStorage[2], flux);
            }
            else
            {
                v_Solve(Fwd, Bwd, flux);
            }
        }

        /**
         * @brief Rotate velocity field to trace normal.
         * 
         * This function performs a rotation of the velocity components provided
         * in inarray so that the first component aligns with the trace normal
         * direction.
         * 
         * In 2D, this is accomplished through the transform:
         * 
         * \f[ (u_x, u_y) = (n_x u_x + n_y u_y, -n_x v_x + n_y v_y) \f]
         * 
         * In 3D, we generate a (non-unique) transformation using
         * RiemannSolver::fromToRotation.
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
                case 1:
                    ASSERTL0(false, "1D not implemented yet.");
                    break;

                case 2:
                {
                    const int vx = (int)velLoc[0];
                    const int vy = (int)velLoc[1];
                    const int nq = normals[0].num_elements();
                    
                    Vmath::Vmul (nq, inarray [vx], 1, normals [0],  1,
                                     outarray[vx], 1);
                    Vmath::Vvtvp(nq, inarray [vy], 1, normals [1],  1,
                                     outarray[vx], 1, outarray[vx], 1);
                    Vmath::Vmul (nq, inarray [vx], 1, normals [1],  1,
                                     outarray[vy], 1);
                    Vmath::Vvtvm(nq, inarray [vy], 1, normals [0],  1,
                                     outarray[vy], 1, outarray[vy], 1);
                    
                    for (int i = 0; i < inarray.num_elements(); ++i)
                    {
                        if (i == vx || i == vy)
                        {
                            continue;
                        }
                        
                        Vmath::Vcopy(inarray [i].num_elements(), inarray[i], 1,
                                     outarray[i], 1);
                    }
                    break;
                }
                
                case 3:
                {
                    const int vx = (int)velLoc[0];
                    const int vy = (int)velLoc[1];
                    const int vz = (int)velLoc[2];
                    const int nq = normals[0].num_elements();

                    // Generate matrices if they don't already exist.
                    if (m_rotMat.num_elements() == 0)
                    {
                        GenerateRotationMatrices();
                    }
                    
                    // Apply rotation matrices.
                    Vmath::Vvtvvtp(nq, inarray [vx], 1, m_rotMat[0],  1,
                                       inarray [vy], 1, m_rotMat[1],  1,
                                       outarray[vx], 1);
                    Vmath::Vvtvp  (nq, inarray [vz], 1, m_rotMat[2],  1,
                                       outarray[vx], 1, outarray[vx], 1);
                    Vmath::Vvtvvtp(nq, inarray [vx], 1, m_rotMat[3],  1,
                                       inarray [vy], 1, m_rotMat[4],  1,
                                       outarray[vy], 1);
                    Vmath::Vvtvp  (nq, inarray [vz], 1, m_rotMat[5],  1,
                                       outarray[vy], 1, outarray[vy], 1);
                    Vmath::Vvtvvtp(nq, inarray [vx], 1, m_rotMat[6],  1,
                                       inarray [vy], 1, m_rotMat[7],  1,
                                       outarray[vz], 1);
                    Vmath::Vvtvp  (nq, inarray [vz], 1, m_rotMat[8],  1,
                                       outarray[vz], 1, outarray[vz], 1);

                    for (int i = 0; i < inarray.num_elements(); ++i)
                    {
                        if (i == vx || i == vy || i == vz)
                        {
                            continue;
                        }

                        Vmath::Vcopy(inarray [i].num_elements(), inarray[i], 1,
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
                    const int vx = (int)velLoc[0];
                    const int vy = (int)velLoc[1];
                    const int nq = normals[0].num_elements();

                    Vmath::Vmul (nq, inarray [vy], 1, normals [1],  1,
                                     outarray[vx], 1);
                    Vmath::Vvtvm(nq, inarray [vx], 1, normals [0],  1,
                                     outarray[vx], 1, outarray[vx], 1);
                    Vmath::Vmul (nq, inarray [vx], 1, normals [1],  1,
                                     outarray[vy], 1);
                    Vmath::Vvtvp(nq, inarray [vy], 1, normals [0],  1,
                                     outarray[vy], 1, outarray[vy], 1);

                    for (int i = 0; i < inarray.num_elements(); ++i)
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
                    const int vx = (int)velLoc[0];
                    const int vy = (int)velLoc[1];
                    const int vz = (int)velLoc[2];
                    const int nq = normals[0].num_elements();
                    
                    Vmath::Vvtvvtp(nq, inarray [vx], 1, m_rotMat[0],  1,
                                       inarray [vy], 1, m_rotMat[3],  1,
                                       outarray[vx], 1);
                    Vmath::Vvtvp  (nq, inarray [vz], 1, m_rotMat[6],  1,
                                       outarray[vx], 1, outarray[vx], 1);
                    Vmath::Vvtvvtp(nq, inarray [vx], 1, m_rotMat[1],  1,
                                       inarray [vy], 1, m_rotMat[4],  1,
                                       outarray[vy], 1);
                    Vmath::Vvtvp  (nq, inarray [vz], 1, m_rotMat[7],  1,
                                       outarray[vy], 1, outarray[vy], 1);
                    Vmath::Vvtvvtp(nq, inarray [vx], 1, m_rotMat[2],  1,
                                       inarray [vy], 1, m_rotMat[5],  1,
                                       outarray[vz], 1);
                    Vmath::Vvtvp  (nq, inarray [vz], 1, m_rotMat[8],  1,
                                       outarray[vz], 1, outarray[vz], 1);

                    for (int i = 0; i < inarray.num_elements(); ++i)
                    {
                        if (i == vx || i == vy || i == vz)
                        {
                            continue;
                        }
                        
                        Vmath::Vcopy(inarray [i].num_elements(), inarray[i], 1,
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
            NekDouble tmp[9];
            const int nq = normals[0].num_elements();
            int i, j;
            xdir[0] = 1.0;

            // Allocate storage for rotation matrices.
            m_rotMat = Array<OneD, Array<OneD, NekDouble> >(9);
            
            for (i = 0; i < 9; ++i)
            {
                m_rotMat[i] = Array<OneD, NekDouble>(nq);
            }

            for (i = 0; i < normals[0].num_elements(); ++i)
            {
                // Generate matrix which takes us from (1,0,0) vector to trace
                // normal.
                tn[0] = normals[0][i];
                tn[1] = normals[1][i];
                tn[2] = normals[2][i];
                FromToRotation(tn, xdir, tmp);
                
                for (j = 0; j < 9; ++j)
                {
                    m_rotMat[j][i] = tmp[j];
                }
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
            NekDouble                    *mat)
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
                        mat[3*i+j] =  - c1 * u[i] * u[j]
                            - c2 * v[i] * v[j]
                            + c3 * v[i] * u[j];
                    }
                    mat[i+3*i] += 1.0;
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
                mat[0] = e + hvx * v[0];
                mat[1] = hvxy - v[2];
                mat[2] = hvxz + v[1];
                mat[3] = hvxy + v[2];
                mat[4] = e + h * v[1] * v[1];
                mat[5] = hvyz - v[0];
                mat[6] = hvxz - v[1];
                mat[7] = hvyz + v[0];
                mat[8] = e + hvz * v[2];
            }
        }
    }
}
