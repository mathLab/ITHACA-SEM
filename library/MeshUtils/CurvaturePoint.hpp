////////////////////////////////////////////////////////////////////////////////
//
//  File: curavturepoint.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_CURAVTUREPOINT_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_CURAVTUREPOINT_H

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
        namespace MeshUtils {
            
            class CurvaturePoint {
                
                public:
                friend class MemoryManager<CurvaturePoint>;
                
                CurvaturePoint(const NekDouble &x,const NekDouble &y,
                               const NekDouble &z,const NekDouble &R,
                               const NekDouble &Nx,const NekDouble &Ny,
                               const NekDouble &Nz) :
                                                    m_x(x),m_y(y),m_z(z),
                                                    m_nx(Nx),m_ny(Ny),m_nz(Nz),
                                                    m_radius(R)
                {
                    m_valid = true;
                }
                
                CurvaturePoint(const NekDouble &x,const NekDouble &y,
                               const NekDouble &z,const NekDouble &Nx,
                               const NekDouble &Ny,const NekDouble &Nz) :
                                                    m_x(x),m_y(y),m_z(z),
                                                    m_nx(Nx),m_ny(Ny),m_nz(Nz)
                {
                    m_delta = -1;
                    m_valid = false;
                }
                
                void Process(const NekDouble &min,
                                                  const NekDouble &max,
                                                  const NekDouble &eps)
                {
                    if(m_valid)
                    {
                        m_delta = 2.0*m_radius*sqrt(eps*(2.0-eps));
                        
                        if(m_delta>max)
                        {
                            m_delta = max;
                        }
                        if(m_delta<min)
                        {
                            m_delta = min;
                        }
                    }
                }
                
                bool IsValid(){return m_valid;}
                NekDouble GetDelta()
                {
                    if(m_valid)
                    {
                        return m_delta;
                    }
                    else
                    {
                        return -1;
                    }
                }
                NekDouble X(){return m_x;}
                NekDouble Y(){return m_y;}
                NekDouble Z(){return m_z;}
                void GetNormal(NekDouble &x, NekDouble &y, NekDouble &z)
                {
                    x = m_nx;
                    y = m_ny;
                    z = m_nz;
                }
                
                private:
                
                
                NekDouble m_x,m_y,m_z;
                NekDouble m_nx, m_ny, m_nz;
                NekDouble m_radius;
                NekDouble m_delta;
                bool m_valid;
            };
            
            typedef boost::shared_ptr<CurvaturePoint> CurvaturePointSharedPtr;
            
        }
}

#endif
