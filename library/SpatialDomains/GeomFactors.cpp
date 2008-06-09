////////////////////////////////////////////////////////////////////////////////
//
//  File: GeomFactors.cpp
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
//  Description:  
//
//
////////////////////////////////////////////////////////////////////////////////
#include "pchSpatialDomains.h"

#include <SpatialDomains/GeomFactors.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        GeomFactors::GeomFactors(void):m_gtype(eRegular),
                                       m_expdim(0), m_coordim(0)
        {
        }

        GeomFactors::GeomFactors(const GeomType gtype,
            const int expdim, const int coordim):m_gtype(gtype),
            m_expdim(expdim), m_coordim(coordim)
        {
        }

        /** \brief One dimensional geometric factors based on one, two or three
        dimensional coordinate description 

        - Calculate the local factors \f$ d \xi/d x_i \f$ 
        as a recipricol derivative \f$ d x_i/d\xi \f$

        - Calcuate the jacobian as \f$ J = \sqrt{\frac{d x_1}{d _\xi}.
        \frac{d x_1}{d _\xi} + \frac{d x_2}{d _\xi}.\frac{d x_2}{d _\xi}  + 
        \frac{d x_3}{d _\xi}\frac{d x_3}{d _\xi}} \f$

        **/

        GeomFactors::GeomFactors(const GeomType gtype, const int coordim, 
            const Array<OneD, const StdRegions::StdExpansion1DSharedPtr>  &Coords):
            m_gtype(gtype), m_coordim(coordim)
        {
            int        i,nquad;
            LibUtilities::PointsType  ptype;

            ASSERTL1(coordim <= 3, "Only understand up to 3 Coordinate for this routine");

            nquad = Coords[0]->GetNumPoints(0);
            m_expdim = nquad;
            ptype = Coords[0]->GetPointsType(0);
            NekDouble fac0,fac1;

            Array<OneD,NekDouble> der[3] = {Array<OneD, NekDouble>(nquad),
                                            Array<OneD, NekDouble>(nquad), 
                                            Array<OneD, NekDouble>(nquad)};

            // Calculate local derivatives usin g physical space storage
            for(i = 0; i < coordim; ++i)
            {
                ASSERTL2(Coords[i]->GetNumPoints(0) == nquad,
                    "Points order are different for each coordinate");
                ASSERTL2(Coords[i]->GetPointsType(0)  == ptype,
                    "Points type are different for each coordinate");
        
                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(), 
                    Coords[i]->UpdatePhys());

                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(), der[i]);
            }

            if(( m_gtype == eRegular)||
               ( m_gtype == eMovingRegular))
            {
                m_jac  = Array<OneD, NekDouble>(1,0.0);
                m_gmat = Array<TwoD, NekDouble>(coordim,1,0.0);
                m_normals = Array<TwoD,NekDouble>(2,coordim,0.0);

                fac0 = fac1 = 0.0;
                for(i = 0; i < coordim; ++i)
                {
                    m_gmat[i][0] = 1.0/der[i][0];
                    m_jac[0] += der[i][0]*der[i][0];

                    m_normals[0][i] = -m_gmat[i][0];
                    m_normals[1][i] =  m_gmat[i][0];

                    fac0 += m_gmat[i][0]*m_gmat[i][0];
                }
                m_jac[0] = sqrt(m_jac[0]);
                
                fac0 = fac1 = sqrt(fac0);
            }
            else
            {
                m_jac     = Array<OneD, NekDouble>(nquad);
                m_gmat    = Array<TwoD, NekDouble>(coordim,nquad);
                m_normals = Array<TwoD,NekDouble>(2,coordim,0.0);

                // invert local derivative for gmat;
                fac0 = fac1 = 0.0;
                for(i = 0; i < coordim; ++i)
                {
                    Vmath::Sdiv(nquad,1.0,&der[i][0],1,&m_gmat[i][0],1);
                    Vmath::Vmul(nquad,&der[i][0],1,&der[i][0],1,&m_jac[0],1);

                    m_normals[0][i] = m_gmat[i][0];
                    m_normals[1][i] = -m_gmat[i][nquad-1];

                    fac0 += m_gmat[i][0]*m_gmat[i][0];
                    fac1 += m_gmat[i][nquad-1]*m_gmat[i][nquad-1];
                }
                Vmath::Vsqrt(nquad,&m_jac[0],1,&m_jac[0],1);

                fac0 = sqrt(fac0);
                fac1 = sqrt(fac1);

            } 

            // normalise normals 
            for(i = 0; i  < coordim; ++i)
            {
                m_normals[0][i] /= fac0;
                m_normals[1][i] /= fac1;
            }

       }


        /**
        \brief Two dimensional geometric factors based on two or three
        dimensional coordinate description

        The geometric factors are evaluated by considering the
        mapping to a coordinate system based on the local tangent
        vectors (which are the local derivatives of the global
        coordinates) and the normal \f$ \bf g \f$ to these two
        tangent vectors. We therefore use the 3 x 3 relationships but
        assume that \f$ \partial x_1/\partial \xi_3 = { g_1},\,
        \partial x_2/\partial \xi_3 = { g_2},\, \partial x_3/\partial
        \xi_3 = { g_3} \f$ i.e.

        \f$ {\bf g }= \left [ \begin{array}{c} g_1 \\ g_2 \\ g_3 \end{array} 
        \right ] = \frac{\partial {\bf x}}{\partial \xi_1} \times 
        \frac{\partial {\bf x}}{\partial \xi_2} =   \left [ \begin{array}{c}
        \frac{\partial x_2}{\partial \xi_1}   \frac{\partial x_3}{\partial \xi_2} - 
        \frac{\partial x_3}{\partial \xi_1}   \frac{\partial x_2}{\partial \xi_2} \\
        \frac{\partial x_3}{\partial \xi_1}   \frac{\partial x_1}{\partial \xi_2} - 
        \frac{\partial x_1}{\partial \xi_1}   \frac{\partial x_3}{\partial \xi_2} \\
        \frac{\partial x_1}{\partial \xi_1}   \frac{\partial x_2}{\partial \xi_2} - 
        \frac{\partial x_2}{\partial \xi_1}   \frac{\partial x_1}{\partial \xi_2}
        \end{array} \right ] \f$
    
        The geometric factors are then given by:

        \f$ \begin{array}{cc}
        \frac{\partial \xi_1}{\partial x_1} = \frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_2} {g_3} - \frac{\partial x_3}{\partial \xi_2} {g_2} \right ) & 
        \frac{\partial \xi_1}{\partial x_2} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} {g_3} - \frac{\partial x_3}{\partial \xi_2} {g_1} \right )\\
        \frac{\partial \xi_1}{\partial x_3} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} {g_2} - \frac{\partial x_2}{\partial \xi_2} {g_1} \right )&
        \frac{\partial \xi_2}{\partial x_1} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_1} {g_3} - \frac{\partial x_3}{\partial \xi_1} {g_2} \right ) \\
        \frac{\partial \xi_2}{\partial x_2} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} {g_3} - \frac{\partial x_3}{\partial \xi_1} {g_1} \right ) &
        \frac{\partial \xi_2}{\partial x_3} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} {g_2} - \frac{\partial x_2}{\partial \xi_1} {g_1} \right )
        \end{array} \f$

        where

        \f$ J_{3D} = 
        {g_3} \left ( 
        \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_2} - 
        \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_2}{\partial \xi_1}\right )
        + {g_2}\left (  
        \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_1}  -
        \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2}\right )
        + {g_1} \left ( 
        \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2} -
        \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_1} 
        \right ) \f$

        and the two-dimensional surface Jacobian  is given by
        \f$ J = \left | \frac{\partial {\bf x}}{\partial \xi_1} \times 
        \frac{\partial {\bf x}}{\partial \xi_2} \right | = \sqrt{J_{3D}} \f$

        **/
        GeomFactors::GeomFactors(const GeomType gtype, const int coordim, 
                                 const Array<OneD, const StdRegions::StdExpansion2DSharedPtr> &Coords)
        {
            int i,j,nquad0,nquad1,nqtot;
            LibUtilities::PointsType  ptype0, ptype1;

            ASSERTL1((coordim >= 2)&&(coordim <= 3),
                "Only understand up to three Coordinate and must have "
                "at least two coordinates for this routine");

            m_gtype = gtype;

            nquad0 = Coords[0]->GetNumPoints(0);
            nquad1 = Coords[0]->GetNumPoints(1);
            ptype0 = Coords[0]->GetPointsType (0);
            ptype1 = Coords[0]->GetPointsType (1);

            nqtot = nquad0*nquad1;

            // setup temp storage
            Array<OneD,NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot), 
                                           Array<OneD, NekDouble>(nqtot)};
            Array<OneD,NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot), 
                                           Array<OneD, NekDouble>(nqtot)};

            // Calculate local derivatives using physical space storage
            for(i = 0; i < coordim; ++i)
            {
                ASSERTL2(Coords[i]->GetNumPoints(0) == nquad0,
                    "Points order are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetNumPoints(1) == nquad1,
                    "Points order are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(0)  == ptype0,
                    "Points type are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointsType(1)  == ptype1,
                    "Points type are different for coordinate 1 ");

                Coords[i]->BwdTrans(Coords[i]->GetCoeffs(), 
                    Coords[i]->UpdatePhys());

                Coords[i]->StdPhysDeriv(Coords[i]->GetPhys(),d1[i],d2[i]);
            }

            if((m_gtype == eRegular)||(m_gtype == eMovingRegular))
            {
                m_jac  = Array<OneD, NekDouble>(1,0.0);
                m_gmat = Array<TwoD, NekDouble>(2*coordim,1,0.0);
                m_normals = Array<TwoD,NekDouble>(coordim,1,0.0);
                
                if(coordim == 2) // assume g = [0,0,1]
                {
                    m_jac[0] = d1[0][0]*d2[1][0] - d2[0][0]*d1[1][0];

                    ASSERTL1(m_jac[0] > 0, "2D Regular Jacobian is not positive");

                    m_gmat[0][0] =  d2[1][0]/m_jac[0]; // d xi_1/d x_1
                    m_gmat[1][0] = -d1[1][0]/m_jac[0]; // d xi_2/d x_1
                    m_gmat[2][0] = -d2[0][0]/m_jac[0]; // d xi_1/d x_2
                    m_gmat[3][0] =  d1[0][0]/m_jac[0]; // d xi_2/d x_2
                }
                else
                {
                    NekDouble g[3];
                    g[0] = d1[1][0]*d2[2][0] - d1[2][0]*d2[1][0];
                    g[1] = d1[2][0]*d2[0][0] - d1[0][0]*d2[2][0];
                    g[2] = d1[0][0]*d2[1][0] - d1[1][0]*d2[0][0];

                    m_jac[0] = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
                    ASSERTL1(m_jac[0] > 0, "Regular Jacobian is not positive");

                    m_gmat[0][0] =  (d2[1][0]*g[2] - d2[2][0]*g[1])/m_jac[0]; // d xi_1/d x_1
                    m_gmat[1][0] = -(d1[1][0]*g[2] - d1[2][0]*g[1])/m_jac[0]; // d xi_2/d x_1
                    m_gmat[2][0] = -(d2[0][0]*g[2] - d2[2][0]*g[0])/m_jac[0]; // d xi_1/d x_2
                    m_gmat[3][0] =  (d1[0][0]*g[2] - d1[2][0]*g[0])/m_jac[0]; // d xi_2/d x_2
                    m_gmat[4][0] =  (d2[0][0]*g[1] - d2[1][0]*g[0])/m_jac[0]; // d xi_1/d x_3
                    m_gmat[5][0] = -(d1[0][0]*g[1] - d1[1][0]*g[0])/m_jac[0]; // d xi_2/d x_3

                    m_jac[0] = sqrt(m_jac[0]);
                }


                // Set up normals
                switch(Coords[0]->DetExpansionType())
                {
                case StdRegions::eTriangle:
                    {
                        Array<OneD,NekDouble> fac(3,0.0);
                        m_normals = Array<TwoD,NekDouble>(3,coordim,0.0);
                        
                        for(i = 0; i < coordim; ++i)
                        {
                            m_normals[0][i] = -m_gmat[2*i+1][0];
                            fac[0] += m_normals[0][i]*m_normals[0][i];
                            m_normals[1][i] =  m_gmat[2*i+1][0] + m_gmat[2*i][0];
                            fac[1] += m_normals[1][i]*m_normals[1][i];
                            m_normals[2][i] = -m_gmat[2*i][0];
                            fac[2] += m_normals[2][i]*m_normals[2][i];
                        }

                        for(i = 0; i < coordim; ++i)
                        {
                            m_normals[0][i] /= sqrt(fac[0]);
                            m_normals[1][i] /= sqrt(fac[1]);
                            m_normals[2][i] /= sqrt(fac[2]);
                        }
                    }
                    break;
                case StdRegions::eQuadrilateral:
                    {
                        Array<OneD,NekDouble> fac(4,0.0);
                        m_normals = Array<TwoD,NekDouble>(4,coordim,0.0);
                        
                        for(i = 0; i < coordim; ++i)
                        {
                            m_normals[0][i] = -m_gmat[2*i+1][0];
                            fac[0] += m_normals[0][i]*m_normals[0][i];
                            m_normals[1][i] =  m_gmat[2*i][0];
                            fac[1] += m_normals[1][i]*m_normals[1][i];
                            m_normals[2][i] =  m_gmat[2*i+1][0];
                            fac[2] += m_normals[2][i]*m_normals[2][i];
                            m_normals[3][i] = -m_gmat[2*i][0];
                            fac[3] += m_normals[3][i]*m_normals[3][i];
                        }

                        for(i = 0; i < coordim; ++i)
                        {
                            m_normals[0][i] /= sqrt(fac[0]);
                            m_normals[1][i] /= sqrt(fac[1]);
                            m_normals[2][i] /= sqrt(fac[2]);
                            m_normals[3][i] /= sqrt(fac[3]);
                        }
                    }
                    break;
                }

            }  
            else
            {
                m_jac  = Array<OneD, NekDouble>(nqtot,0.0);
                m_gmat = Array<TwoD, NekDouble>(2*coordim,nqtot,0.0);

                if(coordim == 2) // assume g = [0,0,1]
                {
                    // set up Jacobian 
                    Vmath::Vmul (nqtot,&d2[0][0],1,&d1[1][0],1,&m_jac[0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&m_jac[0],1,
                                 &m_jac[0],1);

                    ASSERTL1(Vmath::Vmin(nqtot,&m_jac[0],1) > 0, "2D Deformed Jacobian is not positive");
                    
                    Vmath::Vdiv(nqtot,&d2[1][0],1,&m_jac[0],1,&m_gmat[0][0],1); // d xi_1/d x_1
                    Vmath::Vdiv(nqtot,&d1[1][0],1,&m_jac[0],1,&m_gmat[1][0],1); 
                    Vmath::Neg(nqtot,&m_gmat[1][0],1);                          // d xi_2/d x_1
                    Vmath::Vdiv(nqtot,&d2[0][0],1,&m_jac[0],1,&m_gmat[2][0],1); 
                    Vmath::Neg(nqtot,&m_gmat[2][0],1);                          // d xi_1/d x_2
                    Vmath::Vdiv(nqtot,&d1[0][0],1,&m_jac[0],1,&m_gmat[3][0],1); // d xi_2/d x_2
                }
                else
                {
                    Array<OneD,NekDouble> g[3] = {Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot),
                                                  Array<OneD, NekDouble>(nqtot)};
                    // g[0]
                    Vmath::Vmul (nqtot,&d1[2][0],1,&d2[1][0],1,&g[0][0],1);
                    Vmath::Vvtvm(nqtot,&d1[1][0],1,&d2[2][0],1,&g[0][0],1,&g[0][0],1);
                    //g[1]
                    Vmath::Vmul (nqtot,&d1[0][0],1,&d2[2][0],1,&g[1][0],1);
                    Vmath::Vvtvm(nqtot,&d1[2][0],1,&d2[0][0],1,&g[1][0],1,&g[1][0],1);
                    //g[2]
                    Vmath::Vmul (nqtot,&d1[1][0],1,&d2[0][0],1,&g[2][0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&d2[1][0],1,&g[2][0],1,&g[2][0],1);

                    // J_3D
                    Vmath::Vmul (nqtot,&g[0][0],1,&g[0][0],1,&m_jac[0],1);
                    Vmath::Vvtvp(nqtot,&g[1][0],1,&g[1][0],1,&m_jac[0],1,&m_jac[0],1);
                    Vmath::Vvtvp(nqtot,&g[2][0],1,&g[2][0],1,&m_jac[0],1,&m_jac[0],1);

                    // d xi_1/d x_1
                    Vmath::Vmul (nqtot,&d2[2][0],1,&g[1][0],1,&m_gmat[0][0],1);
                    Vmath::Vvtvm(nqtot,&d2[1][0],1,&g[2][0],1,&m_gmat[0][0],1,&m_gmat[0][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[0][0],1,&m_jac[0],1,&m_gmat[0][0],1);

                    // d xi_2/d x_1
                    Vmath::Vmul (nqtot,&d1[1][0],1,&g[2][0],1,&m_gmat[1][0],1);
                    Vmath::Vvtvm(nqtot,&d1[2][0],1,&g[1][0],1,&m_gmat[1][0],1,&m_gmat[1][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[1][0],1,&m_jac[0],1,&m_gmat[1][0],1);

                    // d xi_1/d x_2
                    Vmath::Vmul (nqtot,&d2[0][0],1,&g[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vvtvm(nqtot,&d2[2][0],1,&g[0][0],1,&m_gmat[2][0],1,&m_gmat[2][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[2][0],1,&m_jac[0],1,&m_gmat[2][0],1);

                    // d xi_2/d x_2
                    Vmath::Vmul (nqtot,&d1[2][0],1,&g[0][0],1,&m_gmat[3][0],1);
                    Vmath::Vvtvm(nqtot,&d1[0][0],1,&g[2][0],1,&m_gmat[3][0],1,&m_gmat[3][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[3][0],1,&m_jac[0],1,&m_gmat[3][0],1);

                    // d xi_1/d x_3
                    Vmath::Vmul (nqtot,&d2[1][0],1,&g[0][0],1,&m_gmat[4][0],1);
                    Vmath::Vvtvm(nqtot,&d2[0][0],1,&g[1][0],1,&m_gmat[4][0],1,&m_gmat[4][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[4][0],1,&m_jac[0],1,&m_gmat[4][0],1);

                    // d xi_2/d x_3
                    Vmath::Vmul (nqtot,&d1[0][0],1,&g[1][0],1,&m_gmat[5][0],1);
                    Vmath::Vvtvm(nqtot,&d1[1][0],1,&g[0][0],1,&m_gmat[5][0],1,&m_gmat[5][0],1);
                    Vmath::Vdiv(nqtot,&m_gmat[5][0],1,&m_jac[0],1,&m_gmat[5][0],1);

                    // J = sqrt(J_3D)
                    Vmath::Vsqrt(nqtot,&m_jac[0],1,&m_jac[0],1);
                }

                // Set up normals
                int nquad_m = max(nquad0,nquad1);                
                switch(Coords[0]->DetExpansionType())
                {
                case StdRegions::eTriangle:
                    {
                        m_normals = Array<TwoD,NekDouble>(3,coordim*nquad_m,0.0);

                        for(i = 0; i < coordim; ++i)
                        {
                            for(j = 0; j < nquad0; ++j)
                            {
                                m_normals[0][i*nquad0+j] = -m_gmat[2*i+1][j];
                            }

                            
                            for(j = 0; j < nquad1; ++j)
                            {
                                m_normals[1][i*nquad1+j] = -m_gmat[2*i][nquad0*j + nquad0-1] +
                                    m_gmat[2*i+1][nquad0*j + nquad0-1];
                                m_normals[3][i*nquad1+j] = -m_gmat[2*i][nquad0*j];
                            }
                        }

                        //normalise normal vectors
                        Array<OneD,NekDouble> norm(nquad_m,0.0);

                        //edges 0
                        for(i = 0; i < coordim; ++i)
                        {
                            Vmath::Vvtvp(nquad0,&m_normals[0][i*nquad0],1,
                                         &m_normals[0][i*nquad0],1,&norm[0],1,&norm[0],1);
                        }
                        Vmath::Vsqrt(nquad0,&norm[0],1,&norm[0],1);
                        Vmath::Sdiv(nquad0,1.0,&norm[0],1,&norm[0],1);
                        for(i = 0; i < coordim; ++i)
                        {
                            Vmath::Vmul(nquad0,&m_normals[0][i*nquad0],1,&norm[0],1,
                                        &m_normals[0][i*nquad0],1);
                        }
                        
                        // edge 1 + 2
                        for(j = 1; j < 3; ++j)
                        {
                            Vmath::Zero(nquad1,&norm[0],1);
                            for(i = 0; i < coordim; ++i)
                            {
                                Vmath::Vvtvp(nquad1,&m_normals[j][i*nquad1],1,
                                             &m_normals[j][i*nquad1],1,&norm[0],1,&norm[0],1);
                            }
                            Vmath::Vsqrt(nquad1,&norm[0],1,&norm[0],1);
                            Vmath::Sdiv(nquad1,1.0,&norm[0],1,&norm[0],1);
                            for(i = 0; i < coordim; ++i)
                            {
                                Vmath::Vmul(nquad1,&m_normals[j][i*nquad1],1,&norm[0],1,
                                            &m_normals[j][i*nquad1],1);
                            }
                        }
                    }
                    break;
                case StdRegions::eQuadrilateral:
                    {
                        m_normals = Array<TwoD,NekDouble>(4,coordim*nquad_m,0.0);
                        
                        for(i = 0; i < coordim; ++i)
                        {
                            for(j = 0; j < nquad0; ++j)
                            {
                                m_normals[0][i*nquad0+j] = -m_gmat[2*i+1][j];
                                m_normals[2][i*nquad0+j] =  m_gmat[2*i+1][nquad0*(nquad1-1)+j];
                            }
                            
                            for(j = 0; j < nquad1; ++j)
                            {
                                m_normals[1][i*nquad1+j] =  m_gmat[2*i][nquad0*j + nquad0-1];
                                m_normals[3][i*nquad1+j] = -m_gmat[2*i][nquad0*j];
                            }
                        }
                        
                        //normalise normal vectors
                        Array<OneD,NekDouble> norm(nquad_m,0.0);
                        //edges 0 + 2 
                        for(j = 0; j < 4; j += 2)
                        {
                            Vmath::Zero(nquad0,&norm[0],1);
                            for(i = 0; i < coordim; ++i)
                            {
                                Vmath::Vvtvp(nquad0,&m_normals[j][i*nquad0],1,
                                             &m_normals[j][i*nquad0],1,&norm[0],1,&norm[0],1);
                            }
                            Vmath::Vsqrt(nquad0,&norm[0],1,&norm[0],1);
                            Vmath::Sdiv(nquad0,1.0,&norm[0],1,&norm[0],1);
                            for(i = 0; i < coordim; ++i)
                            {
                                Vmath::Vmul(nquad0,&m_normals[j][i*nquad0],1,&norm[0],1,
                                        &m_normals[j][i*nquad0],1);
                            }
                        }
                        // edge 1 + 3
                        for(j = 1; j < 4; j += 2)
                        {
                            Vmath::Zero(nquad1,&norm[0],1);
                            for(i = 0; i < coordim; ++i)
                            {
                                Vmath::Vvtvp(nquad1,&m_normals[j][i*nquad1],1,
                                             &m_normals[j][i*nquad1],1,&norm[0],1,&norm[0],1);
                            }
                            Vmath::Vsqrt(nquad1,&norm[0],1,&norm[0],1);
                            Vmath::Sdiv(nquad1,1.0,&norm[0],1,&norm[0],1);
                            for(i = 0; i < coordim; ++i)
                            {
                                Vmath::Vmul(nquad1,&m_normals[j][i*nquad1],1,&norm[0],1,
                                            &m_normals[j][i*nquad1],1);
                            }
                        }
                        
                    }
                    break;
                    
                default:
                    break;
                }
            }
        }

        /**
        \brief Three dimensional geometric factors and Jacobian

        The geometric factors are related to the local derivatives of the
        global coordinates by:

        \f$ \begin{array}{lll}
        \frac{\partial \xi_1}{\partial x_1} = \frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_3} - \frac{\partial x_2}{\partial \xi_3} \frac{\partial x_3}{\partial \xi_2} \right ) & 
        \frac{\partial \xi_1}{\partial x_2} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_3}  - \frac{\partial x_3}{\partial \xi_2} \frac{\partial x_1}{\partial \xi_3}  \right )&
        \frac{\partial \xi_1}{\partial x_3} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_2} \frac{\partial x_2}{\partial \xi_3}  - \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_1}{\partial \xi_3}  \right )\\

        \frac{\partial \xi_2}{\partial x_1} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_3} - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_3}  \right ) &
        \frac{\partial \xi_2}{\partial x_2} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_3}  - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_3} \right ) &
        \frac{\partial \xi_2}{\partial x_3} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_3}  - \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_3}  \right )\\

        \frac{\partial \xi_3}{\partial x_1} = \frac{1}{J_{3D}} \left ( \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2} - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_3}  \right ) &
        \frac{\partial \xi_3}{\partial x_2} = -\frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2}  - \frac{\partial x_3}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_2} \right ) &
        \frac{\partial \xi_3}{\partial x_3} = \frac{1}{J_{3D}} \left ( \frac{\partial x_1}{\partial \xi_1} \frac{\partial x_2}{\partial \xi_2}  - \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_1}{\partial \xi_2}  \right )
        \end{array} \f$

        where

        \f$ J_{3D} = 
        \frac{\partial x_1}{\partial \xi_1} \left ( 
        \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_3} - 
        \frac{\partial x_2}{\partial \xi_3} \frac{\partial x_3}{\partial \xi_2}\right )
        - \frac{\partial x_1}{\partial \xi_2} \left (  
        \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_3}  -
        \frac{\partial x_2}{\partial \xi_3} \frac{\partial x_3}{\partial \xi_1}\right )
        + \frac{\partial x_1}{\partial \xi_3} \left ( 
        \frac{\partial x_2}{\partial \xi_1} \frac{\partial x_3}{\partial \xi_2} -
        \frac{\partial x_2}{\partial \xi_2} \frac{\partial x_3}{\partial \xi_1} 
        \right ) \f$

        **/

            GeomFactors(const GeomType gtype, const int coordim,
                        const Array<OneD, const StdRegions::StdExpansion3DSharedPtr> &Coords);
        {

            int        i,nquad0,nquad1,nquad2,nqtot;
            StdRegions::PointsType  ptype0, ptype1, ptype2;
            double      *tmp;

            m_gtype = gtype;
            m_gmat  = new double*[9];

            nquad0 = Coords[0]->GetNumPoints(0);
            nquad1 = Coords[0]->GetNumPoints(1);
            nquad2 = Coords[0]->GetNumPoints(2);
            ptype0 = Coords[0]->GetPointsType (0);
            ptype1 = Coords[0]->GetPointsType (1);
            ptype2 = Coords[0]->GetPointsType (2);

            nqtot = nquad0*nquad1*nquad2;

            Array<OneD,NekDouble> d1[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot), 
                                           Array<OneD, NekDouble>(nqtot)};
            Array<OneD,NekDouble> d2[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot), 
                                           Array<OneD, NekDouble>(nqtot)};
            Array<OneD,NekDouble> d3[3] = {Array<OneD, NekDouble>(nqtot),
                                           Array<OneD, NekDouble>(nqtot), 
                                           Array<OneD, NekDouble>(nqtot)};

            d1[0] = new double [3*nqtot];
            d1[1] = d1[0] + nqtot;
            d1[2] = d1[1] + nqtot;

            d2[0] = new double [3*nqtot];
            d2[1] = d2[0] + nqtot;
            d2[2] = d2[1] + nqtot;

            d3[0] = new double [3*nqtot];
            d3[1] = d3[0] + nqtot;
            d3[2] = d3[1] + nqtot;

            // Calculate local derivatives using physical space storage
            for(i = 0; i < 3; ++i)
            {
                ASSERTL2(Coords[i]->GetNumPoints(0) == nquad0,
                    "Points order are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetNumPoints(1) == nquad1,
                    "Points order are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetNumPoints(2) == nquad2,
                    "Points order are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(0)  == ptype0,
                    "Points type are different for coordinate 0 ");
                ASSERTL2(Coords[i]->GetPointsType(1)  == ptype1,
                    "Points type are different for coordinate 1 ");
                ASSERTL2(Coords[i]->GetPointsType(2)  == ptype2,
                    "Points type are different for coordinate 1 ");

                Coords[i]->StdPhysDeriv(d1[i],d2[i],d3[i]);
            }

            if((m_gtype == eRegular)|| (m_gtype == eMovingRegular))
            {
                m_jac  = Array<OneD, NekDouble>(1,0.0);
                m_gmat = Array<TwoD, NekDouble>(3*coordim,1,0.0);
                m_normals = Array<TwoD,NekDouble>(coordim,1,0.0);

                for(i = 1; i < 9; ++i)
                {
                    m_gmat[i] = m_gmat[i-1]+1;
                }

                m_jac[0] = d1[0][0]*(d2[1][0]*d3[2][0] - d3[1][0]*d2[2][0])
                    - d2[0][0]*(d1[1][0]*d3[2][0] - d3[1][0]*d1[2][0]) 
                    + d3[0][0]*(d1[1][0]*d2[2][0] - d2[1][0]*d1[2][0]);

                // d xi_1/d x_1
                m_gmat[0][0] =  (d2[1][0]*d3[2][0] - d2[2][0]*d2[2][0])/m_jac[0]; 
                // d xi_2/d x_1
                m_gmat[1][0] = -(d1[1][0]*d3[2][0] - d1[2][0]*d3[2][0])/m_jac[0]; 
                // d xi_3/d x_1
                m_gmat[2][0] =  (d1[1][0]*d3[1][0] - d1[2][0]*d3[1][0])/m_jac[0]; 

                // d xi_1/d x_2
                m_gmat[3][0] = -(d2[0][0]*d3[2][0] - d2[2][0]*d3[0][0])/m_jac[0]; 
                // d xi_2/d x_2
                m_gmat[4][0] =  (d1[0][0]*d3[2][0] - d1[2][0]*d3[0][0])/m_jac[0]; 
                // d xi_3/d x_2
                m_gmat[5][0] = -(d1[0][0]*d2[2][0] - d1[2][0]*d2[0][0])/m_jac[0]; 

                // d xi_1/d x_3
                m_gmat[6][0] =  (d2[0][0]*d3[1][0] - d2[1][0]*d3[0][0])/m_jac[0]; 
                // d xi_2/d x_3
                m_gmat[5][0] = -(d1[0][0]*d3[1][0] - d1[1][0]*d3[0][0])/m_jac[0]; 
                // d x7_3/d x_3
                m_gmat[8][0] =  (d1[0][0]*d2[1][0] - d1[1][0]*d2[0][0])/m_jac[0]; 

            }
            else
            {
                m_jac  = Array<OneD, NekDouble>(1,0.0);
                m_gmat = Array<TwoD, NekDouble>(3*coordim,1,0.0);

                for(i = 1; i < 9; ++i)
                {
                    m_gmat[i] = m_gmat[i-1] + nqtot;
                }

                tmp = new double [nqtot]; 

                // Jacobian
                Vmath::Vmul (nqtot,d3[1],1,d2[2],1,tmp,1);
                Vmath::Vvtvm(nqtot,d2[1],1,d3[2],1,tmp,1,tmp,1);
                Vmath::Vmul (nqtot,d1[0],1,tmp,1,m_jac,1);

                Vmath::Vmul (nqtot,d1[1],1,d3[2],1,tmp,1);
                Vmath::Vvtvm(nqtot,d3[1],1,d1[2],1,tmp,1,tmp,1);
                Vmath::Vvtvp(nqtot,d2[0],1,tmp,1,m_jac,1,m_jac,1);

                Vmath::Vmul (nqtot,d2[1],1,d1[2],1,tmp,1);
                Vmath::Vvtvm(nqtot,d1[1],1,d2[2],1,tmp,1,tmp,1);
                Vmath::Vvtvp(nqtot,d3[0],1,tmp,1,m_jac,1,m_jac,1);

                // d xi_1/d x_1
                Vmath::Vmul (nqtot,d2[2],1,d3[1],1,m_gmat[0],1);
                Vmath::Vvtvm(nqtot,d2[1],1,d3[2],1,m_gmat[0],1,m_gmat[0],1);
                Vmath::Vdiv(nqtot,m_gmat[0],1,m_jac,1,m_gmat[0],1);

                // d xi_2/d x_1
                Vmath::Vmul (nqtot,d1[1],1,d3[2],1,m_gmat[1],1);
                Vmath::Vvtvm(nqtot,d1[2],1,d3[1],1,m_gmat[1],1,m_gmat[1],1);
                Vmath::Vdiv(nqtot,m_gmat[1],1,m_jac,1,m_gmat[1],1);

                // d xi_3/d x_1
                Vmath::Vmul (nqtot,d1[1],1,d2[2],1,m_gmat[2],1);
                Vmath::Vvtvm(nqtot,d1[2],1,d3[1],1,m_gmat[2],1,m_gmat[2],1);
                Vmath::Vdiv(nqtot,m_gmat[2],1,m_jac,1,m_gmat[2],1);

                // d xi_1/d x_2
                Vmath::Vmul (nqtot,d2[0],1,d3[2],1,m_gmat[3],1);
                Vmath::Vvtvm(nqtot,d2[2],1,d3[0],1,m_gmat[3],1,m_gmat[3],1);
                Vmath::Vdiv(nqtot,m_gmat[3],1,m_jac,1,m_gmat[3],1);

                // d xi_2/d x_2
                Vmath::Vmul (nqtot,d1[2],1,d3[0],1,m_gmat[4],1);
                Vmath::Vvtvm(nqtot,d1[0],1,d3[2],1,m_gmat[4],1,m_gmat[4],1);
                Vmath::Vdiv(nqtot,m_gmat[4],1,m_jac,1,m_gmat[4],1);

                // d xi_3/d x_2
                Vmath::Vmul (nqtot,d1[0],1,d2[2],1,m_gmat[5],1);
                Vmath::Vvtvm(nqtot,d1[2],1,d3[0],1,m_gmat[5],1,m_gmat[5],1);
                Vmath::Vdiv(nqtot,m_gmat[5],1,m_jac,1,m_gmat[5],1);

                // d xi_1/d x_3
                Vmath::Vmul (nqtot,d2[1],1,d3[0],1,m_gmat[6],1);
                Vmath::Vvtvm(nqtot,d2[0],1,d3[1],1,m_gmat[6],1,m_gmat[6],1);
                Vmath::Vdiv(nqtot,m_gmat[6],1,m_jac,1,m_gmat[6],1);

                // d xi_2/d x_3
                Vmath::Vmul (nqtot,d1[0],1,d3[2],1,m_gmat[7],1);
                Vmath::Vvtvm(nqtot,d1[1],1,d3[0],1,m_gmat[7],1,m_gmat[7],1);
                Vmath::Vdiv(nqtot,m_gmat[7],1,m_jac,1,m_gmat[7],1);

                // d xim_3/d x_3
                Vmath::Vmul (nqtot,d1[1],1,d2[0],1,m_gmat[8],1);
                Vmath::Vvtvm(nqtot,d1[0],1,d2[1],1,m_gmat[8],1,m_gmat[8],1);
                Vmath::Vdiv(nqtot,m_gmat[8],1,m_jac,1,m_gmat[8],1);
            }
        }

        GeomFactors::~GeomFactors(){
        }

        bool operator==(const GeomFactors &lhs, const GeomFactors &rhs)
        {
            return (lhs.m_gmat == rhs.m_gmat &&
                    lhs.m_gtype == rhs.m_gtype &&
                    lhs.m_expdim == rhs.m_expdim &&
                    lhs.m_coordim == rhs.m_coordim);
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: GeomFactors.cpp,v $
// Revision 1.19  2008/05/30 00:33:48  delisi
// Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
//
// Revision 1.18  2008/04/06 22:32:53  bnelson
// Fixed gcc compiler warnings.
//
// Revision 1.17  2008/04/06 06:00:37  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.16  2007/12/17 20:27:19  sherwin
// Added normals to GeomFactors
//
// Revision 1.15  2007/12/06 22:47:15  pvos
// 2D Helmholtz solver updates
//
// Revision 1.14  2007/12/03 21:30:35  sherwin
// Added normal details
//
// Revision 1.13  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.12  2007/07/13 09:02:24  sherwin
// Mods for Helmholtz solver
//
// Revision 1.11  2007/07/10 22:20:55  jfrazier
// Revision of geo fac manager to test for equality.
//
// Revision 1.10  2007/07/10 17:06:30  jfrazier
// Added method and underlying structure to manage geomfactors.
//
// Revision 1.9  2007/06/07 15:54:19  pvos
// Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
// Also made corrections to various ASSERTL2 calls
//
// Revision 1.8  2007/05/31 19:13:12  pvos
// Updated NodalTriExp + LocalRegions/Project2D + some other modifications
//
// Revision 1.7  2007/05/28 21:48:41  sherwin
// Update for 2D functionality
//
// Revision 1.6  2007/05/28 08:35:25  sherwin
// Updated for localregions up to Project1D
//
// Revision 1.5  2007/05/25 17:52:01  jfrazier
// Updated to use new Array classes.
//
// Revision 1.4  2007/04/04 21:49:24  sherwin
// Update for SharedArray
//
// Revision 1.3  2007/03/29 19:23:59  bnelson
// Replaced boost::shared_array with SharedArray
//
// Revision 1.2  2007/03/25 15:48:22  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.1  2007/03/20 09:17:39  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.5  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.4  2007/03/02 12:01:58  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.3  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.2  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.1  2006/05/04 18:58:59  kirby
// *** empty log message ***
//
// Revision 1.18  2006/04/09 02:08:34  jfrazier
// Added precompiled header.
//
// Revision 1.17  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.16  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.15  2006/03/13 18:04:07  sherwin
//
// Corrected silly error in calling new
//
// Revision 1.14  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.13  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.12  2006/02/26 21:19:42  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.11  2006/02/19 01:37:32  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
