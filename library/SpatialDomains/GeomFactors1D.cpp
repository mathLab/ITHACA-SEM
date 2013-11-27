///////////////////////////////////////////////////////////////////////////////
//
// File: GeomFactors1D.cpp
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
// Description: Implementation of 1D geometric factors.
//
///////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/Geometry2D.h>

#include <SpatialDomains/GeomFactors1D.h>
#include <LibUtilities/Foundations/Interp.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         * @class GeomFactors1D
         *
         * Computes and stores geometric factors and pointwise geometry
         * information for a 1D segment.
         */

        /**
         * The argument 'tbasis' contains the information about the quadrature
         * points on which the weighted metric terms should be specified.
         * This constructor
         * - Calculates the local factors \f$ d \xi/d x_i \f$ as a recipricol
         *   derivative \f$ d x_i/d\xi \f$;
         * - Calcuate the jacobian as \f$ J = \sqrt{\frac{d x_1}{d _\xi}.
         *   \frac{d x_1}{d _\xi} + \frac{d x_2}{d _\xi}.\frac{d x_2}{d _\xi}
         *   + \frac{d x_3}{d _\xi}\frac{d x_3}{d _\xi}} \f$
         *
         * @param   gtype       Type of geometry.
         * @param   coordim     Dimension of coordinate system.
         * @param   Coords      ?
         * @param   tbasis      Basis for derivatives
         */
        GeomFactors1D::GeomFactors1D(const GeomType gtype,
                        const int coordim,
                        const Array<OneD, const StdRegions
                                            ::StdExpansion1DSharedPtr> &Coords,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                            &tbasis) :
            GeomFactors(gtype, 1, coordim)
        {
            // Perform sanity checks
            ASSERTL1(coordim == 1 || coordim == 2 || coordim == 3,
                     "The coordinate dimension should be equal to one, two or "
                     "three for one-dimensional elements");
            ASSERTL1(tbasis.num_elements() == 1,
                     "tbasis should be an array of size one");

            for (int i = 0; i < m_coordDim; ++i)
            {
                m_coords[i] = Coords[i];
            }

            // Set the pointskey equal to the pointskey as defined
            // in 'tbasis'
            m_pointsKey[0] = tbasis[0]->GetPointsKey();
        }


        /**
         * Create a copy of an existing GeomFactors1D object.
         */
        GeomFactors1D::GeomFactors1D(const GeomFactors1D& S) :
            GeomFactors(S)
        {
        }


        /**
         *
         */
        GeomFactors1D::~GeomFactors1D()
        {
        }


        void GeomFactors1D::v_ComputeEdgeTangents(
        			const GeometrySharedPtr &geom,
        			const int edge,
        			const LibUtilities::PointsKey &to_key)
        {   	
            int k;
            int nquad= to_key.GetNumPoints();             
            Geometry2DSharedPtr g;
            ASSERTL0(g= boost::dynamic_pointer_cast<Geometry2D>(geom),
                     "FAIL");

            GeomFactorsSharedPtr gf = geom->GetMetricInfo();
            //cannot use m_type here 
            //GeomType gtype = gf->GetGtype();
            GeomType gtype= m_type;                  
            m_tangent =Array<OneD, Array<OneD, NekDouble> >(m_coordDim);
            for( k=0; k< m_coordDim; ++k)
            {
                m_tangent[k] = Array<OneD, NekDouble>(nquad);
            }
	
            int i;

            DerivStorage deriv = GetDeriv(m_pointsKey);
            //FillDeriv(deriv, m_pointsKey);

            NekDouble fac;
            // Regular geometry case
            if((gtype == eRegular)||(gtype == eMovingRegular))
            {

                for(i = 0; i < m_coordDim; ++i)
                {
                        Vmath::Fill(nquad, deriv[0][i][0],m_tangent[i],1);
                }

                // normalise
                fac = 0.0;
                for(i =0 ; i < m_coordDim; ++i)
                {
                    fac += m_tangent[i][0]*m_tangent[i][0];
                }
                fac = 1.0/sqrt(fac);
                for (i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Smul(nquad,fac,m_tangent[i],1,m_tangent[i],1);
                }   
            }
         
            else   // Set up deformed tangents 
            {

           	Array<OneD, NekDouble> jac(m_coordDim*nquad);            	
            	    
            	for(i = 0; i < m_coordDim; ++i)
                {
                	for(int j=0; j<nquad; j++)
                        {
                	m_tangent[i][j] = deriv[0][i][j];
                	}
                }
                //normalise normal vectors
                Array<OneD,NekDouble> work(nquad,0.0);
                for(i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vvtvp(nquad, m_tangent[i],1, m_tangent[i],1,work,1,work,1);
                }

                Vmath::Vsqrt(nquad,work,1,work,1);
                Vmath::Sdiv(nquad,1.0,work,1,work,1);

                for(i = 0; i < m_coordDim; ++i)
                {
                    Vmath::Vmul(nquad, m_tangent[i],1,work,1,m_tangent[i],1);                    
                }                 
            }           
        }
                
        void GeomFactors1D::v_Interp(
                    const PointsKeyArray &map_points,
                    const Array<OneD, const NekDouble> &src,
                    const PointsKeyArray &tpoints,
                    Array<OneD, NekDouble> &tgt) const
        {
            LibUtilities::Interp1D(map_points[0], src, tpoints[0], tgt);
        }

        void GeomFactors1D::v_Adjoint(
                    const Array<TwoD, const NekDouble>& src,
                    Array<TwoD, NekDouble>& tgt) const
        {
            int n = src[0].num_elements();
            //Vmath::Sdiv(n, 1.0, &src[0][0], 1, &tgt[0][0], 1);
            //Vmath::Vcopy(n, &src[0][0], 1, &tgt[0][0], 1);
            Vmath::Fill(n, 1.0, &tgt[0][0], 1);
        }

    }
}


