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
//  Description: Geometric factors base class.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/GeomFactors.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         * @class GeomFactors
         *
         * This class stores the various geometric factors associated with a
         * specific element, necessary for fundamental integration and
         * differentiation operations as well as edge and surface normals.
         * Dimension-specific versions of this class (GeomFactors1D,
         * GeomFactors2D and GeomFactors3D) provide the majority of
         * implementation.
         *
         * Initially, these algorithms are provided with a mapping from the
         * reference region element to the physical element. Practically, this
         * is represented using a corresponding reference region element for
         * each coordinate component. Note that for straight-sided elements,
         * these elements will be of linear order. Curved elements are
         * represented using higher-order coordinate mappings. This geometric
         * order is in contrast to the order of the spectral/hp expansion order
         * on the element.
         *
         * For application of the chain rule during differentiation we require
         * the partial derivatives \f[\frac{\partial \xi_i}{\partial \chi_j}\f]
         * evaluated at the physical points of the expansion basis. We also
         * construct the inverse metric tensor \f$g^{ij}\f$ which, in the case
         * of a domain embedded in a higher-dimensional space, supports the
         * conversion of covariant quantities to contravariant quantities.
         * When the expansion dimension is equal to the coordinate dimension the
         * Jacobian of the mapping \f$\chi_j\f$ is a square matrix and
         * consequently the required terms are the entries of the inverse of the
         * Jacobian. However, in general this is not the case, so we therefore
         * implement the construction of these terms following the derivation
         * in Cantwell, et. al. \cite CaYaKiPeSh13. Given the coordinate maps
         * \f$\chi_i\f$, this comprises of five steps
         * -# Compute the terms of the Jacobian \f$\frac{\partial \chi_i}{\partial \xi_j}\f$.
         * -# Compute the metric tensor \f$g_{ij}=\mathbf{t}_{(i)}\cdot\mathbf{t}_{(j)}\f$.
         * -# Compute the square of the Jacobian determinant \f$g=|\mathbf{g}|\f$.
         * -# Compute the inverse metric tensor \f$g^{ij}\f$.
         * -# Compute the terms \f$\frac{\partial \xi_i}{\partial \chi_j}\f$.
         *
         * @see GeomFactors1D, GeomFactors2D, GeomFactors3D
         */

        /**
         * This constructor is protected since only dimension-specific
         * GeomFactors classes should be instantiated externally.
         * @param   gtype       Specified whether the geometry is regular or
         *                      deformed.
         * @param   expdim      Specifies the dimension of the expansion.
         * @param   coordim     Specifies the dimension of the coordinate
         *                      system
         */
        GeomFactors::GeomFactors(const GeomType gtype,
                                 const int expdim,
                                 const int coordim) :
            m_type(gtype),
            m_expDim(expdim),
            m_coordDim(coordim),
            m_valid(true),
            m_coords(Array<OneD, StdRegions::StdExpansionSharedPtr>(m_coordDim)),
            m_pointsKey(expdim)
        {
        }


        /**
         * @param   S           An instance of a GeomFactors class from which
         *                      to construct a new instance.
         */
        GeomFactors::GeomFactors(const GeomFactors &S) :
            m_type(S.m_type),
            m_expDim(S.m_expDim),
            m_coordDim(S.m_coordDim),
            m_valid(S.m_valid),
            m_coords(S.m_coords),
            m_pointsKey(S.m_pointsKey)
        {
        }


        /**
         *
         */
        GeomFactors::~GeomFactors()
        {
        }


        /**
         * Placeholder function.
         */
        void GeomFactors::v_ComputeEdgeTangents(
            const GeometrySharedPtr       &geom2D,
            const int                      edge,
            const LibUtilities::PointsKey &to_key)
        {
            ASSERTL0(false, "Cannot compute tangents for this geometry.");
        }

        /**
         * Placeholder function.
         */
        void GeomFactors::v_ComputeSurfaceNormals()
        {
            ASSERTL0(false, "Cannot compute surface normal for this geometry.");
        }

        /**
         * Member data equivalence is tested in the following order: shape type,
         * expansion dimension, coordinate dimension, points-keys, determinant
         * of Jacobian matrix, Laplacian coefficients.
         */
        bool operator==(const GeomFactors &lhs, const GeomFactors &rhs)
        {
            if(!(lhs.m_type == rhs.m_type))
            {
                return false;
            }

            if(!(lhs.m_expDim == rhs.m_expDim))
            {
                return false;
            }

            if(!(lhs.m_coordDim == rhs.m_coordDim))
            {
                return false;
            }

            for(int i = 0; i < lhs.m_expDim; i++)
            {
                if(!(lhs.m_pointsKey[i] == rhs.m_pointsKey[i]))
                {
                    return false;
                }
            }

            if (!(lhs.m_jac == rhs.m_jac))
            {
                return false;
            }

            if (!(lhs.m_gmat == rhs.m_gmat))
            {
                return false;
            }

            return true;
        }

        void GeomFactors::FillDeriv(
                DerivStorage &deriv,
                const Array<OneD, const LibUtilities::PointsKey>
                                         &tpoints) const
        {
            ASSERTL1(tpoints.num_elements() == m_expDim,
                     "Dimension of target basis does not match expansion basis.");

            int i = 0;
            int nqtot_map = 1;
            int nqtot_tbasis = 1;
            deriv = DerivStorage(m_expDim);
            DerivStorage d_map = DerivStorage(m_expDim);
            Array<OneD, LibUtilities::PointsKey> map_points(m_expDim);
            for (i = 0; i < m_expDim; ++i)
            {
                map_points[i] =m_coords[0]->GetBasis(i)->GetPointsKey();
                nqtot_map *= map_points[i].GetNumPoints();
                nqtot_tbasis *= tpoints[i].GetNumPoints();
                deriv[i] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
                d_map[i] = Array<OneD, Array<OneD,NekDouble> >(m_coordDim);
            }

            // Calculate local derivatives
            for(int i = 0; i < m_coordDim; ++i)
            {
                for (int j = 0; j < m_expDim; ++j)
                {
                    d_map[j][i] = Array<OneD,NekDouble>(nqtot_map);
                    deriv[j][i] = Array<OneD,NekDouble>(nqtot_tbasis);
                }

                // Transform from coefficient space to physical space
                m_coords[i]->BwdTrans(m_coords[i]->GetCoeffs(),
                                    m_coords[i]->UpdatePhys());
                // Take the derivative (calculated at the points as specified
                // in 'Coords')
                switch (m_expDim)
                {
                    case 1:
                        m_coords[i]->StdPhysDeriv(m_coords[i]->GetPhys(),
                                                d_map[0][i]);
break;
                    case 2:
                        m_coords[i]->StdPhysDeriv(m_coords[i]->GetPhys(),
                                                d_map[0][i],
                                                d_map[1][i]);
break;
                    case 3:
                        m_coords[i]->StdPhysDeriv(m_coords[i]->GetPhys(),
                                                d_map[0][i],
                                                d_map[1][i],
                                                d_map[2][i]);
break;
                }
            }

            v_Interp(map_points, d_map, tpoints, deriv);
        }

    }; //end of namespace
}; //end of namespace
