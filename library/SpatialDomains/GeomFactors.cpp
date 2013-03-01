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
         */

        /**
         * @brief Constructs a GeomFactors base object.
         */
        GeomFactors::GeomFactors(const GeomType gtype,
                                 const int expdim,
                                 const int coordim,
                                 const bool UseQuadMet,
                                 const bool UseLaplMet):
            m_type(gtype),
            m_expDim(expdim),
            m_coordDim(coordim),
            m_coords(Array<OneD, StdRegions::StdExpansionSharedPtr>(m_coordDim)),
            m_isUsingQuadMetrics(UseQuadMet),
            m_isUsingLaplMetrics(UseLaplMet),
            m_pointsKey(expdim)
        {
        }

        /**
         * @brief Copies an existing GeomFactors object.
         */
        GeomFactors::GeomFactors(const GeomFactors &S) :
            m_type(S.m_type),
            m_expDim(S.m_expDim),
            m_coordDim(S.m_coordDim),
            m_coords(S.m_coords),
            m_isUsingQuadMetrics(S.m_isUsingQuadMetrics),
            m_isUsingLaplMetrics(S.m_isUsingLaplMetrics),
            m_pointsKey(S.m_pointsKey)
        {
        }


        /**
         * @brief Destructor.
         */
        GeomFactors::~GeomFactors()
        {
        }


        /**
         * @brief Given a set of vectors, supplied as an array of components,
         * normalise each vector individually.
         * 
         * @param   array       Array of vector components, \a array[i][j] with
         *                      @f$1\leq i\leq m_coordDim@f$ and
         *                      @f$0 \leq j \leq n@f$ with @f$n@f$ the number
         *                      of vectors.
         */
        void GeomFactors::VectorNormalise(
            Array<OneD, Array<OneD, NekDouble> > &array)
        {
            int ndim = array.num_elements();
            ASSERTL0(ndim > 0, "Number of components must be > 0.");
            for (int i = 1; i < ndim; ++i)
            {
                ASSERTL0(array[i].num_elements() == array[0].num_elements(),
                         "Array size mismatch in coordinates.");
            }

            int nq = array[0].num_elements();
            NekDouble Tol = 0.0000000001;
            Array<OneD, NekDouble> norm (nq, 0.0);

            // Compute the norm of each vector.
            for (int i = 0; i < ndim; ++i)
            {
                Vmath::Vvtvp(nq, array[i],  1,
                                 array[i],  1,
                                 norm,      1,
                                 norm,      1);
            }
            Vmath::Vsqrt(nq, norm, 1, norm, 1);

            // Set all norms < tol to 1.0
            for (int i = 0; i < nq; ++i)
            {
                if(abs(norm[i]) < Tol)
                {
                    norm[i] = 1.0;
                }
            }

            // Normalise the vectors by the norm
            for (int i = 0; i < ndim; ++i)
            {
                Vmath::Vdiv(nq, array[i], 1, norm, 1, array[i], 1);
            }
        }


        /**
         * @brief Computes the vector cross-product in 3D of \a v1 and \a v2,
         * storing the result in \a v3.
         * 
         * @param   v1          First input vector.
         * @param   v2          Second input vector.
         * @param   v3          Output vector computed to be orthogonal to
         *                      both \a v1 and \a v2.
         */
        void GeomFactors::VectorCrossProd(
            const Array<OneD, const Array<OneD, NekDouble> > &v1,
            const Array<OneD, const Array<OneD, NekDouble> > &v2,
                  Array<OneD,       Array<OneD, NekDouble> > &v3)
        {
            ASSERTL0(v1.num_elements() == 3,
                     "Input 1 has dimension not equal to 3.");
            ASSERTL0(v2.num_elements() == 3,
                     "Input 2 has dimension not equal to 3.");
            ASSERTL0(v3.num_elements() == 3,
                     "Output vector has dimension not equal to 3.");

            int nq = v1[0].num_elements();
            Array<OneD, NekDouble> temp(nq);

            Vmath::Vmul (nq, v1[2], 1, v2[1], 1, temp, 1);
            Vmath::Vvtvm(nq, v1[1], 1, v2[2], 1, temp, 1, v3[0], 1);

            Vmath::Vmul (nq, v1[0], 1, v2[2], 1, temp, 1);
            Vmath::Vvtvm(nq, v1[2], 1, v2[0], 1, temp, 1, v3[1], 1);

            Vmath::Vmul (nq, v1[1], 1, v2[0], 1, temp, 1);
            Vmath::Vvtvm(nq, v1[0], 1, v2[1], 1, temp, 1, v3[2], 1);
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
         * Placeholder function.
         */
        void GeomFactors::v_ComputeTangents()
        {
            ASSERTL0(false, "Cannot compute tangents for this geometry.");
        }

        /**
         * Placeholder function.
         */
        void GeomFactors::v_SetUpQuadratureMetrics(
                        LibUtilities::ShapeType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                        &tbasis)
        {
            ASSERTL0(false, "Quadrature Metrics not implemented.");
        }

        /**
         * Placeholder function.
         */
        void GeomFactors::v_SetUpLaplacianMetrics(
                        LibUtilities::ShapeType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                        &tbasis)
        {
            ASSERTL0(false, "Laplacian Metrics not implemented.");
        }

        /**
         * @brief Establishes if two GeomFactors objects are equal.
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

            if(!(lhs.m_isUsingQuadMetrics == rhs.m_isUsingQuadMetrics))
            {
                return false;
            }

            if(!(lhs.m_isUsingLaplMetrics == rhs.m_isUsingLaplMetrics))
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
    }; //end of namespace
}; //end of namespace
