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
#include <SpatialDomains/Geometry2D.h>

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
         * Constructs a GeomFactors base object.
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
         * Copies an existing GeomFactors object.
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
         *
         */
        GeomFactors::~GeomFactors()
        {
        }


        /**
         * Given a set of vectors, supplied as an array of components,
         * normalise each vector individually.
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
         * Computes the vector cross-product in 3D of \a v1 and \a v2, storing
         * the result in \a v3.
         * @param   v1          First input vector.
         * @param   v2          Second input vector.
         * @param   v3          Output vector computed to be orthogonal to
         *                      both \a v1 and \a v2.
         */
        void GeomFactors::VectorCrossProd(
                        const Array<OneD, const Array<OneD, NekDouble> > &v1,
                        const Array<OneD, const Array<OneD, NekDouble> > &v2,
                              Array<OneD, Array<OneD, NekDouble> > &v3)
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



		void GeomFactors::v_ComputeEdgeTangents(
	 	 		const GeometrySharedPtr &geom2D,
	 	 		const int edge,
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
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                        &tbasis)
        {
            ASSERTL0(false, "Quadrature Metrics not implemented.");
        }


        /**
         * Placeholder function.
         */
        void GeomFactors::v_SetUpLaplacianMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                        &tbasis)
        {
            ASSERTL0(false, "Laplacian Metrics not implemented.");
        }


        /**
         * Establishes if two GeomFactors objects are equal.
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

//
// $Log: GeomFactors.cpp,v $
// Revision 1.47  2009/08/19 14:13:34  claes
// Removed Gauss-Kronrod parts
//
// Revision 1.46  2009/07/08 17:24:52  sehunchun
// Delete SetUpTanBasis and SetUp SurfaceNormal only when coordim == 3
//
// Revision 1.45  2009/07/08 11:15:51  sehunchun
// Adding Setup function fo Surface Normal and GetSurfaceNormal to obtain Surface Normal vector for a given 2D manifold
//
// Revision 1.44  2009/07/03 15:33:09  sehunchun
// Introducing m_tanbasis for tangential basis of 2D manfiold
//
// Revision 1.43  2009/06/15 01:59:21  claes
// Gauss-Kronrod updates
//
// Revision 1.42  2009/05/15 14:38:41  pvos
// Changed check for regular quads so that it also includes parallellograms
//
// Revision 1.41  2009/05/01 13:23:21  pvos
// Fixed various bugs
//
// Revision 1.40  2009/04/20 16:13:23  sherwin
// Modified Import and Write functions and redefined how Expansion is used
//
// Revision 1.39  2009/04/04 00:29:37  rcantao
// Made a few checks on J3D and companion. Remarks included.
//
// Revision 1.38  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.37  2009/01/01 02:33:29  ehan
// cleaned up the code
//
// Revision 1.36  2008/12/18 14:08:58  pvos
// NekConstants update
//
// Revision 1.35  2008/12/17 16:57:20  pvos
// Performance updates
//
// Revision 1.34  2008/12/16 14:09:07  pvos
// Performance updates
//
// Revision 1.33  2008/12/03 23:42:14  ehan
// Fixed some error for 3D Geomfactors.
//
// Revision 1.32  2008/11/25 23:02:15  ehan
// Fixed level 1 assertion violation: GeomFactors(..) only works for 1D and 2D, but not it works for 3D.
//
// Revision 1.31  2008/11/24 18:33:34  ehan
// Added 3D routines for GeomFactors()
//
// Revision 1.30  2008/09/23 22:07:32  ehan
// Modified 3D GeomFactor constructor
//
// Revision 1.29  2008/09/23 18:19:56  pvos
// Updates for working ProjectContField3D demo
//
// Revision 1.28  2008/09/20 11:32:24  ehan
// Rearranged 3D geomfactors
//
// Revision 1.27  2008/09/15 10:21:32  ehan
// Fixed some errors.
//
// Revision 1.26  2008/09/09 14:18:02  sherwin
// Removed m_normals from GeomFactor. Added GenNormals2D and additional copy type constructor
//
// Revision 1.25  2008/07/28 22:25:43  ehan
// Fixed error for deformed quads.
//
// Revision 1.24  2008/07/19 21:20:36  sherwin
// Changed normal orientation to anticlockwise
//
// Revision 1.23  2008/07/17 19:25:55  ehan
// Added 3D GeomFactors(..) .
//
// Revision 1.22  2008/07/09 11:40:25  sherwin
// Fixed the initialisation of m_expdim
//
// Revision 1.21  2008/06/13 18:07:30  ehan
// Commented out the function GeomFactors(..)
//
// Revision 1.20  2008/06/09 21:34:28  jfrazier
// Added some code for 3d.
//
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
