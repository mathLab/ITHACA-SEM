////////////////////////////////////////////////////////////////////////////////
//
//  File:  GeomFactors.h
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
//  Description: Geometric Factors base class
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS_H

#include <LibUtilities/Foundations/Points.h>

#include <SpatialDomains/SpatialDomains.hpp>
#include <StdRegions/StdExpansion1D.h>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdExpansion3D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        // Forward declarations and useful typedefs
        class GeomFactors;
        class GeomFactors2D;
        typedef boost::shared_ptr<Geometry> GeometrySharedPtr;

        bool operator==(const GeomFactors &lhs, const GeomFactors &rhs);

        /// Pointer to a GeomFactors object.
        typedef boost::shared_ptr<GeomFactors>      GeomFactorsSharedPtr;
        /// A vector of GeomFactor pointers.
        typedef std::vector< GeomFactorsSharedPtr > GeomFactorsVector;
        /// Iterator for the GeomFactorsVector.
        typedef GeomFactorsVector::iterator GeomFactorsVectorIter;

        /// Describes the principle direction for tangents on the domain.
        enum GeomTangents
        {
            eTangentX,          ///< X coordinate direction.
            eTangentY,          ///< Y coordinate direction.
            eTangentZ,          ///< Z coordinate direction.
            eTangentCircular,   ///< Circular around the centre of domain.
            eLOCAL,             ///< No Principal direction.
            SIZE_GeomTangents
        };

        /// Session file names associated with tangent principle directions.
        const char* const GeomTangentsMap[] =
        {
            "TangentX",
            "TangentY",
            "TangentZ",
            "TangentCircular",
            "LOCAL",
        };

        /// Calculation and storage of geometric factors.
        class GeomFactors
        {
        public:

            friend bool operator==( const GeomFactors &lhs,
                                    const GeomFactors &rhs);
            friend bool operator<(  const GeomFactors &lhs,
                                    const GeomFactors &rhs);

            /// Destructor.
            virtual ~GeomFactors();

            /// Return the type of geometry.
            inline GeomType GetGtype();

            /// Return the Jacobian
            inline const Array<OneD, const NekDouble> &GetJac() const;

            /// Return the G matrix.
            inline const Array<TwoD, const NekDouble> &GetGmat() const;

            /// Return the number of dimensions of the coordinate system.
            inline int GetCoordim() const;

            /// Flag indicating if quadrature metrics are in use.
            inline bool IsUsingQuadMetrics() const;

            /// Flag indicating if Laplacian metrics are in use.
            inline bool IsUsingLaplMetrics() const;

            /// Flag indicating if Tangents have been computed.
            inline bool IsUsingTangents() const;

            /// Set up quadrature metrics
            inline void SetUpQuadratureMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                       &tbasis);

            /// Set up Laplacian metrics
            inline void SetUpLaplacianMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                       &tbasis);

            /// Set up Tangents
            inline void SetUpTangents();

            /// Retrieve the quadrature metrics.
            inline const Array<OneD, const NekDouble> &GetQuadratureMetrics()
                        const;

            /// Retrieve the Laplacian metrics.
            inline const Array<TwoD, const NekDouble> &GetLaplacianMetrics()
                        const;

            /// Indicates if the Laplacian metric with index \a indx is zero.
            inline bool LaplacianMetricIsZero(const int indx) const;

            /// Computes the edge normals from a 2D element
            inline void ComputeNormals(
                            const GeometrySharedPtr &geom2D,
                            const int edge,
                            const LibUtilities::PointsKey &to_key);

            /// Computes the edge normals for 1D geometries only.
            inline void ComputeEdgeNormals(
                            const int edge,
                            const LibUtilities::PointsKey &to_key,
                            Array<OneD, Array<OneD, NekDouble> > &output) const;

            /// Set tangent orientation
            inline void SetTangentOrientation(std::string conn);

            /// Set tangent circular orientation centre.
            inline void SetTangentCircularCentre(
                            Array<OneD,NekDouble> &centre);

            /// Returns the normal vectors evaluated at each quadrature point.
            inline const Array<OneD, const Array<OneD, NekDouble> >
                                                            &GetNormal() const;

            /// Returns a single tangent vector.
            inline const Array<OneD, const Array<OneD, NekDouble> >
                                                            &GetTangent(int i);

            /// Set the G-matrix data.
            inline void ResetGmat(const Array<OneD, const NekDouble> &ndata,
                                  const int nq, const int expdim,
                                  const int coordim);

            /// Set the Jacobian data.
            inline void ResetJac(int nq,
                        const Array<OneD, const NekDouble> &ndata);

            /// Normalises a set of vectors.
            static void VectorNormalise(
                        Array<OneD, Array<OneD, NekDouble> > &array);

            /// Computes the cross-product between sets of vectors.
            static void VectorCrossProd(
                        const Array<OneD, const Array<OneD, NekDouble> > &in1,
                        const Array<OneD, const Array<OneD, NekDouble> > &in2,
                              Array<OneD, Array<OneD, NekDouble> > &out);

        protected:
            /// Type of geometry (e.g. eRegular, eDeformed, eMovingRegular).
            GeomType m_type;
            /// Dimension of expansion.
            int m_expDim;
            /// Dimension of coordinate system.
            int m_coordDim;
            /// Stores information about the expansion.
            Array<OneD, StdRegions::StdExpansionSharedPtr> m_coords;
            /// Use Quadrature metrics
            bool m_isUsingQuadMetrics;
            /// Use Laplacian metrics
            bool m_isUsingLaplMetrics;
            /// Principle tangent direction.
            enum GeomTangents m_tangentDir;
            /// Principle tangent circular dir coords
            Array<OneD,NekDouble> m_tangentDirCentre;

            /// Jacobian. If geometry is regular, or moving regular, this is
            /// just an array of one element - the value of the Jacobian across
            /// the whole element. If deformed, the array has the size of the
            /// number of quadrature points and contains the Jacobian evaluated
            /// at each of those points.
            Array<OneD,NekDouble> m_jac;

            ///
            Array<OneD,NekDouble> m_weightedjac;

            /// Array of size coordim x nquad which holds the inverse of the
            /// derivative of the local map in each direction at each
            /// quadrature point.
            Array<TwoD,NekDouble> m_gmat;

            ///
            Array<TwoD,NekDouble> m_laplacianmetrics;

            ///
            Array<OneD,bool>      m_laplacianMetricIsZero;

            /// Array of size coordim which stores a key describing the
            /// location of the quadrature points in each dimension.
            Array<OneD,LibUtilities::PointsKey> m_pointsKey;

            /// Array of derivatives of size (m_expDim)x(mCoordim)x(nq)
            Array<OneD,Array<OneD,Array<OneD,NekDouble> > > m_deriv;

            /// Array of size (m_coordDim-1)x(m_coordDim x nq).
            Array<OneD, Array<OneD, Array<OneD,NekDouble> > > m_tangents;

            /// Array of size (coordim)x(nquad) which holds the components of
            /// the normal vector at each quadrature point. The array is
            /// populated as \a m_coordDim consecutive blocks of size equal to
            /// the number of quadrature points. Each block holds a component
            /// of the normal vectors.
            Array<OneD, Array<OneD,NekDouble> > m_normal;

            /// Instance for a specific expansion/coordinate dimension without
            /// the generation of any factors. This constructor is protected
            /// since only dimension-specific GeomFactors classes should be
            /// instantiated externally.
            GeomFactors(const GeomType gtype,
                        const int expdim,
                        const int coordim,
                        const bool UseQuadMet,
                        const bool UseLaplMet);

            /// Copy constructor.
            GeomFactors(const GeomFactors &S);

        private:
            /// (1D only) Compute normals based on a 2D element.
            virtual void v_ComputeNormals(
                        const GeometrySharedPtr &geom2D,
                        const int edge,
                        const LibUtilities::PointsKey &to_key);

            /// (2D only) Compute the outward normals for a given edge.
            virtual void v_ComputeEdgeNormals(
                        const int edge,
                        const LibUtilities::PointsKey &to_key,
                        Array<OneD, Array<OneD, NekDouble> > &output) const;

            /// Set up surface normals
            virtual void v_ComputeSurfaceNormals();

            /// Set up the tangent vectors
            virtual void v_ComputeTangents();

            /// Set up quadrature metrics
            virtual void v_SetUpQuadratureMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                       &tbasis);

            /// Set up Laplacian metrics
            virtual void v_SetUpLaplacianMetrics(
                        StdRegions::ExpansionType shape,
                        const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                       &tbasis);

        };

        /// Return the type of geometry.
        inline GeomType GeomFactors::GetGtype()
        {
            return m_type;
        }

        /// Return the Jacobian
        inline const Array<OneD, const NekDouble> &GeomFactors::GetJac() const
        {
            return m_jac;
        }

        /// Return the G matrix.
        inline const Array<TwoD, const NekDouble> &GeomFactors::GetGmat() const
        {
            return m_gmat;
        }

        /// Return the number of dimensions of the coordinate system.
        inline int GeomFactors::GetCoordim() const
        {
            return m_coordDim;
        }

        /// Flag indicating if quadrature metrics are in use.
        inline bool GeomFactors::IsUsingQuadMetrics() const
        {
            return m_isUsingQuadMetrics;
        }

        /// Flag indicating if Laplacian metrics are in use.
        inline bool GeomFactors::IsUsingLaplMetrics() const
        {
            return m_isUsingLaplMetrics;
        }

        /// Flag indicating if Tangents are in use.
        inline bool GeomFactors::IsUsingTangents() const
        {
            return (m_tangents.num_elements() != 0);
        }

        /// Set up quadrature metrics
        inline void GeomFactors::SetUpQuadratureMetrics(
                    StdRegions::ExpansionType shape,
                    const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                    &tbasis)
        {
            v_SetUpQuadratureMetrics(shape, tbasis);
        }

        /// Set up Laplacian metrics
        inline void GeomFactors::SetUpLaplacianMetrics(
                    StdRegions::ExpansionType shape,
                    const Array<OneD, const LibUtilities::BasisSharedPtr>
                                                                    &tbasis)
        {
            v_SetUpLaplacianMetrics(shape, tbasis);
        }

        /// Set up Tangents
        inline void GeomFactors::SetUpTangents()
        {
            cout << "GeomFactors: Setting up tangents" << endl;
            v_ComputeTangents();
        }

        /// Retrieve the quadrature metrics.
        inline const Array<OneD, const NekDouble>
                                    &GeomFactors::GetQuadratureMetrics() const
        {
            ASSERTL0(m_isUsingQuadMetrics,
                     "This metric has not been set up for this type of "
                     "expansion");
            return m_weightedjac;
        }

        /// Retrieve the Laplacian metrics.
        inline const Array<TwoD, const NekDouble>
                                    &GeomFactors::GetLaplacianMetrics() const
        {
            ASSERTL0(m_isUsingLaplMetrics,
                     "This metric has not been set up for this type of "
                     "expansion");
            return m_laplacianmetrics;
        }

        /// Indicates if the Laplacian metric with index \a indx is zero.
        inline bool GeomFactors::LaplacianMetricIsZero(const int indx) const
        {
            ASSERTL0(m_isUsingLaplMetrics,
                     "This metric has not been set up for this type of "
                     "expansion");
            return m_laplacianMetricIsZero[indx];
        }

        /// Computes the edge normals from a 2D element
        inline void GeomFactors::ComputeNormals(
                        const GeometrySharedPtr &geom2D,
                        const int edge,
                        const LibUtilities::PointsKey &to_key)
        {
            v_ComputeNormals(geom2D, edge, to_key);
        }

        /// Computes the edge normals for 1D geometries only.
        inline void GeomFactors::ComputeEdgeNormals(
                        const int edge,
                        const LibUtilities::PointsKey &to_key,
                        Array<OneD, Array<OneD, NekDouble> > &output) const
        {
            v_ComputeEdgeNormals(edge, to_key, output);
        }

        /// Set tangent orientation
        inline void GeomFactors::SetTangentOrientation(std::string conn)
        {
            if (conn == "TangentX")         m_tangentDir = eTangentX;
            if (conn == "TangentY")         m_tangentDir = eTangentY;
            if (conn == "TangentZ")         m_tangentDir = eTangentZ;
            if (conn == "TangentCircular")  m_tangentDir = eTangentCircular;
            if (conn == "LOCAL")            m_tangentDir = eLOCAL;
        }

        /**
         * Sets the centre point for circular tangent vectors.
         * @param   centre      Array holding coordinates of centre point.
         */
        inline void GeomFactors::SetTangentCircularCentre(
                    Array<OneD,NekDouble> &centre)
        {
            m_tangentDirCentre = centre;
        }

        /// Returns the normal vectors evaluated at each quadrature point.
        inline const Array<OneD, const Array<OneD, NekDouble> >
                                                &GeomFactors::GetNormal() const
        {
            return m_normal;
        }

        /// Returns a single tangent vector.
        inline const Array<OneD, const Array<OneD, NekDouble> >
                                                &GeomFactors::GetTangent(int i)
        {
            ASSERTL0(i < m_expDim,
                     "Index must be less than expansion dimension.");
            if (m_tangents.num_elements() == 0) {
                v_ComputeTangents();
            }
            return m_tangents[i];
        }

        /// Set the G-matrix data.
        inline void GeomFactors::ResetGmat(
                        const Array<OneD, const NekDouble> &ndata,
                        const int nq, const int expdim,
                        const int coordim)
        {
            m_gmat = Array<TwoD,NekDouble>(m_expDim * m_coordDim, nq,
                                           ndata.data());
        }

        /// Set the Jacobian data.
        inline void GeomFactors::ResetJac(int nq,
                        const Array<OneD, const NekDouble> &ndata)
        {
            m_jac = Array<OneD, NekDouble>(nq, ndata.data());
        }
    } //end of namespace
} //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GeomFactors_H

//
// $Log: GeomFactors.h,v $
// Revision 1.32  2010/01/07 16:00:18  sehunchun
// Generalizing TangentCircular ...
//
// Revision 1.31  2010/01/06 14:53:03  cantwell
// Added geominfo parameters TangentCentre{X,Y} for circular tangent vectors.
//
// Revision 1.30  2009/12/15 18:09:02  cantwell
// Split GeomFactors into 1D, 2D and 3D
// Added generation of tangential basis into GeomFactors
// Updated ADR2DManifold solver to use GeomFactors for tangents
// Added <GEOMINFO> XML session section support in MeshGraph
// Fixed const-correctness in VmathArray
// Cleaned up LocalRegions code to generate GeomFactors
// Removed GenSegExp
// Temporary fix to SubStructuredGraph
// Documentation for GlobalLinSys and GlobalMatrix classes
//
// Revision 1.29  2009/07/08 17:24:52  sehunchun
// Delete SetUpTanBasis and SetUp SurfaceNormal only when coordim == 3
//
// Revision 1.28  2009/07/08 11:15:51  sehunchun
// Adding Setup function fo Surface Normal and GetSurfaceNormal to obtain Surface Normal vector for a given 2D manifold
//
// Revision 1.27  2009/07/03 15:33:09  sehunchun
// Introducing m_tanbasis for tangential basis of 2D manfiold
//
// Revision 1.26  2009/07/02 13:32:24  sehunchun
// *** empty log message ***
//
// Revision 1.25  2009/05/15 14:38:41  pvos
// Changed check for regular quads so that it also includes parallellograms
//
// Revision 1.24  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.23  2008/12/17 16:57:20  pvos
// Performance updates
//
// Revision 1.22  2008/12/16 14:09:07  pvos
// Performance updates
//
// Revision 1.21  2008/11/24 18:33:10  ehan
// Added 3D routines for GeomFactors()
//
// Revision 1.20  2008/09/09 22:46:51  ehan
// Fixed error; extra qualification ‘Nektar::SpatialDomains::GeomFactors::’ on member ‘GenNormals2D’
//
// Revision 1.19  2008/09/09 14:18:02  sherwin
// Removed m_normals from GeomFactor. Added GenNormals2D and additional copy type constructor
//
// Revision 1.18  2008/07/17 19:27:22  ehan
// Added 3D GeomFactors(..).
//
// Revision 1.17  2008/06/09 21:34:28  jfrazier
// Added some code for 3d.
//
// Revision 1.16  2008/04/06 06:00:37  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.15  2007/12/17 20:27:23  sherwin
// Added normals to GeomFactors
//
// Revision 1.14  2007/12/03 21:30:43  sherwin
// Added normal details
//
// Revision 1.13  2007/07/22 23:04:23  bnelson
// Backed out Nektar::ptr.
//
// Revision 1.12  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.11  2007/07/10 22:20:59  jfrazier
// Revision of geo fac manager to test for equality.
//
// Revision 1.10  2007/07/10 17:06:31  jfrazier
// Added method and underlying structure to manage geomfactors.
//
// Revision 1.9  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.8  2007/05/28 08:35:26  sherwin
// Updated for localregions up to Project1D
//
// Revision 1.7  2007/05/25 17:52:02  jfrazier
// Updated to use new Array classes.
//
// Revision 1.6  2007/05/17 18:45:25  jfrazier
// Minor changes to accommodate Array class.
//
// Revision 1.5  2007/04/08 03:34:48  jfrazier
// Updated to compile with SharedArray.  This has not been converted to SharedArray, just made to work with others that have been converted.
//
// Revision 1.4  2007/04/04 21:49:24  sherwin
// Update for SharedArray
//
// Revision 1.3  2007/03/29 19:24:00  bnelson
// Replaced boost::shared_array with SharedArray
//
// Revision 1.2  2007/03/25 15:48:22  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.1  2007/03/20 09:17:39  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.7  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.6  2007/03/02 12:01:59  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.5  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.4  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.3  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.2  2006/05/29 17:05:17  sherwin
// Updated to use shared_ptr around Geom types - added typedef
//
// Revision 1.1  2006/05/04 18:58:59  kirby
// *** empty log message ***
//
// Revision 1.16  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.15  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.14  2006/03/13 18:20:03  sherwin
//
// Fixed error in ResetGmat:
//
// Revision 1.13  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.12  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.11  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.10  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

