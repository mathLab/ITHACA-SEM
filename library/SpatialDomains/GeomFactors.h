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
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS_H

#include <boost/unordered_set.hpp>

#include <LibUtilities/Foundations/Basis.h>
#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <StdRegions/StdExpansion.h>    // for StdExpansionSharedPtr
#include <StdRegions/StdRegions.hpp>


namespace Nektar
{
    namespace LibUtilities
    {
        class PointsKey;
    }

    namespace SpatialDomains
    {
        // Forward declarations and useful typedefs
        class GeomFactors;
        class Geometry;
        typedef boost::shared_ptr<Geometry> GeometrySharedPtr;

        SPATIAL_DOMAINS_EXPORT bool operator==(const GeomFactors &lhs, const GeomFactors &rhs);

        /// Pointer to a GeomFactors object.
        typedef boost::shared_ptr<GeomFactors>      GeomFactorsSharedPtr;
        /// A vector of GeomFactor pointers.
        typedef std::vector< GeomFactorsSharedPtr > GeomFactorsVector;
        /// Iterator for the GeomFactorsVector.
        typedef GeomFactorsVector::iterator GeomFactorsVectorIter;
        /// An unordered set of GeomFactor pointers.
        typedef boost::unordered_set< GeomFactorsSharedPtr > GeomFactorsSet;
        /// Iterator for the GeomFactorsSet
        typedef boost::unordered_set< GeomFactorsSharedPtr >::iterator GeomFactorsSetIter;
        
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
            SPATIAL_DOMAINS_EXPORT friend bool operator==( 
                const GeomFactors &lhs,
                const GeomFactors &rhs);
            SPATIAL_DOMAINS_EXPORT friend bool operator<(
                const GeomFactors &lhs,
                const GeomFactors &rhs);

            SPATIAL_DOMAINS_EXPORT virtual ~GeomFactors();

            /// Return the type of geometry.
            inline GeomType GetGtype();

            /// Return the Jacobian
            inline const Array<OneD, const NekDouble> &GetJac() const;

            /// Return the G matrix.
            inline const Array<TwoD, const NekDouble> &GetGmat() const;

            /// Return the G matrix.
            inline const Array<
                OneD, const Array<OneD, Array<OneD, NekDouble> > > 
                    &GetDeriv() const;

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

            /// Computes the edge tangents from 1D element
            inline void ComputeEdgeTangents(
                const GeometrySharedPtr       &geom2D,
                const int                      edge,
                const LibUtilities::PointsKey &to_key);
            
            /// Set tangent orientation
            inline void SetTangentOrientation(std::string conn);

            /// Set tangent circular orientation centre.
            inline void SetTangentCircularCentre(
                            Array<OneD,NekDouble> &centre);
            
            /// Returns the tangent vectors evaluated at each quadrature point
            /// for 1D elements.  The tangent vectors are set using the
            /// function ComputeEdgeTangents.
            inline const Array<OneD, const Array<OneD, NekDouble> >
                &GetEdgeTangent() const;

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
            SPATIAL_DOMAINS_EXPORT static void VectorNormalise(
                        Array<OneD, Array<OneD, NekDouble> > &array);

            /// Computes the cross-product between sets of vectors.
            SPATIAL_DOMAINS_EXPORT static void VectorCrossProd(
                        const Array<OneD, const Array<OneD, NekDouble> > &in1,
                        const Array<OneD, const Array<OneD, NekDouble> > &in2,
                              Array<OneD, Array<OneD, NekDouble> > &out);

            /// Returns the LibUtilities::PointsKey object associated with a
            /// co-ordinate direction.
            SPATIAL_DOMAINS_EXPORT const LibUtilities::PointsKey 
                &GetPointsKey(unsigned int i) const
            {
                ASSERTL1(i < m_pointsKey.num_elements(), "PointsKey out of range.");
                return m_pointsKey[i];
            }

            /// Computes a hash of this GeomFactors element, based on the
            /// type, expansion/co-ordinate dimensions, metric, Jacobian and
            /// the geometric factors themselves.
            size_t GetHash() const
            {
                size_t hash = 0;
                boost::hash_combine(hash, (int)m_type);
                boost::hash_combine(hash, m_expDim);
                boost::hash_combine(hash, m_coordDim);
                boost::hash_combine(hash, m_isUsingQuadMetrics);
                boost::hash_combine(hash, m_isUsingLaplMetrics);
                boost::hash_range(hash, m_jac.begin(), m_jac.end());
                for (int i = 0; i < m_gmat.GetRows(); ++i)
                {
                    boost::hash_range(hash, m_gmat[i].begin(), m_gmat[i].end());
                }
                return hash;
            }

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

            /// Array of size (coordim)x(nquad) which holds the components of
            /// the tangent vector at each quadrature point. The array is
            /// populated as \a m_coordDim consecutive blocks of size equal to
            /// the number of quadrature points. Each block holds a component
            /// of the tangent vectors.
            Array<OneD, Array<OneD,NekDouble> > m_tangent;
            
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
            /// (2D only) Compute tangents based on a 1D element.
            virtual void v_ComputeEdgeTangents(
            	    	const GeometrySharedPtr &geom2D,
            	    	const int edge,
            	    	const LibUtilities::PointsKey &to_key);
            
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

        /// A hash functor for geometric factors. Utilises
        /// GeomFactors::GetHash.
        struct GeomFactorsHash : std::unary_function<GeomFactorsSharedPtr, std::size_t>
        {
            std::size_t operator()(GeomFactorsSharedPtr const& p) const
            {
                return p->GetHash();
            }
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

        /// Return the G matrix.
        inline const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >
            &GeomFactors::GetDeriv() const
        {
            return m_deriv;
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


        /// Computes the edge tangents from a 1D element
        inline void GeomFactors::ComputeEdgeTangents(
            const GeometrySharedPtr &geom2D,
            const int edge,
            const LibUtilities::PointsKey &to_key)
        {
            v_ComputeEdgeTangents(geom2D, edge, to_key);
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

        /// Returns the tangent vectors evaluated at each quadrature point for
        /// 1D elements.  The tangent vectors are set using the function
        /// ComputeEdgeTangents.
        inline const Array<OneD, const Array<OneD, NekDouble> >
            &GeomFactors::GetEdgeTangent() const
        {
            ASSERTL0(m_tangent.num_elements()>0," tangent vectors are not computed for this line");	
            return m_tangent;
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

#endif
