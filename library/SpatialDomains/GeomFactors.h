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
#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
namespace SpatialDomains
{
    // Forward declarations
    class GeomFactors;
    class Geometry;

    typedef boost::shared_ptr<Geometry> GeometrySharedPtr;

    /// Equivalence test for GeomFactors objects
    SPATIAL_DOMAINS_EXPORT bool operator==(const GeomFactors &lhs,
                                           const GeomFactors &rhs);

    /// Pointer to a GeomFactors object.
    typedef boost::shared_ptr<GeomFactors>      GeomFactorsSharedPtr;
    /// A vector of GeomFactor pointers.
    typedef std::vector< GeomFactorsSharedPtr > GeomFactorsVector;
    /// Iterator for the GeomFactorsVector.
    typedef GeomFactorsVector::iterator         GeomFactorsVectorIter;
    /// An unordered set of GeomFactor pointers.
    typedef boost::unordered_set< GeomFactorsSharedPtr >
                                                GeomFactorsSet;
    /// Iterator for the GeomFactorsSet
    typedef boost::unordered_set< GeomFactorsSharedPtr >::iterator
                                                GeomFactorsSetIter;
    /// Storage type for derivative of mapping.
    typedef Array<OneD, Array<OneD, Array<OneD,NekDouble> > >
                                                DerivStorage;

    /// Calculation and storage of geometric factors associated with the
    /// mapping from StdRegions reference elements to a given LocalRegions
    /// physical element in the mesh.
    class GeomFactors
    {
    public:
            /// Constructor for GeomFactors class.
            GeomFactors(const GeomType                              gtype,
                        const int                                   coordim,
                        const StdRegions::StdExpansionSharedPtr    &xmap,
                        const Array<OneD, Array<OneD, NekDouble> > &coords);

            /// Copy constructor.
            GeomFactors(const GeomFactors &S);

            /// Tests if two GeomFactors classes are equal.
            SPATIAL_DOMAINS_EXPORT friend bool operator==(
                const GeomFactors &lhs,
                const GeomFactors &rhs);

            /// Destructor.
            SPATIAL_DOMAINS_EXPORT ~GeomFactors();

            /// Return the derivative of the mapping with respect to the
            /// reference coordinates,
            /// \f$\frac{\partial \chi_i}{\partial \xi_j}\f$.
            inline DerivStorage GetDeriv(
                    const LibUtilities::PointsKeyVector &keyTgt);

            /// Return the Jacobian of the mapping and cache the result.
            inline const Array<OneD, const NekDouble> GetJac(
                    const LibUtilities::PointsKeyVector &keyTgt);

            /// Return the Laplacian coefficients \f$g_{ij}\f$.
            inline const Array<TwoD, const NekDouble> GetGmat(
                    const LibUtilities::PointsKeyVector &keyTgt);

            /// Return the derivative of the reference coordinates with respect
            /// to the mapping, \f$\frac{\partial \xi_i}{\partial \chi_j}\f$.
            inline const Array<TwoD, const NekDouble> GetDerivFactors(
                    const LibUtilities::PointsKeyVector &keyTgt);

            /// Returns whether the geometry is regular or deformed.
            inline GeomType GetGtype();

            /// Determine if element is valid and not self-intersecting.
            inline bool IsValid() const;

            /// Return the number of dimensions of the coordinate system.
            inline int GetCoordim() const;

            /// Computes a hash of this GeomFactors element.
            inline size_t GetHash();

    protected:
            /// Type of geometry (e.g. eRegular, eDeformed, eMovingRegular).
            GeomType m_type;
            /// Dimension of expansion.
            int m_expDim;
            /// Dimension of coordinate system.
            int m_coordDim;
            /// Validity of element (Jacobian positive)
            bool m_valid;
            /// Stores information about the expansion.
            StdRegions::StdExpansionSharedPtr m_xmap;
            /// Stores coordinates of the geometry.
            Array<OneD, Array<OneD, NekDouble> > m_coords;
            /// Jacobian vector cache
            std::map<LibUtilities::PointsKeyVector, Array<OneD, NekDouble> >
                                                m_jacCache;
            /// DerivFactors vector cache
            std::map<LibUtilities::PointsKeyVector, Array<TwoD, NekDouble> >
                                                m_derivFactorCache;
            /// Return the Xmap;
            inline StdRegions::StdExpansionSharedPtr &GetXmap();

        private:
            /// Tests if the element is valid and not self-intersecting.
            void CheckIfValid();

            SPATIAL_DOMAINS_EXPORT DerivStorage ComputeDeriv(
                    const LibUtilities::PointsKeyVector &keyTgt) const;

            /// Return the Jacobian of the mapping and cache the result.
            SPATIAL_DOMAINS_EXPORT Array<OneD, NekDouble> ComputeJac(
                    const LibUtilities::PointsKeyVector &keyTgt) const;

            /// Computes the Laplacian coefficients \f$g_{ij}\f$.
            SPATIAL_DOMAINS_EXPORT Array<TwoD, NekDouble> ComputeGmat(
                    const LibUtilities::PointsKeyVector &keyTgt) const;

            /// Return the derivative of the reference coordinates with respect
            /// to the mapping, \f$\frac{\partial \xi_i}{\partial \chi_j}\f$.
            SPATIAL_DOMAINS_EXPORT Array<TwoD, NekDouble> ComputeDerivFactors(
                    const LibUtilities::PointsKeyVector &keyTgt) const;

            /// Perform interpolation of data between two point
            /// distributions.
            void Interp(
                    const LibUtilities::PointsKeyVector &src_points,
                    const Array<OneD, const NekDouble> &src,
                    const LibUtilities::PointsKeyVector &tgt_points,
                    Array<OneD, NekDouble> &tgt) const;

            /// Compute the transpose of the cofactors matrix
            void Adjoint(
                    const Array<TwoD, const NekDouble>& src,
                    Array<TwoD, NekDouble>& tgt) const;

    };


    /// A hash functor for geometric factors. Utilises
    /// GeomFactors::GetHash.
    struct GeomFactorsHash : std::unary_function<GeomFactorsSharedPtr,
                                                 std::size_t>
    {
        std::size_t operator()(GeomFactorsSharedPtr const& p) const
        {
            return p->GetHash();
        }
    };


    /**
     * @param   tpoints     Target point distributions.
     * @returns             Derivative evaluated at target point
     *                      distributions.
     * @see                 GeomFactors::ComputeDeriv
     */
    inline DerivStorage GeomFactors::GetDeriv(
            const LibUtilities::PointsKeyVector &tpoints)
    {
        return ComputeDeriv(tpoints);
    }


    /**
     * Returns cached value if available, otherwise computes Jacobian and
     * stores result in cache.
     *
     * @param   keyTgt      Target point distributions.
     * @returns             Jacobian evaluated at target point
     *                      distributions.
     * @see                 GeomFactors::ComputeJac
     */
    inline const Array<OneD, const NekDouble> GeomFactors::GetJac(
            const LibUtilities::PointsKeyVector &keyTgt)
    {
        std::map<LibUtilities::PointsKeyVector,
            Array<OneD, NekDouble> >::const_iterator x;
        
        if ((x = m_jacCache.find(keyTgt)) != m_jacCache.end())
        {
            return x->second;
        }

        m_jacCache[keyTgt] = ComputeJac(keyTgt);

        return m_jacCache[keyTgt];

    }

    /**
     * @param   keyTgt      Target point distributions.
     * @returns             Inverse metric tensor evaluated at target point
     *                      distributions.
     * @see                 GeomFactors::ComputeGmat
     */
    inline const Array<TwoD, const NekDouble> GeomFactors::GetGmat(
            const LibUtilities::PointsKeyVector &keyTgt)
    {
        return ComputeGmat(keyTgt);
    }


    /**
     * Returns cached value if available, otherwise computes derivative
     * factors and stores result in cache.
     *
     * @param   keyTgt      Target point distributions.
     * @returns             Derivative factors evaluated at target point
     *                      distributions.
     * @see                 GeomFactors::ComputeDerivFactors
     */
    inline const Array<TwoD, const NekDouble> GeomFactors::GetDerivFactors(
            const LibUtilities::PointsKeyVector &keyTgt)
    {
        std::map<LibUtilities::PointsKeyVector,
                 Array<TwoD, NekDouble> >::const_iterator x;

        if ((x = m_derivFactorCache.find(keyTgt)) != m_derivFactorCache.end())
        {
            return x->second;
        }

        m_derivFactorCache[keyTgt] = ComputeDerivFactors(keyTgt);

        return m_derivFactorCache[keyTgt];

    }


    /**
     * A geometric shape is considered regular if it has constant geometric
     * information, and deformed if this information changes throughout the
     * shape.
     * @returns             The type of geometry.
     * @see GeomType
     */
    inline GeomType GeomFactors::GetGtype()
    {
        return m_type;
    }


    /**
     * This is greater than or equal to the expansion dimension.
     * @returns             The dimension of the coordinate system.
     */
    inline int GeomFactors::GetCoordim() const
    {
        return m_coordDim;
    }

    /**
     * The validity test is performed by testing if the Jacobian is
     * negative at any point in the shape.
     * @returns             True if the element is not self-intersecting.
     */
    inline bool GeomFactors::IsValid() const
    {
        return m_valid;
    }

    /**
     * The hash is computed from the geometry type, expansion dimension,
     * coordinate dimension and Jacobian.
     * @returns             Hash of this GeomFactors object.
     */
    inline size_t GeomFactors::GetHash()
    {
        LibUtilities::PointsKeyVector ptsKeys = m_xmap->GetPointsKeys();
        const Array<OneD, const NekDouble> jac = GetJac(ptsKeys);

        size_t hash = 0;
        boost::hash_combine(hash, (int)m_type);
        boost::hash_combine(hash, m_expDim);
        boost::hash_combine(hash, m_coordDim);
        if (m_type == eDeformed)
        {
            boost::hash_range(hash, jac.begin(), jac.end());
        }
        else
        {
            boost::hash_combine(hash, jac[0]);
        }
        return hash;
    }

    StdRegions::StdExpansionSharedPtr &GeomFactors::GetXmap(void)
    {
        return m_xmap;
    }

} //end of namespace
} //end of namespace

#endif
