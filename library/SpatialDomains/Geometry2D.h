////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry2D.h
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
//  Description: 2D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion2D.h>  // for StdExpansion2DSharedPtr, etc

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry1D;
        class Geometry2D;
        class VertexComponent;
        class SegGeom;

        // shorthand for boost pointer
        typedef boost::shared_ptr<VertexComponent> VertexComponentSharedPtr;
        typedef boost::shared_ptr<Geometry1D> Geometry1DSharedPtr;
        typedef boost::shared_ptr<Geometry2D> Geometry2DSharedPtr;
        typedef boost::shared_ptr<SegGeom>    SegGeomSharedPtr;
        typedef std::vector< Geometry2DSharedPtr > Geometry2DVector;
        typedef std::vector< Geometry2DSharedPtr >::iterator Geometry2DVectorIter;

        /// 2D geometry information
        class Geometry2D: public Geometry
        {
        public:
            SPATIAL_DOMAINS_EXPORT Geometry2D();
            SPATIAL_DOMAINS_EXPORT Geometry2D(const int coordim);
            SPATIAL_DOMAINS_EXPORT virtual ~Geometry2D();

            //---------------------------------------
            // Helper functions
            //---------------------------------------

            SPATIAL_DOMAINS_EXPORT int GetFid() const;
            SPATIAL_DOMAINS_EXPORT const VertexComponentSharedPtr
                        GetVertex(int i) const;
            SPATIAL_DOMAINS_EXPORT const Geometry1DSharedPtr
                        GetEdge(int i) const;
            SPATIAL_DOMAINS_EXPORT const Geometry2DSharedPtr
                        GetFace(int i) const;

            SPATIAL_DOMAINS_EXPORT int WhichEdge(SegGeomSharedPtr edge);
            SPATIAL_DOMAINS_EXPORT int WhichFace(Geometry2DSharedPtr face);

            SPATIAL_DOMAINS_EXPORT StdRegions::StdExpansion2DSharedPtr
                        GetXmap(const int i);
            SPATIAL_DOMAINS_EXPORT StdRegions::StdExpansion2DSharedPtr
                        operator[](const int i) const;

            SPATIAL_DOMAINS_EXPORT const LibUtilities::BasisSharedPtr
                    GetEdgeBasis(const int i, const int j);

            //---------------------------------------
            // Orientation functions
            //---------------------------------------

            SPATIAL_DOMAINS_EXPORT StdRegions::Orientation
                        GetFaceOrient(const int i) const;
            SPATIAL_DOMAINS_EXPORT StdRegions::Orientation
                        GetCartesianEorient(const int i) const;


        protected:

            Array<OneD, StdRegions::StdExpansion2DSharedPtr> m_xmap;

            void NewtonIterationForLocCoord(const Array<OneD, const NekDouble> &coords, 
                                       Array<OneD,NekDouble> &Lcoords);

        private:
            //---------------------------------------
            // Helper functions
            //---------------------------------------
            
            using Geometry::v_GetFid;

            virtual int                         v_GetShapeDim() const;
            virtual int                         v_GetFid() const;
            virtual int                         v_GetEid(int i) const;
            virtual const VertexComponentSharedPtr v_GetVertex(int i) const;
            virtual const Geometry1DSharedPtr   v_GetEdge(int i) const;
            virtual const Geometry2DSharedPtr   v_GetFace(int i) const;
            virtual StdRegions::Orientation v_GetFaceOrient(const int i) const;
            virtual StdRegions::Orientation v_GetEorient(const int i) const;
            virtual StdRegions::Orientation v_GetCartesianEorient(const int i) const;
            virtual int                         v_WhichEdge(SegGeomSharedPtr edge);
            virtual int                         v_WhichFace(Geometry2DSharedPtr face);

            virtual const LibUtilities::BasisSharedPtr
                            v_GetEdgeBasis(const int i, const int j);
            virtual bool    v_ContainsPoint(
                        const Array<OneD, const NekDouble>& gloCoord,
                        NekDouble tol = 0.0);

        };


    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY2D_H
