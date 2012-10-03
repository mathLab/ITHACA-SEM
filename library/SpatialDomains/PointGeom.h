////////////////////////////////////////////////////////////////////////////////
//
//  File:  PointGeom.h
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
//  Description: Point geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_POINTGEOM_H
#define NEKTAR_SPATIALDOMAINS_POINTGEOM_H


#include <StdRegions/StdRegions.hpp>
//#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/GeomFactors1D.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/MeshComponents.h>

#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdSegExp.h>

#include <LibUtilities/Foundations/Basis.h>

#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
namespace Nektar
{
    namespace SpatialDomains
    {
        class PointGeom;

        // shorthand for boost pointer
        typedef boost::shared_ptr<PointGeom> PointGeomSharedPtr;
        typedef std::vector< PointGeomSharedPtr > PointGeomVector;
        typedef std::vector< PointGeomSharedPtr >::iterator PointGeomVectorIter;
        typedef std::map<int, PointGeomSharedPtr> PointGeomMap;

        class PointGeom: public Geometry1D //???Geometry0D???
        {
            public:
                SPATIAL_DOMAINS_EXPORT PointGeom();

                SPATIAL_DOMAINS_EXPORT ~PointGeom();

                SPATIAL_DOMAINS_EXPORT static StdRegions::Orientation
                            GetPointOrientation(const SegGeom& edge1,
                                                const SegGeom& edge2);

            protected:
                static const int                kNverts = 1;
                SpatialDomains::VertexComponentSharedPtr m_verts[kNverts];

            private:

                virtual VertexComponentSharedPtr v_GetVertex(const int i) const;
                virtual int v_GetVid(int i) const;
                virtual NekDouble v_GetCoord(
                            const int i,
                            const Array<OneD,const NekDouble>& Lcoord);
                virtual void v_GetLocCoords(
                            const Array<OneD,const NekDouble>& coords,
                                  Array<OneD,NekDouble>& Lcoords);
                virtual int v_GetNumVerts() const;
        };


    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_POINTGEOM_H
