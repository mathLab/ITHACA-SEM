////////////////////////////////////////////////////////////////////////////////
//
//  File: PrismGeom.h
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
//  Description: Prismatic geometry definition.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_PRISMGEOM_H
#define NEKTAR_SPATIALDOMAINS_PRISMGEOM_H

#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/GeomFactors3D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class PrismGeom: public LibUtilities::GraphVertexObject, 
                         public Geometry3D
        {
        public:
            SPATIAL_DOMAINS_EXPORT PrismGeom ();
            SPATIAL_DOMAINS_EXPORT PrismGeom(
                const Geometry2DSharedPtr faces[]);
            SPATIAL_DOMAINS_EXPORT ~PrismGeom();
            
            SPATIAL_DOMAINS_EXPORT static const int kNverts  = 6;
            SPATIAL_DOMAINS_EXPORT static const int kNedges  = 9;
            SPATIAL_DOMAINS_EXPORT static const int kNqfaces = 3;
            SPATIAL_DOMAINS_EXPORT static const int kNtfaces = 2;
            SPATIAL_DOMAINS_EXPORT static const int kNfaces  = kNqfaces + kNtfaces;
            SPATIAL_DOMAINS_EXPORT static const std::string XMLElementType;

        protected:
            virtual void v_GetLocCoords(
                const Array<OneD, const NekDouble> &coords,
                      Array<OneD,       NekDouble> &Lcoords);
            virtual int v_GetNumVerts() const;
            virtual int v_GetNumEdges() const;
            
        private:
            void SetUpLocalEdges();
            void SetUpLocalVertices();
            void SetUpEdgeOrientation();
            void SetUpFaceOrientation();
        };

        typedef boost::shared_ptr<PrismGeom> PrismGeomSharedPtr;
        typedef std::vector<PrismGeomSharedPtr> PrismGeomVector;
        typedef std::vector<PrismGeomSharedPtr>::iterator PrismGeomVectorIter;
        typedef std::map<int, PrismGeomSharedPtr> PrismGeomMap;
        typedef std::map<int, PrismGeomSharedPtr>::iterator PrismGeomMapIter;
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_PRISMGEOM_H
