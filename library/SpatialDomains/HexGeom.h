////////////////////////////////////////////////////////////////////////////////
//
//  File: HexGeom.h
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
//  Description: Hexahedral geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_HEXGEOM_H
#define NEKTAR_SPATIALDOMAINS_HEXGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdHexExp.h>

#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class HexGeom : public Geometry3D
        {
        public:
            SPATIAL_DOMAINS_EXPORT HexGeom();
            SPATIAL_DOMAINS_EXPORT HexGeom(
                const QuadGeomSharedPtr faces[]);
            /*
            SPATIAL_DOMAINS_EXPORT HexGeom(
                const QuadGeomSharedPtr faces[], 
                const Array<OneD, StdRegions::StdExpansion3DSharedPtr> &xMap);
            */
            SPATIAL_DOMAINS_EXPORT ~HexGeom();

            SPATIAL_DOMAINS_EXPORT static const int         kNverts  = 8;
            SPATIAL_DOMAINS_EXPORT static const int         kNedges  = 12;
            SPATIAL_DOMAINS_EXPORT static const int         kNqfaces = 6;
            SPATIAL_DOMAINS_EXPORT static const int         kNtfaces = 0;
            SPATIAL_DOMAINS_EXPORT static const int         kNfaces  = kNqfaces + kNtfaces;
            SPATIAL_DOMAINS_EXPORT static const std::string XMLElementType;
            
        protected:
            virtual void v_GenGeomFactors(
                const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis);
            virtual void v_GetLocCoords(
                const Array<OneD, const NekDouble> &coords,
                      Array<OneD,       NekDouble> &Lcoords);
            virtual bool v_ContainsPoint(
                const Array<OneD, const NekDouble> &gloCoord, 
                      NekDouble                     tol = 0.0);
            virtual int v_GetNumVerts() const;
            virtual int v_GetNumEdges() const;
            virtual int v_GetNumFaces() const;
	    virtual int v_GetVertexEdgeMap(
                const int i, const int j) const;
	    virtual int v_GetVertexFaceMap(
                const int i, const int j) const;
	    virtual int v_GetEdgeFaceMap(
                const int i, const int j) const;

        private:
            void SetUpLocalEdges();
            void SetUpLocalVertices();
            void SetUpEdgeOrientation();
            void SetUpFaceOrientation();
        };

        typedef boost::shared_ptr<HexGeom> HexGeomSharedPtr;
        typedef std::vector<HexGeomSharedPtr> HexGeomVector;
        typedef std::vector<HexGeomSharedPtr>::iterator HexGeomVectorIter;
        typedef std::map<int, HexGeomSharedPtr> HexGeomMap;
        typedef std::map<int, HexGeomSharedPtr>::iterator HexGeomMapIter;
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_HEXGEOM_H
