////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry1D.h
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
//  Description:  1D geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <StdRegions/StdExpansion1D.h>  // for StdExpansion1DSharedPtr, etc

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry1D;

        // shorthand for boost pointer
        typedef boost::shared_ptr<Geometry1D> Geometry1DSharedPtr;
        typedef std::vector< Geometry1DSharedPtr > Geometry1DVector;
        typedef std::vector< Geometry1DSharedPtr >::iterator Geometry1DVectorIter;

        /// 1D geometry information
        class Geometry1D: public Geometry
        {
        public:
            SPATIAL_DOMAINS_EXPORT Geometry1D();
            SPATIAL_DOMAINS_EXPORT Geometry1D(const int coordim);
            SPATIAL_DOMAINS_EXPORT virtual ~Geometry1D();

            SPATIAL_DOMAINS_EXPORT const StdRegions::StdExpansion1DSharedPtr&
                        GetXmap(const int i);
            SPATIAL_DOMAINS_EXPORT VertexComponentSharedPtr
                        GetVertex(const int i) const;
            SPATIAL_DOMAINS_EXPORT StdRegions::ExpansionType
                        DetExpansionType() const;
            SPATIAL_DOMAINS_EXPORT void WriteToFile(
                              std::ofstream& outfile,
                        const int dumpVar);
            SPATIAL_DOMAINS_EXPORT int GetEid() const;

        protected:
            using Geometry::v_GetEid;
            
            virtual int v_GetShapeDim() const;
            virtual int v_GetEid() const;
            virtual int v_GetVid(int i) const;

            virtual const StdRegions::StdExpansion1DSharedPtr&
                         v_GetXmap(const int i);
            virtual VertexComponentSharedPtr
                         v_GetVertex(const int i) const;
            virtual StdRegions::ExpansionType
                         v_DetExpansionType() const;
            virtual void v_WriteToFile(
                              std::ofstream& outfile,
                        const int dumpVar);

        };

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H

