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
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY0D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY0D_H

#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

namespace Nektar
{
namespace SpatialDomains
{
class Geometry0D;

// shorthand for boost pointer
typedef std::shared_ptr<Geometry0D> Geometry0DSharedPtr;
typedef std::vector<Geometry0DSharedPtr> Geometry0DVector;
typedef std::vector<Geometry0DSharedPtr>::iterator Geometry0DVectorIter;

/// 1D geometry information
class Geometry0D : public Geometry
{
public:
    SPATIAL_DOMAINS_EXPORT Geometry0D();
    SPATIAL_DOMAINS_EXPORT Geometry0D(const int coordim);
    SPATIAL_DOMAINS_EXPORT virtual ~Geometry0D();

    SPATIAL_DOMAINS_EXPORT static const int kDim = 0;

protected:
    virtual int v_GetShapeDim() const;
};

} // end of namespace
} // end of namespace

#endif // NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H
