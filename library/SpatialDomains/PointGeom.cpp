////////////////////////////////////////////////////////////////////////////////
//
//  File:  PointGeom.cpp
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
//  Description: Point geometry information
//
////////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SegGeom.h>

namespace Nektar
{
namespace SpatialDomains
{
PointGeom::PointGeom() : NekPoint<NekDouble>(0.0, 0.0, 0.0)
{
    m_shapeType = LibUtilities::ePoint;
    m_coordim = 0;
    m_globalID = 0;
}

PointGeom::PointGeom(
    const int coordim, const int vid, NekDouble x, NekDouble y, NekDouble z)
    : NekPoint<NekDouble>(x, y, z)
{
    m_shapeType = LibUtilities::ePoint;
    m_coordim = coordim;
    m_globalID = vid;
}

// copy constructor
PointGeom::PointGeom(const PointGeom &T)
    : Geometry0D(T),
      NekPoint<NekDouble>(T),
      std::enable_shared_from_this<PointGeom>(T)
{
    m_shapeType = T.m_shapeType;
    m_globalID = T.m_globalID;
    m_coordim = T.m_coordim;
}

PointGeom::~PointGeom()
{
}

void PointGeom::GetCoords(NekDouble &x, NekDouble &y, NekDouble &z)
{
    switch (m_coordim)
    {
        case 3:
            z = (*this)(2);
            /* Falls through. */
        case 2:
            y = (*this)(1);
            /* Falls through. */
        case 1:
            x = (*this)(0);
            break;
    }
}

void PointGeom::GetCoords(Array<OneD, NekDouble> &coords)
{
    switch (m_coordim)
    {
        case 3:
            coords[2] = (*this)(2);
            /* Falls through. */
        case 2:
            coords[1] = (*this)(1);
            /* Falls through. */
        case 1:
            coords[0] = (*this)(0);
            break;
    }
}

void PointGeom::UpdatePosition(NekDouble x, NekDouble y, NekDouble z)
{
    (*this)(0) = x;
    (*this)(1) = y;
    (*this)(2) = z;
}

// _this = a + b
void PointGeom::Add(PointGeom &a, PointGeom &b)
{
    (*this)(0) = a[0] + b[0];
    (*this)(1) = a[1] + b[1];
    (*this)(2) = a[2] + b[2];
    m_coordim = std::max(a.GetCoordim(), b.GetCoordim());
}

// _this = a + b
void PointGeom::Sub(PointGeom &a, PointGeom &b)
{
    (*this)(0) = a[0] - b[0];
    (*this)(1) = a[1] - b[1];
    (*this)(2) = a[2] - b[2];
    m_coordim = std::max(a.GetCoordim(), b.GetCoordim());
}

/// \brief _this = a x b
void PointGeom::Mult(PointGeom &a, PointGeom &b)
{
    (*this)(0) = a[1] * b[2] - a[2] * b[1];
    (*this)(1) = a[2] * b[0] - a[0] * b[2];
    (*this)(2) = a[0] * b[1] - a[1] * b[0];
    m_coordim = 3;
}

/// \brief _this  = rotation of a by angle 'angle' around axis dir
void PointGeom::Rotate(PointGeom& a, int dir, NekDouble angle)
{
    switch(dir)
    {
    case 0:
        {
            NekDouble yrot = cos(angle)*a.y() - sin(angle)*a.z();
            NekDouble zrot = sin(angle)*a.y() + cos(angle)*a.z();
            
            (*this)(0) = a.x();
            (*this)(1) = yrot;
            (*this)(2) = zrot;
        }
        break;
    case 1:
        {
            NekDouble zrot = cos(angle)*a.z() - sin(angle)*a.x();
            NekDouble xrot = sin(angle)*a.z() + cos(angle)*a.x();
            
            (*this)(0) = xrot;
            (*this)(1) = a.y(); 
            (*this)(2) = zrot;
        }
        break;
    case 2:
        {
            NekDouble xrot = cos(angle)*a.x() - sin(angle)*a.y();
            NekDouble yrot = sin(angle)*a.x() + cos(angle)*a.y();
            
            (*this)(0) = xrot;
            (*this)(1) = yrot;
            (*this)(2) = a.z(); 
        }
        break;
    }
}
    
///  \brief return distance between this and input a
NekDouble PointGeom::dist(PointGeom &a)
{
    return sqrt((x() - a.x()) * (x() - a.x()) + (y() - a.y()) * (y() - a.y()) +
                (z() - a.z()) * (z() - a.z()));
}

/// \brief retun the dot product between this and input a 
NekDouble PointGeom::dot(PointGeom &a)
{
    return (x() * a.x() + y() * a.y() + z() * a.z());
}

/// Determine equivalence by the ids.  No matter what the position,
/// if the ids are the same, then they are equivalent, and vice versa.
bool operator==(const PointGeom &x, const PointGeom &y)
{
    return (x.m_globalID == y.m_globalID);
}

bool operator==(const PointGeom &x, const PointGeom *y)
{
    return (x.m_globalID == y->m_globalID);
}

bool operator==(const PointGeom *x, const PointGeom &y)
{
    return (x->m_globalID == y.m_globalID);
}

bool operator!=(const PointGeom &x, const PointGeom &y)
{
    return (x.m_globalID != y.m_globalID);
}

bool operator!=(const PointGeom &x, const PointGeom *y)
{
    return (x.m_globalID != y->m_globalID);
}

bool operator!=(const PointGeom *x, const PointGeom &y)
{
    return (x->m_globalID != y.m_globalID);
}

PointGeomSharedPtr PointGeom::v_GetVertex(int i) const
{
    ASSERTL0(i == 0, "Index other than 0 is meaningless.");
    // shared_this_ptr() returns const PointGeom, which cannot be
    // returned.
    return PointGeomSharedPtr(new PointGeom(*this));
}

void PointGeom::v_GenGeomFactors()
{
}

}
}
