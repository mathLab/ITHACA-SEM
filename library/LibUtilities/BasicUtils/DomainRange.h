////////////////////////////////////////////////////////////////////////////////
//
//  File:  DomainRange.h
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
//  Description:  Finds Min and Max X,Y,Z points of a specific domain
//
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_DOMAINRANGE_HPP
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_DOMAINRANGE_HPP

namespace Nektar 
{
namespace LibUtilities
{
    
// set restriction on domain range for post-processing.
struct DomainRange
{
    bool m_doXrange;
    NekDouble m_xmin;
    NekDouble m_xmax;
    bool m_doYrange;
    NekDouble m_ymin;
    NekDouble m_ymax;
    bool m_doZrange;
    NekDouble m_zmin;
    NekDouble m_zmax;

    bool m_checkShape;
    LibUtilities::ShapeType m_shapeType;
};

typedef std::shared_ptr<DomainRange> DomainRangeShPtr;
static DomainRangeShPtr NullDomainRangeShPtr;


}

}
#endif
