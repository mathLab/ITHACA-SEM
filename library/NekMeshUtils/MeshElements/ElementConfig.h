////////////////////////////////////////////////////////////////////////////////
//
//  File: ElementConfig.h
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
//  Description: Mesh elementconfig.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_MESHELEMENTS_ELEMENTCONFIG
#define NEKMESHUTILS_MESHELEMENTS_ELEMENTCONFIG

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief Basic information about an element.
 *
 * ElmtConfig contains four member variables which denote the
 * properties of an element when it is created.
 */
struct ElmtConfig
{
    ElmtConfig(LibUtilities::ShapeType  pE,
               unsigned int             pOrder,
               bool                     pFn,
               bool                     pVn,
               bool                     pReorient = true,
               LibUtilities::PointsType pECt = LibUtilities::ePolyEvenlySpaced,
               LibUtilities::PointsType pFCt = LibUtilities::ePolyEvenlySpaced)
        : m_e(pE), m_faceNodes(pFn), m_volumeNodes(pVn), m_order(pOrder),
          m_reorient(pReorient), m_edgeCurveType(pECt), m_faceCurveType(pFCt)
    {
    }

    ElmtConfig(ElmtConfig const &p)
        : m_e(p.m_e), m_faceNodes(p.m_faceNodes),
          m_volumeNodes(p.m_volumeNodes), m_order(p.m_order),
          m_reorient(p.m_reorient), m_edgeCurveType(p.m_edgeCurveType),
          m_faceCurveType(p.m_faceCurveType)
    {
    }

    ElmtConfig()
    {
    }

    ElmtConfig& operator=(const ElmtConfig &) = default;

    /// Element type (e.g. triangle, quad, etc).
    LibUtilities::ShapeType m_e;
    /// Denotes whether the element contains face nodes. For 2D elements, if
    /// this is true then the element contains interior nodes.
    bool m_faceNodes;
    /// Denotes whether the element contains volume (i.e. interior) nodes. These
    /// are not supported by either the mesh converter or Nektar++ but are
    /// included for completeness and are required for some output modules
    /// (e.g. Gmsh).
    bool m_volumeNodes;
    /// Order of the element.
    unsigned int m_order;
    /// Denotes whether the element needs to be re-orientated for a spectral
    /// element framework.
    bool m_reorient;
    /// Distribution of points in edges.
    LibUtilities::PointsType m_edgeCurveType;
    /// Distribution of points in faces.
    LibUtilities::PointsType m_faceCurveType;
};

NEKMESHUTILS_EXPORT bool operator==(ElmtConfig const &c1, ElmtConfig const &c2);

}
}
#endif
