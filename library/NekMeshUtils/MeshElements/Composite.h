////////////////////////////////////////////////////////////////////////////////
//
//  File: Composite.h
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
//  Description: Mesh Composite Object.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_MESHELEMENTS_COMPOSITE
#define NEKMESHUTILS_MESHELEMENTS_COMPOSITE

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>
#include <NekMeshUtils/MeshElements/Element.h>

namespace Nektar
{
namespace NekMeshUtils
{
/**
 * @brief A composite is a collection of elements.
 *
 * All elements should be of the same type, i.e. have the same tag.
 */
class Composite
{
public:
    NEKMESHUTILS_EXPORT Composite() : m_reorder(true)
    {
    }

    /// ID of composite.
    unsigned int m_id;
    /// Element type tag.
    std::string m_tag;
    /// boundary label
    std::string m_label;
    /// Determines whether items can be reordered.
    bool m_reorder;
    /// List of elements in this composite.
    std::vector<ElementSharedPtr> m_items;
};

/// Shared pointer to a composite.
typedef std::shared_ptr<Composite> CompositeSharedPtr;
/// Container of composites; key is the composite id, value is the
/// composite.
typedef std::map<unsigned int, CompositeSharedPtr> CompositeMap;
}
}

#endif
