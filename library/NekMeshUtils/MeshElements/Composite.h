////////////////////////////////////////////////////////////////////////////////
//
//  File: Face.h
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
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_MESHELEMENTS_COMPOSITE
#define NekMeshUtils_MESHELEMENTS_COMPOSITE

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

    /**
     * @brief Generate a Nektar++ string describing the composite.
     *
     * The list of composites may include individual element IDs or ranges of
     * element IDs.
     */
    NEKMESHUTILS_EXPORT std::string GetXmlString(bool doSort = true)
    {
        std::stringstream st;
        std::vector<ElementSharedPtr>::iterator it;
        bool range = false;
        int vId    = m_items[0]->GetId();
        int prevId = vId;

        st << " " << m_tag << "[" << vId;

        for (it = m_items.begin() + 1; it != m_items.end(); ++it)
        {
            // store previous element ID and get current one
            prevId = vId;
            vId    = (*it)->GetId();

            // continue an already started range
            if (prevId > -1 && vId == prevId + 1)
            {
                range = true;
                // if this is the last element, it's the end of a range, so
                // write
                if (*it == m_items.back())
                {
                    st << "-" << vId;
                }
                continue;
            }

            // terminate a range, if present
            if (range)
            {
                st << "-" << prevId;
                range = false;
            }

            // write what will be either a single entry or start of new range
            st << "," << vId;
        }
        // terminate
        st << "] ";
        return st.str();
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
typedef boost::shared_ptr<Composite> CompositeSharedPtr;
/// Container of composites; key is the composite id, value is the
/// composite.
typedef std::map<unsigned int, CompositeSharedPtr> CompositeMap;
}
}

#endif
