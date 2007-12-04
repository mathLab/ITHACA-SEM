
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/MeshGraph.h,v $
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH_H

#include <cstdlib>
#include <fstream>

#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/InterfaceComponent.h>
#include <SpatialDomains/Equation.h>

class TiXmlDocument;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum ExpansionType
        {
            eNoExpansionType,
            eModified,
            eOrthogonal,
            eExpansionTypeSize
        };

        // Keep this consistent with the enums in ExpansionType.
        // This is used in the BC file to specify the expansion type.
        const std::string kExpansionTypeStr[] = 
        {
            "NOTYPE",
            "MODIFIED",
            "ORTHOGONAL"
        };

        typedef boost::shared_ptr<InterfaceComponent> SharedInterfaceCompPtr;
        typedef std::vector< VertexComponentSharedPtr >  VertexVector;
        typedef std::list< SharedInterfaceCompPtr > InterfaceCompList;

        typedef boost::shared_ptr< GeometryVector > Composite;
        typedef std::vector< Composite >            CompositeVector;
        typedef std::vector< Composite >::iterator  CompositeVectorIter;

        struct Expansion
        {
            Expansion(GeometrySharedPtr geomShPtr, const Equation &numModesEqn, ExpansionType expansionType):
                m_GeomShPtr(geomShPtr),
                m_NumModesEqn(numModesEqn),
                m_ExpansionType(expansionType)
            {};
            GeometrySharedPtr m_GeomShPtr;
            Equation m_NumModesEqn;
            ExpansionType m_ExpansionType;
        };

        typedef boost::shared_ptr<Expansion> ExpansionShPtr;
        typedef std::vector<ExpansionShPtr> ExpansionVector;
        typedef std::vector<ExpansionShPtr>::iterator ExpansionVectorIter;

        class MeshGraph
        {
        public:
            MeshGraph();
            virtual ~MeshGraph();

            void ReadGeometry(std::string &infilename);
            void ReadGeometry(TiXmlDocument &doc);
            void ReadExpansions(std::string &infilename);
            void ReadExpansions(TiXmlDocument &doc);

            inline VertexComponentSharedPtr GetVertex(int id)
            {
                return m_vertset[id];
            }

            /// \brief Dimension of the mesh (can be a 1D curve in 3D space).
            inline int GetMeshDimension(void) const
            {
                return m_MeshDimension;
            }

            /// \brief Dimension of the mesh (can be a 1D curve in 3D space).
            inline int GetSpaceDimension(void) const
            {
                return m_SpaceDimension;
            }

            void Write(std::string &outfilename);

            GeometrySharedPtr GetCompositeItem(int whichComposite, int whichItem);
            Composite GetComposite(int whichComposite) const
            {
                Composite returnval;
                if (whichComposite >= 0 && whichComposite < int(m_MeshCompositeVector.size()))
                {
                    returnval = m_MeshCompositeVector[whichComposite];
                }

                return returnval;
            }

            const CompositeVector &GetDomain(void) const
            {
                return m_Domain;
            }

            void ReadDomain(TiXmlDocument &doc);
            void GetCompositeList(const std::string &compositeStr, CompositeVector &compositeVector) const;
            LibUtilities::BasisKey GetBasisKey(ExpansionShPtr in, const int flag = 0);

            ExpansionShPtr GetExpansion(GeometrySharedPtr geom)
            {
                ExpansionVectorIter iter;
                ExpansionShPtr returnval;

                for (iter = m_ExpansionVector.begin(); iter!=m_ExpansionVector.end(); ++iter)
                {
                    if ((*iter)->m_GeomShPtr == geom)
                    {
                        returnval = *iter;
                        break;
                    }
                }

                return returnval;
            }

            const ExpansionVector &GetExpansions(void) const
            {
                return m_ExpansionVector;
            }

       protected:
            VertexVector  m_vertset;
            InterfaceCompList m_icomps;

            int m_MeshDimension;
            int m_SpaceDimension;

            CompositeVector m_MeshCompositeVector;
            CompositeVector m_Domain;
            ExpansionVector m_ExpansionVector;
       };
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH_H

//
// $Log: MeshGraph.h,v $
// Revision 1.15  2007/11/07 20:31:04  jfrazier
// Added new expansion list to replace the expansion composite list.
//
// Revision 1.14  2007/09/20 22:25:06  jfrazier
// Added expansion information to meshgraph class.
//
// Revision 1.13  2007/09/03 17:05:01  jfrazier
// Cleanup and addition of composite range in domain specification.
//
// Revision 1.12  2007/08/11 23:38:48  sherwin
// Update for full working version of Helmholtz1D
//
// Revision 1.11  2007/07/25 11:01:57  sherwin
// Added GetDomain methods
//
// Revision 1.10  2007/07/24 16:52:09  jfrazier
// Added domain code.
//
// Revision 1.9  2007/07/22 23:04:23  bnelson
// Backed out Nektar::ptr.
//
// Revision 1.8  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.7  2007/06/07 23:55:24  jfrazier
// Intermediate revisions to add parsing for boundary conditions file.
//
// Revision 1.6  2007/01/18 20:59:28  sherwin
// Before new configuration
//
// Revision 1.5  2006/09/26 23:41:53  jfrazier
// Updated to account for highest level NEKTAR tag and changed the geometry tag to GEOMETRY.
//
// Revision 1.4  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.3  2006/06/01 14:58:53  kirby
// *** empty log message ***
//
// Revision 1.2  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.1  2006/05/04 18:59:02  kirby
// *** empty log message ***
//
// Revision 1.16  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.15  2006/04/04 23:12:37  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.14  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.13  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.12  2006/03/12 07:42:03  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.11  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.10  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
