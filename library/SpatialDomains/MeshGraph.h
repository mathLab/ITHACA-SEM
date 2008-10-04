////////////////////////////////////////////////////////////////////////////////
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

#include <SpatialDomains/InterfaceComponent.h>
#include <SpatialDomains/Equation.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/HexGeom.h>

#include <SpatialDomains/Curve.hpp>


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
            eNodal,
            eExpansionTypeSize
        };

        // Keep this consistent with the enums in ExpansionType.
        // This is used in the BC file to specify the expansion type.
        const std::string kExpansionTypeStr[] = 
        {
            "NOTYPE",
            "MODIFIED",
            "ORTHOGONAL",
            "NODAL"
        };

        typedef boost::shared_ptr<InterfaceComponent> SharedInterfaceCompPtr;
        typedef std::vector< VertexComponentSharedPtr >  VertexVector;
        typedef std::list< SharedInterfaceCompPtr > InterfaceCompList;

        typedef boost::shared_ptr< GeometryVector > Composite;
        typedef std::vector< Composite >            CompositeVector;
        typedef std::vector< Composite >::iterator  CompositeVectorIter;

        struct ElementEdge
        {
            GeometrySharedPtr m_Element;
            int m_EdgeIndx;
        };

        struct ElementFace
        {
            GeometrySharedPtr m_Face;
            int m_FaceIndx;
        };

        typedef boost::shared_ptr<ElementEdge> ElementEdgeSharedPtr;
        typedef std::vector<ElementEdgeSharedPtr> ElementEdgeVector;
        typedef boost::shared_ptr<ElementEdgeVector> ElementEdgeVectorSharedPtr;


        typedef boost::shared_ptr<ElementFace> ElementFaceSharedPtr;
        typedef std::vector<ElementFaceSharedPtr> ElementFaceVector;
        typedef boost::shared_ptr<ElementFaceVector> ElementFaceVectorSharedPtr;

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

            static boost::shared_ptr<MeshGraph> Read(std::string &infilename);
            virtual void ReadGeometry(std::string &infilename);
            virtual void ReadGeometry(TiXmlDocument &doc);
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

            /// \brief Dimension of the space (can be a 1D curve in 3D space).
            inline int GetSpaceDimension(void) const
            {
                return m_SpaceDimension;
            }

            inline const int GetNvertices() const 
            {
                return int(m_vertset.size());
            }
            
            struct FieldDefinitions
            {

            FieldDefinitions(GeometryVector &elements,
                             std::vector<LibUtilities::BasisType> &basis,
                             std::vector<unsigned int> &numModes,
                             std::vector<unsigned int> &numFields) : 
                    m_Elements(elements),
                    m_Basis(basis),
                    m_NumModes(numModes),
                    m_Fields(numFields)
                {};
                GeometryVector m_Elements;
                std::vector<LibUtilities::BasisType> m_Basis;
                std::vector<unsigned int> m_NumModes;
                std::vector<unsigned int> m_Fields;
            };
            
            void Write(std::string &outfilename, FieldDefinitions &fielddefs, std::vector<double> &fielddata);
            void Import(std::string &infilename, std::vector<FieldDefinitions> &fielddefs, std::vector<std::vector<double> > &fielddata);
            
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
            void ReadCurves(TiXmlDocument &doc);
            void ReadCurves(std::string &infilename);
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
            VertexVector            m_vertset;
            InterfaceCompList       m_icomps;

            CurveVector             m_curvededges;
            CurveVector             m_curvedfaces;

            SegGeomVector           m_seggeoms;

            TriGeomVector           m_trigeoms;
            QuadGeomVector          m_quadgeoms;
            TetGeomVector           m_tetgeoms;
            PyrGeomVector           m_pyrgeoms;
            PrismGeomVector         m_prismgeoms;
            HexGeomVector           m_hexgeoms;

            int m_MeshDimension;
            int m_SpaceDimension;

            CompositeVector m_MeshCompositeVector;
            CompositeVector m_Domain;
            ExpansionVector m_ExpansionVector;
            
        };

        
        typedef boost::shared_ptr<MeshGraph> MeshGraphSharedPtr;


    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH_H

//
// $Log: MeshGraph.h,v $
// Revision 1.29  2008/09/09 14:20:30  sherwin
// Updated to handle curved edges (first working version)
//
// Revision 1.28  2008/08/26 02:19:39  ehan
// Added struct element face and related shared pointers.
//
// Revision 1.27  2008/07/09 23:41:36  ehan
// Added edge component and face component to the curve reader.
//
// Revision 1.26  2008/07/08 18:58:34  ehan
// Added curve reader.
//
// Revision 1.25  2008/06/30 19:34:54  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.24  2008/06/12 23:27:57  delisi
// Removed MeshGraph.h include from SegGeom.h, to get rid of circular includes. Now can use typedefs from SegGeom.h instead of repeating it in MeshGraph.h.
//
// Revision 1.23  2008/06/11 23:25:29  ehan
// Fixed error : ‘SegGeomVector’ does not name a type
//
// Revision 1.22  2008/06/09 21:33:04  jfrazier
// Moved segment vector to base MeshGraph class since it is used by all derived types.
//
// Revision 1.21  2008/05/29 21:19:23  delisi
// Added the Write(...) and Import(...) functions which write and read XML files for output.
//
// Revision 1.20  2008/03/18 14:14:49  pvos
// Update for nodal triangular helmholtz solver
//
// Revision 1.19  2007/12/11 21:51:52  jfrazier
// Updated 2d components so elements could be retrieved from edges.
//
// Revision 1.18  2007/12/11 18:59:59  jfrazier
// Updated meshgraph so that a generic read could be performed and the proper type read (based on dimension) will be returned.
//
// Revision 1.17  2007/12/06 22:47:44  pvos
// 2D Helmholtz solver updates
//
// Revision 1.16  2007/12/04 02:54:35  jfrazier
// Removed unused declaration.
//
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
