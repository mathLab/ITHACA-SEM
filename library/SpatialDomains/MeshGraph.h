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
            eGLL_Lagrange,
            eGLL_Lagrange_SEM,
            eExpansionTypeSize
        };

        // Keep this consistent with the enums in ExpansionType.
        // This is used in the BC file to specify the expansion type.
        const std::string kExpansionTypeStr[] =
        {
            "NOTYPE",
            "MODIFIED",
            "ORTHOGONAL",
            "GLL_LAGRANGE",
            "GLL_LAGRANGE_SEM"
        };

        typedef boost::shared_ptr< InterfaceComponent > SharedInterfaceCompPtr;
        typedef std::vector< VertexComponentSharedPtr > VertexVector;
        typedef std::list< SharedInterfaceCompPtr >     InterfaceCompList;

        typedef boost::shared_ptr< GeometryVector >     Composite;
        typedef std::vector< Composite >                CompositeVector;
        typedef std::vector< Composite >::iterator      CompositeVectorIter;

        struct ElementEdge
        {
            GeometrySharedPtr m_Element;
            int m_EdgeIndx;
        };

        struct ElementFace
        {
            GeometrySharedPtr m_Element;
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
            Expansion(GeometrySharedPtr geomShPtr,
                      const LibUtilities::BasisKeyVector basiskeyvec):
                m_geomShPtr(geomShPtr),
                m_basisKeyVector(basiskeyvec)
            {
            }

            GeometrySharedPtr             m_geomShPtr;
            LibUtilities::BasisKeyVector  m_basisKeyVector;
        };

        typedef boost::shared_ptr<Expansion> ExpansionShPtr;
        typedef std::vector<ExpansionShPtr>  ExpansionVector;
        typedef std::vector<ExpansionShPtr>::iterator ExpansionVectorIter;

        static std::vector<NekDouble> NullNekDoubleVector;
        static std::vector<LibUtilities::PointsType> NullPointsTypeVector;
        static std::vector<unsigned int> NullUnsignedIntVector;

        struct FieldDefinitions
        {
        FieldDefinitions(SpatialDomains::GeomShapeType shapeType,
                         std::vector<unsigned int> &elementIDs,// vector[2]
                         std::vector<LibUtilities::BasisType> &basis,
                         bool uniOrder,
                         // UniOrder = vector[dimension] - MixOrder
                         //          = vector[element*dimension]
                         std::vector<unsigned int> &numModes,
                         std::vector<std::string>  &fields,
                         int NumHomoDir = 0,
                         std::vector<NekDouble> &HomoLengths =
                         NullNekDoubleVector,
                         std::vector<LibUtilities::PointsType> &points =
                         NullPointsTypeVector,
                         bool pointsDef = false,
                         std::vector<unsigned int> &numPoints =
                         NullUnsignedIntVector,
                         bool numPointsDef = false):
            m_shapeType(shapeType),
                m_elementIDs(elementIDs),
                m_basis(basis),
                m_numHomogeneousDir(NumHomoDir),
                m_homogeneousLengths(HomoLengths),
                m_points(points),
                m_pointsDef(pointsDef),
                m_uniOrder(uniOrder),
                m_numModes(numModes),
                m_numPoints(numPoints),
                m_numPointsDef(numPointsDef),
                m_fields(fields)
            {
            }

            SpatialDomains::GeomShapeType         m_shapeType;
            std::vector<unsigned int>             m_elementIDs;
            std::vector<LibUtilities::BasisType>  m_basis;
            int                                   m_numHomogeneousDir;
            std::vector<NekDouble>                m_homogeneousLengths;
            /// Define the type of points per direction.
            std::vector<LibUtilities::PointsType> m_points;
            bool                                  m_pointsDef;
            /// Define order of the element group.
            /// * UniOrder: same order for each element
            /// * MixOrder: definition of a different order for each element.
            bool                                  m_uniOrder;
            /// Define number of modes per direction.
            std::vector<unsigned int>             m_numModes;
            std::vector<unsigned int>             m_numPoints;
            bool                                  m_numPointsDef;
            std::vector<std::string>              m_fields;
        };

        typedef boost::shared_ptr<FieldDefinitions> FieldDefinitionsSharedPtr;

        typedef std::map<std::string, std::string> GeomInfoMap;


        /// Base class for a spectral/hp element mesh.
        class MeshGraph
        {
        public:
            MeshGraph();
            virtual ~MeshGraph();

            static boost::shared_ptr<MeshGraph> Read(
                    const std::string& infilename);

            virtual void ReadGeometry(const std::string& infilename);
            virtual void ReadGeometry(TiXmlDocument &doc);

            /// Read geometric information from a file.
            void ReadGeometryInfo(const std::string &infilename);

            /// Read geometric information from an XML document.
            void ReadGeometryInfo(TiXmlDocument &doc);

            void ReadExpansions(const std::string &infilename);

            void ReadExpansions(TiXmlDocument &doc);

            inline VertexComponentSharedPtr GetVertex(int id);

            /// \brief Dimension of the mesh (can be a 1D curve in 3D space).
            inline int GetMeshDimension(void) const;

            /// \brief Dimension of the space (can be a 1D curve in 3D space).
            inline int GetSpaceDimension(void) const;

            inline const int GetNvertices() const;

            int CheckFieldDefinition(
                    const FieldDefinitionsSharedPtr  &fielddefs);

            void Write(const std::string &outFile,
                       std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                       std::vector<std::vector<double> >      &fielddata);

            /// This function imports the input xml file. It defines the fields
            /// and their data.
            void Import(const std::string& infilename,
                    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                    std::vector<std::vector<double> > &fielddata);

            /// This function imports the definition of the fields.
            void ImportFieldDefs(TiXmlDocument &doc,
                    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                    bool expChild);

            /// This function imports the data fileds.
            void ImportFieldData(TiXmlDocument &doc,
                    const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                    std::vector<std::vector<double> > &fielddata);

            GeometrySharedPtr GetCompositeItem(int whichComposite,
                    int whichItem);

            inline Composite GetComposite(int whichComposite) const;

            inline const CompositeVector &GetDomain(void) const;

            void ReadDomain(TiXmlDocument &doc);
            void ReadCurves(TiXmlDocument &doc);
            void ReadCurves(std::string &infilename);
            void GetCompositeList(const std::string &compositeStr,
                    CompositeVector &compositeVector) const;

            inline ExpansionShPtr GetExpansion(GeometrySharedPtr geom);

            /// This function sets the expansion giving the definition of the
            /// field and the quadrature points.
            void SetExpansions(
                    std::vector<SpatialDomains::FieldDefinitionsSharedPtr>
                                                                &fielddef,
                    std::vector< std::vector<LibUtilities::PointsType> >
                                                                &pointstype);

            /// This function sets the expansion giving the definition of the
            /// field. The quadrature points type and number is defined as
            /// default.
            void SetExpansions(
                    std::vector<SpatialDomains::FieldDefinitionsSharedPtr>
                                                                &fielddef);

            /// Sets the basis key for all expansions of the given shape.
            void SetBasisKey(SpatialDomains::GeomShapeType shape,
                             LibUtilities::BasisKeyVector &keys);

            inline const ExpansionVector &GetExpansions(void) const;

            inline const bool CheckForGeomInfo(std::string parameter);

            inline const std::string GetGeomInfo(std::string parameter);

       protected:
            VertexVector            m_vertSet;
            InterfaceCompList       m_iComps;

            CurveVector             m_curvedEdges;
            CurveVector             m_curvedFaces;

            SegGeomVector           m_segGeoms;

            TriGeomVector           m_triGeoms;
            QuadGeomVector          m_quadGeoms;
            TetGeomVector           m_tetGeoms;
            PyrGeomVector           m_pyrGeoms;
            PrismGeomVector         m_prismGeoms;
            HexGeomVector           m_hexGeoms;

            int                     m_meshDimension;
            int                     m_spaceDimension;

            CompositeVector         m_meshCompositeVector;
            CompositeVector         m_domain;
            ExpansionVector         m_expansionVector;

            GeomInfoMap             m_geomInfo;
        };

        typedef boost::shared_ptr<MeshGraph> MeshGraphSharedPtr;

        inline VertexComponentSharedPtr MeshGraph::GetVertex(int id)
        {
            return m_vertSet[id];
        }

        /// \brief Dimension of the mesh (can be a 1D curve in 3D space).
        inline int MeshGraph::GetMeshDimension(void) const
        {
            return m_meshDimension;
        }

        /// \brief Dimension of the space (can be a 1D curve in 3D space).
        inline int MeshGraph::GetSpaceDimension(void) const
        {
            return m_spaceDimension;
        }

        inline const int MeshGraph::GetNvertices() const
        {
            return int(m_vertSet.size());
        }

        inline Composite MeshGraph::GetComposite(int whichComposite) const
        {
            Composite returnval;
            if (whichComposite >= 0
                    && whichComposite < int(m_meshCompositeVector.size()))
            {
                returnval = m_meshCompositeVector[whichComposite];
            }

            return returnval;
        }

        inline const CompositeVector &MeshGraph::GetDomain(void) const
        {
            return m_domain;
        }

        inline ExpansionShPtr MeshGraph::GetExpansion(GeometrySharedPtr geom)
        {
            ExpansionVectorIter iter;
            ExpansionShPtr returnval;

            for (iter = m_expansionVector.begin();
                                    iter!=m_expansionVector.end(); ++iter)
            {
                if ((*iter)->m_geomShPtr == geom)
                {
                    returnval = *iter;
                    break;
                }
            }
            return returnval;
        }

        inline const ExpansionVector &MeshGraph::GetExpansions(void) const
        {
            return m_expansionVector;
        }

        inline const bool MeshGraph::CheckForGeomInfo(std::string parameter)
        {
            return m_geomInfo.find(parameter) != m_geomInfo.end();
        }

        inline const std::string MeshGraph::GetGeomInfo(std::string parameter)
        {
            ASSERTL1(m_geomInfo.find(parameter) != m_geomInfo.end(),
                    "Parameter " + parameter + " does not exist.");
            return m_geomInfo[parameter];
        }
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH_H

//
// $Log: MeshGraph.h,v $
// Revision 1.40  2009/12/16 21:09:13  bnelson
// Updated file read methods to take const std::string& instead of std::string&
//
// Revision 1.39  2009/12/15 18:09:02  cantwell
// Split GeomFactors into 1D, 2D and 3D
// Added generation of tangential basis into GeomFactors
// Updated ADR2DManifold solver to use GeomFactors for tangents
// Added <GEOMINFO> XML session section support in MeshGraph
// Fixed const-correctness in VmathArray
// Cleaned up LocalRegions code to generate GeomFactors
// Removed GenSegExp
// Temporary fix to SubStructuredGraph
// Documentation for GlobalLinSys and GlobalMatrix classes
//
// Revision 1.38  2009/11/22 19:43:32  bnelson
// Updating formatting.
//
// Revision 1.37  2009/11/18 22:31:46  bnelson
// Changed Write parameter list to accept a const string& as a first parameter.
//
// Revision 1.36  2009/09/24 10:57:54  cbiotto
// Updates for variable order expansions
//
// Revision 1.35  2009/08/19 14:13:34  claes
// Removed Gauss-Kronrod parts
//
// Revision 1.34  2009/06/15 01:59:21  claes
// Gauss-Kronrod updates
//
// Revision 1.33  2009/05/01 13:23:21  pvos
// Fixed various bugs
//
// Revision 1.32  2009/04/20 16:13:23  sherwin
// Modified Import and Write functions and redefined how Expansion is used
//
// Revision 1.31  2009/01/12 10:26:59  pvos
// Added input tags for nodal expansions
//
// Revision 1.30  2008/10/04 19:32:47  sherwin
// Added SharedPtr Typedef and replaced MeshDimension with SpaceDimension
//
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
