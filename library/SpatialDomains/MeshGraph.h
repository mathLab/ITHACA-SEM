////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph.h
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
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH_H

#include <boost/unordered_map.hpp>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/HexGeom.h>

#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

class TiXmlDocument;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum ExpansionType
        {
            eNoExpansionType,
            eModified,
            eModifiedQuadPlus1,
            eModifiedQuadPlus2,
            eOrthogonal,
            eGLL_Lagrange,
            eGLL_Lagrange_SEM,
            eGauss_Lagrange,
            eGauss_Lagrange_SEM,
            eFourier,
            eFourierSingleMode,
            eFourierHalfModeRe,
            eFourierHalfModeIm,
            eChebyshev,
            eFourierChebyshev,
            eChebyshevFourier,
            eFourierModified,
            eExpansionTypeSize
        };

        // Keep this consistent with the enums in ExpansionType.
        // This is used in the BC file to specify the expansion type.
        const std::string kExpansionTypeStr[] =
        {
            "NOTYPE",
            "MODIFIED",
            "MODIFIEDQUADPLUS1",
            "MODIFIEDQUADPLUS2",
            "ORTHOGONAL",
            "GLL_LAGRANGE",
            "GLL_LAGRANGE_SEM",
            "GAUSS_LAGRANGE",
            "GAUSS_LAGRANGE_SEM",
            "FOURIER",
            "FOURIERSINGLEMODE",
            "FOURIERHALFMODERE",
            "FOURIERHALFMODEIM",
            "CHEBYSHEV",
            "FOURIER-CHEBYSHEV",
            "CHEBYSHEV-FOURIER",
            "FOURIER-MODIFIED"
        };

        class InterfaceComponent;
        typedef boost::shared_ptr< InterfaceComponent > SharedInterfaceCompPtr;
        typedef std::vector< VertexComponentSharedPtr > VertexVector;
        typedef std::map<int, VertexComponentSharedPtr> VertexMap;
        typedef std::list< SharedInterfaceCompPtr >     InterfaceCompList;

        typedef boost::shared_ptr< GeometryVector >     Composite;
        typedef std::map<int, Composite>                CompositeMap;
        typedef std::map<int, Composite>::iterator      CompositeMapIter;

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
        typedef std::map<int, ExpansionShPtr> ExpansionMap;
        typedef std::map<int, ExpansionShPtr>::iterator ExpansionMapIter;
        typedef std::map<int, ExpansionShPtr>::const_iterator ExpansionMapConstIter;

        typedef boost::shared_ptr<ExpansionMap> ExpansionMapShPtr;
        typedef std::map<std::string, ExpansionMapShPtr> ExpansionMapShPtrMap;

        static std::vector<NekDouble> NullNekDoubleVector;
        static std::vector<LibUtilities::PointsType> NullPointsTypeVector;
        static std::vector<unsigned int> NullUnsignedIntVector;

        struct FieldDefinitions
        {
            FieldDefinitions(SpatialDomains::GeomShapeType shapeType,
                             const std::vector<unsigned int> &elementIDs,// vector[2]
                             const std::vector<LibUtilities::BasisType> &basis,
                             bool uniOrder,
                             // UniOrder = vector[dimension] - MixOrder
                             //          = vector[element*dimension]
                             const std::vector<unsigned int> &numModes,
                             const std::vector<std::string>  &fields,
                             int NumHomoDir = 0,
                             const std::vector<NekDouble> &HomoLengths =
                             NullNekDoubleVector,
							 const std::vector<unsigned int> &HomoZIDs =
                             NullUnsignedIntVector,
							 const std::vector<unsigned int> &HomoYIDs =
                             NullUnsignedIntVector,
                             const std::vector<LibUtilities::PointsType> &points =
                             NullPointsTypeVector,
                             bool pointsDef = false,
                             const std::vector<unsigned int> &numPoints =
                             NullUnsignedIntVector,
                             bool numPointsDef = false):
                m_shapeType(shapeType),
                    m_elementIDs(elementIDs),
                    m_basis(basis),
                    m_numHomogeneousDir(NumHomoDir),
                    m_homogeneousLengths(HomoLengths),
			        m_homogeneousZIDs(HomoZIDs),
			        m_homogeneousYIDs(HomoYIDs),
                    m_points(points),
                    m_pointsDef(pointsDef),
                    m_uniOrder(uniOrder),
                    m_numModes(numModes),
                    m_numPoints(numPoints),
                    m_numPointsDef(numPointsDef),
                    m_fields(fields)
                {
                }

                SpatialDomains::GeomShapeType			m_shapeType;
                std::vector<unsigned int>					m_elementIDs;
                std::vector<LibUtilities::BasisType>	m_basis;
                int										m_numHomogeneousDir;
            std::vector<NekDouble>					m_homogeneousLengths;
            std::vector<unsigned int>					m_homogeneousZIDs;
            std::vector<unsigned int>					m_homogeneousYIDs;
            
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
                SPATIAL_DOMAINS_EXPORT MeshGraph();

                SPATIAL_DOMAINS_EXPORT MeshGraph(
                        unsigned int meshDimension,
                        unsigned int spaceDimension);

                SPATIAL_DOMAINS_EXPORT MeshGraph(
                        const LibUtilities::SessionReaderSharedPtr &pSession);

                SPATIAL_DOMAINS_EXPORT virtual ~MeshGraph();


                /* ---- Mesh Reading routines ---- */
                SPATIAL_DOMAINS_EXPORT static boost::shared_ptr<MeshGraph> Read(
                        const LibUtilities::SessionReaderSharedPtr &pSession);

                /// \todo Remove updated routine
                SPATIAL_DOMAINS_EXPORT static boost::shared_ptr<MeshGraph> Read(
                        const std::string& infilename,
                        bool pReadExpansions = true);


                /// Read will read the meshgraph vertices given a filename.
                SPATIAL_DOMAINS_EXPORT virtual void ReadGeometry(
                        const std::string& infilename);

                /// Read will read the meshgraph vertices given a TiXmlDocument.
                SPATIAL_DOMAINS_EXPORT virtual void ReadGeometry(
                        TiXmlDocument &doc);

                /// Read geometric information from a file.
                SPATIAL_DOMAINS_EXPORT void ReadGeometryInfo(
                        const std::string &infilename);

                /// Read geometric information from an XML document.
                SPATIAL_DOMAINS_EXPORT void ReadGeometryInfo(
                        TiXmlDocument &doc);

                /// Read the expansions given the XML file path.
                SPATIAL_DOMAINS_EXPORT void ReadExpansions(
                        const std::string &infilename);

                /// Read the expansions given the XML document reference.
                SPATIAL_DOMAINS_EXPORT void ReadExpansions(
                        TiXmlDocument &doc);

                SPATIAL_DOMAINS_EXPORT void ReadDomain(
                        TiXmlDocument &doc);

                SPATIAL_DOMAINS_EXPORT void ReadCurves(
                        TiXmlDocument &doc);

                SPATIAL_DOMAINS_EXPORT void ReadCurves(
                        std::string &infilename);


                /* --- FLD handling routines ---- */
                SPATIAL_DOMAINS_EXPORT void Write(
                        const std::string &outFile,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> >      &fielddata);

                /// Imports an FLD file.
                SPATIAL_DOMAINS_EXPORT void Import(
                        const std::string& infilename,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata);

                /// Imports the definition of the fields.
                SPATIAL_DOMAINS_EXPORT void ImportFieldDefs(
                        TiXmlDocument &doc,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        bool expChild);

                /// Imports the data fileds.
                SPATIAL_DOMAINS_EXPORT void ImportFieldData(
                        TiXmlDocument &doc,
                        const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata);

                SPATIAL_DOMAINS_EXPORT int CheckFieldDefinition(
                        const FieldDefinitionsSharedPtr  &fielddefs);


                /* ---- Helper functions ---- */
                /// Dimension of the mesh (can be a 1D curve in 3D space).
                inline int GetMeshDimension() const;

                /// Dimension of the space (can be a 1D curve in 3D space).
                inline int GetSpaceDimension() const;


                /* ---- Composites and Domain ---- */
                inline Composite GetComposite(int whichComposite) const;

                SPATIAL_DOMAINS_EXPORT GeometrySharedPtr GetCompositeItem(
                        int whichComposite,
                        int whichItem);

                SPATIAL_DOMAINS_EXPORT void GetCompositeList(
                        const std::string &compositeStr,
                        CompositeMap &compositeVector) const;

                inline const CompositeMap &GetDomain() const;


                /* ---- Expansions ---- */
                inline const ExpansionMap &GetExpansions();

                SPATIAL_DOMAINS_EXPORT const ExpansionMap &GetExpansions(
                        const std::string variable);

                SPATIAL_DOMAINS_EXPORT ExpansionShPtr GetExpansion(
                                                                   GeometrySharedPtr geom, const std::string variable = "DefaultVar");

                /// Sets expansions given field definitions
                SPATIAL_DOMAINS_EXPORT void SetExpansions(
                        std::vector<SpatialDomains::FieldDefinitionsSharedPtr>
                                                                &fielddef);

                /// Sets expansions given field definition, quadrature points.
                SPATIAL_DOMAINS_EXPORT void SetExpansions(
                        std::vector<SpatialDomains::FieldDefinitionsSharedPtr>
                                                                &fielddef,
                        std::vector< std::vector<LibUtilities::PointsType> >
                                                                &pointstype );

                /// This function sets the expansion #exp in map with entry #variable
                inline void SetExpansions(
                        const std::string variable,
                        ExpansionMapShPtr &exp);

                /// Sets the basis key for all expansions of the given shape.
                SPATIAL_DOMAINS_EXPORT void SetBasisKey(
                        SpatialDomains::GeomShapeType shape,
                        LibUtilities::BasisKeyVector &keys,
                        std::string var = "DefaultVar");

                inline bool  SameExpansions(
                        const std::string var1,
                        const std::string var2);

                inline bool CheckForGeomInfo(std::string parameter);

                inline const std::string GetGeomInfo(std::string parameter);

                SPATIAL_DOMAINS_EXPORT static LibUtilities::BasisKeyVector
                                        DefineBasisKeyFromExpansionType(
                        GeometrySharedPtr in,
                        ExpansionType type,
                        const int order);

                SPATIAL_DOMAINS_EXPORT LibUtilities::BasisKeyVector
                                        DefineBasisKeyFromExpansionTypeHomo(
                        GeometrySharedPtr in,
                        ExpansionType type_x,
                        ExpansionType type_y,
                        ExpansionType type_z,
                        const int nummodes_x,
                        const int nummodes_y,
                        const int nummodes_z);


                /* ---- Manipulation of mesh ---- */
                inline int GetNvertices() const;

                inline VertexComponentSharedPtr GetVertex(int id);
                /// Adds a vertex to the with the next available ID.
                SPATIAL_DOMAINS_EXPORT VertexComponentSharedPtr AddVertex(
                        NekDouble x,
                        NekDouble y,
                        NekDouble z);

                /// \brief Adds an edge between two points.  If curveDefinition is 
                /// null, then the edge is straight, otherwise it is curved according 
                /// to the curveDefinition.
                SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr AddEdge(VertexComponentSharedPtr v0, VertexComponentSharedPtr v1,
                    CurveSharedPtr curveDefinition = CurveSharedPtr());
                SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr GetEdge(unsigned int id) { return m_segGeoms[id]; }

                SPATIAL_DOMAINS_EXPORT TriGeomSharedPtr AddTriangle(SegGeomSharedPtr edges[], StdRegions::Orientation orient[]);
                SPATIAL_DOMAINS_EXPORT QuadGeomSharedPtr AddQuadrilateral(SegGeomSharedPtr edges[], StdRegions::Orientation orient[]);
                SPATIAL_DOMAINS_EXPORT TetGeomSharedPtr AddTetrahedron(TriGeomSharedPtr tfaces[TetGeom::kNtfaces]);
                SPATIAL_DOMAINS_EXPORT PyrGeomSharedPtr AddPyramid(TriGeomSharedPtr tfaces[PyrGeom::kNtfaces],
                    QuadGeomSharedPtr qfaces[PyrGeom::kNqfaces]);
                SPATIAL_DOMAINS_EXPORT PrismGeomSharedPtr AddPrism(TriGeomSharedPtr tfaces[PrismGeom::kNtfaces],
                    QuadGeomSharedPtr qfaces[PrismGeom::kNqfaces]);
                SPATIAL_DOMAINS_EXPORT HexGeomSharedPtr AddHexahedron(QuadGeomSharedPtr qfaces[HexGeom::kNqfaces]);
                // void AddExpansion(ExpansionShPtr expansion) { m_expansions[expansion->m_geomShPtr->GetGlobalID()] = expansion; }
                SPATIAL_DOMAINS_EXPORT const SegGeomMap& GetAllSegGeoms() const { return m_segGeoms; }
                SPATIAL_DOMAINS_EXPORT const TriGeomMap& GetAllTriGeoms() const { return m_triGeoms; }
                SPATIAL_DOMAINS_EXPORT const QuadGeomMap& GetAllQuadGeoms() const { return m_quadGeoms; }
                SPATIAL_DOMAINS_EXPORT const TetGeomMap& GetAllTetGeoms() const { return m_tetGeoms; }
                SPATIAL_DOMAINS_EXPORT const PyrGeomMap& GetAllPyrGeoms() const { return m_pyrGeoms; }
                SPATIAL_DOMAINS_EXPORT const PrismGeomMap& GetAllPrismGeoms() const { return m_prismGeoms; }
                SPATIAL_DOMAINS_EXPORT const HexGeomMap& GetAllHexGeoms() const { return m_hexGeoms; }

                /// Convenience method for ElVis.
                template<typename ElementType>
                const std::map<int, boost::shared_ptr<ElementType> >& GetAllElementsOfType() const;

            protected:
                LibUtilities::SessionReaderSharedPtr    m_session;
                VertexMap                               m_vertSet;
                InterfaceCompList                       m_iComps;

                CurveVector                             m_curvedEdges;
                CurveVector                             m_curvedFaces;

                SegGeomMap                              m_segGeoms;

                TriGeomMap                              m_triGeoms;
                QuadGeomMap                             m_quadGeoms;
                TetGeomMap                              m_tetGeoms;
                PyrGeomMap                              m_pyrGeoms;
                PrismGeomMap                            m_prismGeoms;
                HexGeomMap                              m_hexGeoms;

                int                                     m_meshDimension;
                int                                     m_spaceDimension;
                int                                     m_partition;
                bool                                    m_meshPartitioned;

                CompositeMap                            m_meshComposites;
                CompositeMap                            m_domain;

                ExpansionMapShPtrMap                    m_expansionMapShPtrMap;

                GeomInfoMap                             m_geomInfo;

                ExpansionMapShPtr    SetUpExpansionMap(void);
        };
        typedef boost::shared_ptr<MeshGraph> MeshGraphSharedPtr;


        /**
         *
         */
        inline int MeshGraph::GetMeshDimension(void) const
        {
            return m_meshDimension;
        }


        /**
         *
         */
        inline int MeshGraph::GetSpaceDimension(void) const
        {
            return m_spaceDimension;
        }


        /**
         *
         */
        inline Composite MeshGraph::GetComposite(int whichComposite) const
        {
            Composite returnval;
            ASSERTL0(m_meshComposites.find(whichComposite) != m_meshComposites.end(),
                    "Composite not found.");
            return m_meshComposites.find(whichComposite)->second;
        }


        /**
         *
         */
        inline const CompositeMap &MeshGraph::GetDomain() const
        {
            return m_domain;
        }


        /**
         *
         */
        inline const ExpansionMap &MeshGraph::GetExpansions()
        {
            std::string defstr = "DefaultVar";
            return GetExpansions(defstr);
        }


        /**
         *
         */
        void  MeshGraph::SetExpansions(const std::string variable, ExpansionMapShPtr &exp) 
        {
            if(m_expansionMapShPtrMap.count(variable) != 0)
            {
                ASSERTL0(false,(std::string("Expansion field is already set for variable ") + variable).c_str());
            }
            else
            {
                m_expansionMapShPtrMap[variable] = exp;
            }
        }


        /**
         *
         */
        inline bool MeshGraph::SameExpansions(const std::string var1, const std::string var2) 
        {
            ExpansionMapShPtr expVec1 = m_expansionMapShPtrMap.find(var1)->second;
            ExpansionMapShPtr expVec2 = m_expansionMapShPtrMap.find(var2)->second;

            if(expVec1.get() == expVec2.get())
            {
                return true; 
            }

            return false;
        }


        /**
         *
         */
        inline bool MeshGraph::CheckForGeomInfo(std::string parameter)
        {
            return m_geomInfo.find(parameter) != m_geomInfo.end();
        }


        /**
         *
         */
        inline const std::string MeshGraph::GetGeomInfo(std::string parameter)
        {
            ASSERTL1(m_geomInfo.find(parameter) != m_geomInfo.end(),
                    "Parameter " + parameter + " does not exist.");
            return m_geomInfo[parameter];
        }


        /**
         *
         */
        inline int MeshGraph::GetNvertices() const
        {
            return int(m_vertSet.size());
        }


        /**
         *
         */
        inline VertexComponentSharedPtr MeshGraph::GetVertex(int id)
        {
            VertexComponentSharedPtr returnval;
            VertexMap::iterator x = m_vertSet.find(id);
            ASSERTL0(x != m_vertSet.end(),
                     "Vertex " + boost::lexical_cast<string>(id)
                     + " not found.");
            return x->second;
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<HexGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllHexGeoms();
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<PrismGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllPrismGeoms();
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<TetGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllTetGeoms();
        }


        /**
         *
         */
        template<>
        inline const std::map<int, boost::shared_ptr<PyrGeom> >& MeshGraph::GetAllElementsOfType() const
        {
            return GetAllPyrGeoms();
        }
    };
};

#endif

