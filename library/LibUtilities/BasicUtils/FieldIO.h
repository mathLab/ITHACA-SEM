///////////////////////////////////////////////////////////////////////////////
//
// File FieldIO.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Field IO prototype definitions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIO_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_FIELDIO_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Points.h>
#include <boost/assign/list_of.hpp>
#include <tinyxml.h>

#include <LibUtilities/BasicUtils/NekFactory.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        /// Base class for writing hierarchical data (XML or HDF5)
        class TagWriter
        {
            public:
                /// Create a child node
                virtual boost::shared_ptr<TagWriter> AddChild(
                        const std::string& name) = 0;
                /// Set an attribute on the node
                virtual void SetAttr(const std::string& key,
                        const std::string& val) = 0;
            protected:
                virtual ~TagWriter();
        };
        typedef boost::shared_ptr<TagWriter> TagWriterSharedPtr;

        /// Simple class for writing to XML
        class XmlTagWriter : public TagWriter
        {
            public:
                XmlTagWriter(TiXmlElement* elem);
                TagWriterSharedPtr AddChild(const std::string& name);
                void SetAttr(const std::string& key, const std::string& val);
            private:
                TiXmlElement* m_El;
        };

        static std::vector<NekDouble> NullNekDoubleVector;
        static std::vector<LibUtilities::PointsType> NullPointsTypeVector;
        static std::vector<unsigned int> NullUnsignedIntVector;

        typedef std::map<std::string, std::string> FieldMetaDataMap;
        static FieldMetaDataMap NullFieldMetaDataMap;
        static std::vector<std::vector<NekDouble> > NullVectorNekDoubleVector =
                boost::assign::list_of(NullNekDoubleVector);

        struct FieldDefinitions
        {
                FieldDefinitions(ShapeType shapeType,
                        const std::vector<unsigned int> &elementIDs, // vector[2]
                        const std::vector<LibUtilities::BasisType> &basis,
                        bool uniOrder,
                        const std::vector<unsigned int> &numModes,
                        const std::vector<std::string> &fields, int NumHomoDir =
                                0, const std::vector<NekDouble> &HomoLengths =
                                NullNekDoubleVector,
                        bool homoStrips = false,
                        const std::vector<unsigned int> &HomoSIDs =
                                NullUnsignedIntVector,
                        const std::vector<unsigned int> &HomoZIDs =
                                NullUnsignedIntVector,
                        const std::vector<unsigned int> &HomoYIDs =
                                NullUnsignedIntVector,
                        const std::vector<LibUtilities::PointsType> &points =
                                NullPointsTypeVector, bool pointsDef = false,
                        const std::vector<unsigned int> &numPoints =
                                NullUnsignedIntVector,
                        bool numPointsDef = false) :
                        m_shapeType(shapeType), m_elementIDs(elementIDs), m_basis(
                                basis), m_numHomogeneousDir(NumHomoDir), m_homogeneousLengths(
                                HomoLengths), m_homoStrips(homoStrips), m_homogeneousSIDs(HomoSIDs),
                                m_homogeneousZIDs(HomoZIDs), m_homogeneousYIDs(
                                HomoYIDs), m_points(points), m_pointsDef(
                                pointsDef), m_uniOrder(uniOrder), m_numModes(
                                numModes), m_numPoints(numPoints), m_numPointsDef(
                                numPointsDef), m_fields(fields)
                {
                }

                ShapeType m_shapeType;
                std::vector<unsigned int> m_elementIDs;
                std::vector<LibUtilities::BasisType> m_basis;
                int m_numHomogeneousDir;
                std::vector<NekDouble> m_homogeneousLengths;
                bool m_homoStrips;
                std::vector<unsigned int> m_homogeneousSIDs;
                std::vector<unsigned int> m_homogeneousZIDs;
                std::vector<unsigned int> m_homogeneousYIDs;

                /// Define the type of points per direction.
                std::vector<LibUtilities::PointsType> m_points;
                bool m_pointsDef;
                /// Define order of the element group.
                /// * UniOrder: same order for each element
                /// * MixOrder: definition of a different order for each element.
                bool m_uniOrder;
                /// Define number of modes per direction.
                std::vector<unsigned int> m_numModes;
                std::vector<unsigned int> m_numPoints;
                bool m_numPointsDef;
                std::vector<std::string> m_fields;
        };

        typedef boost::shared_ptr<FieldDefinitions> FieldDefinitionsSharedPtr;

        class DataSource
        {
        };
        typedef boost::shared_ptr<DataSource> DataSourceSharedPtr;

        /// Write a field file in serial only
        LIB_UTILITIES_EXPORT void Write(const std::string &outFile,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap);

        /// Imports an FLD file
        LIB_UTILITIES_EXPORT void Import(const std::string& infilename,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata =
                        NullVectorNekDoubleVector,
                FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
                const Array<OneD, int> ElementiDs = NullInt1DArray);

        // Forward declare
        class FieldIO;

        /// Datatype of the NekFactory used to instantiate classes
        typedef LibUtilities::NekFactory<std::string, FieldIO,
            LibUtilities::CommSharedPtr, bool> FieldIOFactory;

        LIB_UTILITIES_EXPORT FieldIOFactory& GetFieldIOFactory();

        /// Class for operating on FLD files
        class FieldIO : public boost::enable_shared_from_this<FieldIO>
        {
            public:

                /// Constructor
                LIB_UTILITIES_EXPORT
                FieldIO(LibUtilities::CommSharedPtr pComm,
                        bool sharedFilesystem);

                /// Write data in FLD format
                LIB_UTILITIES_EXPORT
                inline void Write(const std::string &outFile,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata,
                        const FieldMetaDataMap &fieldinfomap =
                                NullFieldMetaDataMap);

                /// Imports an FLD file.
                LIB_UTILITIES_EXPORT
                inline void Import(const std::string& infilename,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata =
                                NullVectorNekDoubleVector,
                        FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
                        const Array<OneD, int> ElementiDs = NullInt1DArray);

                /// Imports the definition of the meta data
                LIB_UTILITIES_EXPORT
                DataSourceSharedPtr ImportFieldMetaData(std::string filename,
                        FieldMetaDataMap &fieldmetadatamap);

                /// Imports the definition of the fields.
                LIB_UTILITIES_EXPORT
                virtual void ImportFieldDefs(DataSourceSharedPtr dataSource,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        bool expChild) = 0;

                // Figure out what type of FLD file we have.
                // Collective on comm.
                static const std::string GetFileType(
                        const std::string& filename, CommSharedPtr comm);

                virtual const std::string& GetClassName() const = 0;
            protected:
                /// Communicator to use when writing parallel format
                LibUtilities::CommSharedPtr m_comm;
                bool m_sharedFilesystem;

                LIB_UTILITIES_EXPORT
                void AddInfoTag(TiXmlElement * root,
                        const FieldMetaDataMap &fieldmetadatamap);

                LIB_UTILITIES_EXPORT
                void AddInfoTag(TagWriterSharedPtr root,
                        const FieldMetaDataMap &fieldmetadatamap);

                LIB_UTILITIES_EXPORT
                void GenerateSeqString(const std::vector<unsigned int> &elmtids,
                        std::string &idString);

                LIB_UTILITIES_EXPORT
                int CheckFieldDefinition(
                        const FieldDefinitionsSharedPtr &fielddefs);

                LIB_UTILITIES_EXPORT virtual std::string GetFileEnding() const
                {
                    return "fld";
                }

                LIB_UTILITIES_EXPORT
                virtual void v_Write(const std::string &outFile,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata,
                        const FieldMetaDataMap &fieldinfomap) = 0;

                /// Imports an FLD file.
                LIB_UTILITIES_EXPORT
                virtual void v_Import(const std::string& infilename,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata =
                                NullVectorNekDoubleVector,
                        FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
                        const Array<OneD, int> ElementiDs = NullInt1DArray) = 0;

                LIB_UTILITIES_EXPORT
                virtual DataSourceSharedPtr v_ImportFieldMetaData(
                        std::string filename,
                        FieldMetaDataMap &fieldmetadatamap) = 0;

        };

        typedef boost::shared_ptr<FieldIO> FieldIOSharedPtr;

        inline FieldIOSharedPtr MakeDefaultFieldIO(
                const LibUtilities::SessionReaderSharedPtr session)
        {
            std::string iofmt("Xml");
            if (session->DefinesSolverInfo("FieldIO_Format"))
            {
                iofmt = session->GetSolverInfo("FieldIO_Format");
            }
            return GetFieldIOFactory().CreateInstance(
                iofmt, session->GetComm(),
                session->DefinesCmdLineArgument("shared-filesystem"));
        }
        // Collective on session's communicator
        FieldIOSharedPtr MakeFieldIOForFile(
                const LibUtilities::SessionReaderSharedPtr session,
                const std::string& filename);

        inline void FieldIO::Write(const std::string &outFile,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                const FieldMetaDataMap &fieldinfomap)
        {
            v_Write(outFile, fielddefs, fielddata, fieldinfomap);
        }

        inline void FieldIO::Import(
            const std::string& infilename,
            std::vector<FieldDefinitionsSharedPtr> &fielddefs,
            std::vector<std::vector<NekDouble> > &fielddata,
            FieldMetaDataMap &fieldinfomap,
            const Array<OneD, int> ElementiDs)
        {
            v_Import(infilename, fielddefs, fielddata, fieldinfomap, ElementiDs);
        }

        inline DataSourceSharedPtr FieldIO::ImportFieldMetaData(
                std::string filename, FieldMetaDataMap &fieldmetadatamap)
        {
            return v_ImportFieldMetaData(filename, fieldmetadatamap);
        }
    }
}
#endif
