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

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Points.h>
#include <tinyxml/tinyxml.h>

// These are required for the Write(...) and Import(...) functions.
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        
        static std::vector<NekDouble> NullNekDoubleVector;
        static std::vector<LibUtilities::PointsType> NullPointsTypeVector;
        static std::vector<unsigned int> NullUnsignedIntVector;

        typedef std::map<std::string, NekDouble>  FieldMetaDataMap;
        static  FieldMetaDataMap  NullFieldMetaDataMap;
        

        struct FieldDefinitions
        {
        FieldDefinitions(ShapeType shapeType,
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
            
            ShapeType	                          m_shapeType;
            std::vector<unsigned int>		  m_elementIDs;
            std::vector<LibUtilities::BasisType>  m_basis;
            int					  m_numHomogeneousDir;
            std::vector<NekDouble>		  m_homogeneousLengths;
            std::vector<unsigned int>		  m_homogeneousZIDs;
            std::vector<unsigned int>		  m_homogeneousYIDs;
            
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


            /* --- FLD handling routines ---- */
        LIB_UTILITIES_EXPORT void Write(
                        const std::string &outFile,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> >   &fielddata,
                        FieldMetaDataMap &fieldinfomap  = NullFieldMetaDataMap);

        /// Imports an FLD file.
        LIB_UTILITIES_EXPORT void Import(
                        const std::string& infilename,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata,
                        FieldMetaDataMap &fieldinfomap  = NullFieldMetaDataMap);

        /// Imports the definition of the meta data 
        LIB_UTILITIES_EXPORT void ImportFieldMetaData(std::string filename,
                                               FieldMetaDataMap &fieldmetadatamap);

        /// Imports the definition of the meta data 
        LIB_UTILITIES_EXPORT void ImportFieldMetaData(
                        TiXmlDocument &doc,
                        FieldMetaDataMap &fieldmetadatamap);

        /// Imports the definition of the fields.
        LIB_UTILITIES_EXPORT void ImportFieldDefs(
                        TiXmlDocument &doc,
                        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        bool expChild);

        /// Imports the data fileds.
        LIB_UTILITIES_EXPORT void ImportFieldData(
                        TiXmlDocument &doc,
                        const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                        std::vector<std::vector<NekDouble> > &fielddata);

        LIB_UTILITIES_EXPORT int CheckFieldDefinition( 
                        const FieldDefinitionsSharedPtr  &fielddefs);
    }
}
#endif
