////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIOHdf5.cpp
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
//  Description: I/O routines relating to Fields into HDF
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/FieldIOHdf5.h>
#include <LibUtilities/BasicUtils/H5.h>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif

namespace Nektar
{
    namespace LibUtilities
    {
        namespace H5
        {
            template<>
            inline DataTypeSharedPtr DataTypeTraits<BasisType>::GetType()
            {
                return PredefinedDataType::Native<int>();
            }
        }
        H5TagWriter::H5TagWriter(H5::GroupSharedPtr grp) :
                m_Group(grp)
        {
        }
        TagWriterSharedPtr H5TagWriter::AddChild(const std::string& name)
        {
            H5::GroupSharedPtr child = m_Group->CreateGroup(name);
            return TagWriterSharedPtr(new H5TagWriter(child));
        }
        ;

        void H5TagWriter::SetAttr(const std::string& key,
                const std::string& val)
        {
            m_Group->SetAttribute(key, val);
        }
        class H5DataSource : public DataSource
        {
                H5::FileSharedPtr doc;
            public:
                H5DataSource(const std::string& fn) :
                        doc(H5::File::Open(fn, H5F_ACC_RDONLY))
                {
                }

                H5::FileSharedPtr Get()
                {
                    return doc;
                }
                const H5::FileSharedPtr Get() const
                {
                    return doc;
                }

                static DataSourceSharedPtr create(const std::string& fn)
                {
                    return DataSourceSharedPtr(new H5DataSource(fn));
                }
        };
        typedef boost::shared_ptr<H5DataSource> H5DataSourceSharedPtr;

        std::string FieldIOHdf5::className =
                GetFieldIOFactory().RegisterCreatorFunction("Hdf5",
                        FieldIOHdf5::create,
                        "HDF5-based output of field data.");

        FieldIOHdf5::FieldIOHdf5(LibUtilities::CommSharedPtr pComm) :
                FieldIO(pComm)
        {
        }
        /**
         *
         */
        void FieldIOHdf5::v_Write(const std::string &outFile,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            // Check everything seems sensible
            ASSERTL1(fielddefs.size() == fielddata.size(),
                    "Length of fielddefs and fielddata incompatible");
            for (int f = 0; f < fielddefs.size(); ++f)
            {
                ASSERTL1(fielddata[f].size() > 0,
                        "Fielddata vector must contain at least one value.");

                ASSERTL1(fielddata[f].size() ==
                        fielddefs[f]->m_fields.size() *
                        CheckFieldDefinition(fielddefs[f]),
                        "Invalid size of fielddata vector.");
            }

            // Prepare to write out data. In parallel, we must create directory
            // and determine the full pathname to the file to write out.
            // Any existing file/directory which is in the way is removed.
            std::string filename = SetUpOutput(outFile, fielddefs,
                    fieldmetadatamap);

            // Create the file (partition)
            H5::FileSharedPtr outfile = H5::File::Create(filename,
                    H5F_ACC_TRUNC);

            H5::GroupSharedPtr root = outfile->CreateGroup("NEKTAR");
            TagWriterSharedPtr info_writer(new H5TagWriter(root));
            AddInfoTag(info_writer, fieldmetadatamap);

            for (int f = 0; f < fielddefs.size(); ++f)
            {
                //---------------------------------------------
                // Write ELEMENTS

                // first, give it a name that will be unique in the group
                std::string elemGroupName;
                {
                    std::stringstream nameSS;
                    nameSS << "ELEMENTS" << f;
                    elemGroupName = nameSS.str();
                }

                H5::GroupSharedPtr elemTag = root->CreateGroup(elemGroupName);
                // Write FIELDS
                // i.e. m_fields is a vector<string> (variable length)
                elemTag->SetAttribute("FIELDS", fielddefs[f]->m_fields);

                // Write SHAPE
                // because this isn't just a vector of data, serialise it ourselves for now
                std::string shapeString;
                {
                    std::stringstream shapeStringStream;
                    shapeStringStream
                            << ShapeTypeMap[fielddefs[f]->m_shapeType];
                    if (fielddefs[f]->m_numHomogeneousDir == 1)
                    {
                        shapeStringStream << "-HomogenousExp1D";
                    }
                    else if (fielddefs[f]->m_numHomogeneousDir == 2)
                    {
                        shapeStringStream << "-HomogenousExp2D";
                    }

                    shapeString = shapeStringStream.str();
                }
                elemTag->SetAttribute("SHAPE", shapeString);

                // Write BASIS
                elemTag->SetAttribute("BASIS", fielddefs[f]->m_basis);

                // Write homogeneuous length details
                if (fielddefs[f]->m_numHomogeneousDir)
                {
                    elemTag->SetAttribute("HOMOGENEOUSLENGTHS",
                            fielddefs[f]->m_homogeneousLengths);
                }

                // Write homogeneuous planes/lines details
                if (fielddefs[f]->m_numHomogeneousDir)
                {
                    if (fielddefs[f]->m_homogeneousYIDs.size() > 0)
                    {
                        elemTag->SetAttribute("HOMOGENEOUSYIDS",
                                fielddefs[f]->m_homogeneousYIDs);
                    }

                    if (fielddefs[f]->m_homogeneousZIDs.size() > 0)
                    {
                        elemTag->SetAttribute("HOMOGENEOUSZIDS",
                                fielddefs[f]->m_homogeneousZIDs);
                    }
                }

                // Write NUMMODESPERDIR
                // because this isn't just a vector of data, serialise it by hand for now
                std::string numModesString;
                {
                    std::stringstream numModesStringStream;

                    if (fielddefs[f]->m_uniOrder)
                    {
                        numModesStringStream << "UNIORDER:";
                        // Just dump single definition
                        bool first = true;
                        for (std::vector<int>::size_type i = 0;
                                i < fielddefs[f]->m_basis.size(); i++)
                        {
                            if (!first)
                                numModesStringStream << ",";
                            numModesStringStream << fielddefs[f]->m_numModes[i];
                            first = false;
                        }
                    }
                    else
                    {
                        numModesStringStream << "MIXORDER:";
                        bool first = true;
                        for (std::vector<int>::size_type i = 0;
                                i < fielddefs[f]->m_numModes.size(); i++)
                        {
                            if (!first)
                                numModesStringStream << ",";
                            numModesStringStream << fielddefs[f]->m_numModes[i];
                            first = false;
                        }
                    }

                    numModesString = numModesStringStream.str();
                }
                elemTag->SetAttribute("NUMMODESPERDIR", numModesString);

//                // Configure for compression.
//                // This requires chunking - arbitrarily pick 1024 elements
//                H5::PListSharedPtr zip = H5::PList::DatasetCreate();
//                std::vector < hsize_t > chunk;
//                chunk.push_back(std::min(size_t(1024), fielddata[f].size()));
//                zip->SetChunk(chunk);
//                // Higher compression levels will reduce filesize at the cost
//                // increasing time.
//                zip->SetDeflate(1);

                // Write the actual field data
                elemTag->CreateWriteDataSet("DATA", fielddata[f]);

                // Write ID
                elemTag->CreateWriteDataSet("ID", fielddefs[f]->m_elementIDs);
            }
            // Destruction of the H5::File and Group objects closes them. Yay RAII
        }

        void FieldIOHdf5::v_ImportFile(const std::string& infilename,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                DataSourceSharedPtr dataSource)
        {
            if (!dataSource)
                dataSource = H5DataSource::create(infilename);

            ImportFieldDefs(dataSource, fielddefs, false);
            if (fielddata != NullVectorNekDoubleVector)
            {
                ImportFieldData(dataSource, fielddefs, fielddata);
            }
        }
        void FieldIOHdf5::ImportFieldDefs(DataSourceSharedPtr dataSource,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                bool expChild)
        {
            H5DataSourceSharedPtr h5 = boost::static_pointer_cast < H5DataSource
                    > (dataSource);
            ASSERTL0(expChild == false,
                    "FieldDefs in <EXPANSIONS> tag not implemented");
            // Master tag within which all data is contained.
            H5::GroupSharedPtr master = h5->Get()->OpenGroup("NEKTAR");

            // The XML reader says here "Loop through all nektar tags, finding all of the element tags."
            // I believe that there is only allowed to be one NEKTAR tag?
            // If multiple are allowed, reinstate the loop.
            H5::Group::LinkIterator keyIt = master->begin(), keyEnd =
                    master->end();
            for (; keyIt != keyEnd; ++keyIt)
            {
                const std::string& grpName = *keyIt;

                const std::string prefix("ELEMENTS");
                if (!std::equal(prefix.begin(), prefix.end(), grpName.begin()))
                    continue;

                H5::GroupSharedPtr element = master->OpenGroup(grpName);
                // Extract the attributes.
                std::vector < std::string > Fields;
                // Shape is encoded by hand
                std::string shapeString;
                std::vector<BasisType> basis;
                std::vector<NekDouble> homoLengths;
                // Homogeneous points IDs
                std::vector<unsigned int> homoZIDs;
                std::vector<unsigned int> homoYIDs;
                // encoded by hand
                std::string numModesString;
                std::string numPointsString;
                std::string pointsString;
                bool pointDef = false;
                bool numPointDef = false;

                H5::Group::AttrIterator attrIt = element->attr_begin(),
                        attrEnd = element->attr_end();
                for (; attrIt != attrEnd; ++attrIt)
                {
                    const std::string& attrName = *attrIt;
                    if (attrName == "FIELDS")
                    {
                        element->GetAttribute(attrName, Fields);
                    }
                    else if (attrName == "SHAPE")
                    {
                        element->GetAttribute(attrName, shapeString);
                    }
                    else if (attrName == "BASIS")
                    {
                        element->GetAttribute(attrName, basis);
                    }
                    else if (attrName == "HOMOGENEOUSLENGTHS")
                    {
                        element->GetAttribute(attrName, homoLengths);
                    }
                    else if (attrName == "HOMOGENEOUSZIDS")
                    {
                        element->GetAttribute(attrName, homoZIDs);
                    }
                    else if (attrName == "HOMOGENEOUSYIDS")
                    {
                        element->GetAttribute(attrName, homoYIDs);
                    }
                    else if (attrName == "NUMMODESPERDIR")
                    {
                        element->GetAttribute(attrName, numModesString);
                    }
                    else if (attrName == "POINTSTYPE")
                    {
                        element->GetAttribute(attrName, pointsString);
                        pointDef = true;
                    }
                    else if (attrName == "NUMPOINTSPERDIR")
                    {
                        element->GetAttribute(attrName, numPointsString);
                        numPointDef = true;
                    }
                    else
                    {
                        std::string errstr("Unknown attribute: ");
                        errstr += attrName;
                        ASSERTL1(false, errstr.c_str());
                    }
                }

                // Check to see if homogeneous expansion and if so
                // strip down the shapeString definition
                int numHomoDir = 0;
                size_t loc;
                //---> This finds the first location of  'n'!
                if ((loc = shapeString.find_first_of("-")) != string::npos)
                {
                    if (shapeString.find("Exp1D") != string::npos)
                    {
                        numHomoDir = 1;
                    }
                    else // HomogeneousExp1D
                    {
                        numHomoDir = 2;
                    }

                    shapeString.erase(loc, shapeString.length());
                }

                // Reconstruct the fielddefs.
                // Get the geometrical shape
                ShapeType shape;
                bool valid = false;
                for (unsigned int j = 0; j < SIZE_ShapeType; j++)
                {
                    if (ShapeTypeMap[j] == shapeString)
                    {
                        shape = (ShapeType) j;
                        valid = true;
                        break;
                    }
                }

                ASSERTL0(valid,
                        std::string(
                                "Unable to correctly parse the shape type: ").append(
                                shapeString).c_str());

                // Check the basis is in range
                for (std::vector<BasisType>::const_iterator bIt = basis.begin(),
                        bEnd = basis.end(); bIt != bEnd; ++bIt)
                {
                    BasisType bt = *bIt;
                    ASSERTL0(bt >= 0 && bt < SIZE_BasisType,
                            "Unable to correctly parse the basis types.")
                }

                // Get points type
                std::vector<PointsType> points;

                if (pointDef)
                {
                    std::vector < std::string > pointsStrings;
                    valid = ParseUtils::GenerateOrderedStringVector(
                            pointsString.c_str(), pointsStrings);
                    ASSERTL0(valid,
                            "Unable to correctly parse the points types.");
                    for (std::vector<std::string>::size_type i = 0;
                            i < pointsStrings.size(); i++)
                    {
                        valid = false;
                        for (unsigned int j = 0; j < SIZE_PointsType; j++)
                        {
                            if (kPointsTypeStr[j] == pointsStrings[i])
                            {
                                points.push_back((PointsType) j);
                                valid = true;
                                break;
                            }
                        }

                        ASSERTL0(valid,
                                std::string(
                                        "Unable to correctly parse the points type: ").append(
                                        pointsStrings[i]).c_str());
                    }
                }

                // Get numModes
                std::vector<unsigned int> numModes;
                bool UniOrder = false;

                if (strstr(numModesString.c_str(), "UNIORDER:"))
                {
                    UniOrder = true;
                }

                valid = ParseUtils::GenerateOrderedVector(
                        numModesString.c_str() + 9, numModes);

                ASSERTL0(valid, "Unable to correctly parse the number of modes");

                // Get numPoints
                std::vector<unsigned int> numPoints;
                if (numPointDef)
                {
                    valid = ParseUtils::GenerateOrderedVector(
                            numPointsString.c_str(), numPoints);
                    ASSERTL0(valid,
                            "Unable to correctly parse the number of points.");
                }

                std::vector<unsigned int> elementIds;
                H5::DataSetSharedPtr idData = element->OpenDataSet("ID");
                idData->Read(elementIds);

                FieldDefinitionsSharedPtr fielddef = MemoryManager<
                        FieldDefinitions>::AllocateSharedPtr(shape, elementIds,
                        basis, UniOrder, numModes, Fields, numHomoDir,
                        homoLengths, homoZIDs, homoYIDs, points, pointDef,
                        numPoints, numPointDef);

                fielddefs.push_back(fielddef);
            }

        }

        void FieldIOHdf5::ImportFieldData(DataSourceSharedPtr dataSource,
                const std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata)
        {
            int cntdumps = 0;
            H5DataSourceSharedPtr h5 = boost::static_pointer_cast < H5DataSource
                    > (dataSource);
            // Master tag within which all data is contained.
            H5::GroupSharedPtr master = h5->Get()->OpenGroup("NEKTAR");

            // The XML reader says here "Loop through all nektar tags, finding all of the element tags."
            // I believe that there is only allowed to be one NEKTAR tag?
            // If multiple are allowed, reinstate the loop.
            H5::Group::LinkIterator keyIt = master->begin(), keyEnd =
                    master->end();
            for (; keyIt != keyEnd; ++keyIt)
            {
                const std::string& grpName = *keyIt;

                const std::string prefix("ELEMENTS");
                if (!std::equal(prefix.begin(), prefix.end(), grpName.begin()))
                    continue;

                H5::GroupSharedPtr element = master->OpenGroup(grpName);

                std::vector<NekDouble> elementFieldData;
                H5::DataSetSharedPtr dataset = element->OpenDataSet("DATA");
                dataset->Read(elementFieldData);
                fielddata.push_back(elementFieldData);
                int datasize = CheckFieldDefinition(fielddefs[cntdumps]);
                ASSERTL0(
                        fielddata[cntdumps].size()
                                == datasize
                                        * fielddefs[cntdumps]->m_fields.size(),
                        "Input data is not the same length as header infoarmation");

                cntdumps++;
            }
        }

        DataSourceSharedPtr FieldIOHdf5::v_ImportFieldMetaData(
                std::string filename, FieldMetaDataMap &fieldmetadatamap)
        {
            DataSourceSharedPtr ans = H5DataSource::create(filename);
            v_ImportFieldMetaData(ans, fieldmetadatamap);
            return ans;
        }

        void FieldIOHdf5::v_ImportFieldMetaData(DataSourceSharedPtr dataSource,
                FieldMetaDataMap &fieldmetadatamap)
        {
            H5DataSourceSharedPtr hdf = boost::static_pointer_cast
                    < H5DataSource > (dataSource);

            H5::GroupSharedPtr master = hdf->Get()->OpenGroup("NEKTAR");
            // New metadata format only in HDF
            H5::GroupSharedPtr metadata = master->OpenGroup("Metadata");

            if (metadata)
            {
                H5::Group::AttrIterator param = metadata->attr_begin(), pEnd =
                        metadata->attr_end();
                for (; param != pEnd; ++param)
                {
                    std::string paramString = *param;
                    if (paramString != "Provenance")
                    {
                        std::string paramBodyStr;
                        metadata->GetAttribute(paramString, paramBodyStr);
                        fieldmetadatamap[paramString] = paramBodyStr;
                    }
                }
            }
        }

    }
}
