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

namespace berrc = boost::system::errc;

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

        std::string GetElemGroupName(const int i)
        {
            std::stringstream nameSS;
            nameSS << "ELEMENTS" << i;
            return nameSS.str();
        }

        /**
         *
         */
        void FieldIOHdf5::v_Write(const std::string &outFile,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                const FieldMetaDataMap &fieldmetadatamap)
        {
            size_t nFields = fielddefs.size();
            // Check everything seems sensible
            ASSERTL1(fielddefs.size() == fielddata.size(),
                    "Length of fielddefs and fielddata incompatible");
            for (int f = 0; f < nFields; ++f)
            {
                ASSERTL1(fielddata[f].size() > 0,
                        "Fielddata vector must contain at least one value.");

                ASSERTL1(fielddata[f].size() ==
                        fielddefs[f]->m_fields.size() *
                        CheckFieldDefinition(fielddefs[f]),
                        "Invalid size of fielddata vector.");
            }

            // If the file exists already, rmtree it
            if (m_comm->RemoveExistingFiles())
            {
                if (m_comm->GetRank() == 0)
                {
                    try
                    {
                        fs::path specPath(outFile);
                        fs::remove_all(specPath);
                    } catch (fs::filesystem_error& e)
                    {
                        ASSERTL0(
                                e.code().value()
                                        == berrc::no_such_file_or_directory,
                                "Filesystem error: " + string(e.what()));
                    }
                }
            }

            // So first we're going to create the file structure on rank 0
            // because this involves a lot of metadata operations which are
            // collective and hence slow in parallel.

            // This does involve the other ranks because everyone needs to know
            // where in the global list they're writing their subset
            typedef Array<OneD, unsigned long long> SizeArray;
            SizeArray local_n_elem(nFields);
            SizeArray global_n_elem(nFields);
            SizeArray vals_per_elem(nFields);
            for (int f = 0; f < nFields; ++f)
            {
                // Get the global index at which we're writing our local data
                // and the global number of elements
                local_n_elem[f] = fielddefs[f]->m_elementIDs.size();
                vals_per_elem[f] = fielddata[f].size() / local_n_elem[f];
            }

            bool amRoot = (m_comm->GetRank() == 0);

            SizeArray all_n_elem = m_comm->Gather(0, local_n_elem);
            SizeArray all_idx(all_n_elem.num_elements());

            if (amRoot)
            {
                // Compute global_idx = sum(local_n_elem[rank i] for i in range(0, my rank))
                for (int f = 0; f < nFields; ++f)
                    all_idx[f] = 0;

                for (int p = 1; p < m_comm->GetSize(); ++p)
                {
                    for (int f = 0; f < nFields; ++f)
                    {
                        all_idx[p * nFields + f] =
                                all_idx[(p - 1) * nFields + f]
                                        + all_n_elem[(p - 1) * nFields + f];
                    }
                }
                for (int f = 0; f < nFields; ++f)
                {
                    int i = (m_comm->GetSize() - 1) + f;
                    global_n_elem[f] = all_idx[i] + all_n_elem[i];
                }
            }

            SizeArray global_idx = m_comm->Scatter(0, all_idx);
            ASSERTL0(global_idx.num_elements() == nFields,
                    "Didn't get the right number of field indices")

            // Only one rank touches the file to create the layout.
            if (amRoot)
            {
                H5::FileSharedPtr outfile = H5::File::Create(outFile,
                        H5F_ACC_TRUNC);
                H5::GroupSharedPtr root = outfile->CreateGroup("NEKTAR");
                TagWriterSharedPtr info_writer(new H5TagWriter(root));
                AddInfoTag(info_writer, fieldmetadatamap);

                std::vector < std::vector<unsigned long long> > decompositions(nFields);
                for (int f = 0; f < nFields; ++f)
                    decompositions[f].resize(m_comm->GetSize());
                for (int p = 0; p < m_comm->GetSize(); ++p)
                {
                    for (int f = 0; f < nFields; ++f)
                    {
                        decompositions[f][p] = all_n_elem[p * nFields + f];
                    }
                }

                for (int f = 0; f < nFields; ++f)
                {
                    //---------------------------------------------
                    // Write ELEMENTS

                    // first, give it a name that will be unique in the group
                    std::string elemGroupName = GetElemGroupName(f);

                    H5::GroupSharedPtr elemTag = root->CreateGroup(
                            elemGroupName);
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
                                numModesStringStream
                                        << fielddefs[f]->m_numModes[i];
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
                                numModesStringStream
                                        << fielddefs[f]->m_numModes[i];
                                first = false;
                            }
                        }

                        numModesString = numModesStringStream.str();
                    }
                    elemTag->SetAttribute("NUMMODESPERDIR", numModesString);

                    // DECOMPOSITION
                    elemTag->SetAttribute("DECOMPOSITION", decompositions[f]);

                    // Create empty datasets
                    // ID
                    H5::DataSpaceSharedPtr idSpace = H5::DataSpace::OneD(
                            global_n_elem[f]);
                    H5::DataTypeSharedPtr idType = H5::DataType::OfObject(
                            fielddefs[f]->m_elementIDs[0]);
                    H5::DataSetSharedPtr idDS = elemTag->CreateDataSet("ID",
                            idType, idSpace);

                    // DATA
                    H5::DataSpaceSharedPtr dataSpace = H5::DataSpace::OneD(
                            global_n_elem[f] * vals_per_elem[f]);
                    H5::DataTypeSharedPtr dataType = H5::DataType::OfObject(
                            fielddata[f][0]);
                    H5::DataSetSharedPtr dataDS = elemTag->CreateDataSet("DATA",
                            dataType, dataSpace);

                }
            }
            // RAII => file gets closed automatically

            // Set properties for parallel (if we're in parallel)
            H5::PListSharedPtr parallelProps = H5::PList::Default();
            H5::PListSharedPtr writePL = H5::PList::Default();
            if (m_comm->GetSize() > 1)
            {
                // Use MPI/O to access the file
                parallelProps = H5::PList::FileAccess();
                parallelProps->SetMpio(m_comm);
                // Use collective IO
                writePL = H5::PList::DatasetXfer();
                writePL->SetDxMpioCollective();
            }

            // Reopen the file
            H5::FileSharedPtr outfile = H5::File::Open(outFile, H5F_ACC_RDWR,
                    parallelProps);
            H5::GroupSharedPtr root = outfile->OpenGroup("NEKTAR");

            for (int f = 0; f < nFields; ++f)
            {
                std::string elemGroupName = GetElemGroupName(f);
                H5::GroupSharedPtr elemTag = root->OpenGroup(elemGroupName);
                // TODO: should probably add a check that the attributes are right to make sure we're OK to write

                H5::DataSetSharedPtr dset;
                H5::DataSpaceSharedPtr filespace;

                // Write the ID dataset to the file
                dset = elemTag->OpenDataSet("ID");
                // Select the section in the file
                filespace = dset->GetSpace();
                filespace->SelectRange(global_idx[f], local_n_elem[f]);
                dset->Write(fielddefs[f]->m_elementIDs, filespace, writePL);

                // Write the actual field data
                dset = elemTag->OpenDataSet("DATA");

                // Select the section in the file
                filespace = dset->GetSpace();
                filespace->SelectRange(global_idx[f] * vals_per_elem[f],
                        local_n_elem[f] * vals_per_elem[f]);

                dset->Write(fielddata[f], filespace, writePL);
            }
            // RAII => file gets closed automatically
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
