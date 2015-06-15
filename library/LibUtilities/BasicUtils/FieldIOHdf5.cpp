////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldIO.cpp
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
            H5::File outfile(filename, H5F_ACC_TRUNC);

            H5::GroupSharedPtr root = outfile.CreateGroup("NEKTAR");
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

                // Configure for compression.
                // This requires chunking - arbitrarily pick 1024 elements
                H5::PListSharedPtr zip = H5::PList::DatasetCreate();
                std::vector<hsize_t> chunk;
                chunk.push_back(std::min(size_t(1024), fielddata[f].size()));
                zip->SetChunk(chunk);
                // Higher compression levels will reduce filesize at the cost
                // increasing time.
                zip->SetDeflate(1);

                // Write the actual field data
                elemTag->CreateWriteDataSet("DATA", fielddata[f], zip);

                // Write ID
                elemTag->CreateWriteDataSet("ID", fielddefs[f]->m_elementIDs);
            }
            // Destruction of the H5::File and Group objects closes them. Yay RAII
        }

        void FieldIOHdf5::v_Import(const std::string& infilename,
                std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                std::vector<std::vector<NekDouble> > &fielddata,
                FieldMetaDataMap &fieldinfomap,
                const Array<OneD, int> ElementiDs)
        {
            ASSERTL1(0, "Import of HDF not implemented.");
        }

        void FieldIOHdf5::v_ImportFieldMetaData(std::string filename,
                FieldMetaDataMap &fieldmetadatamap)
        {
            ASSERTL1(0, "Import of HDF metadata not implemented.");
        }
    }
}
