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
#include <LibUtilities/BasicUtils/ShapeType.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Points.h>
#include <tinyxml.h>

#include <LibUtilities/BasicUtils/NekFactory.hpp>

namespace Nektar
{
namespace LibUtilities
{

typedef std::map<std::string, std::string> FieldMetaDataMap;
static FieldMetaDataMap NullFieldMetaDataMap;

/**
 * @brief Base class for writing hierarchical data (XML or HDF5).
 */
class TagWriter
{
public:
    /// Create a child node.
    virtual std::shared_ptr<TagWriter> AddChild(const std::string &name) = 0;
    /// Set an attribute on the node.
    virtual void SetAttr(const std::string &key, const std::string &val) = 0;

protected:
    virtual ~TagWriter() {}
};
typedef std::shared_ptr<TagWriter> TagWriterSharedPtr;

/**
 * @class A simple class encapsulating a data source. This allows us to pass
 * around native file formats in virtual functions without resorting to using
 * the filename.
 */
class DataSource
{
};
typedef std::shared_ptr<DataSource> DataSourceSharedPtr;

/**
 * @brief Metadata that describes the storage properties of field output.
 *
 * The purpose of this struct is to describe the format of binary field data.
 * This can then be used in the library to determine appropriate actions. For
 * example, when restarting a simulation, the information this struct
 * encapsulates can be used to determine whether interpolation is required to a
 * different polynomial order depending on the order of the simulation versus
 * the order of the restart file.
 *
 * We note that some of the parameters here include:
 *
 * - Element shape type and the basis used
 * - Element IDs, which determines the order of data written in a field
 * - The field names (e.g. u for x-velocity) so that multi-field storage can be
 *   agglomerated into one data block. Each field is written in the order
 *   specified here.
 * - Number of modes, including support for variable polynomial order
 * - Homogeneous information (dimension of homogeneity, IDs of planes and/or
 *   strips if they are used)
 */
struct FieldDefinitions
{
    /// Default constructor
    FieldDefinitions() : m_shapeType(eNoShapeType),
                         m_numHomogeneousDir(0),
                         m_homoStrips(false),
                         m_pointsDef(false),
                         m_uniOrder(true),
                         m_numPointsDef(false)
    {
    }

    /// Simple constructor to allocate all internal properties.
    FieldDefinitions(
        ShapeType shapeType,
        const std::vector<unsigned int> &elementIDs, // vector[2]
        const std::vector<LibUtilities::BasisType> &basis,
        bool uniOrder,
        const std::vector<unsigned int> &numModes,
        const std::vector<std::string> &fields,
        int NumHomoDir                            = 0,
        const std::vector<NekDouble> &HomoLengths = NullNekDoubleVector,
        bool homoStrips                           = false,
        const std::vector<unsigned int> &HomoSIDs = NullUnsignedIntVector,
        const std::vector<unsigned int> &HomoZIDs = NullUnsignedIntVector,
        const std::vector<unsigned int> &HomoYIDs = NullUnsignedIntVector,
        const std::vector<LibUtilities::PointsType> &points =
            NullPointsTypeVector,
        bool pointsDef                             = false,
        const std::vector<unsigned int> &numPoints = NullUnsignedIntVector,
        bool numPointsDef = false)
        : m_shapeType(shapeType), m_elementIDs(elementIDs), m_basis(basis),
          m_numHomogeneousDir(NumHomoDir), m_homogeneousLengths(HomoLengths),
          m_homoStrips(homoStrips), m_homogeneousSIDs(HomoSIDs),
          m_homogeneousZIDs(HomoZIDs), m_homogeneousYIDs(HomoYIDs),
          m_points(points), m_pointsDef(pointsDef), m_uniOrder(uniOrder),
          m_numModes(numModes), m_numPoints(numPoints),
          m_numPointsDef(numPointsDef), m_fields(fields)
    {
    }

    /// Shape type of this field data.
    ShapeType                             m_shapeType;
    /// Element IDs of the field data.
    std::vector<unsigned int>             m_elementIDs;
    /// Vector of basis types for each of the coordinate directions.
    std::vector<LibUtilities::BasisType>  m_basis;
    /// Number of homogeneous directions, in the range \f$ 0\leq d \leq 3 \f$.
    int                                   m_numHomogeneousDir;
    /// Spatial lengths of each homogeneous direction.
    std::vector<NekDouble>                m_homogeneousLengths;
    /// Boolean determining whether homogeneous strips are used.
    bool                                  m_homoStrips;
    /// IDs corresponding to homogeneous strip IDs.
    std::vector<unsigned int>             m_homogeneousSIDs;
    /// IDs corresponding to z-direction homogeneous IDs.
    std::vector<unsigned int>             m_homogeneousZIDs;
    /// IDs corresponding to y-direction homogeneous IDs.
    std::vector<unsigned int>             m_homogeneousYIDs;
    /// Define the type of points per direction.
    std::vector<LibUtilities::PointsType> m_points;
    /// Boolean determining whether points have been defined in output.
    bool                                  m_pointsDef;
    /// Define order of the element group.
    /// * UniOrder: same order for each element
    /// * MixOrder: definition of a different order for each element.
    bool                                  m_uniOrder;
    /// Define number of modes per direction.
    std::vector<unsigned int>             m_numModes;
    /// Define number of points per direction.
    std::vector<unsigned int>             m_numPoints;
    /// Boolean determining whether number of points has been defined.
    bool                                  m_numPointsDef;
    /// Vector of field names that this data encapsulates.
    std::vector<std::string>              m_fields;
};

typedef std::shared_ptr<FieldDefinitions> FieldDefinitionsSharedPtr;

LIB_UTILITIES_EXPORT void Write(
    const std::string &outFile,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata,
    const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
    const bool backup = false);
LIB_UTILITIES_EXPORT void Import(
    const std::string &infilename,
    std::vector<FieldDefinitionsSharedPtr> &fielddefs,
    std::vector<std::vector<NekDouble> > &fielddata = NullVectorNekDoubleVector,
    FieldMetaDataMap &fieldinfomap                  = NullFieldMetaDataMap,
    const Array<OneD, int> &ElementIDs = NullInt1DArray);

// Forward declare
class FieldIO;

/// Datatype of the NekFactory used to instantiate classes
typedef LibUtilities::NekFactory<std::string,
                                 FieldIO,
                                 LibUtilities::CommSharedPtr,
                                 bool> FieldIOFactory;

LIB_UTILITIES_EXPORT FieldIOFactory &GetFieldIOFactory();

/**
 * @brief Class for operating on Nektar++ input/output files.
 *
 * Nektar++ input/output of field data can be described as follows:
 *
 *   - The FieldDefinitions class defines the metadata that is associated with
 *     one or more scalar field variables, and determines the storage length.
 *   - Each scalar field is stored in a single contiguous vector which is
 *     written/read by the functions defined in this class.
 *   - Optional metadata can be read/written in a simple key/value pair map
 *     FieldMetaDataMap. This can be used to define, e.g. the timestep of the
 *     simulation.
 *
 * This base class represents the minimum functionality that subclasses need to
 * implement in order to implement the above functionality. Each subclass is
 * free to determine its own file structure and parallel behaviour.
 */
class FieldIO : public std::enable_shared_from_this<FieldIO>
{
public:
    LIB_UTILITIES_EXPORT FieldIO(
        LibUtilities::CommSharedPtr pComm, bool sharedFilesystem);

    LIB_UTILITIES_EXPORT virtual ~FieldIO()
    {
    }

    LIB_UTILITIES_EXPORT inline void Write(
        const std::string &outFile,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata,
        const FieldMetaDataMap &fieldinfomap = NullFieldMetaDataMap,
        const bool backup = false);

    LIB_UTILITIES_EXPORT inline void Import(
        const std::string &infilename,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> > &fielddata =
            NullVectorNekDoubleVector,
        FieldMetaDataMap &fieldinfomap    = NullFieldMetaDataMap,
        const Array<OneD, int> &ElementIDs = NullInt1DArray);

    LIB_UTILITIES_EXPORT DataSourceSharedPtr ImportFieldMetaData(
        const std::string &filename,
        FieldMetaDataMap  &fieldmetadatamap);

    LIB_UTILITIES_EXPORT static const std::string GetFileType(
        const std::string &filename, CommSharedPtr comm);
    LIB_UTILITIES_EXPORT virtual const std::string &GetClassName() const = 0;

    LIB_UTILITIES_EXPORT static std::shared_ptr<FieldIO> CreateDefault(
        const LibUtilities::SessionReaderSharedPtr session);
    LIB_UTILITIES_EXPORT static std::shared_ptr<FieldIO> CreateForFile(
        const LibUtilities::SessionReaderSharedPtr session,
        const std::string &filename);
    LIB_UTILITIES_EXPORT static void AddInfoTag(
        TagWriterSharedPtr      root,
        const FieldMetaDataMap &fieldmetadatamap);

protected:
    /// Communicator to use when writing parallel format
    LibUtilities::CommSharedPtr m_comm;
    /// Boolean dictating whether we are on a shared filesystem.
    bool                        m_sharedFilesystem;

    LIB_UTILITIES_EXPORT int CheckFieldDefinition(
        const FieldDefinitionsSharedPtr &fielddefs);

    /**
     * @brief Helper function that determines default file extension.
     */
    LIB_UTILITIES_EXPORT virtual std::string GetFileEnding() const
    {
        return "fld";
    }

    LIB_UTILITIES_EXPORT std::string SetUpOutput(
        const std::string outname, bool perRank, bool backup = false);

    /// @copydoc FieldIO::Write
    LIB_UTILITIES_EXPORT virtual void v_Write(
        const std::string                      &outFile,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> >   &fielddata,
        const FieldMetaDataMap                 &fieldinfomap,
        const bool                              backup = false) = 0;

    /// @copydoc FieldIO::Import
    LIB_UTILITIES_EXPORT virtual void v_Import(
        const std::string &infilename,
        std::vector<FieldDefinitionsSharedPtr> &fielddefs,
        std::vector<std::vector<NekDouble> >
            &fielddata                     = NullVectorNekDoubleVector,
        FieldMetaDataMap &fieldinfomap     = NullFieldMetaDataMap,
        const Array<OneD, int> &ElementIDs = NullInt1DArray) = 0;

    /// @copydoc FieldIO::ImportFieldMetaData
    LIB_UTILITIES_EXPORT virtual DataSourceSharedPtr v_ImportFieldMetaData(
        const std::string &filename, FieldMetaDataMap &fieldmetadatamap) = 0;
};

typedef std::shared_ptr<FieldIO> FieldIOSharedPtr;

/**
 * @brief Write out the field information to the file @p outFile.
 *
 * @param outFile       Output filename
 * @param fielddefs     Field definitions that define the output
 * @param fielddata     Binary field data that stores the output corresponding
 *                      to @p fielddefs.
 * @param fieldinfomap  Associated field metadata map.
 */
inline void FieldIO::Write(const std::string                      &outFile,
                           std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                           std::vector<std::vector<NekDouble> > &fielddata,
                           const FieldMetaDataMap                 &fieldinfomap,
                           const bool                              backup)
{
    v_Write(outFile, fielddefs, fielddata, fieldinfomap, backup);
}

/**
 * @brief Read field information from the file @p infilename.
 *
 * @param infilename    Input filename (or directory if parallel format)
 * @param fielddefs     On return contains field definitions as read from the
 *                      input.
 * @param fielddata     On return, contains binary field data that stores the
 *                      input corresponding to @p fielddefs.
 * @param fieldinfo     On returnm, contains the associated field metadata map.
 * @param ElementIDs    Element IDs that lie on this processor, which can be
 *                      optionally supplied to avoid reading the entire file on
 *                      each processor.
 */
inline void FieldIO::Import(const std::string                      &infilename,
                            std::vector<FieldDefinitionsSharedPtr> &fielddefs,
                            std::vector<std::vector<NekDouble> >   &fielddata,
                            FieldMetaDataMap                       &fieldinfo,
                            const Array<OneD, int>                 &ElementIDs)
{
    v_Import(infilename, fielddefs, fielddata, fieldinfo, ElementIDs);
}

/**
 * @brief Import the metadata from a field file.
 *
 * @param filename          Input filename.
 * @param fieldmetadatamap  On return contains the field metadata map from @p
 *                          filename.
 */
inline DataSourceSharedPtr FieldIO::ImportFieldMetaData(
    const std::string &filename, FieldMetaDataMap &fieldmetadatamap)
{
    return v_ImportFieldMetaData(filename, fieldmetadatamap);
}

}
}
#endif
