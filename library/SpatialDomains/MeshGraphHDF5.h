////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description:  This file contains the base class specification for the
//                Geometry class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MGHDF5_H
#define NEKTAR_SPATIALDOMAINS_MGHDF5_H

#include <SpatialDomains/MeshEntities.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/BasicUtils/H5.h>

namespace Nektar
{
namespace SpatialDomains
{

class MeshGraphHDF5 : public MeshGraph
{
public:
    MeshGraphHDF5()
    {
    }

    SPATIAL_DOMAINS_EXPORT virtual void WriteGeometry(
        std::string &outfilename,
        bool defaultExp = false,
        const LibUtilities::FieldMetaDataMap &metadata
                                         = LibUtilities::NullFieldMetaDataMap);

    virtual ~MeshGraphHDF5()
    {
    }

    static MeshGraphSharedPtr create()
    {
        return MemoryManager<MeshGraphHDF5>::AllocateSharedPtr();
    }

    static std::string className;

protected:
    SPATIAL_DOMAINS_EXPORT virtual void ReadGeometry(
        DomainRangeShPtr rng,
        bool             fillGraph);
    SPATIAL_DOMAINS_EXPORT virtual void PartitionMesh(
        LibUtilities::SessionReaderSharedPtr session);

private:

    void ReadCurveMap(
        CurveMap                      &curveMap,
        std::string                    dsName,
        const std::unordered_set<int> &readIds);
    void ReadDomain();
    void ReadComposites();

    template<class T>
    void WriteGeometryMap(std::map<int, std::shared_ptr<T>> &geomMap,
                          std::string datasetName);

    template<class T, typename DataType> void ReadGeometryData(
        std::map<int, std::shared_ptr<T>>      &geomMap,
        std::string                             dataSet,
        const std::unordered_set<int>          &readIds,
        std::vector<int>                       &ids,
        std::vector<DataType>                  &geomData);
    template<class T, typename DataType> void FillGeomMap(
        std::map<int, std::shared_ptr<T>> &geomMap,
        const CurveMap                    &curveMap,
        std::vector<int>                  &ids,
        std::vector<DataType>             &geomData);
    template<class T, typename DataType> void ConstructGeomObject(
        std::map<int, std::shared_ptr<T>> &geomMap, int id,
        DataType *data, CurveSharedPtr curve);

    CompositeDescriptor CreateCompositeDescriptor(
        std::unordered_map<int, int> &id2row);

    void WriteCurveMap(CurveMap &curves,
                       std::string dsName,
                       MeshCurvedPts &curvedPts,
                       int &ptOffset,
                       int &newIdx);
    void WriteCurvePoints(MeshCurvedPts &curvedPts);

    void WriteComposites(CompositeMap &comps);
    void WriteDomain(std::vector<CompositeMap> &domain);

    std::string m_hdf5Name;
    LibUtilities::H5::FileSharedPtr m_file;
    LibUtilities::H5::PListSharedPtr m_readPL;
    LibUtilities::H5::GroupSharedPtr m_mesh;
    LibUtilities::H5::GroupSharedPtr m_maps;

    static const unsigned int FORMAT_VERSION;
};

} // end of namespace
} // end of namespace

#endif
