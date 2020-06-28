////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_CADSYSTEM_CADSYSTEM
#define NekMeshUtils_CADSYSTEM_CADSYSTEM

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <LibUtilities/BasicUtils/NekFactory.hpp>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>

#include <NekMeshUtils/CADSystem/CADVert.h>
#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

#include "CADObject.h"

namespace Nektar
{
namespace NekMeshUtils
{

// forward declorators
class CADVert;
typedef std::shared_ptr<CADVert> CADVertSharedPtr;
class CADCurve;
typedef std::shared_ptr<CADCurve> CADCurveSharedPtr;
class CADSurf;
typedef std::shared_ptr<CADSurf> CADSurfSharedPtr;

/**
 * @brief Base class for CAD interface system.
 *
 * A class which can load and interact with CAD for Nektar++.
 * This class contains maps to subclasses surface and curves.
 */
class CADSystem
{
public:
    friend class MemoryManager<CADSystem>;

    /**
     * @brief Default constructor.
     */
    CADSystem(std::string name, std::string engine)
        : m_name(name), m_engine(engine)
    {
    }

    virtual ~CADSystem()
    {
    }

    /**
     * @brief Return the name of the CAD file.
     */
    std::string GetName()
    {
        return m_name;
    }

    /**
     * @brief Set the CAD engine up to handle 2D geometries explicitly.
     */
    void Set2D()
    {
        m_2d = true;
    }

    /**
     * @brief Returns whether the CAD engine is set in 2D mode.
     */
    bool Is2D()
    {
        return m_2d;
    }

    /**
     * @brief Gets a configuration option.
     *
     * @param key  Configuration key.
     *
     * @return The configuration value.
     */
    const std::string &GetConfig(const std::string &key) const
    {
        auto it = m_config.find(key);
        return it->second;
    }

    /**
     * @brief Sets a configuration option @p key to @p value.
     *
     * @param key    Configuration key.
     * @param value  The configuration value.
     */
    void SetConfig(const std::string &key, const std::string &value)
    {
        m_config[key] = value;
    }

    /**
     * @brief Sets the verbose flag for additional output.
     */
    void SetVerbose(bool verbose = true)
    {
        m_verbose = verbose;
    }

    /**
     * @brief Initialises CAD and makes surface, curve and vertex maps.
     *
     * @return true if completed successfully
     */
    NEKMESHUTILS_EXPORT virtual bool LoadCAD() = 0;

    /**
     * @brief Returns bounding box of the domain.
     *
     * Gets the bounding box of the domain by considering the start and end
     * points of each curve in the geometry.
     *
     * @return Array with 6 entries: xmin, xmax, ymin, ymax, zmin and zmax.
     */
    NEKMESHUTILS_EXPORT virtual Array<OneD, NekDouble> GetBoundingBox() = 0;

    /**
     * @brief Get the number of surfaces.
     */
    int GetNumSurf()
    {
        return m_surfs.size();
    }

    /**
     * @brief Get the number of curves.
     */
    int GetNumCurve()
    {
        return m_curves.size();
    }

    /**
     * @brief Gets a curve from the map.
     */
    CADCurveSharedPtr GetCurve(int i)
    {
        auto search = m_curves.find(i);
        ASSERTL1(search != m_curves.end(), "curve does not exist");

        return search->second;
    }

    /**
     * @brief Gets a surface from the map.
     */
    CADSurfSharedPtr GetSurf(int i)
    {
        auto search = m_surfs.find(i);
        ASSERTL1(search != m_surfs.end(), "surface does not exist");

        return search->second;
    }

    /**
     * @brief Gets a vert from the map.
     */
    CADVertSharedPtr GetVert(int i)
    {
        auto search = m_verts.find(i);
        ASSERTL1(search != m_verts.end(), "vert does not exist");

        return search->second;
    }

    /**
     * @brief Gets map of all vertices.
     */
    std::map<int, CADVertSharedPtr> GetVerts()
    {
        return m_verts;
    }

    /**
     * @brief Gets number of vertices
     */
    int GetNumVerts() const
    {
        return m_verts.size();
    }

    /**
     * @brief Gets void points for tetrahedral meshing.
     */
    std::vector<Array<OneD, NekDouble>> GetVoidPoints()
    {
        return m_voidPoints;
    }

    /**
     * @brief Sets void points for tetrahedral meshing.
     */
    void SetVoidPoints(const std::vector<Array<OneD, NekDouble>> &voidPts)
    {
        m_voidPoints = voidPts;
    }

    /**
     * @brief Return the vector of translation from one curve to another to
     * allow for periodic mesh generation in 2D.
     */
    NEKMESHUTILS_EXPORT Array<OneD, NekDouble> GetPeriodicTranslationVector(
        int first, int second);

    /**
     * @brief Return the name of the CAD system engine (e.g. OCE).
     */
    std::string GetEngine() const
    {
        return m_engine;
    }

protected:
    /// Name of CAD file
    std::string m_name;
    /// Name of CAD engine.
    std::string m_engine;
    /// Map of curves
    std::map<int, CADCurveSharedPtr> m_curves;
    /// Map of surfaces
    std::map<int, CADSurfSharedPtr> m_surfs;
    /// Map of vertices
    std::map<int, CADVertSharedPtr> m_verts;
    /// Verbosity
    bool m_verbose = false;
    /// 2D cad flag
    bool m_2d = false;
    /// Configuration options which might be used per-engine, serialised in
    /// strings.
    std::map<std::string, std::string> m_config;
    /// Points contained within volume voids for tetrahedralisation
    std::vector<Array<OneD, NekDouble>> m_voidPoints;

    /**
     * @brief Reports basic properties to screen.
     */
    void Report()
    {
        std::cout << std::endl << "CAD report:" << std::endl;
        std::cout << "\tCAD has: " << m_verts.size() << " verts." << std::endl;
        std::cout << "\tCAD has: " << m_curves.size() << " curves."
                  << std::endl;
        std::cout << "\tCAD has: " << m_surfs.size() << " surfaces."
                  << std::endl;
    }
};

typedef std::shared_ptr<CADSystem> CADSystemSharedPtr;
typedef LibUtilities::NekFactory<std::string, CADSystem, std::string>
    EngineFactory;

EngineFactory &GetEngineFactory();
}
}

#endif
