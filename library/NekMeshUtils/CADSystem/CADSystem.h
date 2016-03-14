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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_CADSYSTEM_CADSYSTEM
#define NekMeshUtils_CADSYSTEM_CADSYSTEM

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <NekMeshUtils/CADSystem/OpenCascade.h>

namespace Nektar
{
namespace NekMeshUtils
{

class CADVert;
typedef boost::shared_ptr<CADVert> CADVertSharedPtr;
class CADCurve;
typedef boost::shared_ptr<CADCurve> CADCurveSharedPtr;
class CADSurf;
typedef boost::shared_ptr<CADSurf> CADSurfSharedPtr;
struct EdgeLoop;

/**
 * @brief Base class for CAD interface system.
 *
 * A class which can load and interact with CAD for Nektar++ using OpenCascade.
 * This class contains maps to subclasses surface and curves.
 */
class CADSystem
{
public:
    friend class MemoryManager<CADSystem>;

    /**
     * @brief Default constructor.
     */
    CADSystem(const std::string &name) : m_name(name)
    {
    }

    /**
     * @brief Return the name of the CAD system.
     */
    std::string GetName();

    /**
     * @brief Initialises CAD and makes surface, curve and vertex maps.
     *
     * @return true if completed successfully
     */
    bool LoadCAD();

    /**
     * @brief Reports basic properties to screen.
     */
    void Report();

    /**
     * @brief Returns bounding box of the domain.
     *
     * Gets the bounding box of the domain by considering the start and end
     * points of each curve in the geometry.
     *
     * @return Array with 6 entries: xmin, xmax, ymin, ymax, zmin and zmax.
     */
    Array<OneD, NekDouble> GetBoundingBox();

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
        std::map<int, CADCurveSharedPtr>::iterator search = m_curves.find(i);
        ASSERTL0(search != m_curves.end(), "curve does not exist");

        return search->second;
    }

    /**
     * @brief Gets a surface from the map.
     */
    CADSurfSharedPtr GetSurf(int i)
    {
        std::map<int, CADSurfSharedPtr>::iterator search = m_surfs.find(i);
        ASSERTL0(search != m_surfs.end(), "surface does not exist");

        return search->second;
    }

    /**
     * @brief Gets map of all vertices
     */
    std::map<int, CADVertSharedPtr> GetVerts()
    {
        return m_verts;
    }

    /**
     * @brief Gets number of vertices
     */
    int GetNumVerts()
    {
        return m_verts.size();
    }

    /**
     * @brief based on location in space, uses opencascade routines to
     * determin if the point is within the domain. This routine is slow
     * and should be used sparingly, it is smart enough to take and form
     * of geometry
     */
    bool InsideShape(Array<OneD, NekDouble> loc);

private:
    /// Function to add curve to CADSystem::m_verts.
    void AddVert(int i, TopoDS_Shape in);
    /// Function to add curve to CADSystem::m_curves.
    void AddCurve(int i, TopoDS_Shape in, int fv, int lv);
    /// Function to add surface to CADSystem::m_surfs.
    void AddSurf(int i, TopoDS_Shape in, std::vector<EdgeLoop> ein);
    /// Name of cad file to be opened, including file extension.
    std::string m_name;
    /// Map of curves
    std::map<int, CADCurveSharedPtr> m_curves;
    /// Map of surfaces
    std::map<int, CADSurfSharedPtr> m_surfs;
    /// Map of vertices
    std::map<int, CADVertSharedPtr> m_verts;
    /// OCC master object
    TopoDS_Shape shape;
};

typedef boost::shared_ptr<CADSystem> CADSystemSharedPtr;
}
}

#endif
