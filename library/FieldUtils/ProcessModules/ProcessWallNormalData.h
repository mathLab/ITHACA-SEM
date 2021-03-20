////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessWallNormalData.h
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
//  Description: Get the wall-normal data at a given origin.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSWALLNORMALDATA
#define FIELDUTILS_PROCESSWALLNORMALDATA

#include "ProcessBoundaryExtract.h"

namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief This processing module calculates the wall shear stress and adds it
 * as an extra-field to the output file, and writes it to a surface output file.
 */
class ProcessWallNormalData : public ProcessBoundaryExtract
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessWallNormalData>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessWallNormalData(FieldSharedPtr f);
    virtual ~ProcessWallNormalData();

    /// Write mesh to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessWallNormalData";
    }

    virtual std::string GetModuleDescription()
    {
        return "Get the wall-normal data at a given origin.";
    }

protected:

private:
    int m_spacedim;

    /**
    * @brief Project a single point along the given direction to a plane
    * @param gloCoord     Global coordinate of the point. size=3.
    * @param projDir      Projection direction, which is also the normal vector 
    *                     of the target plane. size=3, norm=1.
    * @param distToOrig   The distance from the origin (0,0,0) to the target plane.
    * @param projGloCoord The global coordinate of the projecion result.
    */ 
    void ProjectPoint(
        const Array<OneD, const NekDouble > & gloCoord,
        const Array<OneD, const NekDouble > & projDir,
        const NekDouble distToOrig,
        Array<OneD, NekDouble > & projGloCoord);

    /**
    * @brief Project a single point along the given direction to a plane
    * @param pts          Global coordinate of the vertices of the elmt. size=2/3.
    * @param projDir      Projection direction, which is also the normal vector of 
    *                     the target plane. size=3, norm=1.
    * @param distToOrig   The distance from the origin (0,0,0) to the target plane.
    * @param projPts      The global coordinate of the projecion result.
    */ 
    void ProjectVertices(
        const Array<OneD, const Array<OneD, NekDouble> > & pts,
        const Array<OneD, const NekDouble > & projDir,
        const NekDouble distToOrig,
        Array<OneD, Array<OneD, NekDouble> > & projPts);

    /**
    * @brief Determine if the projected point is inside the projected element.
    * @param projGloCoord The global coordinate of the projected single point.
    * @param projPts      The global coordinate of the projected vertices,size=2/3
    * @param projDir      Projection direction, which is used as the reference
    *                     direction. size=3, norm=1. 
    * @param paralTol     Tolerence to check if two vectors are parallel.
    * @param angleTol     Tolerence to check if the total angle is 2*pi.
    * @return             Inside (true) or not (false)
    */ 
    bool isInProjectedArea2D(
        const Array<OneD, const NekDouble > & projGloCoord,
        const Array<OneD, const Array<OneD, NekDouble> > & projPts,
        const NekDouble paralTol = 1.0e-12);

    bool isInProjectedArea3D(
        const Array<OneD, const NekDouble > & projGloCoord,
        const Array<OneD, const Array<OneD, NekDouble> > & projPts,
        const Array<OneD, const NekDouble > & projDir,
        const NekDouble paralTol = 1.0e-12,
        const NekDouble angleTol = 1.0e-6);

    /**
     * @brief Use iteration to get the locCoord. This routine should be used after
     *        we have checked the projected point is inside the projected element.
     * @param bndGeom      Geometry to get the xmap.
     * @param gloCoord     Global coordinate of the point. size=3.
     * @param pts          Global coordinate of the vertices of the elmt. size=2/3.
     * @param dieUse       The main direction(s) used to compute local coordinate
     * @param locCoord     Iteration results for local coordinate(s) 
     * @param dist         Returned distance in physical space if the collapsed 
     *                     locCoord is out of range [-1,1].
     * @param iterTol      Tolerence for iteration.
     * @param iterMax      Maximum iteration steps
     * @return             Converged (true) or not (false)
     */
    bool BisectionForLocCoordOnBndElmt(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & gloCoord,
        const Array<OneD, const Array<OneD, NekDouble> > & pts,
        const Array<OneD, const int > & dirUse,
        Array<OneD, NekDouble > & locCoord,
        const NekDouble iterTol = 1.0e-8,
        const int iterMax = 51);
    
    bool NewtonIterForLocCoordOnBndElmt(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble> & gloCoord,
        const Array<OneD, const Array<OneD, NekDouble> > & pts,
        const Array<OneD, const int > & dirUse,
        Array<OneD, NekDouble> & locCoord,
        NekDouble & dist,
        const NekDouble iterTol = 1.0e-8,
        const int iterMax = 51);
    
    /**
    * @brief Check if a point can be projected onto an oundary element in a given
    *        direction. If yes, give the local coordinates of the projected point.
    *        we have checked the projected point is inside the projected element.
    * @param bndGeom      Pointer to the geometry of the boundary element.
    * @param gloCoord     Global coordinate of the point. size=3.
    * @param projDir      Projection direction, which is used as the reference
    *                     direction in the 3D routine. size=3, norm=1. 
    * @param locCoord     Iteration results for local coordinates (if inside).
    * @param projDist     Projection distance betweem the point to the wall point.
    * @param geomTol      Disntance to check if the wall point is desired.
    * @param iterTol      Tolerence for iteration.
    * @return             Inside (true) or not (false)
    */
    bool BndElmtContainsPoint(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & gloCoord,
        const Array<OneD, const NekDouble > & projDir,
        Array< OneD, NekDouble > & locCoord,
        NekDouble & projDist,
        const NekDouble geomTol = 1.0,
        const NekDouble iterTol = 1.0e-8);

    /**
    * @brief Get the normals for a given locCoord
    * @param bndGeom      Pointer to the geometry of the boundary element.
    * @param locCoord     Iteration results for local coordinates (if inside).
    * @param normals      Wall normal as the result
    */
    void GetNormals(
        SpatialDomains::GeometrySharedPtr bndGeom,
        const Array<OneD, const NekDouble > & locCoord, 
        Array< OneD, NekDouble > & normals);

};
}
}

#endif
