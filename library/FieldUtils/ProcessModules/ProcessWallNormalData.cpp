////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessWallNormalData.cpp
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

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

#include "ProcessWallNormalData.h"


using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessWallNormalData::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "wallNormalData"),
    ProcessWallNormalData::create,
    "Export data in wall-normal direction from a single point on the wall.");

ProcessWallNormalData::ProcessWallNormalData(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
    f->m_writeBndFld = false; // turned on in upstream ProcessBoundaryExtract

    m_config["xorig"]   = ConfigOption(false, "0.5,0.0,0.0",
                          "The point to be projected onto the wall to get the \
                          sampling origin. default=[0.5,0,0]");
    m_config["projDir"] = ConfigOption(false, "0.0,1.0,0.0",
                          "The direction to project the point onto the wall to \
                          get the sampling origin. default=[0,1,0]");
    m_config["maxDist"] = ConfigOption(false, "1.0",
                          "Distance to limit projection distance to find the \
                          desired sampling origin. defalut=1.0");
    m_config["distH"]   = ConfigOption(false, "0.01",
                          "Sampling distance along the wall normal at the \
                          sampling origin. default=0.1");
    m_config["nptsH"]   = ConfigOption(false, "5",
                          "Number of sampling points along the wall normal. \
                          default=5");
    m_config["d"]       = ConfigOption(false, "0.1",
                          "The parameter that controls the sampling points' \
                          distribution. d should be in the range (0,inf). d \
                          in (0,0.95] gives controled points; d in (0.95,inf) \
                          gives evenly spaced points");
}

ProcessWallNormalData::~ProcessWallNormalData()
{
}


/*Note
* This module is used to get field data in the wall-normal direction.
* The input cases can be 2D, 2.5D and 3D. 
* The data will be exported with .pts extension.
*
* The user defined parameters are: bnd, xorig, searchDir, maxDist, distH,
* nptsH, and d.
*  - bnd=0 is the boundary id. This boundary should contain the desired origin.
*  - xorig="x,y,z" are the coordinates of the input sampling origin. They are
*    used as the references to get exact origin on the wall.
*  - projDir="0,1,0" is the projection direction to find the point on the wall
*  - maxDist=1.0 is the distance that constrains the projection on the wall.
*    The projection distance should be smaller than this distance. This
*    parameter is set in case the boundary is wavey so that multiple
*    projections exist at the same time.
*  - distH=0.01 is the sampling depth in the wall-normal direction.
*  - nptsH=11 is the number of sampling points along wall-normal direction.
*  - d=0.1 is a destribution control parameter of the sampling points. It
*    should be in the range (0,inf). d in range (0,0.95] for controlled array.
*    d>0.95 gives evenly spaced array.
*/
void ProcessWallNormalData::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);
    
    // Get dim to store data
    const int nfields   = m_f->m_variables.size();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    m_spacedim          = nCoordDim + m_f->m_numHomogeneousDir;
    const int nBndLcoordDim = nCoordDim - 1;
    const int totVars   = m_spacedim + m_f->m_variables.size();


    // Initialize the sampling parameters
    std::vector<NekDouble> xorig, searchDir;
    ASSERTL0(ParseUtils::GenerateVector(m_config["xorig"].as<string>(), xorig),
             "Failed to interpret origin coordinates string");
    ASSERTL0(ParseUtils::GenerateVector(m_config["projDir"].as<string>(), searchDir),
             "Failed to interpret search direction string");
    const NekDouble maxDist = m_config["maxDist"].as<NekDouble>();
    const NekDouble distH   = m_config["distH"].as<NekDouble>();
    const int       nptsH   = m_config["nptsH"].as<int>();
    const NekDouble delta   = m_config["d"].as<NekDouble>();


    Array<OneD, NekDouble> orig(3), projDir(3); // gloCoord of the origin
    for (int i=0; i<3; ++i)
    {
        orig[i]    = (xorig.size()>=(i+1))     ? xorig[i]     : 0.0;
        projDir[i] = (searchDir.size()>=(i+1)) ? searchDir[i] : 0.0;
    }
    if (nCoordDim==2)
    {
        projDir[2] = 0.0;
    }
    Vmath::Smul(3, 1.0/sqrt(Vmath::Dot(3, projDir, 1, projDir, 1)), 
                projDir, 1, projDir, 1);

    // Update z according to the closest plane in 2.5D cases
    if (m_f->m_numHomogeneousDir==1)
    {
        int nPlanes    = m_f->m_exp[0]->GetHomogeneousBasis()->GetZ().size();
        NekDouble lHom = m_f->m_exp[0]->GetHomoLen();
        if (orig[2]<0.0 || orig[2]>lHom)
        {
            orig[2] = 0.0;
        }
        else
        {
            NekDouble dZ = lHom / nPlanes;
            NekDouble zTmp, zCur=0.0, distTmp, distCur = 999.0;
            for(int i=0; i<=nPlanes; ++i)
            {
                zTmp    = dZ * i;
                distTmp = fabs(orig[2] - zTmp); 
                if(distTmp < distCur)
                {
                    distCur = distTmp;
                    zCur    = zTmp;
                    if(i==nPlanes)
                    {
                        zCur = 0.0;
                    }
                }
            }
            orig[2] = zCur;
        }
    }


    // Get the bnd id
    SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                           m_f->m_exp[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection bregions =
        bcs.GetBoundaryRegions();
    map<int, int> BndRegionMap;
    int cnt = 0;
    for (auto &breg_it : bregions)
    {
        BndRegionMap[breg_it.first] = cnt++;
    }
    int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[0]];

    // Get expansion list for boundary
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(nfields); 
    for (int i = 0; i < nfields; ++i) {
        BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(bnd);
    }
    

    //-------------------------------------------------------------------------
    // Find the element that contains the origin, and give the precise origin
    SpatialDomains::GeometrySharedPtr bndGeom;
    Array<OneD, NekDouble> locCoord(nBndLcoordDim, -999.0);
    NekDouble projDist;
    int elmtid;
    bool isInside = false;

    // Search and get precise locCoord
    const int nElmts = BndExp[0]->GetNumElmts();
    for (elmtid=0; elmtid<nElmts; ++elmtid) //nElmts
    {    
        bndGeom  = BndExp[0]->GetExp(elmtid)->GetGeom();         
        isInside = BndElmtContainsPoint(bndGeom, orig, projDir, locCoord,
                                        projDist, maxDist);
        if (isInside) 
        {
            break;
        }

    }
    ASSERTL0(isInside, "Failed to find the sampling origin on the boundary."); 
 

    // Then Update the precise sampling position
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    const int npts = bndXmap->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);

    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i);              // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
        orig[i] = bndXmap->PhysEvaluate(locCoord, pts[i]); // Update wall point
    }
 
    
    // Get outward-pointing normal vectors for all quadrature points on bnd
    // Use these normals as diretion references
    Array<OneD, Array<OneD, NekDouble> > normalsQ; 
    m_f->m_exp[0]->GetBoundaryNormals(bnd, normalsQ);

    // Get the normals in the element
    // Set point key to get point id
    Array<OneD, LibUtilities::PointsKey> from_key(nBndLcoordDim);
    int from_nPtsPerElmt = 1;
    for (int i=0; i<nBndLcoordDim; ++i)
    {
        from_key[i] = BndExp[0]->GetExp(elmtid)->GetBasis(i)->GetPointsKey();
        from_nPtsPerElmt *= from_key[i].GetNumPoints();
    }

    //Get ref normal
    Array< OneD, NekDouble > normalsRef(3, 0.0);
    int refId = elmtid*from_nPtsPerElmt + 0; // the 1st normal in the bnd elmt
    for(int i=0;i<m_spacedim; ++i)
    {
        normalsRef[i] = normalsQ[i][refId];
    }
    
    // Get the precise normals and correct the direction to be inward
    Array< OneD, NekDouble > normals(3, 0.0);
    GetNormals(bndGeom, locCoord, normals);

    if(Vmath::Dot(3, normals, normalsRef) > 0.0)
    {
        Vmath::Neg(3, normals, 1);
    }
    
    // Output the info if -v is set
    if (m_f->m_verbose)
    {
        cout << "------ wallNormalData module ------\n";
        cout << "Input point:\n";
        cout << "  - [Px,Py,Pz] = [" << xorig[0] << ", " << xorig[1];
        if (xorig.size()>=3)
        {
            cout << ", " << xorig[2] << "]\n";
        }
        else
        {
            cout << ", " << 0.0 << "]\n";
        }
        cout << "Projection direction:\n";
        cout << "  - [vx,vy,vz] = ["
             << projDir[0] << ", " << projDir[1] << ", " << projDir[2] << "]\n";
        cout << "Sampling origin on the wall:\n";
        cout << "  - [Ox,Oy,Oz] = ["
             << orig[0] << ", " << orig[1] << ", " << orig[2] << "]\n";
        cout << "Normals at the origin:\n";
        cout << "  - [nx,ny,nz] = ["
             << normals[0] << ", " << normals[1] << ", " << normals[2] << "]\n";
        cout << "Ref normals (at quadrature points in the projected element):\n";
        for (int i=0; i<from_nPtsPerElmt; ++i)
        {
            cout << "  - " << i << ": [nx,ny,nz] = ["
                 << -normalsQ[0][elmtid*from_nPtsPerElmt+i] << ", "
                 << -normalsQ[1][elmtid*from_nPtsPerElmt+i];
            if (m_spacedim==3)
            {
                cout << ", " << -normalsQ[2][elmtid*from_nPtsPerElmt+i] << "]\n";
            }
            else
            {
                cout << "]\n";
            }
        }
        cout << endl;
    }

    
    //-------------------------------------------------------------------------
    // Set the depth of sampling array
    Array<OneD, NekDouble> h(nptsH, 0.0);
    if (delta > 0 && delta <= 0.95)
    {
        // Use expression in Agrawal's paper:
        // h = 1-tanh((1-ksi)*atanh(sqrt(1-delta)))/sqrt(1-delta), ksi in [0,1]
        NekDouble tmp1;
        const NekDouble tmp2 = 1.0/(static_cast<NekDouble>(nptsH)-1.0);
        const NekDouble tmp3 = sqrt(1.0-delta);
        const NekDouble tmp4 = atanh(tmp3);
        const NekDouble tmp5 = 1.0/tmp3;
        for (int i=1; i<nptsH; ++i)
        {
            tmp1 = 1.0 - i * tmp2; // tmp1 = 1-ksi
            h[i] = 1 - tanh(tmp1*tmp4)*tmp5;
            h[i] *= distH; // physical distance in normal direction
        }
    }
    else if (delta > 0.95)
    {
        // Evenly spaced array
        const NekDouble tmp1 = 1.0/(static_cast<NekDouble>(nptsH)-1.0);
        for (int i=1; i<nptsH; ++i)
        {
            h[i] = i * tmp1;
            h[i] *= distH; // physical distance in normal direction
        }

    }
    else{
        ASSERTL0(false, "Input error. Delta needs to be greater than 0.0.");
    }

    // Set pts coordinates and interpoate the data
    Array<OneD, Array<OneD, NekDouble> > ptsH(totVars);
    for (int i=0; i<totVars; ++i)
    {
        ptsH[i] = Array<OneD, NekDouble>(nptsH, 0.0);
    }

    for(int i=0; i<m_spacedim; ++i)
    {
        for(int j=0; j<nptsH; ++j)
        {
            ptsH[i][j] = orig[i] + h[j] * normals[i]; // x0+dist*nx 
        }
    }

    m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
        m_spacedim, m_f->m_variables, ptsH, LibUtilities::NullPtsInfoMap);

    Interpolator interp;
    interp.Interpolate(m_f->m_exp, m_f->m_fieldPts, NekConstants::kNekUnsetDouble);
}



/**
 * @brief Project a single point along the given direction to a plane
 * @param gloCoord     Global coordinate of the point. size=3.
 * @param projDir      Projection direction, which is also the normal vector of 
 *                     the target plane. size=3, norm=1.
 * @param distToOrig   The distance from the origin (0,0,0) to the target plane
 * @param projGloCoord The global coordinate of the projecion result.
 */ 
void ProcessWallNormalData::ProjectPoint(
    const Array<OneD, const NekDouble > & gloCoord,
    const Array<OneD, const NekDouble > & projDir,
    const NekDouble distToOrig,
    Array<OneD, NekDouble > & projGloCoord)
{
    NekDouble tmp1;
    Array<OneD, NekDouble > tmp2(3);
    
    tmp1 = Vmath::Dot(3, gloCoord, 1, projDir, 1);        // |xn| = x dot n
    Vmath::Smul(3, distToOrig-tmp1, projDir, 1, tmp2, 1); // d_xn = (dist-|xn|) n
    Vmath::Vadd(3, gloCoord, 1, tmp2, 1, projGloCoord,1); // x' = x + d_xn
}


/**
 * @brief Project the vertices of the elmt along the given direction to a plane
 * @param pts          Global coordinate of the vertices of the elmt. size=2/3.
 * @param projDir      Projection direction, which is also the normal vector of 
 *                     the target plane. size=3, norm=1.
 * @param distToOrig   The distance from the origin (0,0,0) to the target plane.
 * @param projPts      The global coordinate of the projecion result.
 */ 
void ProcessWallNormalData::ProjectVertices(
    const Array<OneD, const Array<OneD, NekDouble> > & pts,
    const Array<OneD, const NekDouble > & projDir,
    const NekDouble distToOrig,
    Array<OneD, Array<OneD, NekDouble> > & projPts)
{
    const int nCoordDim = pts.size();   // 2 for 2.5D cases
    const int npts      = pts[0].size();

    NekDouble tmp1;
    Array<OneD, NekDouble > singlePnt(nCoordDim), tmp2(nCoordDim);

    for (int i=0; i<npts; ++i)
    {
        // Get a point
        for (int j=0; j<nCoordDim; ++j)
        {
            singlePnt[j] = pts[j][i];
        }
        
        // Projection
        tmp1 = Vmath::Dot(nCoordDim, singlePnt, 1, projDir, 1);       
        Vmath::Smul(nCoordDim, distToOrig-tmp1, projDir, 1, tmp2, 1);
        Vmath::Vadd(nCoordDim, singlePnt, 1, tmp2, 1, singlePnt,1);

        // Save to projPts
        for (int j=0; j<nCoordDim; ++j)
        {
            projPts[j][i] = singlePnt[j];
        }
    }

}


/**
 * @brief Determine if the projected point is inside the projected element.
 * @param projGloCoord The global coordinate of the projected single point.
 * @param projPts      The global coordinate of the projected vertices,size=2/3
 * @param projDir      Projection direction, which is used as the reference
 *                     direction in the 3D routine. size=3, norm=1. 
 * @param paralTol     Tolerence to check if two vectors are parallel.
 * @param angleTol     Tolerence to check if the total angle is 2*pi.
 * @return             Inside (true) or not (false)
 */
bool ProcessWallNormalData::isInProjectedArea2D(
    const Array<OneD, const NekDouble > & projGloCoord,
    const Array<OneD, const Array<OneD, NekDouble> > & projPts,
    const NekDouble paralTol)
{
    const NekDouble npts = projPts[0].size();
    
    Array<OneD, NekDouble> vec1(2, 0.0), vec2(2, 0.0);

    vec1[0] = projPts[0][0]      - projGloCoord[0];
    vec1[1] = projPts[1][0]      - projGloCoord[1];
    vec2[0] = projPts[0][npts-1] - projGloCoord[0];
    vec2[1] = projPts[1][npts-1] - projGloCoord[1];
    
    Vmath::Smul(2, 1.0/sqrt(Vmath::Dot(2, vec1, 1, vec1, 1)),
                vec1, 1, vec1, 1);
    Vmath::Smul(2, 1.0/sqrt(Vmath::Dot(2, vec2, 1, vec2, 1)),
                vec2, 1, vec2, 1);

    if( (fabs(vec1[0]+vec2[0]) + fabs(vec1[1]+vec2[1])) < paralTol )
    {
        // On the line segement, true
        return true;
    }
    else
    {
        return false;
    }

}


bool ProcessWallNormalData::isInProjectedArea3D(
    const Array<OneD, const NekDouble > & projGloCoord,
    const Array<OneD, const Array<OneD, NekDouble> > & projPts,
    const Array<OneD, const NekDouble > & projDir,
    const NekDouble paralTol,
    const NekDouble angleTol)
{

    //Generate a polygen (quad), re-order the points on the edge
    const NekDouble npts  = projPts[0].size();
    const int nptsEdge    = sqrt(npts);  // num of points on edge for geo representation
    const int nptsPolygon = 4*nptsEdge-4;
    Array<OneD, Array<OneD, NekDouble> > ptsPolygon(3);
    for (int i=0; i<3; ++i) 
    {
        ptsPolygon[i] = Array<OneD, NekDouble>(nptsPolygon);

        for(int j=0; j<nptsEdge; ++j)
        {
            ptsPolygon[i][j] = projPts[i][j];
        }
        for(int j=0; j<nptsEdge-2; ++j)
        {
            ptsPolygon[i][nptsEdge+j] = projPts[i][(j+2)*nptsEdge-1];
        }
        for(int j=0; j<nptsEdge; ++j)
        {
            ptsPolygon[i][2*nptsEdge-2+j] = projPts[i][npts-1-j];
        }
        for(int j=0; j<nptsEdge-2; ++j)
        {
            ptsPolygon[i][3*nptsEdge-2+j] = projPts[i][nptsEdge*(nptsEdge-j-2)];
        }
    }
        

    // Determine relation using angle method (sum=2*pi)
    NekDouble angleCos, angleAbs, angleSign, angleSum = 0.0;
    Array<OneD, NekDouble> vec1(3, 0.0), vec2(3, 0.0), vec3(3, 0.0);
    int id1, id2;

    for(int i=0; i<nptsPolygon; ++i)
    {
        id1 = i;
        id2 = (id1==(nptsPolygon-1)) ? 0 : (id1+1);

        for (int j=0; j<3; ++j)
        {
            vec1[j] = ptsPolygon[j][id1] - projGloCoord[j];
            vec2[j] = ptsPolygon[j][id2] - projGloCoord[j];
        }

        Vmath::Smul(3, 1.0/sqrt(Vmath::Dot(3, vec1, 1, vec1, 1)),
                    vec1, 1, vec1, 1);
        Vmath::Smul(3, 1.0/sqrt(Vmath::Dot(3, vec2, 1, vec2, 1)),
                    vec2, 1, vec2, 1);


        if( ( fabs(vec1[0]+vec2[0]) + fabs(vec1[1]+vec2[1]) +
              fabs(vec1[2]+vec2[2]) ) < paralTol )
        {
            // On the line segement, true
            // This branch is used to aviod angle=pi but angleSign=-1 (not 1)
            return true;
        }
        else
        {
            // Off the line segement, compute angle
            // vec3 = vec1 x vec2, not scaled
            vec3[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
            vec3[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
            vec3[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

            // Use the projDir as reference direction (positive)
            angleSign = Vmath::Dot(3, vec3, 1, projDir, 1) > 0.0 ? 1.0 : -1.0;
            angleCos  = Vmath::Dot(3, vec1, 1, vec2, 1);

            // Avoid error in using acos()
            if(angleCos>1.0)
            {
                angleCos=1.0;
            }
            else if(angleCos<-1.0)
            {
                angleCos=-1.0;
            }

            angleAbs = acos(angleCos);
            angleSum += angleSign * angleAbs;
        }   
    }
    
    // Check
    angleSum = fabs(angleSum);
    if( fabs(angleSum-2.0*M_PI) < angleTol )
    {
        return true;
    }
    else
    {
        return false;
    }

}


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
bool ProcessWallNormalData::BisectionForLocCoordOnBndElmt(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & gloCoord,
    const Array<OneD, const Array<OneD, NekDouble> > & pts,
    const Array<OneD, const int > & dirUse,
    Array<OneD, NekDouble > & locCoord,
    const NekDouble iterTol,
    const int iterMax)
{
    // Initial settings
    Array<OneD, NekDouble> etaLR(2); // range [-1,1]
    etaLR[0] = -1.0;                 // left
    etaLR[1] =  1.0;                 // right
    NekDouble tmpL, tmpR;            // tmp values for L/R

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    locCoord[0] = -2.0;
    int cnt = 0;
    bool isConverge = false;
    while(cnt<iterMax)
    {
        tmpL = bndXmap->PhysEvaluate(etaLR  , pts[dirUse[0]]);
        tmpR = bndXmap->PhysEvaluate(etaLR+1, pts[dirUse[0]]);

        if (fabs(gloCoord[dirUse[0]]-tmpL) >= fabs(gloCoord[dirUse[0]]-tmpR))
        {
            etaLR[0] = 0.5 * (etaLR[0]+etaLR[1]);
        }
        else
        {
            etaLR[1] = 0.5 * (etaLR[0]+etaLR[1]);
        }

        if ((etaLR[1]-etaLR[0]) < iterTol)
        {
            locCoord[0] = 0.5 * (etaLR[0]+etaLR[1]);
            isConverge = true;
            break;
        }

        ++cnt;
    }

    // Warning if failed
    if(cnt >= iterMax)
    {
        WARNINGL1(false, "Bisection iteration is not converged");
    }

    return isConverge;
}


bool ProcessWallNormalData::NewtonIterForLocCoordOnBndElmt(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble> & gloCoord,
    const Array<OneD, const Array<OneD, NekDouble> > & pts,
    const Array<OneD, const int > & dirUse,
    Array<OneD, NekDouble> & locCoord,
    NekDouble & dist,
    const NekDouble iterTol,
    const int iterMax)
{

    const NekDouble LcoordDiv = 15.0;

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();

    Array<OneD, const NekDouble> Jac =
        bndGeom->GetMetricInfo()->GetJac(bndXmap->GetPointsKeys());
    NekDouble scaledTol = Vmath::Vsum(Jac.size(), Jac, 1) /
                          ((NekDouble)Jac.size());
    scaledTol *= iterTol;


    // Set the gloCoord used to compute locCoord
    const int dir1 = dirUse[0]; 
    const int dir2 = dirUse[1];

    Array<OneD, NekDouble> Dx1D1(pts[dir1].size());
    Array<OneD, NekDouble> Dx1D2(pts[dir1].size());
    Array<OneD, NekDouble> Dx2D1(pts[dir1].size());
    Array<OneD, NekDouble> Dx2D2(pts[dir1].size());

    // Ideally this will be stored in m_geomfactors
    bndXmap->PhysDeriv(pts[dir1], Dx1D1, Dx1D2);
    bndXmap->PhysDeriv(pts[dir2], Dx2D1, Dx2D2);
       
    // Initial the locCoord, in [-1,1]
    locCoord[0] = 0.0;
    locCoord[1] = 0.0;

    NekDouble x1map, x2map, F1, F2;
    NekDouble derx1_1, derx1_2, derx2_1, derx2_2, jac;
    NekDouble resid;
    int cnt = 0;
    bool isConverge = false;
    

    F1 = F2 = 2000; // Starting value of Function
    while (cnt++ < iterMax)
    {
        x1map = bndXmap->PhysEvaluate(locCoord, pts[dir1]);
        x2map = bndXmap->PhysEvaluate(locCoord, pts[dir2]);

        F1 = gloCoord[dir1] - x1map;
        F2 = gloCoord[dir2] - x2map;

        if (F1 * F1 + F2 * F2 < scaledTol)
        {
            resid = sqrt(F1 * F1 + F2 * F2);
            isConverge = true;
            break;
        }

        // Interpolate derivative metric at locCoord
        derx1_1 = bndXmap->PhysEvaluate(locCoord, Dx1D1);
        derx1_2 = bndXmap->PhysEvaluate(locCoord, Dx1D2);
        derx2_1 = bndXmap->PhysEvaluate(locCoord, Dx2D1);
        derx2_2 = bndXmap->PhysEvaluate(locCoord, Dx2D2);

        jac = derx2_2 * derx1_1 - derx2_1 * derx1_2;
        
        // Use analytical inverse of derivitives which are
        // also similar to those of metric factors.
        locCoord[0] = locCoord[0] + ( derx2_2 * (gloCoord[dir1] - x1map) -
                                      derx1_2 * (gloCoord[dir2] - x2map)) / jac;

        locCoord[1] = locCoord[1] + (-derx2_1 * (gloCoord[dir1] - x1map) +
                                      derx1_1 * (gloCoord[dir2] - x2map)) / jac;

        
        // locCoord have diverged so stop iteration
        if( !(std::isfinite(locCoord[0]) && std::isfinite(locCoord[1])) )
        {
            dist = 1e16;
            std::ostringstream ss;
            ss << "nan or inf found in NewtonIterForLocCoordOnProjBndElmt in element "
               << bndGeom->GetGlobalID();
            WARNINGL1(false, ss.str());
            return false;
        }
        if (fabs(locCoord[0]) > LcoordDiv || fabs(locCoord[1]) > LcoordDiv)
        {
            break; 
        }
    }

    // Check distance for collapsed coordinate 
    Array<OneD, NekDouble> eta(2);
    bndXmap->LocCoordToLocCollapsed(locCoord, eta);

    if(bndGeom->ClampLocCoords(eta, 0.0))
    {
        // calculate the global point corresponding to locCoord
        x1map = bndXmap->PhysEvaluate(eta, pts[dir1]);
        x2map = bndXmap->PhysEvaluate(eta, pts[dir2]);

        F1 = gloCoord[dir1] - x1map;
        F2 = gloCoord[dir2] - x2map;

        dist = sqrt(F1 * F1 + F2 * F2);
    }
    else
    {
        dist = 0.0;
    }

    // Warning if failed
    if (cnt >= iterMax)
    {
        Array<OneD, NekDouble> collCoords(2);
        bndXmap->LocCoordToLocCollapsed(locCoord, collCoords);

        // if coordinate is inside element dump error!
        if ((collCoords[0] >= -1.0 && collCoords[0] <= 1.0) &&
            (collCoords[1] >= -1.0 && collCoords[1] <= 1.0))
        {
            std::ostringstream ss;

            ss << "Reached MaxIterations (" << iterMax
               << ") in Newton iteration ";
            ss << "Init value (" << setprecision(4) << 0 << "," << 0
               << ","
               << ") ";
            ss << "Fin  value (" << locCoord[0] << "," << locCoord[1] << ","
               << ") ";
            ss << "Resid = " << resid << " Tolerance = " << sqrt(scaledTol);

            WARNINGL1(cnt < iterMax, ss.str());
        }
    }

    return isConverge;

}


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
 * @param maxDist      Disntance to check if the wall point is desired.
 * @param iterTol      Tolerence for iteration.
 * @return             Inside (true) or not (false)
 */
bool ProcessWallNormalData::BndElmtContainsPoint(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & gloCoord,
    const Array<OneD, const NekDouble > & projDir,
    Array< OneD, NekDouble > & locCoord,
    NekDouble & projDist,
    const NekDouble maxDist,
    const NekDouble iterTol)
{
    // Get variables
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    const int npts      = bndXmap->GetTotPoints();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);     // =2 for 2.5D cases

    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, NekDouble> >       projPts(nCoordDim);
    Array<OneD, NekDouble >                    projGloCoord(3, 0.0);
    
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i]     = Array<OneD, NekDouble>(npts);
        projPts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    // Project the point and vertices of the element in the input direction
    ProjectPoint(gloCoord, projDir, 0.0, projGloCoord);
    ProjectVertices(pts, projDir, 0.0, projPts);


    // Set the main direction(s) and the minor direction
    // The gloCoord for minor direction will not be used for locCoord iteration
    // dirUse[0] is the main dir for 2D/2.5D cases, dirUse[1]/[2] is the minor one
    // dirUse[0] and dirUse[1] are the main dir for 3D cases, dirUse[2] is hte minor one
    int dirMaxId = 0; // id to get the dir with largest projDir component
    for (int i=1; i<nCoordDim; ++i)
    {
        if (fabs(projDir[i])>fabs(projDir[dirMaxId]))
        {
            dirMaxId = i;
        }
    }

    Array<OneD, int > dirUse(3, 0); 
    if (nCoordDim==2)
    {
        // 2D or 2.5D cases
        if (dirMaxId==0)
        {
            dirUse[0] = 1;
            dirUse[1] = 0;
            dirUse[2] = 2;
        }
        else 
        {
            dirUse[0] = 0;
            dirUse[1] = 1;
            dirUse[2] = 2;
        }
    }
    else
    {
        // 3D cases
        if (dirMaxId==0)
        {
            dirUse[0] = 1;
            dirUse[1] = 2;
            dirUse[2] = 0;
        }
        else if (dirMaxId==1)
        {
            dirUse[0] = 2;
            dirUse[1] = 0;
            dirUse[2] = 1;
        }
        else
        {
            dirUse[0] = 0;
            dirUse[1] = 1;
            dirUse[2] = 2;
        }

    }


    // Check if the projected point is in the projected elmt
    // If yes, then compute the locCoord and check if desired point is found
    if(nCoordDim==2)
    {
        if (isInProjectedArea2D(projGloCoord, projPts, 1.0e-12))
        {
            bool isConverge, isDesired;
            
            isConverge = BisectionForLocCoordOnBndElmt(bndGeom, projGloCoord,
                             projPts, dirUse, locCoord, iterTol);

            Array<OneD, NekDouble > tmp(2, 0.0);
            tmp[0] = bndXmap->PhysEvaluate(locCoord, pts[0]) - gloCoord[0];
            tmp[1] = bndXmap->PhysEvaluate(locCoord, pts[1]) - gloCoord[1];
            projDist = Vmath::Dot(2, tmp, 1, projDir, 1);  // can be negative
            
            isDesired = (projDist > 0.0) && (projDist < maxDist);

            return isConverge && isDesired;
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        if (isInProjectedArea3D(projGloCoord, projPts, projDir, 1.0e-12, 1.0e-6))
        {
            NekDouble dist; 
            bool isConverge, isDesired;

            isConverge = NewtonIterForLocCoordOnBndElmt(bndGeom, projGloCoord,
                             projPts, dirUse, locCoord, dist, iterTol);
            
            if (dist>iterTol)
            {
                std::ostringstream ss;
                ss << "Collapsed locCoord out of range.\n"
                   << "Newton iteration gives the distance: " << dist;
                WARNINGL1(false, ss.str());
            }
            
            Array<OneD, NekDouble > tmp(3, 0.0);
            tmp[0] = bndXmap->PhysEvaluate(locCoord, pts[0]) - gloCoord[0];
            tmp[1] = bndXmap->PhysEvaluate(locCoord, pts[1]) - gloCoord[1];
            tmp[2] = bndXmap->PhysEvaluate(locCoord, pts[2]) - gloCoord[2];
            projDist = Vmath::Dot(3, tmp, 1, projDir, 1);  // can be negative

            isDesired = (projDist > 0.0) && (projDist < maxDist);

            return isConverge && isDesired;
        }
        else
        {
            return false;
        }
        
    }
    
}


/**
 * @brief Get the normals for a given locCoord
 * @param bndGeom      Pointer to the geometry of the boundary element.
 * @param locCoord     Iteration results for local coordinates (if inside).
 * @param normals      Wall normal as the result
 */
void ProcessWallNormalData::GetNormals(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & locCoord, 
    Array< OneD, NekDouble > & normals)
{
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);     // =2 for 2.5D cases
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    Array<OneD, Array<OneD, NekDouble> >        pts(m_spacedim);
    Array<OneD, Array<OneD, const NekDouble> >  bndCoeffs(m_spacedim);
    int npts = bndXmap->GetTotPoints();

    // Get pts in the element
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    // Get the outward-pointing normals according to the given locCoord
    if(nCoordDim==2)
    {
        Array<OneD, NekDouble> DxD1(pts[0].size());
        Array<OneD, NekDouble> DyD1(pts[0].size());

        bndXmap->PhysDeriv(pts[0], DxD1);
        bndXmap->PhysDeriv(pts[1], DyD1);

        NekDouble dxd1, dyd1, norm;
        dxd1 = bndXmap->PhysEvaluate(locCoord, DxD1);
        dyd1 = bndXmap->PhysEvaluate(locCoord, DyD1);
        norm = sqrt(dxd1*dxd1 + dyd1*dyd1);

        normals[0] =  dyd1/norm;
        normals[1] = -dxd1/norm;
        normals[2] =  0.0;
    }
    else
    {
        Array<OneD, NekDouble> DxD1(pts[0].size());
        Array<OneD, NekDouble> DxD2(pts[0].size());
        Array<OneD, NekDouble> DyD1(pts[0].size());
        Array<OneD, NekDouble> DyD2(pts[0].size());
        Array<OneD, NekDouble> DzD1(pts[0].size());
        Array<OneD, NekDouble> DzD2(pts[0].size());

        bndXmap->PhysDeriv(pts[0], DxD1, DxD2);
        bndXmap->PhysDeriv(pts[1], DyD1, DyD2);
        bndXmap->PhysDeriv(pts[2], DzD1, DzD2);

        NekDouble dxd1, dyd1, dzd1, dxd2, dyd2, dzd2;
        dxd1 = bndXmap->PhysEvaluate(locCoord, DxD1);
        dxd2 = bndXmap->PhysEvaluate(locCoord, DxD2);
        dyd1 = bndXmap->PhysEvaluate(locCoord, DyD1);
        dyd2 = bndXmap->PhysEvaluate(locCoord, DyD2);
        dzd1 = bndXmap->PhysEvaluate(locCoord, DzD1);
        dzd2 = bndXmap->PhysEvaluate(locCoord, DzD2);

        NekDouble n1, n2, n3, norm;
        n1 = dyd1*dzd2 - dyd2*dzd1;
        n2 = dzd1*dxd2 - dzd2*dxd1;
        n3 = dxd1*dyd2 - dxd2*dyd1;
        norm = sqrt(n1*n1 + n2*n2 + n3*n3);

        normals[0] = n1/norm;
        normals[1] = n2/norm;
        normals[2] = n3/norm;
    }

}



}
}
