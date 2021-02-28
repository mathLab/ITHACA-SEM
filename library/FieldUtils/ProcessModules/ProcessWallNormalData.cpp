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

#include "ProcessWallNormalData.h"

#define _DEBUG_

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
    m_config["x"]    = ConfigOption(false, "0.1", 
                       "x-coordinate of the sample origin.");
    m_config["y"]    = ConfigOption(false, "0.0",
                       "y-coordinate of the sample origin.");
    m_config["z"]    = ConfigOption(false, "0.0",
                       "z-coordinate of the sample origin.");
    m_config["tol"]  = ConfigOption(false, "0.1", 
                       "Relative tolerence to find the exact origin.");
    m_config["useY"] = ConfigOption(true, "0", 
                       "Flag to set using y-coordinate to get the origin.");
    m_config["h"]    = ConfigOption(false, "0.01", 
                       "Sampling distance along the wall normals.");
    m_config["nh"]   = ConfigOption(false, "5", 
                       "Number of sampling points along the wall normals.");
    m_config["d"]    = ConfigOption(false, "0.1", 
                       "Points distribution control in h direction, in (0,1)");
}

ProcessWallNormalData::~ProcessWallNormalData()
{
}


/*Note
* This module is used to get field data in the wall-normal direction.
* The input cases can be 2D, 2.5D and 3D. 
* The data will be exported with .pts extension.
*
* The user defined parameters are: bnd, x, y, z, tol, useY, h, nh, d
* bnd is the boundary id. This boundary should contain the desired origin.
* (x,y,z) are the coordinates of the input sampling origin. They are used as
* references to get exact origin on the wall. By default, the projection in the
* x-direction is used. If flag useY is set to be 1, y-direction projection will
* be used instead. The definations of the x/y/z-direction are as follows:
*   - The x-direction is set to be the chordwise direction
*   - The y-direction is set to be the vertical direction
*   - The z-direction is set to be the spanwise direction
* tol is the relative tolerence to find the exact origin. It is the absolute 
* tolerence over the averaged boundary element edge length
* h is the sampling depth in the wall-normal direction
* nh is the number of sampling points along h.
* d is a destribution control parameter of the sampling points. It should be in
* the range (0,1). d=0.99 gives approximately evenly spaced array.
*/
void ProcessWallNormalData::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    // Initialize the sampling parameters
    Array<OneD, NekDouble> orig(3); // gloCoord of the origin
    orig[0]                   = m_config["x"].as<NekDouble>();
    orig[1]                   = m_config["y"].as<NekDouble>();
    orig[2]                   = m_config["z"].as<NekDouble>();
    const NekDouble relTol    = m_config["tol"].as<NekDouble>();
    const bool      isUseY    = m_config["useY"].as<bool>();
    const NekDouble distanceH = m_config["h"].as<NekDouble>();
    const int       nptsH     = m_config["nh"].as<int>();
    const NekDouble delta     = m_config["d"].as<NekDouble>();  
   
    // Get dim to store data
    const int nfields   = m_f->m_variables.size();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    m_spacedim          = nCoordDim + m_f->m_numHomogeneousDir;
    const int nBndLcoordDim = 
        (m_spacedim==3 && m_f->m_numHomogeneousDir==0) ? 2 : 1;
    const int totVars   = m_spacedim + m_f->m_variables.size();

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
    StdRegions::StdExpansionSharedPtr bndXmap;
    Array<OneD, const NekDouble> Jac;
    NekDouble scaledTol;
    int npts;
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    Array<OneD, NekDouble> locCoord(nCoordDim-1, -999.0);
    NekDouble resid;
    int elmtid;
    bool isInside = false;

    // Search and get precise Lcoord
    const int nElmts = BndExp[0]->GetNumElmts();
    for (elmtid=0; elmtid<nElmts; ++elmtid) //nElmts
    {    
        bndGeom = BndExp[0]->GetExp(elmtid)->GetGeom(); 
        bndXmap = bndGeom->GetXmap();
        
        // Get the scaled tol by averaged bnd elemt edge length 
        Jac = bndGeom->GetMetricInfo()->GetJac(bndXmap->GetPointsKeys());
        scaledTol = Vmath::Vsum(Jac.size(), Jac, 1)/((NekDouble)Jac.size());
        scaledTol = pow(Jac[0], 
                        1.0/(static_cast<NekDouble>(nCoordDim)-1.0))*2.0;
        scaledTol *= relTol;

        // Get the point
        npts = bndXmap->GetTotPoints();
        for (int i=0; i<nCoordDim; ++i) 
        {
            pts[i] = Array<OneD, NekDouble>(npts);
            bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
            bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
        }

        // Key routine
        isInside = BndElmtContainsPoint(bndGeom, orig, locCoord, isUseY, scaledTol, resid); 

        if (isInside) 
        {
            break;
        }

    }
    ASSERTL0(isInside, "Failed to find the sampling origin on the boundary."); 

    // Update the precise sampling position
    // x is precise, update y, or vice versa
    if(!isUseY)
    {
        orig[1] = bndXmap->PhysEvaluate(locCoord, pts[1]);                 
    }
    else
    {
        orig[0] = bndXmap->PhysEvaluate(locCoord, pts[0]);
    }

    // Update z according to the closest plane in 2.5D cases
    if(m_f->m_numHomogeneousDir==1)
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
    int refId = elmtid*from_nPtsPerElmt; // the 1st normal in the bnd element
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
    
#ifdef _DEBUG_
    cout << "Ref normals:" << endl;
    for(int i=0;i<from_nPtsPerElmt;++i){
        cout << i << ", Q [nx,ny,nz] = [" 
             << -normalsQ[0][elmtid*from_nPtsPerElmt+i] << ", "
             << -normalsQ[1][elmtid*from_nPtsPerElmt+i] << "]" << endl;
    }
    cout << "Final normals:" << endl;
    cout << "[nx,ny,nz] = [" << normals[0] << ", " << normals[1] << ", " << normals[2] << "]" << endl;
    cout << "Final origin coordinates:" << endl;
    cout << "[Ox,Oy,Oz] = [" << orig[0] << ", " << orig[1] << ", " << orig[2] << "]" << endl;   
#endif

    
    //-------------------------------------------------------------------------
    // Set the sampling array
    // Expression in Agrawal's paper:
    // h = 1- tanh((1-ksi)*atanh(sqrt(1-delta)))/sqrt(1-delta), ksi in [0,1]
    Array<OneD, NekDouble> h(nptsH);
    NekDouble tmp1;
    const NekDouble tmp2 = 1.0/(static_cast<NekDouble>(nptsH)-1.0);
    const NekDouble tmp3 = sqrt(1.0-delta);
    const NekDouble tmp4 = atanh(tmp3);
    const NekDouble tmp5 = 1.0/tmp3;
    for (int i=0; i<nptsH; ++i)
    {
        tmp1 = 1.0 - i * tmp2; // tmp1 = 1-ksi
        h[i] = 1 - tanh(tmp1*tmp4)*tmp5;
        h[i] *= distanceH; // physical distance in normal direction
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



/*
* Check if the point is in the projected area or not
* bndGeom   geometry of single element. It is unnecessary if
            this routine is moved to Geometry class 
* gloCoord  global coordinates of the point to be checked
* projDir   projection direction, the coordinate to be computed
*/
bool ProcessWallNormalData::isInProjectedArea2D(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & gloCoord,
    const int projDir)
{
    // Get the variables 
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    const int npts      = bndXmap->GetTotPoints();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    // Set direction
    int dir1 = 0; // Main(1) and minor(2) coordinates
    if(projDir==0)
    {
        dir1 = 1;
    }
    else if(projDir==1)
    {
        dir1 = 0;
    }
    else{
        ASSERTL0(false, "Incorrect minor (projection) direction."); 
    }

    // Set range
    NekDouble valueL, valueR;
    if(pts[dir1][0] <= pts[dir1][npts-1])
    {
        valueL = pts[dir1][0];
        valueR = pts[dir1][npts-1];
    }
    else
    {
        valueL = pts[dir1][npts-1];
        valueR = pts[dir1][0];
    }

    // Check
    if(valueL <= gloCoord[dir1] && gloCoord[dir1] <= valueR)
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool ProcessWallNormalData::isInProjectedArea3D(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & gloCoord,
    const int projDir)
{
    // Get the variables 
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    const int npts      = bndXmap->GetTotPoints();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    // Set direction
    int dir1 = 2, dir2 = 0; // Main(1,2) and minor(3) coordinates
    if(projDir==0)
    {
        dir1 = 1;
        dir2 = 2;
    }
    else if(projDir==1)
    {
        dir1 = 2;
        dir2 = 0;
    }
    else if(projDir==2)
    {
        dir1 = 0;
        dir2 = 1;
    }
    else{
        ASSERTL0(false, "Incorrect minor (projection) direction."); 
    }


    //Generate a polygen (quad), re-order the points on the edge
    const int nptsEdge    = sqrt(npts);  // num of points on edge for geo representation
    const int nptsPolygon = 4*nptsEdge-4;
    Array<OneD, Array<OneD, NekDouble> > ptsPolygon(3);
    for (int i=0; i<nCoordDim; ++i) 
    {
        ptsPolygon[i] = Array<OneD, NekDouble>(nptsPolygon);

        for(int j=0; j<nptsEdge; ++j)
        {
            ptsPolygon[i][j] = pts[i][j];
        }
        for(int j=0; j<nptsEdge-2; ++j)
        {
            ptsPolygon[i][nptsEdge+j] = pts[i][(j+2)*nptsEdge-1];
        }
        for(int j=0; j<nptsEdge; ++j)
        {
            ptsPolygon[i][2*nptsEdge-2+j] = pts[i][npts-1-j];
        }
        for(int j=0; j<nptsEdge-2; ++j)
        {
            ptsPolygon[i][3*nptsEdge-2+j] = pts[i][nptsEdge*(nptsEdge-j-2)];
        }
    }
        
    
    // Determine relation using angle method (sum=2*pi)
    const NekDouble paralTol = 1.0e-12;
    const NekDouble angleTol = 1.0e-6;
    NekDouble angleCos, angleAbs, angleSign, angleSum = 0.0;
    NekDouble norm1, norm2;
    Array<OneD, NekDouble> vec1(2, 0.0), vec2(2, 0.0);
    int id1, id2;

    for(int i=0; i<nptsPolygon; ++i)
    {
        id1 = i;
        id2 = (id1==(nptsPolygon-1)) ? 0 : id1+1;

        vec1[0] = ptsPolygon[dir1][id1] - gloCoord[dir1];
        vec1[1] = ptsPolygon[dir2][id1] - gloCoord[dir2];
        vec2[0] = ptsPolygon[dir1][id2] - gloCoord[dir1];
        vec2[1] = ptsPolygon[dir2][id2] - gloCoord[dir2];

        norm1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1]);
        norm2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1]);
        vec1[0] /= norm1;
        vec1[1] /= norm1;
        vec2[0] /= norm2;
        vec2[1] /= norm2;

        if( (fabs(vec1[0]+vec2[0]) + fabs(vec1[1]+vec2[1])) < paralTol )
        {
            // On the line segement, true
            return true;
        }
        else
        {
            // Off the line segement, compute angle
            angleSign = ((vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0) ? 1.0 : -1.0;
            angleCos  = vec1[0]*vec2[0]+vec1[1]*vec2[1];

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


/*
* Compute the Lcoord using the coord
* bndGeom   geometry of single element. It is unnecessary if
            this routine is moved to Geometry class 
* gloCoord  global coordinates of the point
* pts       global coordinates of the geometry feature points
* locCoord  local coordinates to be computed
* projDir   projection direction, the coordinate to be computed
*/
bool ProcessWallNormalData::BisectionForLocCoordOnBndElmt(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble > & gloCoord,
    const Array<OneD, const Array<OneD, NekDouble> > & pts,
    const int projDir,
    Array< OneD, NekDouble > & locCoord,
    NekDouble & dist)
{
    const int MaxIterations = 51;     
    const NekDouble iterTol = 1.e-8;
  
    // Set direction
    int dir1 = 0; // Main(1) and minor(2) coordinates
    if(projDir==0)
    {
        dir1 = 1;
    }
    else if(projDir==1)
    {
        dir1 = 0;
    }
    else{
        ASSERTL0(false, "Incorrect minor (projection) direction."); 
    }

    // Bisection iteration
    Array<OneD, NekDouble> etaLR(2); // range [-1,1]
    etaLR[0] = -1.0; // left
    etaLR[1] =  1.0; // right
    NekDouble tmpL, tmpR;
    int cnt = 0;

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    locCoord[0] = -2.0;
    dist = 999.0;
    bool isConverge = false;
    while(cnt<MaxIterations)
    {
        tmpL = bndXmap->PhysEvaluate(etaLR  , pts[dir1]);
        tmpR = bndXmap->PhysEvaluate(etaLR+1, pts[dir1]);

        if (fabs(gloCoord[dir1]-tmpL) >= fabs(gloCoord[dir1]-tmpR))
        {
            etaLR[0] = 0.5 * (etaLR[0]+etaLR[1]);
        }
        else
        {
            etaLR[1] = 0.5 * (etaLR[0]+etaLR[1]);
        }

        if ( (etaLR[1]-etaLR[0]) < iterTol )
        {
            locCoord[0] = 0.5 * (etaLR[0]+etaLR[1]);
            dist = fabs(0.5*(tmpL+tmpR)-gloCoord[dir1]);
            isConverge = true;
            break;
        }

        ++cnt;
    }

    // Warning if failed
    if(cnt >= MaxIterations)
    {
        WARNINGL1(false, "Bisection iteration is not converged");
    }

    return isConverge;

}


bool ProcessWallNormalData::NewtonIterationForLocCoordOnBndElmt(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array<OneD, const NekDouble> &coords,
    const Array<OneD, const Array<OneD, NekDouble> > &pts,
    const int projDir,
    Array<OneD, NekDouble> &Lcoords,
    NekDouble &dist)
{ 
    // Maximum iterations for convergence
    const int MaxIterations = 51;
    // |x-xp|^2 < EPSILON  error    tolerance
    const NekDouble Tol = 1.e-8;
    // |r,s|    > LcoordDIV stop   the search
    const NekDouble LcoordDiv = 15.0;

    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();

    Array<OneD, const NekDouble> Jac =
        bndGeom->GetMetricInfo()->GetJac(bndXmap->GetPointsKeys());
    NekDouble ScaledTol = Vmath::Vsum(Jac.size(), Jac, 1) /
                          ((NekDouble)Jac.size());
    ScaledTol *= Tol;



    int dir1 = 2, dir2 = 0;
    if(projDir==0)
    {
        dir1 = 1;
        dir2 = 2;
    }
    else if(projDir==1)
    {
        dir1 = 2;
        dir2 = 0;   
    }
    else if(projDir==2)
    {
        dir1 = 0;
        dir2 = 1;
    }
    else{
        ASSERTL0(false, "The projection direction needs to be 0 or 1 or 2.");  
    }
    

    Array<OneD, NekDouble> Dx1D1(pts[dir1].size());
    Array<OneD, NekDouble> Dx1D2(pts[dir1].size());
    Array<OneD, NekDouble> Dx2D1(pts[dir1].size());
    Array<OneD, NekDouble> Dx2D2(pts[dir1].size());

    // Ideally this will be stored in m_geomfactors
    bndXmap->PhysDeriv(pts[dir1], Dx1D1, Dx1D2);
    bndXmap->PhysDeriv(pts[dir2], Dx2D1, Dx2D2);
       
    // Initial the locCoord
    Lcoords[0] = 0.0;
    Lcoords[1] = 0.0;

    NekDouble x1map, x2map, F1, F2;
    NekDouble derx1_1, derx1_2, derx2_1, derx2_2, jac;
    bool isConverge = false;
    int cnt = 0;

    F1 = F2 = 2000; // Starting value of Function
    NekDouble resid;
    while (cnt++ < MaxIterations)
    {
        x1map = bndXmap->PhysEvaluate(Lcoords, pts[dir1]);
        x2map = bndXmap->PhysEvaluate(Lcoords, pts[dir2]);

        F1 = coords[dir1] - x1map;
        F2 = coords[dir2] - x2map;

        if (F1 * F1 + F2 * F2 < ScaledTol)
        {
            resid = sqrt(F1 * F1 + F2 * F2);
            isConverge = true;
            break;
        }

        // Interpolate derivative metric at Lcoords
        derx1_1 = bndXmap->PhysEvaluate(Lcoords, Dx1D1);
        derx1_2 = bndXmap->PhysEvaluate(Lcoords, Dx1D2);
        derx2_1 = bndXmap->PhysEvaluate(Lcoords, Dx2D1);
        derx2_2 = bndXmap->PhysEvaluate(Lcoords, Dx2D2);

        jac = derx2_2 * derx1_1 - derx2_1 * derx1_2;
        
        // use analytical inverse of derivitives which are
        // also similar to those of metric factors.
        Lcoords[0] =
            Lcoords[0] +
            ( derx2_2 * (coords[dir1] - x1map) - derx1_2 * (coords[dir2] - x2map)) / jac;

        Lcoords[1] =
            Lcoords[1] +
            (-derx2_1 * (coords[dir1] - x1map) + derx1_1 * (coords[dir2] - x2map)) / jac;

        
        if( !(std::isfinite(Lcoords[0]) && std::isfinite(Lcoords[1])) )
        {
            dist = 1e16;
            std::ostringstream ss;
            ss << "nan or inf found in NewtonIterationForLocCoord in element "
               << bndGeom->GetGlobalID();
            WARNINGL1(false, ss.str());
            return false;
        }
        if (fabs(Lcoords[0]) > LcoordDiv || fabs(Lcoords[1]) > LcoordDiv)
        {
            break; // lcoords have diverged so stop iteration
        }
    }

    Array<OneD, DNekMatSharedPtr> I(2);
    Array<OneD, NekDouble> eta(2);

    bndXmap->LocCoordToLocCollapsed(Lcoords, eta);
    if(bndGeom->ClampLocCoords(eta, 0.0))
    {
        I[0] = bndXmap->GetBasis(0)->GetI(eta);
        I[1] = bndXmap->GetBasis(1)->GetI(eta + 1);
        // calculate the global point corresponding to Lcoords
        x1map = bndXmap->PhysEvaluate(I, pts[dir1]);
        x2map = bndXmap->PhysEvaluate(I, pts[dir2]);
        F1 = coords[dir1] - x1map;
        F2 = coords[dir2] - x2map;
        dist = sqrt(F1 * F1 + F2 * F2);
    }
    else
    {
        dist = 0.0;
    }

    // Warning if failed
    if (cnt >= MaxIterations)
    {
        Array<OneD, NekDouble> collCoords(2);
        bndXmap->LocCoordToLocCollapsed(Lcoords, collCoords);

        // if coordinate is inside element dump error!
        if ((collCoords[0] >= -1.0 && collCoords[0] <= 1.0) &&
            (collCoords[1] >= -1.0 && collCoords[1] <= 1.0))
        {
            std::ostringstream ss;

            ss << "Reached MaxIterations (" << MaxIterations
               << ") in Newton iteration ";
            ss << "Init value (" << setprecision(4) << 0 << "," << 0
               << ","
               << ") ";
            ss << "Fin  value (" << Lcoords[0] << "," << Lcoords[1] << ","
               << ") ";
            ss << "Resid = " << resid << " Tolerance = " << sqrt(ScaledTol);

            WARNINGL1(cnt < MaxIterations, ss.str());
        }
    }

    return isConverge;

}


/*
* Check if the boundary element contain an input point
* Find the Lcoord if it is inside
*/
bool ProcessWallNormalData::BndElmtContainsPoint(
    SpatialDomains::GeometrySharedPtr bndGeom,
    const Array< OneD, const NekDouble > & gloCoord,
    Array< OneD, NekDouble > & locCoord,
    const bool isUseY, 
    const NekDouble geomTol,
    NekDouble & dist)
{
    // Get variables
    StdRegions::StdExpansionSharedPtr bndXmap = bndGeom->GetXmap();
    const int npts      = bndXmap->GetTotPoints();
    const int nCoordDim = m_f->m_exp[0]->GetCoordim(0);     // =2 for 2.5D cases
    Array<OneD, Array<OneD, NekDouble> >       pts(nCoordDim);
    Array<OneD, Array<OneD, const NekDouble> > bndCoeffs(nCoordDim);
    
    for (int i=0; i<nCoordDim; ++i) 
    {
        pts[i] = Array<OneD, NekDouble>(npts);
        bndCoeffs[i] = bndGeom->GetCoeffs(i); // 0/1/2 for x/y/z
        bndXmap->BwdTrans(bndCoeffs[i], pts[i]);
    }

    const int dir1 = isUseY ? 1 : 0; // Use x or y to determine the location
    const int dir2 = 1 - dir1;       // The other Lcoord to be computed


    if(nCoordDim==2)
    {
        if( isInProjectedArea2D(bndGeom, gloCoord, dir2) )
        {
            NekDouble tmp;
            bool isConverge, isDesired;
            
            isConverge = BisectionForLocCoordOnBndElmt(bndGeom, gloCoord, pts, dir2, locCoord, dist);

            tmp = bndXmap->PhysEvaluate(locCoord, pts[dir2]);
            if (fabs(gloCoord[dir2]-tmp) < geomTol)
            {
                isDesired = true;
            }

            return isConverge && isDesired;

         }
         else{
             return false;
         }
    }
    else
    {
        if( isInProjectedArea3D(bndGeom, gloCoord, dir2) )
         {
            NekDouble tmp;
            bool isConverge, isDesired;
            
            isConverge = NewtonIterationForLocCoordOnBndElmt(bndGeom, gloCoord, pts, dir2, locCoord, dist);

            tmp = bndXmap->PhysEvaluate(locCoord, pts[dir2]);
            if (fabs(gloCoord[dir2]-tmp) < geomTol)
            {
                isDesired = true;
            }

            return isConverge && isDesired;
            
         }
         else{
             return false;
         }
    }

}


/*
* Get the normals for a given locCoord
*/
void ProcessWallNormalData::GetNormals(
    SpatialDomains::GeometrySharedPtr bndGeom,
    Array< OneD, NekDouble > & locCoord, 
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
