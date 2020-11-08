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
//  Description: Export data in the wall normal direction along the surface.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include "ProcessWallNormalData.h"

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessWallNormalData::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "wnd"),
    ProcessWallNormalData::create,
    "Export data in the wall normal direction along the surface.");

ProcessWallNormalData::ProcessWallNormalData(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
    f->m_writeBndFld = false; // turned on in upstream ProcessBoundaryExtract

    m_config["h"]  = ConfigOption(false, "0.01", 
                     "Sampling distance along the wall normals.");
    m_config["nh"] = ConfigOption(false, "3", 
                     "Number of sampling points along the wall normals.");
    m_config["ne"] = ConfigOption(false, "3", 
                     "Number of sampling points on each edge of an element.");
    m_config["d"]  = ConfigOption(false, "0.1", 
                     "Points distribution control in h direction, in (0,1)");
}

ProcessWallNormalData::~ProcessWallNormalData()
{
}


/* To heapify (max-haeading) a subtree rooted with node rootId
   A[0~5][0~currentLength-1] is the data to be heapified
   A.size()=6; A[0~5] for x/y/z/nx/ny/nz respectively
   The heap is adjusted regards to A[0] (x-value)
*/
void ProcessWallNormalData::Heapify_max(Array<OneD, Array<OneD, NekDouble> > A, 
                                 const int curLen, 
                                 const int rootId)
{
    const int dataDim = A.size();
    int maxId = rootId;   // Initialize the maximum as root 
    const int leftId  = 2 * rootId + 1; // left  child = 2*i + 1 
    const int rightId = 2 * rootId + 2; // right child = 2*i + 2 

    // Check if left child exists and it is larger than root
    // Then check if right child exits and it is larger than largest  
    if (leftId  < curLen && A[0][leftId]  > A[0][maxId]) { maxId = leftId;  } 
    if (rightId < curLen && A[0][rightId] > A[0][maxId]) { maxId = rightId; }

    // If largest is not the root, swap values at [maxId] and [rootId]
    // then recursively heapify the affected sub-tree, rooted at [maxId] 
    if (maxId != rootId)
    {
        for (int j=0; j<dataDim; ++j) { std::swap(A[j][rootId], A[j][maxId]); }
        Heapify_max(A, curLen, maxId);
    }

}

/* Sort the array using heap*/
void ProcessWallNormalData::HeapSort(Array<OneD, Array<OneD, NekDouble> > A)
{
    const int dataDim = A.size();
    const int totLen  = A[0].size();
    
    // Build max heap, starting from the last non-leaf node
    for (int i = floor(totLen / 2) - 1; i >= 0; --i)
    {
        Heapify_max(A, totLen, i);
    }

    // Move current root to end [0]->[currentlength-1]
    // and adjust the reduced heap, startig from root
    for (int curLen = totLen; curLen > 1; --curLen)
    {
        for (int j=0; j<dataDim; ++j) { std::swap(A[j][0], A[j][curLen-1]); }

        Heapify_max(A, curLen-1, 0);
    }

}

/* Clean the repeated points in the array
   put the repeated points in the end 
   return the non-repeated array length
*/
int ProcessWallNormalData::CleanRepeatedPts(Array<OneD, Array<OneD, NekDouble> > A)
{
    const int dataDim = A.size();
    const int totLen = A[0].size();
    int       newLen = totLen; // initialize the new length as the total length
    for (int i=1; i<newLen; ++i) {
        // For each i, check indax smaller than i
        for (int j=i-1; j>=0; --j) {

            if (abs(A[0][i]-A[0][j]) < NekConstants::kNekZeroTol &&
                abs(A[1][i]-A[1][j]) < NekConstants::kNekZeroTol &&
                abs(A[2][i]-A[2][j]) < NekConstants::kNekZeroTol)
            {
                // repeated point found, [i]==[j]
                // move element [i+1] to [totLen-1] forword
                for (int k=0; k<dataDim; ++k)
                {
                    for (int t=i+1; t<totLen; ++t)
                    {
                        A[k][t-1] = A[k][t];
                    }
                    A[k][totLen-1] = A[k][j];
                }

                --newLen; // key origins -- 
                --i;      // check the same i again
            }

        }
    }

    return (newLen);
}


void ProcessWallNormalData::WriteDataInPts(const std::string &outFile, 
    const Array<OneD, Array<OneD, Array<OneD, NekDouble> > > data, 
    const int len)
{

    // Step 1 - settings
    bool is2D;  // true=2D, false=3D
    bool isInc; // true=incompressible. false=compressible
    int nfields = m_f->m_variables.size();

    if      (boost::iequals(m_f->m_variables[0], "u"  )) { isInc = true;  }
    else if (boost::iequals(m_f->m_variables[0], "rho")) { isInc = false; }
    else { cout << "Other types of field. Might be an issue." << endl; }
 
    if (isInc==true)
    {
        if      (nfields==3) { is2D = true;  }
        else if (nfields==4) { is2D = false; }
        else  { cout << "Incorrect dimension." << endl; }
    }
    else
    {
        if      (nfields==4) { is2D = true; }
        else if (nfields==5) { is2D = false; }
        else  { cout << "Incorrect dimension." << endl; }
    }


    // Step 2 - outpput
    std::ofstream ptsFile;
    ptsFile.open(outFile.c_str());
    ptsFile << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << endl;
    ptsFile << "<NEKTAR>" << endl;
    ptsFile << "  <POINTS ";
    ptsFile << "DIM=\"" << (is2D ? 2 : 3) << "\" ";
    ptsFile << "FIELDS=\"";
    for (int i=0; i<nfields; ++i)
    {
        ptsFile << m_f->m_variables[i];
        if (i < nfields-1)
        {
            ptsFile << ",";
        }
    }
    ptsFile << "\">" << endl;

    for (size_t i = 0; i<len; ++i)
    {
        for (size_t j=0; j<data[0].size(); ++j)
        {
            ptsFile << "    ";
            ptsFile << data[i][j][0] << " " << data[i][j][1] << " ";
            if (is2D==false)
            { 
                ptsFile << data[i][j][2] << " ";
            }
            for (int k=0; k<nfields; ++k)
            {
                ptsFile << data[i][j][3+k];
                if (k < nfields-1) {
                    ptsFile << " "; 
                }
            }
            ptsFile << endl;
        }
    }

    ptsFile << "  </POINTS>" << endl;
    ptsFile << "</NEKTAR>" << endl;
    ptsFile.close();
}


void ProcessWallNormalData::CreateFieldPts(
    const Array<OneD, Array<OneD, Array<OneD, NekDouble> > > data, 
    const int len)
{
    
    int npts_h  = data[0].size(); // number of points in normal direction
    int npts    = len * npts_h;   // len is number of sampling origins
    int totvars = m_spacedim + m_f->m_variables.size();
    Array<OneD, Array<OneD, NekDouble> > pts(totvars);

    for (int i = 0; i < totvars; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(npts);
    }

    int idx_o, idx_n, idx_v;
    for (int i = 0; i < npts; ++i)
    {
        idx_n = i%npts_h;         // idx for normal direction
        idx_o = (i-idx_n)/npts_h; // idx for origin
        idx_v = 0;                // idx for variable
        for (int j = 0; j < totvars; ++j)
        {
            pts[j][i] = data[idx_o][idx_n][idx_v];
            ++idx_v;
            if(m_spacedim==2 && j==1){++idx_v;}
        }
    }

    m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(m_spacedim, m_f->m_variables, pts, LibUtilities::NullPtsInfoMap);

}


// $NEK/build-g/dist/bin/FieldConvert-g -m wnd:bnd=0:h=0.005:nh=21:ne=11 mesh.xml session.xml field.fld/ data.pts

void ProcessWallNormalData::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    // Step 1 - Get sampling parameters
    const NekDouble distance_h = m_config["h"].as<NekDouble>(); //0.005
    const int       npts_h     = m_config["nh"].as<int>();      //21
    const int to_nPtsPerEdge   = m_config["ne"].as<int>();      //6
    const NekDouble delta      = m_config["d"].as<NekDouble>(); //0.1


    // Step 2 - Get data dimension
    int nfields   = m_f->m_variables.size();
    int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    m_spacedim    = nCoordDim + m_f->m_numHomogeneousDir;


    // Step 3 - Get boundary info
    // Create map of boundary ids for partitioned domains
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

    // Get boundary id
    // m_f->m_bndRegionsToWrite.size() is the number of input bnd
    // eg. =3 if bnd=0,1,2; =1 if bnd=0
    // ref. ProcedureWSS.cpp
    int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[0]];
    
    // Get expansion list for boundary and the number of points
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(nfields); 
    for (int i = 0; i < nfields; ++i) {
        BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(bnd);
    }
    int nqb = BndExp[0]->GetTotPoints(); // points for all HomModesZ planes
    
    // Get inward-pointing wall-normal vectors for all points on bnd
    Array<OneD, Array<OneD, NekDouble> > normals; 
    m_f->m_exp[0]->GetBoundaryNormals(bnd, normals);
    for (int i = 0; i < m_spacedim; ++i) {
        Vmath::Neg(nqb, normals[i], 1);
    }
 
    /*// Get coordinates for all points on the bnd (not used, just to check)
    Array<OneD, Array<OneD, NekDouble> > xyz_bnd(3);
    for (int i=0; i<3; ++i) {
        xyz_bnd[i] = Array<OneD, NekDouble>(nqb, 0.0);
    }
    BndExp[0]->GetCoords(xyz_bnd[0], xyz_bnd[1], xyz_bnd[2]);
    */


    // Step 4 - set sampling origins on the wall
    // dimension for parametric coordinate system and 
    // dimension for physical space that this cooridinate located in
    // eg. dim_para=1 dim_phys=2 for 2.5D case
    const int dim_para = BndExp[0]->GetExp(0)->GetNumBases(); 
    const int dim_phys = nCoordDim; 

    // set point keys
    LibUtilities::PointsType to_pointstype = 
                                 LibUtilities::PointsType::ePolyEvenlySpaced;
    LibUtilities::PointsKey  to_key(to_nPtsPerEdge, to_pointstype);
    Array<OneD, LibUtilities::PointsKey> from_key(dim_para);
    int from_nPtsPerElmt = 1;
    int to_nPtsPerElmt   = 1; 
    for (int i=0; i<dim_para; ++i)
    {
        from_key[i] = BndExp[0]->GetExp(0)->GetBasis(i)->GetPointsKey();
        from_nPtsPerElmt *= from_key[i].GetNumPoints();
        to_nPtsPerElmt   *= to_nPtsPerEdge;
    }

    // declare arrays to save points
    Array<OneD, Array<OneD, NekDouble> > from_ptsInElmt(3);
    Array<OneD, Array<OneD, NekDouble> > to_ptsInElmt(3);
    Array<OneD, Array<OneD, NekDouble> > to_normalsInElmt(3);
    for (int i=0; i<3; ++i)
    {
        from_ptsInElmt[i]   = Array<OneD, NekDouble>(from_nPtsPerElmt, 0.0);
        to_ptsInElmt[i]     = Array<OneD, NekDouble>(to_nPtsPerElmt, 0.0);
        to_normalsInElmt[i] = Array<OneD, NekDouble>(to_nPtsPerElmt, 0.0);
    }

    const int nElmts = BndExp[0]->GetNumElmts(); //42
    const int nOrigs = to_nPtsPerElmt * nElmts;  // num of origins on the wall
    Array<OneD, Array<OneD, NekDouble> > origs(6); // 6 for x/y/z/nx/ny/nz
    for (int i=0; i<6; ++i) {
        origs[i] = Array<OneD, NekDouble>(nOrigs, 0.0); 
    }
    
    // loop the element on the bnd
    int ptr = 0; // pointed to head to copy the arrays
    for (int i=0; i<nElmts; ++i)
    {
        // obtain the points in the element
        BndExp[0]->GetExp(i)->GetCoords(from_ptsInElmt[0], 
                                        from_ptsInElmt[1], 
                                        from_ptsInElmt[2]);
        // interp x/y/z and nx/ny/nz
        // dim_phys determins the times (xy/xyz) to loop
        // dim_para determins the interpolation function (Interp1D/Interp2D) to ues
        // ref: Expansion::v_GetCoords in Expansion.cpp
        for (int j=0; j<dim_phys; ++j)
        {
            //LibUtilities::Interp1D(from_key, &xyz_bnd[j][i*from_nPtsPerElmt], to_key, &to_ptsInElmt[j][0]); //alternative code
            //LibUtilities::Interp1D(from_key, &from_ptsInElmt[j][0], to_key, &to_ptsInElmt[j][0]); //x/y/z
            //LibUtilities::Interp1D(from_key, &normals[j][i*from_nPtsPerElmt], to_key, &to_normalsInElmt[j][0]);
            switch (dim_para)
            {
                case 1:
                {
                    // ksi - line segement element on bnd
                    LibUtilities::Interp1D(
                        from_key[0], &from_ptsInElmt[j][0], 
                        to_key,      &to_ptsInElmt[j][0]); //x/y/z
                    LibUtilities::Interp1D(
                        from_key[0], &normals[j][i*from_nPtsPerElmt], 
                        to_key,      &to_normalsInElmt[j][0]); // nx/ny/nz
                    break;
                }
                case 2:
                {
                    // ksi and eta - tri or quad element on bnd
                    LibUtilities::Interp2D(
                        from_key[0], from_key[1], 
                        &from_ptsInElmt[j][0], 
                        to_key,      to_key,      
                        &to_ptsInElmt[j][0]); //x/y/z
                    LibUtilities::Interp2D(
                        from_key[0], from_key[1],
                        &normals[j][i*from_nPtsPerElmt], 
                        to_key,      to_key,      
                        &to_normalsInElmt[j][0]);
                    break;
                }
                default:
                {
                    ASSERTL0(false, "Dimension error.");
                }
            }

            // save the interpolated results
            Vmath::Vcopy(to_nPtsPerElmt, &to_ptsInElmt[j][0],
                         1, &origs[j][ptr],   1); // copy coordinates
            Vmath::Vcopy(to_nPtsPerElmt, &to_normalsInElmt[j][0],
                         1, &origs[j+3][ptr], 1); // copy normals
        }
        ptr += to_nPtsPerElmt;

    }
 
    // sort array regards to x-value and remove repeated origin points
    HeapSort(origs);
    int nOrigs_diff = CleanRepeatedPts(origs); // number of different origins
     

    // Step 5 - set parametric distance of sampling points
    // Expression in Agrawal's paper:
    // h = 1- tanh((1-ksi)*atanh(sqrt(1-delta)))/sqrt(1-delta), ksi in [0,1]
    Array<OneD, NekDouble> h(npts_h);
    NekDouble tmp1;
    const NekDouble tmp2 = 1.0/(static_cast<NekDouble>(npts_h)-1.0);
    const NekDouble tmp3 = sqrt(1.0-delta);
    const NekDouble tmp4 = atanh(tmp3);
    const NekDouble tmp5 = 1.0/tmp3;
    for (int i=0; i<npts_h; ++i)
    {
        tmp1 = 1.0 - i * tmp2; // tmp1 = 1-ksi
        h[i] = 1 - tanh(tmp1*tmp4)*tmp5;
    }


    //=========================================================================
    // Set sampling array and interpolate
    // if 1 ---- final format
    // if 0 ---- just to debug

#if 1   
    // Step 6 - create fieldPts and interpolate
    int npts    = nOrigs_diff * npts_h;   // total number of points
    int totVars = m_spacedim  + m_f->m_variables.size();

    Array<OneD, Array<OneD, NekDouble> > pts(totVars);
    for (int i=0; i<totVars; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(npts, 0.0);
    }

    int idx=0;
    for (int i=0; i<nOrigs_diff; ++i)
    {
        for (int j=0; j<npts_h; ++j)
        {
            idx = i*npts_h + j;     // 2nd index in pts
            tmp1 = distance_h*h[j]; // physical distance in normal direction

            for (int k=0; k<m_spacedim; ++k)
            {
                pts[k][idx] = origs[k][i] + tmp1*origs[k+3][i]; // x0+dist*nx 
            }
        }
    }

    m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
        m_spacedim, m_f->m_variables, pts, LibUtilities::NullPtsInfoMap);

    Interpolator interp;
    interp.Interpolate(m_f->m_exp, m_f->m_fieldPts, NekConstants::kNekUnsetDouble);
    
    //=========================================================================
#else
    // Alternative code for Step 6 (last part)
    // declare the data array and fill in the coordinates
    // data[originId][normalId][variableId]
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > data(nOrigs_diff);
    for (int i=0; i<nOrigs_diff; ++i) {
        data[i] = Array<OneD, Array<OneD, NekDouble> >(npts_h);
        for (int j=0; j<npts_h; ++j) {
            data[i][j] = Array<OneD, NekDouble>(3+BndExp.size(), 0.0); //x/y/z+flow variables
            
            tmp1 = distance_h*h[j]; // physical distance in normal direction
            for (int k=0; k<3; ++k) {
                data[i][j][k] = origs[k][i] + tmp1*origs[k+3][i]; // x0+dist*nx
            }
        }
    }


    // Step 6 - interpolate the variables for each point
    Array<OneD, NekDouble> Lcoords(nCoordDim, 0.0); // 2.5D=2, 3D=3
    Array<OneD, NekDouble> coords(3);               // directly use 3
   
    for (int i=0; i<nOrigs_diff; ++i) {
        for (int j=0; j<npts_h; ++j) {

            // Get donor element and local coordinates
            Vmath::Vcopy(3, &data[i][j][0], 1, &coords[0], 1);
            int elmtid = -1;
            elmtid = m_f->m_exp[0]->GetExpIndex(coords, Lcoords, 
                     NekConstants::kGeomFactorsTol, false, elmtid); 

            // Homogeneous case, need to find the right plane
            int targetPlane = -1;
            if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D) {
                int nPlanes    = m_f->m_exp[0]->GetHomogeneousBasis()->GetZ().size();
                NekDouble lHom = m_f->m_exp[0]->GetHomoLen();
                targetPlane = std::round((coords[2]*nPlanes)/lHom);

                // Reset from last plane to plane 0 (periodic bnd)
                if(targetPlane == nPlanes) {targetPlane = 0;}
            }

            // limit Lcoords to [-1,1]
            for (int k = 0; k < nCoordDim; ++k) {
                Lcoords[k] = std::max(Lcoords[k], -1.0);
                Lcoords[k] = std::min(Lcoords[k],  1.0);
            }

            //-------------
            // interpolate the value for each field
            // interpolation routine already there 
            int offset;
            NekDouble value;
            if (elmtid >= 0) {
                // Get offset
                offset = m_f->m_exp[0]->GetPhys_Offset(elmtid); 
        
                // interpolate each field
                for (int f = 0; f < m_f->m_exp.size(); ++f) {
                    // interpolate a field
                    if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D){
                        // 2.5D case, interpolate on the target plane
                        auto planeExp = m_f->m_exp[f]->GetPlane(targetPlane);
                        value         = planeExp->GetExp(elmtid)->StdPhysEvaluate(
                                        Lcoords, planeExp->GetPhys() + offset);
                    }
                    else {
                        // 2D/3D and other cases
                        value = m_f->m_exp[f]->GetExp(elmtid)->StdPhysEvaluate(
                                Lcoords, m_f->m_exp[f]->GetPhys() + offset);
                    }

                    // Check and save
                    if ((boost::math::isnan)(value)) {
                        ASSERTL0(false, "NaN for interpolation.");
                    }    
                    else {
                        data[i][j][3+f] = value;
                    }

                } // loop f for diferent fields
            }
            else {
                ASSERTL0(false, "Incorrect Id for donor element.");
            }

        } // loop j for normal array
    } // loop i for orgins 
    
    // Export data in pts format directly in the current module
    // can be used to check the output file by OutputPts module
    WriteDataInPts("outTest.pts", data, nOrigs_diff);

    // create m_fieldPts for OutputPts module to write
    //CreateFieldPts(data, nOrigs_diff);

    // interpolate the pts object
    //Interpolator interp;
    //interp.Interpolate(m_f->m_exp, m_f->m_fieldPts, -99999.0);
#endif   
    //========================================================================= 

}



}
}
