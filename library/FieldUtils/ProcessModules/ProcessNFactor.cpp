////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessNFactor.cpp
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

#include "ProcessNFactor.h"

#include <LibUtilities/Foundations/Interp.h>


#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessNFactor::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "nf"),
    ProcessNFactor::create,
    "Export data in the wall normal direction along the surface.");

ProcessNFactor::ProcessNFactor(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
}

ProcessNFactor::~ProcessNFactor()
{
}


/* To heapify (max-haeading) a subtree rooted with node i
   A[0~5][0~currentLength-1] is the data to be heapified
   A.size()=6; A[0~5] for x/y/z/nx/ny/nz respectively
   The heap is adjusted regards to A[0] (x-value)
*/
void heapify_max(Array<OneD, Array<OneD, NekDouble> > A, 
                 const int currentLength, 
                 const int i) {
    
    int largest = i;   // Initialize the largest as root 
    const int l = 2 * i + 1; // left  child = 2*i + 1 
    const int r = 2 * i + 2; // right child = 2*i + 2 

    // Check if left child exists and it is larger than root
    // Then check if right child exits and it is larger than largest  
    if (l < currentLength && A[0][l] > A[0][largest]) { largest = l; } 
    if (r < currentLength && A[0][r] > A[0][largest]) { largest = r; }

    // If largest is not the root, swap values at [largest] and [i]
    // then recursively heapify the affected sub-tree 
    if (largest != i) {
        for (int j=0; j<6; ++j) { std::swap(A[j][i], A[j][largest]); }

        heapify_max(A, currentLength, largest);
    }
}




void ProcessNFactor::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);


    // Input paramaters (move to other routines later)
    

    int i;
    int nfields = m_f->m_variables.size();
    int expdim  = m_f->m_graph->GetSpaceDimension();
    m_spacedim  = expdim + m_f->m_numHomogeneousDir;

    std::cout<<"Inside the N-factor module!" << std::endl;
    std::cout<< nfields<< ", " << expdim << ", " << m_spacedim <<std::endl;
    std::cout<< m_f->m_exp[0]->GetNumElmts() <<std::endl;

    // Declare arrays
    // This part needs to be updated to suitable for compressible cases. [!]
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(m_spacedim + 1); //uvw+p
    Array<OneD, MultiRegions::ExpListSharedPtr> BndElmtExp(nfields);


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

    // m_f->m_bndRegionsToWrite.size() is the number of input bnd
    // eg. =3 if bnd=0,1,2; =1 if bnd=0
    int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[0]];
    cout << "bnd = " << bnd << endl;

    // Get expansion list for boundary and for elements containing this
    // bnd
    for (i = 0; i < (m_spacedim + 1); ++i) {
        BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(bnd);
    }
    for (i = 0; i < nfields; i++) {
        m_f->m_exp[i]->GetBndElmtExpansion(bnd, BndElmtExp[i]);
    }


    // Get number of points in expansions
    // For 2.5D cases, nqb = numElem * numPoints * HomModesZ
    // where numPoints = numModes + 1 = P + 2 by default
    // and nqe = nqb * numPoints. Why is this?
    int nqb = BndExp[0]->GetTotPoints();
    int nqe = BndElmtExp[0]->GetTotPoints(); // seems to be not used [!]

    // Get inward-pointing wall-normal vectors 
    Array<OneD, Array<OneD, NekDouble> > normals; 
    m_f->m_exp[0]->GetBoundaryNormals(bnd, normals);
    // Reverse normals, to get correct orientation for the body
    // normals[i][j], where i is direction (x/y/z) varying from 0 to 2;
    // j is the point varying from 0 to (nqb-1)
    for (i = 0; i < m_spacedim; ++i) {
        Vmath::Neg(nqb, normals[i], 1);
    }
    
    cout <<"normals1 "<< normals[1][nqb-2]<<" "<< normals[1][nqb-1] << endl; 
    cout << "nqb = " << nqb << ", nqe = " << nqe <<endl;


    Array<OneD, Array<OneD, NekDouble> > xyz_bnd(3);
    for (int i=0; i<3; ++i) {
        xyz_bnd[i] = Array<OneD, NekDouble>(nqb, 0.0);
    }
    BndExp[0]->GetCoords(xyz_bnd[0],xyz_bnd[1],xyz_bnd[2]);

    for (int i=0;i<nqb/4;++i){   // 0 ~ nqb/4, where 4 if for HomModesZ=4
        cout << i << " - " <<xyz_bnd[0][i] <<", "<<xyz_bnd[1][i]<<", "<<xyz_bnd[2][i]<<endl;
    }


    //=========================================================================
    //------------------------------ Set origins ------------------------------
    // set dimensions
    const int dim_para = BndExp[0]->GetExp(0)->GetNumBases(); // dimension for parametric coordinate system, eg. =1
    const int dim_phys = BndExp[0]->GetCoordim(0); // dimension for the physical space that the parametric coordinate system located on, eg. =2
    cout << "dim_para = " << dim_para <<", dim_coor = " << dim_phys <<endl;

    // input parameters
    // use prefix const later
    Array<OneD, NekDouble> range_x(2); //output range, 0-lower and 1-upper limit
    range_x[0] = 0.189;
    range_x[1] = 0.401;

    // set point key
    const int to_nPtsPerElmt = 6; // number of points per element, key parameter, better to be odd
    LibUtilities::PointsKey from_key = BndExp[0]->GetExp(0)->GetBasis(0)->GetPointsKey();
    LibUtilities::PointsKey to_key( to_nPtsPerElmt, LibUtilities::PointsType::ePolyEvenlySpaced ); //[!] important!
    const int from_nPtsPerElmt = from_key.GetNumPoints();

    // declare arrays to save points
    Array<OneD, Array<OneD, NekDouble> > from_ptsInElmt(3); // 3 for 3D,//offset=0,8,16,...,328
    Array<OneD, Array<OneD, NekDouble> > to_ptsInElmt(3);
    Array<OneD, Array<OneD, NekDouble> > to_normalsInElmt(3);
    for (int i=0; i<3; ++i) {
        from_ptsInElmt[i]   = Array<OneD, NekDouble>(from_nPtsPerElmt, 0.0);
        to_ptsInElmt[i]     = Array<OneD, NekDouble>(to_nPtsPerElmt, 0.0);
        to_normalsInElmt[i] = Array<OneD, NekDouble>(to_nPtsPerElmt, 0.0);
    }

    const int nElmts = BndExp[0]->GetNumElmts(); //42
    const int nOrigs = to_nPtsPerElmt * nElmts;
    Array<OneD, Array<OneD, NekDouble> > origs(6); // samping origins (have same points), 6 for x/y/z/nx/ny/nz
    for (int i=0; i<6; ++i) {
        origs[i] = Array<OneD, NekDouble>(nOrigs, 0.0); 
    }
    
    int ptr = 0;
    // loop the element on the bnd
    for ( int i = 0; i < nElmts; ++i ) { //i < nElmts

        // obtain the points in the element
        BndExp[0]->GetExp(i)->GetCoords( from_ptsInElmt[0], from_ptsInElmt[1], from_ptsInElmt[2] ); 

        // skip some elements, needs further improved
        if (from_ptsInElmt[0][0]<range_x[0] || from_ptsInElmt[0][from_nPtsPerElmt-1]>range_x[1]) { continue; } 

        // interp x/y/z and nx/ny/nz
        // needs to be further improved for cases with different dimensions
        // dim_phys determins times (xy/xyz) to loop
        // dim_para determins functions (Interp1D/Interp2D) to ues
        // ref: Expansion::v_GetCoords in Expansion.cpp
        for (int j = 0; j < dim_phys; ++j ) {
            LibUtilities::Interp1D(from_key, &from_ptsInElmt[j][0], to_key, &to_ptsInElmt[j][0]); //x/y/z
            //LibUtilities::Interp1D(from_key, &xyz_bnd[j][i*from_nPtsPerElmt], to_key, &to_ptsInElmt[j][0]); //alternative code
            LibUtilities::Interp1D(from_key, &normals[j][i*from_nPtsPerElmt], to_key, &to_normalsInElmt[j][0]);

            // save the interpolated results
            Vmath::Vcopy( to_nPtsPerElmt, &to_ptsInElmt[j][0],     1, &origs[j][ptr],   1); // copy coordinates
            Vmath::Vcopy( to_nPtsPerElmt, &to_normalsInElmt[j][0], 1, &origs[j+3][ptr], 1); // copy coordinates
        }
        ptr = ptr + to_nPtsPerElmt;

    }
 

    // build max heap and sort
    // Build max heap, starting from the last non-leaf node
    for (int i = floor(nOrigs / 2) - 1; i >= 0; --i) {
        heapify_max(origs, nOrigs, i);
    }
    // Move current root to end [0]->[length-1]
    // and adjust the reduced heap, startig from root
    for (int length = nOrigs; length > 1; --length) {
        for (int j=0; j<6; ++j) { std::swap(origs[j][0], origs[j][length-1]); }

        heapify_max(origs, length-1, 0);
    }

    // remove repeated origin points 
    int nOrigs_new = nOrigs; // new length
    for (int i=1; i<nOrigs; ++i) {
        // For each i, check indax smaller than i
        for (int j=i-1; j>=0; --j) {
            
            if ( abs(origs[0][i]-origs[0][j])>NekConstants::kNekZeroTol ) {
                break; // The array has already sorted regards to x
            }
            else {
                if (abs(origs[1][i]-origs[1][i-1])<NekConstants::kNekZeroTol &&
                    abs(origs[2][i]-origs[2][i-1])<NekConstants::kNekZeroTol) {
                    
                    // repeated points found, [i]==[j]
                    // move [i+1] to [nOrigs-1] forword
                    for (int k=0; k<6; ++k) {
                        for (int t=i+1; t<nOrigs; ++t) {
                            origs[k][t-1] = origs[k][t];
                        }
                        origs[k][nOrigs-1] = origs[k][j];
                    }

                    --nOrigs_new; // key origins -- 
                    --i;         // check the same i again
                }
            }

        }
    }

    cout << "len1 = "<< nOrigs <<", len2 = " << nOrigs_new << endl;

    for (int j=0; j<nOrigs_new; ++j){
        cout <<"-array_3- " << origs[0][j] <<", "<< origs[1][j]<<", "<< origs[2][j] <<", "
                            << origs[3][j] <<", "<< origs[4][j]<<", "<< origs[5][j] <<endl;
    }  




    /*
    Array<OneD, NekDouble> test_tmpGeom(8); // physical points
    Array<OneD, NekDouble> test_tmpGeom2(15);
    test_tmpGeom[0] = 0.365;
    test_tmpGeom[1] = 0.36532;
    test_tmpGeom[2] = 0.36602;
    test_tmpGeom[3] = 0.366976;
    test_tmpGeom[4] = 0.368023;
    test_tmpGeom[5] = 0.368979;
    test_tmpGeom[6] = 0.369679;
    test_tmpGeom[7] = 0.37;
    LibUtilities::Interp1D(from_key, &test_tmpGeom[0], to_key, &test_tmpGeom2[0]);
    for (int i=0; i<test_tmpGeom2.size(); ++i){
        cout <<"---x2--- = " << test_tmpGeom2[i] <<endl;
    }
    */
    //-------------------------------------------------------------------------
    //=========================================================================

    // Sampling setting
    const NekDouble distance_n = 0.005; // from wall to wall + H in normal direction
    const NekDouble x0 = 0.19;
    const NekDouble x1 = 0.4;
    const NekInt npts_n     = 21;    // npts in wall normal direction, use npts points for export
    const NekInt npts_x     = 3;     // npts in x direction
    cout << distance_n << ", " << npts_n<<", "<< x1-x0 << npts_x << endl;

    // Sampling points
    // h = 1- tanh((1-ksi)*atanh(sqrt(1-delta)))/sqrt(1-delta), ksi in [0,1]
    // from Agrawal's paper
    Array<OneD, NekDouble> h(npts_n);
    const NekDouble delta = 0.1;
    NekDouble tmp1;
    const NekDouble tmp2 = 1.0/(static_cast<NekDouble>(npts_n)-1.0); // 1/(npts-1)
    const NekDouble tmp3 = sqrt(1.0-delta);
    const NekDouble tmp4 = atanh(tmp3);
    const NekDouble tmp5 = 1.0/tmp3;
    for (int i=0;i<npts_n;++i){
        tmp1 = 1.0 - i * tmp2; // tmp1 = 1-ksi, ksi = i/(npts_n-1) belonging to [0,1]
        h[i] = 1 - tanh(tmp1*tmp4)*tmp5;
        //cout << i<<" - ksi = "<<1-tmp1<<", h = "<< h[i] <<endl;
    }

    // Get the sampled y/z on the wall according to given x
    // Call NekMesh or OpenCascade routine to get y/z based on x array
    // or directly interpolate using the available x/y/z/normals data on the bnd
    // temporarily use the following point, assume it's already on the wall
    NekDouble x_middle = (x1+x0)/2;
    NekDouble y_middle = 0.0;
    NekDouble z_middle = 0.0;
    NekDouble normals_x = 0.0;
    NekDouble normals_y = 1.0;
    NekDouble normals_z = 0.0;


    // Get target x/y/z array
    Array<OneD, NekDouble> x_target(npts_n, x_middle);
    Array<OneD, NekDouble> y_target(npts_n, y_middle);
    Array<OneD, NekDouble> z_target(npts_n, z_middle);

    for (int i=0;i<npts_n;++i){
        x_target[i] = x_target[i] + distance_n * h[i] * normals_x;
        y_target[i] = y_target[i] + distance_n * h[i] * normals_y;
        z_target[i] = z_target[i] + distance_n * h[i] * normals_z;
    }

    for(int i=0;i<npts_n;++i){
         cout << "#"<<i<<"-- "<<x_target[i]<<", "<<y_target[i]<<", "<<z_target[i]<<endl;
    }
   

    // Get variables' value
    // Mimic the interpolation procedure in 
    // /disk_two/Nek_Test/nektar++/library/FieldUtils/Interpolator.cpp
    int nCoordDim = m_f->m_exp[0]->GetCoordim(0); // tell the difference between nCoorDim and Dimension of the space
    cout << "Dimension = " << nCoordDim <<endl;

    Array<OneD, NekDouble> Lcoords(nCoordDim, 0.0); 
    Array<OneD, NekDouble> coords(3);
    coords[0] = x_target[20]; //0
    coords[1] = y_target[20]; //
    coords[2] = z_target[20]; //

    cout << "Point: "<<coords[0]<<", "<< coords[1]<<", "<<coords[2]<<endl;

    // Get donor element and local coordinates
    int elmtid = -1;
    elmtid = m_f->m_exp[0]->GetExpIndex(
            coords, Lcoords, NekConstants::kGeomFactorsTol, false, elmtid); // Get donor elmt 
    cout <<"elmtid = "<<elmtid<<endl;
    cout <<"Local coord: "<< Lcoords[0] <<", "<< Lcoords[1]<<endl;//", "<< Lcoords[2]<<endl;

    // Homogeneous case, need to find the right plane
    int targetPlane = -1;
    if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D) {

        cout << "e3DH1D" << endl;

        int nPlanes    = m_f->m_exp[0]->GetHomogeneousBasis()->GetZ().size();
        NekDouble lHom = m_f->m_exp[0]->GetHomoLen();
        targetPlane = std::round((coords[2]*nPlanes)/lHom);

        cout << "Plane = " << targetPlane << endl;

        // Reset from last plane to plane 0 (same physical result)
        if(targetPlane == nPlanes) {
            targetPlane = 0;
        }
    }

    // limit Lcoords to avoid warnings, ref Interpolator.cpp
    for (int j = 0; j < 2; ++j) {
        Lcoords[j] = std::max(Lcoords[j], -1.0);
        Lcoords[j] = std::min(Lcoords[j], 1.0);
    }
    
    // interpolate the value for each field
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
                // 1D/2D/3D and other cases [?]
                value = m_f->m_exp[f]->GetExp(elmtid)->StdPhysEvaluate(
                    Lcoords, m_f->m_exp[f]->GetPhys() + offset);
            }

            // Check and save
            if ((boost::math::isnan)(value)){
                ASSERTL0(false, "new value is not a number");
            }    
            else {
                // Save the value
                //m_ptsOutField->SetPointVal(m_ptsOutField->GetDim() + f, i, value);
                cout <<value <<endl;
            }
        }
    }
    else {
        cout << "Not coded yet."<<endl;
    }






    //----------
    /*
    if (m_spacedim == 2)
    {
        m_f->m_variables[0] = "Shear_x";
        m_f->m_variables[1] = "Shear_y";
        m_f->m_variables[2] = "Shear_mag";
    }
    else
    {
        m_f->m_variables[0] = "Shear_x";
        m_f->m_variables[1] = "Shear_y";
        m_f->m_variables[2] = "Shear_z";
        m_f->m_variables[3] = "Shear_mag";
    }
    */
}



}
}
