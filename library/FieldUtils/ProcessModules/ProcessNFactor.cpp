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
    // here is the place to register new options!!!
    // have fun
    // h,nh,ne,x1,x2,y1,y2,z1,z2
    //================================================
    m_config["newConfig"] = ConfigOption(false, "NotSet", "test the new option");

}

ProcessNFactor::~ProcessNFactor()
{
}


/* To heapify (max-haeading) a subtree rooted with node rootId
   A[0~5][0~currentLength-1] is the data to be heapified
   A.size()=6; A[0~5] for x/y/z/nx/ny/nz respectively
   The heap is adjusted regards to A[0] (x-value)
*/
void ProcessNFactor::Heapify_max(Array<OneD, Array<OneD, NekDouble> > A, 
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
    if (maxId != rootId) {
        for (int j=0; j<dataDim; ++j) { std::swap(A[j][rootId], A[j][maxId]); }
        Heapify_max(A, curLen, maxId);
    }

}

/* Sort the array using heap*/
void ProcessNFactor::HeapSort(Array<OneD, Array<OneD, NekDouble> > A)
{
    const int dataDim = A.size();
    const int totLen  = A[0].size();
    
    // Build max heap, starting from the last non-leaf node
    for (int i = floor(totLen / 2) - 1; i >= 0; --i) {
        Heapify_max(A, totLen, i);
    }
    // Move current root to end [0]->[currentlength-1]
    // and adjust the reduced heap, startig from root
    for (int curLen = totLen; curLen > 1; --curLen) {
        for (int j=0; j<dataDim; ++j) { std::swap(A[j][0], A[j][curLen-1]); }

        Heapify_max(A, curLen-1, 0);
    }

}

/* clean the repeated points in the array
   put the repeated points in the end 
   return the array length without repeated points
*/
int ProcessNFactor::CleanRepeatedPts(Array<OneD, Array<OneD, NekDouble> > A)
{
    const int dataDim = A.size();
    const int totLen = A[0].size();
    int       newLen = totLen; // initialize the new length as the total length
    for (int i=1; i<newLen; ++i) {
        // For each i, check indax smaller than i
        for (int j=i-1; j>=0; --j) {
            
            if ( abs(A[0][i]-A[0][j]) > NekConstants::kNekZeroTol ) {
                break; // The array has already sorted regards to x
            }
            else {
                if ( abs(A[1][i]-A[1][i-1]) < NekConstants::kNekZeroTol &&
                     abs(A[2][i]-A[2][i-1]) < NekConstants::kNekZeroTol ) {
                    
                    // repeated points found, [i]==[j]
                    // move [i+1] to [totLen-1] forword
                    for (int k=0; k<dataDim; ++k) {
                        for (int t=i+1; t<totLen; ++t) {
                            A[k][t-1] = A[k][t];
                        }
                        A[k][totLen-1] = A[k][j];
                    }

                    --newLen; // key origins -- 
                    --i;      // check the same i again
                }
            }

        }
    }

    return (newLen);
}


void ProcessNFactor::WriteDataInPts(const std::string &outFile, 
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
 
    if (isInc==true){
        if      (nfields==3) { is2D = true;  }
        else if (nfields==4) { is2D = false; }
        else  { cout << "Incorrect dimension." << endl; }
    }
    else {
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
    for (int i=0; i<nfields; ++i){
        ptsFile << m_f->m_variables[i];
        if (i < nfields-1) {
            ptsFile << ",";
        }
    }
    ptsFile << "\">" << endl;

    for (size_t i = 0; i<len; ++i){
        for (size_t j=0; j<data[0].size(); ++j) {
            ptsFile << "    ";
            ptsFile << data[i][j][0] << " " << data[i][j][1] << " ";
            if (is2D==false) { ptsFile << data[i][j][2] << " "; }
            for (int k=0; k<nfields; ++k){
                ptsFile << data[i][j][3+k];
                if (k < nfields-1) { ptsFile << " "; }
            }
            ptsFile << endl;
        }
    }

    ptsFile << "  </POINTS>" << endl;
    ptsFile << "</NEKTAR>" << endl;
    ptsFile.close();
}


void ProcessNFactor::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    //---------- Input paramaters (move to upstream routines later) -----------
    // sampling origins setting
    const int to_nPtsPerEdge = 6; // needed number of points per element
    Array<OneD, NekDouble> boundingBox(6); // use origins inside the box
    boundingBox[0] = 0.189;                               // x_lower
    boundingBox[1] = 0.401;                               // x_upper
    boundingBox[2] = -abs(NekConstants::kNekUnsetDouble); // y_lower
    boundingBox[3] =  abs(NekConstants::kNekUnsetDouble); // y_upper
    boundingBox[4] = -abs(NekConstants::kNekUnsetDouble); // z_lower
    boundingBox[5] =  abs(NekConstants::kNekUnsetDouble); // z_upper

    // Sampling setting
    const NekDouble distance_n = 0.005; // from wall to wall + H in normal direction
    const int       npts_n     = 21;    // npts in wall normal direction, use npts points for export
    const NekDouble delta      = 0.1;   // needs to be smaller than 1.
    //-------------------------------------------------------------------------

    // Step 1 - get data type
    int nfields   = m_f->m_variables.size();
    int nCoordDim = m_f->m_exp[0]->GetCoordim(0);
    m_spacedim    = nCoordDim + m_f->m_numHomogeneousDir;

    //string str_inc = "u";
    //string str_com = "rho";
    //if      (m_f->m_variables[0].compare(str_inc)==0){cout <<"incompressible"<<endl;}
    //else if (m_f->m_variables[0].compare(str_com)==0){cout <<"compressible"<<endl;}
    //else {ASSERTL0(false, "Other type of fields, not coded yet.");}
    //if(boost::iequals(m_f->m_variables[0], "u"))
    //if(boost::iequals(m_f->m_variables[0], "rho") && boost::iequals(m_f->m_variables[1], "rhou")


    // Step 2 - get boundary info
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

    //=========================================================================
    // Step 3 - set origins for sampling
    // only support 1D interpolation now
    // set dimensions bndElmtDim_para bndElmtDim_phys
    const int dim_para = BndExp[0]->GetExp(0)->GetNumBases(); // dimension for parametric coordinate system, eg. =1
    const int dim_phys = nCoordDim; // dimension for the physical space that the parametric coordinate system located on, eg. =2
    
    // set point key
    //const int to_nPtsPerElmt = 6; // needed number of points per element
    //LibUtilities::PointsKey from_key = BndExp[0]->GetExp(0)->GetBasis(0)->GetPointsKey();
    //const int from_nPtsPerElmt = from_key.GetNumPoints();
    LibUtilities::PointsType to_pointstype = LibUtilities::PointsType::ePolyEvenlySpaced;
    LibUtilities::PointsKey  to_key(to_nPtsPerEdge, to_pointstype);
    Array<OneD, LibUtilities::PointsKey> from_key(dim_para);
    int from_nPtsPerElmt = 1;
    int to_nPtsPerElmt   = 1; 
    for (int i=0; i<dim_para; ++i) {
        from_key[i] = BndExp[0]->GetExp(0)->GetBasis(i)->GetPointsKey();
        from_nPtsPerElmt *= from_key[i].GetNumPoints();
        to_nPtsPerElmt   *= to_nPtsPerEdge;
    }

    // declare arrays to save points
    Array<OneD, Array<OneD, NekDouble> > from_ptsInElmt(3);
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
    
    // loop the element on the bnd
    int ptr = 0; // pointed to head to copy the arrays
    for ( int i = 0; i < nElmts; ++i ) { //i < nElmts

        // obtain the points in the element
        BndExp[0]->GetExp(i)->GetCoords(from_ptsInElmt[0], from_ptsInElmt[1], from_ptsInElmt[2]); 

        // skip some elements, needs further improved
        if (from_ptsInElmt[0][0]<boundingBox[0] ||
            from_ptsInElmt[0][from_nPtsPerElmt-1]>boundingBox[1]) { continue; } 

        // interp x/y/z and nx/ny/nz
        // dim_phys determins the times (xy/xyz) to loop
        // dim_para determins the interpolation function (Interp1D/Interp2D) to ues
        // ref: Expansion::v_GetCoords in Expansion.cpp
        for (int j = 0; j < dim_phys; ++j ) {
            //LibUtilities::Interp1D(from_key, &xyz_bnd[j][i*from_nPtsPerElmt], to_key, &to_ptsInElmt[j][0]); //alternative code
            //LibUtilities::Interp1D(from_key, &from_ptsInElmt[j][0], to_key, &to_ptsInElmt[j][0]); //x/y/z
            //LibUtilities::Interp1D(from_key, &normals[j][i*from_nPtsPerElmt], to_key, &to_normalsInElmt[j][0]);
            switch (dim_para) {
                case 1: {
                    // ksi - line segement element on bnd
                    LibUtilities::Interp1D(
                        from_key[0], &from_ptsInElmt[j][0], 
                        to_key,      &to_ptsInElmt[j][0]); //x/y/z
                    LibUtilities::Interp1D(
                        from_key[0], &normals[j][i*from_nPtsPerElmt], 
                        to_key,      &to_normalsInElmt[j][0]); // nx/ny/nz
                    break;
                }
                case 2: {
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
                default: {
                    ASSERTL0(false, "Dimension error.");
                }
            }

            // save the interpolated results
            Vmath::Vcopy(to_nPtsPerElmt, &to_ptsInElmt[j][0],
                         1, &origs[j][ptr],   1); // copy coordinates
            Vmath::Vcopy(to_nPtsPerElmt, &to_normalsInElmt[j][0],
                         1, &origs[j+3][ptr], 1); // copy coordinates
        }
        ptr += to_nPtsPerElmt;

    }
 
    // sort array regards to x-value and remove repeated origin points
    HeapSort(origs);
    int nOrigs_new = CleanRepeatedPts(origs);
     
    //-------------------------------------------------------------------------
    //=========================================================================

    // Step 4 - set sampling points
    // Expression in Agrawal's paper:
    // h = 1- tanh((1-ksi)*atanh(sqrt(1-delta)))/sqrt(1-delta), ksi in [0,1]
    Array<OneD, NekDouble> h(npts_n);
    NekDouble tmp1;
    const NekDouble tmp2 = 1.0/(static_cast<NekDouble>(npts_n)-1.0); // 1/(npts-1)
    const NekDouble tmp3 = sqrt(1.0-delta);
    const NekDouble tmp4 = atanh(tmp3);
    const NekDouble tmp5 = 1.0/tmp3;
    for (int i=0; i<npts_n; ++i){
        tmp1 = 1.0 - i * tmp2; // tmp1 = 1-ksi, ksi = i/(npts_n-1) belonging to [0,1]
        h[i] = 1 - tanh(tmp1*tmp4)*tmp5;
        //cout << i <<" - ksi = "<<1-tmp1<<", h = "<< h[i] <<endl;
    }

    // declare the data array and fill in the coordinates
    // data[originId][normalId][variableId]
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > data(nOrigs_new);
    for (int i=0; i<nOrigs_new; ++i) {
        data[i] = Array<OneD, Array<OneD, NekDouble> >(npts_n);
        for (int j=0; j<npts_n; ++j) {
            data[i][j] = Array<OneD, NekDouble>(3+BndExp.size(), 0.0); //x/y/z+flow variables
            
            tmp1 = distance_n*h[j]; // physical distance in normal direction
            for (int k=0; k<3; ++k) {
                data[i][j][k] = origs[k][i] + tmp1*origs[k+3][i]; // x0+dist*nx
            }
        }
    }

    // Step 5 - interpolate the variables for each point
    Array<OneD, NekDouble> Lcoords(nCoordDim, 0.0); // 2.5D=2, 3D=3
    Array<OneD, NekDouble> coords(3);               // directly use 3
   
    for (int i=0; i<nOrigs_new; ++i) {
        for (int j=0; j<npts_n; ++j) {

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
    

    // Export data in pts format
    WriteDataInPts("111.pts", data, nOrigs_new);  //nOrigs_new

    
    //=========================================================================
    //vector<NekDouble> values;
    //bool test = ParseUtils::GenerateVector(m_config["plane"].as<string>(), values);

    /*
     * Process list of modules. Each element of the vector of module
     * strings can be in the following form:
     *
     * modname:arg1=a:arg2=b:arg3=c:arg4:arg5=asd
     * 
     * wallNormalDat:bnd=0:ne=11:h=0.001:nh=21 mesh.xml session.xml field.fld data.pts
     * 
     *
     * where the only required argument is 'modname', specifing the
     * name of the module to load.
     */
    
    vector<string> inout   = vm["input-file"].as<vector<string> >();
    vector<string> modcmds = vm["module"].as<vector<string> >();


    cout << "input.size = " << inout.size() << endl;
    for (int i=0; i<inout.size(); ++i){
        cout << "inout " << i << " = " << inout[i] << endl;
    }

    cout << "modcmds.size =" << modcmds.size() << endl;
    for (int i=0; i<modcmds.size(); ++i){
        cout << "modcmds " << i << " = " << modcmds[i] << endl;
    }


    cout << "---" << endl;
    vector<string> tmpStr1;
    boost::split(tmpStr1, modcmds[0], boost::is_any_of(":"));
    cout << "tmpStr1.size = " << tmpStr1.size() << endl;
    for (int i=0; i<tmpStr1.size(); ++i){
        cout << "tmpStr1 " << i << " = " << tmpStr1[i] << endl;
        vector<string> tmpStr2;
        boost::split(tmpStr2, tmpStr1[i], boost::is_any_of("="));
        for (int j=0; j<tmpStr2.size(); ++j){
            cout << "tmpStr2 " << j << " = " << tmpStr2[j] << endl;
        }
    }


    cout << "----" << endl;
    if ( m_config["bnd"].as<string>().compare("0") == 0 ) {
        cout << "0000" << endl;
    }
    else {
        cout << "1111" << endl;
    }

    if ( m_config["addnormals"].as<string>().compare("2") == 0 ) {
        cout << "0000" << endl;
    }
    else {
        cout << "1111" << endl;
    }
    if ( m_config["newConfig"].as<string>().compare("333") == 0 ) {
        cout << "0000" << endl;
    }
    else {
        cout << "1111" << endl;
    }









    //---------------------------------------------------------------
    // To do list (done)
    // 1. Will the following 3 expressions always the same?
    //    m_f->m_graph->GetSpaceDimension()
    //    m_f->m_exp[0]->GetCoordim(0)
    //    BndExp[0]->GetCoordim(0)
    // 2. PointTypes for quad bnd mesh and tri mesh on the boundary?
    //    PointsKey only accept 1 pointtype but there are 2 directions.
    // 3. Does this routine give the parametric dimension?
    //    BndExp[0]->GetExp(0)->GetNumBases()

    // ans:

    //---------------output some  middle results --------------------
    /*
    std::cout << "dim_para = " << dim_para <<", dim_coor = " << dim_phys <<std::endl;
    std::cout<< m_f->m_exp[0]->GetNumElmts() <<std::endl;
    std::cout<< m_f->m_variables[0]<<", "<<m_f->m_variables[1]<<", "
             << m_f->m_variables[2]<<", "<<m_f->m_variables[3]<<endl;
    cout << "m_f->m_exp.size() = " << m_f->m_exp.size() <<endl;

    cout << "bnd = " << bnd << endl;
    cout << "nqb = " << nqb << endl;
    cout << "normals_y = " << normals[1][nqb-2]<<" "<< normals[1][nqb-1] << endl; 
    cout << "Dimension = " << nCoordDim <<endl;
    */
    /*
    for (int i=0;i<nqb/4;++i){   // 0 ~ nqb/4, where 4 if for HomModesZ=4
        cout << i << " - " <<xyz_bnd[0][i] <<", "<<xyz_bnd[1][i]<<", "<<xyz_bnd[2][i]<<endl;
    }*/
    /*
    cout << "len1 = "<< nOrigs <<", len2 = " << nOrigs_new << endl;
    for (int j=0; j<nOrigs_new; ++j){
        cout <<"-array_3- " << origs[0][j] <<", "<< origs[1][j]<<", "<< origs[2][j] <<", "
                            << origs[3][j] <<", "<< origs[4][j]<<", "<< origs[5][j] <<endl;
    } 

    //---------------------------------------------------------------
    int i=0;
    cout << "======Result check======\nInput an origin index:"<<endl;
    cin >> i;
    while (i>=0 && i<nOrigs_new) {
        
        for (int j=0; j<npts_n; ++j) {
            cout << "#"<<j<<" - [" <<data[i][j][0] <<", "<<data[i][j][1]<<", "<< data[i][j][2]<<"]\n     "
                 <<data[i][j][3] <<", "<<data[i][j][4]<<", "<< data[i][j][5]<<", "<< data[i][j][6]<<endl;

        }
        system("pause");
        cout << "Dumped.\nInput a new origin index:"<<endl;
        cin >> i;
    }
    //---------------------------------------------------------------
    */

}



}
}
