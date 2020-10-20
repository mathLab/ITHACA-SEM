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
//  Description: Computes the N-factor along the surface.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include "ProcessNFactor.h"
#include <NekMesh/CADSystem/CADCurve.h>


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
    "Computes the N-factor along the surface.");

ProcessNFactor::ProcessNFactor(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
}

ProcessNFactor::~ProcessNFactor()
{
}

void ProcessNFactor::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    int i;
    int nfields = m_f->m_variables.size();
    int expdim  = m_f->m_graph->GetSpaceDimension();
    m_spacedim  = expdim + m_f->m_numHomogeneousDir;

    std::cout<<"Inside the N-factor module!" << std::endl;
    std::cout<< nfields<< ", " << expdim << ", " << m_spacedim <<std::endl;
    std::cout<< m_f->m_exp[0]->GetNumElmts() <<std::endl;

    // Declare arrays
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(m_spacedim + 1); //[?]
    Array<OneD, MultiRegions::ExpListSharedPtr> BndElmtExp(nfields);


    // Create map of boundary ids for partitioned domains [?]
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
    cout << bnd << endl;

    // Get expansion list for boundary and for elements containing this
    // bnd
    // But why there is m_exp[m_spacedim] [?]
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
    int nqe = BndElmtExp[0]->GetTotPoints();

    // Get inward-pointing wall-normal vectors 
    Array<OneD, Array<OneD, NekDouble> > normals; 
    m_f->m_exp[0]->GetBoundaryNormals(bnd, normals);
    // Reverse normals, to get correct orientation for the body
    // normals[i][j], where i is direction varying from 0 to 2;
    // j is the point varying from 0 to (nqb-1)
    for (i = 0; i < m_spacedim; ++i)
    {
        Vmath::Neg(nqb, normals[i], 1);
    }
    
    cout <<"normals1 "<< normals[1][nqb-2]<<" "<< normals[1][nqb-1] << endl; 
    cout << "nqb = " << nqb << ", nqe = " << nqe <<endl;
    
    
    Array<OneD, NekDouble> x_bnd(nqb);
    Array<OneD, NekDouble> y_bnd(nqb);
    Array<OneD, NekDouble> z_bnd(nqb);
    cout << x_bnd.size() <<", "<< y_bnd.size() << ", " << z_bnd.size() << endl;
    // m_f->m_exp[0]->GetCoords(x,y,z); // is the coordinates of the whole field 
    //cout << m_f->m_exp.size() << endl;
    //cout << m_f->m_exp[0]->GetCoordim(0) <<endl; // output  = 2
    BndExp[0]->GetCoords(x_bnd,y_bnd,z_bnd);
    // cout << BndExp.size() << endl; // 4 for u,v,w,p
    // cout << BndExp[0]->GetCoordim(0) << endl; // 2 for expansion in the x-y plane
    for (int i=8;i<16;++i){   // 0 ~ nqb/4, where 4 if for HomModesZ=4
        cout << i << " - " <<x_bnd[i] <<", "<<y_bnd[i]<<", "<<z_bnd[i]<<endl;
    }

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
        cout << i<<" - ksi = "<<1-tmp1<<", h = "<< h[i] <<endl;
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
    Array<OneD, NekDouble> Lcoords(2, 0.0); 
    Array<OneD, NekDouble> coords(2);
    coords[0] = x_target[0];
    coords[1] = y_target[0];
    //coords[2] = z_target[0];

    cout << "Point: "<<coords[0]<<", "<< coords[1]<<endl;//", "<<coords[2]<<endl;

    // Get donor element
    int elmtid = -1;
    elmtid = m_f->m_exp[0]->GetExpIndex(
            coords, Lcoords, NekConstants::kGeomFactorsTol, false, elmtid); // Get donor elmt 
    cout <<"elmtid = "<<elmtid<<endl;
    
    // limit Lcoords to avoid warnings, ref Interpolator.cpp
    for (int j = 0; j < 2; ++j) {
        Lcoords[j] = std::max(Lcoords[j], -1.0);
        Lcoords[j] = std::min(Lcoords[j], 1.0);
    }

    if (elmtid >= 0) {


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
