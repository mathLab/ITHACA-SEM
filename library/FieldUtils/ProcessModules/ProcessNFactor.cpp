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
    
    
    Array<OneD, NekDouble> x(nqb);
    Array<OneD, NekDouble> y(nqb);
    Array<OneD, NekDouble> z(nqb);
    cout << x.size() <<", "<< y.size() <<endl;
    //Array<OneD, NekDouble> y(nqb);
    //cout << m_f->m_exp.size() << endl;
    //cout << m_f->m_exp[0]->GetCoordim(0) <<endl;
    m_f->m_exp[0]->GetCoords(x,y,z);
    

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
