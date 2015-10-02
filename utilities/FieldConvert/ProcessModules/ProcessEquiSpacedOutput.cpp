////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessEquiSpacedOutput.cpp
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
//  Description: Set up fields as interpolation to equispaced output
//
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
using namespace std;

#include "ProcessEquiSpacedOutput.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessEquiSpacedOutput::className =
    GetModuleFactory().RegisterCreatorFunction(
           ModuleKey(eProcessModule, "equispacedoutput"),
           ProcessEquiSpacedOutput::create,
           "Write data as equi-spaced output using simplices to represent the data for connecting points");


ProcessEquiSpacedOutput::ProcessEquiSpacedOutput(FieldSharedPtr f)
    : ProcessModule(f)
{
    f->m_setUpEquiSpacedFields = true;

    m_config["tetonly"] = ConfigOption(true, "NotSet",
                                "Only process tetrahedral elements");

    m_config["modalenergy"] = ConfigOption(true,"NotSet","Write output as modal energy");

}

ProcessEquiSpacedOutput::~ProcessEquiSpacedOutput()
{
}

void ProcessEquiSpacedOutput::Process(po::variables_map &vm)
{
    SetupEquiSpacedField();
}


void ProcessEquiSpacedOutput::SetupEquiSpacedField(void)
{

    if(m_f->m_verbose)
    {
        cout << "Interpolating fields to equispaced" << endl;
    }

    int coordim  = m_f->m_exp[0]->GetCoordim(0);
    int shapedim = m_f->m_exp[0]->GetExp(0)->GetShapeDimension();
    int npts     = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > coords(3);

    int nel = m_f->m_exp[0]->GetExpSize();

    // set up the number of points in each element
    int newpoints = 0;
    int newtotpoints = 0;

    Array<OneD,int> conn;
    int prevNcoeffs = 0;
    int prevNpoints = 0;
    int cnt = 0;

    // identify face 1 connectivity for prisms
    map<int,StdRegions::Orientation > face0orient;
    set<int> prismorient;
    LocalRegions::ExpansionSharedPtr e;

    // prepare PtsField
    vector<std::string> fieldNames;
    vector<int> ppe;
    vector<Array<OneD, int> > ptsConn;
    int nfields;

    for(int i = 0; i < nel; ++i)
    {
        e = m_f->m_exp[0]->GetExp(i);
        if(e->DetShapeType() == LibUtilities::ePrism)
        {
            StdRegions::Orientation forient = e->GetForient(0);
            int fid = e->GetGeom()->GetFid(0);
            if(face0orient.count(fid))
            { // face 1 meeting face 1 so reverse this id
                prismorient.insert(i);
            }
            else
            {
                // just store if Dir 1 is fwd or bwd
                if((forient == StdRegions::eDir1BwdDir1_Dir2FwdDir2) ||
                   (forient == StdRegions::eDir1BwdDir1_Dir2BwdDir2) ||
                   (forient == StdRegions::eDir1BwdDir2_Dir2FwdDir1) ||
                   (forient == StdRegions::eDir1BwdDir2_Dir2BwdDir1))
                {
                    face0orient[fid] = StdRegions::eBwd;
                }
                else
                {
                    face0orient[fid] = StdRegions::eFwd;
                }
            }
        }
    }

    for(int i = 0; i < nel; ++i)
    {
        e = m_f->m_exp[0]->GetExp(i);
        if(e->DetShapeType() == LibUtilities::ePrism)
        {
            int fid = e->GetGeom()->GetFid(2);
            // check to see if face 2 meets face 1
            if(face0orient.count(fid))
            {
                // check to see how face 2 is orientated
                StdRegions::Orientation forient2 = e->GetForient(2);
                StdRegions::Orientation forient0 = face0orient[fid];

                // If dir 1 or forient2 is bwd then check agains
                // face 1 value
                if((forient2 == StdRegions::eDir1BwdDir1_Dir2FwdDir2) ||
                   (forient2 == StdRegions::eDir1BwdDir1_Dir2BwdDir2) ||
                   (forient2 == StdRegions::eDir1BwdDir2_Dir2FwdDir1) ||
                   (forient2 == StdRegions::eDir1BwdDir2_Dir2BwdDir1))
                {
                    if(forient0 == StdRegions::eFwd)
                    {
                        prismorient.insert(i);
                    }
                }
                else
                {
                    if(forient0 == StdRegions::eBwd)
                    {
                        prismorient.insert(i);
                    }
                }
            }
        }
    }

    for(int i = 0; i < nel; ++i)
    {
        e = m_f->m_exp[0]->GetExp(i);
        if(m_config["tetonly"].m_beenSet)
        {
            if(m_f->m_exp[0]->GetExp(i)->DetShapeType() !=
                    LibUtilities::eTetrahedron)
            {
                continue;
            }
        }

        switch(e->DetShapeType())
        {
        case LibUtilities::eSegment:
            {
                int npoints0 = e->GetBasis(0)->GetNumPoints();

                newpoints = LibUtilities::StdSegData::
                                    getNumberOfCoefficients(npoints0);
            }
            break;
        case LibUtilities::eTriangle:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np = max(np0,np1);
                newpoints = LibUtilities::StdTriData::
                                    getNumberOfCoefficients(np,np);
            }
            break;
        case LibUtilities::eQuadrilateral:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np = max(np0,np1);

                newpoints = LibUtilities::StdQuadData::
                                    getNumberOfCoefficients(np,np);
            }
            break;
        case LibUtilities::eTetrahedron:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np = max(np0,max(np1,np2));

                newpoints = LibUtilities::StdTetData::
                                    getNumberOfCoefficients(np,np,np);
            }
            break;
        case LibUtilities::ePrism:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np = max(np0,max(np1,np2));

                newpoints = LibUtilities::StdPrismData::
                                    getNumberOfCoefficients(np,np,np);
            }
            break;
        case LibUtilities::ePyramid:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np = max(np0,max(np1,np2));

                newpoints = LibUtilities::StdPyrData::
                                    getNumberOfCoefficients(np,np,np);
            }
            break;
        case LibUtilities::eHexahedron:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np = max(np0,max(np1,np2));

                newpoints = LibUtilities::StdPyrData::
                                    getNumberOfCoefficients(np,np,np);
            }
            break;
        default:
            {
                ASSERTL0(false,"Points not known");
            }
        }

        ppe.push_back(newpoints);
        newtotpoints += newpoints;


        if(e->DetShapeType() == LibUtilities::ePrism)
        {
            bool standard = true;

            if(prismorient.count(i))
            {
                standard = false; // reverse direction
            }

            e->GetSimplexEquiSpacedConnectivity(conn,standard);
        }
        else
        {

            if((prevNcoeffs != e->GetNcoeffs()) ||
               (prevNpoints != e->GetTotPoints()))
            {
                prevNcoeffs = e->GetNcoeffs();
                prevNpoints = e->GetTotPoints();

                e->GetSimplexEquiSpacedConnectivity(conn);
            }
        }
        Array<OneD, int> newconn(conn.num_elements());
        for(int j = 0; j < conn.num_elements(); ++j)
        {
            newconn[j] = conn[j] + cnt;
        }

        ptsConn.push_back(newconn);
        cnt += newpoints;
    }

    if(m_f->m_fielddef.size())
    {
        nfields = m_f->m_exp.size();
    }
    else // just the mesh points
    {
        nfields = 0;
    }

    Array<OneD, Array<OneD, NekDouble> > pts(nfields + coordim);

    for(int i = 0; i < nfields + coordim; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(newtotpoints);
    }

    // Interpolate coordinates
    for(int i = 0; i < coordim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(npts);
    }

    for(int i = coordim; i < 3; ++i)
    {
        coords[i] = NullNekDouble1DArray;
    }

    m_f->m_exp[0]->GetCoords(coords[0],coords[1],coords[2]);

    int nq1 = m_f->m_exp[0]->GetTotPoints();

    Array<OneD, NekDouble> x1(nq1);
    Array<OneD, NekDouble> y1(nq1);
    Array<OneD, NekDouble> z1(nq1);

    m_f->m_exp[0]->GetCoords(x1, y1, z1);


    Array<OneD, NekDouble> tmp;

    for(int n = 0; n < coordim; ++n)
    {
        cnt = 0;
        int cnt1 = 0;
        for(int i = 0; i < nel; ++i)
        {
            m_f->m_exp[0]->GetExp(i)->PhysInterpToSimplexEquiSpaced(
                                        coords[n] + cnt,
                                        tmp = pts[n] + cnt1);
            cnt1 += ppe[i];
            cnt  += m_f->m_exp[0]->GetExp(i)->GetTotPoints();
        }
    }

    if(m_f->m_fielddef.size())
    {
        ASSERTL0(m_f->m_fielddef[0]->m_fields.size() == m_f->m_exp.size(),
                 "More expansion defined than fields");

        for(int n = 0; n < m_f->m_exp.size(); ++n)
        {
            cnt = 0;
            int cnt1 = 0;

            if(m_config["modalenergy"].m_beenSet)
            {
                Array<OneD, const NekDouble> phys = m_f->m_exp[n]->GetPhys();
                for(int i = 0; i < nel; ++i)
                {
                    GenOrthoModes(i,phys+cnt,tmp = pts[coordim + n] + cnt1);
                    cnt1 += ppe[i];
                    cnt  += m_f->m_exp[0]->GetExp(i)->GetTotPoints();
                }
            }
            else
            {
                Array<OneD, const NekDouble> phys = m_f->m_exp[n]->GetPhys();
                for(int i = 0; i < nel; ++i)
                {
                    m_f->m_exp[0]->GetExp(i)->PhysInterpToSimplexEquiSpaced(
                            phys + cnt,
                            tmp = pts[coordim + n] + cnt1);
                    cnt1 += ppe[i];
                    cnt  += m_f->m_exp[0]->GetExp(i)->GetTotPoints();
                }
            }

            // Set up Variable string.
            fieldNames.push_back(m_f->m_fielddef[0]->m_fields[n]);
        }
    }

    m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(coordim, fieldNames, pts);
    if (shapedim == 2)
    {
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsTriBlock);
    }
    else if (shapedim == 3)
    {
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsTetBlock);
    }
    m_f->m_fieldPts->SetConnectivity(ptsConn);
}


    void ProcessEquiSpacedOutput::GenOrthoModes(
            int n,
            const Array<OneD,const NekDouble> &phys,
                  Array<OneD, NekDouble> &coeffs)
    {
        LocalRegions::ExpansionSharedPtr e;
        e = m_f->m_exp[0]->GetExp(n);

        switch(e->DetShapeType())
        {
            case LibUtilities::eTriangle:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np = max(np0,np1);

                // to ensure points are correctly projected to np need
                // to increase the order slightly of coordinates
                LibUtilities::PointsKey pa(np+1,e->GetPointsType(0));
                LibUtilities::PointsKey pb(np,e->GetPointsType(1));
                Array<OneD, NekDouble> tophys(np*(np+1));

                LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,np,pa);
                LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B,np,pb);
                StdRegions::StdTriExp OrthoExp(Ba,Bb);

                // interpolate points to new phys points!
                LibUtilities::Interp2D(e->GetBasis(0)->GetBasisKey(),
                                       e->GetBasis(1)->GetBasisKey(),
                                       phys,Ba,Bb,tophys);

                OrthoExp.FwdTrans(tophys,coeffs);
                break;
            }
            case LibUtilities::eQuadrilateral:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np = max(np0,np1);

                LibUtilities::PointsKey pa(np+1,e->GetPointsType(0));
                LibUtilities::PointsKey pb(np+1,e->GetPointsType(1));
                Array<OneD, NekDouble> tophys((np+1)*(np+1));

                LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,np,pa);
                LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A,np,pb);
                StdRegions::StdQuadExp OrthoExp(Ba,Bb);

                // interpolate points to new phys points!
                LibUtilities::Interp2D(e->GetBasis(0)->GetBasisKey(),
                                       e->GetBasis(1)->GetBasisKey(),
                                       phys,Ba,Bb,tophys);

                OrthoExp.FwdTrans(phys,coeffs);
                break;
            }
            default:
                ASSERTL0(false,"Shape needs setting up");
                break;
        }
    }

}
}


