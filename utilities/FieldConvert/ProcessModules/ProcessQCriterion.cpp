////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessQCriterion.cpp
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
//  Description: Computes Q Criterion field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessQCriterion.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessQCriterion::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "QCriterion"),
        ProcessQCriterion::create, "Computes Q-Criterion.");

ProcessQCriterion::ProcessQCriterion(FieldSharedPtr f)
    : ProcessModule(f)
{
}

ProcessQCriterion::~ProcessQCriterion()
{
}

void ProcessQCriterion::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessQCriterion: Calculating Q Criterion..." << endl;
    }

    int i, j, s;
    int expdim    = m_f->m_graph->GetMeshDimension();
    int spacedim  = expdim;
    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
        (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
    {
        spacedim = 3;
    }
    int nfields = m_f->m_fielddef[0]->m_fields.size();
    if (spacedim == 1 || spacedim == 2)
    {
        cerr << "\n Error: ProcessQCriterion must be computed for a 3D"
        " (or quasi-3D) case. \n" << endl;
    }

    //For calculating Q-Criterion only 1 field must be added
    int addfields = 1;

    int npoints = m_f->m_exp[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble> > grad(nfields * nfields);

    Array<OneD, Array<OneD, NekDouble> > omega(nfields * nfields);
    Array<OneD, Array<OneD, NekDouble> > S    (nfields * nfields);

    Array<OneD, Array<OneD, NekDouble> > outfield (addfields);
    Array<OneD, Array<OneD, NekDouble> > outfield1(addfields);
    Array<OneD, Array<OneD, NekDouble> > outfield2(addfields);
    Array<OneD, Array<OneD, NekDouble> > outfield3(addfields);

    int nstrips;

    m_f->m_session->LoadParameter("Strip_Z",nstrips,1);

    m_f->m_exp.resize(nfields*nstrips);

    for (i = 0; i < nfields*nfields; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    for (i = 0; i < addfields; ++i)
    {
        //Will store the Q-Criterion
        outfield[i] = Array<OneD, NekDouble>(npoints);
        outfield1[i] = Array<OneD, NekDouble>(npoints);
        outfield2[i] = Array<OneD, NekDouble>(npoints);
        outfield3[i] = Array<OneD, NekDouble>(npoints);

        omega[i] = Array<OneD, NekDouble>(npoints);
        S[i] = Array<OneD, NekDouble>(npoints);
    }

    vector<MultiRegions::ExpListSharedPtr> Exp(nstrips*addfields);

    for(s = 0; s < nstrips; ++s) //homogeneous strip varient
    {
        for (i = 0; i < nfields; ++i)
        {
            m_f->m_exp[s*nfields+i]->PhysDeriv(m_f->m_exp[s*nfields+i]->GetPhys(),
                                               grad[i*nfields],
                                               grad[i*nfields+1],
                                               grad[i*nfields+2]);
        }

        // W_x = Wy - Vz
        Vmath::Vsub(npoints, grad[2 * nfields + 1], 1,
                    grad[1 * nfields + 2], 1,
                    outfield1[0],          1);
        // W_x^2
        Vmath::Vmul(npoints, outfield1[0],          1,
                    outfield1[0],          1,
                    outfield1[0],          1);

        // W_y = Uz - Wx
        Vmath::Vsub(npoints, grad[0 * nfields + 2], 1,
                    grad[2 * nfields + 0], 1,
                    outfield2[0],          1);
        // W_y^2
        Vmath::Vmul(npoints, outfield2[0],          1,
                    outfield2[0],          1,
                    outfield2[0],          1);

        // W_z = Vx - Uy
        Vmath::Vsub(npoints, grad[1 * nfields + 0], 1,
                    grad[0 * nfields + 1], 1,
                    outfield3[0],          1);
        // W_z^2
        Vmath::Vmul(npoints, outfield3[0],          1,
                    outfield3[0],          1,
                    outfield3[0],          1);

        // add fields omega = 0.5*(W_x^2 + W_y^2 + W_z^2)

        NekDouble fac = 0.5;
        Vmath::Vadd(npoints, &outfield1[0][0],      1,
                    &outfield2[0][0],      1,
                    &omega[0][0],          1);
        Vmath::Vadd(npoints, &omega[0][0],          1,
                    &outfield3[0][0],      1,
                    &omega[0][0],          1);

        for (int k = 0; k < addfields; ++k)
        {
            Vmath::Smul(npoints, fac, &omega[k][0], 1, &omega[k][0], 1);
        }

        Vmath::Zero(npoints, &outfield1[0][0], 1);
        Vmath::Zero(npoints, &outfield2[0][0], 1);
        Vmath::Zero(npoints, &outfield3[0][0], 1);

        Vmath::Vmul(npoints, grad[0 * nfields + 0], 1,
                    grad[0 * nfields + 0], 1,
                    outfield1[0],          1);
        Vmath::Vmul(npoints, grad[1 * nfields + 1], 1,
                    grad[1 * nfields + 1], 1,
                    outfield2[0],          1);
        Vmath::Vmul(npoints, grad[2 * nfields + 2], 1,
                    grad[2 * nfields + 2], 1,
                    outfield3[0],          1);

        Vmath::Vadd(npoints, &outfield1[0][0], 1,
                    &outfield2[0][0], 1,
                    &S[0][0], 1);
        Vmath::Vadd(npoints, &S[0][0], 1,
                    &outfield3[0][0], 1,
                    &S[0][0], 1);

        // W_y + V_z
        Vmath::Vadd(npoints, grad[2 * nfields + 1], 1,
                    grad[1 * nfields + 2], 1,
                    outfield1[0], 1);
        Vmath::Vmul(npoints, &outfield1[0][0], 1,
                    &outfield1[0][0], 1,
                    &outfield1[0][0], 1);

        // U_z + W_x
        Vmath::Vadd(npoints, grad[0 * nfields + 2], 1,
                    grad[2 * nfields + 0], 1,
                    outfield2[0], 1);
        Vmath::Vmul(npoints, &outfield2[0][0], 1,
                    &outfield2[0][0], 1,
                    &outfield2[0][0], 1);

        // V_x + U_y
        Vmath::Vadd(npoints, grad[1 * nfields + 0], 1,
                    grad[0 * nfields + 1], 1,
                    outfield3[0], 1);
        Vmath::Vmul(npoints, &outfield3[0][0], 1,
                    &outfield3[0][0], 1,
                    &outfield3[0][0], 1);

        Vmath::Vadd(npoints, &outfield1[0][0], 1,
                    &outfield2[0][0], 1,
                    &outfield2[0][0], 1);
        Vmath::Vadd(npoints, &outfield2[0][0], 1,
                    &outfield3[0][0], 1,
                    &outfield3[0][0], 1);

        for (int k = 0; k < addfields; ++k)
        {
            Vmath::Smul(npoints, fac, &outfield3[k][0], 1,
                        &outfield3[k][0], 1);
        }

        Vmath::Vadd(npoints, &outfield3[0][0], 1, &S[0][0], 1, &S[0][0], 1);
        Vmath::Vsub(npoints, omega[0], 1, S[0], 1, outfield[0], 1);

        for (int k = 0; k < addfields; ++k)
        {
            Vmath::Smul(npoints, fac, &outfield[k][0], 1,
                        &outfield[k][0], 1);
        }


        for (i = 0; i < addfields; ++i)
        {
            int n = s*addfields + i;
            Exp[n] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
            Exp[n]->UpdatePhys() = outfield[i];
            Exp[n]->FwdTrans(outfield[i],
                             Exp[n]->UpdateCoeffs());
        }
    }

    vector<MultiRegions::ExpListSharedPtr>::iterator it;

    for(s = 0; s < nstrips; ++s)
    {
        for(i = 0; i < addfields; ++i)
        {
            it = m_f->m_exp.begin()+s*(nfields+addfields)+nfields+i;
            m_f->m_exp.insert(it, Exp[s*addfields+i]);
        }
    }

    vector<string> outname;
    outname.push_back("Q");

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(s = 0; s < nstrips; ++s) //homogeneous strip varient
    {
        for (j = 0; j < nfields + addfields; ++j)
        {
            for (i = 0; i < FieldDef.size()/nstrips; ++i)
            {
                int n = s * FieldDef.size()/nstrips + i;

                if (j >= nfields)
                {
                    FieldDef[n]->m_fields.push_back(outname[j-nfields]);
                }
                else
                {
                    FieldDef[n]->m_fields.push_back(
                        m_f->m_fielddef[0]->m_fields[j]);
                }
                m_f->m_exp[s*(nfields + addfields)+j]->AppendFieldData(FieldDef[n], FieldData[n]);
            }
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

}
}
