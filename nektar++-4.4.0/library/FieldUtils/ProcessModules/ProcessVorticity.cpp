////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVorticity.cpp
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
//  Description: Computes vorticity field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessMapping.h"
#include "ProcessVorticity.h"
#include <GlobalMapping/Mapping.h>

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessVorticity::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "vorticity"),
        ProcessVorticity::create,
        "Computes vorticity field.");

ProcessVorticity::ProcessVorticity(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessVorticity::~ProcessVorticity()
{
}

void ProcessVorticity::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "ProcessVorticity: Calculating vorticity..." << endl;
        }
    }

    int i, j, s;
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = expdim;
    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
        (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
    {
        spacedim = 3;
    }
    int nfields = m_f->m_fielddef[0]->m_fields.size();
    if (spacedim == 1)
    {
        ASSERTL0(false, "Error: Vorticity for a 1D problem cannot "
                        "be computed")
    }
    int addfields = (spacedim == 2) ? 1 : 3;

    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(spacedim * spacedim);
    Array<OneD, Array<OneD, NekDouble> > outfield(addfields);

    int nstrips;

    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    m_f->m_exp.resize(nfields * nstrips);

    for (i = 0; i < spacedim * spacedim; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    for (i = 0; i < addfields; ++i)
    {
        outfield[i] = Array<OneD, NekDouble>(npoints);
    }

    Array<OneD, Array<OneD, NekDouble> > tmp(spacedim);
    for (int i = 0; i < spacedim; i++)
    {
        tmp[i] = Array<OneD, NekDouble>(npoints);
    }

    vector<MultiRegions::ExpListSharedPtr> Exp(nstrips * addfields);

    // Get mapping
    GlobalMapping::MappingSharedPtr mapping = ProcessMapping::GetMapping(m_f);

    for (s = 0; s < nstrips; ++s) // homogeneous strip varient
    {
        // Get velocity and convert to Cartesian system,
        //      if it is still in transformed system
        Array<OneD, Array<OneD, NekDouble> > vel(spacedim);
        if (m_f->m_fieldMetaDataMap.count("MappingCartesianVel"))
        {
            if (m_f->m_fieldMetaDataMap["MappingCartesianVel"] == "False")
            {
                // Initialize arrays and copy velocity
                for (int i = 0; i < spacedim; ++i)
                {
                    vel[i] = Array<OneD, NekDouble>(npoints);
                    if (m_f->m_exp[0]->GetWaveSpace())
                    {
                        m_f->m_exp[0]->HomogeneousBwdTrans(
                            m_f->m_exp[s * nfields + i]->GetPhys(), vel[i]);
                    }
                    else
                    {
                        Vmath::Vcopy(npoints,
                                     m_f->m_exp[s * nfields + i]->GetPhys(), 1,
                                     vel[i], 1);
                    }
                }
                // Convert velocity to cartesian system
                mapping->ContravarToCartesian(vel, vel);
                // Convert back to wavespace if necessary
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    for (int i = 0; i < spacedim; ++i)
                    {
                        m_f->m_exp[0]->HomogeneousFwdTrans(vel[i], vel[i]);
                    }
                }
            }
            else
            {
                for (int i = 0; i < spacedim; ++i)
                {
                    vel[i] = Array<OneD, NekDouble>(npoints);
                    Vmath::Vcopy(npoints,
                                 m_f->m_exp[s * nfields + i]->GetPhys(), 1,
                                 vel[i], 1);
                }
            }
        }
        else
        {
            for (int i = 0; i < spacedim; ++i)
            {
                vel[i] = Array<OneD, NekDouble>(npoints);
                Vmath::Vcopy(npoints, m_f->m_exp[s * nfields + i]->GetPhys(), 1,
                             vel[i], 1);
            }
        }

        // Calculate Gradient & Vorticity
        if (spacedim == 2)
        {
            for (i = 0; i < spacedim; ++i)
            {
                m_f->m_exp[s * nfields + i]->PhysDeriv(vel[i], tmp[0], tmp[1]);
                mapping->CovarToCartesian(tmp, tmp);
                for (int j = 0; j < spacedim; j++)
                {
                    Vmath::Vcopy(npoints, tmp[j], 1, grad[i * spacedim + j], 1);
                }
            }
            // W_z = Vx - Uy
            Vmath::Vsub(npoints, grad[1 * spacedim + 0], 1,
                        grad[0 * spacedim + 1], 1, outfield[0], 1);
        }
        else
        {
            for (i = 0; i < spacedim; ++i)
            {
                m_f->m_exp[s * nfields + i]->PhysDeriv(vel[i], tmp[0], tmp[1],
                                                       tmp[2]);
                mapping->CovarToCartesian(tmp, tmp);
                for (int j = 0; j < spacedim; j++)
                {
                    Vmath::Vcopy(npoints, tmp[j], 1, grad[i * spacedim + j], 1);
                }
            }

            // W_x = Wy - Vz
            Vmath::Vsub(npoints, grad[2 * spacedim + 1], 1,
                        grad[1 * spacedim + 2], 1, outfield[0], 1);
            // W_y = Uz - Wx
            Vmath::Vsub(npoints, grad[0 * spacedim + 2], 1,
                        grad[2 * spacedim + 0], 1, outfield[1], 1);
            // W_z = Vx - Uy
            Vmath::Vsub(npoints, grad[1 * spacedim + 0], 1,
                        grad[0 * spacedim + 1], 1, outfield[2], 1);
        }

        for (i = 0; i < addfields; ++i)
        {
            int n = s * addfields + i;
            Exp[n] =
                m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
            Vmath::Vcopy(npoints, outfield[i], 1, Exp[n]->UpdatePhys(), 1);
            Exp[n]->FwdTrans_IterPerExp(outfield[i], Exp[n]->UpdateCoeffs());
        }
    }

    vector<MultiRegions::ExpListSharedPtr>::iterator it;
    for (s = 0; s < nstrips; ++s)
    {
        for (i = 0; i < addfields; ++i)
        {
            it = m_f->m_exp.begin() + s * (nfields + addfields) + nfields + i;
            m_f->m_exp.insert(it, Exp[s * addfields + i]);
        }
    }

    vector<string> outname;
    if (addfields == 1)
    {
        outname.push_back("W_z");
    }
    else
    {
        outname.push_back("W_x");
        outname.push_back("W_y");
        outname.push_back("W_z");
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (s = 0; s < nstrips; ++s) // homogeneous strip varient
    {
        for (j = 0; j < nfields + addfields; ++j)
        {
            for (i = 0; i < FieldDef.size() / nstrips; ++i)
            {
                int n = s * FieldDef.size() / nstrips + i;

                if (j >= nfields)
                {
                    FieldDef[n]->m_fields.push_back(outname[j - nfields]);
                }
                else
                {
                    FieldDef[n]->m_fields.push_back(
                        m_f->m_fielddef[0]->m_fields[j]);
                }
                m_f->m_exp[s * (nfields + addfields) + j]->AppendFieldData(
                    FieldDef[n], FieldData[n]);
            }
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}
}
}
