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

#include <boost/core/ignore_unused.hpp>

#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessMapping.h"
#include "ProcessVorticity.h"

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
    boost::ignore_unused(vm);

    int i, s;
    int expdim   = m_f->m_graph->GetMeshDimension();
    m_spacedim = expdim;
    if ((m_f->m_numHomogeneousDir) == 1 || (m_f->m_numHomogeneousDir) == 2)
    {
        m_spacedim = 3;
    }
    int nfields = m_f->m_variables.size();
    ASSERTL0(m_spacedim != 1,
             "Error: Vorticity for a 1D problem cannot be computed");
    int addfields = (m_spacedim == 2) ? 1 : 3;

    // Append field names
    if (addfields == 1)
    {
        m_f->m_variables.push_back("W_z");
    }
    else
    {
        m_f->m_variables.push_back("W_x");
        m_f->m_variables.push_back("W_y");
        m_f->m_variables.push_back("W_z");
    }

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }
    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(m_spacedim * m_spacedim);
    Array<OneD, Array<OneD, NekDouble> > outfield(addfields);

    int nstrips;

    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    for (i = 0; i < m_spacedim * m_spacedim; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    for (i = 0; i < addfields; ++i)
    {
        outfield[i] = Array<OneD, NekDouble>(npoints);
    }

    Array<OneD, Array<OneD, NekDouble> > tmp(m_spacedim);
    for (int i = 0; i < m_spacedim; i++)
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
        Array<OneD, Array<OneD, NekDouble> > vel(m_spacedim);
        GetVelocity(vel, s);
        if (m_f->m_fieldMetaDataMap.count("MappingCartesianVel"))
        {
            if (m_f->m_fieldMetaDataMap["MappingCartesianVel"] == "False")
            {
                // Initialize arrays and copy velocity
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        m_f->m_exp[0]->HomogeneousBwdTrans(vel[i], vel[i]);
                    }
                }
                // Convert velocity to cartesian system
                mapping->ContravarToCartesian(vel, vel);
                // Convert back to wavespace if necessary
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    for (int i = 0; i < m_spacedim; ++i)
                    {
                        m_f->m_exp[0]->HomogeneousFwdTrans(vel[i], vel[i]);
                    }
                }
            }
        }

        // Calculate Gradient & Vorticity
        if (m_spacedim == 2)
        {
            for (i = 0; i < m_spacedim; ++i)
            {
                m_f->m_exp[s * nfields + i]->PhysDeriv(vel[i], tmp[0], tmp[1]);
                mapping->CovarToCartesian(tmp, tmp);
                for (int j = 0; j < m_spacedim; j++)
                {
                    Vmath::Vcopy(npoints, tmp[j], 1, grad[i * m_spacedim + j], 1);
                }
            }
            // W_z = Vx - Uy
            Vmath::Vsub(npoints, grad[1 * m_spacedim + 0], 1,
                        grad[0 * m_spacedim + 1], 1, outfield[0], 1);
        }
        else
        {
            for (i = 0; i < m_spacedim; ++i)
            {
                m_f->m_exp[s * nfields + i]->PhysDeriv(vel[i], tmp[0], tmp[1],
                                                       tmp[2]);
                mapping->CovarToCartesian(tmp, tmp);
                for (int j = 0; j < m_spacedim; j++)
                {
                    Vmath::Vcopy(npoints, tmp[j], 1, grad[i * m_spacedim + j], 1);
                }
            }

            // W_x = Wy - Vz
            Vmath::Vsub(npoints, grad[2 * m_spacedim + 1], 1,
                        grad[1 * m_spacedim + 2], 1, outfield[0], 1);
            // W_y = Uz - Wx
            Vmath::Vsub(npoints, grad[0 * m_spacedim + 2], 1,
                        grad[2 * m_spacedim + 0], 1, outfield[1], 1);
            // W_z = Vx - Uy
            Vmath::Vsub(npoints, grad[1 * m_spacedim + 0], 1,
                        grad[0 * m_spacedim + 1], 1, outfield[2], 1);
        }

        for (i = 0; i < addfields; ++i)
        {
            int n = s * addfields + i;
            Exp[n] =
                m_f->AppendExpList(m_f->m_numHomogeneousDir);
            Vmath::Vcopy(npoints, outfield[i], 1, Exp[n]->UpdatePhys(), 1);
            Exp[n]->FwdTrans_IterPerExp(outfield[i], Exp[n]->UpdateCoeffs());
        }
    }

    for (s = 0; s < nstrips; ++s)
    {
        for (i = 0; i < addfields; ++i)
        {
            m_f->m_exp.insert(
                m_f->m_exp.begin() + s * (nfields + addfields) + nfields + i,
                Exp[s * addfields + i]);
        }
    }
}

void ProcessVorticity::GetVelocity( Array<OneD, Array<OneD, NekDouble> > &vel,
                                    int strip)
{
    int nfields = m_f->m_variables.size();
    int npoints = m_f->m_exp[0]->GetNpoints();
    if(boost::iequals(m_f->m_variables[0], "u"))
    {
        // IncNavierStokesSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints,
                         m_f->m_exp[strip * nfields + i]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
    else if(boost::iequals(m_f->m_variables[0], "rho") &&
            boost::iequals(m_f->m_variables[1], "rhou"))
    {
        // CompressibleFlowSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vdiv(npoints,
                         m_f->m_exp[strip * nfields + i + 1]->GetPhys(), 1,
                         m_f->m_exp[strip * nfields + 0    ]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
    else
    {
        // Unknown
        ASSERTL0(false, "Could not identify velocity for ProcessVorticity");
    }
}

}
}
