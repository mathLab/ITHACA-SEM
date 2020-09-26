////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessGrad.cpp
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
//  Description: Computes gradient of fields.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessGrad.h"
#include "ProcessMapping.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessGrad::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "gradient"),
    ProcessGrad::create,
    "Computes gradient of fields.");

ProcessGrad::ProcessGrad(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessGrad::~ProcessGrad()
{
}

void ProcessGrad::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int i, j;
    int expdim    = m_f->m_graph->GetMeshDimension();
    int spacedim  = m_f->m_numHomogeneousDir + expdim;
    int nfields   = m_f->m_variables.size();
    int addfields = nfields * spacedim;

    for (i = 0; i < nfields; ++i)
    {
        if (spacedim == 1)
        {
            m_f->m_variables.push_back(m_f->m_variables[i] + "_x");
        }
        else if (spacedim == 2)
        {
            m_f->m_variables.push_back(m_f->m_variables[i] + "_x");
            m_f->m_variables.push_back(m_f->m_variables[i] + "_y");
        }
        else if (spacedim == 3)
        {
            m_f->m_variables.push_back(m_f->m_variables[i] + "_x");
            m_f->m_variables.push_back(m_f->m_variables[i] + "_y");
            m_f->m_variables.push_back(m_f->m_variables[i] + "_z");
        }
    }

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(addfields);
    m_f->m_exp.resize(nfields + addfields);

    for (i = 0; i < addfields; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    Array<OneD, Array<OneD, NekDouble> > tmp(spacedim);
    for (int i = 0; i < spacedim; i++)
    {
        tmp[i] = Array<OneD, NekDouble>(npoints);
    }

    // Get mapping
    GlobalMapping::MappingSharedPtr mapping = ProcessMapping::GetMapping(m_f);

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
                    m_f->m_exp[0]->HomogeneousBwdTrans(m_f->m_exp[i]->GetPhys(),
                                                       vel[i]);
                }
                else
                {
                    Vmath::Vcopy(npoints, m_f->m_exp[i]->GetPhys(), 1, vel[i],
                                 1);
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
                Vmath::Vcopy(npoints, m_f->m_exp[i]->GetPhys(), 1, vel[i], 1);
            }
        }
    }
    else
    {
        for (int i = 0; i < spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints, m_f->m_exp[i]->GetPhys(), 1, vel[i], 1);
        }
    }

    // Calculate Gradient
    for (i = 0; i < nfields; ++i)
    {
        for (j = 0; j < spacedim; ++j)
        {
            if (i < spacedim)
            {
                m_f->m_exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                         vel[i], tmp[j]);
            }
            else
            {
                m_f->m_exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                         m_f->m_exp[i]->GetPhys(), tmp[j]);
            }
        }
        mapping->CovarToCartesian(tmp, tmp);
        for (int j = 0; j < spacedim; j++)
        {
            Vmath::Vcopy(npoints, tmp[j], 1, grad[i * spacedim + j], 1);
        }
    }

    for (i = 0; i < addfields; ++i)
    {
        m_f->m_exp[nfields + i] =
            m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, grad[i], 1, m_f->m_exp[nfields + i]->UpdatePhys(),
                     1);
        m_f->m_exp[nfields + i]->FwdTrans_IterPerExp(
            grad[i], m_f->m_exp[nfields + i]->UpdateCoeffs());
    }
}
}
}
