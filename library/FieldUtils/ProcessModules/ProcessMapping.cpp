////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessMapping.cpp
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
//  Description: Add mapping coordinates to field
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessMapping.h"

namespace Nektar
{
namespace FieldUtils
{
ModuleKey ProcessMapping::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "mapping"),
        ProcessMapping::create,
        "Add mapping coordinates to output file.");

ProcessMapping::ProcessMapping(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessMapping::~ProcessMapping()
{
}

void ProcessMapping::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Determine dimensions of mesh, solution, etc...
    int npoints  = m_f->m_exp[0]->GetNpoints();
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = expdim;
    if ((m_f->m_numHomogeneousDir) == 1 || (m_f->m_numHomogeneousDir) == 2)
    {
        spacedim = 3;
    }
    int nfields   = m_f->m_variables.size();
    int addfields = expdim;

    string fieldNames[3] = {"xCoord", "yCoord", "zCoord"};
    for (int i = 0; i < addfields; ++i)
    {
        m_f->m_variables.push_back(fieldNames[i]);
    }

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    m_f->m_exp.resize(nfields + addfields);

    // Load mapping
    GlobalMapping::MappingSharedPtr mapping = GetMapping(m_f);

    // Convert velocity to Cartesian system
    if (m_f->m_fieldMetaDataMap.count("MappingCartesianVel"))
    {
        if (m_f->m_fieldMetaDataMap["MappingCartesianVel"] == "False")
        {
            m_f->m_fieldMetaDataMap["MappingCartesianVel"] = "True";

            Array<OneD, Array<OneD, NekDouble> > vel(spacedim);
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
            // Copy result back
            for (int i = 0; i < spacedim; ++i)
            {
                if (m_f->m_exp[0]->GetWaveSpace())
                {
                    m_f->m_exp[0]->HomogeneousFwdTrans(
                        vel[i], m_f->m_exp[i]->UpdatePhys());
                }
                else
                {
                    Vmath::Vcopy(npoints, vel[i], 1,
                                 m_f->m_exp[i]->UpdatePhys(), 1);
                }
                m_f->m_exp[i]->FwdTrans_IterPerExp(
                    m_f->m_exp[i]->GetPhys(), m_f->m_exp[i]->UpdateCoeffs());
            }
        }
    }

    // Get coordinates from mapping
    Array<OneD, Array<OneD, NekDouble> > coords(3);
    mapping->GetCartesianCoordinates(coords[0], coords[1], coords[2]);

    // Add new information to m_f
    for (int i = 0; i < addfields; ++i)
    {
        m_f->m_exp[nfields + i] =
            m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, coords[i], 1,
                     m_f->m_exp[nfields + i]->UpdatePhys(), 1);
        m_f->m_exp[nfields + i]->FwdTrans_IterPerExp(
            coords[i], m_f->m_exp[nfields + i]->UpdateCoeffs());
    }
}

GlobalMapping::MappingSharedPtr ProcessMapping::GetMapping(FieldSharedPtr f)
{
    // Create mapping object
    Array<OneD, MultiRegions::ExpListSharedPtr> field(1);
    field[0] = f->m_exp[0];
    GlobalMapping::MappingSharedPtr mapping =
        GlobalMapping::Mapping::Load(f->m_session, field);

    // Get time from metadata
    NekDouble time;
    if (f->m_fieldMetaDataMap.count("Time"))
    {
        string s_time = f->m_fieldMetaDataMap["Time"];
        time          = atof(s_time.c_str());
    }
    else
    {
        time = 0.0;
    }

    // Get field information
    int npoints  = f->m_exp[0]->GetNpoints();
    int expdim   = f->m_graph->GetMeshDimension();
    int spacedim = expdim + f->m_numHomogeneousDir;

    // Declare coordinates storage
    Array<OneD, Array<OneD, NekDouble> > coords_new(3);
    Array<OneD, Array<OneD, NekDouble> > coords_vel(3);
    for (int i = 0; i < 3; i++)
    {
        coords_new[i] = Array<OneD, NekDouble>(npoints);
        coords_vel[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    string fieldNames[3]    = {"x", "y", "z"};
    string velFieldNames[3] = {"vx", "vy", "vz"};

    // Evaluate coordinates and coordinates velocity
    if (f->m_fieldMetaDataMap.count("MappingType"))
    {
        if (f->m_fieldMetaDataMap["MappingType"] == "Expression")
        {
            // Get name of the functions
            string funcName;
            string velFuncName;
            if (f->m_fieldMetaDataMap.count("MappingExpression"))
            {
                funcName = f->m_fieldMetaDataMap["MappingExpression"];
            }
            else
            {
                funcName = "";
            }
            if (f->m_fieldMetaDataMap.count("MappingVelExpression"))
            {
                velFuncName = f->m_fieldMetaDataMap["MappingVelExpression"];
            }
            else
            {
                velFuncName = "";
            }

            // Get original coordinates (in case some of them are not changed)
            Array<OneD, Array<OneD, NekDouble> > coords(3);
            for (int i = 0; i < 3; i++)
            {
                coords[i] = Array<OneD, NekDouble>(npoints);
            }
            f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);

            // Load coordinates
            std::string s_FieldStr;
            for (int i = 0; i < 3; i++)
            {
                s_FieldStr = fieldNames[i];
                if (f->m_session->DefinesFunction(funcName, s_FieldStr))
                {
                    LibUtilities::EquationSharedPtr ffunc =
                        f->m_session->GetFunction(funcName, s_FieldStr);
                    ffunc->Evaluate(coords[0], coords[1], coords[2], time,
                                    coords_new[i]);
                }
                else
                {
                    // This coordinate is not defined, so use (x^i)' = x^i
                    Vmath::Vcopy(npoints, coords[i], 1, coords_new[i], 1);
                }
            }
            // Load velocities
            if (f->m_session->DefinesFunction(velFuncName))
            {
                for (int i = 0; i < 3; i++)
                {
                    s_FieldStr = velFieldNames[i];
                    if (f->m_session->DefinesFunction(velFuncName, s_FieldStr))
                    {
                        LibUtilities::EquationSharedPtr ffunc =
                            f->m_session->GetFunction(velFuncName, s_FieldStr);
                        ffunc->Evaluate(coords[0], coords[1], coords[2], time,
                                        coords_vel[i]);
                    }
                }
            }

            // Update mapping with coordinates
            mapping->SetFromFunction(false);
            mapping->UpdateMapping(time, coords_new, coords_vel);
        }
        else if (f->m_fieldMetaDataMap["MappingType"] == "File")
        {
            ASSERTL0(f->m_fieldMetaDataMap.count("FileName"),
                     "FileName parameter for Mapping missing in field file.");
            string fileName = f->m_fieldMetaDataMap["FileName"];
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;

            f->FieldIOForFile(fileName)->Import(fileName, FieldDef, FieldData);

            for (int j = 0; j < spacedim; ++j)
            {
                int ncoeffs = f->m_exp[0]->GetNcoeffs();
                Array<OneD, NekDouble> fieldcoeffs(ncoeffs, 0.0);
                for (int i = 0; i < FieldData.size(); ++i)
                {
                    f->m_exp[j]->ExtractDataToCoeffs(
                        FieldDef[i], FieldData[i], fieldNames[j], fieldcoeffs);
                }
                bool wavespace = f->m_exp[0]->GetWaveSpace();
                f->m_exp[0]->SetWaveSpace(false);

                f->m_exp[0]->BwdTrans(fieldcoeffs, coords_new[j]);

                // Load coordinate velocity
                if (std::find(FieldDef[0]->m_fields.begin(),
                              FieldDef[0]->m_fields.end(),
                              velFieldNames[j]) != FieldDef[0]->m_fields.end())
                {
                    for (int i = 0; i < FieldData.size(); ++i)
                    {
                        f->m_exp[j]->ExtractDataToCoeffs(
                            FieldDef[i], FieldData[i], velFieldNames[j],
                            fieldcoeffs);
                    }
                    f->m_exp[0]->BwdTrans(fieldcoeffs, coords_vel[j]);
                }
                f->m_exp[0]->SetWaveSpace(wavespace);
            }
            // Update mapping with coordinates
            mapping->SetFromFunction(false);
            mapping->UpdateMapping(time, coords_new, coords_vel);
        }
    }
    else
    {
        // Use trivial mapping
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        Array<OneD, Array<OneD, NekDouble> > coords_vel(3);
        for (int i = 0; i < 3; i++)
        {
            coords[i]     = Array<OneD, NekDouble>(npoints);
            coords_vel[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }
        f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);
        mapping->SetFromFunction(false);
        mapping->UpdateMapping(time, coords, coords_vel);
    }

    return mapping;
}
}
}
