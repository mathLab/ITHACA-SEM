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
//  Description: Add mapping coordinates to field
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessMapping.h"
#include <GlobalMapping/Mapping.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessMapping::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "mapping"),
        ProcessMapping::create, "Add mapping coordinates to output file.");

ProcessMapping::ProcessMapping(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessMapping::~ProcessMapping()
{
}

void ProcessMapping::Process(po::variables_map &vm)
{
    ASSERTL0(m_f->m_fieldMetaDataMap.count("MappingType"),
            "Failed to get mapping information from fld file.");
    
    // Get time from metadata
    string s_time = m_f->m_fieldMetaDataMap["Time"];
    NekDouble time = atof(s_time.c_str());
    
    // Determine dimensions of mesh, solution, etc...
    int npoints = m_f->m_exp[0]->GetNpoints();
    int expdim    = m_f->m_graph->GetMeshDimension();
    int spacedim  = expdim;
    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
        (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
    {
        spacedim = 3;
    }
    int nfields = m_f->m_fielddef[0]->m_fields.size();
    int addfields = spacedim;
    m_f->m_exp.resize(nfields+addfields);
    
    // Declare coordinates storage
    Array<OneD, Array<OneD, NekDouble> > coords_new(3);
    for (int i = 0; i < 3; i++)
    {
        coords_new[i]  = Array<OneD, NekDouble> (npoints);
    }    
    string fieldNames[3] = {"x", "y", "z"};
    
    // Evaluate coordinates
    if (m_f->m_fieldMetaDataMap["MappingType"] == "Expression")
    {            
        // Get name of the function
        string funcName = m_f->m_fieldMetaDataMap["MappingExpression"];
        ASSERTL0(m_f->m_session->DefinesFunction(funcName),
                        "Function '" + funcName + "' not defined.");
        
        // Get original coordinates (in case some of them are not changed)
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        for (int i = 0; i < 3; i++)
        {
            coords[i]      = Array<OneD, NekDouble> (npoints);
        }       
        m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);
        
        // Load coordinates
        std::string s_FieldStr; 
        for(int i = 0; i < 3; i++)
        {
            s_FieldStr = fieldNames[i];
            if ( m_f->m_session->DefinesFunction(funcName, s_FieldStr))
            {
                LibUtilities::EquationSharedPtr ffunc =
                        m_f->m_session->GetFunction(funcName, s_FieldStr);
                ffunc->Evaluate(coords[0], coords[1], coords[2], 
                                            time, coords_new[i]);
            }
            else
            {
                // This coordinate is not defined, so use (x^i)' = x^i
                Vmath::Vcopy(npoints, coords[i], 1, coords_new[i], 1);
            }
        }                
    }
    else
    {
        string fileName = m_f->m_fieldMetaDataMap["FileName"];
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;    
        std::vector<std::vector<NekDouble> > FieldData;
        
        m_f->m_fld->Import(fileName,
                            FieldDef,
                            FieldData);

        for (int j = 0; j < spacedim; ++j)
        {
            int ncoeffs = m_f->m_exp[0]->GetNcoeffs();
            Array<OneD, NekDouble> fieldcoeffs(ncoeffs,0.0);
            for (int i = 0; i < FieldData.size(); ++i)
            {
                m_f->m_exp[j]->ExtractDataToCoeffs(FieldDef[i],
                                                    FieldData[i],
                                                    fieldNames[j],
                                                           fieldcoeffs);
            }
            bool wavespace = m_f->m_exp[0]->GetWaveSpace();
            m_f->m_exp[0]->SetWaveSpace(false);
            
            m_f->m_exp[0]->BwdTrans(fieldcoeffs,
                                    coords_new[j]);
            
            m_f->m_exp[0]->SetWaveSpace(wavespace);
        }             
    }
    
    // Convert velocity to Cartesian system
    if (m_f->m_fieldMetaDataMap["MappingCorrectVel"] != "True")
    {
        m_f->m_fieldMetaDataMap["MappingCorrectVel"] = "True";
        
        Array<OneD, MultiRegions::ExpListSharedPtr>  field(1);
        field[0] = m_f->m_exp[0];
        GlobalMapping::MappingSharedPtr mapping = 
                                GlobalMapping::Mapping::Load(m_f->m_session,
                                field);
        // Update mapping with coordinates
        //     the coordinates velocity don't affect the transformation,
        //      so they can be set to zero
        Array<OneD, Array<OneD, NekDouble> > coords_vel(3);
        for (int i = 0; i < 3; i++)
        {
            coords_vel[i]  = Array<OneD, NekDouble> (npoints,0.0);
        }        
        mapping->UpdateMapping(time, false, coords_new,coords_vel);
        
        Array<OneD, Array<OneD, NekDouble> > vel (spacedim);
        Array<OneD, Array<OneD, NekDouble> > velCart (spacedim);
        // Initialize arrays and copy velocity
        for ( int i =0; i<spacedim; ++i )
        {
            vel[i] = Array<OneD, NekDouble> (npoints);
            velCart[i] = Array<OneD, NekDouble> (npoints);                    
            Vmath::Vcopy(npoints, m_f->m_exp[i]->GetPhys(), 1, vel[i], 1);
        }
        // Convert velocity to cartesian system
        mapping->ContravarToCartesian(vel, velCart);            
        // Copy result back
        for ( int i =0; i<spacedim; ++i )
        {
            Vmath::Vcopy(npoints, velCart[i], 1, m_f->m_exp[i]->UpdatePhys(), 1);
            m_f->m_exp[i]->FwdTrans_IterPerExp(m_f->m_exp[i]->GetPhys(),
                                                m_f->m_exp[i]->UpdateCoeffs());
        }        
    }
    
    // Add new information to m_f
    vector<string > outname;
    for (int i = 0; i < addfields; ++i)
    {
        m_f->m_exp[nfields + i] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
        m_f->m_exp[nfields + i]->UpdatePhys() = coords_new[i];
        m_f->m_exp[nfields + i]->FwdTrans_IterPerExp(coords_new[i],
                            m_f->m_exp[nfields + i]->UpdateCoeffs());
        outname.push_back(fieldNames[i]);
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (int j = 0; j < nfields + addfields; ++j)
    {
        for (int i = 0; i < FieldDef.size(); ++i)
        {   
            if (j >= nfields)
            {
                FieldDef[i]->m_fields.push_back(outname[j-nfields]);
            }
            else
            {
                FieldDef[i]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[j]);
            }
            m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;

}

}
}
