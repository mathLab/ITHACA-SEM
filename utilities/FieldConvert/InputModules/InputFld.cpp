////////////////////////////////////////////////////////////////////////////////
//
//  File: InputFld.cpp
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
//  Description: Reads a Nektar++ FLD file.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "InputFld.h"


static std::string npts = LibUtilities::SessionReader::RegisterCmdLineArgument(
                "NumberOfPoints","n","Define number of points to dump output");

namespace Nektar
{
namespace Utilities
{

ModuleKey InputFld::m_className[4] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "fld"), InputFld::create,
        "Reads Fld file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "chk"), InputFld::create,
        "Reads checkpoint file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "rst"), InputFld::create,
        "Reads restart file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "bse"), InputFld::create,
        "Reads stability base-flow file.")
};


/**
 * @brief Set up InputFld object.
 *
 */
InputFld::InputFld(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("fld");
    m_allowedFiles.insert("chk");
    m_allowedFiles.insert("rst");
    m_allowedFiles.insert("bse");
}


/**
 *
 */
InputFld::~InputFld()
{
}


/**
 *
 */
void InputFld::Process(po::variables_map &vm)
{
    if(m_f->m_verbose)
    {
        cout << "Processing input fld file" << endl;
    }

    int i,j;
    string fldending;
    //Determine appropriate field input
    if(m_f->m_inputfiles.count("fld") != 0)
    {
        fldending = "fld";
    }
    else if(m_f->m_inputfiles.count("chk") != 0)
    {
        fldending = "chk";
    }
    else if (m_f->m_inputfiles.count("rst") != 0)
    {
        fldending = "rst";
    }
    else if (m_f->m_inputfiles.count("bse") != 0)
    {
        fldending = "bse";
    }
    else
    {
        ASSERTL0(false,"no input file found");
    }

    if(!m_f->m_fld)
    {
        if(m_f->m_session)
        {
            m_f->m_fld = MemoryManager<LibUtilities::FieldIO>
                ::AllocateSharedPtr(m_f->m_session->GetComm());
        }
        else // serial communicator
        {
            LibUtilities::CommSharedPtr c = LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
            m_f->m_fld = MemoryManager<LibUtilities::FieldIO>
                ::AllocateSharedPtr(c);
        }
    }


    if(m_f->m_graph)  // all for restricted expansion defintion when loading field
    {
        // currently load all field (possibly could read data from expansion list
        // but it is re-arranged in expansion)

        const SpatialDomains::ExpansionMap &expansions = m_f->m_graph->GetExpansions();

        // if Range has been speficied it is possible to have a
        // partition which is empty so ccheck this and return if
        // no elements present.
        if(!expansions.size())
        {
            return;
        }

        m_f->m_exp.resize(1);

        Array<OneD,int> ElementGIDs(expansions.size());
        SpatialDomains::ExpansionMap::const_iterator expIt;

        i = 0;
        for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
        {
            ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
        }

        m_f->m_fielddef.clear();
        m_f->m_data.clear();

        m_f->m_fld->Import(m_f->m_inputfiles[fldending][0],
                           m_f->m_fielddef,
                           m_f->m_data,
                           m_f->m_fieldMetaDataMap,
                           ElementGIDs);
    }
    else // load all data.
    {
        m_f->m_fld->Import(m_f->m_inputfiles[fldending][0],
                           m_f->m_fielddef,
                           m_f->m_data,
                           m_f->m_fieldMetaDataMap);
    }


    // if m_exp defined presume we want to load all field  into expansions
    if(m_f->m_exp.size())
    {
        int nfields,nstrips;

        m_f->m_session->LoadParameter("Strip_Z",nstrips,1);

        if(vm.count("useSessionVariables"))
        {
            nfields = m_f->m_session->GetVariables().size();
        }
        else
        {
            nfields = m_f->m_fielddef[0]->m_fields.size();
        }


        m_f->m_exp.resize(nfields*nstrips);

        vector<string> vars = m_f->m_session->GetVariables();

        // declare other fields;
        for (int s = 0; s < nstrips; ++s) //homogeneous strip varient
        {
            for (i = 0; i < nfields; ++i)
            {
                if(i < vars.size())
                {
                     m_f->m_exp[s*nfields+i] = m_f->AppendExpList(
                             m_f->m_fielddef[0]->m_numHomogeneousDir,
                             vars[i]);
                }
                else
                {
                    if(vars.size())
                    {
                        m_f->m_exp[s*nfields+i] =
                            m_f->AppendExpList(
                                m_f->m_fielddef[0]->m_numHomogeneousDir,
                                vars[0]);
                    }
                    else
                    {
                        m_f->m_exp[s*nfields+i] =
                            m_f->AppendExpList(
                                m_f->m_fielddef[0]->m_numHomogeneousDir);
                    }
                }
            }
        }

        for(int s = 0; s < nstrips; ++s) //homogeneous strip varient
        {
            for (j = 0; j < nfields; ++j)
            {
                for (i = 0; i < m_f->m_data.size()/nstrips; ++i)
                {
                    m_f->m_exp[s*nfields+j]->
                              ExtractDataToCoeffs(m_f->m_fielddef[i*nstrips+s],
                    m_f->m_data[i*nstrips+s],
                    m_f->m_fielddef[i*nstrips+s]
                            ->m_fields[j],
                    m_f->m_exp[s*nfields+j]->UpdateCoeffs());
                }
                m_f->m_exp[s*nfields+j]->BwdTrans(m_f->m_exp[s*nfields+j]->GetCoeffs(),
                            m_f->m_exp[s*nfields+j]->UpdatePhys());
            }
        }

        // if range is defined reset up output field in case or
        // reducing fld definition
        if(vm.count("range"))
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = m_f->m_exp[0]->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

            for (j = 0; j < nfields; ++j)
            {
                for (i = 0; i < FieldDef.size(); ++i)
                {
                    FieldDef[i]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[j]);
                    m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
                }
            }
            m_f->m_fielddef = FieldDef;
            m_f->m_data     = FieldData;
        }
    }
}

}
}
