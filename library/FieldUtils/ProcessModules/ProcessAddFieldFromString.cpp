///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessAddFieldFromString.cpp
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
//  Description: Add a new field from a string based on existing variable
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>
using namespace std;

#include "ProcessAddFieldFromString.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessAddFieldFromString::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "addfieldfromstring"),
        ProcessAddFieldFromString::create,
        "Create a new field from the existing field as specified by a string"
        " using a required qrgument of the form fieldstr=\"x + y + u\" ");

ProcessAddFieldFromString::ProcessAddFieldFromString(FieldSharedPtr f)
    : ProcessModule(f)
{

    m_config["fieldstr"] =
        ConfigOption(false, "NotSet", "string of new field to be added (required)");
    m_config["fieldname"] =
        ConfigOption(false, "newfield", "name for isocontour if fieldstr "
                                      "specified, default is newfield (optional)");
}

ProcessAddFieldFromString::~ProcessAddFieldFromString(void)
{
}

void ProcessAddFieldFromString::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "ProcessAddFieldFromString: Processing new field..." << endl;
        }
    }
    
    if (m_config["fieldstr"].m_beenSet) // generate field of interest
    {
        ASSERTL0(m_f->m_fielddef.size(),"Expect a field file to be provided with this method");
        
        int nstrips;
        m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);
        ASSERTL0(nstrips == 1,"Routine is currently only setup for non-strip files");

        int npoints = m_f->m_exp[0]->GetTotPoints();
        int nfields = m_f->m_exp.size();
        
        // set up field name if provided otherwise called "newfield" from
        // default
        string fieldName = m_config["fieldname"].as<string>();
        
        m_f->m_exp.resize(nfields+1);
        m_f->m_exp[nfields] = m_f->AppendExpList(m_f->m_fielddef[0]->
                                                 m_numHomogeneousDir,fieldName);
        
        // get hold of an array of coordinates and the physical fields values
        string varstr; 
        vector<Array<OneD, const NekDouble> > interpfields;
        
        int Dim = m_f->m_graph->GetMeshDimension();
        if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
            (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
        {
            Dim = 3;
        }
        
        switch(Dim)
        {
        case 1:
        {
            varstr += "x ";
            Array<OneD, NekDouble> x(npoints);
            m_f->m_exp[0]->GetCoords(x);
            interpfields.push_back(x);
        }
        break;
        case 2:
        {
            varstr += "x y";
            
            Array<OneD, NekDouble> x(npoints);
            Array<OneD, NekDouble> y(npoints);
            m_f->m_exp[0]->GetCoords(x,y);
            interpfields.push_back(x);
                interpfields.push_back(y);
        }
        break;
        case 3:
        {
            varstr += "x y z";
            
            Array<OneD, NekDouble> x(npoints);
            Array<OneD, NekDouble> y(npoints);
            Array<OneD, NekDouble> z(npoints);
            m_f->m_exp[0]->GetCoords(x,y,z);
            interpfields.push_back(x);
            interpfields.push_back(y);
            interpfields.push_back(z);
        }
        }

        // add on the field values
        for (int i = 0; i < nfields; ++i)
        {
            varstr += " " + m_f->m_fielddef[0]->m_fields[i];
            m_f->m_exp[i]->BwdTrans(m_f->m_exp[i]->GetCoeffs(),
                                    m_f->m_exp[i]->UpdatePhys());
            interpfields.push_back(m_f->m_exp[i]->GetPhys());
        }

        // evaluate new function
        LibUtilities::AnalyticExpressionEvaluator strEval;
        int ExprId      = -1;
        std::string fieldstr = m_config["fieldstr"].as<string>();
        ExprId          = strEval.DefineFunction(varstr.c_str(), fieldstr);
        
        strEval.Evaluate(ExprId, interpfields, m_f->m_exp[nfields]->UpdatePhys());
        
        m_f->m_exp[nfields]->FwdTrans_IterPerExp(m_f->m_exp[nfields]->GetPhys(),m_f->m_exp[nfields]->UpdateCoeffs());


        // add in new field to fielddef
        for (int i = 0; i < m_f->m_fielddef.size(); ++i)
        {
            m_f->m_fielddef[i]->m_fields.push_back(fieldName);

            m_f->m_exp[nfields]->AppendFieldData(m_f->m_fielddef[i],
                                                 m_f->m_data[i]);
        }
    }
    else
    {
        ASSERTL0(false, "fieldstrmust be specified");
    }
}
}
}
