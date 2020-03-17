////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInnerProduct.cpp
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
//  Description: Compute inner product between two fields.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessInnerProduct.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessInnerProduct::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "innerproduct"),
        ProcessInnerProduct::create,
        "take inner product between two fields and return value.");

ProcessInnerProduct::ProcessInnerProduct(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["fromfld"] = ConfigOption(
        false, "NotSet", "Fld file form which to interpolate field");
    m_config["fields"] =
        ConfigOption(false, "All", "field id's to be used in inner product");
    m_config["multifldids"] = ConfigOption(
        false, "NotSet", "Take inner product of multiple field fields with "
                         "ids given in string. i.e. file_0.chk file_1.chk ...");
    m_config["allfromflds"] =
        ConfigOption(true, "0", "Take inner product between all fromflds, "
                                     "requires multifldids to be set");
}

ProcessInnerProduct::~ProcessInnerProduct()
{
}

void ProcessInnerProduct::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    string fromfld           = m_config["fromfld"].as<string>();
    FieldSharedPtr fromField = std::shared_ptr<Field>(new Field());

    ASSERTL0(m_config["fromfld"].as<string>() != "NotSet",
             "The config parameter "
             "fromfld needs to be defined");

    // Set up ElementGIDs in case of parallel processing
    Array<OneD, int> ElementGIDs(m_f->m_exp[0]->GetExpSize());
    for (int i = 0; i < m_f->m_exp[0]->GetExpSize(); ++i)
    {
        ElementGIDs[i] = m_f->m_exp[0]->GetExp(i)->GetGeom()->GetGlobalID();
    }

    int nfields = m_f->m_variables.size();
    int nphys   = m_f->m_exp[0]->GetTotPoints();
    NekDouble totiprod;
    string fields = m_config["fields"].as<string>();
    vector<unsigned int> processFields;
    string multifldidsstr = m_config["multifldids"].as<string>();
    vector<unsigned int> multiFldIds;
    vector<string> fromfiles;
    bool allfromflds = m_config["allfromflds"].as<bool>();

    if (fields.compare("All") == 0)
    {
        for (int i = 0; i < nfields; ++i)
        {
            processFields.push_back(i);
        }
    }
    else
    {
        ASSERTL0(ParseUtils::GenerateVector(fields, processFields),
                 "Failed to interpret field string in module innerproduct");
    }

    if (multifldidsstr.compare("NotSet") == 0)
    {
        fromfiles.push_back(fromfld);
    }
    else
    {
        ASSERTL0(
            ParseUtils::GenerateSeqVector(multifldidsstr, multiFldIds),
            "Failed to interpret multifldids string in module innerproduct");
        int end        = fromfld.find_first_of('.', 0);
        string endstr  = fromfld.substr(end, fromfld.size());
        string bodystr = fromfld.substr(0, end);
        for (int i = 0; i < multiFldIds.size(); ++i)
        {
            string infile = bodystr + "_" +
                            boost::lexical_cast<string>(multiFldIds[i]) +
                            endstr;
            fromfiles.push_back(infile);
        }
    }

    Array<OneD, Array<OneD, NekDouble> > SaveFld(processFields.size());
    for (int j = 0; j < processFields.size(); ++j)
    {
        int fid    = processFields[j];
        SaveFld[j] = Array<OneD, NekDouble>(nphys);
        m_f->m_exp[fid]->BwdTrans(m_f->m_exp[fid]->GetCoeffs(), SaveFld[j]);
    }

    if (allfromflds == false)
    {

        for (int f = 0; f < fromfiles.size(); ++f)
        {
            m_f->FieldIOForFile(fromfiles[f])->Import(
                fromfiles[f], fromField->m_fielddef, fromField->m_data,
                LibUtilities::NullFieldMetaDataMap, ElementGIDs);

            totiprod = IProduct(processFields, fromField, SaveFld);

            if (m_f->m_comm->GetRank() == 0)
            {
                cout << "Inner Product WRT " << fromfiles[f] << " : "
                     << totiprod << endl;
            }
        }
    }
    else // evaluate all from fields, first by loading them all up and then
         // calling IProduct
    {

        // Load all from fields.
        Array<OneD, FieldSharedPtr> allFromField(fromfiles.size());
        for (int i = 0; i < fromfiles.size(); ++i)
        {
            allFromField[i] = std::shared_ptr<Field>(new Field());

            m_f->FieldIOForFile(fromfiles[i])->Import(
                fromfiles[i], allFromField[i]->m_fielddef,
                allFromField[i]->m_data, LibUtilities::NullFieldMetaDataMap,
                ElementGIDs);
        }

        for (int g = 0; g < fromfiles.size(); ++g)
        {
            for (int j = 0; j < processFields.size(); ++j)
            {
                int fid = processFields[j];

                // load new field
                for (int i = 0; i < allFromField[g]->m_data.size(); ++i)
                {
                    m_f->m_exp[fid]->ExtractDataToCoeffs(
                        allFromField[g]->m_fielddef[i],
                        allFromField[g]->m_data[i],
                        m_f->m_variables[fid],
                        m_f->m_exp[fid]->UpdateCoeffs());
                }

                m_f->m_exp[fid]->BwdTrans(m_f->m_exp[fid]->GetCoeffs(),
                                          SaveFld[j]);
            }

            // take inner product from this g field with all other above
            for (int f = g; f < fromfiles.size(); ++f)
            {
                totiprod = IProduct(processFields, allFromField[f], SaveFld);

                if (m_f->m_comm->GetRank() == 0)
                {
                    cout << "Inner Product of " << fromfiles[g] << " WRT "
                         << fromfiles[f] << " : " << totiprod << endl;
                }
            }
        }
    }
}

NekDouble ProcessInnerProduct::IProduct(
    vector<unsigned int> &processFields,
    FieldSharedPtr &fromField,
    Array<OneD, const Array<OneD, NekDouble> > &SaveFld)
{
    int nphys          = m_f->m_exp[0]->GetTotPoints();
    NekDouble totiprod = 0.0;

    for (int j = 0; j < processFields.size(); ++j)
    {
        int fid = processFields[j];

        // load new field
        for (int i = 0; i < fromField->m_data.size(); ++i)
        {
            m_f->m_exp[fid]->ExtractDataToCoeffs(
                fromField->m_fielddef[i], fromField->m_data[i],
                m_f->m_variables[fid],
                m_f->m_exp[fid]->UpdateCoeffs());
        }

        m_f->m_exp[fid]->BwdTrans(m_f->m_exp[fid]->GetCoeffs(),
                                  m_f->m_exp[fid]->UpdatePhys());

        Vmath::Vmul(nphys, SaveFld[j], 1, m_f->m_exp[fid]->GetPhys(), 1,
                    m_f->m_exp[fid]->UpdatePhys(), 1);

        NekDouble iprod =
            m_f->m_exp[fid]->PhysIntegral(m_f->m_exp[fid]->UpdatePhys());

        // put in parallel summation
        m_f->m_comm->AllReduce(iprod, Nektar::LibUtilities::ReduceSum);

        totiprod += iprod;
    }
    return totiprod;
}
}
}
