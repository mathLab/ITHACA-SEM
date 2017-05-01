////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFld.cpp
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
//  Description: FLD file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputFld.h"
#include <LibUtilities/BasicUtils/FileSystem.h>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputFld::m_className[2] = {
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "fld"),
                                               OutputFld::create,
                                               "Writes a Fld file."),
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "chk"),
                                               OutputFld::create,
                                               "Writes a Fld file."),
};

OutputFld::OutputFld(FieldSharedPtr f) : OutputFileBase(f)
{
    m_config["format"] = ConfigOption(
        false, "Xml", "Output format of field file");
}

OutputFld::~OutputFld()
{
}

void OutputFld::OutputFromPts(po::variables_map &vm)
{
    //ASSERTL0(false, "OutputFld can't write using Pts information.");
    OutputFromExp(vm);
}

void OutputFld::OutputFromExp(po::variables_map &vm)
{
    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    // Set up communicator and FieldIO object.
    LibUtilities::CommSharedPtr c = m_f->m_session ? m_f->m_session->GetComm() :
        LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
    LibUtilities::FieldIOSharedPtr fld =
        LibUtilities::GetFieldIOFactory().CreateInstance(
            m_config["format"].as<string>(), c, true);

    if (m_f->m_writeBndFld)
    {
        if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
        {
            cout << "\t OutputFld: Writing boundary file(s): ";
            for (int i = 0; i < m_f->m_bndRegionsToWrite.size(); ++i)
            {
                cout << m_f->m_bndRegionsToWrite[i];
                if (i < m_f->m_bndRegionsToWrite.size() - 1)
                {
                    cout << ", ";
                }
            }
            cout << endl;
        }

        // Extract data to boundaryconditions
        if (m_f->m_fldToBnd)
        {
            for (int i = 0; i < m_f->m_exp.size(); ++i)
            {
                m_f->m_exp[i]->FillBndCondFromField();
            }
        }

        int nfields = m_f->m_exp.size();
        Array<OneD, Array<OneD, const MultiRegions::ExpListSharedPtr> > BndExp(
            nfields);
        for (int i = 0; i < nfields; ++i)
        {
            BndExp[i] = m_f->m_exp[i]->GetBndCondExpansions();
        }

        // get hold of partition boundary regions so we can match it to desired
        // region extraction
        SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                               m_f->m_exp[0]->GetGraph());
        const SpatialDomains::BoundaryRegionCollection bregions =
            bcs.GetBoundaryRegions();
        SpatialDomains::BoundaryRegionCollection::const_iterator breg_it;
        map<int, int> BndRegionMap;
        int cnt = 0;
        for (breg_it = bregions.begin(); breg_it != bregions.end();
             ++breg_it, ++cnt)
        {
            BndRegionMap[breg_it->first] = cnt;
        }

        // find ending of output file and insert _b1, _b2
        int dot     = filename.find_last_of('.') + 1;
        string ext  = filename.substr(dot, filename.length() - dot);
        string name = filename.substr(0, dot - 1);

        for (int i = 0; i < m_f->m_bndRegionsToWrite.size(); ++i)
        {
            string outname =
                name + "_b" +
                boost::lexical_cast<string>(m_f->m_bndRegionsToWrite[i]) + "." +
                ext;

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;

            if (BndRegionMap.count(m_f->m_bndRegionsToWrite[i]) == 1)
            {
                int Border = BndRegionMap[m_f->m_bndRegionsToWrite[i]];

                FieldDef = BndExp[0][Border]->GetFieldDefinitions();
                FieldData.resize(FieldDef.size());

                for (int j = 0; j < nfields; ++j)
                {
                    for (int k = 0; k < FieldDef.size(); ++k)
                    {
                        BndExp[j][Border]->AppendFieldData(FieldDef[k],
                                                           FieldData[k]);
                        FieldDef[k]->m_fields.push_back(m_f->m_variables[j]);
                    }
                }

                if (m_f->m_addNormals)
                {
                    int normdim       = m_f->m_graph->GetMeshDimension();
                    string normstr[3] = {"Norm_x", "Norm_y", "Norm_z"};

                    // Get normals
                    Array<OneD, Array<OneD, NekDouble> > NormPhys;
                    m_f->m_exp[0]->GetBoundaryNormals(Border, NormPhys);

                    // add normal coefficients to list to be dumped
                    for (int j = 0; j < normdim; ++j)
                    {
                        BndExp[0][Border]->FwdTrans(
                            NormPhys[j], BndExp[0][Border]->UpdateCoeffs());

                        for (int k = 0; k < FieldDef.size(); ++k)
                        {
                            BndExp[0][Border]->AppendFieldData(FieldDef[k],
                                                               FieldData[k]);
                            FieldDef[k]->m_fields.push_back(normstr[j]);
                        }
                    }
                }

                // output error for regression checking.
                if (vm.count("error"))
                {
                    int rank = m_f->m_session->GetComm()->GetRank();

                    for (int j = 0; j < nfields; ++j)
                    {
                        BndExp[j][Border]->BwdTrans(
                            BndExp[j][Border]->GetCoeffs(),
                            BndExp[j][Border]->UpdatePhys());

                        // Note currently these calls will
                        // hange since not all partitions will
                        // call error.
                        NekDouble l2err =
                            BndExp[j][Border]->L2(BndExp[j][Border]->GetPhys());

                        NekDouble linferr = BndExp[j][Border]->Linf(
                            BndExp[j][Border]->GetPhys());

                        if (rank == 0)
                        {
                            cout << "L 2 error (variable "
                                 << FieldDef[0]->m_fields[j] << ") : " << l2err
                                 << endl;

                            cout << "L inf error (variable "
                                 << FieldDef[0]->m_fields[j]
                                 << ") : " << linferr << endl;
                        }
                    }
                }
            }

            fld->Write(outname, FieldDef, FieldData, m_f->m_fieldMetaDataMap);
        }
    }
    else
    {
        int i, j, s;
        int nfields = m_f->m_variables.size();
        int nstrips;
        m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
            m_f->m_exp[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
        for (s = 0; s < nstrips; ++s)
        {
            for (j = 0; j < nfields; ++j)
            {
                for (i = 0; i < FieldDef.size() / nstrips; ++i)
                {
                    int n = s * FieldDef.size() / nstrips + i;

                    FieldDef[n]->m_fields.push_back(m_f->m_variables[j]);
                    m_f->m_exp[s * nfields + j]->AppendFieldData(
                        FieldDef[n], FieldData[n]);
                }
            }
        }
        fld->Write(filename, FieldDef, FieldData, m_f->m_fieldMetaDataMap);

        // output error for regression checking.
        if (vm.count("error"))
        {
            int rank = m_f->m_session->GetComm()->GetRank();

            for (int j = 0; j < m_f->m_exp.size(); ++j)
            {
                if (m_f->m_exp[j]->GetPhysState() == false)
                {
                    m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(),
                                            m_f->m_exp[j]->UpdatePhys());
                }

                NekDouble l2err = m_f->m_exp[j]->L2(m_f->m_exp[j]->GetPhys());

                NekDouble linferr =
                    m_f->m_exp[j]->Linf(m_f->m_exp[j]->GetPhys());
                if (rank == 0)
                {
                    cout << "L 2 error (variable "
                         << m_f->m_variables[j] << ") : " << l2err
                         << endl;

                    cout << "L inf error (variable "
                         << m_f->m_variables[j] << ") : " << linferr
                         << endl;
                }
            }
        }
    }
}

void OutputFld::OutputFromData(po::variables_map &vm)
{
    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();
    // Set up communicator and FieldIO object.
    LibUtilities::CommSharedPtr c = m_f->m_session ? m_f->m_session->GetComm() :
        LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
    LibUtilities::FieldIOSharedPtr fld =
        LibUtilities::GetFieldIOFactory().CreateInstance(
            m_config["format"].as<string>(), c, true);

    fld->Write(filename, m_f->m_fielddef, m_f->m_data,
                   m_f->m_fieldMetaDataMap);
}

}
}
