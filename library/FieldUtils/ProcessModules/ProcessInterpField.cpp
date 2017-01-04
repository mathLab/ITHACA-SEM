////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInterpField.cpp
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
//  Description: Interpolate one field to another.
//
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>
using namespace std;

#include "ProcessInterpField.h"

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessInterpField::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "interpfield"),
        ProcessInterpField::create,
        "Interpolates one field to another, requires fromxml, "
        "fromfld to be defined");

ProcessInterpField::ProcessInterpField(FieldSharedPtr f) : ProcessModule(f)
{

    m_config["fromxml"] = ConfigOption(
        false, "NotSet", "Xml file form which to interpolate field");
    m_config["fromfld"] = ConfigOption(
        false, "NotSet", "Fld file form which to interpolate field");

    m_config["clamptolowervalue"] =
        ConfigOption(false, "-10000000", "Lower bound for interpolation value");
    m_config["clamptouppervalue"] =
        ConfigOption(false, "10000000", "Upper bound for interpolation value");
    m_config["defaultvalue"] =
        ConfigOption(false, "0", "Default value if point is outside domain");
}

ProcessInterpField::~ProcessInterpField()
{
}

void ProcessInterpField::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "ProcessInterpField: Interpolating field..." << endl;
        }
    }

    m_fromField = boost::shared_ptr<Field>(new Field());

    std::vector<std::string> files;

    // set up session file for from field
    ParseUtils::GenerateOrderedStringVector(m_config["fromxml"].as<string>().c_str(), files);
    m_fromField->m_session =
        LibUtilities::SessionReader::CreateInstance(0, 0, files);

    // Set up range based on min and max of local parallel partition
    SpatialDomains::DomainRangeShPtr rng =
        MemoryManager<SpatialDomains::DomainRange>::AllocateSharedPtr();

    int coordim = m_f->m_exp[0]->GetCoordim(0);
    int npts    = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > coords(3);

    for (int i = 0; i < coordim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(npts);
    }

    for (int i = coordim; i < 3; ++i)
    {
        coords[i] = NullNekDouble1DArray;
    }

    m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);

    rng->m_checkShape = false;
    switch (coordim)
    {
        case 3:
            rng->m_doZrange = true;
            rng->m_zmin     = Vmath::Vmin(npts, coords[2], 1);
            rng->m_zmax     = Vmath::Vmax(npts, coords[2], 1);
        case 2:
            rng->m_doYrange = true;
            rng->m_ymin     = Vmath::Vmin(npts, coords[1], 1);
            rng->m_ymax     = Vmath::Vmax(npts, coords[1], 1);
        case 1:
            rng->m_doXrange = true;
            rng->m_xmin     = Vmath::Vmin(npts, coords[0], 1);
            rng->m_xmax     = Vmath::Vmax(npts, coords[0], 1);
            break;
        default:
            ASSERTL0(false, "too many values specfied in range");
    }

    // setup rng parameters.
    m_fromField->m_graph =
        SpatialDomains::MeshGraph::Read(m_fromField->m_session, rng);

    // Read in local from field partitions
    const SpatialDomains::ExpansionMap &expansions =
        m_fromField->m_graph->GetExpansions();

    // check for case where no elements are specified on this
    // parallel partition
    if (!expansions.size())
    {
        return;
    }

    Array<OneD, int> ElementGIDs(expansions.size());
    SpatialDomains::ExpansionMap::const_iterator expIt;

    int i = 0;
    for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
    {
        ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
    }

    string fromfld = m_config["fromfld"].as<string>();
    m_f->FieldIOForFile(fromfld)->Import(
        fromfld, m_fromField->m_fielddef, m_fromField->m_data,
        LibUtilities::NullFieldMetaDataMap, ElementGIDs);

    int NumHomogeneousDir = m_fromField->m_fielddef[0]->m_numHomogeneousDir;

    //----------------------------------------------
    // Set up Expansion information to use mode order from field
    m_fromField->m_graph->SetExpansions(m_fromField->m_fielddef);

    int nfields = m_fromField->m_fielddef[0]->m_fields.size();

    m_fromField->m_exp.resize(nfields);
    m_fromField->m_exp[0] =
        m_fromField->SetUpFirstExpList(NumHomogeneousDir, true);

    m_f->m_exp.resize(nfields);

    // declare auxiliary fields.
    for (i = 1; i < nfields; ++i)
    {
        m_f->m_exp[i]         = m_f->AppendExpList(NumHomogeneousDir);
        m_fromField->m_exp[i] = m_fromField->AppendExpList(NumHomogeneousDir);
    }

    // load field into expansion in fromfield.
    for (int j = 0; j < nfields; ++j)
    {
        for (i = 0; i < m_fromField->m_fielddef.size(); i++)
        {
            m_fromField->m_exp[j]->ExtractDataToCoeffs(
                m_fromField->m_fielddef[i], m_fromField->m_data[i],
                m_fromField->m_fielddef[0]->m_fields[j],
                m_fromField->m_exp[j]->UpdateCoeffs());
        }
        m_fromField->m_exp[j]->BwdTrans(m_fromField->m_exp[j]->GetCoeffs(),
                                        m_fromField->m_exp[j]->UpdatePhys());
    }

    int nq1 = m_f->m_exp[0]->GetTotPoints();

    Array<OneD, NekDouble> x1(nq1);
    Array<OneD, NekDouble> y1(nq1);
    Array<OneD, NekDouble> z1(nq1);

    if (coordim == 2)
    {
        m_f->m_exp[0]->GetCoords(x1, y1);
    }
    else if (coordim == 3)
    {
        m_f->m_exp[0]->GetCoords(x1, y1, z1);
    }

    NekDouble clamp_low = m_config["clamptolowervalue"].as<NekDouble>();
    NekDouble clamp_up  = m_config["clamptouppervalue"].as<NekDouble>();
    NekDouble def_value = m_config["defaultvalue"].as<NekDouble>();

    for (int i = 0; i < nfields; i++)
    {
        for (int j = 0; j < nq1; ++j)
        {
            m_f->m_exp[i]->UpdatePhys()[j] = def_value;
        }
    }

    Interpolator interp;
    if (m_f->m_comm->GetRank() == 0)
    {
        interp.SetProgressCallback(&ProcessInterpField::PrintProgressbar, this);
    }
    interp.Interpolate(m_fromField->m_exp, m_f->m_exp);
    if (m_f->m_comm->GetRank() == 0)
    {
        cout << endl;
    }

    for (int i = 0; i < nfields; ++i)
    {
        for (int j = 0; j < nq1; ++j)
        {
            if (m_f->m_exp[i]->GetPhys()[j] > clamp_up)
            {
                m_f->m_exp[i]->UpdatePhys()[j] = clamp_up;
            }
            else if (m_f->m_exp[i]->GetPhys()[j] < clamp_low)
            {
                m_f->m_exp[i]->UpdatePhys()[j] = clamp_low;
            }
        }
    }

    // put field into field data for output
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (int j = 0; j < nfields; ++j)
    {
        m_f->m_exp[j]->FwdTrans(m_f->m_exp[j]->GetPhys(),
                                m_f->m_exp[j]->UpdateCoeffs());
        for (i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back(
                m_fromField->m_fielddef[0]->m_fields[j]);
            m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

void ProcessInterpField::PrintProgressbar(const int position,
                                          const int goal) const
{
    LibUtilities::PrintProgressbar(position, goal, "Interpolating");
}
}
}
