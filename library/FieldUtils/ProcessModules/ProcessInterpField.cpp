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

#include <boost/core/ignore_unused.hpp>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessInterpField.h"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

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
        false, "NotSet", "Xml file from which to interpolate field");
    m_config["fromfld"] = ConfigOption(
        false, "NotSet", "Fld file from which to interpolate field");

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
    boost::ignore_unused(vm);

    FieldSharedPtr fromField = std::shared_ptr<Field>(new Field());

    std::vector<std::string> files;

    // set up session file for from field
    char *argv[] = { const_cast<char *>("FieldConvert"), nullptr };
    ParseUtils::GenerateVector(m_config["fromxml"].as<string>(), files);
    fromField->m_session =
        LibUtilities::SessionReader::CreateInstance(
            1, argv, files,
            LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0));

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
            /* Falls through. */
        case 2:
            rng->m_doYrange = true;
            rng->m_ymin     = Vmath::Vmin(npts, coords[1], 1);
            rng->m_ymax     = Vmath::Vmax(npts, coords[1], 1);
            /* Falls through. */
        case 1:
            rng->m_doXrange = true;
            rng->m_xmin     = Vmath::Vmin(npts, coords[0], 1);
            rng->m_xmax     = Vmath::Vmax(npts, coords[0], 1);
            break;
        default:
            NEKERROR(ErrorUtil::efatal, "coordim should be <= 3");
    }

    // setup rng parameters.
    fromField->m_graph =
        SpatialDomains::MeshGraph::Read(fromField->m_session, rng);

    // Read in local from field partitions
    const SpatialDomains::ExpansionMap &expansions =
        fromField->m_graph->GetExpansions();

    // check for case where no elements are specified on this
    // parallel partition
    if (!expansions.size())
    {
        return;
    }

    Array<OneD, int> ElementGIDs(expansions.size());

    int i = 0;
    for (auto &expIt : expansions)
    {
        ElementGIDs[i++] = expIt.second->m_geomShPtr->GetGlobalID();
    }

    string fromfld = m_config["fromfld"].as<string>();
    m_f->FieldIOForFile(fromfld)->Import(
        fromfld, fromField->m_fielddef, fromField->m_data,
        LibUtilities::NullFieldMetaDataMap, ElementGIDs);

    int NumHomogeneousDir = fromField->m_fielddef[0]->m_numHomogeneousDir;

    //----------------------------------------------
    // Set up Expansion information to use mode order from field
    fromField->m_graph->SetExpansions(fromField->m_fielddef);

    int nfields = fromField->m_fielddef[0]->m_fields.size();

    fromField->m_exp.resize(nfields);
    fromField->m_exp[0] =
        fromField->SetUpFirstExpList(NumHomogeneousDir, true);

    m_f->m_exp.resize(nfields);

    // declare auxiliary fields.
    for (i = 1; i < nfields; ++i)
    {
        m_f->m_exp[i]       = m_f->AppendExpList(NumHomogeneousDir);
        fromField->m_exp[i] = fromField->AppendExpList(NumHomogeneousDir);
    }

    // load field into expansion in fromfield.
    for (int j = 0; j < nfields; ++j)
    {
        for (i = 0; i < fromField->m_fielddef.size(); i++)
        {
            fromField->m_exp[j]->ExtractDataToCoeffs(
                fromField->m_fielddef[i], fromField->m_data[i],
                fromField->m_fielddef[0]->m_fields[j],
                fromField->m_exp[j]->UpdateCoeffs());
        }
        fromField->m_exp[j]->BwdTrans(fromField->m_exp[j]->GetCoeffs(),
                                        fromField->m_exp[j]->UpdatePhys());
    }

    int nq1 = m_f->m_exp[0]->GetTotPoints();

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
    if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
    {
        interp.SetProgressCallback(&ProcessInterpField::PrintProgressbar, this);
    }
    interp.Interpolate(fromField->m_exp, m_f->m_exp);
    if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
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
        m_f->m_exp[i]->FwdTrans_IterPerExp(
                    m_f->m_exp[i]->GetPhys(), m_f->m_exp[i]->UpdateCoeffs());
    }
    // save field names
    m_f->m_variables = fromField->m_fielddef[0]->m_fields;
}

void ProcessInterpField::PrintProgressbar(const int position,
                                          const int goal) const
{
    LibUtilities::PrintProgressbar(position, goal, "Interpolating");
}
}
}
