////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInterpPoints.cpp
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
//  Description: Interpolate  field to a series of specified points.
//
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
using namespace std;

#include "ProcessInterpPoints.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessInterpPoints::className =
    GetModuleFactory().RegisterCreatorFunction(
       ModuleKey(eProcessModule, "interppoints"),
       ProcessInterpPoints::create,
       "Interpolates a set of points to another, requires fromfld and "
       "fromxml to be defined, a line or plane of points can be "
       "defined");


ProcessInterpPoints::ProcessInterpPoints(FieldSharedPtr f) : ProcessModule(f)
{

    m_config["fromxml"] =
            ConfigOption(false, "NotSet",
                         "Xml file from which to interpolate field");

    ASSERTL0(m_config["fromxml"].as<string>().compare("NotSet") != 0,
             "Need to specify fromxml=file.xml");

    m_config["fromfld"] =
            ConfigOption(false,"NotSet",
                         "Fld file from which to interpolate field");

    ASSERTL0(m_config["fromfld"].as<string>().compare("NotSet") != 0,
             "Need to specify fromfld=file.fld ");

    m_config["clamptolowervalue"] =
            ConfigOption(false, "-10000000",
                         "Lower bound for interpolation value");
    m_config["clamptouppervalue"] =
            ConfigOption(false, "10000000",
                         "Upper bound for interpolation value");
    m_config["defaultvalue"] =
            ConfigOption(false, "0",
                         "Default value if point is outside domain");
    m_config["line"] =
            ConfigOption(false, "NotSet",
                         "Specify a line of N points using "
                         "line=N,x0,y0,z0,z1,y1,z1");

    m_config["plane"] =
            ConfigOption(false, "NotSet",
                         "Specify a plane of N1 x N1 points using "
                         "plane=N1,N2,x0,y0,z0,z1,y1,z1,x2,y2,z2,x3,"
                         "y3,z3");
}

ProcessInterpPoints::~ProcessInterpPoints()
{
}

void ProcessInterpPoints::Process(po::variables_map &vm)
{

    if(m_f->m_verbose)
    {
        cout << "Processing point interpolation" << endl;
    }


    // Check for command line point specification if no .pts file specified
    if(m_f->m_fieldPts == LibUtilities::NullPtsField)
    {
        if(m_config["line"].as<string>().compare("NotSet") != 0)
        {
            string help = m_config["line"].as<string>();
            vector<NekDouble> values;
            ASSERTL0(ParseUtils::GenerateUnOrderedVector(
                        m_config["line"].as<string>().c_str(),values),
                     "Failed to interpret line string");

            ASSERTL0(values.size() > 2,
                     "line string should contain 2Dim+1 values "
                     "N,x0,y0,z0,x1,y1,z1");

            double tmp;
            ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N is not an integer");
            ASSERTL0(values[0] > 1, "N is not a valid number");
           
            int dim = (values.size()-1)/2;
            int npts = values[0];
            Array<OneD, Array<OneD, NekDouble> > pts(dim);

            for(int i = 0; i < dim; ++i)
            {
                pts[i] = Array<OneD,NekDouble>(npts);
            }

            for(int i = 0; i < npts; ++i)
            {
                pts[0][i] = values[1]
                        + i/((NekDouble)(npts-1))*(values[dim+1] - values[1]);
                if(dim > 1)
                {
                    pts[1][i] = values[2]
                        + i/((NekDouble)(npts-1))*(values[dim+2] - values[2]);

                    if(dim > 2)
                    {
                        pts[2][i] = values[3]
                        + i/((NekDouble)(npts-1))*(values[dim+3] - values[3]);
                    }
                }
            }

            vector<int> ppe;
            ppe.push_back(npts);
            m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim, pts);
            m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsLine);
            m_f->m_fieldPts->SetPointsPerEdge(ppe);
       
        }
        else if(m_config["plane"].as<string>().compare("NotSet") != 0)
        {
            string help = m_config["plane"].as<string>();
            vector<NekDouble> values;
            ASSERTL0(ParseUtils::GenerateUnOrderedVector(
                             m_config["plane"].as<string>().c_str(),values),
                     "Failed to interpret plane string");

            ASSERTL0(values.size() > 9,
                     "plane string should contain 4Dim+2 values "
                     "N1,N2,x0,y0,x1,y1,x2,y2,x3,y3");

            double tmp;
            ASSERTL0(std::modf(values[0], &tmp) == 0.0, "N1 is not an integer");
            ASSERTL0(std::modf(values[1], &tmp) == 0.0, "N2 is not an integer");
            
            ASSERTL0(values[0] > 1, "N1 is not a valid number");
            ASSERTL0(values[1] > 1, "N2 is not a valid number");
            
            int dim = (values.size()-2)/4;

            int npts1 = values[0];
            int npts2 = values[1];

            Array<OneD, Array<OneD, NekDouble> > pts(dim);

            for(int i = 0; i < dim; ++i)
            {
                pts[i] = Array<OneD,NekDouble>(npts1*npts2);
            }

            for(int j = 0; j < npts2; ++j)
            {
                for(int i = 0; i < npts1; ++i)
                {
                    pts[0][i+j*npts1] =
                        (values[2] + i/((NekDouble)(npts1-1))*(values[dim+2] - values[2]))*(1.0-j/((NekDouble)(npts2-1))) +
                        (values[3*dim+2] + i/((NekDouble)(npts1-1))*(values[2*dim+2] - values[3*dim+2]))*(j/((NekDouble)(npts2-1)));
                    pts[1][i+j*npts1] =
                        (values[3] + i/((NekDouble)(npts1-1))*(values[dim+3] - values[3]))*(1.0-j/((NekDouble)(npts2-1))) +
                        (values[3*dim+3] + i/((NekDouble)(npts1-1))*(values[2*dim+3] - values[3*dim+3]))*(j/((NekDouble)(npts2-1)));
                    
                    if(dim > 2)
                    {
                        pts[2][i+j*npts1] =
                            (values[4] + i/((NekDouble)(npts1-1))*(values[dim+4] - values[4]))*(1.0-j/((NekDouble)(npts2-1))) +
                            (values[3*dim+4] + i/((NekDouble)(npts1-1))*(values[2*dim+4] - values[3*dim+4]))*(j/((NekDouble)(npts2-1)));
                    }
                }
            }

            vector<int> ppe;
            ppe.push_back(npts1);
            ppe.push_back(npts2);
            m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(dim, pts);
            m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsPlane);
            m_f->m_fieldPts->SetPointsPerEdge(ppe);
           
        }
    }


    FieldSharedPtr fromField = boost::shared_ptr<Field>(new Field());

    std::vector<std::string> files;
    // set up session file for from field
    files.push_back(m_config["fromxml"].as<string>());
    fromField->m_session = LibUtilities::SessionReader::
        CreateInstance(0, 0, files);

    // Set up range based on min and max of local parallel partition
    SpatialDomains::DomainRangeShPtr rng = MemoryManager<SpatialDomains::DomainRange>::AllocateSharedPtr();

    int coordim = m_f->m_fieldPts->GetDim();
    int npts    = m_f->m_fieldPts->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > pts;
    m_f->m_fieldPts->GetPts(pts);

    rng->m_checkShape   = false;
    switch(coordim)
    {
    case 3:
        rng->m_doZrange = true;
        rng->m_zmin = Vmath::Vmin(npts, pts[2],1);
        rng->m_zmax = Vmath::Vmax(npts, pts[2],1);
    case 2:
        rng->m_doYrange = true;
        rng->m_ymin = Vmath::Vmin(npts, pts[1],1);
        rng->m_ymax = Vmath::Vmax(npts, pts[1],1);
    case 1:
        rng->m_doXrange = true;
        rng->m_xmin = Vmath::Vmin(npts, pts[0],1);
        rng->m_xmax = Vmath::Vmax(npts, pts[0],1);
        break;
    default:
        ASSERTL0(false,"too many values specfied in range");
    }

    // setup rng parameters.
    fromField->m_graph = SpatialDomains::MeshGraph::Read(fromField->m_session,rng);

    // Read in local from field partitions
    const SpatialDomains::ExpansionMap &expansions = fromField->m_graph->GetExpansions();

    fromField->m_fld = MemoryManager<LibUtilities::FieldIO>
        ::AllocateSharedPtr(fromField->m_session->GetComm());

    Array<OneD,int> ElementGIDs(expansions.size());
    SpatialDomains::ExpansionMap::const_iterator expIt;

    int i = 0;
    for (expIt = expansions.begin(); expIt != expansions.end();
         ++expIt)
    {
        ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
    }

    string fromfld = m_config["fromfld"].as<string>();
    fromField->m_fld->Import(fromfld,fromField->m_fielddef,
                             fromField->m_data,
                             LibUtilities::NullFieldMetaDataMap,
                             ElementGIDs);

    int NumHomogeneousDir = fromField->m_fielddef[0]->m_numHomogeneousDir;

    //----------------------------------------------
    // Set up Expansion information to use mode order from field
    fromField->m_graph->SetExpansions(fromField->m_fielddef);

    int nfields = fromField->m_fielddef[0]->m_fields.size();

    fromField->m_exp.resize(nfields);
    fromField->m_exp[0] = fromField->SetUpFirstExpList(NumHomogeneousDir,true);

    m_f->m_exp.resize(nfields);

    // declare auxiliary fields.
    for(i = 1; i < nfields; ++i)
    {
        fromField->m_exp[i] = fromField->AppendExpList(NumHomogeneousDir);
    }

    // load field into expansion in fromfield.
    for(int j = 0; j < nfields; ++j)
    {
        for (i = 0; i < fromField->m_fielddef.size(); i++)
        {
            fromField->m_exp[j]->ExtractDataToCoeffs(
                                                     fromField->m_fielddef[i],
                                                     fromField->m_data[i],
                                                     fromField->m_fielddef[0]->m_fields[j],
                                                     fromField->m_exp[j]->UpdateCoeffs());
        }
        fromField->m_exp[j]->BwdTrans(fromField->m_exp[j]->GetCoeffs(),
                                      fromField->m_exp[j]->UpdatePhys());

        Array< OneD, NekDouble > newPts(m_f->m_fieldPts->GetNpoints());
        m_f->m_fieldPts->AddField(newPts, fromField->m_fielddef[0]->m_fields[j]);
    }

    if(fromField->m_session->GetComm()->GetRank() == 0)
    {
        cout << "Interpolating [" << flush;
    }

    NekDouble clamp_low = m_config["clamptolowervalue"].as<NekDouble>();
    NekDouble clamp_up  = m_config["clamptouppervalue"].as<NekDouble>();
    NekDouble def_value = m_config["defaultvalue"].as<NekDouble>();

    InterpolateFieldToPts(fromField->m_exp, pts,
                          clamp_low, clamp_up, def_value);

    if(fromField->m_session->GetComm()->GetRank() == 0)
    {
        cout << "]" << endl;
    }

}

void ProcessInterpPoints::InterpolateFieldToPts(
                                                vector<MultiRegions::ExpListSharedPtr> &field0,
                                                Array<OneD, Array<OneD, NekDouble> >   &pts,
                                                NekDouble                              clamp_low,
                                                NekDouble                              clamp_up,
                                                NekDouble                              def_value)
{
    int expdim = pts.num_elements();

    Array<OneD, NekDouble> coords(expdim), Lcoords(expdim);
    int nq1 = pts[0].num_elements();
    int elmtid, offset;
    int r, f;
    int intpts = 0;

    // resize data field
    m_f->m_data.resize(field0.size());

    for (f = 0; f < field0.size(); ++f)
    {
        m_f->m_data[f].resize(nq1);
    }

    for (r = 0; r < nq1; r++)
    {
        coords[0] = pts[0][r];
        coords[1] = pts[1][r];
        if (expdim == 3)
        {
            coords[2] = pts[2][r];
        }

        // Obtain Element and LocalCoordinate to interpolate
        elmtid = field0[0]->GetExpIndex(coords, Lcoords, 1e-3);

        if(elmtid >= 0)
        {
            offset = field0[0]->GetPhys_Offset(field0[0]->
                                               GetOffset_Elmt_Id(elmtid));

            for (f = 0; f < field0.size(); ++f)
            {
                NekDouble value;
                value = field0[f]->GetExp(elmtid)->
                    StdPhysEvaluate(Lcoords, field0[f]->GetPhys() +offset);

                if ((boost::math::isnan)(value))
                {
                    ASSERTL0(false, "new value is not a number");
                }
                else
                {
                    value = (value > clamp_up)? clamp_up :
                        ((value < clamp_low)? clamp_low :
                         value);

                    m_f->m_data[f][r] = value;
                }
            }
        }
        else
        {
            for (f = 0; f < field0.size(); ++f)
            {
                m_f->m_data[f][r] = def_value;
            }
        }

        if (intpts%1000 == 0)
        {
            cout <<"." << flush;
        }
        intpts ++;
    }
}

}
}
