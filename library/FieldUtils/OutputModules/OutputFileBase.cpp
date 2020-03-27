////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFileBase.cpp
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
//  Description: Base class for outputting to a file
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputFileBase.h"
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <boost/format.hpp>
#include <iomanip>

namespace Nektar
{
namespace FieldUtils
{

OutputFileBase::OutputFileBase(FieldSharedPtr f) : OutputModule(f)
{
    m_requireEquiSpaced            = false;
    m_config["writemultiplefiles"] = ConfigOption(
        true, "0",
        "Write multiple files in parallel or when using nparts option");
}

OutputFileBase::~OutputFileBase()
{
}

void OutputFileBase::Process(po::variables_map &vm)
{
    string filename = m_config["outfile"].as<string>();

    if (m_f->m_fieldPts != LibUtilities::NullPtsField)
    {
        ASSERTL0(!m_f->m_writeBndFld, "Boundary can't be obtained from pts.");
        if (WriteFile(filename, vm))
        {
            OutputFromPts(vm);

            if (vm.count("error"))
            {
                PrintErrorFromPts();
            }
        }
    }
    else if (m_f->m_exp.size())
    {
        // reset expansion definition to use equispaced points if required.
        if (m_requireEquiSpaced && (vm.count("noequispaced") == 0) &&
            m_f->m_exp[0]->GetNumElmts() != 0)
        {
            ConvertExpToEquispaced(vm);
        }

        if (m_f->m_writeBndFld)
        {
            if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
            {
                cout << "\t" << GetModuleName()
                     << ": Writing boundary file(s): ";
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

            int nfields = m_f->m_variables.size();
            int normdim = m_f->m_graph->GetMeshDimension();

            // Prepare for normals output
            if (m_f->m_addNormals)
            {
                // Prepare for creating expansions for normals
                m_f->m_exp.resize(nfields + normdim);
                ;

                // Include normal name in m_variables
                string normstr[3] = {"Norm_x", "Norm_y", "Norm_z"};
                for (int j = 0; j < normdim; ++j)
                {
                    m_f->m_exp[nfields + j] =
                        m_f->AppendExpList(m_f->m_numHomogeneousDir);
                    m_f->m_variables.push_back(normstr[j]);
                }
            }

            // Move m_exp to a new expansion vector
            vector<MultiRegions::ExpListSharedPtr> exp(m_f->m_exp.size());
            exp.swap(m_f->m_exp);

            Array<OneD, Array<OneD, const MultiRegions::ExpListSharedPtr>>
                BndExp(exp.size());
            for (int i = 0; i < exp.size(); ++i)
            {
                BndExp[i] = exp[i]->GetBndCondExpansions();
            }

            // get hold of partition boundary regions so we can match it to
            // desired region extraction
            SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                                   exp[0]->GetGraph());
            const SpatialDomains::BoundaryRegionCollection bregions =
                bcs.GetBoundaryRegions();
            map<int, int> BndRegionMap;
            map<int, LibUtilities::CommSharedPtr> BndRegionComm;
            int cnt = 0;
            for (auto &breg_it : bregions)
            {
                BndRegionMap[breg_it.first] = cnt++;
                BndRegionComm[breg_it.first] =
                    bcs.GetBoundaryCommunicators()[breg_it.first];
            }

            // find ending of output file and insert _b1, _b2
            int dot     = filename.find_last_of('.') + 1;
            string ext  = filename.substr(dot, filename.length() - dot);
            string name = filename.substr(0, dot - 1);

            // Store temporary communicator
            LibUtilities::CommSharedPtr tmpComm = m_f->m_comm;

            for (int i = 0; i < m_f->m_bndRegionsToWrite.size(); ++i)
            {
                string outname =
                    name + "_b" +
                    boost::lexical_cast<string>(m_f->m_bndRegionsToWrite[i]) +
                    "." + ext;

                if (!WriteFile(outname, vm))
                {
                    continue;
                }
                RegisterConfig("outfile", outname);

                if (BndRegionMap.count(m_f->m_bndRegionsToWrite[i]) == 1)
                {
                    m_f->m_comm = BndRegionComm[m_f->m_bndRegionsToWrite[i]];

                    int Border = BndRegionMap[m_f->m_bndRegionsToWrite[i]];

                    for (int j = 0; j < exp.size(); ++j)
                    {
                        m_f->m_exp[j] = BndExp[j][Border];
                        m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(),
                                                m_f->m_exp[j]->UpdatePhys());
                    }

                    if (m_f->m_addNormals)
                    {
                        // Get normals
                        Array<OneD, Array<OneD, NekDouble>> NormPhys;
                        exp[0]->GetBoundaryNormals(Border, NormPhys);

                        // add normal coefficients to expansions
                        for (int j = 0; j < normdim; ++j)
                        {
                            m_f->m_exp[nfields + j] =
                                BndExp[nfields + j][Border];
                            Vmath::Vcopy(
                                m_f->m_exp[nfields + j]->GetTotPoints(),
                                NormPhys[j], 1,
                                m_f->m_exp[nfields + j]->UpdatePhys(), 1);
                            m_f->m_exp[nfields + j]->FwdTrans_IterPerExp(
                                m_f->m_exp[nfields + j]->GetPhys(),
                                m_f->m_exp[nfields + j]->UpdateCoeffs());
                        }
                    }
                    OutputFromExp(vm);
                    // output error for regression checking.
                    if (vm.count("error"))
                    {
                        PrintErrorFromExp();
                    }

                    // Reset communicator
                    m_f->m_comm = tmpComm;
                }

                // put outfile back to filename in case of nparts option
                RegisterConfig("outfile", filename);
            }
            // Restore m_exp
            exp.swap(m_f->m_exp);
        }
        else
        {
            if (WriteFile(filename, vm))
            {
                OutputFromExp(vm);
                // output error for regression checking.
                if (vm.count("error"))
                {
                    PrintErrorFromExp();
                }
            }
        }
    }
    else if (m_f->m_data.size())
    {
        ASSERTL0(!m_f->m_writeBndFld, "Boundary extraction requires xml file.");
        if (WriteFile(filename, vm))
        {
            OutputFromData(vm);
        }
    }
}

bool OutputFileBase::WriteFile(std::string &filename, po::variables_map &vm)
{
    // Get path to file. If procid was defined, get the full name
    //     to avoid checking files from other partitions
    fs::path outFile;
    if (vm.count("nparts"))
    {
        outFile = GetFullOutName(filename, vm);
    }
    else
    {
        outFile = GetPath(filename, vm);
    }

    LibUtilities::CommSharedPtr comm;
    if (m_f->m_comm)
    {
        comm = m_f->m_comm;
    }
    else
    {
        comm = LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
    }

    int count = fs::exists(outFile) ? 1 : 0;
    comm->AllReduce(count, LibUtilities::ReduceSum);

    int writeFile = 1;
    if (count && (vm.count("forceoutput") == 0))
    {
        if (vm.count("nparts") == 0) // do not do check if --nparts is enabled.
        {

            writeFile = 0; // set to zero for reduce all to be correct.

            if (comm->TreatAsRankZero())
            {
                string answer;
                cout << "Did you wish to overwrite " << outFile << " (y/n)? ";
                getline(cin, answer);
                if (answer.compare("y") == 0)
                {
                    writeFile = 1;
                }
                else
                {
                    cout << "Not writing file " << filename
                         << " because it already exists" << endl;
                }
            }
            comm->AllReduce(writeFile, LibUtilities::ReduceSum);
        }
    }
    return (writeFile == 0) ? false : true;
}

void OutputFileBase::ConvertExpToEquispaced(po::variables_map &vm)
{
    // Information to create new expansion
    int numFields   = m_f->m_exp.size();
    m_f->m_fielddef = m_f->m_exp[0]->GetFieldDefinitions();

    // Set points to equispaced
    int nPointsNew = 0;
    if (vm.count("output-points"))
    {
        nPointsNew = vm["output-points"].as<int>();
    }
    m_f->m_graph->SetExpansionsToEvenlySpacedPoints(nPointsNew);

    // Save original expansion
    vector<MultiRegions::ExpListSharedPtr> expOld = m_f->m_exp;
    // Create new expansion
    m_f->m_exp[0] = m_f->SetUpFirstExpList(m_f->m_numHomogeneousDir, true);
    for (int i = 1; i < numFields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    }
    // Extract result to new expansion
    for (int i = 0; i < numFields; ++i)
    {
        m_f->m_exp[i]->ExtractCoeffsToCoeffs(expOld[i], expOld[i]->GetCoeffs(),
                                             m_f->m_exp[i]->UpdateCoeffs());
        m_f->m_exp[i]->BwdTrans(m_f->m_exp[i]->GetCoeffs(),
                                m_f->m_exp[i]->UpdatePhys());
    }
    // Extract boundary expansion if needed
    if (m_f->m_writeBndFld)
    {
        Array<OneD, const MultiRegions::ExpListSharedPtr> BndExpOld;
        MultiRegions::ExpListSharedPtr BndExp;
        for (int i = 0; i < numFields; ++i)
        {
            BndExpOld = expOld[i]->GetBndCondExpansions();
            for (int j = 0; j < BndExpOld.size(); ++j)
            {
                BndExp = m_f->m_exp[i]->UpdateBndCondExpansion(j);

                BndExp->ExtractCoeffsToCoeffs(BndExpOld[j],
                                              BndExpOld[j]->GetCoeffs(),
                                              BndExp->UpdateCoeffs());
            }
        }
    }
}

void OutputFileBase::PrintErrorFromPts()
{
    int coordim             = m_f->m_fieldPts->GetDim();
    std::string coordVars[] = {"x", "y", "z"};

    vector<string> variables = m_f->m_variables;
    variables.insert(variables.begin(), coordVars, coordVars + coordim);
    // Get fields and coordinates
    Array<OneD, Array<OneD, NekDouble>> fields(variables.size());

    // We can just grab everything from points. This should be a
    // reference, not a copy.
    m_f->m_fieldPts->GetPts(fields);
    for (int i = 0; i < fields.size(); ++i)
    {
        // calculate L2 and Linf value
        int npts = fields[i].size();

        NekDouble l2err   = 0.0;
        NekDouble linferr = 0.0;
        for (int j = 0; j < npts; ++j)
        {
            l2err += fields[i][j] * fields[i][j];
            linferr = max(linferr, fabs(fields[i][j]));
        }

        m_f->m_comm->AllReduce(l2err, LibUtilities::ReduceSum);
        m_f->m_comm->AllReduce(npts, LibUtilities::ReduceSum);
        m_f->m_comm->AllReduce(linferr, LibUtilities::ReduceMax);

        l2err /= npts;
        l2err = sqrt(l2err);

        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "L 2 error (variable " << variables[i] << ") : " << l2err
                 << endl;

            cout << "L inf error (variable " << variables[i]
                 << ") : " << linferr << endl;
        }
    }
}

void OutputFileBase::PrintErrorFromExp()
{
    int coordim =
        m_f->m_exp[0]->GetExp(0)->GetCoordim() + m_f->m_numHomogeneousDir;
    int totpoints           = m_f->m_exp[0]->GetTotPoints();
    std::string coordVars[] = {"x", "y", "z"};

    // Set up storage for coordinates
    Array<OneD, Array<OneD, NekDouble>> coords(coordim);
    for (int i = 0; i < coordim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(totpoints);
    }

    // Get coordinates
    if (coordim == 1)
    {
        m_f->m_exp[0]->GetCoords(coords[0]);
    }
    else if (coordim == 2)
    {
        m_f->m_exp[0]->GetCoords(coords[0], coords[1]);
    }
    else
    {
        m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);
    }

    for (int j = 0; j < coordim; ++j)
    {
        NekDouble l2err   = m_f->m_exp[0]->L2(coords[j]);
        NekDouble linferr = m_f->m_exp[0]->Linf(coords[j]);

        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "L 2 error (variable " << coordVars[j] << ") : " << l2err
                 << endl;

            cout << "L inf error (variable " << coordVars[j]
                 << ") : " << linferr << endl;
        }
    }

    for (int j = 0; j < m_f->m_exp.size(); ++j)
    {
        NekDouble l2err   = m_f->m_exp[j]->L2(m_f->m_exp[j]->GetPhys());
        NekDouble linferr = m_f->m_exp[j]->Linf(m_f->m_exp[j]->GetPhys());

        if (m_f->m_comm->TreatAsRankZero() && m_f->m_variables.size() > 0)
        {
            cout << "L 2 error (variable " << m_f->m_variables[j]
                 << ") : " << l2err << endl;

            cout << "L inf error (variable " << m_f->m_variables[j]
                 << ") : " << linferr << endl;
        }
    }
}

} // namespace FieldUtils
} // namespace Nektar
