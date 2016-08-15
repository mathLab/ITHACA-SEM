///////////////////////////////////////////////////////////////////////////////
//
// File FilterFieldConvert.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Base clase for filters performing operations on samples
//              of the field.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Filters/FilterFieldConvert.h>
#include <boost/program_options.hpp>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterFieldConvert::className =
        GetFilterFactory().RegisterCreatorFunction(
                "FieldConvert", FilterFieldConvert::create);

FilterFieldConvert::FilterFieldConvert(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams)
    : Filter(pSession)
{
    ParamMap::const_iterator it;

    // OutputFile
    it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        std::stringstream outname;
        outname << m_session->GetSessionName() << ".fld";
        m_outputFile = outname.str();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        if ( it->second.find_last_of('.') != string::npos)
        {
            m_outputFile = it->second;
        }
        else
        {
            std::stringstream outname;
            outname << it->second << ".fld";
            m_outputFile = outname.str();
        }
    }

    // Restart file
    it = pParams.find("RestartFile");
    if (it == pParams.end())
    {
        m_restartFile = "";
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'RestartFile'.");
        if ( it->second.find_last_of('.') != string::npos)
        {
            m_restartFile = it->second;
        }
        else
        {
            std::stringstream outname;
            outname << it->second << ".fld";
            m_restartFile = outname.str();
        }
    }

    // SampleFrequency
    it = pParams.find("SampleFrequency");
    if (it == pParams.end())
    {
        m_sampleFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_sampleFrequency = floor(equ.Evaluate());
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = m_session->GetParameter("NumSteps");
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_outputFrequency = floor(equ.Evaluate());
    }

    m_numSamples  = 0;
    m_index       = 0;
    m_outputIndex = 0;

    //
    // FieldConvert modules
    //
    m_f = boost::shared_ptr<Field>(new Field());
    vector<string>          modcmds;
    // Process modules
    std::stringstream moduleStream;
    it = pParams.find("Modules");
    if (it != pParams.end())
    {
        moduleStream.str(it->second);
    }
    while (!moduleStream.fail())
    {
        std::string sMod;
        moduleStream >> sMod;
        if (!moduleStream.fail())
        {
            modcmds.push_back(sMod);
        }
    }
    // Output module
    modcmds.push_back(m_outputFile);
    // Create modules
    CreateModules(modcmds);
}

FilterFieldConvert::~FilterFieldConvert()
{
}

void FilterFieldConvert::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    v_FillVariablesName(pFields);

    int ncoeff = pFields[0]->GetNcoeffs();
    // m_variables need to be filled by a derived class
    m_outFields.resize(m_variables.size());

    for (int n = 0; n < m_variables.size(); ++n)
    {
        m_outFields[n] = Array<OneD, NekDouble>(ncoeff, 0.0);
    }

    m_fieldMetaData["InitialTime"] = boost::lexical_cast<std::string>(time);

    // Fill some parameters of m_f
    m_f->m_session = m_session;
    m_f->m_graph = pFields[0]->GetGraph();
    m_f->m_comm = m_f->m_session->GetComm();

    // Load restart file if necessary
    if (m_restartFile != "")
    {
        // Load file
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
        std::vector<std::vector<NekDouble> > fieldData;
        LibUtilities::FieldMetaDataMap fieldMetaData;
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateForFile(m_session, m_restartFile);
        fld->Import(m_restartFile, fieldDef, fieldData, fieldMetaData);

        // Extract fields to output
        for (int j = 0; j < m_variables.size(); ++j)
        {
            for (int i = 0; i < fieldData.size(); ++i)
            {
                pFields[0]->ExtractDataToCoeffs(
                    fieldDef[i],
                    fieldData[i],
                    m_variables[j],
                    m_outFields[j]);
            }
        }

        // Load information for numSamples
        if (fieldMetaData.count("NumberOfFieldDumps"))
        {
            m_numSamples = atoi(fieldMetaData["NumberOfFieldDumps"].c_str());
        }
        else
        {
            m_numSamples = 1;
        }

        // Divide by scale
        NekDouble scale = v_GetScale();
        for (int n = 0; n < m_outFields.size(); ++n)
        {
            Vmath::Smul(m_outFields[n].num_elements(),
                        1.0/scale,
                        m_outFields[n],
                        1,
                        m_outFields[n],
                        1);
        }
    }
}

void FilterFieldConvert::v_FillVariablesName(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    int nfield = pFields.num_elements();
    m_variables.resize(pFields.num_elements());
    for (int n = 0; n < nfield; ++n)
    {
        m_variables[n] = pFields[n]->GetSession()->GetVariable(n);
    }
}

void FilterFieldConvert::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_index++;
    if (m_index % m_sampleFrequency > 0)
    {
        return;
    }

    m_numSamples++;
    v_ProcessSample(pFields, time);

    if (m_index % m_outputFrequency == 0)
    {
        m_fieldMetaData["FinalTime"] = boost::lexical_cast<std::string>(time);
        v_PrepareOutput(pFields, time);
        OutputField(pFields, ++m_outputIndex);
    }
}

void FilterFieldConvert::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_fieldMetaData["FinalTime"] = boost::lexical_cast<std::string>(time);
    v_PrepareOutput(pFields, time);
    OutputField(pFields);
}

void FilterFieldConvert::v_PrepareOutput(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    for(int n = 0; n < pFields.num_elements(); ++n)
    {
        Vmath::Vcopy(m_outFields[n].num_elements(),
                    pFields[n]->GetCoeffs(),
                    1,
                    m_outFields[n],
                    1);
    }
}

void FilterFieldConvert::OutputField(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, int dump)
{
    NekDouble scale = v_GetScale();
    for (int n = 0; n < m_outFields.size(); ++n)
    {
        Vmath::Smul(m_outFields[n].num_elements(),
                    scale,
                    m_outFields[n],
                    1,
                    m_outFields[n],
                    1);
    }

    CreateFields(pFields);

    // Determine new file name
    std::stringstream outname;
    int    dot    = m_outputFile.find_last_of('.');
    string name   = m_outputFile.substr(0, dot);
    string ext    = m_outputFile.substr(dot, m_outputFile.length() - dot);
    std::string suffix = v_GetFileSuffix();
    if (dump == -1) // final dump
    {
        outname << name << suffix << ext;
    }
    else
    {
        outname << name << "_" << dump << suffix << ext;
    }
    m_modules[m_modules.size()-1]->RegisterConfig("outfile", outname.str());

    // Prevent checking before overwritting
    po::options_description desc("Available options");
        desc.add_options()
            ("forceoutput,f",
                "Force the output to be written without any checks");
    po::variables_map vm;
    vm.insert(std::make_pair("forceoutput", po::variable_value()));
    // Run field process.
    for (int i = 0; i < m_modules.size(); ++i)
    {
        m_modules[i]->Process(vm);
        cout.flush();
    }

    // Empty m_f to save memory
    ClearFields();

    if (dump != -1) // not final dump so rescale
    {
        for (int n = 0; n < m_outFields.size(); ++n)
        {
            Vmath::Smul(m_outFields[n].num_elements(),
                        1.0 / scale,
                        m_outFields[n],
                        1,
                        m_outFields[n],
                        1);
        }
    }
}

bool FilterFieldConvert::v_IsTimeDependent()
{
    return true;
}

void FilterFieldConvert::CreateModules( vector<string> &modcmds)
{
    for (int i = 0; i < modcmds.size(); ++i)
    {
        // First split each command by the colon separator.
        vector<string> tmp1;
        ModuleKey module;
        int offset = 1;

        boost::split(tmp1, modcmds[i], boost::is_any_of(":"));

        if (i == modcmds.size() - 1)
        {
            module.first = eOutputModule;

            // If no colon detected, automatically detect mesh type from
            // file extension. Otherwise override and use tmp1[1] as the
            // module to load. This also allows us to pass options to
            // input/output modules. So, for example, to override
            // filename.xml to be read as vtk, you use:
            //
            // filename.xml:vtk:opt1=arg1:opt2=arg2
            if (tmp1.size() == 1)
            {
                int    dot    = tmp1[0].find_last_of('.') + 1;
                string ext    = tmp1[0].substr(dot, tmp1[0].length() - dot);

                module.second = ext;
                tmp1.push_back(string("outfile=") + tmp1[0]);
            }
            else
            {
                module.second = tmp1[1];
                tmp1.push_back(string("outfile=") + tmp1[0]);
                offset++;
            }
        }
        else
        {
            module.first  = eProcessModule;
            module.second = tmp1[0];
        }

        // Create modules
        ModuleSharedPtr mod;
        mod = GetModuleFactory().CreateInstance(module, m_f);
        m_modules.push_back(mod);

        // Set options for this module.
        for (int j = offset; j < tmp1.size(); ++j)
        {
            vector<string> tmp2;
            boost::split(tmp2, tmp1[j], boost::is_any_of("="));

            if (tmp2.size() == 1)
            {
                mod->RegisterConfig(tmp2[0], "1");
            }
            else if (tmp2.size() == 2)
            {
                mod->RegisterConfig(tmp2[0], tmp2[1]);
            }
            else
            {
                cerr << "ERROR: Invalid module configuration: format is "
                     << "either :arg or :arg=val" << endl;
                abort();
            }
        }

        // Ensure configuration options have been set.
        mod->SetDefaults();
    }

    bool RequiresEquiSpaced = false;
    for (int i = 0; i < m_modules.size(); ++i)
    {
        if(m_modules[i]->GetRequireEquiSpaced())
        {
            RequiresEquiSpaced = true;
        }
    }
    if (RequiresEquiSpaced)
    {
        for (int i = 0; i < m_modules.size(); ++i)
        {
            m_modules[i]->SetRequireEquiSpaced(true);
        }
    }
}

void FilterFieldConvert::CreateFields(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    m_f->m_fieldMetaDataMap = m_fieldMetaData;
    m_f->m_fieldPts = LibUtilities::NullPtsField;
    // Create m_f->m_exp
    int NumHomogeneousDir = 0;
    if (pFields[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        NumHomogeneousDir = 1;
    }
    else if (pFields[0]->GetExpType() == MultiRegions::e3DH2D)
    {
        NumHomogeneousDir = 2;
    }

    m_f->m_exp.resize(m_variables.size());
    m_f->m_exp[0] = pFields[0];
    for (int n = 0; n < m_variables.size(); ++n)
    {
        m_f->m_exp[n] = m_f->AppendExpList(
                            NumHomogeneousDir, m_variables[0]);
        m_f->m_exp[n]->SetWaveSpace(false);
        Vmath::Vcopy( m_outFields[n].num_elements(),
                      m_outFields[n], 1,
                      m_f->m_exp[n]->UpdateCoeffs(), 1);
        m_f->m_exp[n]->BwdTrans( m_f->m_exp[n]->GetCoeffs(),
                                 m_f->m_exp[n]->UpdatePhys());
    }
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        pFields[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    for (int n = 0; n < m_outFields.size(); ++n)
    {
        for (int i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back(m_variables[n]);
            m_f->m_exp[n]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }
    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

void FilterFieldConvert::ClearFields()
{
    m_f->m_fieldPts = LibUtilities::NullPtsField;
    m_f->m_exp.clear();
    m_f->m_fielddef = std::vector<LibUtilities::FieldDefinitionsSharedPtr>();
    m_f->m_data = std::vector<std::vector<NekDouble> > ();
}

}
}
