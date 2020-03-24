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
#include <boost/core/ignore_unused.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace Nektar
{
namespace SolverUtils
{
std::string FilterFieldConvert::className =
        GetFilterFactory().RegisterCreatorFunction(
                "FieldConvert", FilterFieldConvert::create);

FilterFieldConvert::FilterFieldConvert(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    m_dt = m_session->GetParameter("TimeStep");

    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        std::stringstream outname;
        outname << m_session->GetSessionName() << ".fld";
        m_outputFile = outname.str();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        if ( it->second.find_last_of('.') != std::string::npos)
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
        if ( it->second.find_last_of('.') != std::string::npos)
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

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = m_session->GetParameter("NumSteps");
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // The base class can use SampleFrequency = OutputFrequency
    //    (Derived classes need to override this if needed)
    m_sampleFrequency = m_outputFrequency;

    // Phase sampling option
    it = pParams.find("PhaseAverage");
    if (it == pParams.end())
    {
        m_phaseSample = false;
    }
    else
    {
        std::string sOption = it->second.c_str();
        m_phaseSample = (boost::iequals(sOption, "true")) ||
                   (boost::iequals(sOption, "yes"));
    }

    if(m_phaseSample)
    {
        auto itPeriod = pParams.find("PhaseAveragePeriod");
        auto itPhase = pParams.find("PhaseAveragePhase");

        // Error if only one of the required params for PhaseAverage is present
        ASSERTL0((itPeriod != pParams.end() && itPhase != pParams.end()),
            "The phase sampling feature requires both 'PhaseAveragePeriod' and "
            "'PhaseAveragePhase' to be set.");

        LibUtilities::Equation equPeriod(
            m_session->GetInterpreter(), itPeriod->second);
        m_phaseSamplePeriod = equPeriod.Evaluate();

        LibUtilities::Equation equPhase(
            m_session->GetInterpreter(), itPhase->second);
        m_phaseSamplePhase = equPhase.Evaluate();

        // Check that phase and period are within required limits
        ASSERTL0(m_phaseSamplePeriod > 0,
                 "PhaseAveragePeriod must be greater than 0.");
        ASSERTL0(m_phaseSamplePhase >= 0 && m_phaseSamplePhase <= 1,
                 "PhaseAveragePhase must be between 0 and 1.");

        // Load sampling frequency, overriding the previous value
        it = pParams.find("SampleFrequency");
        if (it == pParams.end())
        {
            m_sampleFrequency = 1;
        }
        else
        {
            LibUtilities::Equation equ(
                m_session->GetInterpreter(), it->second);
            m_sampleFrequency = round(equ.Evaluate());
        }

        // Compute tolerance within which sampling occurs.
        m_phaseTolerance = m_dt * m_sampleFrequency /
            (m_phaseSamplePeriod * 2);

        // Display worst case scenario sampling tolerance for exact phase, if
        // verbose option is active
        if (m_session->GetComm()->GetRank() == 0 &&
            m_session->DefinesCmdLineArgument("verbose"))
        {
            std::cout << "Phase sampling activated with period "
                      << m_phaseSamplePeriod << " and phase "
                      << m_phaseSamplePhase << "." << std::endl
                      << "Sampling within a tolerance of "
                      << std::setprecision(6) << m_phaseTolerance << "."
                      << std::endl;
        }
    }

    m_numSamples  = 0;
    m_index       = 0;
    m_outputIndex = 0;

    //
    // FieldConvert modules
    //
    m_f = std::shared_ptr<Field>(new Field());
    std::vector<std::string> modcmds;
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
    // Strip options from m_outputFile
    std::vector<std::string> tmp;
    boost::split(tmp, m_outputFile, boost::is_any_of(":"));
    m_outputFile = tmp[0];
}

FilterFieldConvert::~FilterFieldConvert()
{
}

void FilterFieldConvert::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    v_FillVariablesName(pFields);

    // m_variables need to be filled by a derived class
    m_outFields.resize(m_variables.size());
    int nfield;

    for (int n = 0; n < m_variables.size(); ++n)
    {
        // if n >= pFields.size() assum we have used n=0 field
        nfield = (n < pFields.size())? n: 0;

        m_outFields[n] =
                Array<OneD, NekDouble>(pFields[nfield]->GetNcoeffs(), 0.0);
    }

    m_fieldMetaData["InitialTime"] = boost::lexical_cast<std::string>(time);

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
        int nfield = -1, k;
        for (int j = 0; j < m_variables.size(); ++j)
        {
            // see if m_variables is part of pFields definition and if
            // so use that field for extract
            for(k = 0; k < pFields.size(); ++k)
            {
                if(pFields[k]->GetSession()->GetVariable(k)
                   == m_variables[j])
                {
                    nfield = k;
                    break;
                }
            }
            if(nfield == -1)
            {
                nfield = 0;
            }

            for (int i = 0; i < fieldData.size(); ++i)
            {
                pFields[nfield]->ExtractDataToCoeffs(
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

        if(fieldMetaData.count("InitialTime"))
        {
            m_fieldMetaData["InitialTime"] = fieldMetaData["InitialTime"];
        }

        // Divide by scale
        NekDouble scale = v_GetScale();
        for (int n = 0; n < m_outFields.size(); ++n)
        {
            Vmath::Smul(m_outFields[n].size(),
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
    int nfield = pFields.size();
    m_variables.resize(pFields.size());
    for (int n = 0; n < nfield; ++n)
    {
        m_variables[n] = pFields[n]->GetSession()->GetVariable(n);
    }

    // Need to create a dummy coeffs vector to get extra variables names...
    std::vector<Array<OneD, NekDouble> > coeffs(nfield);
    for (int n = 0; n < nfield; ++n)
    {
        coeffs[n] = pFields[n]->GetCoeffs();
    }
    // Get extra variables
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");
    equ->ExtraFldOutput(coeffs, m_variables);
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

    // Append extra fields
    int nfield = pFields.size();
    std::vector<Array<OneD, NekDouble> > coeffs(nfield);
    for (int n = 0; n < nfield; ++n)
    {
        coeffs[n] = pFields[n]->GetCoeffs();
    }
    std::vector<std::string> variables = m_variables;
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");
    equ->ExtraFldOutput(coeffs, variables);

    if(m_phaseSample)
    {
        // The sample is added to the filter only if the current time
        // corresponds to the correct phase. Introducing M as number of
        // cycles and N nondimensional phase (between 0 and 1):
        // t = M * m_phaseSamplePeriod + N * m_phaseSamplePeriod
        int currentCycle       = floor(time / m_phaseSamplePeriod);
        NekDouble currentPhase = time / m_phaseSamplePeriod - currentCycle;

        // Evaluate phase relative to the requested value.
        NekDouble relativePhase = fabs(m_phaseSamplePhase - currentPhase);

        // Check if relative phase is within required tolerance and sample.
        // Care must be taken to handle the cases at phase 0 as the sample might
        // have to be taken at the very end of the previous cycle instead.
        if (relativePhase < m_phaseTolerance ||
            fabs(relativePhase-1) < m_phaseTolerance)
        {
            m_numSamples++;
            v_ProcessSample(pFields, coeffs, time);
            if (m_session->GetComm()->GetRank() == 0 &&
                m_session->DefinesCmdLineArgument("verbose"))
            {
                std::cout << "Sample: " << std::setw(8) << std::left
                          << m_numSamples << "Phase: " << std::setw(8)
                          << std::left << currentPhase << std::endl;
            }
        }
    }
    else
    {
        m_numSamples++;
        v_ProcessSample(pFields, coeffs, time);
    }

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

void FilterFieldConvert::v_ProcessSample(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
          std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);

    for(int n = 0; n < m_outFields.size(); ++n)
    {
        Vmath::Vcopy(m_outFields[n].size(),
                    fieldcoeffs[n],
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
        Vmath::Smul(m_outFields[n].size(),
                    scale,
                    m_outFields[n],
                    1,
                    m_outFields[n],
                    1);
    }

    CreateFields(pFields);

    // Determine new file name
    std::stringstream outname;
    int         dot    = m_outputFile.find_last_of('.');
    std::string name   = m_outputFile.substr(0, dot);
    std::string ext    = m_outputFile.substr(dot, m_outputFile.length() - dot);
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

    // Prevent checking before overwriting
    po::options_description desc("Available options");
        desc.add_options()
            ("forceoutput,f",
                "Force the output to be written without any checks");
    po::variables_map vm;
    vm.insert(std::make_pair("forceoutput", po::variable_value()));

    // Run field process.
    for (int n = 0; n < SIZE_ModulePriority; ++n)
    {
        ModulePriority priority = static_cast<ModulePriority>(n);
        for (int i = 0; i < m_modules.size(); ++i)
        {
            if(m_modules[i]->GetModulePriority() == priority)
            {
                m_modules[i]->Process(vm);
            }
        }
    }

    // Empty m_f to save memory
    m_f->ClearField();

    if (dump != -1) // not final dump so rescale
    {
        for (int n = 0; n < m_outFields.size(); ++n)
        {
            Vmath::Smul(m_outFields[n].size(),
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

void FilterFieldConvert::CreateModules(std::vector<std::string> &modcmds)
{
    for (int i = 0; i < modcmds.size(); ++i)
    {
        // First split each command by the colon separator.
        std::vector<std::string> tmp1;
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
                int         dot = tmp1[0].find_last_of('.') + 1;
                std::string ext = tmp1[0].substr(dot, tmp1[0].length() - dot);

                module.second = ext;
                tmp1.push_back(std::string("outfile=") + tmp1[0]);
            }
            else
            {
                module.second = tmp1[1];
                tmp1.push_back(std::string("outfile=") + tmp1[0]);
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
            std::vector<std::string> tmp2;
            boost::split(tmp2, tmp1[j], boost::is_any_of("="));

            if (tmp2.size() == 1)
            {
                mod->RegisterConfig(tmp2[0]);
            }
            else if (tmp2.size() == 2)
            {
                mod->RegisterConfig(tmp2[0], tmp2[1]);
            }
            else
            {
                std::cerr << "ERROR: Invalid module configuration: format is "
                          << "either :arg or :arg=val" << std::endl;
                abort();
            }
        }

        // Ensure configuration options have been set.
        mod->SetDefaults();
    }

    // Include equispaced output if needed
    Array< OneD, int>  modulesCount(SIZE_ModulePriority,0);
    for (int i = 0; i < m_modules.size(); ++i)
    {
        ++modulesCount[m_modules[i]->GetModulePriority()];
    }
    if( modulesCount[eModifyPts] != 0 &&
        modulesCount[eCreatePts] == 0 &&
        modulesCount[eConvertExpToPts] == 0)
    {
        ModuleKey               module;
        ModuleSharedPtr         mod;
        module.first  = eProcessModule;
        module.second = std::string("equispacedoutput");
        mod = GetModuleFactory().CreateInstance(module, m_f);
        m_modules.insert(m_modules.end()-1, mod);
        mod->SetDefaults();
    }

    // Check if modules provided are compatible
    CheckModules(m_modules);
}

void FilterFieldConvert::CreateFields(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    // Fill some parameters of m_f
    m_f->m_session = m_session;
    m_f->m_graph = pFields[0]->GetGraph();
    m_f->m_comm = m_f->m_session->GetComm();
    m_f->m_fieldMetaDataMap = m_fieldMetaData;
    m_f->m_fieldPts = LibUtilities::NullPtsField;
    // Create m_f->m_exp
    m_f->m_numHomogeneousDir = 0;
    if (pFields[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        m_f->m_numHomogeneousDir = 1;
    }
    else if (pFields[0]->GetExpType() == MultiRegions::e3DH2D)
    {
        m_f->m_numHomogeneousDir = 2;
    }

    m_f->m_exp.resize(m_variables.size());
    m_f->m_exp[0] = pFields[0];
    int nfield;
    for (int n = 0; n < m_variables.size(); ++n)
    {
        // if n >= pFields.size() assume we have used n=0 field
        nfield = (n < pFields.size())? n: 0;

        m_f->m_exp[n] = m_f->AppendExpList(
                            m_f->m_numHomogeneousDir, m_variables[0]);
        m_f->m_exp[n]->SetWaveSpace(false);

        ASSERTL1(pFields[nfield]->GetNcoeffs() == m_outFields[n].size(),
                 "pFields[nfield] does not have the "
                 "same number of coefficients as m_outFields[n]");

        m_f->m_exp[n]->ExtractCoeffsToCoeffs(pFields[nfield], m_outFields[n],
                                             m_f->m_exp[n]->UpdateCoeffs());

        m_f->m_exp[n]->BwdTrans( m_f->m_exp[n]->GetCoeffs(),
                                 m_f->m_exp[n]->UpdatePhys());
    }
    m_f->m_variables= m_variables;
}

// This function checks validity conditions for the list of modules provided
void FilterFieldConvert::CheckModules(std::vector<ModuleSharedPtr> &modules)
{
    // Count number of modules by priority
    Array< OneD, int>  modulesCount(SIZE_ModulePriority,0);
    for (int i = 0; i < modules.size(); ++i)
    {
        ++modulesCount[modules[i]->GetModulePriority()];
    }

    // FilterFieldConvert already starts with m_exp, so anything before
    //    eModifyExp is not valid, and also eCreatePts
    if( modulesCount[eCreateGraph] != 0     ||
        modulesCount[eCreateFieldData] != 0 ||
        modulesCount[eModifyFieldData] != 0 ||
        modulesCount[eCreateExp] != 0       ||
        modulesCount[eFillExp] != 0         ||
        modulesCount[eCreatePts] != 0)
    {
        std::stringstream ss;
        ss << "Module(s): ";
        for (int i = 0; i < modules.size(); ++i)
        {
            if(modules[i]->GetModulePriority() == eCreateGraph     ||
               modules[i]->GetModulePriority() == eCreateFieldData ||
               modules[i]->GetModulePriority() == eModifyFieldData ||
               modules[i]->GetModulePriority() == eCreateExp       ||
               modules[i]->GetModulePriority() == eFillExp         ||
               modules[i]->GetModulePriority() == eCreatePts)
            {
                ss << modules[i]->GetModuleName()<<" ";
            }
        }
        ss << "not compatible with FilterFieldConvert.";
        ASSERTL0(false, ss.str());
    }

    // Modules of type eConvertExpToPts are not compatible with eBndExtraction
    if( modulesCount[eConvertExpToPts] != 0 &&
        modulesCount[eBndExtraction]   != 0)
    {
        std::stringstream ss;
        ss << "Module(s): ";
        for (int i = 0; i < modules.size(); ++i)
        {
            if(modules[i]->GetModulePriority() == eBndExtraction)
            {
                ss << modules[i]->GetModuleName()<<" ";
            }
        }
        ss << "is not compatible with module(s): ";
        for (int i = 0; i < modules.size(); ++i)
        {
            if(modules[i]->GetModulePriority() == eConvertExpToPts)
            {
                ss << modules[i]->GetModuleName()<<" ";
            }
        }
        ss << ".";
        ASSERTL0(false, ss.str());
    }

}

}
}
