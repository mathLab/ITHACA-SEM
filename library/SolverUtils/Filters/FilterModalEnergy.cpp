///////////////////////////////////////////////////////////////////////////////
//
// File FilterModalEnergy.cpp
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
// Description: Output values of the modal energy
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SolverUtils/Filters/FilterModalEnergy.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterModalEnergy::className = GetFilterFactory().
    RegisterCreatorFunction("ModalEnergy", FilterModalEnergy::create);

/**
 *  Constructor.
 */
FilterModalEnergy::FilterModalEnergy(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams) :
    Filter(pSession, pEquation)
{
    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    if (!(m_outputFile.length() >= 4
          && m_outputFile.substr(m_outputFile.length() - 4) == ".mdl"))
    {
        m_outputFile += ".mdl";
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }


    m_session->MatchSolverInfo("Homogeneous", "1D", m_isHomogeneous1D, false);
    m_session->MatchSolverInfo("Homogeneous", "2D", m_isHomogeneous2D, false);
    m_session->MatchSolverInfo("CalculatePerturbationEnergy", "True",
                               m_PertEnergy, false);
    m_session->LoadParameter  ("NumQuadPointsError", m_NumQuadPointsError, 0);
    m_EqTypeStr = m_session->GetSolverInfo("EQTYPE");

    // OutputPlane
    if(m_isHomogeneous1D || m_isHomogeneous2D)
    {
        m_session->LoadParameter("LZ", m_LhomZ);

        it = pParams.find("OutputPlane");
        if (it == pParams.end())
        {
            m_outputPlane = 0;
        }
        else
        {
            LibUtilities::Equation equ(
                m_session->GetInterpreter(), it->second);
            m_outputPlane = round(equ.Evaluate());
        }
    }

    m_fld = LibUtilities::FieldIO::CreateDefault(pSession);
}

/**
 *  Destructor.
 */
FilterModalEnergy::~FilterModalEnergy()
{

}

/**
 *  Initialize the parallel communication and the output stream.
 */
void FilterModalEnergy::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    if (vComm->GetRank() == 0)
    {
        // Open output stream
        bool adaptive;
        m_session->MatchSolverInfo("Driver", "Adaptive",
                                    adaptive, false);
        if (adaptive)
        {
            m_outputStream.open(m_outputFile.c_str(), ofstream::app);
        }
        else
        {
            m_outputStream.open(m_outputFile.c_str());
        }
        if(m_isHomogeneous1D)
        {
            m_outputStream << "# Time,  Fourier Mode, Energy ";
            m_outputStream << endl;
        }
        else
        {
            m_outputStream << "# Time, Energy ";
            m_outputStream << endl;
        }

    }

    m_index = 0;
    v_Update(pFields, time);
}


/**
 *  Update the modal energy every m_outputFrequency.
 */
void FilterModalEnergy::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }

    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    // Homogeneous 1D implementation
    if (m_isHomogeneous1D)
    {
        int colrank = vComm->GetColumnComm()->GetRank();
        int nproc   = vComm->GetColumnComm()->GetSize();
        m_npointsZ  = (m_session->GetParameter("HomModesZ"));
        int locsize = m_npointsZ/nproc/2;

        Array<OneD, NekDouble> energy    (locsize, 0.0);
        Array<OneD, NekDouble> energy_tmp(locsize, 0.0);
        Array<OneD, NekDouble> tmp;

        // Calculate the energy of the perturbation for stability
        // analysis
        if (m_PertEnergy)
        {
            // Compressible Flow Solver
            if (m_EqTypeStr=="EulerCFE"   ||
                m_EqTypeStr=="EulerADCFE" ||
                m_EqTypeStr=="NavierStokesCFE")
            {
                ASSERTL0(false, "Stability analysis module not "
                         "implemented for the Compressible Flow "
                         "Solver. Please remove the function BaseFlow "
                         "from your .xml file");
            }
            // Incompressible Navier-Stokes Solver
            else
            {
                SpatialDomains::MeshGraphSharedPtr graphShrPtr =
                    SpatialDomains::MeshGraph::Read(m_session);
                SetUpBaseFields(graphShrPtr);
                string file = m_session->
                    GetFunctionFilename("BaseFlow", 0);
                ImportFldBase(file);

                for (int i = 0; i < pFields.size()-1; ++i)
                {
                    Vmath::Vsub(pFields[i]->GetNcoeffs(),
                                pFields[i]->GetCoeffs(), 1,
                                m_base [i]->GetCoeffs(), 1,
                                pFields[i]->UpdateCoeffs(), 1);

                    energy_tmp = pFields[i]->HomogeneousEnergy();
                    Vmath::Vadd(locsize, energy_tmp, 1,
                                energy, 1, energy, 1);

                    Vmath::Vadd(pFields[i]->GetNcoeffs(),
                                pFields[i]->GetCoeffs(), 1,
                                m_base[i]->GetCoeffs(), 1,
                                pFields[i]->UpdateCoeffs(), 1);
                }
            }
        }
        // Calculate the modal energy for general simulation
        else
        {
            // Compressible Flow Solver
            if (m_EqTypeStr=="EulerCFE"   ||
                m_EqTypeStr=="EulerADCFE" ||
                m_EqTypeStr=="NavierStokesCFE")
            {
                // Extracting kinetic energy
                for (int i = 1; i < pFields.size()-1; ++i)
                {
                    energy_tmp = pFields[i]->HomogeneousEnergy();
                    Vmath::Vadd(locsize, energy_tmp, 1,
                                energy, 1, energy, 1);
                }
            }
            // Incompressible Navier-Stokes Solver
            else
            {
                // Extracting kinetic energy
                for (int i = 0; i < pFields.size()-1; ++i)
                {
                    energy_tmp = pFields[i]->HomogeneousEnergy();
                    Vmath::Vadd(locsize, energy_tmp, 1,
                                energy, 1, energy, 1);
                }
            }
        }

        // Send to root process
        if (colrank == 0)
        {
            int j, m = 0;

            for (j = 0; j < energy.size(); ++j, ++m)
            {
                m_outputStream << setw(10) << time
                               << setw(5)  << m
                               << setw(18) << energy[j] << endl;
            }

            for (int i = 1; i < nproc; ++i)
            {
                vComm->GetColumnComm()->Recv(i, energy);

                for (j = 0; j < energy.size(); ++j, ++m)
                {
                    m_outputStream << setw(10) << time
                                   << setw(5)  << m
                                   << setw(18) << energy[j] << endl;
                }
            }
        }
        else
        {
            vComm->GetColumnComm()->Send(0, energy);
        }
    }
    // Homogeneous 2D implementation
    else if (m_isHomogeneous2D)
    {
        ASSERTL0(false, "3D Homogeneous 2D energy "
                 "dumping not implemented yet");
    }
    // General implementation
    else
    {
        // Compressible Flow Solver
        if (m_EqTypeStr=="EulerCFE"   ||
            m_EqTypeStr=="EulerADCFE" ||
            m_EqTypeStr=="NavierStokesCFE")
        {
            // Total energy
            NekDouble energy = 0.0;
            for (int i = 1; i < pFields.size()-1; ++i)
            {
                pFields[i]->SetPhysState(true);
                NekDouble norm = L2Error(pFields, i, time);
                energy += norm * norm;
            }

            m_outputStream << setprecision(6) << time;
            m_outputStream.width(25);
            m_outputStream << setprecision(8) << 0.5*energy;
            m_outputStream << endl;
        }
        // Incompressible Navier-Stokes Solver
        else
        {
            // Kinetic energy
            NekDouble energy = 0.0;
            for (int i = 0; i < pFields.size()-1; ++i)
            {
                pFields[i]->SetPhysState(true);
                NekDouble norm = L2Error(pFields, i, time);
                energy += norm * norm;
            }
            m_outputStream << setprecision(6) << time;
            m_outputStream.width(25);
            m_outputStream << setprecision(8) << 0.5*energy;
            m_outputStream << endl;
        }
    }
}

/**
 *  Close the output stream.
 */
void FilterModalEnergy::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(time);

    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}

/**
 *  Calculate the L2 norm of a given field for calculating the
 *  modal energy.
 */
NekDouble FilterModalEnergy::L2Error(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    unsigned int field,
    const NekDouble &time)
{
    boost::ignore_unused(time);

    NekDouble L2error = -1.0;
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    if (m_NumQuadPointsError == 0)
    {
        if (pFields[field]->GetPhysState() == false)
        {
            pFields[field]->BwdTrans(pFields[field]->GetCoeffs(),
                                     pFields[field]->UpdatePhys());
        }
    }

    L2error = pFields[field]->L2(pFields[field]->GetPhys());
    return L2error;
}

/**
 *  Setup the base fields in case of stability analyses.
 */
void FilterModalEnergy::SetUpBaseFields(
    SpatialDomains::MeshGraphSharedPtr &graphShrPtr)
{
    int i;
    int m_expdim  = graphShrPtr->GetMeshDimension();

    //definition of the projection tipe:
    if(m_session->DefinesSolverInfo("PROJECTION"))
    {
        std::string ProjectStr = m_session->GetSolverInfo("PROJECTION");

        if ((ProjectStr == "Continuous") ||
            (ProjectStr == "Galerkin")   ||
            (ProjectStr == "CONTINUOUS") ||
            (ProjectStr == "GALERKIN"))
        {
            m_projectionType = MultiRegions::eGalerkin;
        }
        else if ((ProjectStr == "MixedCGDG") ||
                 (ProjectStr == "Mixed_CG_Discontinuous"))
        {
            m_projectionType = MultiRegions::eMixed_CG_Discontinuous;
        }
        else if(ProjectStr == "DisContinuous")
        {
            m_projectionType = MultiRegions::eDiscontinuous;
        }
        else
        {
            ASSERTL0(false, "PROJECTION value not recognised");
        }
    }
    else
    {
        cerr << "Projection type not specified in SOLVERINFO,"
                        "defaulting to continuous Galerkin" << endl;
        m_projectionType = MultiRegions::eGalerkin;
    }

    if (m_session->DefinesSolverInfo("ModeType"))
    {
        m_session->MatchSolverInfo("ModeType", "SingleMode",
                                   m_SingleMode, false);
        m_session->MatchSolverInfo("ModeType", "HalfMode",
                                   m_HalfMode, false);
        m_session->MatchSolverInfo("ModeType", "MultipleModes",
                                   m_MultipleModes, false);
    }

    m_session->MatchSolverInfo("USEFFT","FFTW", m_useFFT, false);
    m_session->MatchSolverInfo("DEALIASING", "True",
                               m_homogen_dealiasing, false);

    // Stability Analysis flags
    if (m_session->DefinesSolverInfo("ModeType"))
    {
        if (m_SingleMode)
        {
            m_npointsZ = 2;
        }
        else if (m_HalfMode)
        {
            m_npointsZ = 1;
        }
        else if (m_MultipleModes)
        {
            m_npointsZ = m_session->GetParameter("HomModesZ");
        }
        else
        {
            ASSERTL0(false, "SolverInfo ModeType not valid");
        }
    }
    else
    {
        m_npointsZ = m_session->GetParameter("HomModesZ");
    }

    if (m_projectionType == MultiRegions::eGalerkin ||
        m_projectionType == MultiRegions::eMixed_CG_Discontinuous)
    {
        switch (m_expdim)
        {
            case 1:
            {
                for(i = 0; i < m_base.size(); i++)
                {
                    m_base[i] = MemoryManager<MultiRegions::ContField1D>
                        ::AllocateSharedPtr(m_session, graphShrPtr,
                                            m_session->GetVariable(0));
                }
            }
            break;
            case 2:
            {
                if (m_isHomogeneous1D)
                {
                    if (m_SingleMode)
                    {
                        const LibUtilities::PointsKey PkeyZ(
                            m_npointsZ,
                            LibUtilities::eFourierSingleModeSpaced);
                        const LibUtilities::BasisKey  BkeyZ(
                            LibUtilities::eFourier,
                            m_npointsZ, PkeyZ);

                        for (i = 0 ; i < m_base.size(); i++)
                        {
                            m_base[i] = MemoryManager<MultiRegions::
                                ContField3DHomogeneous1D>::
                                    AllocateSharedPtr(
                                        m_session, BkeyZ, m_LhomZ,
                                        m_useFFT, m_homogen_dealiasing,
                                        graphShrPtr,
                                        m_session->GetVariable(i));

                            m_base[i]->SetWaveSpace(true);
                        }
                    }
                    else if (m_HalfMode)
                    {
                        //1 plane field (half mode expansion)
                        const LibUtilities::PointsKey PkeyZ(
                            m_npointsZ,
                            LibUtilities::eFourierSingleModeSpaced);
                        const LibUtilities::BasisKey  BkeyZ(
                            LibUtilities::eFourierHalfModeRe,
                            m_npointsZ,PkeyZ);

                        for (i = 0 ; i < m_base.size(); i++)
                        {
                            m_base[i] = MemoryManager<MultiRegions::
                                ContField3DHomogeneous1D>::
                                    AllocateSharedPtr(
                                        m_session, BkeyZ, m_LhomZ,
                                        m_useFFT, m_homogen_dealiasing,
                                        graphShrPtr,
                                        m_session->GetVariable(i));

                            m_base[i]->SetWaveSpace(true);
                        }
                    }
                    else
                    {
                        const LibUtilities::PointsKey PkeyZ(
                            m_npointsZ,
                            LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyZ(
                            LibUtilities::eFourier, m_npointsZ, PkeyZ);

                        for (i = 0 ; i < m_base.size(); i++)
                        {
                            m_base[i] = MemoryManager<MultiRegions::
                                ContField3DHomogeneous1D>::
                                    AllocateSharedPtr(
                                        m_session, BkeyZ, m_LhomZ,
                                        m_useFFT, m_homogen_dealiasing,
                                        graphShrPtr,
                                        m_session->GetVariable(i));

                            m_base[i]->SetWaveSpace(false);
                        }
                    }
                }
                else
                {
                    i = 0;
                    MultiRegions::ContField2DSharedPtr firstbase =
                        MemoryManager<MultiRegions::ContField2D>::
                            AllocateSharedPtr(
                                m_session,graphShrPtr,
                                m_session->GetVariable(i));

                    m_base[0] = firstbase;

                    for (i = 1 ; i < m_base.size(); i++)
                    {
                        m_base[i] = MemoryManager<MultiRegions::
                            ContField2D>::AllocateSharedPtr(
                                *firstbase, graphShrPtr,
                                m_session->GetVariable(i));
                    }
                }
            }
            break;
            case 3:
            {
                MultiRegions::ContField3DSharedPtr firstbase =
                    MemoryManager<MultiRegions::ContField3D>::
                        AllocateSharedPtr(m_session, graphShrPtr,
                                          m_session->GetVariable(0));
                m_base[0] = firstbase;
                for (i = 1 ; i < m_base.size(); i++)
                {
                    m_base[i] = MemoryManager<MultiRegions::
                    ContField3D>::AllocateSharedPtr(
                        *firstbase, graphShrPtr,
                        m_session->GetVariable(0));
                }
            }
            break;
            default:
                NEKERROR(ErrorUtil::efatal,
                         "Expansion dimension not recognised");
                break;
        }
    }
    else
    {
        switch (m_expdim)
        {
            case 1:
            {
                // need to use zero for variable as may be more base
                // flows than variables
                for (i = 0 ; i < m_base.size(); i++)
                {
                    m_base[i] = MemoryManager<MultiRegions::
                        DisContField1D>::AllocateSharedPtr(
                            m_session, graphShrPtr,
                            m_session->GetVariable(0));
                }
                break;
            }
            case 2:
            {
                for (i = 0 ; i < m_base.size(); i++)
                {
                    m_base[i] = MemoryManager<MultiRegions::
                        DisContField2D>::AllocateSharedPtr(
                            m_session, graphShrPtr,
                            m_session->GetVariable(0));
                }
                break;
            }
            case 3:
                NEKERROR(ErrorUtil::efatal, "3D not set up");
                break;
            default:
                NEKERROR(ErrorUtil::efatal,
                         "Expansion dimension not recognised");
                break;
        }
    }
}

/**
 *  Import the base flow fld file.
 */
void FilterModalEnergy::ImportFldBase(
    std::string pInfile)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::vector<NekDouble> > FieldData;

    // Get Homogeneous
    m_fld->Import(pInfile,FieldDef,FieldData);

    int nvar = m_session->GetVariables().size();
    if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
    }
    // Copy FieldData into m_fields
    for (int j = 0; j < nvar; ++j)
    {
        for (int i = 0; i < FieldDef.size(); ++i)
        {
            bool flag =
                FieldDef[i]->m_fields[j] == m_session->GetVariable(j);

            ASSERTL0(flag, (std::string("Order of ") + pInfile
                        + std::string(" data and that defined in "
                            "m_boundaryconditions differs")).c_str());

            m_base[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                           FieldDef[i]->m_fields[j],
                                           m_base[j]->UpdateCoeffs());
        }
    }
}

/**
 *  Flag for time-dependent flows.
 */
bool FilterModalEnergy::v_IsTimeDependent()
{
    return true;
}
}
}
