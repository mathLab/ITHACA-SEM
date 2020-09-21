///////////////////////////////////////////////////////////////////////////////
//
// File DriverAdaptive.cpp
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
// Description: Driver class for adaptive solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <SolverUtils/DriverAdaptive.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdHexExp.h>
#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

string DriverAdaptive::className = GetDriverFactory().RegisterCreatorFunction(
    "Adaptive", DriverAdaptive::create);
string DriverAdaptive::driverLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("Driver", "Adaptive", 0);

/**
 *
 */
DriverAdaptive::DriverAdaptive(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : Driver(pSession, pGraph)
{
}

/**
 *
 */
DriverAdaptive::~DriverAdaptive()
{
}

/**
 *
 */
void DriverAdaptive::v_InitObject(ostream &out)
{
    Driver::v_InitObject(out);
}

void DriverAdaptive::v_Execute(ostream &out)
{
    time_t starttime, endtime;
    NekDouble CPUtime;

    m_equ[0]->PrintSummary(out);

    // First run using original order
    time(&starttime);
    m_equ[0]->DoInitialise();

    // Obtain initial time in case a restart was used
    NekDouble startTime = m_equ[0]->GetFinalTime();
    m_equ[0]->DoSolve();

    // Load session parameters and solver info
    bool isHomogeneous1D;
    int nRuns, minP, maxP, sensorVar;
    NekDouble lowerTol, upperTol;
    m_session->LoadParameter  ("NumRuns",                nRuns,     1);
    m_session->LoadParameter  ("AdaptiveMinModes",       minP,      4);
    m_session->LoadParameter  ("AdaptiveMaxModes",       maxP,      12);
    m_session->LoadParameter  ("AdaptiveLowerTolerance", lowerTol,  1e-8);
    m_session->LoadParameter  ("AdaptiveUpperTolerance", upperTol,  1e-6);
    m_session->LoadParameter  ("AdaptiveSensorVariable", sensorVar, 0);
    m_session->MatchSolverInfo("Homogeneous", "1D", isHomogeneous1D, false);

    // Get number of elements and planes
    int nExp, nPlanes;
    if (isHomogeneous1D)
    {
        nExp    = m_equ[0]->UpdateFields()[0]->GetPlane(0)->GetExpSize();
        nPlanes = m_equ[0]->UpdateFields()[0]->GetZIDs().size();
    }
    else
    {
        nExp    = m_equ[0]->UpdateFields()[0]->GetExpSize();
        nPlanes = 1;
    }
    int  expdim   = m_equ[0]->UpdateFields()[0]->GetGraph()->GetMeshDimension();

    int       nFields  = m_equ[0]->UpdateFields().size();
    int       numSteps = m_session->GetParameter("NumSteps");
    NekDouble period   = m_session->GetParameter("TimeStep") * numSteps;

    Array<OneD, NekDouble> coeffs, phys, physReduced, tmpArray;

    // Get mapping
    GlobalMapping::MappingSharedPtr mapping;
    mapping = GlobalMapping::Mapping::Load(m_session,
                                           m_equ[0]->UpdateFields());

    // Adaptive loop
    Array<OneD, int> P(expdim);
    Array<OneD, int> numPoints(expdim);
    Array<OneD, LibUtilities::PointsKey> ptsKey(expdim);
    LocalRegions::ExpansionSharedPtr Exp;
    for (int i = 1; i < nRuns; i++)
    {
        // Get field expansions
        Array<OneD, MultiRegions::ExpListSharedPtr> fields =
            m_equ[0]->UpdateFields();

        // Determine the change to be applied in the order
        map<int, int> deltaP;
        int offset = 0;
        for (int n = 0; n < nExp; n++)
        {
            offset = fields[sensorVar]->GetPhys_Offset(n);
            Exp    = fields[sensorVar]->GetExp(n);

            for( int k = 0; k < expdim; ++k)
            {
                P[k]         = Exp->GetBasis(k)->GetNumModes();
                numPoints[k] = Exp->GetBasis(k)->GetNumPoints();
                ptsKey[k]    = LibUtilities::PointsKey (numPoints[k],
                               Exp->GetBasis(k)->GetPointsType());
            }

            // Declare orthogonal basis.
            StdRegions::StdExpansionSharedPtr OrthoExp;
            switch (Exp->GetGeom()->GetShapeType())
            {
                case LibUtilities::eQuadrilateral:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - 1,
                                              ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - 1,
                                              ptsKey[1]);
                    OrthoExp = MemoryManager<
                        StdRegions::StdQuadExp>::AllocateSharedPtr(Ba, Bb);
                    break;
                }
                case LibUtilities::eTriangle:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - 1,
                                              ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B, P[1] - 1,
                                              ptsKey[1]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(
                            Ba, Bb);
                    break;
                }
                case LibUtilities::eTetrahedron:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - 1,
                                              ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B, P[1] - 1,
                                              ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_C, P[2] - 1,
                                              ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                case LibUtilities::ePyramid:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - 1,
                                              ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - 1,
                                              ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_C, P[2] - 1,
                                              ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdPyrExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                case LibUtilities::ePrism:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - 1,
                                              ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - 1,
                                              ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_B, P[2] - 1,
                                              ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                case LibUtilities::eHexahedron:
                {
                    LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, P[0] - 1,
                                              ptsKey[0]);
                    LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, P[1] - 1,
                                              ptsKey[1]);
                    LibUtilities::BasisKey Bc(LibUtilities::eOrtho_A, P[2] - 1,
                                              ptsKey[2]);
                    OrthoExp =
                        MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(
                            Ba, Bb, Bc);
                    break;
                }
                default:
                    ASSERTL0(false, "Shape not supported.");
                    break;
            }

            int nq = OrthoExp->GetTotPoints();

            NekDouble error = 0;
            NekDouble fac   = 0;
            NekDouble tmp   = 0;

            coeffs      = Array<OneD, NekDouble>(OrthoExp->GetNcoeffs());
            physReduced = Array<OneD, NekDouble>(OrthoExp->GetTotPoints());
            tmpArray    = Array<OneD, NekDouble>(OrthoExp->GetTotPoints(), 0.0);

            // Refinement based only on one variable
            for (int plane = 0; plane < nPlanes; plane++)
            {
                if (isHomogeneous1D)
                {
                    phys =
                        fields[sensorVar]->GetPlane(plane)->GetPhys() + offset;
                }
                else
                {
                    phys = fields[sensorVar]->GetPhys() + offset;
                }

                // Project solution to lower order
                OrthoExp->FwdTrans(phys, coeffs);
                OrthoExp->BwdTrans(coeffs, physReduced);

                // Calculate error =||phys-physReduced||^2 / ||phys||^2
                Vmath::Vsub(nq, phys,     1, physReduced, 1, tmpArray, 1);
                Vmath::Vmul(nq, tmpArray, 1, tmpArray,    1, tmpArray, 1);
                tmp = Exp->Integral(tmpArray);

                Vmath::Vmul(nq, phys, 1, phys, 1, tmpArray, 1);
                fac = Exp->Integral(tmpArray);

                tmp = abs(tmp / fac);

                if (tmp != tmp)
                {
                    ASSERTL0(false,
                             "Adaptive procedure encountered NaN value.");
                }

                // Get maximum value along planes
                error = (tmp > error) ? tmp : error;
            }

            // Combine planes from different processes
            m_session->GetComm()->GetColumnComm()->AllReduce(
                error, LibUtilities::ReduceMax);

            // Override tolerances if function is defined
            if (m_session->DefinesFunction("AdaptiveLowerTolerance"))
            {
                int nq = Exp->GetTotPoints();

                // Obtain points from the the element
                Array<OneD, NekDouble> xc0, xc1, xc2;
                xc0 = Array<OneD, NekDouble>(nq, 0.0);
                xc1 = Array<OneD, NekDouble>(nq, 0.0);
                xc2 = Array<OneD, NekDouble>(nq, 0.0);
                Exp->GetCoords(xc0, xc1, xc2);

                // Evaluate function from session file
                Array<OneD, NekDouble> tolerance(nq, 0.0);
                LibUtilities::EquationSharedPtr ffunc =
                    m_session->GetFunction("AdaptiveLowerTolerance", 0);
                ffunc->Evaluate(xc0, xc1, xc2, tolerance);
                lowerTol = Vmath::Vsum(nq, tolerance, 1) / nq;
            }

            if (m_session->DefinesFunction("AdaptiveUpperTolerance"))
            {
                int nq = Exp->GetTotPoints();

                // Obtain points from the the element
                Array<OneD, NekDouble> xc0, xc1, xc2;
                xc0 = Array<OneD, NekDouble>(nq, 0.0);
                xc1 = Array<OneD, NekDouble>(nq, 0.0);
                xc2 = Array<OneD, NekDouble>(nq, 0.0);
                Exp->GetCoords(xc0, xc1, xc2);

                // Evaluate function from session file
                Array<OneD, NekDouble> tolerance(nq, 0.0);
                LibUtilities::EquationSharedPtr ffunc =
                    m_session->GetFunction("AdaptiveUpperTolerance", 0);
                ffunc->Evaluate(xc0, xc1, xc2, tolerance);
                upperTol = Vmath::Vsum(nq, tolerance, 1) / nq;
            }

            // Determine what to do with the polynomial order
            if ((error > upperTol) && (P[0] < maxP))
            {
                deltaP[Exp->GetGeom()->GetGlobalID()] = 1;
            }
            else if ((error < lowerTol) && P[0] > minP)
            {
                deltaP[Exp->GetGeom()->GetGlobalID()] = -1;
            }
            else
            {
                deltaP[Exp->GetGeom()->GetGlobalID()] = 0;
            }
        }

        // Write new expansion section to the session reader and re-read graph.
        ReplaceExpansion(fields, deltaP);
        m_graph->ReadExpansions();

        // Reset GlobalLinSys Manager to avoid using too much memory
        //
        // @todo This could be made better by replacing individual matrices
        //       within the linear system.
        if (LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
                                     MultiRegions::GlobalLinSys>::
                PoolCreated(std::string("GlobalLinSys")))
        {
            LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,
                                     MultiRegions::GlobalLinSys>::
                ClearManager(std::string("GlobalLinSys"));
        }

        int chkNumber = m_equ[0]->GetCheckpointNumber();
        int chkSteps  = m_equ[0]->GetCheckpointSteps();

        // Initialise driver again
        Driver::v_InitObject(out);

        // Update mapping (must be before m_equ[0]->DoInitialise();)
        mapping->ReplaceField(m_equ[0]->UpdateFields());

        // Set chkSteps to zero to avoid writing initial condition
        m_equ[0]->SetCheckpointSteps(0);

        // Initialise equation
        m_equ[0]->DoInitialise();
        m_equ[0]->SetInitialStep(i * numSteps);
        m_equ[0]->SetSteps(i * numSteps + numSteps);
        m_equ[0]->SetTime(startTime + i * period);
        m_equ[0]->SetBoundaryConditions(startTime + i * period);
        m_equ[0]->SetCheckpointNumber(chkNumber);
        m_equ[0]->SetCheckpointSteps(chkSteps);

        // Project solution to new expansion
        for (int n = 0; n < nFields; n++)
        {
            m_equ[0]->UpdateFields()[n]->ExtractCoeffsToCoeffs(
                fields[n], fields[n]->GetCoeffs(),
                m_equ[0]->UpdateFields()[n]->UpdateCoeffs());
            m_equ[0]->UpdateFields()[n]->BwdTrans_IterPerExp(
                m_equ[0]->UpdateFields()[n]->GetCoeffs(),
                m_equ[0]->UpdateFields()[n]->UpdatePhys());
        }

        // Solve equation
        m_equ[0]->DoSolve();
    }

    time(&endtime);

    m_equ[0]->Output();

    if (m_comm->GetRank() == 0)
    {
        CPUtime = difftime(endtime, starttime);
        cout << "-------------------------------------------" << endl;
        cout << "Total Computation Time = " << CPUtime << "s" << endl;
        cout << "-------------------------------------------" << endl;
    }

    // Evaluate and output computation time and solution accuracy.
    // The specific format of the error output is essential for the
    // regression tests to work.

    // Evaluate L2 Error
    for (int i = 0; i < m_equ[0]->GetNvariables(); ++i)
    {
        Array<OneD, NekDouble> exactsoln(m_equ[0]->GetTotPoints(), 0.0);

        // Evaluate "ExactSolution" function, or zero array
        m_equ[0]->EvaluateExactSolution(i, exactsoln, m_equ[0]->GetFinalTime());

        NekDouble vL2Error   = m_equ[0]->L2Error  (i, exactsoln);
        NekDouble vLinfError = m_equ[0]->LinfError(i, exactsoln);

        if (m_comm->GetRank() == 0)
        {
            out << "L 2 error (variable " << m_equ[0]->GetVariable(i)
                << ") : " << vL2Error << endl;
            out << "L inf error (variable " << m_equ[0]->GetVariable(i)
                << ") : " << vLinfError << endl;
        }
    }
}

/**
 * @brief Update EXPANSIONS tag inside XML schema to reflect new polynomial
 * order distribution.
 *
 * @param fields  Input fields.
 * @param deltaP  Map of polynomial order expansions
 */
void DriverAdaptive::ReplaceExpansion(
    Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    map<int, int>                                deltaP)
{
    int nExp, nDim;
    int  expdim   = m_equ[0]->UpdateFields()[0]->GetGraph()->GetMeshDimension();

    // Get field definitions
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs =
        fields[0]->GetFieldDefinitions();

    if (fielddefs[0]->m_numHomogeneousDir == 1)
    {
        nDim = expdim+1;
    }
    else
    {
        nDim = expdim;
    }

    // Add variables to field definition
    for (int i = 0; i < fielddefs.size(); ++i)
    {
        for (int j = 0; j < fields.size(); ++j)
        {
            fielddefs[i]->m_fields.push_back(m_session->GetVariable(j));
        }
    }

    // Get tinyxml objects
    TiXmlElement *exp_tag = m_session->GetElement("NEKTAR/EXPANSIONS");

    // Clear current expansions
    exp_tag->Clear();

    // Write new expansion information
    for (int f = 0; f < fielddefs.size(); ++f)
    {
        nExp = fielddefs[f]->m_elementIDs.size();

        //---------------------------------------------
        // Write ELEMENTS
        TiXmlElement *elemTag = new TiXmlElement("ELEMENTS");
        exp_tag->LinkEndChild(elemTag);

        // Write FIELDS
        std::string fieldsString;
        {
            std::stringstream fieldsStringStream;
            bool first = true;
            for (std::vector<int>::size_type i = 0;
                 i < fielddefs[f]->m_fields.size(); i++)
            {
                if (!first)
                {
                    fieldsStringStream << ",";
                }
                fieldsStringStream << fielddefs[f]->m_fields[i];
                first = false;
            }
            fieldsString = fieldsStringStream.str();
        }
        elemTag->SetAttribute("FIELDS", fieldsString);

        // Write SHAPE
        std::string shapeString;
        {
            std::stringstream shapeStringStream;
            shapeStringStream
                << LibUtilities::ShapeTypeMap[fielddefs[f]->m_shapeType];
            if (fielddefs[f]->m_numHomogeneousDir == 1)
            {
                shapeStringStream << "-HomogenousExp1D";
            }
            else if (fielddefs[f]->m_numHomogeneousDir == 2)
            {
                shapeStringStream << "-HomogenousExp2D";
            }

            shapeString = shapeStringStream.str();
        }
        elemTag->SetAttribute("SHAPE", shapeString);

        // Write BASIS
        std::string basisString;
        {
            std::stringstream basisStringStream;
            bool first = true;
            for (std::vector<LibUtilities::BasisType>::size_type i = 0;
                 i < fielddefs[f]->m_basis.size(); i++)
            {
                if (!first)
                {
                    basisStringStream << ",";
                }
                basisStringStream
                    << LibUtilities::BasisTypeMap[fielddefs[f]->m_basis[i]];
                first = false;
            }
            basisString = basisStringStream.str();
        }
        elemTag->SetAttribute("BASIS", basisString);

        // Write homogeneuous length details
        if (fielddefs[f]->m_numHomogeneousDir)
        {
            std::string homoLenString;
            {
                std::stringstream homoLenStringStream;
                bool first = true;
                for (int i = 0; i < fielddefs[f]->m_numHomogeneousDir; ++i)
                {
                    if (!first)
                    {
                        homoLenStringStream << ",";
                    }
                    homoLenStringStream
                        << fielddefs[f]->m_homogeneousLengths[i];
                    first = false;
                }
                homoLenString = homoLenStringStream.str();
            }
            elemTag->SetAttribute("HOMOGENEOUSLENGTHS", homoLenString);
        }

        // Write homogeneuous planes/lines details
        if (fielddefs[f]->m_numHomogeneousDir)
        {
            if (fielddefs[f]->m_homogeneousYIDs.size() > 0)
            {
                std::string homoYIDsString;
                {
                    std::stringstream homoYIDsStringStream;
                    bool first = true;
                    for (int i = 0; i < fielddefs[f]->m_homogeneousYIDs.size();
                         i++)
                    {
                        if (!first)
                        {
                            homoYIDsStringStream << ",";
                        }
                        homoYIDsStringStream
                            << fielddefs[f]->m_homogeneousYIDs[i];
                        first = false;
                    }
                    homoYIDsString = homoYIDsStringStream.str();
                }
                elemTag->SetAttribute("HOMOGENEOUSYIDS", homoYIDsString);
            }

            if (fielddefs[f]->m_homogeneousZIDs.size() > 0)
            {
                std::string homoZIDsString;
                {
                    std::stringstream homoZIDsStringStream;
                    bool first = true;
                    for (int i = 0; i < fielddefs[f]->m_homogeneousZIDs.size();
                         i++)
                    {
                        if (!first)
                        {
                            homoZIDsStringStream << ",";
                        }
                        homoZIDsStringStream
                            << fielddefs[f]->m_homogeneousZIDs[i];
                        first = false;
                    }
                    homoZIDsString = homoZIDsStringStream.str();
                }
                elemTag->SetAttribute("HOMOGENEOUSZIDS", homoZIDsString);
            }
        }

        // Write NUMMODESPERDIR
        std::string numModesString;
        {
            std::stringstream numModesStringStream;

            numModesStringStream << "MIXORDER:";
            bool first = true;
            int eId;
            Array<OneD, int> order(nDim, 0);
            for (int n = 0; n < nExp; n++)
            {
                eId = fielddefs[f]->m_elementIDs[n];

                for (int i = 0; i < expdim; i++)
                {
                    order[i] = deltaP[eId];
                }

                for (int i = 0; i < nDim; i++)
                {
                    if (fielddefs[f]->m_uniOrder)
                    {
                        order[i] += fielddefs[f]->m_numModes[i];
                    }
                    else
                    {
                        order[i] += fielddefs[f]->m_numModes[n * nDim + i];
                    }

                    if (!first)
                    {
                        numModesStringStream << ",";
                    }

                    numModesStringStream << order[i];
                    first = false;
                }
            }
            numModesString = numModesStringStream.str();
        }
        elemTag->SetAttribute("NUMMODESPERDIR", numModesString);

        // Write IDs. Should ideally look at ways of compressing this stream if
        // just sequential.
        elemTag->SetAttribute(
            "ID", ParseUtils::GenerateSeqString(fielddefs[f]->m_elementIDs));
    }
}

}
}
