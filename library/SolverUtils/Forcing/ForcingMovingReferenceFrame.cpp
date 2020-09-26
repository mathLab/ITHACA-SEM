///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.cpp
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
// Description: Allows for a moving frame of reference, through adding c * du/dx
// to the body force, where c is the frame velocity vector
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Forcing/ForcingMovingReferenceFrame.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

std::string ForcingMovingReferenceFrame::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("MovingReferenceFrame",
                                ForcingMovingReferenceFrame::create,
                                "Moving Frame");
std::string ForcingMovingReferenceFrame::classNameField = GetForcingFactory().
        RegisterCreatorFunction("Field",
                                ForcingMovingReferenceFrame::create,
                                "Field Forcing");

/**
 * @brief
 * @param pSession
 * @param pEquation
 */
ForcingMovingReferenceFrame::ForcingMovingReferenceFrame(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation)
        : Forcing(pSession, pEquation)
{
}


/**
 * @brief Initialise the forcing module
 * @param pFields
 * @param pNumForcingFields
 * @param pForce
 */
void ForcingMovingReferenceFrame::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields,
        const TiXmlElement *pForce)
{
    boost::ignore_unused(pNumForcingFields);

    int  npoints = pFields[0]->GetNpoints();
    int  expdim  = pFields[0]->GetGraph()->GetMeshDimension();
    bool isH1d;
    bool isH2d;
    m_session->MatchSolverInfo("Homogeneous", "1D", isH1d, false);
    m_session->MatchSolverInfo("Homogeneous", "2D", isH2d, false);
    m_spacedim    = expdim + (isH1d ? 1 : 0) + (isH2d ? 2 : 0);
    m_NumVariable = m_spacedim;

    const TiXmlElement *funcNameElmt = pForce->FirstChildElement(
            "FRAMEVELOCITY");
    ASSERTL0(funcNameElmt, "Requires FRAMEVELOCITY tag specifying function "
                           "name which prescribes velocity of the moving "
                           "frame.");

    m_funcName = funcNameElmt->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName),
             "Function '" + m_funcName + "' not defined.");

    m_session->MatchSolverInfo("ModeType", "SingleMode", m_SingleMode, false);
    m_session->MatchSolverInfo("ModeType", "HalfMode",   m_HalfMode,   false);
    bool homogeneous = pFields[0]->GetExpType() == MultiRegions::e3DH1D ||
                       pFields[0]->GetExpType() == MultiRegions::e3DH2D;
    m_transform = (m_SingleMode || m_HalfMode || homogeneous);

    m_homogen_dealiasing = m_session->DefinesSolverInfo("dealiasing");

    m_Forcing       = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
    m_frameVelocity = Array<OneD, NekDouble>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_Forcing[i] = Array<OneD, NekDouble>(npoints, 0.0);

        std::string s_FieldStr = m_session->GetVariable(i);

        ASSERTL0(m_session->DefinesFunction(m_funcName, s_FieldStr),
                 "Variable '" + s_FieldStr + "' not defined.");
        auto ep = m_session->GetFunction(m_funcName, s_FieldStr);
        m_frameVelocity[i] = ep->Evaluate();
    }

    m_grad = Array<OneD, Array<OneD, NekDouble>>(m_spacedim * m_spacedim);
    for (int i = 0; i < m_spacedim * m_spacedim; ++i)
    {
        m_grad[i] = Array<OneD, NekDouble>(npoints);
    }

    Update(pFields, 0.0);
}


/**
 * @brief Calculates the gradient of the velocity fields.
 * @param pFields
 */
void ForcingMovingReferenceFrame::CalculateGradient(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    int  npoints = pFields[0]->GetNpoints();
    auto tmp     = Array<OneD, NekDouble>(npoints);

    // Evaluate Grad(u)
    switch (m_spacedim)
    {
        case 1:
            // du/dx
            pFields[0]->PhysDeriv(pFields[0]->GetPhys(), m_grad[0]);
            break;
        case 2:
            // du/dx du/dy
            pFields[0]->PhysDeriv(pFields[0]->GetPhys(), m_grad[0],
                                  m_grad[1]);
            // dv/dx dv/dy
            pFields[1]->PhysDeriv(pFields[1]->GetPhys(), m_grad[2],
                                  m_grad[3]);
            break;
        case 3:
            if (pFields[0]->GetWaveSpace() == true &&
                pFields[0]->GetExpType()   == MultiRegions::e3DH1D)
            {
                // take d/dx, d/dy  gradients in physical Fourier space
                pFields[0]->PhysDeriv(pFields[0]->GetPhys(),
                                      m_grad[0], m_grad[1]);
                pFields[1]->PhysDeriv(pFields[1]->GetPhys(),
                                      m_grad[3], m_grad[4]);

                pFields[0]->HomogeneousBwdTrans(pFields[2]->GetPhys(),
                                                tmp);
                pFields[0]->PhysDeriv(tmp, m_grad[6],
                                      m_grad[7]);

                // Take d/dz derivative using wave space field
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                      pFields[0]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[2]);
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                      pFields[1]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[5]);
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                      pFields[2]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[8]);
            }
            else if (pFields[0]->GetWaveSpace() == true &&
                     pFields[0]->GetExpType()   == MultiRegions::e3DH2D)
            {
                // take d/dx,  gradients in physical Fourier space
                pFields[0]->PhysDeriv(pFields[0]->GetPhys(),
                                      m_grad[0]);
                pFields[1]->PhysDeriv(pFields[1]->GetPhys(),
                                      m_grad[3]);

                pFields[0]->HomogeneousBwdTrans(pFields[2]->GetPhys(),
                                                tmp);
                pFields[0]->PhysDeriv(tmp, m_grad[6]);

                // Take d/dy derivative using wave space field
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],
                                      pFields[0]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[1]);
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],
                                      pFields[1]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[4]);
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],
                                      pFields[2]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[7]);

                // Take d/dz derivative using wave space field
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                      pFields[0]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[2]);
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                      pFields[1]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[5]);
                pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                      pFields[2]->GetPhys(), tmp);
                pFields[0]->HomogeneousBwdTrans(tmp, m_grad[8]);
            }
            else
            {
                pFields[0]->PhysDeriv(pFields[0]->GetPhys(),
                                      m_grad[0], m_grad[1],
                                      m_grad[2]);
                pFields[1]->PhysDeriv(pFields[1]->GetPhys(),
                                      m_grad[3], m_grad[4],
                                      m_grad[5]);
                pFields[2]->PhysDeriv(pFields[2]->GetPhys(),
                                      m_grad[6], m_grad[7],
                                      m_grad[8]);
            }
            break;
        default:
            ASSERTL0(false, "dimension unknown");
    }

}


/**
 * @brief Updates the forcing array with the current required forcing.
 * @param pFields
 * @param time
 */
void ForcingMovingReferenceFrame::Update(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    boost::ignore_unused(time);

    int npoints = pFields[0]->GetNpoints();
    Array<OneD, NekDouble> tmp(npoints, 0.0);

    CalculateGradient(pFields);

    switch (m_spacedim)
    {
        case 1:
            // - cu * du/dx
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[0],
                        m_grad[0], 1, m_Forcing[0], 1);
            break;
        case 2:
            // - cu * du/dx - cv * du/dy
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[0],
                        m_grad[0], 1, m_Forcing[0], 1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[1],
                        m_grad[1], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[0], 1, tmp, 1, m_Forcing[0],
                        1);
            // - cu * dv/dx - cv * dv/dy
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[0],
                        m_grad[2], 1, m_Forcing[1], 1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[1],
                        m_grad[3], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[1], 1, tmp, 1, m_Forcing[1],
                        1);
            break;
        case 3:
            //How do I deal with Homogeneous expansion?
            // - cu * du/dx - cv * du/dy - cz * du/dz
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[0],
                        m_grad[0], 1, m_Forcing[0], 1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[1],
                        m_grad[1], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[0], 1, tmp, 1, m_Forcing[0],
                        1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[2],
                        m_grad[2], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[0], 1, tmp, 1, m_Forcing[0],
                        1);
            // - cu * dv/dx - cv * dv/dy  - cz * dv/dz
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[0],
                        m_grad[3], 1, m_Forcing[1], 1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[1],
                        m_grad[4], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[1], 1, tmp, 1, m_Forcing[1],
                        1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[2],
                        m_grad[5], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[1], 1, tmp, 1, m_Forcing[1],
                        1);
            // - cu * dw/dx - cv * dw/dy  - cz * dw/dz
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[0],
                        m_grad[6], 1, m_Forcing[2], 1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[1],
                        m_grad[7], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[2], 1, tmp, 1, m_Forcing[2],
                        1);
            Vmath::Smul(npoints, -1.0 * m_frameVelocity[2],
                        m_grad[8], 1, tmp, 1);
            Vmath::Vadd(npoints, m_Forcing[2], 1, tmp, 1, m_Forcing[2],
                        1);
            break;
        default:
            ASSERTL0(false, "dimension unknown");
    }

    // If singleMode or halfMode, transform the forcing term to be in
    // physical space in the plane, but Fourier space in the homogeneous
    // direction
    if (m_transform)
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            pFields[0]->HomogeneousFwdTrans(m_Forcing[i], m_Forcing[i]);
        }
    }
}


/**
 * @brief Apply the forcing term
 * @param fields
 * @param inarray
 * @param outarray
 * @param time
 */
void ForcingMovingReferenceFrame::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> > &inarray,
              Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble &time)
{
    boost::ignore_unused(inarray);

    Update(fields, time);

    for (int i = 0; i < m_spacedim; i++)
    {
        Vmath::Vadd(outarray[i].size(), outarray[i], 1,
                    m_Forcing[i], 1, outarray[i], 1);
    }
}


}
}

