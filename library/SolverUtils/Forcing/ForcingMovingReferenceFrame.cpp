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
// Description: Solving the absolute flow in a moving body frame,
// by adding (U0 + Omega X (x - x0)) . grad u - Omega X u
// as the body force.
// U0 is the translational velocity of the body frame.
// Omega is the angular velocity.
// x0 is the rotation pivot in the body frame.
// All vectors use the basis of the body frame.
// Translational motion is allowed for all dimensions.
// Rotation is not allowed for 1D, 2DH1D, 3DH2D.
// Rotation in z direction is allowed for 2D and 3DH1D.
// Rotation in 3 directions are allowed for 3D.
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
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isH1d, false);
    m_session->MatchSolverInfo("Homogeneous", "2D", m_isH2d, false);
    bool singleMode, halfMode;
    m_session->MatchSolverInfo("ModeType", "SingleMode", singleMode, false);
    m_session->MatchSolverInfo("ModeType", "HalfMode",   halfMode,   false);
    if (singleMode || halfMode)
    {
        m_isH1d = false;
    }

    int  npoints = pFields[0]->GetNpoints();
    int  expdim = m_isH2d ? 1 : pFields[0]->GetGraph()
                                ->GetMeshDimension();
    m_spacedim    = expdim + (m_isH1d ? 1 : 0) + (m_isH2d ? 2 : 0);
    m_NumVariable = m_spacedim;

    m_hasPlane0 = true;
    if (m_isH1d)
    {
       m_hasPlane0 = pFields[0]->GetZIDs()[0] == 0;
    }

    const TiXmlElement *funcNameElmt = pForce->FirstChildElement(
            "FRAMEVELOCITY");
    ASSERTL0(funcNameElmt, "Requires FRAMEVELOCITY tag specifying function "
                           "name which prescribes velocity of the moving "
                           "frame.");

    m_funcName = funcNameElmt->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName),
             "Function '" + m_funcName + "' not defined.");

    for (int i = 0; i < 6; ++i)
    {
        m_frameVelocity[i] = 0.;
    }

    m_hasRotation = false;
    //frame linear velocity
    for (int i = 0; i < m_spacedim; ++i)
    {
        std::string s_FieldStr = m_session->GetVariable(i);

        if (m_session->DefinesFunction(m_funcName, s_FieldStr))
        {
            m_frameFunction[i] = m_session->
                GetFunction(m_funcName, s_FieldStr);
        }
    }

    if (expdim==1)
    {
        return;
    }
    //frame angular velocity
    std::vector<std::string> angularVar = {"Omega_x", "Omega_y", "Omega_z"};
    for (int i = (expdim==2 ? 2 : 0); i < 3; ++i)
    {
        std::string s_FieldStr = angularVar[i];

        if (m_session->DefinesFunction(m_funcName, s_FieldStr))
        {
            m_hasRotation = true;
            m_frameFunction[m_spacedim + i] = m_session->
                GetFunction(m_funcName, s_FieldStr);
        }
    }
    if (m_hasRotation)
    {
        m_coords = Array<OneD, Array<OneD, NekDouble>>(3);
        for(int j=0; j<m_spacedim; ++j)
        {
            m_coords[j] = Array<OneD, NekDouble>(npoints);
        }
        for(int j=m_spacedim; j<3; ++j)
        {
            m_coords[j] = NullNekDouble1DArray;
        }
        pFields[0]->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
        std::vector<std::string> pivotVar = {"x0", "y0", "z0"};
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string s_FieldStr = pivotVar[i];

            if (m_session->DefinesFunction(m_funcName, s_FieldStr))
            {
                NekDouble x0 = m_session->GetFunction(m_funcName,
                    s_FieldStr)->Evaluate();
                Vmath::Sadd(npoints, -x0, m_coords[i], 1, m_coords[i], 1);
            }
        }
    }
}

/**
 * @brief Updates the forcing array with the current required forcing.
 * @param pFields
 * @param time
 */
void ForcingMovingReferenceFrame::Update(const NekDouble &time)
{
    for (auto it : m_frameFunction)
    {
        m_frameVelocity[it.first] = it.second->Evaluate(0., 0., 0., time);
    }
}

/**
 * @brief Adds the body force, -Omega X u.
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
    // If there is no rotation, body force is zero,
    // nothing needs to be done here.
    if (!m_hasRotation)
    {
        return;
    }
    Update(time);
    int  npoints = fields[0]->GetNpoints();
    addRotation(npoints, outarray, -1., inarray, outarray);
}

/**
 * @brief outarray = inarray0 + angVelScale Omega x inarray1
 */
void ForcingMovingReferenceFrame::addRotation(
            int npoints,
            const Array<OneD, Array<OneD, NekDouble> > &inarray0,
            NekDouble angVelScale,
            const Array<OneD, Array<OneD, NekDouble> > &inarray1,
            Array<OneD, Array<OneD, NekDouble> > &outarray)
{
    ASSERTL0(&inarray1!=&outarray, "inarray1 and outarray "
                        "should have different addresses.");

    if ((m_spacedim>=2 && &inarray0 != &outarray) ||
        (m_spacedim>=2 && m_frameFunction.count(m_spacedim+2)) )
    {
        Vmath::Svtvp(npoints,
                    -m_frameVelocity[m_spacedim+2]*angVelScale,
                    inarray1[1], 1,
                    inarray0[0], 1,
                    outarray[0], 1);
        Vmath::Svtvp(npoints,
                    +m_frameVelocity[m_spacedim+2]*angVelScale,
                    inarray1[0], 1,
                    inarray0[1], 1,
                    outarray[1], 1);
    }

    if ((m_spacedim==3 && &inarray0 != &outarray) ||
        (m_spacedim==3 && m_frameFunction.count(m_spacedim+0)) )
    {
        Vmath::Svtvp(npoints,
                    +m_frameVelocity[m_spacedim+0]*angVelScale,
                    inarray1[1], 1,
                    inarray0[2], 1,
                    outarray[2], 1);
    }

    if (m_spacedim==3 && m_frameFunction.count(m_spacedim+0))
    {
        Vmath::Svtvp(npoints,
                    -m_frameVelocity[m_spacedim+0]*angVelScale,
                    inarray1[2], 1,
                    outarray[1], 1,
                    outarray[1], 1);
    }

    if (m_spacedim==3 && m_frameFunction.count(m_spacedim+1))
    {
        Vmath::Svtvp(npoints,
                    +m_frameVelocity[m_spacedim+1]*angVelScale,
                    inarray1[2], 1,
                    outarray[0], 1,
                    outarray[0], 1);
        Vmath::Svtvp(npoints,
                    -m_frameVelocity[m_spacedim+1]*angVelScale,
                    inarray1[0], 1,
                    outarray[2], 1,
                    outarray[2], 1);
    }
}

void ForcingMovingReferenceFrame::v_PreApply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
{
    Update(time);
    int npoints = fields[0]->GetNpoints();
    if (m_isH2d && fields[0]->GetWaveSpace())
    {
        for (int i=0; i<m_spacedim; ++i)
        {
            if (m_frameFunction.count(i))
            {
                Array<OneD, NekDouble> tmpphys(npoints, -m_frameVelocity[i]);
                Array<OneD, NekDouble> tmpcoef(npoints);
                fields[0]->HomogeneousFwdTrans(tmpphys, tmpcoef);
                Vmath::Vadd(npoints, tmpcoef, 1, inarray[i], 1, outarray[i], 1);
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
    }
    else
    {
        int npoints0 = npoints;
        if (m_isH1d && fields[0]->GetWaveSpace())
        {
            npoints0 = m_hasPlane0 ? fields[0]->GetPlane(0)->GetNpoints() : 0;
        }
        for (int i=0; i<m_spacedim; ++i)
        {
            if (m_frameFunction.count(i))
            {
                Vmath::Sadd(npoints0, -m_frameVelocity[i], inarray[i], 1, outarray[i], 1);
                if (&inarray != &outarray && npoints != npoints0)
                {
                    Array<OneD, NekDouble> tmp = outarray[i]+npoints0;
                    Vmath::Vcopy(npoints - npoints0, inarray[i] + npoints0, 1, tmp, 1);
                }
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
        if (m_hasRotation)
        {
            addRotation(npoints0, outarray, -1., m_coords, outarray);
        }
    }
}

}
}