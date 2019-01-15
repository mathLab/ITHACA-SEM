///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingBody.cpp
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
// Description: Allows for a moving frame of reference, through adding c * du/dx to the body force, where c is the frame velocity vector
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingMovingFrame.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingMovingFrame::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("MovingFrame",
                                ForcingMovingFrame::create,
                                "Moving Frame");
    std::string ForcingMovingFrame::classNameField = GetForcingFactory().
        RegisterCreatorFunction("Field",
                                ForcingMovingFrame::create,
                                "Field Forcing");

    ForcingMovingFrame::ForcingMovingFrame(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation)
    {
    }

    void ForcingMovingFrame::v_InitObject(
                                   const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                                   const unsigned int& pNumForcingFields,
                                   const TiXmlElement* pForce)
    {
        m_NumVariable = pNumForcingFields;
        int npoints = pFields[0]->GetNpoints();

        const TiXmlElement* funcNameElmt = pForce->FirstChildElement("FRAMEVELOCITY");
        if(!funcNameElmt)
        {
            ASSERTL0(funcNameElmt, "Requires FRAMEVELOCITY tag "
                     "specifying function name which prescribes velocity of the moving frame.");
        }

        m_funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName), "Function '" + m_funcName + "' not defined.");

        bool singleMode, halfMode;
        m_session->MatchSolverInfo("ModeType","SingleMode",singleMode,false);
        m_session->MatchSolverInfo("ModeType","HalfMode",  halfMode,  false);
        bool homogeneous = pFields[0]->GetExpType() == MultiRegions::e3DH1D ||
                           pFields[0]->GetExpType() == MultiRegions::e3DH2D;
        m_transform = (singleMode || halfMode || homogeneous);

        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (npoints, 0.0);
        }

        m_frameVelocity = Array<OneD, NekDouble>(m_NumVariable); // Assume as many directions as variables
        for (int i = 0; i < m_NumVariable; ++i)
        {
            std::string  s_FieldStr   = m_session->GetVariable(i);
            auto ep = m_session->GetFunction(m_funcName, s_FieldStr);
            ASSERTL0(m_session->DefinesFunction(m_funcName, s_FieldStr),
                     "Variable '" + s_FieldStr + "' not defined.");
            m_frameVelocity[i] = ep->Evaluate();
        }

        //allocate space for gradient
        // Skip in case of empty partition
        if (pFields[0]->GetNumElmts() == 0)
        {
            return;
        }

        int expdim = pFields[0]->GetGraph()->GetMeshDimension();
        bool isH1d;
        bool isH2d;
        m_session->MatchSolverInfo("Homogeneous", "1D", isH1d, false);
        m_session->MatchSolverInfo("Homogeneous", "2D", isH2d, false);
        m_spacedim = expdim + (isH1d ? 1 : 0) + (isH2d ? 2 : 0);

        m_gradient = Array<OneD, Array<OneD, NekDouble>>(m_spacedim * m_spacedim);

        for (int i = 0; i < m_spacedim * m_spacedim; ++i)
        {
            m_gradient[i] = Array<OneD, NekDouble>(npoints);
        }

        Update(pFields, 0.0);
    }

    void ForcingMovingFrame::CalculateGradient(const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields)
    {
        // Evaluate Grad(u)
        switch(m_spacedim)
        {
            case 1:
                // du/dx
                pFields[0]->PhysDeriv(pFields[0]->GetPhys(), m_gradient[0]);
                break;
            case 2:
                // du/dx du/dy
                pFields[0]->PhysDeriv(pFields[0]->GetPhys(), m_gradient[0], m_gradient[1]);
                // dv/dx dv/dy
                pFields[1]->PhysDeriv(pFields[1]->GetPhys(), m_gradient[2], m_gradient[3]);
                break;
            case 3:
                //How do I deal with Homogeneous expansion?
                pFields[0]->PhysDeriv(pFields[0]->GetPhys(), m_gradient[0], m_gradient[1], m_gradient[2]);
                pFields[1]->PhysDeriv(pFields[1]->GetPhys(), m_gradient[3], m_gradient[4], m_gradient[5]);
                pFields[2]->PhysDeriv(pFields[2]->GetPhys(), m_gradient[6], m_gradient[7], m_gradient[8]);
                break;
            default:
                ASSERTL0(false,"dimension unknown");
        }
    }


    void ForcingMovingFrame::Update(
            const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const NekDouble &time)
    {
        int npoints = pFields[0]->GetNpoints();
        CalculateGradient( pFields );
        Array<OneD, NekDouble> tmp(npoints, 0.0);
        switch(m_spacedim)
        {
            case 1:
                // - cu * du/dx
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[0],m_gradient[0],1,m_Forcing[0], 1);
                break;
            case 2:
                // - cu * du/dx - cv * du/dy
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[0],m_gradient[0],1,m_Forcing[0], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[1],m_gradient[1],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[0], 1, tmp, 1, m_Forcing[0], 1);
                // - cu * dv/dx - cv * dv/dy
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[0],m_gradient[2],1,m_Forcing[1], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[1],m_gradient[3],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[1], 1, tmp, 1, m_Forcing[1], 1);
                break;
            case 3:
                //How do I deal with Homogeneous expansion?
                // - cu * du/dx - cv * du/dy - cz * du/dz
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[0],m_gradient[0],1,m_Forcing[0], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[1],m_gradient[1],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[0], 1, tmp, 1, m_Forcing[0], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[2],m_gradient[2],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[0], 1, tmp, 1, m_Forcing[0], 1);
                // - cu * dv/dx - cv * dv/dy  - cz * dv/dz
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[0],m_gradient[3],1,m_Forcing[1], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[1],m_gradient[4],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[1], 1, tmp, 1, m_Forcing[1], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[2],m_gradient[5],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[1], 1, tmp, 1, m_Forcing[1], 1);
                // - cu * dw/dx - cv * dw/dy  - cz * dw/dz
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[0],m_gradient[6],1,m_Forcing[2], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[1],m_gradient[7],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[2], 1, tmp, 1, m_Forcing[2], 1);
                Vmath::Smul(npoints,-1.0 * m_frameVelocity[2],m_gradient[8],1,tmp, 1);
                Vmath::Vadd(npoints, m_Forcing[2], 1, tmp, 1, m_Forcing[2], 1);
                break;
            default:
                ASSERTL0(false,"dimension unknown");
        }

        // If singleMode or halfMode, transform the forcing term to be in
        // physical space in the plane, but Fourier space in the homogeneous
        // direction
        if (m_transform)
        {
            for (int i = 0; i < m_NumVariable; ++i)
            {
                pFields[0]->HomogeneousFwdTrans(m_Forcing[i], m_Forcing[i]);
            }
        }
    }


    void ForcingMovingFrame::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        Update(fields, time);

        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }
    }

}
}

