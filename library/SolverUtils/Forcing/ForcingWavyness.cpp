///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingWavyness.cpp
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
// Description: Implementation of a forcing term to simulate a wavy body
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingWavyness.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace SolverUtils
{

std::string ForcingWavyness::className =
        GetForcingFactory().RegisterCreatorFunction("Wavyness",
                ForcingWavyness::create, "Wavyness Forcing");

/**
 *
 */
ForcingWavyness::ForcingWavyness(
        const LibUtilities::SessionReaderSharedPtr &pSession)
    : Forcing(pSession)
{
}


/**
 *
 */
void ForcingWavyness::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int                                &pNumForcingFields,
        const TiXmlElement                                *pForce)
{
    // Just 3D homogenous 1D problems can use this techinque
    ASSERTL0(pFields[0]->GetExpType() == MultiRegions::e3DH1D,
             "Wavyness implemented just for 3D Homogenous 1D expansions.");

    m_session->LoadParameter("Kinvis", m_kinvis);
    // forcing size (it must be 3)
    m_NumVariable       = pNumForcingFields;
    m_Forcing           = Array<OneD, Array<OneD, NekDouble> >(m_NumVariable);
    m_wavyGeom = Array<OneD, Array<OneD, NekDouble> >(6);
    int phystot         = pFields[0]->GetTotPoints();

    // Allocation of forcing and geometry memory
    for (int i = 0; i < m_NumVariable; ++i)
    {
        m_Forcing[i] = Array<OneD, NekDouble>(phystot, 0.0);
    }
    for (int i = 0; i < m_wavyGeom.num_elements(); i++)
    {
        m_wavyGeom[i] = Array<OneD, NekDouble>(phystot, 0.0);
    }

    const TiXmlElement* funcNameElmt;
    funcNameElmt = pForce->FirstChildElement("WAVYNESS");
    ASSERTL0(funcNameElmt, "Requires WAVYNESS tag, specifying function "
            "name which prescribes wavy function.");

    string funcName = funcNameElmt->GetText();
    ASSERTL0(m_session->DefinesFunction(funcName),
            "Function '" + funcName + "' not defined.");

    std::string s_FieldStr = m_session->GetVariable(0);
    ASSERTL0(m_session->DefinesFunction(funcName, s_FieldStr),
            "Variable '" + s_FieldStr + "' not defined.");

    EvaluateFunction(pFields, m_session, s_FieldStr, m_wavyGeom[0],
            funcName);

    // Calculate derivatives of transformation
    for (int i = 1; i < 4; i++)
    {
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                m_wavyGeom[i - 1], m_wavyGeom[i]);
    }

    Vmath::Vmul(phystot, m_wavyGeom[1], 1, m_wavyGeom[1], 1,
                         m_wavyGeom[4], 1);
    Vmath::Vmul(phystot, m_wavyGeom[1], 1, m_wavyGeom[2], 1,
                         m_wavyGeom[5], 1);
}


/**
 *
 */
void ForcingWavyness::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        const Array<OneD, Array<OneD, NekDouble> >          &inarray,
              Array<OneD, Array<OneD, NekDouble> >          &outarray,
        const NekDouble                                     &time)
{
    //calcualte the forcing components Ax,Ay,Az and put them in m_Forcing
    CalculateForcing(fields);

    // Apply forcing terms
    for (int i = 0; i < m_NumVariable; i++)
    {
        Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1, m_Forcing[i], 1,
                                                outarray[i], 1);
    }
}


/**
 *
 */
void ForcingWavyness::CalculateForcing(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    int np = fields[0]->GetNpoints();

    Array<OneD, NekDouble> U, V, W, P, tmp1, tmp2, tmp3, Wz, Wzz, Px;

    U = Array<OneD, NekDouble>(np);
    V = Array<OneD, NekDouble>(np);
    W = Array<OneD, NekDouble>(np);
    P = Array<OneD, NekDouble>(np);

    tmp1 = Array<OneD, NekDouble>(np);
    tmp2 = Array<OneD, NekDouble>(np);
    tmp3 = Array<OneD, NekDouble>(np);

    Wz = Array<OneD, NekDouble>(np);
    Wzz = Array<OneD, NekDouble>(np);
    Px = Array<OneD, NekDouble>(np);

    fields[0]->HomogeneousBwdTrans(fields[0]->GetPhys(), U);
    fields[0]->HomogeneousBwdTrans(fields[1]->GetPhys(), V);
    fields[0]->HomogeneousBwdTrans(fields[2]->GetPhys(), W);
    fields[0]->HomogeneousBwdTrans(fields[3]->GetPhys(), P);
    //-------------------------------------------------------------------------
    // Ax calculation
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], P, Px); // Px
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[3]->GetPhys(), tmp2); // Pz
    // Pz in physical space
    fields[0]->HomogeneousBwdTrans(tmp2, tmp3);
    // Pz * Xz
    Vmath::Vmul(np, tmp3,         1, m_wavyGeom[1], 1, tmp3,         1);
    // Px * Xz^2
    Vmath::Vmul(np, Px,           1, m_wavyGeom[4], 1, tmp1,         1);
    // W^2
    Vmath::Vmul(np, W,            1, W,             1, tmp2,         1);
    // Xzz * W^2
    Vmath::Vmul(np, tmp2,         1, m_wavyGeom[2], 1, tmp2,         1);
    // Pz * Xz - Px * Xz^2
    Vmath::Vsub(np, tmp3,         1, tmp1,          1, m_Forcing[0], 1);
    // A0 = Pz * Xz - Px * Xz^2 - Xzz * W^2
    Vmath::Vsub(np, m_Forcing[0], 1, tmp2,          1, m_Forcing[0], 1);

    // here part to be multiplied by 1/Re we use P to store it, since we dont't
    // need it anymore
    // W * Xzzz
    Vmath::Vmul(np, W,            1, m_wavyGeom[3], 1, P,            1);
    // P = W * Xzzz + Xz * Xzz
    Vmath::Vadd(np, P,            1, m_wavyGeom[5], 1, P,            1);
    // Wz
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[2]->GetPhys(), tmp1);
    // Wzz
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], tmp1, tmp2);
    //Wz in physical space
    fields[0]->HomogeneousBwdTrans(tmp1, Wz);
    //Wzz in physical space
    fields[0]->HomogeneousBwdTrans(tmp2, Wzz);
    // Wzz * Xz
    Vmath::Vmul(np, Wzz,          1, m_wavyGeom[1], 1, tmp1,         1);
    // P = W * Xzzz + Xz * Xzz - Wzz * Xz
    Vmath::Vsub(np, P,            1, tmp1,          1, P,            1);
    // Ux
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], U, tmp1);
    // Ux * Xzz
    Vmath::Vmul(np, tmp1,         1, m_wavyGeom[2], 1, tmp2,         1);
    // P = W * Xzzz + Xz * Xzz - Wzz * Xz - Ux * Xzz
    Vmath::Vsub(np, P,            1, tmp2,          1, P,            1);
    // Uxx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2);
    // Uxx * Xz^2
    Vmath::Vmul(np, tmp2,         1, m_wavyGeom[4], 1, tmp2,         1);
    // P = W * Xzzz + Xz * Xzz - Wzz * Xz - Ux * Xzz + Uxx * Xz^2
    Vmath::Vadd(np, P,            1, tmp2,          1, P,            1);
    // Uz
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[0]->GetPhys(), tmp1);
    // Uzx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2);
    // Uzx in physical space
    fields[0]->HomogeneousBwdTrans(tmp2, tmp1);
    // Uzx * Xz
    Vmath::Vmul(np, tmp1,         1, m_wavyGeom[1], 1, tmp1,         1);
    // 2 * Uzx * Xz
    Vmath::Smul(np, 2.0,             tmp1,          1, tmp2,         1);
    // P = W * Xzzz + Xz * Xzz - Wzz * Xz - Ux * Xzz + Uxx * Xz^2 - 2 * Uzx * Xz
    Vmath::Vsub(np, P,            1, tmp2,          1, P,            1);
    // *1/Re
    Vmath::Smul(np, m_kinvis,        P,             1, tmp1,         1);
    Vmath::Vadd(np, tmp1,         1, m_Forcing[0],  1, tmp2,         1);
    // back to Fourier Space
    fields[0]->HomogeneousFwdTrans(tmp2, m_Forcing[0]);

    //-------------------------------------------------------------------------
    // Ay calucaltion

    // Vx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], V, tmp1);
    // Vxx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2);
    // Vx * Xzz
    Vmath::Vmul(np, tmp1,         1, m_wavyGeom[2], 1, tmp1,         1);
    // Vxx * Xz^2
    Vmath::Vmul(np, tmp2,         1, m_wavyGeom[4], 1, tmp2,         1);
    // Vxx * Xz^2 - Vx* Xzz
    Vmath::Vsub(np, tmp2,         1, tmp1,          1, m_Forcing[1], 1);
    //Vz
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[1]->GetPhys(), tmp3);
    //Vzx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp3, tmp2);
    // Vzx physical space
    fields[0]->HomogeneousBwdTrans(tmp2, tmp1);
    // Vzx * Xz
    Vmath::Vmul(np, tmp1,         1, m_wavyGeom[1], 1, tmp2,         1);
    // 2 * Vzx * Xz
    Vmath::Smul(np, 2.0,             tmp2,          1, tmp1,         1);
    // Vxx * Xz^2 - Vx* Xzz - Vxz * Xz
    Vmath::Vsub(np, m_Forcing[1], 1, tmp1,          1, tmp3,         1);
    // * 1/Re
    Vmath::Smul(np, m_kinvis,        tmp3,          1, tmp1,         1);
    // back to Fourier Space
    fields[0]->HomogeneousFwdTrans(tmp1, m_Forcing[1]);

    //-------------------------------------------------------------------------
    // Az calculation
    // Px * Xz
    Vmath::Vmul(np, Px,           1, m_wavyGeom[1], 1, P,            1);
    // Wzz - Xzz
    Vmath::Vsub(np, Wzz,          1, m_wavyGeom[2], 1, Wzz,          1);
    // Wx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], W, tmp1);
    // Wx * Xzz
    Vmath::Vmul(np, tmp1,         1, m_wavyGeom[2], 1, tmp2,         1);
    // Wzz - Xzz - Wx * Xzz
    Vmath::Vsub(np, Wzz,          1, tmp2,          1, Wzz,          1);
    // Wxx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2);
    // Wxx * Xz^2
    Vmath::Vmul(np, tmp2,         1, m_wavyGeom[4], 1, tmp2,         1);
    // Wzz - Xzz - Wx * Xzz + Wxx * Xz^2
    Vmath::Vadd(np, Wzz,          1, tmp2,          1, Wzz,          1);
    // Wzx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Wz, tmp1);
    // Wzx * Xz
    Vmath::Vmul(np, tmp1,         1, m_wavyGeom[1], 1, tmp1,         1);
    // 2 * Vzx * Xz
    Vmath::Smul(np, 2.0,             tmp1,          1, tmp1,         1);
    // Wzz - Xzz - Wx * Xzz + Wxx * Xz^2 - 2 * Vzx * Xz
    Vmath::Vsub(np, Wzz,          1, tmp1,          1, Wzz,          1);
    // * 1/Re
    Vmath::Smul(np, m_kinvis,        Wzz,           1, tmp1,         1);
    Vmath::Vadd(np, tmp1,         1, P,             1, tmp2,         1);
    // back to Fourier Space
    fields[0]->HomogeneousFwdTrans(tmp2, m_Forcing[2]);
}

}
}
