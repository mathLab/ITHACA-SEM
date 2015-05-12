///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingBody.cpp
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
// Description: Moving Body m_forcing (movement of a body in a domain is achieved
// via a m_forcing term which is the results of a coordinate system motion)
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingMovingBody.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
std::string ForcingMovingBody::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("MovingBody",
                                    ForcingMovingBody::create,
                                    "Moving Body Forcing");

ForcingMovingBody::ForcingMovingBody(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : Forcing(pSession)
{
}

void ForcingMovingBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pfields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pforce)
{
    // Just 3D homogenous 1D problems can use this techinque
    ASSERTL0(pfields[0]->GetExpType()==MultiRegions::e3DH1D,
             "Moving body implemented just for 3D Homogenous 1D expansions.");

    // At this point we know in the xml file where those quantities
    // are declared (equation or file) - via a function name which is now
    // stored in funcNameD etc. We need now to fill in with this info the
    // m_zta and m_eta vectors (actuallythey are matrices) Array to control if
    // the motion is determined by an equation or is from a file.(not Nektar++)
    // check if we need to load a file or we have an equation
    CheckIsFromFile(pforce);

    // Initialise movingbody filter
    InitialiseFilter(m_session, pfields, pforce);

    // Initialise the cable model
    InitialiseCableModel(m_session, pfields);

    m_pressureCalls = 0;

    m_zta = Array<OneD, Array< OneD, NekDouble> > (10);
    m_eta = Array<OneD, Array< OneD, NekDouble> > (10);
    // What are this bi-dimensional vectors ------------------------------------------
    // m_zta[0] = m_zta                     |  m_eta[0] = m_eta                      |
    // m_zta[1] = d(m_zta)/dt               |  m_eta[1] = d(m_eta)/dt                |
    // m_zta[2] = dd(m_zta)/ddtt            |  m_eta[2] = dd(m_eta)/ddtt             |
    // m_zta[3] = d(m_zta)/dz               |  m_eta[3] = d(m_eta)/dz                |
    // m_zta[4] = dd(m_zta)/ddzz            |  m_eta[4] = dd(m_eta)/ddzz             |
    // m_zta[5] = ddd(m_zta)/dddzzz         |  m_eta[5] = ddd(m_eta)/dddzzz          |
    // m_zta[6] = dd(m_zta)/ddtz            |  m_eta[6] = dd(m_eta)/ddtz             |
    // m_zta[7] = ddd(m_zta)/dddtzz         |  m_eta[7] = ddd(m_eta)/dddtzz          |
    // m_zta[8] = (d(m_zta)/dz)(d(m_eta)/dz)|  m_eta[8] = (d(m_eta)/dz)(d(m_zta)/dz) |
    // m_zta[9] = (d(zm_zta)/dz)^2          |  m_eta[9] = (d(m_eta)/dz)^2            |
    //--------------------------------------------------------------------------------
    int phystot = pfields[0]->GetTotPoints();
    for(int i = 0; i < m_zta.num_elements(); i++)
    {
        m_zta[i] = Array<OneD, NekDouble>(phystot,0.0);
        m_eta[i] = Array<OneD, NekDouble>(phystot,0.0);
    }

    // m_forcing size (it must be 3: Ax, Ay and Az)
    // total number of structural variables(Disp, Vel and Accel) must be 3
    m_forcing = Array<OneD, Array<OneD, NekDouble> >(3);
    // create the storage space for the body motion description

    for(int i = 0; i < 3; i++)
    {
        m_forcing[i] = Array<OneD, NekDouble> (phystot, 0.0);
    }
 
}

void ForcingMovingBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  fields,
        const Array<OneD, Array<OneD, NekDouble> >&         inarray,
              Array<OneD, Array<OneD, NekDouble> >&         outarray,
        const NekDouble&                                    time)
{
    // Fill in m_zta and m_eta with all the required values
    UpdateMotion(fields,inarray,time);

    //calcualte the m_forcing components Ax,Ay,Az and put them in m_forcing
    CalculateForcing(fields);
    // Apply m_forcing terms
    for (int i = 0; i < 3; i++)
    {
        Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                               m_forcing[i], 1, outarray[i], 1);
    }

	// Mapping boundary conditions
    MappingBndConditions(fields, inarray, time);
}

/**
 *
 */
void ForcingMovingBody::UpdateMotion(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  pfields,
        const Array<OneD, Array<OneD, NekDouble> >       &  fields,
              NekDouble                                     time)
{
    // Update the forces from the calculation of fluid field, which is
    // implemented in the movingbody filter
    m_MovBodyfilter->UpdateForce(m_session, pfields, m_Aeroforces, time);

    // for "free" type, the cable vibrates both in streamwise and crossflow
    // dimections, for "constrained" type, the cable only vibrates in crossflow
    // dimection, and for "forced" type, the calbe vibrates specifically along
    // a given function or file
    std::string vibtype = m_session->GetSolverInfo("VibrationType"); 
    if(vibtype == "Free" || vibtype == "FREE" ||
       vibtype == "Constrained" || vibtype == "CONSTRAINED")
    {
        // For free vibration case, displacements, velocities and acceleartions
        // are obtained through solving structure dynamic model
        EvaluateStructDynModel(pfields, time);
    }
    else if(vibtype == "Forced" || vibtype == "FORCED")
    {
        // For forced vibration case, load from specified file or function
        int cnt = 0;
        for(int j = 0; j < m_funcName.num_elements(); j++)
        {
            if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
            {
                ASSERTL0(false, "Motion loading from file needs specific "
                                "implementation: Work in Progress!");
            }
            else
            {
                EvaluateFunction(pfields, m_session, m_motion[0], m_zta[j],
                                 m_funcName[j], time);
                EvaluateFunction(pfields, m_session, m_motion[1], m_eta[j],
                                 m_funcName[j], time);
                cnt = cnt + 2;
            }
        }

        for(int var = 0; var < 3; var++)
        {
            for(int plane = 0; plane < m_np; plane++)
            {
                int n = pfields[0]->GetPlane(plane)->GetTotPoints();
                int offset  = plane * n;
                int xoffset = var * m_np+plane;
                int yoffset = 3*m_np+xoffset;

                m_MotionVars[xoffset] = m_zta[var][offset];
                m_MotionVars[yoffset] = m_eta[var][offset];
            }
        }
    }
    else
    {
        ASSERTL0(false, 
                 "Unrecogized vibration type for cable's dynamic model");
    }

    // Pass the variables of the cable's motion to the movingbody filter
    m_MovBodyfilter->UpdateMotion(m_session, pfields, m_MotionVars, time);

    // Now we need to calcualte all the required z-derivatives
    bool OriginalWaveSpace = pfields[0]->GetWaveSpace();
    pfields[0]->SetWaveSpace(false);

    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zta[0], m_zta[3]); //d(m_zta)/dz
    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zta[3], m_zta[4]); //dd(m_zta)/ddzz
    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zta[4], m_zta[5]); //ddd(m_zta)/dddzzz

    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[0], m_eta[3]); //d(m_eta)/dz
    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[3], m_eta[4]); //dd(m_eta)/ddzz
    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[4], m_eta[5]); //ddd(m_eta)/dddzzz

    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zta[1], m_zta[6]); //dd(m_zta)/ddtz
    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zta[6], m_zta[7]); //ddd(m_zta)/ddtzz

    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[1], m_eta[6]); //dd(m_eta)/ddtz
    pfields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[6], m_eta[7]); //ddd(m_eta)/ddtzz

    int NumPoints = pfields[0]->GetTotPoints();

    // (d(m_zta)/dz)(d(m_eta)/dz)
    Vmath::Vmul(NumPoints, m_zta[3], 1, m_eta[3], 1, m_zta[8], 1);
    // (d(m_eta)/dz)(d(m_zta)/dz) // not really needed
    Vmath::Vmul(NumPoints, m_eta[3], 1, m_zta[3], 1, m_eta[8], 1);

    // (d(m_zta)/dz)^2
    Vmath::Vmul(NumPoints, m_zta[3], 1, m_zta[3], 1, m_zta[9], 1);
    // (d(m_eta)/dz)^2
    Vmath::Vmul(NumPoints, m_eta[3], 1, m_eta[3], 1, m_eta[9], 1);

    pfields[0]->SetWaveSpace(OriginalWaveSpace);
}


/**
 *
 */
void ForcingMovingBody::CalculateForcing(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{

    int nPointsTot = fields[0]->GetNpoints();
    Array<OneD, NekDouble> U,V,W;
    Array<OneD, NekDouble> Wt,Wx,Wy,Wz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz;
    Array<OneD, NekDouble> Px,Py,Pz;
    Array<OneD, NekDouble> tmp,tmp1,tmp2,tmp3;
    Array<OneD, NekDouble> Fx,Fy,Fz;
    U    = Array<OneD, NekDouble> (nPointsTot);
    V    = Array<OneD, NekDouble> (nPointsTot);
    W    = Array<OneD, NekDouble> (nPointsTot);
    Wt   = Array<OneD, NekDouble> (nPointsTot);
    Wx   = Array<OneD, NekDouble> (nPointsTot);
    Wy   = Array<OneD, NekDouble> (nPointsTot);
    Wz   = Array<OneD, NekDouble> (nPointsTot);
    Wxx  = Array<OneD, NekDouble> (nPointsTot);
    Wxy  = Array<OneD, NekDouble> (nPointsTot);
    Wyy  = Array<OneD, NekDouble> (nPointsTot);
    Wxz  = Array<OneD, NekDouble> (nPointsTot);
    Wyz  = Array<OneD, NekDouble> (nPointsTot);
    Wzz  = Array<OneD, NekDouble> (nPointsTot);
    Px   = Array<OneD, NekDouble> (nPointsTot);
    Py   = Array<OneD, NekDouble> (nPointsTot);
    Pz   = Array<OneD, NekDouble> (nPointsTot);
    Fx   = Array<OneD, NekDouble> (nPointsTot,0.0);
    Fy   = Array<OneD, NekDouble> (nPointsTot,0.0);
    Fz   = Array<OneD, NekDouble> (nPointsTot,0.0);
    tmp  = Array<OneD, NekDouble> (nPointsTot,0.0);
    tmp1 = Array<OneD, NekDouble> (nPointsTot);
    tmp2 = Array<OneD, NekDouble> (nPointsTot);
    tmp3 = Array<OneD, NekDouble> (nPointsTot);
    fields[0]->HomogeneousBwdTrans(fields[0]->GetPhys(),U);
    fields[0]->HomogeneousBwdTrans(fields[1]->GetPhys(),V);
    fields[0]->HomogeneousBwdTrans(fields[2]->GetPhys(),W);
    fields[0]->HomogeneousBwdTrans(fields[3]->GetPhys(),tmp1); // pressure

    //-------------------------------------------------------------------------
    // Setting the pressure derivatives
    //-------------------------------------------------------------------------
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],tmp1,Px); // Px
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp1,Py); // Py

    //-------------------------------------------------------------------------
    // Setting the z-component velocity derivatives
    //-------------------------------------------------------------------------
    EvaluateAccelaration(W,Wt,nPointsTot); //Wt
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], W, Wx); // Wx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], W, Wy); // Wy
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[2]->GetPhys(),tmp1); // Wz
    fields[0]->HomogeneousBwdTrans(tmp1, Wz); // Wz in physical space
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Wx, Wxx); // Wxx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Wx, Wxy); // Wxy
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Wy, Wyy); // Wyy
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Wz, Wxz); // Wxz
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Wz, Wyz); // Wyz
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], tmp1, tmp2); // Wzz
    fields[0]->HomogeneousBwdTrans(tmp2, Wzz); // Wzz in physical space
    //-------------------------------------------------------------------------
    // Setting the z-component of the accelaration term:
    //   tmp = (Wt + U * Wx + V * Wy + W * Wz)
    //-------------------------------------------------------------------------
    Vmath::Vadd(nPointsTot, tmp, 1, Wt,   1, tmp,  1);
    Vmath::Vmul(nPointsTot, U,   1, Wx,   1, tmp1, 1);
    Vmath::Vadd(nPointsTot, tmp, 1, tmp1, 1, tmp,  1);
    Vmath::Vmul(nPointsTot, V,   1, Wy,   1, tmp1, 1);
    Vmath::Vadd(nPointsTot, tmp, 1, tmp1, 1, tmp,  1);
    Vmath::Vmul(nPointsTot, W,   1, Wz,   1, tmp1, 1);
    Vmath::Vadd(nPointsTot, tmp, 1, tmp1, 1, tmp,  1);

    //-------------------------------------------------------------------------
    // x-component of the m_forcing - accelaration component
    //-------------------------------------------------------------------------
    Vmath::Vsub(nPointsTot, Fx, 1, m_zta[2], 1, Fx,   1);
    Vmath::Vmul(nPointsTot, m_zta[3], 1, tmp,1, tmp1, 1);
    Vmath::Vsub(nPointsTot, Fx, 1, tmp1,   1, Fx,   1);

    Vmath::Vmul(nPointsTot, m_zta[6], 1, W,  1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,    tmp1,  1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Fx, 1, tmp2,   1, Fx,   1);
    // W^2 - we reuse it later
    Vmath::Vmul(nPointsTot, W,   1, W,      1, tmp3, 1);
    Vmath::Vmul(nPointsTot, m_zta[4], 1, tmp3,1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Fx, 1, tmp2,    1, Fx,   1);

    //-------------------------------------------------------------------------
    // y-component of the m_forcing - accelaration component
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_eta[4],  1, tmp3, 1, tmp2, 1); // reusing W^2
    Vmath::Vsub(nPointsTot, Fy, 1, tmp2,      1, Fy,   1);

    Vmath::Vmul(nPointsTot, m_eta[3],  1, tmp,  1, tmp1, 1);
    Vmath::Vsub(nPointsTot, Fy, 1, tmp1,      1, Fy,   1);

    Vmath::Vsub(nPointsTot, Fy,   1, m_eta[2],  1, Fy,   1);
    Vmath::Vmul(nPointsTot, m_eta[6],  1, W,    1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,    tmp1,     1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Fy, 1, tmp2,      1, Fy,   1);

    //-------------------------------------------------------------------------
    // z-component of the m_forcing - accelaration component
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_zta[3], 1, Px,  1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3], 1, Py,  1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Fz,   1, tmp1,  1, Fz,   1);
    Vmath::Vadd(nPointsTot, Fz,   1, tmp2,  1, Fz,   1);

    //-------------------------------------------------------------------------
    // Note: Now we use Px,Py,Pz to store the viscous component of the m_forcing
    // before we multiply them by m_kinvis = 1/Re. Since we build up on them we
    // need to set the entries to zero.
    //-------------------------------------------------------------------------
    Vmath::Zero(nPointsTot,Px,1);
    Vmath::Zero(nPointsTot,Py,1);
    Vmath::Zero(nPointsTot,Pz,1);

    //-------------------------------------------------------------------------
    // x-component of the m_forcing - viscous component1:  (U_z^'z^' - U_zz + ZETA_tz^'z^')
    //-------------------------------------------------------------------------
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[0]->GetPhys(), tmp1); // Uz
    fields[0]->HomogeneousBwdTrans(tmp1, tmp3); // Uz in physical space
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp3, tmp1); // Uzx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp3, tmp2); // Uzy

    Vmath::Vmul(nPointsTot, m_zta[3], 1, tmp1,      1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp1,      1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3], 1, tmp2,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp2,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp1,      1, Px,   1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp2,      1, Px,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], U,    tmp1); // Ux
    Vmath::Vmul(nPointsTot, m_zta[4], 1, tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp2,      1, Px,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2); // Uxx
    Vmath::Vmul(nPointsTot, m_zta[9], 1, tmp2,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,     1, tmp3,      1, Px,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Uxy
    Vmath::Vmul(nPointsTot, m_zta[8], 1, tmp2,      1, tmp3, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp3,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,     1, tmp3,      1, Px,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], U,    tmp1); // Uy
    Vmath::Vmul(nPointsTot, m_eta[4],  1, tmp1,     1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,      1, tmp2,     1, Px,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Uyy
    Vmath::Vmul(nPointsTot, m_eta[9],  1, tmp2,     1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,      1, tmp3,     1, Px,   1);
    Vmath::Vadd(nPointsTot, Px,      1, m_zta[7],   1, Px,   1);

    //-------------------------------------------------------------------------
    // x-component of the m_forcing - viscous component2:
    //      ((ZETA_z * W)_z^'z^' + ZETA_z * (W_xx + W_yy))
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_zta[5], 1, W,         1, tmp1, 1);
    Vmath::Vadd(nPointsTot, Px,     1, tmp1,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[3], 1, Wx,        1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_zta[4], 1, tmp1,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 3.0,       tmp2,      1, tmp3, 1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp3,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_eta[3], 1, Wy,        1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_zta[4], 1, tmp1,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp2,      1, tmp3, 1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp3,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[3], 1, Wy,        1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[4], 1, tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[4], 1, Wz,        1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp1,      1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Px,     1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[9], 1, Wxy,       1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3], 1, tmp1,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp2,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,     1, tmp3,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[9], 1, Wxz,       1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[8], 1, Wyz,       1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,       tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,     1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[9], 1, m_zta[3],    1, tmp1, 1);
    Vmath::Vmul(nPointsTot, tmp1,1, Wxx,          1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Px,  1, tmp2,         1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zta[3], 1, Wyy,       1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[9], 1, tmp1,      1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Px,     1, tmp2,      1, Px,   1);

    Vmath::Vadd(nPointsTot, Wxx,    1, Wyy,       1, tmp1, 1);
    Vmath::Vadd(nPointsTot, Wzz,    1, tmp1,      1, tmp2, 1);
    Vmath::Vmul(nPointsTot, m_zta[3], 1, tmp2,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,     1, tmp3,      1, Px,   1);

    Vmath::Smul(nPointsTot, m_kinvis,  Px,        1, Px,   1); //* 1/Re

    //-------------------------------------------------------------------------
    // y-component of the m_forcing - viscous component1:
    //  (V_z^'z^' - V_zz + ETA_tz^'z^')
    //-------------------------------------------------------------------------
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[1]->GetPhys(), tmp1); // Vz
    fields[0]->HomogeneousBwdTrans(tmp1, tmp3); // Vz in physical space
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp3, tmp1); // Vzx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp3, tmp2); // Vzy
    Vmath::Vmul(nPointsTot, m_zta[3],  1,  tmp1,   1,  tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,         tmp1,   1,  tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3],  1,  tmp2,   1,  tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,         tmp2,   1,  tmp2, 1);
    Vmath::Vsub(nPointsTot, Py,      1,  tmp1,   1,  Py,   1);
    Vmath::Vsub(nPointsTot, Py,      1,  tmp2,   1,  Py,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], V,    tmp1); // Vx
    Vmath::Vmul(nPointsTot, m_zta[4],  1,  tmp1,   1,  tmp2, 1);
    Vmath::Vsub(nPointsTot, Py,      1,  tmp2,   1,  Py,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2); // Vxx
    Vmath::Vmul(nPointsTot, m_zta[9],  1,  tmp2,   1,  tmp3, 1);
    Vmath::Vadd(nPointsTot, Py,      1,  tmp3,   1,  Py,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Vxy
    Vmath::Vmul(nPointsTot, m_zta[8],  1,  tmp2,   1,  tmp3, 1);
    Vmath::Smul(nPointsTot, 2.0,         tmp3,   1,  tmp3, 1);
    Vmath::Vadd(nPointsTot, Py,      1,  tmp3,   1,  Py,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], V,    tmp1); // Vy
    Vmath::Vmul(nPointsTot, m_eta[4],   1,  tmp1,  1,  tmp2, 1);
    Vmath::Vsub(nPointsTot, Py,       1,  tmp2,  1,  Py,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Vyy
    Vmath::Vmul(nPointsTot, m_eta[9],   1,  tmp2,  1,  tmp3, 1);
    Vmath::Vadd(nPointsTot, Py,       1,  tmp3,  1,  Py,   1);
    Vmath::Vadd(nPointsTot, m_eta[7],   1,  Py,    1,  Py,   1);

    //-------------------------------------------------------------------------
    // y-component of the m_forcing - viscous component2:
    //   ((ETA_z * W)_z^'z^' + ETA_z * (W_xx + W_yy))
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_eta[5],   1,  W,     1,  tmp1, 1);
    Vmath::Vadd(nPointsTot, Py,       1,  tmp1,  1,  Py,   1);

    Vmath::Vmul(nPointsTot, m_zta[3],   1,  Wx,    1,  tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[4],   1,  tmp1,  1,  tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp2,  1,  tmp3, 1);
    Vmath::Vsub(nPointsTot, Py,       1,  tmp3,  1,  Py,   1);

    Vmath::Vmul(nPointsTot, m_zta[4],  1,  Wx,     1,  tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3],  1,  tmp1,   1,  tmp2, 1);
    Vmath::Vsub(nPointsTot, Py,      1,  tmp2,   1,  Py,   1);

    Vmath::Vmul(nPointsTot, m_eta[3],   1,  Wy,    1,  tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[4],   1,  tmp1,  1,  tmp2, 1);
    Vmath::Smul(nPointsTot, 3.0,          tmp2,  1,  tmp3, 1);
    Vmath::Vsub(nPointsTot, Py,       1,  tmp3,  1,  Py,   1);

    Vmath::Vmul(nPointsTot, m_eta[4],   1,  Wz,    1,  tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,     tmp1,       1,  tmp2, 1);
    Vmath::Vadd(nPointsTot, Py,  1,  tmp2,       1,  Py,   1);

    Vmath::Vmul(nPointsTot, m_eta[9],   1,  Wxy, 1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_zta[3],  1,  tmp1, 1,  tmp2,   1);
    Vmath::Smul(nPointsTot, 2.0,   tmp2,       1,  tmp3,   1);
    Vmath::Vadd(nPointsTot, Py, 1, tmp3,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_zta[8],  1,  Wxz,  1,  tmp1,   1);
    Vmath::Smul(nPointsTot, 2.0,   tmp1,       1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Py, 1,  tmp2,      1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[9],   1,  Wyz, 1,  tmp1,   1);
    Vmath::Smul(nPointsTot, 2.0,   tmp1,       1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Py, 1,  tmp2,      1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[3],   1,  Wxx, 1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_zta[9],  1,  tmp1, 1,  tmp2,   1);
    Vmath::Vadd(nPointsTot, Py, 1,  tmp2,      1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[9], 1,  m_eta[3], 1,  tmp1,  1);
    Vmath::Vmul(nPointsTot, tmp1,   1,  Wyy,    1,  tmp2,  1);
    Vmath::Vadd(nPointsTot, Py,     1,  tmp2,   1,  Py,    1);

    Vmath::Vadd(nPointsTot, Wxx,    1,  Wyy,    1,  tmp1,  1);
    Vmath::Vadd(nPointsTot, Wzz,    1,  tmp1,   1,  tmp2,  1);
    Vmath::Vmul(nPointsTot, m_eta[3], 1,  tmp2,   1,  tmp3,  1);
    Vmath::Vadd(nPointsTot, Py,     1,  tmp3,   1,  Py,    1);

    Vmath::Smul(nPointsTot, m_kinvis,   Py,     1,  Py,    1); //* 1/Re

    //-------------------------------------------------------------------------
    // z-component of the m_forcing - viscous component: (W_z^'z^' - W_zz)
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_zta[3],  1,  Wxz,   1,  tmp1,  1);
    Vmath::Smul(nPointsTot, 2.0,    tmp1,       1,  tmp1,  1);
    Vmath::Vmul(nPointsTot, m_eta[3],   1,  Wyz,  1,  tmp2,  1);
    Vmath::Smul(nPointsTot, 2.0,    tmp2,       1,  tmp2,  1);
    Vmath::Vsub(nPointsTot, Pz, 1,  tmp1,       1,  Pz,    1);
    Vmath::Vsub(nPointsTot, Pz, 1,  tmp2,       1,  Pz,    1);

    Vmath::Vmul(nPointsTot, m_zta[4],  1,  Wx,    1,  tmp2,  1);
    Vmath::Vsub(nPointsTot, Pz,     1,  tmp2,   1,  Pz,    1);
    Vmath::Vmul(nPointsTot, m_zta[9],  1,  Wxx,   1,  tmp3,  1);
    Vmath::Vadd(nPointsTot, Pz,     1,  tmp3,   1,  Pz,    1);
    Vmath::Vmul(nPointsTot, m_zta[8],  1,  Wxy,   1,  tmp3,  1);
    Vmath::Smul(nPointsTot, 2.0,    tmp3,       1,  tmp3,  1);
    Vmath::Vadd(nPointsTot, Pz, 1,  tmp3,       1,  Pz,    1);

    Vmath::Vmul(nPointsTot, m_eta[4], 1,  Wy,     1,  tmp2,  1);
    Vmath::Vsub(nPointsTot, Pz,     1,  tmp2,   1,  Pz,    1);
    Vmath::Vmul(nPointsTot, m_eta[9],   1,  Wyy,  1,  tmp3,  1);
    Vmath::Vadd(nPointsTot, Pz,     1,  tmp3,   1,  Pz,    1);

    Vmath::Smul(nPointsTot, m_kinvis,   Pz,     1,  Pz,    1); //* 1/Re

    //-------------------------------------------------------------------------
    // adding viscous and pressure components and transfroming back to wave
    // space
    //-------------------------------------------------------------------------
    Vmath::Vadd(nPointsTot, Fx,     1,  Px,     1,  Fx,    1);
    Vmath::Vadd(nPointsTot, Fy,     1,  Py,     1,  Fy,    1);
    Vmath::Vadd(nPointsTot, Fz,     1,  Pz,     1,  Fz,    1);
    fields[0]->HomogeneousFwdTrans(Fx, m_forcing[0]);
    fields[0]->HomogeneousFwdTrans(Fy, m_forcing[1]);
    fields[0]->HomogeneousFwdTrans(Fz, m_forcing[2]);
}


/**
 *
 */
void ForcingMovingBody::TensionedCableModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields,
              Array<OneD, NekDouble> &AeroForces,
              Array<OneD, NekDouble> &CableMotions)
{  
    std::string supptype = m_session->GetSolverInfo("SupportType");

    bool homostrip;
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);
    // homostrip == false indicates that a single fourier transformation is
    // implemented for the motion of cable, so it is matched with that carried
    // out in fluid fields; on the other hand, if homostrip == true, there is
    // a mismatch between the structure and fluid fields in terms of fourier
    // transformation, then the routines such as FFTFwdTrans/
    // FFTBwdTrans can not be used dimectly for cable's solution.

    int npts;

    if(!homostrip) //full resolutions
    {
        npts = m_session->GetParameter("HomModesZ");
    }
    else
    {
        m_session->LoadParameter("Strip_Z", npts);
    }

    Array<OneD, Array<OneD, NekDouble> > fft_i(4);
    Array<OneD, Array<OneD, NekDouble> > fft_o(4);

    for(int i = 0; i < 4; i++)
    {
        fft_i[i] = Array<OneD, NekDouble>(npts, 0.0);
        fft_o[i] = Array<OneD, NekDouble>(npts, 0.0);
    }

    LibUtilities::CommSharedPtr vcomm = pfields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    int nproc   = vcomm->GetColumnComm()->GetSize();
    
    if(!homostrip) //full resolutions
    {
        Array<OneD, NekDouble> tmp(4*m_np);
        
        for(int n = 0; n < m_np; n++)
        {
            tmp[n] = AeroForces[n];
            for(int j = 0; j < 3; ++j)
            {
                tmp[(j+1)*m_np + n] = 
                    CableMotions[j*m_np + n];
            }
        }
        // Send to root process.
        if(colrank == 0)
        {
            Vmath::Vcopy(m_np, AeroForces, 1, fft_i[0], 1);

            for (int var = 0; var < 3; var++)
            {
                Vmath::Vcopy(m_np, CableMotions+var*m_np, 1, fft_i[var+1], 1);
            }

            for (int i = 1; i < nproc; ++i)
            {
                vcomm->GetColumnComm()->Recv(i, tmp);
                for(int n = 0; n < m_np; n++)
                {
                    for(int var = 0; var < 4; ++var)
                    {
                        fft_i[var][i*m_np + n] = tmp[var*m_np + n];
                    }
                }
            }
        }
        else
        {
            vcomm->GetColumnComm()->Send(0, tmp);
        }
    }
    else //strip modeling
    {
        Array<OneD, NekDouble> tmp(4);

        tmp[0] = AeroForces[0];
        for(int j = 0; j < 3; ++j)
        {
            tmp[j+1] = CableMotions[j*m_np];
        }

        // Send to root process.
        if(colrank == 0)
        {
            for(int j = 0 ; j < 4; ++j)
            {
                fft_i[j][0] = tmp[j];
            }

            for(int i = 1; i < npts; ++i)
            {
                vcomm->GetColumnComm()->Recv(i, tmp);

                for(int j = 0 ; j < 4; ++j)
                {
                    fft_i[j][i] = tmp[j];
                }
            }
        }
        else
        {
            for(int i = 1; i < npts; ++i)
            {
                if(colrank == i)
                {
                    vcomm->GetColumnComm()->Send(0, tmp);
                }
            }
        }
    }
    
    if(colrank == 0)
    {
        // Implement Fourier transformation of the motion variables
        if(supptype == "FREE-FREE" || 
                        supptype == "Free-Free")
        {
            for(int j = 0 ; j < 4; ++j)
            {
                m_FFT->FFTFwdTrans(fft_i[j], fft_o[j]);
            }
        }
        else if(supptype == "PINNED-PINNED" || 
                            supptype == "Pinned-Pinned")
        {
            //TODO:
            int N = fft_i[0].num_elements();

            for(int var = 0; var < 4; var++)
            {
                for(int i = 0; i < N; i++)
                {
                    fft_o[var][i] = 0;

                    for(int k = 0; k < N; k++)
                    {
                        fft_o[var][i] +=
                            fft_i[var][k]*
                                sin(M_PI/(N)*(k+1/2)*(i+1));
                    }
                }
            }
        }
        else
        {
            ASSERTL0(false,
                        "Unrecognized support type for cable's motion");
        }

        // solve the ODE in the wave space
        for(int i = 0; i < npts; ++i)
        {
            int nrows = 3;

            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            tmp0 = Array<OneD, NekDouble> (3,0.0);
            tmp1 = Array<OneD, NekDouble> (3,0.0);
            tmp2 = Array<OneD, NekDouble> (3,0.0);

            for(int var = 0; var < 3; var++)
            {
                tmp0[var] = fft_o[var+1][i];
            }

            tmp2[0] = fft_o[0][i];

            Blas::Dgemv('N', nrows, nrows, 1.0,
                        &(m_CoeffMat_B[i]->GetPtr())[0],
                        nrows, &tmp0[0], 1,
                        0.0,   &tmp1[0], 1);
            Blas::Dgemv('N', nrows, nrows, 1.0/m_structrho,
                        &(m_CoeffMat_A[i]->GetPtr())[0],
                        nrows, &tmp2[0], 1,
                        1.0,   &tmp1[0], 1);

            for(int var = 0; var < 3; var++)
            {
                fft_i[var][i] = tmp1[var];
            }
        }

        // get physical coeffients via Backward fourier transformation of wave
        // coefficients
        if(supptype == "FREE-FREE" || 
                        supptype == "Free-Free")
        {
            for(int var = 0; var < 3; var++)
            {
                m_FFT->FFTBwdTrans(fft_i[var], fft_o[var]);
            }
        }
        else if(supptype == "PINNED-PINNED" || 
                            supptype == "Pinned-Pinned")
        {
            //TODO:
            int N = fft_i[0].num_elements();

            for(int var = 0; var < 3; var++)
            {
                for(int i = 0; i < N; i++)
                {
                    fft_o[var][i] = 0;

                    for(int k = 0; k < N; k++)
                    {
                        fft_o[var][i] +=
                            fft_i[var][k]*
                                sin(M_PI/(N)*(k+1)*(i+1/2))*2/N;
                    }
                }
            }
        }
        else
        {
            ASSERTL0(false,
                        "Unrecognized support type for cable's motion");
        }
    }

    // send physical coeffients to all planes of each processor
    if(!homostrip)//full resolutions
    {
        Array<OneD, NekDouble> tmp(3*m_np);

        if(colrank != 0)
        {
            vcomm->GetColumnComm()->Recv(0, tmp);
            Vmath::Vcopy(3*m_np, tmp, 1, 
                            CableMotions, 1);
        }
        else
        {
            for (int i = 1; i < nproc; ++i)
            {
                for(int j = 0; j < 3; j++)
                {
                    for(int n = 0; n < m_np; n++)
                    {
                        tmp[j*m_np+n] = fft_o[j][i*m_np+n];
                    }
                }
                vcomm->GetColumnComm()->Send(i, tmp);
            }

            for(int j = 0; j < 3; j++)
            {
                for(int n = 0; n < m_np; n++)
                {
                    tmp[j*m_np+n] = fft_o[j][n];
                }
                Vmath::Vcopy(3*m_np, tmp, 1,  
                            CableMotions, 1);
            }
        }
    }
    else //strip modeling
    {
        Array<OneD, NekDouble> tmp(4);

        if(colrank != 0)
        {
            for (int j = 1; j < nproc/npts; j++)
            {
                if(colrank == j*npts)
                {
                    vcomm->GetColumnComm()->Recv(0, tmp);

                    for(int plane = 0; plane < m_np; plane++)
                    {
                        for(int var = 0; var < 3; var++)
                        {
                            CableMotions[var*m_np+plane]= tmp[var];
                        }
                    }
                }
            }

            for(int i = 1; i < npts; i++)
            {
                for (int j = 0; j < nproc/npts; j++)
                {
                    if(colrank == i+j*npts)
                    {
                        vcomm->GetColumnComm()->Recv(0, tmp);

                        for(int plane = 0; plane < m_np; plane++)
                        {
                            for(int var = 0; var < 3; var++)
                            {
                                CableMotions[var*m_np+plane] = tmp[var];
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for(int j = 0; j < 3; ++j)
            {
                tmp[j] = fft_o[j][0];
            }

            for (int j = 1; j < nproc/npts; j++)
            {
                vcomm->GetColumnComm()->Send(j*npts, tmp);
            }

            for(int plane = 0; plane < m_np; plane++)
            {
                for(int var = 0; var < 3; var++)
                {
                    CableMotions[var*m_np+plane] = tmp[var];
                }
            }

            for(int i = 1; i < npts; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    tmp[j] = fft_o[j][i];
                }

                for (int j = 0; j < nproc/npts; j++)
                {
                    vcomm->GetColumnComm()->Send(i+j*npts, tmp);
                }
            }
        }
    }
}


/**
 *
 */
void ForcingMovingBody::EvaluateStructDynModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields,
              NekDouble  time)
{
    //Get the hydrodynamic forces from the fluid solver
    Array<OneD, Array <OneD, NekDouble> > fces(m_motion.num_elements());

    fces[0] = Array <OneD, NekDouble> (m_np,0.0);
    fces[1] = Array <OneD, NekDouble> (m_np,0.0);

    for(int plane = 0; plane < m_np; plane++)
    {
        fces[0][plane] = m_Aeroforces[plane];
        fces[1][plane] = m_Aeroforces[plane+m_np];
    }

    // Fictitious mass method used to stablize the explicit coupling at
    // relatively lower mass ratio
    bool fictmass;
    m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                fictmass, false);
    if(fictmass)
    {
        NekDouble fictrho, fictdamp;
        m_session->LoadParameter("FictMass", fictrho);
        m_session->LoadParameter("FictDamp", fictdamp);

        static NekDouble Betaq_Coeffs[2][2] = 
                            {{1.0,  0.0},{2.0, -1.0}};

        // only consider second order approximation for fictitious variables
        int  intOrder= 2;
        int  nint    = min(m_movingBodyCalls+1,intOrder);
        int  nlevels = m_fV[0].num_elements();

        for(int i = 0; i < m_motion.num_elements(); ++i)
        {
            RollOver(m_fV[i]);
            RollOver(m_fA[i]);

            int Voffset = i*3*m_np+m_np;
            int Aoffset = i*3*m_np+2*m_np;

            Vmath::Vcopy(m_np, m_MotionVars+Voffset, 1, m_fV[i][0], 1);
            Vmath::Vcopy(m_np, m_MotionVars+Aoffset, 1, m_fA[i][0], 1);

            // Extrapolate to n+1
            Vmath::Smul(m_np, 
                        Betaq_Coeffs[nint-1][nint-1],
                        m_fV[i][nint-1],    1,
                        m_fV[i][nlevels-1], 1);
            Vmath::Smul(m_np, 
                        Betaq_Coeffs[nint-1][nint-1],
                        m_fA[i][nint-1],    1,
                        m_fA[i][nlevels-1], 1);

            for(int n = 0; n < nint-1; ++n)
            {
                Vmath::Svtvp(m_np, 
                             Betaq_Coeffs[nint-1][n],
                             m_fV[i][n],1,m_fV[i][nlevels-1],1,
                             m_fV[i][nlevels-1],1);
                Vmath::Svtvp(m_np, 
                             Betaq_Coeffs[nint-1][n],
                             m_fA[i][n],1,m_fA[i][nlevels-1],1,
                             m_fA[i][nlevels-1],1);
            }

            // Add the fictitious forces on the RHS of the equation
            Vmath::Svtvp(m_np, fictdamp,m_fV[i][nlevels-1],1,
                         fces[i],1,fces[i],1);
            Vmath::Svtvp(m_np, fictrho, m_fA[i][nlevels-1],1,
                         fces[i],1,fces[i],1);
        }
    }

    // Tensioned cable model is evaluated in wave space
    for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
    {
        int offset = cn*3*m_np;
        Array<OneD, NekDouble> tmp(3*m_np);
        TensionedCableModel(pfields, fces[cn], 
                            tmp = m_MotionVars+offset);
    }

    // Set the m_forcing term based on the motion of the cable
    for(int var = 0; var < 3; var++)
    {
        for(int plane = 0; plane < m_np; plane++)
        {
            int n = pfields[0]->GetPlane(plane)->GetTotPoints();

            Array<OneD, NekDouble> tmp;

            int offset  = plane * n;
            int xoffset = var * m_np+plane;
            int yoffset = 3*m_np + xoffset;

            Vmath::Fill(n, m_MotionVars[xoffset], tmp = m_zta[var] + offset, 1);
            Vmath::Fill(n, m_MotionVars[yoffset], tmp = m_eta[var] + offset, 1);
        }
    }
}

/**
 *
 */
void ForcingMovingBody::MappingBndConditions(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields,
            const Array<OneD, Array<OneD, NekDouble> >        &fields,
                  NekDouble  time) 
{
    Array<OneD, Array<OneD, NekDouble> > BndV (
                                        m_motion.num_elements());
    for (int i = 0; i < m_motion.num_elements(); ++i)
    {
        BndV[i] = Array<OneD, NekDouble>(m_np, 0.0);
    }

    for(int i = 0; i < m_motion.num_elements(); ++i)
    {
        int offset = i*3*m_np+m_np;
        Vmath::Vcopy(m_np, m_MotionVars+offset, 1, BndV[i], 1);
    }

    Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExps;
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds;

    const Array<OneD, const NekDouble> z
                            = pfields[0]->GetHomogeneousBasis()->GetZ();
    int nbnd = pfields[0]->GetBndCondExpansions().num_elements();

    for (int i = 0; i < nbnd; ++i)
    {
        for (int plane = 0; plane < m_np; plane++)
        {
            for ( int dim = 0; dim < m_motion.num_elements(); dim++)
            {
                bndCondExps = pfields[dim]->GetPlane(plane)
                                            ->GetBndCondExpansions();
                bndConds    = pfields[dim]->GetPlane(plane)
                                            ->GetBndConditions();

                if (bndConds[i]->GetUserDefined() ==
                                            SpatialDomains::eMovingBody)
                {
                    int npoints = bndCondExps[i]->GetNpoints();
                    Array<OneD, NekDouble> x0(npoints, 0.0);
                    Array<OneD, NekDouble> x1(npoints, 0.0);
                    Array<OneD, NekDouble> x2(npoints, 0.0);
                    Array<OneD, NekDouble> tmp(npoints,0.0);

                    NekDouble local_z = z[pfields[0]->GetTransposition()
                                                    ->GetPlaneID(plane)];
                    NekDouble x2_in   = 0.5*m_lhom*(1.0+local_z);
                    // Homogeneous input case for x2.
                    if (x2_in == NekConstants::kNekUnsetDouble)
                    {
                        bndCondExps[i]->GetCoords(x0,x1,x2);
                    }
                    else
                    {
                        bndCondExps[i]->GetCoords(x0, x1, x2);
                        Vmath::Fill(npoints, x2_in, x2, 1);
                    }

                    LibUtilities::Equation condition =
                        boost::static_pointer_cast<
                            SpatialDomains::DirichletBoundaryCondition>
                                (bndConds[i])->
                                    m_dirichletCondition;

                    condition.Evaluate(x0, x1, x2, time,
                                    bndCondExps[i]->UpdatePhys());
                    Vmath::Fill(npoints, BndV[dim][plane], tmp, 1);
                    Vmath::Vsub(npoints, bndCondExps[i]->UpdatePhys(), 1,
                                         tmp,                          1,
                                         bndCondExps[i]->UpdatePhys(), 1);
                }
            }
        }
    }
    
    // Set up high order outflow conditions 
    // for movingbody module
    CalcOutflowBCs(pfields, fields);
}


/**
 *
 **/
void ForcingMovingBody::CalcOutflowBCs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields,
        const Array<OneD, Array<OneD, NekDouble> >        &fields)
    {
        static bool init = true;
        static bool noHOBC = false;

        if(noHOBC == true)
        {
           return;
        }
        
        if(init) // set up storage for boundary velocity at outflow
        {
            init = false;

            m_PBndConds = pfields[3]->GetBndConditions();
            m_PBndExp   = pfields[3]->GetBndCondExpansions();

           int totbndpts = 0;
            for(int n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(m_PBndConds[n]->GetUserDefined()
                        == SpatialDomains::eHighOutflow)
                {
                    totbndpts += m_PBndExp[n]->GetTotPoints();
                }
            }
            
            if(totbndpts == 0)
            {
                noHOBC = true;
                return;
            }
            
            m_curl_dim = 3;
            m_bnd_dim  = 2;
            m_intSteps = 2;

            // find the maximum values of points for pressure BC evaluation
            m_pressureBCsMaxPts = 0;
            m_pressureBCsElmtMaxPts = 0;
            // Set up mapping from pressure boundary condition to pressure element
            // dm_etails.
            pfields[3]->GetBoundaryToElmtMap(m_pressureBCtoElmtID,
                                             m_pressureBCtoTraceID);
            int cnt, n;
            for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i)
                {
                    m_pressureBCsMaxPts = max(m_pressureBCsMaxPts,
                                    m_PBndExp[n]->GetExp(i)->GetTotPoints());
                    m_pressureBCsElmtMaxPts = max(m_pressureBCsElmtMaxPts,
                                    pfields[0]->GetExp(m_pressureBCtoElmtID[cnt++])
                                                                ->GetTotPoints());
                }
            }
 
            m_outflowVel = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (m_curl_dim);
            m_PhyoutfVel = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (m_curl_dim);
            m_movbodyVel = Array<OneD, Array<OneD, NekDouble> >(m_curl_dim);
            m_movbodyDev = Array<OneD, Array<OneD, NekDouble> >(m_curl_dim);
            m_ub         = Array<OneD, Array<OneD, NekDouble> >(m_curl_dim);
            m_nonlinearterm_phys   = Array<OneD, NekDouble> (totbndpts,0.0);
            for(int i = 0; i < m_curl_dim; ++i)
            {
                m_outflowVel[i] = Array<OneD, Array<OneD, NekDouble> >(m_curl_dim);
                m_PhyoutfVel[i] = Array<OneD, Array<OneD, NekDouble> >(m_curl_dim);
                for(int j = 0; j < m_curl_dim; ++j)
                {             
                    // currently just set up for 2nd order extrapolation 
                    m_outflowVel[i][j] = Array<OneD, NekDouble>(totbndpts,0.0);
                    m_PhyoutfVel[i][j] = Array<OneD, NekDouble>(totbndpts,0.0);
                }
                m_movbodyVel[i] = Array<OneD, NekDouble>(totbndpts,0.0);
                m_movbodyDev[i] = Array<OneD, NekDouble>(totbndpts,0.0);
                m_ub[i]         = Array<OneD, NekDouble>(totbndpts,0.0);
            }

            Array<OneD, unsigned int> planes;
            planes = pfields[0]->GetZIDs();
            int num_planes = planes.num_elements();
            m_expsize_per_plane = 
					Array<OneD, unsigned int> (m_PBndConds.num_elements());
            for(int n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                int exp_size = m_PBndExp[n]->GetExpSize();
                m_expsize_per_plane[n] = exp_size/num_planes;
            }
            m_totexps_per_plane = 0;
            for(int n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                m_totexps_per_plane += m_PBndExp[n]->GetExpSize()/num_planes;
            }
        }
        
        StdRegions::StdExpansionSharedPtr Bc; 
        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >
                                                        UBndConds(m_curl_dim);
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                                                        UBndExp(m_curl_dim);

        for (int i = 0; i < m_curl_dim; ++i)
        {
            UBndConds[i] = pfields[i]->GetBndConditions();
            UBndExp[i]   = pfields[i]->GetBndCondExpansions();
        }

        Array<OneD, Array<OneD, NekDouble> > BndValues(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > MovValues(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > BndElmt  (m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > ngradu   (m_curl_dim);
        Array<OneD, NekDouble > gradtmp (m_pressureBCsElmtMaxPts),
                                fgradtmp(m_pressureBCsElmtMaxPts);
        
        ngradu[0] = Array<OneD, NekDouble>(m_curl_dim*m_pressureBCsMaxPts);
        for(int i = 0; i < m_curl_dim; ++i)
        {
            BndElmt[i]   = Array<OneD, NekDouble> (m_pressureBCsElmtMaxPts,
                                                   0.0);
            BndValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
            MovValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
            ngradu[i]  = ngradu[0] + i*m_pressureBCsMaxPts;
            RollOver(m_outflowVel[i]);
            RollOver(m_PhyoutfVel[i]);
        }
            
        m_pressureCalls++;
        int  nint = min(m_pressureCalls,m_intSteps);

        StdRegions::StdExpansionSharedPtr elmt;
        Array<OneD, NekDouble> PBCvals, UBCvals;
        Array<OneD, Array<OneD, NekDouble> > normals;

        NekDouble U0,delta;
        m_session->LoadParameter("U0_HighOrderBC",U0,1.0);
        m_session->LoadParameter("Delta_HighOrderBC",delta,1/20.0);

        int cnt;
        int count = 0;
        for(int n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // Do outflow boundary conditions if they exist
            if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHighOutflow)
            {

                int cnt_exp   = 0; int cnt_plane = 0;
                int veloffset = 0;
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i, cnt_exp++)
                {
                    // count the expansion list in each plane for e3DH1D case
                    if(cnt_exp == m_expsize_per_plane[n])
                    {
                        cnt_exp = 0; cnt_plane++;
                    }
                    int cnt = cnt_plane * m_totexps_per_plane + cnt_exp + count;

                    // find element and edge of this expansion.
                    Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion>
                        (m_PBndExp[n]->GetExp(i));

                    int elmtid = m_pressureBCtoElmtID[cnt];
                    elmt       = pfields[0]->GetExp(elmtid);
                    int offset = pfields[0]->GetPhys_Offset(elmtid);

                    int boundary = m_pressureBCtoTraceID[cnt];

                    // Determine extrapolated U,V values
                    int nq = elmt->GetTotPoints();
                    int nbc = m_PBndExp[n]->GetExp(i)->GetTotPoints();
                    // currently just using first order approximation here.
                    // previously have obtained value from m_integrationSoln
                    Array<OneD, NekDouble> veltmp;

                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        Vmath::Vcopy(nq, &fields[j][offset], 1,
                                         &BndElmt[j][0],     1);
                        elmt->GetTracePhysVals(boundary,Bc,BndElmt[j],
                                     veltmp = m_outflowVel[j][0] + veloffset);
                    }

                    for(int j = 0; j < m_bnd_dim; ++j)
                    {
                        if(j == 0)
                        {
                            Vmath::Vcopy(nq, &m_zta[1][offset], 1, &BndElmt[0][0], 1);
                        }
                        if(j == 1)
                        {
                            Vmath::Vcopy(nq, &m_eta[1][offset], 1, &BndElmt[1][0], 1);
                        }
                        // extract the values of velocity of movingbody on outflow boundaries
                        elmt->GetTracePhysVals(boundary,Bc,BndElmt[j],
                                         veltmp = m_movbodyVel[j] + veloffset);
                        if(j == 0)
                        {
                            Vmath::Vcopy(nq, &m_zta[3][offset], 1, &BndElmt[0][0], 1);
                        }
                        if(j == 1)
                        {
                            Vmath::Vcopy(nq, &m_eta[3][offset], 1, &BndElmt[1][0], 1);
                        }
                        // extract the values of d(m_zta)/dz and d(m_eta)/dz on outflow boundaries
                        elmt->GetTracePhysVals(boundary,Bc,BndElmt[j],
                                         veltmp = m_movbodyDev[j] + veloffset);
                    }

                    veloffset += nbc;
                }

                // for velocity on the outflow boundary in e3DH1D,
                // we need to make a backward fourier transformation
                // to get the physical coeffs at the outflow BCs.
                for(int j = 0; j < m_curl_dim; ++j)
                {
                    m_PBndExp[n]->HomogeneousBwdTrans(
                            m_outflowVel[j][0],
                            m_PhyoutfVel[j][0]);
                }
                
                cnt_plane = 0; cnt_exp = 0;
                veloffset = 0;
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i, cnt_exp++)
                {   
                    // count the expansion list for each plane for e3DH1D
                    if(cnt_exp == m_expsize_per_plane[n])
                    {
                        cnt_exp = 0; cnt_plane++;
                    }
                    cnt = cnt_plane * m_totexps_per_plane + cnt_exp + count;

                    int elmtid = m_pressureBCtoElmtID[cnt];
                    elmt       = pfields[0]->GetExp(elmtid);
                    int nbc = m_PBndExp[n]->GetExp(i)->GetTotPoints();

                    Array<OneD, NekDouble>  veltmp(nbc,0.0),
                                            normDotu(nbc,0.0), utot(nbc,0.0);
                    int boundary = m_pressureBCtoTraceID[cnt];
                    normals=elmt->GetSurfaceNormal(boundary);

                    // extrapolate velocity
                    if(nint <= 1)
                    {
                        for(int j = 0; j < m_curl_dim; ++j)
                        {
                            Vmath::Vcopy(nbc,
                                    veltmp = m_PhyoutfVel[j][0] +veloffset, 1,
                                    BndValues[j],                           1);
                        }
                    }
                    else // only set up for 2nd order extrapolation
                    {
                        for(int j = 0; j < m_curl_dim; ++j)
                        {
                            Vmath::Smul(nbc, 2.0,
                                    veltmp = m_PhyoutfVel[j][0] + veloffset, 1,
                                    BndValues[j],                            1);
                            Vmath::Svtvp(nbc, -1.0,
                                    veltmp = m_PhyoutfVel[j][1] + veloffset, 1,
                                    BndValues[j],                            1,
                                    BndValues[j],                            1);
                        }    
                    }

                    // Set up |u|^2, n.u in physical space
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        Vmath::Vvtvp(nbc, BndValues[j], 1, BndValues[j], 1,
                                          utot,         1, utot,         1);
                    }
                    for(int j = 0; j < m_bnd_dim; ++j)
                    {
                        Vmath::Vvtvp(nbc, normals[j],   1, BndValues[j], 1,
                                          normDotu,     1, normDotu,     1);
                    }
                        
                    int Offset = m_PBndExp[n]->GetPhys_Offset(i);

                    for(int k = 0; k < nbc; ++k)
                    {
                        // calculate the nonlinear term (kinetic energy
                        //  multiplies step function) in physical space
                        NekDouble fac = 0.5*(1.0-tanh(normDotu[k]/(U0*delta)));
                        m_nonlinearterm_phys[k + Offset] =  0.5 * utot[k] * fac;
                    }


                    // Include the constribution of movingbody effects
                    for(int j = 0; j < m_bnd_dim; ++j)
                    {
                        Vmath::Vvtvp(nbc,veltmp = m_PhyoutfVel[2][0] + veloffset, 1,
                                         veltmp = m_movbodyDev[j]    + veloffset, 1,
                                         veltmp = m_movbodyVel[j]    + veloffset, 1,
                                         MovValues[j],                            1);
                        Vmath::Vadd(nbc, MovValues[j], 1,
                                         BndValues[j], 1,
                                         BndValues[j], 1);
                    }


                    Vmath::Zero(nbc,utot,1);
                    Vmath::Zero(nbc,normDotu,1);

                    // Set up |u|^2, n.u in physical space
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        Vmath::Vvtvp(nbc, BndValues[j], 1, BndValues[j], 1,
                                          utot,         1, utot,         1);
                    }
                    for(int j = 0; j < m_bnd_dim; ++j)
                    {
                        Vmath::Vvtvp(nbc, normals[j],   1, BndValues[j], 1,
                                          normDotu,     1, normDotu,     1);
                    }

                    for(int k = 0; k < nbc; ++k)
                    {
                        // calculate the nonlinear term (kinetic energy
                        //  multiplies step function) in physical space
                        NekDouble fac = 0.5*(1.0-tanh(normDotu[k]/(U0*delta)));
                        m_nonlinearterm_phys[k + Offset] += -0.5 * utot[k] * fac;
                    }
                    veloffset += nbc;
                }

                // m_forcing size (it must be 3: Ax, Ay and Az)
                // total number of structural variables(Disp, Vel and Accel) must be 3
                Array<OneD, Array<OneD, NekDouble> > Fields(3);
                // create the storage space for the body motion description
                int phystot = pfields[0]->GetTotPoints();
                for(int i = 0; i < 3; i++)
                {
                    Fields[i] = Array<OneD, NekDouble> (phystot, 0.0);
                }
                for(int j = 0; j < m_curl_dim; ++j)
                {
                    pfields[0]->HomogeneousBwdTrans(fields[j],Fields[j]);
                }
    
                veloffset = 0;
                cnt_exp = 0; cnt_plane = 0; //only useful in e3DH1D case
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i)
                {
                    // count the expansion list for e3DH1D
                    if(cnt_exp == m_expsize_per_plane[n])
                    {
                        cnt_exp = 0; cnt_plane++;
                    }
                    cnt = cnt_plane * m_totexps_per_plane + cnt_exp + count;
                    cnt_exp++;

                    // find element and edge of this expansion. 
                    Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                        (m_PBndExp[n]->GetExp(i));

                    int elmtid = m_pressureBCtoElmtID[cnt];
                    elmt       = pfields[0]->GetExp(elmtid);
                    int offset = pfields[0]->GetPhys_Offset(elmtid);

                    // Determine extrapolated U,V values
                    int nq = elmt->GetTotPoints();

                    // currently just using first order approximation here. 
                    // previously have obtained value from m_integrationSoln
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        Vmath::Vcopy(nq, &Fields[j][offset], 1,
                                         &BndElmt[j][0],     1);
                    }

                    int nbc      = m_PBndExp[n]->GetExp(i)->GetTotPoints();
                    int boundary = m_pressureBCtoTraceID[cnt];

                    Array<OneD, NekDouble>  ptmp(nbc,0.0),devtmp(nbc,0.0),
                                            divU(nbc,0.0);

                    normals=elmt->GetSurfaceNormal(boundary);
                    Vmath::Zero(m_bnd_dim*m_pressureBCsMaxPts,ngradu[0],1);

                    for (int j = 0; j < m_bnd_dim; j++)
                    {
                        // Calculate gradient terms:
                        // nx(w_x*d(m_zta)/dz)+ny(w_y*d(m_zta)/dz),
                        // nx(w_x*d(m_eta)/dz)+ny(w_y*d(m_eta)/dz)
                        for (int k = 0; k< m_bnd_dim; k++)
                        {
                            elmt->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                            BndElmt[2], gradtmp);
                            elmt->GetTracePhysVals(boundary, Bc, gradtmp,
                                                   fgradtmp);
                            Vmath::Vmul(nbc, devtmp = m_movbodyDev[j] + veloffset, 1,
                                             fgradtmp, 1, fgradtmp, 1);
                            Vmath::Vvtvp(nbc,normals[k], 1, fgradtmp, 1,
                                             ngradu[j],  1, ngradu[j],1);
                        }
                    }

                    // Set up (n.grad(u) . n) for
                    // pressure condition
                    for(int j = 0; j < m_bnd_dim; ++j)
                    {
                        Vmath::Vvtvp(nbc, normals[j],   1, ngradu[j],    1,
                                          ptmp,         1, ptmp,         1);
                    }
                    int p_offset = m_PBndExp[n]->GetPhys_Offset(i);

                    for(int k = 0; k < nbc; ++k)
                    {
                        // Set up Dirichlet pressure condition and
                        // store in ptmp (m_UBndCoeffs contains Fourier Coeffs of the
                        // function from the input file )

                        ptmp[k] = - m_kinvis * ptmp[k] 
                                  - m_nonlinearterm_phys[k + p_offset];
                    }

                    // set up pressure boundary condition
                    PBCvals = m_PBndExp[n]->UpdatePhys()
                            + m_PBndExp[n]->GetPhys_Offset(i);
                    Vmath::Vcopy(nbc, ptmp, 1, PBCvals, 1);
                    for(int j = 0; j < m_curl_dim; j++)
                    {
                        int u_offset = UBndExp[j][n]->GetPhys_Offset(i);

                        for(int k = 0; k < nbc; ++k)
                        {
                            if(j < m_curl_dim-1)
                            {
                                m_ub[j][k + u_offset] = - m_nonlinearterm_phys[k + u_offset] 
                                                                  * normals[j][k] / m_kinvis
                                                                              - ngradu[j][k];
                            }
                            else
                            {
                                m_ub[j][k + u_offset] = - ngradu[j][k];
                            }
                        }
                    }

                    veloffset += nbc;
                }

                //
                for(int j = 0; j < m_curl_dim; ++j)
                {
                    m_PBndExp[n]->HomogeneousFwdTrans(
                            m_ub[j],
                            m_ub[j]);
                }

                // Now set up Velocity conditions.
                for(int j = 0; j < m_curl_dim; j++)
                {
                    if(UBndConds[j][n]->GetUserDefined()
                                        == SpatialDomains::eHighOutflow)
                    {
                        for(int i = 0; i < UBndExp[0][n]->GetExpSize();
                                       ++i)
                        {
                            int u_offset = UBndExp[j][n]->GetPhys_Offset(i);
                            UBCvals = UBndExp[j][n]->UpdateCoeffs()
                                    + UBndExp[j][n]->GetCoeff_Offset(i);
                            Bc->IProductWRTBase(m_ub[j]+u_offset,UBCvals);
                        }
                    }
                }
            }
            else
            {
                count  += m_expsize_per_plane[n];
            }
        }
    }
    
/**
 *
 */
void ForcingMovingBody::InitialiseCableModel(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields)
{
    // storage of spanwise-velocity for current and previous time levels
    // used to calculate dw/dt in m_forcing term
    int nPointsTot = pfields[0]->GetNpoints();
    m_W = Array<OneD, Array<OneD, NekDouble> > (2);
    for(int n = 0; n < 2; ++n)
    {
        m_W[n] = Array<OneD, NekDouble> (nPointsTot, 0.0);
    }

    m_session->LoadParameter("Kinvis",m_kinvis);
    m_session->LoadParameter("TimeStep", m_timestep, 0.01);

    m_movingBodyCalls = 0;

    Array<OneD, unsigned int> ZIDs;
    ZIDs = pfields[0]->GetZIDs();
    m_np = ZIDs.num_elements();
    
    m_Aeroforces = Array<OneD, NekDouble>(2*m_np,0.0);
    m_MotionVars = Array<OneD, NekDouble>(6*m_np,0.0);

    std::string vibtype = m_session->GetSolverInfo("VibrationType");

    if(vibtype == "Constrained" || vibtype == "CONSTRAINED")
    {
        m_vdim = 1;
    }
    else if (vibtype == "Free" || vibtype == "FREE")
    {
        m_vdim = 2;
    }
    else if (vibtype == "Forced" || vibtype == "FORCED")
    {
        return;
    }

    LibUtilities::CommSharedPtr vcomm = pfields[0]->GetComm();

    bool homostrip;
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);

    if(!homostrip)
    {
        m_session->LoadParameter("LZ", m_lhom);
        int nplanes = m_session->GetParameter("HomModesZ");
        m_FFT = 
            LibUtilities::GetNektarFFTFactory().CreateInstance(
                                            "NekFFTW", nplanes);
    }
    else
    {
        int nstrips;

        NekDouble DistStrip;

        m_session->LoadParameter("DistStrip", DistStrip);
        m_session->LoadParameter("Strip_Z", nstrips);
        m_lhom = nstrips * DistStrip;
        m_FFT = 
            LibUtilities::GetNektarFFTFactory().CreateInstance(
                                            "NekFFTW", nstrips);
    }

    // load the structural dynamic parameters from xml file
    m_session->LoadParameter("StructRho",  m_structrho);
    m_session->LoadParameter("StructDamp", m_structdamp, 0.0);

    // Identify whether the fictitious mass method is active for explicit
    // coupling of fluid solver and structural dynamics solver
    bool fictmass;
    m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                fictmass, false);
    if(fictmass)
    {
        NekDouble fictrho, fictdamp;
        m_session->LoadParameter("FictMass", fictrho);
        m_session->LoadParameter("FictDamp", fictdamp);
        m_structrho  += fictrho;
        m_structdamp += fictdamp;

        // Storage array of Struct Velocity and Acceleration used for
        // extrapolation of fictitious force
        m_fV = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                            m_motion.num_elements());
        m_fA = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                            m_motion.num_elements());
        for (int i = 0; i < m_motion.num_elements(); ++i)
        {
            m_fV[i] = Array<OneD, Array<OneD, NekDouble> > (2);
            m_fA[i] = Array<OneD, Array<OneD, NekDouble> > (2);

            for(int n = 0; n < 2; ++n)
            {
                m_fV[i][n] = Array<OneD, NekDouble>(m_np, 0.0);
                m_fA[i][n] = Array<OneD, NekDouble>(m_np, 0.0);
            }
        }
    }

    // Setting the coefficient matrices for solving structural dynamic ODEs
    SetDynEqCoeffMatrix(pfields);

    // Set initial condition for cable's motion
    int cnt = 0;

    for(int j = 0; j < m_funcName.num_elements(); j++)
    {
                
        // loading from the specified files through inputstream
        if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
        {
    
            std::ifstream inputStream;

            int nzpoints;

            if(homostrip)
            {
                int nstrips;
                m_session->LoadParameter("Strip_Z", nstrips);
                nzpoints = nstrips;
            }
            else
            { 
                nzpoints = pfields[0]->GetHomogeneousBasis()->GetNumModes();
            }

 
            if (vcomm->GetRank() == 0)
            {

                Array<OneD, NekDouble> motion_x(3*nzpoints, 0.0);
                Array<OneD, NekDouble> motion_y(3*nzpoints, 0.0);

                Array<OneD, std::string> tmp(9);

                std::string filename = m_session->GetFunctionFilename(
                    m_funcName[j], m_motion[0]);
                
                // Open intputstream for cable motions
                inputStream.open(filename.c_str());

                // Import the head string from the file
                for(int n = 0; n < tmp.num_elements(); n++)
                {
                    inputStream >> tmp[n];
                }

                NekDouble time, z_cds;
                // Import the motion variables from the file
                for (int n = 0; n < nzpoints; n++)
                {
                    inputStream >> setprecision(6) >> time;
                    inputStream >> setprecision(6) >> z_cds;
                    inputStream >> setprecision(8) >> motion_x[n];
                    inputStream >> setprecision(8) >> motion_x[n+nzpoints];
                    inputStream >> setprecision(8) >> motion_x[n+2*nzpoints];
                    inputStream >> setprecision(8) >> motion_y[n];
                    inputStream >> setprecision(8) >> motion_y[n+nzpoints];
                    inputStream >> setprecision(8) >> motion_y[n+2*nzpoints];
                }

                // Close inputstream for cable motions
                inputStream.close();
                
                if(!homostrip)
                {
                    int nproc = vcomm->GetColumnComm()->GetSize();
                    for (int i = 1; i < nproc; ++i)
                    {
                        for(int plane = 0; plane < m_np; plane++)
                        {
                            for (int var = 0; var < 3; var++)
                            {
                                int xoffset = var*m_np+plane;
                                int yoffset = 3*m_np+xoffset;
                                m_MotionVars[xoffset] = 
                                            motion_x[plane+i*m_np+var*nzpoints];
                                m_MotionVars[yoffset] = 
                                            motion_y[plane+i*m_np+var*nzpoints];
                            }
                        }

                        vcomm->GetColumnComm()->Send(i, m_MotionVars);
                    }

                    for(int plane = 0; plane < m_np; plane++)
                    {
                        for (int var = 0; var < 3; var++)
                        {
                            int xoffset = var*m_np+plane;
                            int yoffset = 3*m_np+xoffset;
                            m_MotionVars[xoffset] = motion_x[plane+var*nzpoints];
                            m_MotionVars[yoffset] = motion_y[plane+var*nzpoints];
                        }
                    }
                }
                else //for strip model, make sure that each processor gets its motion values correctly
                {
                    int nstrips;
                    m_session->LoadParameter("Strip_Z", nstrips);

                    int nproc = vcomm->GetColumnComm()->GetSize();
                    
                    for (int i = 1; i < nstrips; ++i)
                    {
                        for(int plane = 0; plane < m_np; plane++)
                        {
                            for (int var = 0; var < 3; var++)
                            {
                                int xoffset = var*m_np+plane;
                                int yoffset = 3*m_np+xoffset;
                                m_MotionVars[xoffset] = motion_x[i+var*nstrips];
                                m_MotionVars[yoffset] = motion_y[i+var*nstrips];
                            }
                        }

                        for (int j = 0; j < nproc/nstrips; j++)
                        {
                            vcomm->GetColumnComm()->Send(i+j*nstrips, m_MotionVars);
                        }
                    }
                    
                    for(int plane = 0; plane < m_np; plane++)
                    {
                        for (int var = 0; var < 3; var++)
                        {
                            int xoffset = var*m_np+plane;
                            int yoffset = 3*m_np+xoffset;
                            m_MotionVars[xoffset] = motion_x[var*nstrips];
                            m_MotionVars[yoffset] = motion_y[var*nstrips];
                        }
                    }

                    for (int j = 1; j < nproc/nstrips; j++)
                    {
                        vcomm->GetColumnComm()->Send(j*nstrips, m_MotionVars);
                    }
                }
            }
            else
            {
                if(!homostrip)
                {
                    int colrank = vcomm->GetColumnComm()->GetRank();
                    int nproc   = vcomm->GetColumnComm()->GetSize();

                    for (int j = 1; j < nproc; j++)
                    {
                        if(colrank == j)
                        {
                            vcomm->GetColumnComm()->Recv(0, m_MotionVars);
                        }
                    }
                }
                else //for strip model
                {
                    int nstrips;
                    m_session->LoadParameter("Strip_Z", nstrips);

                    int colrank = vcomm->GetColumnComm()->GetRank();
                    int nproc   = vcomm->GetColumnComm()->GetSize();

                    for(int i = 1; i < nstrips; i++)
                    {
                        for (int j = 0; j < nproc/nstrips; j++)
                        {
                            if(colrank == i+j*nstrips)
                            {
                                vcomm->GetColumnComm()->Recv(0, m_MotionVars);
                            }
                        }
                    }

                    for (int j = 1; j < nproc/nstrips; j++)
                    {
                        if(colrank == j*nstrips)
                        {
                            vcomm->GetColumnComm()->Recv(0, m_MotionVars);
                        }
                    }
                }
            }

            cnt = cnt + 2;
        }
        else //Evaluate from the functions specified in xml file
        {
            Array<OneD, NekDouble> x0(m_np, 0.0);
            Array<OneD, NekDouble> x1(m_np, 0.0);
            Array<OneD, NekDouble> x2(m_np, 0.0);
            
            if(!homostrip)
            {
                int nzpoints = pfields[0]->GetHomogeneousBasis()->GetNumModes();
                Array<OneD, NekDouble> z_coords(nzpoints,0.0);
                Array<OneD, const NekDouble> pts
                                    = pfields[0]->GetHomogeneousBasis()->GetZ();

                Vmath::Smul(nzpoints,m_lhom/2.0,pts,1,z_coords,1);
                Vmath::Sadd(nzpoints,m_lhom/2.0,z_coords,1,z_coords,1);

                for (int plane = 0; plane < m_np; plane++)
                {
                    x2[plane] = z_coords[ZIDs[plane]];
                }
            }
            else
            {   
                int nstrips;
                m_session->LoadParameter("Strip_Z", nstrips);

                ASSERTL0(m_session->DefinesSolverInfo("USEFFT"),
                            "Fourier transformation of cable motion is currently "
                            "implemented only for FFTW module.");
                
                NekDouble DistStrip;
                m_session->LoadParameter("DistStrip", DistStrip);

                int colrank = vcomm->GetColumnComm()->GetRank();
                int nproc   = vcomm->GetColumnComm()->GetSize();
                
                for(int i = 0; i < nstrips; ++i)
                {
                    for(int j = 0; j < nproc/nstrips; j++)
                    {
                        if (colrank == i+j*nstrips)
                        {
                            for (int plane = 0; plane < m_np; plane++)
                            {
                                x2[plane] = i*DistStrip;
                            }
                        }
                    }
                }
            }

            int xoffset = j*m_np;
            int yoffset = 3*m_np+xoffset;

            Array<OneD, NekDouble> tmp(m_np,0.0);
            LibUtilities::EquationSharedPtr ffunc0,ffunc1;

            ffunc0 = m_session->GetFunction(m_funcName[j], m_motion[0]);
            ffunc1 = m_session->GetFunction(m_funcName[j], m_motion[1]);

            ffunc0->Evaluate(x0, x1, x2, 0.0, tmp = m_MotionVars+xoffset);
            ffunc1->Evaluate(x0, x1, x2, 0.0, tmp = m_MotionVars+yoffset);
            cnt = cnt + 2;
        }
    }
}


/**
 *
 */
void ForcingMovingBody::SetDynEqCoeffMatrix(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields)
{
    int nplanes;

    bool homostrip;
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);

    if(!homostrip)
    {
        nplanes = m_session->GetParameter("HomModesZ"); 
    }
    else
    {
        m_session->LoadParameter("Strip_Z", nplanes);
    }

    m_CoeffMat_A = Array<OneD, DNekMatSharedPtr>(nplanes);
    m_CoeffMat_B = Array<OneD, DNekMatSharedPtr>(nplanes);

    NekDouble tmp1, tmp2, tmp3;
    NekDouble tmp4, tmp5, tmp6, tmp7;

    // load the structural dynamic parameters from xml file
    NekDouble cabletension;
    NekDouble bendingstiff;
    NekDouble structstiff;
    m_session->LoadParameter("StructStiff",  structstiff,  0.0);
    m_session->LoadParameter("CableTension", cabletension, 0.0);
    m_session->LoadParameter("BendingStiff", bendingstiff, 0.0);

    tmp1 =   m_timestep * m_timestep;
    tmp2 =  structstiff / m_structrho;
    tmp3 = m_structdamp / m_structrho;
    tmp4 = cabletension / m_structrho;
    tmp5 = bendingstiff / m_structrho;

    // solve the ODE in the wave space for cable motion to obtain disp, vel and
    // accel

    std::string supptype = m_session->GetSolverInfo("SupportType");

    for(int plane = 0; plane < nplanes; plane++)
    {
        int nel = 3;
        m_CoeffMat_A[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);
        m_CoeffMat_B[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);

        unsigned int K;
        NekDouble beta;

        if(supptype == "FREE-FREE" || 
                    supptype == "Free-Free")
        {
            K = plane/2;
            beta = 2.0 * M_PI/m_lhom;
        }
        else if(supptype == "PINNED-PINNED" || 
                        supptype == "Pinned-Pinned")
        {   
            K = plane+1;
            beta = M_PI/m_lhom;
        }
        else
        {
            ASSERTL0(false,
                        "Unrecognized support type for cable's motion");
        }

        tmp6 = beta * K;
        tmp6 = tmp6 * tmp6;
        tmp7 = tmp6 * tmp6;
        tmp7 = tmp2 + tmp4 * tmp6 + tmp5 * tmp7;

        (*m_CoeffMat_A[plane])(0,0) = tmp7;
        (*m_CoeffMat_A[plane])(0,1) = tmp3;
        (*m_CoeffMat_A[plane])(0,2) = 1.0;
        (*m_CoeffMat_A[plane])(1,0) = 0.0;
        (*m_CoeffMat_A[plane])(1,1) = 1.0;
        (*m_CoeffMat_A[plane])(1,2) =-m_timestep/2.0;
        (*m_CoeffMat_A[plane])(2,0) = 1.0;
        (*m_CoeffMat_A[plane])(2,1) = 0.0;
        (*m_CoeffMat_A[plane])(2,2) =-tmp1/4.0;
        (*m_CoeffMat_B[plane])(0,0) = 0.0;
        (*m_CoeffMat_B[plane])(0,1) = 0.0;
        (*m_CoeffMat_B[plane])(0,2) = 0.0;
        (*m_CoeffMat_B[plane])(1,0) = 0.0;
        (*m_CoeffMat_B[plane])(1,1) = 1.0;
        (*m_CoeffMat_B[plane])(1,2) = m_timestep/2.0;
        (*m_CoeffMat_B[plane])(2,0) = 1.0;
        (*m_CoeffMat_B[plane])(2,1) = m_timestep;
        (*m_CoeffMat_B[plane])(2,2) = tmp1/4.0;

        m_CoeffMat_A[plane]->Invert();
        (*m_CoeffMat_B[plane]) =
            (*m_CoeffMat_A[plane]) * (*m_CoeffMat_B[plane]);
    }
}


/**
 * Function to roll time-level storages to the next step layout.
 * The stored data associated with the oldest time-level
 * (not required anymore) are moved to the top, where they will
 * be overwritten as the solution process progresses.
 */
void ForcingMovingBody::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
{
    int nlevels = input.num_elements();

    Array<OneD, NekDouble> tmp;
    tmp = input[nlevels-1];

    for(int n = nlevels-1; n > 0; --n)
    {
        input[n] = input[n-1];
    }

    input[0] = tmp;
}


/**
 *
 */
void ForcingMovingBody::EvaluateAccelaration(
        const Array<OneD, NekDouble> &input,
              Array<OneD, NekDouble> &output,
              int npoints)
{
    m_movingBodyCalls++;

    // Rotate acceleration term
    RollOver(m_W);

    Vmath::Vcopy(npoints, input, 1, m_W[0], 1);

    //Calculate acceleration term at level n based on previous steps
    if (m_movingBodyCalls > 1)
    {
        Vmath::Smul(npoints,
                    1.0,
                    m_W[0], 1,
                    output, 1);
        Vmath::Svtvp(npoints,
                    -1.0,
                    m_W[1], 1,
                    output,   1,
                    output,   1);
        Vmath::Smul(npoints,
                    1.0/m_timestep,
                    output, 1,
                    output, 1);
    }
}


/**
 *
 */
void ForcingMovingBody::CheckIsFromFile(const TiXmlElement* pforce)
{

    m_funcName = Array<OneD, std::string> (3);
    m_motion = Array<OneD, std::string> (2);
    m_motion[0] = "x";
    m_motion[1] = "y";

    m_IsFromFile = Array<OneD, bool> (6);
    // Loading the x-dispalcement (m_zta) and the y-displacement (m_eta)
    // Those two variables are bith functions of z and t and the may come
    // from an equation (forced vibration) or from another solver which, given
    // the aerodynamic forces at the previous step, calculates the 
    // displacements.

    //Get the body displacement: m_zta and m_eta
    const TiXmlElement* funcNameElmt_D
                    = pforce->FirstChildElement("DISPLACEMENTS");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body displacement as d(z,t).");

    m_funcName[0] = funcNameElmt_D->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[0]),
             "Function '" + m_funcName[0] + "' not defined.");

    //Get the body velocity of movement: d(m_zta)/dt and d(m_eta)/dt
    const TiXmlElement* funcNameElmt_V
                    = pforce->FirstChildElement("VELOCITIES");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body velocity of movement as v(z,t).");

    m_funcName[1] = funcNameElmt_V->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[1]),
             "Function '" + m_funcName[1] + "' not defined.");


    //Get the body acceleration: dd(m_zta)/ddt and dd(m_eta)/ddt
    const TiXmlElement* funcNameElmt_A
                    = pforce->FirstChildElement("ACCELERATIONS");
    ASSERTL0(funcNameElmt_A,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body acceleration as a(z,t).");

    m_funcName[2] = funcNameElmt_A->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[2]),
             "Function '" + m_funcName[2] + "' not defined.");

    LibUtilities::FunctionType vType;

    // Check Displacement x
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[0]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[0] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[0] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Displacement y
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[1]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[1] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[1] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in y must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Velocity x
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[2] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[2] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in x must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Velocity y
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[3] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[3] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in y must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Acceleration x
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[4] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[4] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Acceleration y
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[5] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[5] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in y must be specified via an "
                        "equation <E> or a file <F>");
    }
}

/**
 *
 */
void ForcingMovingBody::InitialiseFilter(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields,
        const TiXmlElement* pforce)
{
    // Get the outputfile name, output frequency and 
    // the boundary's ID for the cable's wall
    std::string typeStr = pforce->Attribute("TYPE");
    std::map<std::string, std::string> vParams;

    const TiXmlElement *param = pforce->FirstChildElement("PARAM");
    while (param)
    {
        ASSERTL0(param->Attribute("NAME"),
                 "Missing attribute 'NAME' for parameter in filter "
                 + typeStr + "'.");
        std::string nameStr = param->Attribute("NAME");

        ASSERTL0(param->GetText(), "Empty value string for param.");
        std::string valueStr = param->GetText();

        vParams[nameStr] = valueStr;

        param = param->NextSiblingElement("PARAM");
    }

    // Creat the filter for MovingBody, where we performed the calculation of
    // fluid forces and write both the aerodynamic forces and motion variables
    // into the output files
    m_MovBodyfilter = MemoryManager<FilterMovingBody>::
                                    AllocateSharedPtr(pSession, vParams);

    // Initialise the object of MovingBody filter
    m_MovBodyfilter->Initialise(pfields, 0.0);

}

}
