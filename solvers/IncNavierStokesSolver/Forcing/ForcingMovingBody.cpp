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
// Description: Moving Body forcing (movement of a body in a domain is achieved
// via a forcing term which is the results of a coordinate system motion)
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingMovingBody.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{

NekDouble ForcingMovingBody::StifflyStable_Betaq_Coeffs[3][3] = {
    { 1.0,  0.0, 0.0},{ 2.0, -1.0, 0.0},{ 3.0, -3.0, 1.0}};
NekDouble ForcingMovingBody::StifflyStable_Alpha_Coeffs[3][3] = {
    { 1.0,  0.0, 0.0},{ 2.0, -0.5, 0.0},{ 3.0, -1.5, 1.0/3.0}};
NekDouble ForcingMovingBody::StifflyStable_Gamma0_Coeffs[3] = {
      1.0,  1.5, 11.0/6.0};

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
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
    // Just 3D homogenous 1D problems can use this techinque
    ASSERTL0(pFields[0]->GetExpType()==MultiRegions::e3DH1D,
             "Moving body implemented just for 3D Homogenous 1D expansions.");

    // only a first order approximation is considered for dw/dt in forcing term
    m_intSteps        = 1;
    // forcing size (it must be 3: Ax, Ay and Az)
    // total number of structural variables(Disp, Vel and Accel) must be 3
    m_Forcing         = Array<OneD, Array<OneD, NekDouble> > (3);
    m_funcName        = Array<OneD, std::string> (3);
    m_motion          = Array<OneD, std::string> (2);
    m_motion[0]       = "x";
    m_motion[1]       = "y";
    m_movingBodyCalls =  0 ;
    m_IsFromFile      = Array<OneD, bool> (6);

    m_session->LoadParameter("Kinvis",m_kinvis);
    m_session->LoadParameter("TimeStep", m_timestep, 0.01);

    // storage of spanwise-velocity for current and previous time levels
    // used to calculate dw/dt in forcing term
    int nPointsTot = pFields[0]->GetNpoints();
    m_W    = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);
    m_W[0] = Array<OneD, NekDouble> (nPointsTot, 0.0);
    for(int n = 0; n < m_intSteps; ++n)
    {
        m_W[n+1] = Array<OneD, NekDouble> (nPointsTot, 0.0);
    }

    // Loading the x-dispalcement (m_zeta) and the y-displacement (m_eta)
    // Those two variables are bith functions of z and t and the may come
    // from an equation (forced vibration) or from another solver which, given
    // the aerodynamic forces at the previous step, calculates the 
    // displacements.

    //Get the body displacement: m_zeta and m_eta
    const TiXmlElement* funcNameElmt_D 
                    = pForce->FirstChildElement("DISPLACEMENTS");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body displacement as d(z,t).");

    m_funcName[0] = funcNameElmt_D->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[0]),
             "Function '" + m_funcName[0] + "' not defined.");

    //Get the body velocity of movement: d(m_zeta)/dt and d(m_eta)/dt
    const TiXmlElement* funcNameElmt_V 
                    = pForce->FirstChildElement("VELOCITIES");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body velocity of movement as v(z,t).");

    m_funcName[1] = funcNameElmt_V->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[1]),
             "Function '" + m_funcName[1] + "' not defined.");


    //Get the body acceleration: dd(m_zeta)/ddt and dd(m_eta)/ddt
    const TiXmlElement* funcNameElmt_A 
                    = pForce->FirstChildElement("ACCELERATIONS");
    ASSERTL0(funcNameElmt_A,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body acceleration as a(z,t).");

    m_funcName[2] = funcNameElmt_A->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[2]),
             "Function '" + m_funcName[2] + "' not defined.");

    // At this point we know in the xml file where those quantities
    // are declared (equation or file) - via a function name which is now
    // stored in funcNameD etc. We need now to fill in with this info the
    // m_zeta and m_eta vectors (actuallythey are matrices) Array to control if
    // the motion is determined by an equation or is from a file.(not Nektar++)
    // check if we need to load a file or we have an equation
    CheckIsFromFile();

    // Get the outputfile name, output frequency and 
    // the boundary's ID for the cable's wall
    std::string typeStr = pForce->Attribute("TYPE");
    std::map<std::string, std::string> vParams;

    const TiXmlElement *param = pForce->FirstChildElement("PARAM");
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
    m_filter = MemoryManager<FilterMovingBody>::
                                    AllocateSharedPtr(m_session, vParams);

    // Initialise the object of MovingBody filter
    m_filter->Initialise(pFields, 0.0);

    // create the storage space for the body motion description
    int phystot = pFields[0]->GetTotPoints();

    for(int i = 0; i < 3; i++)
    {
        m_Forcing[i] = Array<OneD, NekDouble> (phystot, 0.0);
    }

    m_zeta = Array<OneD, Array< OneD, NekDouble> >(10);
    m_eta  = Array<OneD, Array< OneD, NekDouble> >(10);

    // What are this bi-dimensional vectors -----------------------------------
    // m_zeta[0] = zeta                 |  m_eta[0] = eta                      |
    // m_zeta[1] = d(zeta)/dt           |  m_eta[1] = d(eta)/dt                |
    // m_zeta[2] = dd(zeta)/ddtt        |  m_eta[2] = dd(eta)/ddtt             |
    // m_zeta[3] = d(zeta)/dz           |  m_eta[3] = d(eta)/dz                |
    // m_zeta[4] = dd(zeta)/ddzz        |  m_eta[4] = dd(eta)/ddzz             |
    // m_zeta[5] = ddd(zeta)/dddzzz     |  m_eta[5] = ddd(eta)/dddzzz          |
    // m_zeta[6] = dd(zeta)/ddtz        |  m_eta[6] = dd(eta)/ddtz             |
    // m_zeta[7] = ddd(zeta)/dddtzz     |  m_eta[7] = ddd(eta)/dddtzz          |
    // m_zeta[8] = (d(zeta)/dz)(d(eta)/dz)|m_eta[8] = (d(eta)/dz)(d(zeta)/dz)  |
    // m_zeta[9] = (d(zeta)/dz)^2       |  m_eta[9] = (d(eta)/dz)^2            |
    //-------------------------------------------------------------------------

    for(int i = 0; i < m_zeta.num_elements(); i++)
    {
        m_zeta[i] = Array<OneD, NekDouble>(phystot,0.0);
        m_eta[i]  = Array<OneD, NekDouble>(phystot,0.0);
    }

    // Initialise the cable model
    InitialiseCableModel(m_session, pFields);
}

void ForcingMovingBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  fields,
        const Array<OneD, Array<OneD, NekDouble> >&         inarray,
              Array<OneD, Array<OneD, NekDouble> >&         outarray,
        const NekDouble&                                    time)
{
    // Fill in m_zeta and m_eta with all the required values
    UpdateMotion(fields,time);
    //calcualte the forcing components Ax,Ay,Az and put them in m_Forcing
    CalculateForcing(fields);
    // Apply forcing terms
    for (int i = 0; i < 3; i++)
    {
        Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                                               m_Forcing[i], 1, outarray[i], 1);
    }
}

/**
 *
 */
void ForcingMovingBody::UpdateMotion(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  pFields,
              NekDouble                                     time)
{
    // Update the forces from the calculation of fluid field, which is
    // implemented in the movingbody filter
    m_filter->UpdateForce(m_session, pFields, m_Aeroforces, time);

    // for "free" type, the cable vibrates both in streamwise and crossflow
    // directions, for "constrained" type, the cable only vibrates in crossflow
    // direction, and for "forced" type, the calbe vibrates specifically along
    // a given function or file
    if(m_vibrationtype == "Free" || m_vibrationtype == "FREE" ||
        m_vibrationtype == "Constrained" || m_vibrationtype == "CONSTRAINED")
    {
        // For free vibration case, displacements, velocities and acceleartions
        // are obtained through solving structure dynamic model
        EvaluateStructDynModel(pFields, time);
    }
    else if(m_vibrationtype == "Forced" || m_vibrationtype == "FORCED")
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
                EvaluateFunction(pFields, m_session, m_motion[0], m_zeta[j],
                                 m_funcName[j], time);
                EvaluateFunction(pFields, m_session, m_motion[1], m_eta[j],
                                 m_funcName[j], time);
                cnt = cnt + 2;
            }
        }

        for(int var = 0; var < 3; var++)
        {
            for(int plane = 0; plane < m_np; plane++)
            {
                int NumPhyPoints = pFields[0]->GetPlane(plane)->GetTotPoints();
                int offset       = plane * NumPhyPoints;
                int xoffset      = var * m_np+plane;
                int yoffset      = m_vsize+xoffset;

                m_MotionVars[xoffset] = m_zeta[var][offset];
                m_MotionVars[yoffset] = m_eta [var][offset];
            }
        }
    }
    else
    {
        ASSERTL0(false, "Unrecogized vibration type for cable's dynamic model");
    }

    // Pass the variables of the cable's motion to the movingbody filter
    m_filter->UpdateMotion(m_session, pFields, m_MotionVars, time);

    // Now we need to calcualte all the required z-derivatives
    bool OriginalWaveSpace = pFields[0]->GetWaveSpace();
    pFields[0]->SetWaveSpace(false);

    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zeta[0], m_zeta[3]); //d(zeta)/dz
    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zeta[3], m_zeta[4]); //dd(zeta)/ddzz
    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zeta[4], m_zeta[5]); //ddd(zeta)/dddzzz

    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[0], m_eta[3]); //d(eta)/dz
    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[3], m_eta[4]); //dd(eta)/ddzz
    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[4], m_eta[5]); //ddd(eta)/dddzzz

    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zeta[1], m_zeta[6]); //dd(zeta)/ddtz
    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_zeta[6], m_zeta[7]); //ddd(zeta)/ddtzz

    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[1], m_eta[6]); //dd(eta)/ddtz
    pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                            m_eta[6], m_eta[7]); //ddd(eta)/ddtzz

    int NumPoints = pFields[0]->GetTotPoints();

    // (d(zeta)/dz)(d(eta)/dz)
    Vmath::Vmul(NumPoints, m_zeta[3], 1, m_eta[3],  1, m_zeta[8], 1);
    // (d(eta)/dz)(d(zeta)/dz) // not really needed
    Vmath::Vmul(NumPoints, m_eta[3],  1, m_zeta[3], 1, m_eta[8],  1);

    // (d(zeta)/dz)^2
    Vmath::Vmul(NumPoints, m_zeta[3], 1, m_zeta[3], 1, m_zeta[9], 1);
    // (d(eta)/dz)^2
    Vmath::Vmul(NumPoints, m_eta[3],  1, m_eta[3],  1, m_eta[9],  1);

    pFields[0]->SetWaveSpace(OriginalWaveSpace);
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
    // x-component of the forcing - accelaration component
    //-------------------------------------------------------------------------
    Vmath::Vsub(nPointsTot, Fx,        1, m_zeta[2], 1, Fx,   1);
    Vmath::Vmul(nPointsTot, m_zeta[3], 1, tmp,       1, tmp1, 1);
    Vmath::Vsub(nPointsTot, Fx,        1, tmp1,      1, Fx,   1);

    Vmath::Vmul(nPointsTot, m_zeta[6], 1, W,         1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Fx,        1, tmp2,      1, Fx,   1);
    // W^2 - we reuse it later
    Vmath::Vmul(nPointsTot, W,         1, W,         1, tmp3, 1);
    Vmath::Vmul(nPointsTot, m_zeta[4], 1, tmp3,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Fx,        1, tmp2,      1, Fx,   1);

    //-------------------------------------------------------------------------
    // y-component of the forcing - accelaration component
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_eta[4],  1, tmp3,      1, tmp2, 1); // reusing W^2
    Vmath::Vsub(nPointsTot, Fy,        1, tmp2,      1, Fy,   1);

    Vmath::Vmul(nPointsTot, m_eta[3],  1, tmp,       1, tmp1, 1);
    Vmath::Vsub(nPointsTot, Fy,        1, tmp1,      1, Fy,   1);

    Vmath::Vsub(nPointsTot, Fy,        1, m_eta[2],  1, Fy,   1);
    Vmath::Vmul(nPointsTot, m_eta[6],  1, W,         1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Fy,        1, tmp2,      1, Fy,   1);

    //-------------------------------------------------------------------------
    // z-component of the forcing - accelaration component
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_zeta[3], 1, Px,        1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3],  1, Py,        1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Fz,        1, tmp1,      1, Fz,   1);
    Vmath::Vadd(nPointsTot, Fz,        1, tmp2,      1, Fz,   1);

    //-------------------------------------------------------------------------
    // Note: Now we use Px,Py,Pz to store the viscous component of the forcing
    // before we multiply them by m_kinvis = 1/Re. Since we build up on them we
    // need to set the entries to zero.
    //-------------------------------------------------------------------------
    Vmath::Zero(nPointsTot,Px,1);
    Vmath::Zero(nPointsTot,Py,1);
    Vmath::Zero(nPointsTot,Pz,1);

    //-------------------------------------------------------------------------
    // x-component of the forcing - viscous component1:  (U_z^'z^' - U_zz + ZETA_tz^'z^')
    //-------------------------------------------------------------------------
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[0]->GetPhys(), tmp1); // Uz
    fields[0]->HomogeneousBwdTrans(tmp1, tmp3); // Uz in physical space
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp3, tmp1); // Uzx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp3, tmp2); // Uzy

    Vmath::Vmul(nPointsTot, m_zeta[3], 1, tmp1,      1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp1,      1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3],  1, tmp2,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp2,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp1,      1, Px,   1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp2,      1, Px,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], U, tmp1); // Ux
    Vmath::Vmul(nPointsTot, m_zeta[4], 1, tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp2,      1, Px,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2); // Uxx
    Vmath::Vmul(nPointsTot, m_zeta[9], 1, tmp2,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp3,      1, Px,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Uxy
    Vmath::Vmul(nPointsTot, m_zeta[8], 1, tmp2,      1, tmp3, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp3,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp3,      1, Px,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], U, tmp1); // Uy
    Vmath::Vmul(nPointsTot, m_eta[4],  1, tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp2,      1, Px,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Uyy
    Vmath::Vmul(nPointsTot, m_eta[9],  1, tmp2,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp3,      1, Px,   1);
    Vmath::Vadd(nPointsTot, Px,        1, m_zeta[7], 1, Px,   1);

    //-------------------------------------------------------------------------
    // x-component of the forcing - viscous component2:
    //      ((ZETA_z * W)_z^'z^' + ZETA_z * (W_xx + W_yy))
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_zeta[5], 1, W,         1, tmp1, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp1,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[3], 1, Wx,        1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_zeta[4], 1, tmp1,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 3.0,          tmp2,      1, tmp3, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp3,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_eta[3],  1, Wy,        1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_zeta[4], 1, tmp1,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp2,      1, tmp3, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp3,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[3], 1, Wy,        1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[4],  1, tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[4], 1, Wz,        1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp1,      1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[9], 1, Wxy,       1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3],  1, tmp1,      1, tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp2,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp3,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[9], 1, Wxz,       1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[8], 1, Wyz,       1, tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,          tmp1,      1, tmp2, 1);
    Vmath::Vsub(nPointsTot, Px,        1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[9], 1, m_zeta[3], 1, tmp1, 1);
    Vmath::Vmul(nPointsTot, tmp1,      1, Wxx,       1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp2,      1, Px,   1);

    Vmath::Vmul(nPointsTot, m_zeta[3], 1, Wyy,       1, tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[9],  1, tmp1,      1, tmp2, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp2,      1, Px,   1);

    Vmath::Vadd(nPointsTot, Wxx,       1, Wyy,       1, tmp1, 1);
    Vmath::Vadd(nPointsTot, Wzz,       1, tmp1,      1, tmp2, 1);
    Vmath::Vmul(nPointsTot, m_zeta[3], 1, tmp2,      1, tmp3, 1);
    Vmath::Vadd(nPointsTot, Px,        1, tmp3,      1, Px,   1);

    Vmath::Smul(nPointsTot, m_kinvis,     Px,        1, Px,   1); //* 1/Re

    //-------------------------------------------------------------------------
    // y-component of the forcing - viscous component1:
    //  (V_z^'z^' - V_zz + ETA_tz^'z^')
    //-------------------------------------------------------------------------
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                         fields[1]->GetPhys(), tmp1); // Vz
    fields[0]->HomogeneousBwdTrans(tmp1, tmp3); // Vz in physical space
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp3, tmp1); // Vzx
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp3, tmp2); // Vzy
    Vmath::Vmul(nPointsTot, m_zeta[3],  1,  tmp1,   1,  tmp1, 1);
    Vmath::Smul(nPointsTot, 2.0,            tmp1,   1,  tmp1, 1);
    Vmath::Vmul(nPointsTot, m_eta[3],   1,  tmp2,   1,  tmp2, 1);
    Vmath::Smul(nPointsTot, 2.0,            tmp2,   1,  tmp2, 1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp1,   1,  Py,   1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp2,   1,  Py,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], V, tmp1); // Vx
    Vmath::Vmul(nPointsTot, m_zeta[4],  1,  tmp1,   1,  tmp2, 1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp2,   1,  Py,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], tmp1, tmp2); // Vxx
    Vmath::Vmul(nPointsTot, m_zeta[9],  1,  tmp2,   1,  tmp3, 1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp3,   1,  Py,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Vxy
    Vmath::Vmul(nPointsTot, m_zeta[8],  1,  tmp2,   1,  tmp3, 1);
    Vmath::Smul(nPointsTot, 2.0,            tmp3,   1,  tmp3, 1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp3,   1,  Py,   1);

    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], V, tmp1); // Vy
    Vmath::Vmul(nPointsTot, m_eta[4],   1,  tmp1,   1,  tmp2, 1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp2,   1,  Py,   1);
    fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], tmp1, tmp2); // Vyy
    Vmath::Vmul(nPointsTot, m_eta[9],   1,  tmp2,   1,  tmp3, 1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp3,   1,  Py,   1);
    Vmath::Vadd(nPointsTot, m_eta[7],   1,  Py,     1,  Py,   1);

    //-------------------------------------------------------------------------
    // y-component of the forcing - viscous component2:
    //   ((ETA_z * W)_z^'z^' + ETA_z * (W_xx + W_yy))
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_eta[5],   1,  W,          1,  tmp1,   1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp1,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_zeta[3],  1,  Wx,         1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_eta[4],   1,  tmp1,       1,  tmp2,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp2,       1,  tmp3,   1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp3,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_zeta[4],  1,  Wx,         1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_eta[3],   1,  tmp1,       1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp2,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[3],   1,  Wy,         1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_eta[4],   1,  tmp1,       1,  tmp2,   1);
    Vmath::Smul(nPointsTot, 3.0,            tmp2,       1,  tmp3,   1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp3,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[4],   1,  Wz,         1,  tmp1,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp1,       1,  tmp2,   1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp2,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[9],   1,  Wxy,        1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_zeta[3],  1,  tmp1,       1,  tmp2,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp2,       1,  tmp3,   1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp3,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_zeta[8],  1,  Wxz,        1,  tmp1,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp1,       1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp2,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[9],   1,  Wyz,        1,  tmp1,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp1,       1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Py,         1,  tmp2,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[3],   1,  Wxx,        1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_zeta[9],  1,  tmp1,       1,  tmp2,   1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp2,       1,  Py,     1);

    Vmath::Vmul(nPointsTot, m_eta[9],   1,  m_eta[3],   1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, tmp1,       1,  Wyy,        1,  tmp2,   1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp2,       1,  Py,     1);

    Vmath::Vadd(nPointsTot, Wxx,        1,  Wyy,        1,  tmp1,   1);
    Vmath::Vadd(nPointsTot, Wzz,        1,  tmp1,       1,  tmp2,   1);
    Vmath::Vmul(nPointsTot, m_eta[3],   1,  tmp2,       1,  tmp3,   1);
    Vmath::Vadd(nPointsTot, Py,         1,  tmp3,       1,  Py,     1);

    Vmath::Smul(nPointsTot, m_kinvis,       Py,         1,  Py,     1); //* 1/Re

    //-------------------------------------------------------------------------
    // z-component of the forcing - viscous component: (W_z^'z^' - W_zz)
    //-------------------------------------------------------------------------
    Vmath::Vmul(nPointsTot, m_zeta[3],  1,  Wxz,        1,  tmp1,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp1,       1,  tmp1,   1);
    Vmath::Vmul(nPointsTot, m_eta[3],   1,  Wyz,        1,  tmp2,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp2,       1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Pz,         1,  tmp1,       1,  Pz,     1);
    Vmath::Vsub(nPointsTot, Pz,         1,  tmp2,       1,  Pz,     1);

    Vmath::Vmul(nPointsTot, m_zeta[4],  1,  Wx,         1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Pz,         1,  tmp2,       1,  Pz,     1);
    Vmath::Vmul(nPointsTot, m_zeta[9],  1,  Wxx,        1,  tmp3,   1);
    Vmath::Vadd(nPointsTot, Pz,         1,  tmp3,       1,  Pz,     1);
    Vmath::Vmul(nPointsTot, m_zeta[8],  1,  Wxy,        1,  tmp3,   1);
    Vmath::Smul(nPointsTot, 2.0,            tmp3,       1,  tmp3,   1);
    Vmath::Vadd(nPointsTot, Pz,         1,  tmp3,       1,  Pz,     1);

    Vmath::Vmul(nPointsTot, m_eta[4],   1,  Wy,         1,  tmp2,   1);
    Vmath::Vsub(nPointsTot, Pz,         1,  tmp2,       1,  Pz,     1);
    Vmath::Vmul(nPointsTot, m_eta[9],   1,  Wyy,        1,  tmp3,   1);
    Vmath::Vadd(nPointsTot, Pz,         1,  tmp3,       1,  Pz,     1);

    Vmath::Smul(nPointsTot, m_kinvis,       Pz,         1,  Pz,     1); //* 1/Re

    //-------------------------------------------------------------------------
    // adding viscous and pressure components and transfroming back to wave
    // space
    //-------------------------------------------------------------------------
    Vmath::Vadd(nPointsTot, Fx,         1,  Px,         1,  Fx,     1);
    Vmath::Vadd(nPointsTot, Fy,         1,  Py,         1,  Fy,     1);
    Vmath::Vadd(nPointsTot, Fz,         1,  Pz,         1,  Fz,     1);
    fields[0]->HomogeneousFwdTrans(Fx, m_Forcing[0]);
    fields[0]->HomogeneousFwdTrans(Fy, m_Forcing[1]);
    fields[0]->HomogeneousFwdTrans(Fz, m_Forcing[2]);
}


/**
 *
 */
void ForcingMovingBody::TensionedCableModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble> &AeroForces,
              Array<OneD, NekDouble> &CableMotions)
{   
    // m_homostrip == false indicates that a single fourier transformation is
    // implemented for the motion of cable, so it is matched with that carried
    // out in fluid fields; on the other hand, if m_homostrip == true, there is
    // a mismatch between the structure and fluid fields in terms of fourier
    // transformation, then the routines such as FFTFwdTrans/
    // FFTBwdTrans can not be used directly for cable's solution.

    int npts;

    if(!m_homostrip) //full resolutions
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

    int colrank = m_comm->GetColumnComm()->GetRank();
    int nproc   = m_comm->GetColumnComm()->GetSize();
    
    if(!m_homostrip) //full resolutions
    {
        Array <OneD, NekDouble> tmp;
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
                m_comm->GetColumnComm()->Recv(i, AeroForces);
                m_comm->GetColumnComm()->Recv(i, CableMotions);

                Vmath::Vcopy(m_np, AeroForces, 1, 
                        tmp = fft_i[0]+i*m_np, 1);
                for(int j = 0; j < 3; j++)
                {
                    Vmath::Vcopy(m_np, CableMotions+j*m_np, 1, 
                                            tmp = fft_i[j+1]+i*m_np, 1);
                }
            }
        }
        else
        {
            m_comm->GetColumnComm()->Send(0, AeroForces);
            m_comm->GetColumnComm()->Send(0, CableMotions);
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
                m_comm->GetColumnComm()->Recv(i, tmp);

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
                    m_comm->GetColumnComm()->Send(0, tmp);
                }
            }
        }
    }
    
    if(colrank == 0)
    {
        // Implement Fourier transformation of the motion variables
        if(m_supporttype == "FREE-FREE" || m_supporttype == "Free-Free")
        {
            for(int j = 0 ; j < 4; ++j)
            {
                m_FFT->FFTFwdTrans(fft_i[j], fft_o[j]);
            }
        }
        else if(m_supporttype == "PINNED-PINNED" || m_supporttype == "Pinned-Pinned")
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
        if(m_supporttype == "FREE-FREE" || m_supporttype == "Free-Free")
        {
            for(int var = 0; var < 3; var++)
            {
                m_FFT->FFTBwdTrans(fft_i[var], fft_o[var]);
            }
        }
        else if(m_supporttype == "PINNED-PINNED" || m_supporttype == "Pinned-Pinned")
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
    if(!m_homostrip)//full resolutions
    {
        Array<OneD, NekDouble> tmp;

        if(colrank != 0)
        {
            m_comm->GetColumnComm()->Recv(0, CableMotions);
        }
        else
        {
            for (int i = 1; i < nproc; ++i)
            {
                for(int j = 0; j < 3; j++)
                {
                    Vmath::Vcopy(m_np, fft_o[j]+i*m_np, 1, 
                                      tmp = CableMotions+j*m_np, 1);
                }
                m_comm->GetColumnComm()->Send(i, CableMotions);
            }

            for (int var = 0; var < 3; var++)
            {
                Vmath::Vcopy(m_np, fft_o[var], 1, 
                                    tmp = CableMotions+var*m_np, 1);
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
                    m_comm->GetColumnComm()->Recv(0, tmp);

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
                        m_comm->GetColumnComm()->Recv(0, tmp);

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
                m_comm->GetColumnComm()->Send(j*npts, tmp);
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
                    m_comm->GetColumnComm()->Send(i+j*npts, tmp);
                }
            }
        }
    }
}


/**
 *
 */
void ForcingMovingBody::EvaluateStructDynModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time )
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
    if(m_FictitiousMass)
    {
        // only consider second order approximation for fictitious variables
        int  intOrder= 2;
        int  nint    = min(m_movingBodyCalls+1,intOrder);
        int  nlevels = m_fV[0].num_elements();

        for(int i = 0; i < m_motion.num_elements(); ++i)
        {
            RollOver(m_fV[i]);
            RollOver(m_fA[i]);

            int Voffset = i*m_vsize+m_np;
            int Aoffset = i*m_vsize+2*m_np;

            Vmath::Vcopy(m_np, m_MotionVars+Voffset, 1, m_fV[i][0], 1);
            Vmath::Vcopy(m_np, m_MotionVars+Aoffset, 1, m_fA[i][0], 1);

            // Extrapolate to n+1
            Vmath::Smul(m_np, 
                        StifflyStable_Betaq_Coeffs[nint-1][nint-1],
                        m_fV[i][nint-1],    1,
                        m_fV[i][nlevels-1], 1);
            Vmath::Smul(m_np, 
                        StifflyStable_Betaq_Coeffs[nint-1][nint-1],
                        m_fA[i][nint-1],    1,
                        m_fA[i][nlevels-1], 1);

            for(int n = 0; n < nint-1; ++n)
            {
                Vmath::Svtvp(m_np, 
                             StifflyStable_Betaq_Coeffs[nint-1][n],
                             m_fV[i][n],1,m_fV[i][nlevels-1],1,
                             m_fV[i][nlevels-1],1);
                Vmath::Svtvp(m_np, 
                             StifflyStable_Betaq_Coeffs[nint-1][n],
                             m_fA[i][n],1,m_fA[i][nlevels-1],1,
                             m_fA[i][nlevels-1],1);
            }

            // Add the fictitious forces on the RHS of the equation
            Vmath::Svtvp(m_np, m_fictdamp,m_fV[i][nlevels-1],1,
                         fces[i],1,fces[i],1);
            Vmath::Svtvp(m_np, m_fictrho, m_fA[i][nlevels-1],1,
                         fces[i],1,fces[i],1);
        }
    }

    for(int n = 0, cn = 1; n < m_NumD; n++, cn--)
    {
        int offset = cn*m_vsize;

        Array<OneD, NekDouble> tmp(m_vsize);

        TensionedCableModel(pFields, fces[cn], tmp = m_MotionVars+offset);
    }

    // Set the forcing term based on the motion of the cable
    for(int var = 0; var < 3; var++)
    {
        for(int plane = 0; plane < m_np; plane++)
        {
            int NumPhyPoints = pFields[0]->GetPlane(plane)->GetTotPoints();

            Array<OneD, NekDouble> tmp;

            int offset  = plane * NumPhyPoints;
            int xoffset = var * m_np+plane;
            int yoffset = m_vsize + xoffset;

            Vmath::Fill(NumPhyPoints, m_MotionVars[xoffset], tmp = m_zeta[var] + offset, 1);
            Vmath::Fill(NumPhyPoints, m_MotionVars[yoffset], tmp =  m_eta[var] + offset, 1);
        }
    }

    // velocity of moving body is used to set the boundary conditions
    // according to the coordinate system mapping. for forced vibration,
    // however, that can be simply set by prior according to the moving function
    // in .xml file

    for(int i = 0; i < m_motion.num_elements(); ++i)
    {
        //RollOver(m_BndV[i]);

        int offset = i*m_vsize+m_np;
        Vmath::Vcopy(m_np, m_MotionVars+offset, 1, m_BndV[i][0], 1);

        // TODO: Extrapolate to n+1
        // TODO: Second order make the coupling unstable
        /*Vmath::Smul(m_np, StifflyStable_Betaq_Coeffs[nint-1][nint-1],
                         m_BndV[i][nint-1],    1,
                         m_BndV[i][nlevels-1], 1);

        for(int n = 0; n < nint-1; ++n)
        {
            Vmath::Svtvp(m_np, StifflyStable_Betaq_Coeffs[nint-1][n],
                         m_BndV[i][n],1,m_BndV[i][nlevels-1],1,
                         m_BndV[i][nlevels-1],1);
        }*/
    }

    Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExps;
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds;

    const Array<OneD, const NekDouble> z
                            = pFields[0]->GetHomogeneousBasis()->GetZ();
    int nbnd = pFields[0]->GetBndCondExpansions().num_elements();

    for (int i = 0; i < nbnd; ++i)
    {
        for (int plane = 0; plane < m_np; plane++)
        {
            for ( int dir = 0; dir < m_motion.num_elements(); dir++)
            {
                bndCondExps = pFields[dir]->GetPlane(plane)
                                            ->GetBndCondExpansions();
                bndConds    = pFields[dir]->GetPlane(plane)
                                            ->GetBndConditions();

                if (bndConds[i]->GetUserDefined() ==
                                            SpatialDomains::eMovingBody)
                {
                    int npoints = bndCondExps[i]->GetNpoints();
                    Array<OneD, NekDouble> x0(npoints, 0.0);
                    Array<OneD, NekDouble> x1(npoints, 0.0);
                    Array<OneD, NekDouble> x2(npoints, 0.0);
                    Array<OneD, NekDouble> tmp(npoints,0.0);

                    NekDouble local_z = z[pFields[0]->GetTransposition()
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
                    Vmath::Fill(npoints, m_BndV[dir][0][plane], tmp, 1);
                    Vmath::Vsub(npoints, bndCondExps[i]->UpdatePhys(), 1,
                                         tmp,                          1,
                                         bndCondExps[i]->UpdatePhys(), 1);
                }
            }
        }
    }
}


/**
 *
 */
void ForcingMovingBody::InitialiseCableModel(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    Array<OneD, unsigned int> ZIDs;
    ZIDs           = pFields[0]->GetZIDs();
    m_np  = ZIDs.num_elements();
    m_vsize = m_np*3;
    m_Aeroforces   = Array<OneD, NekDouble>(2*m_np,0.0);
    m_MotionVars   = Array<OneD, NekDouble>(2*m_vsize,0.0);

    m_vibrationtype = m_session->GetSolverInfo("VibrationType");

    if(m_vibrationtype == "Constrained" || m_vibrationtype == "CONSTRAINED")
    {
        m_NumD = 1;
    }
    else if (m_vibrationtype == "Free" || m_vibrationtype == "FREE")
    {
        m_NumD = 2;
    }
    else if (m_vibrationtype == "Forced" || m_vibrationtype == "FORCED")
    {
        return;
    }

    m_supporttype = m_session->GetSolverInfo("SupportType");

    m_comm        = pFields[0]->GetComm();

    m_session->MatchSolverInfo("HomoStrip","True",m_homostrip,false);

    if(!m_homostrip)
    {
        m_session->LoadParameter("LZ", m_lhom);
        int nplanes = m_session->GetParameter("HomModesZ");
        m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance(
                                                            "NekFFTW", nplanes);
    }
    else
    {
        int nstrips;

        NekDouble DistStrip;

        m_session->LoadParameter("DistStrip", DistStrip);
        m_session->LoadParameter("Strip_Z", nstrips);
        m_lhom = nstrips * DistStrip;
        m_FFT  = LibUtilities::GetNektarFFTFactory().CreateInstance(
                                                            "NekFFTW", nstrips);
    }

    // load the structural dynamic parameters from xml file
    m_session->LoadParameter("StructRho",    m_structrho);
    m_session->LoadParameter("StructDamp",   m_structdamp,   0.0);
    m_session->LoadParameter("StructStiff",  m_structstiff,  0.0);
    m_session->LoadParameter("CableTension", m_cabletension, 0.0);
    m_session->LoadParameter("BendingStiff", m_bendingstiff, 0.0);

    // only consider second order interpolation for fictitious variables
    int intOrder = 2;
    m_BndV = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                    m_motion.num_elements());
    for (int i = 0; i < m_motion.num_elements(); ++i)
    {
        m_BndV[i] = Array<OneD, Array<OneD, NekDouble> > (intOrder);
        for(int n = 0; n < intOrder; ++n)
        {
            m_BndV[i][n] = Array<OneD, NekDouble>(m_np, 0.0);
        }
    }

    // Identify whether the fictitious mass method is active for explicit
    // coupling of fluid solver and structural dynamics solver
    m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                m_FictitiousMass,      false);
    if(m_FictitiousMass)
    {
        m_session->LoadParameter("FictMass", m_fictrho, 1.5/m_structrho);
        m_session->LoadParameter("FictDamp", m_fictdamp, 0.1);
        m_structrho  += m_fictrho;
        m_structdamp += m_fictdamp;

        // Storage array of Struct Velocity and Acceleration used for
        // extrapolation of fictitious force
        m_fV = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                    m_motion.num_elements());
        m_fA = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (
                                                    m_motion.num_elements());
        for (int i = 0; i < m_motion.num_elements(); ++i)
        {
            m_fV[i] = Array<OneD, Array<OneD, NekDouble> > (intOrder);
            m_fA[i] = Array<OneD, Array<OneD, NekDouble> > (intOrder);

            for(int n = 0; n < intOrder; ++n)
            {
                m_fV[i][n] = Array<OneD, NekDouble>(m_np, 0.0);
                m_fA[i][n] = Array<OneD, NekDouble>(m_np, 0.0);
            }
        }
    }

    // Setting the coefficient matrices for solving structural dynamic ODEs
    SetDynEqCoeffMatrix(pFields);

    // Set initial condition for cable's motion
    int cnt = 0;

    for(int j = 0; j < m_funcName.num_elements(); j++)
    {
                
        // loading from the specified files through inputstream
        if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
        {
    
            std::ifstream inputStream;

            int nzpoints;

            if(m_homostrip)
            {
                int nstrips;
                m_session->LoadParameter("Strip_Z", nstrips);
                nzpoints = nstrips;
            }
            else
            { 
                nzpoints = pFields[0]->GetHomogeneousBasis()->GetNumModes();
            }

 
            if (m_comm->GetRank() == 0)
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
                
                if(!m_homostrip)
                {
                    int nproc = m_comm->GetColumnComm()->GetSize();
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

                        m_comm->GetColumnComm()->Send(i, m_MotionVars);
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

                    int nproc = m_comm->GetColumnComm()->GetSize();
                    
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
                            m_comm->GetColumnComm()->Send(i+j*nstrips, m_MotionVars);
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
                        m_comm->GetColumnComm()->Send(j*nstrips, m_MotionVars);
                    }
                }
            }
            else
            {
                if(!m_homostrip)
                {
                    int colrank = m_comm->GetColumnComm()->GetRank();
                    int nproc   = m_comm->GetColumnComm()->GetSize();

                    for (int j = 1; j < nproc; j++)
                    {
                        if(colrank == j)
                        {
                            m_comm->GetColumnComm()->Recv(0, m_MotionVars);
                        }
                    }
                }
                else //for strip model
                {
                    int nstrips;
                    m_session->LoadParameter("Strip_Z", nstrips);

                    int colrank = m_comm->GetColumnComm()->GetRank();
                    int nproc   = m_comm->GetColumnComm()->GetSize();

                    for(int i = 1; i < nstrips; i++)
                    {
                        for (int j = 0; j < nproc/nstrips; j++)
                        {
                            if(colrank == i+j*nstrips)
                            {
                                m_comm->GetColumnComm()->Recv(0, m_MotionVars);
                            }
                        }
                    }

                    for (int j = 1; j < nproc/nstrips; j++)
                    {
                        if(colrank == j*nstrips)
                        {
                            m_comm->GetColumnComm()->Recv(0, m_MotionVars);
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
            
            if(!m_homostrip)
            {
                int nzpoints = pFields[0]->GetHomogeneousBasis()->GetNumModes();
                Array<OneD, NekDouble> z_coords(nzpoints,0.0);
                Array<OneD, const NekDouble> pts
                                    = pFields[0]->GetHomogeneousBasis()->GetZ();

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

                int colrank = m_comm->GetColumnComm()->GetRank();
                int nproc   = m_comm->GetColumnComm()->GetSize();
                
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

            int Xoffset = j*m_np;
            int Yoffset = m_vsize+Xoffset;

            Array<OneD, NekDouble> tmp(m_np,0.0);
            LibUtilities::EquationSharedPtr ffunc0,ffunc1;

            ffunc0 = m_session->GetFunction(m_funcName[j], m_motion[0]);
            ffunc1 = m_session->GetFunction(m_funcName[j], m_motion[1]);

            ffunc0->Evaluate(x0, x1, x2, 0.0, tmp = m_MotionVars+Xoffset);
            ffunc1->Evaluate(x0, x1, x2, 0.0, tmp = m_MotionVars+Yoffset);
            cnt = cnt + 2;
        }
    }
}


/**
 *
 */
void ForcingMovingBody::SetDynEqCoeffMatrix(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    int nplanes;

    if(!m_homostrip)
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
    tmp1 =     m_timestep * m_timestep;
    tmp2 =  m_structstiff / m_structrho;
    tmp3 =   m_structdamp / m_structrho;
    tmp4 = m_cabletension / m_structrho;
    tmp5 = m_bendingstiff / m_structrho;

    // solve the ODE in the wave space for cable motion to obtain disp, vel and
    // accel
    for(int plane = 0; plane < nplanes; plane++)
    {
        int nel = 3;
        m_CoeffMat_A[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);
        m_CoeffMat_B[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);

        unsigned int K;
        NekDouble beta;

        if(m_supporttype == "FREE-FREE" || 
                    m_supporttype == "Free-Free")
        {
            K = plane/2;
            beta = 2.0 * M_PI/m_lhom;
        }
        else if(m_supporttype == "PINNED-PINNED" || 
                        m_supporttype == "Pinned-Pinned")
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

    int acc_order = 0;

    // Rotate acceleration term
    RollOver(m_W);

    Vmath::Vcopy(npoints, input, 1, m_W[0], 1);

    //Calculate acceleration term at level n based on previous steps
    if (m_movingBodyCalls > 2)
    {
        acc_order = min(m_movingBodyCalls-2,m_intSteps);
        Vmath::Smul(npoints,
                    StifflyStable_Gamma0_Coeffs[acc_order-1],
                    m_W[0], 1,
                    output, 1);
        for(int i = 0; i < acc_order; i++)
        {
            Vmath::Svtvp(npoints,
                         -1*StifflyStable_Alpha_Coeffs[acc_order-1][i],
                         m_W[i+1], 1,
                         output,   1,
                         output,   1);
        }
        Vmath::Smul(npoints,
                    1.0/m_timestep,
                    output, 1,
                    output, 1);
    }
}


/**
 *
 */
void ForcingMovingBody::CheckIsFromFile()
{
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

}
