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
// Description: Moving Body forcing (movement of a body in a domain is achieved via
// a forcing term which is the results of a coordinate system motion)
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingMovingBody.h>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
namespace SolverUtils
{
    std::string ForcingMovingBody::className = GetForcingFactory().
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
		ASSERTL0(pFields[0]->GetExpType()==MultiRegions::e3DH1D,"Moving body implemented just for"
                 "3D Homogenous 1D expansions.");

        int nPointsTot    = pFields[0]->GetNpoints();

        m_NumVariable     = pNumForcingFields; // forcing size (it must be 3)
	m_Forcing         = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        m_funcName        = Array<OneD, std::string> (3);
        m_motion          = Array<OneD, std::string> (2);
        m_motion[0]       = "x";
        m_motion[1]       = "y";
        m_movingBodyCalls = 0;
        m_IsFromFile      = Array<OneD, bool> (6);        
        m_session->LoadParameter("Kinvis",m_kinvis);
        m_session->LoadParameter("TimeStep", m_timestep, 0.01); 
        
        m_StifflyStable_Alpha_Coeffs = Array<OneD, Array<OneD, NekDouble> > (3);
        m_StifflyStable_Alpha_Coeffs[0] = Array<OneD, NekDouble> (3);
        m_StifflyStable_Alpha_Coeffs[1] = Array<OneD, NekDouble> (3);
        m_StifflyStable_Alpha_Coeffs[2] = Array<OneD, NekDouble> (3);
        m_StifflyStable_Gamma0_Coeffs = Array<OneD, NekDouble> (3);
        m_StifflyStable_Alpha_Coeffs[0][0] = 1.0;
        m_StifflyStable_Alpha_Coeffs[0][1] = 0.0;
        m_StifflyStable_Alpha_Coeffs[0][2] = 0.0;
        m_StifflyStable_Alpha_Coeffs[1][0] = 2.0;
        m_StifflyStable_Alpha_Coeffs[1][1] =-0.5;
        m_StifflyStable_Alpha_Coeffs[1][2] = 0.0;
        m_StifflyStable_Alpha_Coeffs[2][0] = 3.0;
        m_StifflyStable_Alpha_Coeffs[2][1] =-1.5;
        m_StifflyStable_Alpha_Coeffs[2][2] = 1.0/3.0;
        m_StifflyStable_Gamma0_Coeffs[0]   = 1.0;
        m_StifflyStable_Gamma0_Coeffs[1]   = 1.5;
        m_StifflyStable_Gamma0_Coeffs[2]   = 11.0/6.0;
       
        if (m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"))
        {
            std::string TimeInMethod = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
            if ((TimeInMethod == "IMEXORDER1") || (TimeInMethod == "IMEXOrder1"))
            {
                m_intSteps = 1;
            }
            else if ((TimeInMethod == "IMEXORDER2") || (TimeInMethod == "IMEXOrder2"))
            {
                m_intSteps = 2;
            }
            else if ((TimeInMethod == "IMEXORDER3") || (TimeInMethod == "IMEXOrder3"))
            {
                m_intSteps = 3;
            }
            else
            {
                ASSERTL0(false, "Unrecognised time integration method for evaluation of MovingBody forcing");
            }
        }
        
        m_acceleration    = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);
        m_acceleration[0] = Array<OneD, NekDouble>(nPointsTot, 0.0);
       
        for(int n = 0; n < m_intSteps; ++n)
        {
            m_acceleration[n+1] = Array<OneD, NekDouble>(nPointsTot, 0.0);
        }        

        // Loading the x-dispalcement (m_zeta) and the y-displacement (m_eta)
		// Those two variables are bith functions of z and t and the may come
		// from an equation (forced vibration) or from another solver which, given
		// the aerodynamic forces at the previous step, calculates the displacements.
        
        //Get the body displacement: m_zeta and m_eta
        const TiXmlElement* funcNameElmt_D = pForce->FirstChildElement("DISPLACEMENTS");
        ASSERTL0(funcNameElmt_D,"MOVINGBODYFORCE tag has to specify a function name which " 
                                "prescribes the body displacement as d(z,t).");

        m_funcName[0] = funcNameElmt_D->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName[0]),
                 "Function '" + m_funcName[0] + "' not defined.");
        
        //Get the body velocity of movement: d(m_zeta)/dt and d(m_eta)/dt
        const TiXmlElement* funcNameElmt_V = pForce->FirstChildElement("VELOCITIES");
        ASSERTL0(funcNameElmt_D,"MOVINGBODYFORCE tag has to specify a function name which " 
                 "prescribes the body velocity of movement as v(z,t).");
        
        m_funcName[1] = funcNameElmt_V->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName[1]),
                 "Function '" + m_funcName[1] + "' not defined.");
        
        
        //Get the body acceleration: dd(m_zeta)/ddt and dd(m_eta)/ddt
        const TiXmlElement* funcNameElmt_A = pForce->FirstChildElement("ACCELERATIONS");
        ASSERTL0(funcNameElmt_A,"MOVINGBODYFORCE tag has to specify a function name which " 
                 "prescribes the body acceleration as a(z,t).");
        
        m_funcName[2] = funcNameElmt_A->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName[2]),
                 "Function '" + m_funcName[2] + "' not defined.");
        
        // At this point we know in the xml file where those quantities
        // are declared (equation or file) - via a function name which is now stored in funcNameD etc.
        // We need now to fill in with this info the m_zeta and m_eta vectors (actuallythey are matrices)
        // Array to control if the motion is determined by an equation or is from a file.(not Nektar++)
        // check if we need to load a file or we have an equation
        CheckIsFromFile();
        
        // create the storage space for the body motion description
        int phystot = pFields[0]->GetTotPoints();
        
        for(int i = 0; i < m_NumVariable; i++)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (phystot, 0.0);
        }
        
        m_zeta = Array<OneD, Array< OneD, NekDouble> >(10);
        m_eta  = Array<OneD, Array< OneD, NekDouble> >(10);
        
        // What are this bi-dimensional vectors ----------------------------------------
        // m_zeta[0] = zeta                     |  m_eta[0] = eta                      |
        // m_zeta[1] = d(zeta)/dt               |  m_eta[1] = d(eta)/dt                |
        // m_zeta[2] = dd(zeta)/ddtt            |  m_eta[2] = dd(eta)/ddtt             |
        // m_zeta[3] = d(zeta)/dz               |  m_eta[3] = d(eta)/dz                |
        // m_zeta[4] = dd(zeta)/ddzz            |  m_eta[4] = dd(eta)/ddzz             |
        // m_zeta[5] = ddd(zeta)/dddzzz         |  m_eta[5] = ddd(eta)/dddzzz          |
        // m_zeta[6] = dd(zeta)/ddtz            |  m_eta[6] = dd(eta)/ddtz             |
        // m_zeta[7] = ddd(zeta)/dddtzz         |  m_eta[7] = ddd(eta)/dddtzz          |
        // m_zeta[8] = (d(zeta)/dz)(d(eta)/dz)  |  m_eta[8] = (d(eta)/dz)(d(zeta)/dz)  |
        // m_zeta[9] = (d(zeta)/dz)^2           |  m_eta[9] = (d(eta)/dz)^2            |
        //------------------------------------------------------------------------------
        
        for(int i = 0; i < m_zeta.num_elements(); i++)
        {
            m_zeta[i] = Array<OneD, NekDouble>(phystot,0.0);
            m_eta[i]  = Array<OneD, NekDouble>(phystot,0.0);
        }
    }

    void ForcingMovingBody::v_Apply(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                    const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                                    NekDouble time)
    {
	// Fill in m_zeta and m_eta with all the required values       
        UpdateMotion(fields,time);
		
	//calcualte the forcing components Ax,Ay,Az and put them in m_Forcing
	CalculateForcing(fields);
		
	// Apply forcing terms
        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }
    }
				 				 
    void ForcingMovingBody::UpdateMotion(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                                         NekDouble time)
    {
        // Loading displacements, velocities, acceleration
        int cnt = 0;
        for(int j = 0; j < m_funcName.num_elements();j++)
        {
            if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
            {
               ASSERTL0(false,"Motion loading from file needs specific implementation: Work in Progress!"); 
            }
            else 
            {
                EvaluateFunction(pFields, m_session, m_motion[0],m_zeta[j],m_funcName[j],time);
                EvaluateFunction(pFields, m_session, m_motion[1],m_eta[j],m_funcName[j],time);
                cnt = cnt + 2;
            }
        }
        
        // Now we need to calcualte all the required z-derivatives
        bool OriginalWaveSpace = pFields[0]->GetWaveSpace();
        pFields[0]->SetWaveSpace(false);
        
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_zeta[0],m_zeta[3]); //d(zeta)/dz
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_zeta[3],m_zeta[4]); //dd(zeta)/ddzz
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_zeta[4],m_zeta[5]); //ddd(zeta)/dddzzz
        
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_eta[0],m_eta[3]); //d(eta)/dz
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_eta[3],m_eta[4]); //dd(eta)/ddzz
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_eta[4],m_eta[5]); //ddd(eta)/dddzzz
        
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_zeta[1],m_zeta[6]); //dd(zeta)/ddtz
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_zeta[6],m_zeta[7]); //ddd(zeta)/ddtzz
        
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_eta[1],m_eta[6]); //dd(eta)/ddtz
        pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],m_eta[6],m_eta[7]); //ddd(eta)/ddtzz
        
        int NumPoints = pFields[0]->GetTotPoints();
        
        Vmath::Vmul(NumPoints,m_zeta[3],1,m_eta[3],1,m_zeta[8],1); //(d(zeta)/dz)(d(eta)/dz)
        Vmath::Vmul(NumPoints,m_eta[3],1,m_zeta[3],1,m_eta[8],1); //(d(eta)/dz)(d(zeta)/dz) // not really needed
        
        Vmath::Vmul(NumPoints,m_zeta[3],1,m_zeta[3],1,m_zeta[9],1); //(d(zeta)/dz)^2
        Vmath::Vmul(NumPoints,m_eta[3],1,m_eta[3],1,m_eta[9],1); //(d(eta)/dz)^2
                    
        pFields[0]->SetWaveSpace(OriginalWaveSpace);
    }
	
    void ForcingMovingBody::CalculateForcing(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
    {

        int nPointsTot = fields[0]->GetNpoints();
        Array<OneD, NekDouble> U,V,W;
        Array<OneD, NekDouble> Wt,Wx,Wy,Wz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz;
        Array<OneD, NekDouble> Px,Py,Pz;
        Array<OneD, NekDouble> tmp,tmp1,tmp2,tmp3;
        Array<OneD, NekDouble> Fx,Fy,Fz;
        U = Array<OneD, NekDouble> (nPointsTot);
        V = Array<OneD, NekDouble> (nPointsTot);
        W = Array<OneD, NekDouble> (nPointsTot);
        Wt = Array<OneD, NekDouble> (nPointsTot);
        Wx = Array<OneD, NekDouble> (nPointsTot);
        Wy = Array<OneD, NekDouble> (nPointsTot);
        Wz = Array<OneD, NekDouble> (nPointsTot);
        Wxx = Array<OneD, NekDouble> (nPointsTot);
        Wxy = Array<OneD, NekDouble> (nPointsTot);
        Wyy = Array<OneD, NekDouble> (nPointsTot);
        Wxz = Array<OneD, NekDouble> (nPointsTot);
        Wyz = Array<OneD, NekDouble> (nPointsTot);
        Wzz = Array<OneD, NekDouble> (nPointsTot);
        Px = Array<OneD, NekDouble> (nPointsTot);
        Py = Array<OneD, NekDouble> (nPointsTot);
        Pz = Array<OneD, NekDouble> (nPointsTot);
        Fx = Array<OneD, NekDouble> (nPointsTot,0.0);
        Fy = Array<OneD, NekDouble> (nPointsTot,0.0);
        Fz = Array<OneD, NekDouble> (nPointsTot,0.0);
        tmp  = Array<OneD, NekDouble> (nPointsTot,0.0);
        tmp1 = Array<OneD, NekDouble> (nPointsTot);
        tmp2 = Array<OneD, NekDouble> (nPointsTot);
        tmp3 = Array<OneD, NekDouble> (nPointsTot);
        fields[0]->HomogeneousBwdTrans(fields[0]->GetPhys(),U);
        fields[0]->HomogeneousBwdTrans(fields[1]->GetPhys(),V);
        fields[0]->HomogeneousBwdTrans(fields[2]->GetPhys(),W);
        fields[0]->HomogeneousBwdTrans(fields[3]->GetPhys(),tmp1); // pressure

        //-------------------------------------------------------------------------------------------------
        // Setting the pressure derivatives
        //-------------------------------------------------------------------------------------------------
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],tmp1,Px); // Px
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp1,Py); // Py
     
        //-------------------------------------------------------------------------------------------------
        // Setting the z-component velocity derivatives
        //-------------------------------------------------------------------------------------------------
        EvaluateAccelaration(W,Wt,nPointsTot); //Wt
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],W,Wx); // Wx
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],W,Wy); // Wy
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],fields[2]->GetPhys(),tmp1); // Wz
        fields[0]->HomogeneousBwdTrans(tmp1,Wz); // Wz in physical space
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],Wx,Wxx); // Wxx
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],Wx,Wxy); // Wxy
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],Wy,Wyy); // Wyy
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],Wz,Wxz); // Wxz
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],Wz,Wyz); // Wyz
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],tmp1,tmp2); // Wzz
        fields[0]->HomogeneousBwdTrans(tmp2,Wzz); // Wzz in physical space
        //-------------------------------------------------------------------------------------------------
        // Setting the z-component of the accelaration term: tmp = (Wt + U * Wx + V * Wy + W * Wz) 
        //-------------------------------------------------------------------------------------------------
        Vmath::Vadd(nPointsTot,tmp,1,Wt,1,tmp,1); 
        Vmath::Vmul(nPointsTot,U,1,Wx,1,tmp1,1);
        Vmath::Vadd(nPointsTot,tmp,1,tmp1,1,tmp,1);
        Vmath::Vmul(nPointsTot,V,1,Wy,1,tmp1,1);
        Vmath::Vadd(nPointsTot,tmp,1,tmp1,1,tmp,1);
        Vmath::Vmul(nPointsTot,W,1,Wz,1,tmp1,1);
        Vmath::Vadd(nPointsTot,tmp,1,tmp1,1,tmp,1);

        //-------------------------------------------------------------------------------------------------
        // x-component of the forcing - accelaration component
        //-------------------------------------------------------------------------------------------------
        Vmath::Vsub(nPointsTot,Fx,1,m_zeta[2],1,Fx,1);
        
        Vmath::Vmul(nPointsTot,m_zeta[3],1,tmp,1,tmp1,1);
        Vmath::Vsub(nPointsTot,Fx,1,tmp1,1,Fx,1);     
 
        Vmath::Vmul(nPointsTot,m_zeta[6],1,W,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Fx,1,tmp2,1,Fx,1);
        Vmath::Vmul(nPointsTot,W,1,W,1,tmp3,1); //W^2 - we reuse it later
        Vmath::Vmul(nPointsTot,m_zeta[4],1,tmp3,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Fx,1,tmp2,1,Fx,1);

        //-------------------------------------------------------------------------------------------------
        // y-component of the forcing - accelaration component
        //-------------------------------------------------------------------------------------------------
        Vmath::Vmul(nPointsTot,m_eta[4],1,tmp3,1,tmp2,1); // reusing W^2
        Vmath::Vsub(nPointsTot,Fy,1,tmp2,1,Fy,1);

        Vmath::Vmul(nPointsTot,m_eta[3],1,tmp,1,tmp1,1);
        Vmath::Vsub(nPointsTot,Fy,1,tmp1,1,Fy,1);

        Vmath::Vsub(nPointsTot,Fy,1,m_eta[2],1,Fy,1);
        Vmath::Vmul(nPointsTot,m_eta[6],1,W,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Fy,1,tmp2,1,Fy,1);

        //-------------------------------------------------------------------------------------------------
        // z-component of the forcing - accelaration component
        //-------------------------------------------------------------------------------------------------
        Vmath::Vmul(nPointsTot,m_zeta[3],1,Px,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[3],1,Py,1,tmp2,1);
        Vmath::Vadd(nPointsTot,Fz,1,tmp1,1,Fz,1);
        Vmath::Vadd(nPointsTot,Fz,1,tmp2,1,Fz,1);

        //-------------------------------------------------------------------------------------------------
        // Note: Now we use Px,Py,Pz to store the viscous component of the forcing before we multiply
        // them by m_kinvis = 1/Re. Since we build up on them we need to set the entries to zero.
        //-------------------------------------------------------------------------------------------------
        Vmath::Zero(nPointsTot,Px,1);
        Vmath::Zero(nPointsTot,Py,1);
        Vmath::Zero(nPointsTot,Pz,1);

        //-------------------------------------------------------------------------------------------------
        // x-component of the forcing - viscous component1:  (U_z^'z^' - U_zz + ZETA_tz^'z^') 
        //-------------------------------------------------------------------------------------------------
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],fields[0]->GetPhys(),tmp1); // Uz
        fields[0]->HomogeneousBwdTrans(tmp1,tmp3); // Uz in physical space
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],tmp3,tmp1); // Uzx
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp3,tmp2); // Uzy
        
        Vmath::Vmul(nPointsTot,m_zeta[3],1,tmp1,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[3],1,tmp2,1,tmp2,1);
        Vmath::Smul(nPointsTot,2.0,tmp2,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp1,1,Px,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp2,1,Px,1);

        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],U,tmp1); // Ux
        Vmath::Vmul(nPointsTot,m_zeta[4],1,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp2,1,Px,1);
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],tmp1,tmp2); // Uxx
        Vmath::Vmul(nPointsTot,m_zeta[9],1,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp3,1,Px,1);
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp1,tmp2); // Uxy
        Vmath::Vmul(nPointsTot,m_zeta[8],1,tmp2,1,tmp3,1);
        Vmath::Smul(nPointsTot,2.0,tmp3,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp3,1,Px,1);

        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],U,tmp1); // Uy
        Vmath::Vmul(nPointsTot,m_eta[4],1,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp2,1,Px,1);
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp1,tmp2); // Uyy
        Vmath::Vmul(nPointsTot,m_eta[9],1,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp3,1,Px,1);
        Vmath::Vadd(nPointsTot,Px,1,m_zeta[7],1,Px,1);

        //------------------------------------------------------------------------------------------------
        // x-component of the forcing - viscous component2: ((ZETA_z * W)_z^'z^' + ZETA_z * (W_xx + W_yy))
        //------------------------------------------------------------------------------------------------
        Vmath::Vmul(nPointsTot,m_zeta[5],1,W,1,tmp1,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp1,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[3],1,Wx,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_zeta[4],1,tmp1,1,tmp2,1);
        Vmath::Smul(nPointsTot,3.0,tmp2,1,tmp3,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp3,1,Px,1);
        
        Vmath::Vmul(nPointsTot,m_eta[3],1,Wy,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_zeta[4],1,tmp1,1,tmp2,1);
        Vmath::Smul(nPointsTot,2.0,tmp2,1,tmp3,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp3,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[3],1,Wy,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[4],1,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp2,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[4],1,Wz,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp2,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[9],1,Wxy,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[3],1,tmp1,1,tmp2,1);
        Vmath::Smul(nPointsTot,2.0,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp3,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[9],1,Wxz,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp2,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[8],1,Wyz,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Px,1,tmp2,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[9],1,m_zeta[3],1,tmp1,1);
        Vmath::Vmul(nPointsTot,tmp1,1,Wxx,1,tmp2,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp2,1,Px,1);

        Vmath::Vmul(nPointsTot,m_zeta[3],1,Wyy,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[9],1,tmp1,1,tmp2,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp2,1,Px,1);

        Vmath::Vadd(nPointsTot,Wxx,1,Wyy,1,tmp1,1);
        Vmath::Vadd(nPointsTot,Wzz,1,tmp1,1,tmp2,1);
        Vmath::Vmul(nPointsTot,m_zeta[3],1,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Px,1,tmp3,1,Px,1);

        Vmath::Smul(nPointsTot,m_kinvis,Px,1,Px,1); //* 1/Re

        //-------------------------------------------------------------------------------------------------
        // y-component of the forcing - viscous component1: (V_z^'z^' - V_zz + ETA_tz^'z^')
        //-------------------------------------------------------------------------------------------------
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],fields[1]->GetPhys(),tmp1); // Vz
        fields[0]->HomogeneousBwdTrans(tmp1,tmp3); // Vz in physical space
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],tmp3,tmp1); // Vzx
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp3,tmp2); // Vzy
        Vmath::Vmul(nPointsTot,m_zeta[3],1,tmp1,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[3],1,tmp2,1,tmp2,1);
        Vmath::Smul(nPointsTot,2.0,tmp2,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp1,1,Py,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp2,1,Py,1);
       
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],V,tmp1); // Vx
        Vmath::Vmul(nPointsTot,m_zeta[4],1,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp2,1,Py,1);
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],tmp1,tmp2); // Vxx
        Vmath::Vmul(nPointsTot,m_zeta[9],1,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp3,1,Py,1);
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp1,tmp2); // Vxy
        Vmath::Vmul(nPointsTot,m_zeta[8],1,tmp2,1,tmp3,1);
        Vmath::Smul(nPointsTot,2.0,tmp3,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp3,1,Py,1);
       
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],V,tmp1); // Vy
        Vmath::Vmul(nPointsTot,m_eta[4],1,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp2,1,Py,1);
        fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],tmp1,tmp2); // Vyy
        Vmath::Vmul(nPointsTot,m_eta[9],1,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp3,1,Py,1);
        Vmath::Vadd(nPointsTot,Py,1,m_eta[7],1,Py,1);
       
        //------------------------------------------------------------------------------------------------
        // y-component of the forcing - viscous component2: ((ETA_z * W)_z^'z^' + ETA_z * (W_xx + W_yy))
        //------------------------------------------------------------------------------------------------
        Vmath::Vmul(nPointsTot,m_eta[5],1,W,1,tmp1,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp1,1,Py,1);
       
        Vmath::Vmul(nPointsTot,m_zeta[3],1,Wx,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[4],1,tmp1,1,tmp2,1);
        Vmath::Smul(nPointsTot,2.0,tmp2,1,tmp3,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp3,1,Py,1);

        Vmath::Vmul(nPointsTot,m_zeta[4],1,Wx,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[3],1,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp2,1,Py,1);

        Vmath::Vmul(nPointsTot,m_eta[3],1,Wy,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[4],1,tmp1,1,tmp2,1);
        Vmath::Smul(nPointsTot,3.0,tmp2,1,tmp3,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp3,1,Py,1);

        Vmath::Vmul(nPointsTot,m_eta[4],1,Wz,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp2,1,Py,1);

        Vmath::Vmul(nPointsTot,m_eta[9],1,Wxy,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_zeta[3],1,tmp1,1,tmp2,1);
        Vmath::Smul(nPointsTot,2.0,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp3,1,Py,1);

        Vmath::Vmul(nPointsTot,m_zeta[8],1,Wxz,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp2,1,Py,1);

        Vmath::Vmul(nPointsTot,m_eta[9],1,Wyz,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Py,1,tmp2,1,Py,1);

        Vmath::Vmul(nPointsTot,m_eta[3],1,Wxx,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_zeta[9],1,tmp1,1,tmp2,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp2,1,Py,1);

        Vmath::Vmul(nPointsTot,m_eta[9],1,m_eta[3],1,tmp1,1);
        Vmath::Vmul(nPointsTot,tmp1,1,Wyy,1,tmp2,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp2,1,Py,1);

        Vmath::Vadd(nPointsTot,Wxx,1,Wyy,1,tmp1,1);
        Vmath::Vadd(nPointsTot,Wzz,1,tmp1,1,tmp2,1);
        Vmath::Vmul(nPointsTot,m_eta[3],1,tmp2,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Py,1,tmp3,1,Py,1);
        
        Vmath::Smul(nPointsTot,m_kinvis,Py,1,Py,1); //* 1/Re

        //-------------------------------------------------------------------------------------------------
        // z-component of the forcing - viscous component: (W_z^'z^' - W_zz)
        //-------------------------------------------------------------------------------------------------
        Vmath::Vmul(nPointsTot,m_zeta[3],1,Wxz,1,tmp1,1);
        Vmath::Smul(nPointsTot,2.0,tmp1,1,tmp1,1);
        Vmath::Vmul(nPointsTot,m_eta[3],1,Wyz,1,tmp2,1);
        Vmath::Smul(nPointsTot,2.0,tmp2,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Pz,1,tmp1,1,Pz,1);
        Vmath::Vsub(nPointsTot,Pz,1,tmp2,1,Pz,1);

        Vmath::Vmul(nPointsTot,m_zeta[4],1,Wx,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Pz,1,tmp2,1,Pz,1);
        Vmath::Vmul(nPointsTot,m_zeta[9],1,Wxx,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Pz,1,tmp3,1,Pz,1);
        Vmath::Vmul(nPointsTot,m_zeta[8],1,Wxy,1,tmp3,1);
        Vmath::Smul(nPointsTot,2.0,tmp3,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Pz,1,tmp3,1,Pz,1);

        Vmath::Vmul(nPointsTot,m_eta[4],1,Wy,1,tmp2,1);
        Vmath::Vsub(nPointsTot,Pz,1,tmp2,1,Pz,1);
        Vmath::Vmul(nPointsTot,m_eta[9],1,Wyy,1,tmp3,1);
        Vmath::Vadd(nPointsTot,Pz,1,tmp3,1,Pz,1);

        Vmath::Smul(nPointsTot,m_kinvis,Pz,1,Pz,1); //* 1/Re

        //-------------------------------------------------------------------------------------------------
        // adding viscous and pressure components and transfroming back to wave space
        //-------------------------------------------------------------------------------------------------
        Vmath::Vadd(nPointsTot,Fx,1,Px,1,Fx,1);
        Vmath::Vadd(nPointsTot,Fy,1,Py,1,Fy,1);
        Vmath::Vadd(nPointsTot,Fz,1,Pz,1,Fz,1);
        fields[0]->HomogeneousFwdTrans(Fx,m_Forcing[0]);
        fields[0]->HomogeneousFwdTrans(Fy,m_Forcing[1]);
        fields[0]->HomogeneousFwdTrans(Fz,m_Forcing[2]);
    }
    
    /** 
     * Function to roll time-level storages to the next step layout.
     * The stored data associated with the oldest time-level 
     * (not required anymore) are moved to the top, where they will
     * be overwritten as the solution process progresses.
     */
    void ForcingMovingBody::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
    {
        int  nlevels = input.num_elements();

        Array<OneD, NekDouble> tmp;

        tmp = input[nlevels-1];
     
        for(int n = nlevels-1; n > 0; --n)
        {
            input[n] = input[n-1];
        }
    
        input[0] = tmp;
    }    

    void ForcingMovingBody::EvaluateAccelaration(const Array<OneD, NekDouble> &input, Array<OneD, NekDouble> &output, int npoints)
    {

        m_movingBodyCalls++;

        int acc_order = 0;

        // Rotate acceleration term
        RollOver(m_acceleration);

        Vmath::Vcopy(npoints, input, 1, m_acceleration[0], 1);

        //Calculate acceleration term at level n based on previous steps
        if (m_movingBodyCalls > 2)
        {
            acc_order = min(m_movingBodyCalls-2,m_intSteps);
            Vmath::Smul(npoints,
                       m_StifflyStable_Gamma0_Coeffs[acc_order-1],
                       m_acceleration[0], 1,
                       output,            1);
            for(int i = 0; i < acc_order; i++)
            {
                Vmath::Svtvp(npoints,
                            -1*m_StifflyStable_Alpha_Coeffs[acc_order-1][i],
                            m_acceleration[i+1], 1,
                            output,              1,
                            output,              1);
            }
            Vmath::Smul(npoints,
                       1.0/m_timestep,
                       output, 1,
                       output, 1);
        }
    }

    void ForcingMovingBody::CheckIsFromFile()
    {
        LibUtilities::FunctionType vType;
        
        // Check Displacement x
        vType = m_session->GetFunctionType(m_funcName[0],m_motion[0]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {m_IsFromFile[0] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {m_IsFromFile[0] = true;}
        else 
        {ASSERTL0(false,"The displacements in x must be specified via an equation <E> or a file <F>");}
        
        // Check Displacement y
        vType = m_session->GetFunctionType(m_funcName[0],m_motion[1]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {m_IsFromFile[1] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {m_IsFromFile[1] = true;}
        else 
        {ASSERTL0(false,"The displacements in y must be specified via an equation <E> or a file <F>");}
        
        // Check Velocity x
        vType = m_session->GetFunctionType(m_funcName[1],m_motion[0]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {m_IsFromFile[2] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {m_IsFromFile[2] = true;}
        else 
        {ASSERTL0(false,"The velocities in x must be specified via an equation <E> or a file <F>");}
        
        // Check Velocity y
        vType = m_session->GetFunctionType(m_funcName[1],m_motion[1]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {m_IsFromFile[3] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {m_IsFromFile[3] = true;}
        else 
        {ASSERTL0(false,"The velocities in y must be specified via an equation <E> or a file <F>");}
        
        // Check Acceleration x
        vType = m_session->GetFunctionType(m_funcName[2],m_motion[0]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {m_IsFromFile[4] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {m_IsFromFile[4] = true;}
        else 
        {ASSERTL0(false,"The accelerations in x must be specified via an equation <E> or a file <F>");}

        // Check Acceleration y
        vType = m_session->GetFunctionType(m_funcName[2],m_motion[1]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {m_IsFromFile[5] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {m_IsFromFile[5] = true;}
        else 
        {ASSERTL0(false,"The accelerations in y must be specified via an equation <E> or a file <F>");}
    }
 
}
}
