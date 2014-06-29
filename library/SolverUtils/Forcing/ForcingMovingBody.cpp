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
        
        Array<OneD, unsigned int> ZIDs;
        ZIDs = pFields[0]->GetZIDs();
        m_NumLocPlanes    = ZIDs.num_elements();
        m_NumVariable     = pNumForcingFields; // forcing size (it must be 3)
		m_Forcing         = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        m_funcName        = Array<OneD, std::string> (3);
        m_motion          = Array<OneD, std::string> (2);
        m_motion[0]       = "x";
        m_motion[1]       = "y";
        m_movingBodyCalls = 0;
        m_IsFromFile      = Array<OneD, bool> (6);
		m_NumStructVars   = 3; //variables number must be 3 including displacement, velocity and acceleration
        m_session->LoadParameter("LZ", m_lhomz);        
        m_session->LoadParameter("Kinvis",m_kinvis);
        m_session->LoadParameter("TimeStep", m_timestep, 0.01); 

		//Loading strcutural parameters from input file
        m_session->LoadParameter("StructRho",m_structrho); 
        m_session->LoadParameter("StructDamp",m_structdamp, 0.0);
        m_session->LoadParameter("StructStiff",m_structstiff, 0.0); 
		m_session->LoadParameter("CableTension",m_cabletension, 0.0);
		m_session->LoadParameter("BendingStiff",m_bendingstiff, 0.0);

        m_Motion_x = Array<OneD, Array<OneD, NekDouble> > (m_NumStructVars);
        m_Motion_y = Array<OneD, Array<OneD, NekDouble> > (m_NumStructVars);
        for(int i = 0; i < m_Motion_x.num_elements(); i++)
        {
            m_Motion_x[i] = Array<OneD, NekDouble>(m_NumLocPlanes,0.0);
            m_Motion_y[i] = Array<OneD, NekDouble>(m_NumLocPlanes,0.0);
        }

		//Setting the coefficient matrices for solving structural dynamic ODEs
		SetDynEqCoeffMatrix(pFields);

        //Setting the finite difference coefficients used to calculate accelerration term dw/dt
		SetBWFinitDiffCoeffs();
		
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
       
		// other choises of higher order would make the algorithm to become unstable, 
		// work on this issue is in progress 

        m_acceleration = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);
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

        // Set initial condition for cable motion
        int nzpoints = pFields[0]->GetHomogeneousBasis()->GetNumModes();
        Array<OneD, NekDouble> z_coords(nzpoints,0.0);
        Array<OneD, const NekDouble> pts = pFields[0]->GetHomogeneousBasis()->GetZ();

        Vmath::Smul(nzpoints,m_lhomz/2.0,pts,1,z_coords,1);
       	Vmath::Sadd(nzpoints,m_lhomz/2.0,z_coords,1,z_coords,1);
        Array<OneD, NekDouble> x0(m_NumLocPlanes, 0.0);
        Array<OneD, NekDouble> x1(m_NumLocPlanes, 0.0);
        Array<OneD, NekDouble> x2(m_NumLocPlanes, 0.0);

		for (int plane = 0; plane < m_NumLocPlanes; plane++)
		{
        	x2[plane] = z_coords[ZIDs[plane]];
		}

        for(int j = 0; j < m_funcName.num_elements(); j++)
        {         
            LibUtilities::EquationSharedPtr ffunc0,ffunc1;
            ffunc0 = m_session->GetFunction(m_funcName[j], m_motion[0]);
			ffunc1 = m_session->GetFunction(m_funcName[j], m_motion[1]);
            ffunc0 ->Evaluate(x0, x1, x2, 0.0, m_Motion_x[j]);
            ffunc1 ->Evaluate(x0, x1, x2, 0.0, m_Motion_y[j]);
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
		
        // For fluid structure interaction case, displacements, velocities and acceleartions
        // are obtained through solving structure dynamic model
        EvaluateStructDynModel(pFields, time);
	
		/*// For forced vibration case, displacements, velocities, acceleration 
		// are loading from specified file or function	
		int cnt = 0;
		for(int j = 0; j < m_funcName.num_elements(); j++)
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
		}*/
		
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
    
   	void ForcingMovingBody::TensionedCableModel(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                        							 Array<OneD, NekDouble> &Force,
                        							 Array<OneD, Array<OneD, NekDouble> > &Motion)
    {
    	//
    	Array<OneD, NekDouble> ForceCoeffs;
        ForceCoeffs = Array<OneD, NekDouble> (m_NumLocPlanes);
    	pFields[0]->HomogeneousFwdTrans(Force, ForceCoeffs);

        Array<OneD, Array<OneD, NekDouble> > MotionCoeffs;
        MotionCoeffs = Array<OneD, Array<OneD, NekDouble> > (m_NumStructVars);
        for (int var = 0; var < m_NumStructVars; var++)
        {
            MotionCoeffs[var] = Array<OneD, NekDouble> (m_NumLocPlanes,0.0);
        	pFields[0]->HomogeneousFwdTrans(Motion[var], MotionCoeffs[var]);
        }
    	// solve the ODE in the wave space for cable motion to obtain disp, vel and accel
    	for(int plane = 0; plane < m_NumLocPlanes; plane++)
    	{
        	int nrows = m_NumStructVars;
            Array<OneD, NekDouble> tmp1,tmp2;
            tmp1 = Array<OneD, NekDouble> (m_NumStructVars,0.0);
            tmp2 = Array<OneD, NekDouble> (m_NumStructVars,0.0);

            for(int var = 0; var < m_NumStructVars; var++)
            {
                tmp1[var] = MotionCoeffs[var][plane];
            }
        	tmp2[0] = ForceCoeffs[plane];

        	Blas::Dgemv('N', nrows, nrows, 1.0, 
						&(m_CoeffMat_B[plane]->GetPtr())[0],
                        nrows, &tmp1[0], 1, 
						0.0,   &tmp1[0], 1);
            Blas::Dgemv('N', nrows, nrows, 1.0/m_structrho, 
						&(m_CoeffMat_A[plane]->GetPtr())[0],
                        nrows, &tmp2[0], 1, 
						1.0,   &tmp1[0], 1);

            for(int var = 0; var < m_NumStructVars; var++)
            {
            	MotionCoeffs[var][plane] = tmp1[var];
        	}
    	}

        for(int var = 0; var < m_NumStructVars; var++)
        {
        	pFields[0]->HomogeneousBwdTrans(MotionCoeffs[var], Motion[var]);
    	}
    }

    void ForcingMovingBody::EvaluateStructDynModel(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                                                   const NekDouble &time )
    {
        //Get the hydrodynamic forces from the fluid solver
		Array<OneD, unsigned int> ZIDs;
        ZIDs = pFields[0]->GetZIDs();
		Array <OneD, NekDouble> tmp;
	
		int ndirmotion = m_motion.num_elements();
        Array<OneD, Array <OneD, NekDouble> > Hydrofces(ndirmotion);
		Hydrofces[0] = Array <OneD, NekDouble> (m_NumLocPlanes);
        Hydrofces[1] = Array <OneD, NekDouble> (m_NumLocPlanes);
        for(int plane = 0; plane < m_NumLocPlanes; plane++)
        {
            pFields[0]->GetPlane(plane)->GetMovBodyForces(tmp);
		 	Hydrofces[0][plane] = tmp[0];
            Hydrofces[1][plane] = tmp[1];
        }
		
        TensionedCableModel(pFields, Hydrofces[0], m_Motion_x);
        TensionedCableModel(pFields, Hydrofces[1], m_Motion_y);
		
		//Set the forcing term based on the motion of the cable
        for(int plane = 0 ; plane < m_NumLocPlanes; plane++)
        {
            int NumPhyPoints = pFields[0]->GetPlane(plane)->GetTotPoints();
            for(int var = 0; var < m_NumStructVars; var++)
            {
                for(int point = 0; point < NumPhyPoints; point++)
                {
                    int j;
                    j = point + plane * NumPhyPoints;
                    m_zeta[var][j] = m_Motion_x[var][plane];
                    m_eta[var][j]  = m_Motion_y[var][plane];
                }
            }
        }
	
		
        for(int plane = 0; plane < m_NumLocPlanes; plane++)
        {
			Array<OneD, NekDouble> tmp0(3);
			Array<OneD, NekDouble> tmp1(3);
			for (int var = 0; var < m_NumStructVars; var++)
			{
				tmp0[var] = m_Motion_x[var][plane];
				tmp1[var] = m_Motion_y[var][plane];
			}
            pFields[0]->GetPlane(plane)->SetMovBodyMotionVars(tmp0);
            pFields[1]->GetPlane(plane)->SetMovBodyMotionVars(tmp1);
        }
	
        // velocity of moving body is used to set the boundary conditions
		// according to the coordinate system mapping
        Array<OneD, Array<OneD, const MultiRegions::ExpListSharedPtr> > bndCondExps;
        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> > bndConds;
        bndCondExps = Array<OneD, Array<OneD, const MultiRegions::ExpListSharedPtr> >(ndirmotion);
        bndConds = Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >(ndirmotion);
        for (int plane = 0; plane < m_NumLocPlanes; plane++)
        {
            Array<OneD, Array<OneD, NekDouble> > tmp(2);
			tmp[0] = Array<OneD, NekDouble>(3,0.0);
			tmp[1] = Array<OneD, NekDouble>(3,0.0);
			tmp[0][1] = m_Motion_x[1][plane];
			tmp[1][1] = m_Motion_y[1][plane];
			 
            for (int dir = 0; dir < ndirmotion; dir++)
            {
                bndCondExps[dir] = pFields[dir]->GetPlane(plane)->GetBndCondExpansions();
                bndConds[dir] = pFields[dir]->GetPlane(plane)->GetBndConditions();
                int nbnd = bndCondExps[dir].num_elements();
                for (int i = 0; i < nbnd; ++i)
                {
                    if (bndConds[dir][i]->GetUserDefined() ==
                        SpatialDomains::eFluidStructInt &&
                        bndConds[dir][i]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                    {
						bndCondExps[dir][i]->SetMovBodyMotionVars(tmp[dir]);
                    }
                }
            }
        }
    }

    void ForcingMovingBody::SetDynEqCoeffMatrix(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
    {
        NekDouble tmp1, tmp2, tmp3;
        NekDouble tmp4, tmp5, tmp6, tmp7;
        NekDouble beta =  2.0 * M_PI/m_lhomz;
        tmp1 =     m_timestep * m_timestep;
        tmp2 =  m_structstiff / m_structrho;
        tmp3 =   m_structdamp / m_structrho;
        tmp4 = m_cabletension / m_structrho;
        tmp5 = m_bendingstiff / m_structrho;

        m_CoeffMat_A = Array<OneD, DNekMatSharedPtr>(m_NumLocPlanes);
        m_CoeffMat_B = Array<OneD, DNekMatSharedPtr>(m_NumLocPlanes);

        // solve the ODE in the wave space for cable motion to obtain disp, vel and accel
        for(int plane = 0; plane < m_NumLocPlanes; plane++)
        {
            int nel = m_NumStructVars;
            m_CoeffMat_A[plane] = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);
            m_CoeffMat_B[plane] = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);

            Array<OneD, unsigned int> ZIDs;
            ZIDs = pFields[0]->GetZIDs();
            unsigned int n_mode = ZIDs[plane]/2;
            tmp6 = beta * boost::lexical_cast<NekDouble>(n_mode);
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

    void ForcingMovingBody::SetBWFinitDiffCoeffs()
    {
        // backward finite difference coefficient for calculation
        // of spanwise acceleration dw/dt
        m_Alpha       = Array<OneD, Array<OneD, NekDouble> > (3);
        m_Alpha[0]    = Array<OneD, NekDouble> (3);
        m_Alpha[1]    = Array<OneD, NekDouble> (3);
        m_Alpha[2]    = Array<OneD, NekDouble> (3);
        m_Gamma0      = Array<OneD, NekDouble> (3);
        m_Alpha[0][0] = 1.0;
        m_Alpha[0][1] = 0.0;
        m_Alpha[0][2] = 0.0;
        m_Alpha[1][0] = 2.0;
        m_Alpha[1][1] =-0.5;
        m_Alpha[1][2] = 0.0;
        m_Alpha[2][0] = 3.0;
        m_Alpha[2][1] =-1.5;
        m_Alpha[2][2] = 1.0/3.0;
        m_Gamma0[0]   = 1.0;
        m_Gamma0[1]   = 1.5;
        m_Gamma0[2]   = 11.0/6.0;
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

    void ForcingMovingBody::EvaluateAccelaration(const Array<OneD, NekDouble> &input,
                                                       Array<OneD, NekDouble> &output,
                                                   int npoints)
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
                        m_Gamma0[acc_order-1],
                        m_acceleration[0], 1,
                        output,            1);
            for(int i = 0; i < acc_order; i++)
            {
                Vmath::Svtvp(npoints,
                             -1*m_Alpha[acc_order-1][i],
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
