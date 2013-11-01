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
                                                        ForcingBody::create,
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
		// forcing size (it must be 3)
        m_NumVariable  = pNumForcingFields;
		m_Forcing      = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        m_funcName     = Array<OneD, std::string> (3);
        m_motion       = Array<OneD, std::string> (2);
        m_motion[0]    = "x";
        m_motion[1]    = "y";
        m_IsFromFile   = Array<OneD, bool> (6);
        m_session->LoadParameter("Kinvis",m_kinvis);
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
        int NumPoints = pFields[0]->GetTotPoints();
        
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
	
	void ForcingMovingBody::CalculateForcing(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                             const Array<OneD, Array<OneD, NekDouble> > &inarray)
	{
        //////////////// Work in progress
		
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
