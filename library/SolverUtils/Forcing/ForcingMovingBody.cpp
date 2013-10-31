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
        m_NumVariable = pNumForcingFields;
		m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        
        // Loading the x-dispalcement (m_zeta) and the y-displacement (m_eta)
		// Those two variables are bith functions of z and t and the may come
		// from an equation (forced vibration) or from another solver which, given
		// the aerodynamic forces at the previous step, calculates the displacements.
        
        //Get the body displacement: m_zeta and m_eta
        const TiXmlElement* funcNameElmt_D = pForce->FirstChildElement("DISPLACEMENTS");
        ASSERTL0(funcNameElmt_D,"MOVINGBODYFORCE tag has to specify a function name which " 
                                "prescribes the body displacement as d(z,t).");

        funcName[0] = funcNameElmt_D->GetText();
        ASSERTL0(m_session->DefinesFunction(funcName[0]),
                 "Function '" + funcName[0] + "' not defined.");
        
        //Get the body velocity of movement: d(m_zeta)/dt and d(m_eta)/dt
        const TiXmlElement* funcNameElmt_V = pForce->FirstChildElement("VELOCITIES");
        ASSERTL0(funcNameElmt_D,"MOVINGBODYFORCE tag has to specify a function name which " 
                 "prescribes the body velocity of movement as v(z,t).");
        
        funcName[1] = funcNameElmt_V->GetText();
        ASSERTL0(m_session->DefinesFunction(funcName[1]),
                 "Function '" + funcName[1] + "' not defined.");
        
        
        //Get the body acceleration: dd(m_zeta)/ddt and dd(m_eta)/ddt
        const TiXmlElement* funcNameElmt_A = pForce->FirstChildElement("ACCELERATIONS");
        ASSERTL0(funcNameElmt_A,"MOVINGBODYFORCE tag has to specify a function name which " 
                 "prescribes the body acceleration as a(z,t).");
        
        funcName[2] = funcNameElmt_A->GetText();
        ASSERTL0(m_session->DefinesFunction(funcName[2]),
                 "Function '" + funcName[2] + "' not defined.");
        
        // At this point we know in the xml file where those quantities
        // are declared (equation or file) - via a function name which is now stored in funcNameD etc.
        // We need now to fill in with this info the m_zeta and m_eta vectors (actuallythey are matrices)
        // Array to control if the motion is determined by an equation or is from a file.(not Nektar++)
        IsFromFile[3] = {false,false,false};
        // Defining motion directions
        motion[2] = {"x","y"};
        // check if we need to load a file or we have an equation
        CheckIsFromFile();
        
        int NumPoints = pFields->GetTotPoints();
        
        m_zeta = Array<OneD, Array< OneD, NekDouble> >(10);
        m_eta  = Array<OneD, Array< OneD, NekDouble> >(10);
        
        // What are this bi-dimensional vectors ----------------------------------------
        // m_zeta[0] = zeta                     |  m_eta[0] = eta                      |
        // m_zeta[1] = d(zeta)/dt               |  m_eta[1] = d(eta)/dt                |
        // m_zeta[2] = dd(zeta)/ddtt            |  m_eta[2] = dd(eta)/ddtt             |
        // m_zeta[3] = d(zeta)/dz               |  m_eta[3] = d(eta)/dz                |
        // m_zeta[4] = dd(zeta)/ddzz            |  m_eta[4] = dd(eta)/ddzz             |
        // m_zeta[5] = ddd(zeta)/dddzzz         |  m_eta[5] = ddd(eta)/dddzzz          |
        // m_zeta[6] = ddd(zeta)/dddtzz         |  m_eta[6] = ddd(eta)/dddtzz          |
        // m_zeta[7] = (d(zeta)/dz)^2           |  m_eta[7] = (d(eta)/dz)^2            |
        // m_zeta[8] = (d(zeta)/dz)(d(eta)/dz)  |  m_eta[8] = (d(eta)/dz)(d(zeta)/dz)  |
        // m_zeta[9] = dd(zeta)/ddtz            |  m_eta[9] = dd(eta)/ddtz             |
        //------------------------------------------------------------------------------
        
        for(int i = 0; i < m_zeta.num_elements(); i++)
        {
            m_zeta[i] = Array<OneD, NekDouble>(phystot,0.0);
            m_eta[i]  = Array<OneD, NekDouble>(phystot,0.0);
        }
    }

    void ForcingMovingBody::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
                  Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
		// Fill in m_zeta and m_eta with all the required values
        UpdateMotion(fields);
		
		//calcualte the forcing components Ax,Ay,Az
		CalculateForcing(fields);
		
		// Apply forcing terms
        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }
    }
				 

				 
	void ForcingMovingBody::UpdateMotion()
	{
					 
	}
	
	void ForcingMovingBody::CalculateForcing(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
	{
		
	}
        
    void ForcingMovingBody::CheckIsFromFile()
    {
        // Note: we just check the first component - 'x' (m_zeta) - of each time derivative
        LibUtilities::FunctionType vType;
        
        // Check Displacement
        vType = m_session->GetFunctionType(funcName[0],motion[0]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {IsFromFile[0] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {IsFromFile[0] = true;}
        else 
        {ASSERTL0(false,"The displacements must be specified via and equation <E> or a file <F>");}
        
        // Check Velocity
        vType = m_session->GetFunctionType(funcName[1],motion[0]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {IsFromFile[0] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {IsFromFile[0] = true;}
        else 
        {ASSERTL0(false,"The velocities must be specified via and equation <E> or a file <F>");}
        
        // Check Acceleration
        vType = m_session->GetFunctionType(funcName[2],motion[0]);
        if(vType == LibUtilities::eFunctionTypeExpression)
        {IsFromFile[0] = false;}
        else if(vType == LibUtilities::eFunctionTypeFile)
        {IsFromFile[0] = true;}
        else 
        {ASSERTL0(false,"The accelerations must be specified via and equation <E> or a file <F>");}
    }

}
}
