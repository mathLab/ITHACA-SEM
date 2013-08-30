///////////////////////////////////////////////////////////////////////////////
//
// File: Extrapolate.cpp
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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    ExtrapolateFactory& GetExtrapolateFactory()
    {
        typedef Loki::SingletonHolder<ExtrapolateFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy > Type;
        return Type::Instance();
    }

    Extrapolate::Extrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        Array<OneD, int> pVel)
        : m_session(pSession),
          m_fields(pFields),
          m_velocity(pVel),
          m_comm(pSession->GetComm())
    {

    }

    Extrapolate::~Extrapolate()
    {
    }

        
    std::string Extrapolate::def =
        LibUtilities::SessionReader::RegisterDefaultSolverInfo(
            "NoSubStepping", "No SubStepping");
	
	
	/**
	 * Curl Curl routine - dimension dependent
	 */
	void Extrapolate::CurlCurl(const MultiRegionss::ExpListSharedPtr &pField,
							   const Array<OneD, Array<OneD, NekDouble> > &Vel,
							   Array<OneD, Array<OneD, NekDouble> > &Q,
							   const int j)
	{
		m_elmt = pField->GetExp(m_HBCdata[j].m_globalElmtID);
		
		Array<OneD,NekDouble> Vx(m_pressureBCsMaxPts);
		Array<OneD,NekDouble> Uy(m_pressureBCsMaxPts);
		
		switch(pField->GetExpType())
		{
			case MultiRegions::e2D:
			{
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[1],Vx);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vel[2],Uy);  
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Uy,1,Dummy,1);
				
				m_elmt->PhysDeriv(Dummy,Q[1],Q[0]);
				
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,-1.0,Q[1],1,Q[1],1);
			}
				break;
				
			case MultiRegions::e3DH1D:
			{
				Array<OneD,NekDouble> Wz(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Dummy1(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> DUmmy2(m_pressureBCsMaxPts);
				
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[1],Vx);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vel[0],Uy);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_wavenumber[j],Vel[2],1,Wz,1);
				
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vx,Dummy1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Uy,Dummy2);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Dummy1,1,Dummy2,1,Q[0],1);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],Vel[0],1,Dummy1,1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Wz,Dummy2);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Q[0],1,Dummy1,1,Q[0],1);
				Vmath::Vadd(m_HBCdata[j].m_ptsInElmt,Q[0],1,Dummy2,1,Q[0],1);
				
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Wz,Dummy1);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],Vel[1],1,Dummy2,1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Dummy1,1,Dummy2,1,Q[1],1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vx,Dummy1);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Uy,Dummy2);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Q[1],1,Dummy1,1,Q[1],1);
				Vmath::Vadd(m_HBCdata[j].m_ptsInElmt,Q[1],1,Dummy2,1,Q[1],1);			
			}
				break;
				
			case MultiRegions::e3DH2D:
			{
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[2],Wx);
				m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[0],Vel[1],Vx);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],Vel[0],Uy);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],Vel[0],1,Uy,1);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[2],Vel[0],Uz);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_wavenumber[j],Vel[0],1,Uz,1);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Wz,1,Wx,1,qy,1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Uy,1,qz,1);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[1],qz,Uy);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],qz,1,Uy,1);
				
				//m_elmt->PhysDeriv(MultiRegionss::DirCartesianMap[2],qy,Uz);
				Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_wavenumber[j],qy,1,Uz,1);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Uy,1,Uz,1,Q[0],1);
			}
				break;
				
			case MultiRegions::e3D:
			{
				Array<OneD,NekDouble> Dummy(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Vz(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Uz(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Wx(m_pressureBCsMaxPts);
				Array<OneD,NekDouble> Wy(m_pressureBCsMaxPts);
				
				m_elmt->PhysDeriv(Vel[0],Dummy,Uy,Uz);
				m_elmt->PhysDeriv(Vel[1],Vx,Dummy,Vz);
				m_elmt->PhysDeriv(Vel[2],Wx,Wy,Dummy);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Wy,1,Vz,1,Q[0],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Uz,1,Wx,1,Q[1],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Uy,1,Q[2],1);
				
				m_elmt->PhysDeriv(Q[0],Dummy,Wy,Vx);
				m_elmt->PhysDeriv(Q[1],Wx,Dummy,Uz);
				m_elmt->PhysDeriv(Q[2],Vz,Uy,Dummy);
				
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Uy,1,Uz,1,Q[0],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vx,1,Vz,1,Q[1],1);
				Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Wx,1,Wy,1,Q[2],1);
			}
				break;
		}
	}
	
	
	/** 
	 * Function to roll time-level storages to the next step layout.
	 * The stored data associated with the oldest time-level 
	 * (not required anymore) are moved to the top, where they will
	 * be overwritten as the solution process progresses.
	 */
	void Extrapolate::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
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
	
	
	/**
	 * Map to directly locate HOPBCs position and offsets in all scenarios
	 */
	void Extrapolate::GenerateHOPBCMap(const MultiRegionss::ExpListSharedPtr &pField,
									const int HOPBCnumber)
	{
		m_HBCdata = Array<OneD, HBCInfo>(HOPBCnumber);
		
        switch(pField->GetExpType())
        {
			case MultiRegions::e2D:
			case MultiRegions::e3D:
			{
				int coeff_count = 0;
				int exp_size;
				int j=0;
				int cnt = 0;
				for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
				{
					exp_size = m_PBndExp[n]->GetExpSize();
                    
					if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
					{
						for(int i = 0; i < exp_size; ++i,cnt++)
						{
							m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];   
							m_elmt      = pField->GetExp(m_HBCdata[j].m_globalElmtID);
							m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();         
							m_HBCdata[j].m_physOffset = pField->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
							m_HBCdata[j].m_bndElmtOffset = i;       
							m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];      
							m_HBCdata[j].m_bndryElmtID = n;
							m_HBCdata[j].m_coeffOffset = coeff_count;
							coeff_count += m_elmt->GetEdgeNcoeffs(m_HBCdata[j].m_elmtTraceID);
							j = j+1;
						}
					}
					else // setting if just standard BC no High order
					{
						cnt += exp_size;
					}
				}
			}
				break;
				
			case MultiRegions::e3DH1D:
			{
				Array<OneD, unsigned int> planes;
				planes = pField->GetZIDs();
				int num_planes = planes.num_elements();            
				int num_elm_per_plane = (pField->GetExpSize())/num_planes;
				
				m_wavenumber      = Array<OneD, NekDouble>(HOPBCnumber);
				m_negWavenumberSq = Array<OneD, NekDouble>(HOPBCnumber);
				
				int coeff_count = 0;
				int exp_size, exp_size_per_plane;
				int j=0;
				int K;
				NekDouble sign = -1.0;
				int cnt = 0;
				
				for(int k = 0; k < num_planes; k++)
				{
					K = planes[k]/2;
					for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
					{
						exp_size = m_PBndExp[n]->GetExpSize();
						exp_size_per_plane = exp_size/num_planes;
						
						if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
						{
							for(int i = 0; i < exp_size_per_plane; ++i,cnt++)
							{
								m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];   
								m_elmt      = pField->GetExp(m_HBCdata[j].m_globalElmtID);
								m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();         
								m_HBCdata[j].m_physOffset = pField->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
								m_HBCdata[j].m_bndElmtOffset = i+k*exp_size_per_plane;       
								m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];      
								m_HBCdata[j].m_bndryElmtID = n;
								m_HBCdata[j].m_coeffOffset = coeff_count;
								coeff_count += m_elmt->GetEdgeNcoeffs(m_HBCdata[j].m_elmtTraceID);
								
								if(m_SingleMode)
								{
									m_wavenumber[j]      = -2*M_PI/m_LhomZ;       
									m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
								}
								else if(m_HalfMode || m_MultipleModes)
								{
									m_wavenumber[j]      = 2*M_PI/m_LhomZ;       
									m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
								}
								else
								{
									m_wavenumber[j]     = 2*M_PI*sign*(NekDouble(K))/m_LhomZ; 
									m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
								}
								
								int assElmtID;
								
								if(k%2==0)
								{
									if(m_HalfMode)
									{
										assElmtID = m_HBCdata[j].m_globalElmtID;
										
									}
									else
									{
										assElmtID = m_HBCdata[j].m_globalElmtID + num_elm_per_plane;
									}
								}
								else 
								{
									assElmtID = m_HBCdata[j].m_globalElmtID - num_elm_per_plane;
								}
								
								m_HBCdata[j].m_assPhysOffset = pField->GetPhys_Offset(assElmtID);
								
								j = j+1;
							}
						}
						else // setting if just standard BC no High order
						{
							cnt += exp_size_per_plane;
						}
					}
					sign = -1.0*sign;
				}
			}
				break;
				
			case MultiRegions::e3DH2D:
			{
				int cnt = 0;
				int exp_size, exp_size_per_line;
				int j=0;
				
				for(int k1 = 0; k1 < m_npointsZ; k1++)
				{
					for(int k2 = 0; k2 < m_npointsY; k2++)
					{
						for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
						{
							exp_size = m_PBndExp[n]->GetExpSize();
							
							exp_size_per_line = exp_size/(m_npointsZ*m_npointsY);
							
							if(m_PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
							{
								for(int i = 0; i < exp_size_per_line; ++i,cnt++)
								{
									// find element and edge of this expansion. 
									// calculate curl x curl v;
									m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];
									m_elmt      = pField->GetExp(m_HBCdata[j].m_globalElmtID);
									m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();
									m_HBCdata[j].m_physOffset = pField->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
									m_HBCdata[j].m_bndElmtOffset = i+(k1*m_npointsY+k2)*exp_size_per_line;
									m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];                
									m_HBCdata[j].m_bndryElmtID = n;
									m_wavenumber[j] = 2*M_PI*sign*(NekDouble(k1))/m_LhomZ;
									m_negWavenumberSq[j] = 2*M_PI*sign*(NekDouble(k2))/m_LhomY;
								}
							}
							else
							{
								cnt += exp_size_per_line;
							}
						}
					}
				}
			}
				break;
		}
    }
}

