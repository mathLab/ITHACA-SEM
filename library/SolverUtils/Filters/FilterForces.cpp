///////////////////////////////////////////////////////////////////////////////
//
// File FilterForces.cpp
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
// Description: Output values of aerodynamic forces during time-stepping
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <SolverUtils/Filters/FilterForces.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterForces::className = GetFilterFactory().RegisterCreatorFunction("Forces", FilterForces::create);

        /**
         *
         */
        FilterForces::FilterForces(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            if (pParams.find("OutputFile") == pParams.end())
            {
                m_outputFile = m_session->GetSessionName();
            }
            else
            {
                ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                         "Missing parameter 'OutputFile'.");
                m_outputFile = pParams.find("OutputFile")->second;
            }
            if (!(m_outputFile.length() >= 4 
                  && m_outputFile.substr(m_outputFile.length() - 4) == ".frc"))
            {
                m_outputFile += ".frc";
            }

            if (pParams.find("OutputFrequency") == pParams.end())
            {
                m_outputFrequency = 1;
            }
            else
            {
                m_outputFrequency = atoi(pParams.find("OutputFrequency")->second.c_str());
            }


            m_session->MatchSolverInfo("Homogeneous","1D",m_isHomogeneous1D,false);
            
            if(m_isHomogeneous1D)
            {
                if (pParams.find("OutputPlane") == pParams.end())
                {
                    m_outputPlane = 0;
                }
                else
                {
                    m_outputPlane = atoi(pParams.find("OutputPlane")->second.c_str());
                }
            }

        }


        /**
         *
         */
        FilterForces::~FilterForces()
        {

        }


        /**
         *
         */
        void FilterForces::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
			LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
			
			/*
			int n, cnt, elmtid, nq, offset, nt, boundary, nfq;
			nt = pFields[0]->GetNpoints();
			int dim = pFields[0]->GetGraph()->GetSpaceDimension();
            LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
            int vRank = vComm->GetRank();
			
			
			// Set up mapping from Boundary condition to element details.
			StdRegions::StdExpansionSharedPtr elmt;
			StdRegions::StdExpansion1DSharedPtr bc;
			Array<OneD, int> BoundarytoElmtID;
			Array<OneD, int> BoundarytoTraceID;
			Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
			Array<OneD, const NekDouble> P(nt);
			Array<OneD, NekDouble> outvalues;
			NekDouble D,L;
			
			//initialisation of drag and lift
			D=0.0;
			L=0.0;
			
			//boundary map and BndExp
			pFields[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);
			BndExp = pFields[0]->GetBndCondExpansions();
			

			if(m_isHomogeneous1D)
			{
				
			m_homogeneousBasis=pFields[0]->GetHomogeneousBasis();				
			Array<OneD, const NekDouble> w = m_homogeneousBasis->GetW();
			int num_planes=m_session->GetParameter("HomModesZ");

			cout << "plot " << pFields[2]->GetPlane(0)->GetPhys()[22] << endl;
				cout<< "exp size " << BndExp[0]->GetPlane(0)->GetExpSize() << endl;
			
			// loop over the types of boundary conditions
			for(cnt = n = 0; n < BndExp.num_elements(); ++n)
			{ 
				cout << "Bnd Expansion n. " << n << endl;
				//bluff body condition
				if((pFields[0]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet && 
					pFields[0]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eBluffBody)&&	
				   (pFields[1]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet &&
					pFields[1]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eBluffBody)&&
				   (pFields[2]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet &&
					pFields[2]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eBluffBody)&&
				   (pFields[3]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eNeumann && 
					pFields[3]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eHigh))
				{
					for(int i = 0; i < BndExp[n]->GetPlane(0)->GetExpSize(); ++i, cnt++)
					{
						
						cout << "Inside " << i << endl;
						// find element and face of this expansion.
						elmtid = BoundarytoElmtID[cnt];
						cout << "elmID" << endl;
						elmt   = pFields[0]->GetPlane(0)->GetExp(elmtid);
						cout << "GetExp(elmID)" << endl;
						nq     = elmt->GetTotPoints();
						offset = pFields[0]->GetPlane(0)->GetPhys_Offset(elmtid);
						cout << "GetPhys_offest" << endl;

						//identify boundary of element
						boundary = BoundarytoTraceID[cnt];
						//Extract pressure field
						P = pFields[2]->GetPlane(0)->GetPhys() + offset;
						
						cout << "p calculated" << endl;
						// Get face 2D expansion from element expansion
						bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(i));
						
						cout << "bc done" << endl;
						//number of points on the boundary 
						int nbc = bc->GetTotPoints();
						Array<OneD, NekDouble> Pb(nbc);
						Array<OneD, NekDouble>  drag(nbc);
						Array<OneD, NekDouble>  lift(nbc);

						
						//identify boundary of element looking at.
						boundary = BoundarytoTraceID[cnt];
						
						elmt->GetEdgePhysVals(boundary,bc,P,Pb);
												
						const Array<OneD, Array<OneD, NekDouble> > &normals
						= elmt->GetEdgeNormal(boundary);
						
						Vmath::Vmul(nbc,Pb,1,normals[0],1,drag,1);
						Vmath::Vmul(nbc,Pb,1,normals[1],1,lift,1);
						
						D=D+bc->Integral(drag);
						L=L+bc->Integral(lift);
						cout << "Done" << endl;
					}
				
										
				}
				else 
				{
					cout << "advance " << endl;
					cnt += BndExp[n]->GetExpSize();
				}
	
			}
*/
			if (vComm->GetRank() == 0)
			{


			// Open output stream
			m_outputStream.open(m_outputFile.c_str());
			m_outputStream << "# Time, drag and lift:" ;
			m_outputStream << endl;

				
			/*	
			m_outputStream.width(8);
			m_outputStream << setprecision(6) << time;
			
			m_outputStream.width(25);
			m_outputStream << setprecision(16) << D;
			
			m_outputStream.width(25);
			m_outputStream << setprecision(16) << L;
			
			m_outputStream << endl;
			
			}*/
			}
		
            v_Update(pFields, time);
        }


        /**
         *
         */
        void FilterForces::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            // Only output every m_outputFrequency.
            if ((m_index++) % m_outputFrequency)
            {
                return;
            }

			int n, cnt, elmtid, nq, offset, nt, boundary;
			nt = pFields[0]->GetNpoints();
			int dim = pFields.num_elements()-1;
			
			// Set up mapping from Boundary condition to element details.
			StdRegions::StdExpansionSharedPtr elmt;
			StdRegions::StdExpansion1DSharedPtr bc;
			Array<OneD, int> BoundarytoElmtID;
			Array<OneD, int> BoundarytoTraceID;
			Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
	
			Array<OneD, const NekDouble> P(nt);
			Array<OneD, const NekDouble> U(nt);
			Array<OneD, const NekDouble> V(nt);
			Array<OneD, const NekDouble> W(nt);

			
			Array<OneD, Array<OneD, NekDouble> > gradU(dim);
			Array<OneD, Array<OneD, NekDouble> > gradV(dim);
			Array<OneD, Array<OneD, NekDouble> > gradW(dim);
			
			Array<OneD, Array<OneD, NekDouble> > fgradU(dim);
			Array<OneD, Array<OneD, NekDouble> > fgradV(dim);
			Array<OneD, Array<OneD, NekDouble> > fgradW(dim);


			/*
			Array<OneD, NekDouble > Ux(nt);
			Array<OneD, NekDouble > Uy(nt);
			Array<OneD, NekDouble > Uz(nt);

			Array<OneD, NekDouble > Vx(nt);
			Array<OneD, NekDouble > Vy(nt);
			Array<OneD, NekDouble > Vz(nt);

			
			Array<OneD, NekDouble > Wx(nt);
			Array<OneD, NekDouble > Wy(nt);
			Array<OneD, NekDouble > Wz(nt);
			*/

			Array<OneD, NekDouble> outvalues;
			LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

			NekDouble D,L;
			D=0.0;
			L=0.0;
			
			
			NekDouble rho=(m_session->DefinesParameter("rho")) ? (m_session->GetParameter("rho")):1;
			NekDouble kinvis=rho*m_session->GetParameter("Kinvis");

			if(m_isHomogeneous1D)
			{
				pFields[0]->GetPlane(0)->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);
				BndExp = pFields[0]->GetPlane(0)->GetBndCondExpansions();
				
				m_homogeneousBasis=pFields[0]->GetHomogeneousBasis();				
				Array<OneD, const NekDouble> w = m_homogeneousBasis->GetW();
		
			// loop over the types of boundary conditions
			for(cnt = n = 0; n < BndExp.num_elements(); ++n)
			{ 
				//condition for selecting the boundary
				if((pFields[0]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet && 
					pFields[0]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eWall_Forces)&&	
				   (pFields[1]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet &&
					pFields[1]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eWall_Forces)&&
				   (pFields[2]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet &&
					pFields[2]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eWall_Forces)&&
				   (pFields[3]->GetPlane(0)->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eNeumann && 
					pFields[3]->GetPlane(0)->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eHigh))
				{
					cout << "Inside" << endl;

					for(int i = 0; i <  BndExp[n]->GetExpSize(); ++i, cnt++)
					{
						// find element and face of this expansion.
						elmtid = BoundarytoElmtID[cnt];
						elmt   = pFields[0]->GetPlane(0)->GetExp(elmtid);
						nq     = elmt->GetTotPoints();
						offset = pFields[0]->GetPlane(0)->GetPhys_Offset(elmtid);
						
						// Initialise local arrays for the velocity gradients
						// size of total number of quadrature points for each element (hence local).
						for(int j = 0; j < dim; ++j)
						{
							gradU[j] = Array<OneD, NekDouble>(nq);
							gradV[j] = Array<OneD, NekDouble>(nq);
							gradW[j] = Array<OneD, NekDouble>(nq);

						}
						
						
						//identify boundary of element
						boundary = BoundarytoTraceID[cnt];
						
						//Extract  field
						U = pFields[0]->GetPlane(0)->GetPhys() + offset;
						V = pFields[1]->GetPlane(0)->GetPhys() + offset;
						W = pFields[2]->GetPlane(0)->GetPhys() + offset;
						P = pFields[3]->GetPlane(0)->GetPhys() + offset;

						elmt->PhysDeriv(U,gradU[0],gradU[1],gradU[2]);
						elmt->PhysDeriv(V,gradV[0],gradV[1],gradV[2]);
						elmt->PhysDeriv(W,gradW[0],gradW[1],gradW[2]);

						// Get face 2D expansion from element expansion
						bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(i));
						
						//number of points on the boundary 
						int nbc = bc->GetTotPoints();
						Array<OneD, NekDouble> Pb(nbc);
						Array<OneD, NekDouble> Ub(nbc);

						Array<OneD, NekDouble>  drag(nbc);
						Array<OneD, NekDouble>  lift(nbc);
						Array<OneD, NekDouble>  temp(nbc);
						Array<OneD, NekDouble>  temp2(nbc);

						//identify boundary of element looking at.
						boundary = BoundarytoTraceID[cnt];
						
						elmt->GetEdgePhysVals(boundary,bc,P,Pb);
						elmt->GetEdgePhysVals(boundary,bc,U,Ub);

						
						for(int j = 0; j < dim; ++j)
						{
							fgradU[j] = Array<OneD, NekDouble>(nbc);
							fgradV[j] = Array<OneD, NekDouble>(nbc);
							fgradW[j] = Array<OneD, NekDouble>(nbc);

							
						} 
						
						for(int j = 0; j < dim; ++j)
						{
							elmt->GetEdgePhysVals(boundary,bc,gradU[j],fgradU[j]);
						    elmt->GetEdgePhysVals(boundary,bc,gradV[j],fgradV[j]);
							elmt->GetEdgePhysVals(boundary,bc,gradW[j],fgradW[j]);

						}
						
						const Array<OneD, Array<OneD, NekDouble> > &normals
						= elmt->GetEdgeNormal(boundary);
						
						
						//tangential viscous tractive forces
						
						//====drag terms=========
						//-(du/dz+dw/dx)*nz (nz=1)
						Vmath::Vadd(nbc,fgradU[2],1,fgradW[0],1,temp,1);
						Vmath::Neg(nbc,temp,1);
						//-(du/dy+dv/dx)*ny
						Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,drag,1);
						Vmath::Neg(nbc,drag,1);
						Vmath::Vmul(nbc,drag,1,normals[1],1,drag,1);
						//-2*du/dx*nx-(du/dy+dv/dx)*ny
						Vmath::Smul(nbc,-2.0,fgradU[0],1,fgradU[0],1);
						Vmath::Vmul(nbc,fgradU[0],1,normals[0],1,temp2,1);
						
						//take back gradU as before
						Vmath::Smul(nbc,-0.5,fgradU[0],1,fgradU[0],1);

						Vmath::Vadd(nbc,temp,1,temp2,1,temp,1);
						Vmath::Vadd(nbc,temp,1,drag,1,drag,1);
						//multiply by viscosity
						Vmath::Smul(nbc,kinvis,drag,1,drag,1);
						//zero temporary storage vector
						Vmath::Zero(nbc,temp,0);
						Vmath::Zero(nbc,temp2,0);

						
						//========lift terms==============
						//-(dv/dz+dw/dy)*ny
						Vmath::Vadd(nbc,fgradV[2],1,fgradW[1],1,temp,1);
						Vmath::Neg(nbc,temp,1);
						//-(du/dy+dv/dx)*nx
						Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,lift,1);
						Vmath::Neg(nbc,lift,1);
						Vmath::Vmul(nbc,lift,1,normals[0],1,lift,1);
						//-2*dv/dy*nx-(du/dy+dv/dx)
						Vmath::Smul(nbc,-2.0,fgradV[1],1,fgradV[1],1);
						Vmath::Vmul(nbc,fgradV[1],1,normals[1],1,temp2,1);
						
						Vmath::Smul(nbc,-0.5,fgradV[1],1,fgradV[1],1);

						
						Vmath::Vadd(nbc,temp,1,temp2,1,temp,1);
						Vmath::Vadd(nbc,temp,1,lift,1,lift,1);
						//multiply by viscosity
						Vmath::Smul(nbc,kinvis,lift,1,lift,1);
						
						//normal tractive forces added
						Vmath::Vvtvp(nbc,Pb,1,normals[0],1,drag,1,drag,1);
						Vmath::Vvtvp(nbc,Pb,1,normals[1],1,lift,1,lift,1);
						
						D=D+bc->Integral(drag);
						L=L+bc->Integral(lift);
						
						cout << "L " << L << endl;
						cout << "D " << D << endl;


					}
					
				}
				else 
				{
					
					cout << "not bb" << endl;
					cnt += BndExp[n]->GetExpSize();
				}
			}
					
		} //m_homogeneous
		else {
			//2D case
			cout << "2D case" << endl;
			pFields[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);
			BndExp = pFields[0]->GetBndCondExpansions();

			// loop over the types of boundary conditions
			for(cnt = n = 0; n < BndExp.num_elements(); ++n)
			{ 
				
				if((pFields[0]->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet && 
					pFields[0]->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eWall_Forces)&&	
				   (pFields[1]->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eDirichlet &&
					pFields[1]->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eWall_Forces)&&
				   (pFields[2]->GetBndConditions()[n]->GetBoundaryConditionType()==SpatialDomains::eNeumann && 
					pFields[2]->GetBndConditions()[n]->GetUserDefined()==SpatialDomains::eHigh))	
				{
					for(int i = 0; i <  BndExp[n]->GetExpSize(); ++i, cnt++)
					{
						
						// find element and face of this expansion.
						elmtid = BoundarytoElmtID[cnt];
						elmt   = pFields[0]->GetExp(elmtid);
						nq     = elmt->GetTotPoints();
						offset = pFields[0]->GetPhys_Offset(elmtid);
						
						// Initialise local arrays for the velocity gradients
						// size of total number of quadrature points for each element (hence local).
						for(int j = 0; j < dim; ++j)
						{
							gradU[j] = Array<OneD, NekDouble>(nq);
							gradV[j] = Array<OneD, NekDouble>(nq);

						}
						
						//identify boundary of element
						boundary = BoundarytoTraceID[cnt];
						
						//Extract fields
						U = pFields[0]->GetPhys() + offset;
						V = pFields[1]->GetPhys() + offset;
						P = pFields[2]->GetPhys() + offset;

						elmt->PhysDeriv(U,gradU[0],gradU[1]);
						elmt->PhysDeriv(V,gradV[0],gradV[1]);
						
						
						// Get face 2D expansion from element expansion
						bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (BndExp[n]->GetExp(i));
						

						
						//number of points on the boundary 
						int nbc = bc->GetTotPoints();
						Array<OneD, NekDouble> Pb(nbc);
						
						Array<OneD, NekDouble>  drag(nbc);
						Array<OneD, NekDouble>  lift(nbc);
						Array<OneD, NekDouble>  temp(nbc);

						
						//identify boundary of element looking at.
						boundary = BoundarytoTraceID[cnt];
						
						elmt->GetEdgePhysVals(boundary,bc,P,Pb);
						
						for(int j = 0; j < dim; ++j)
						{
							fgradU[j] = Array<OneD, NekDouble>(nbc);
							fgradV[j] = Array<OneD, NekDouble>(nbc);
							
						} 
						
						for(int j = 0; j < dim; ++j)
						{
							elmt->GetEdgePhysVals(boundary,bc,gradU[j],fgradU[j]);
						    elmt->GetEdgePhysVals(boundary,bc,gradV[j],fgradV[j]);
						}
						
						
						const Array<OneD, Array<OneD, NekDouble> > &normals
						= elmt->GetEdgeNormal(boundary);
						
						
						//tangential viscous tractive forces
						
						//====drag terms=========
						
						//-(du/dy+dv/dx)*ny
						Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,drag,1);
						Vmath::Neg(nbc,drag,1);
						Vmath::Vmul(nbc,drag,1,normals[1],1,drag,1);
						//-2*du/dx*nx-(du/dy+dv/dx)*ny
						Vmath::Smul(nbc,-2.0,fgradU[0],1,fgradU[0],1);
						Vmath::Vmul(nbc,fgradU[0],1,normals[0],1,temp,1);
						Vmath::Vadd(nbc,temp,1,drag,1,drag,1);
						//multiply by viscosity
						Vmath::Smul(nbc,kinvis,drag,1,drag,1);
						
						//=======lift terms=========
						
						//-(du/dy+dv/dx)*ny
						Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,lift,1);
						Vmath::Neg(nbc,lift,1);
						Vmath::Vmul(nbc,lift,1,normals[0],1,lift,1);
						//-2*dv/dy*nx-(du/dy+dv/dx)
						Vmath::Svtvp(nbc,-2.0,fgradV[1],1,lift,1,lift,1);
						Vmath::Smul(nbc,kinvis,lift,1,lift,1);
						
						
						//normal tractive forces 
						Vmath::Vvtvp(nbc,Pb,1,normals[0],1,drag,1,drag,1);
						Vmath::Vvtvp(nbc,Pb,1,normals[1],1,lift,1,lift,1);
						

						
						D=D+bc->Integral(drag);
						L=L+bc->Integral(lift);
						cout << "L " << L << endl;

					}
				}
				else 
				{
					cnt += BndExp[n]->GetExpSize();
				}
				
			}
				
		}

			cout << "L OUT" << L << endl;
			
				
			if (vComm->GetRank() == 0)
            {
				m_outputStream.width(8);
				m_outputStream << setprecision(6) << time;

				m_outputStream.width(25);
				m_outputStream << setprecision(16) << D;
                    
				m_outputStream.width(25);
				m_outputStream << setprecision(16) << L;
					
				m_outputStream << endl;
			}
                
            
        }


        /**
         *
         */
        void FilterForces::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            if (pFields[0]->GetComm()->GetRank() == 0)
            {
                m_outputStream.close();
            }
        }


        /**
         *
         */
        bool FilterForces::v_IsTimeDependent()
        {
            return true;
        }
    }
}
