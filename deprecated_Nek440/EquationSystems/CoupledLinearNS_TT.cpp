///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLinearNS.cpp
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
// Description: Coupled  Solver for the Linearised Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <boost/algorithm/string.hpp>

#include <LibUtilities/TimeIntegration/TimeIntegrationWrapper.h>
#include "CoupledLinearNS_TT.h"
#include "CoupledLinearNS_trafoP.h"
#include <LibUtilities/BasicUtils/Timer.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include "../Eigen/Dense"
#include <time.h> 

using namespace std;

namespace Nektar
{

    string CoupledLinearNS_TT::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("CoupledLinearisedNS_TT", CoupledLinearNS_TT::create);


    CoupledLinearNS_TT::CoupledLinearNS_TT(const LibUtilities::SessionReaderSharedPtr &pSession):
        UnsteadySystem(pSession),
        CoupledLinearNS(pSession),
        m_zeroMode(false)
    {



    }

    void CoupledLinearNS_TT::v_InitObject()
    {
        IncNavierStokes::v_InitObject();

        int  i;
        int  expdim = m_graph->GetMeshDimension();

        // Get Expansion list for orthogonal expansion at p-2
        const SpatialDomains::ExpansionMap &pressure_exp = GenPressureExp(m_graph->GetExpansions("u"));

        m_nConvectiveFields = m_fields.num_elements();
        if(boost::iequals(m_boundaryConditions->GetVariable(m_nConvectiveFields-1), "p"))
        {
            ASSERTL0(false,"Last field is defined as pressure but this is not suitable for this solver, please remove this field as it is implicitly defined");
        }
        // Decide how to declare explist for pressure. 
        if(expdim == 2)
        {
            int nz; 

            if(m_HomogeneousType == eHomogeneous1D)
            {
                ASSERTL0(m_fields.num_elements() > 2,"Expect to have three at least three components of velocity variables");
                LibUtilities::BasisKey Homo1DKey = m_fields[0]->GetHomogeneousBasis()->GetBasisKey();

                m_pressure = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(m_session, Homo1DKey, m_LhomZ, m_useFFT,m_homogen_dealiasing, pressure_exp);

                ASSERTL1(m_npointsZ%2==0,"Non binary number of planes have been specified");
                nz = m_npointsZ/2;                

            }
            else
            {
                //m_pressure2 = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(m_session, pressure_exp);
                m_pressure = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(m_session, pressure_exp);
                nz = 1;
            }

            Array<OneD, MultiRegions::ExpListSharedPtr> velocity(m_velocity.num_elements());
            for(i =0 ; i < m_velocity.num_elements(); ++i)
            {
                velocity[i] = m_fields[m_velocity[i]];
            }

            // Set up Array of mappings

            m_locToGloMap = Array<OneD, CoupledLocalToGlobalC0ContMapSharedPtr> (nz);
	//    AssemblyMapCGSharedPtr my_m_locToGloMap = MemoryManager<AssemblyMapCG>::AllocateSharedPtr(m_session);

            if(m_singleMode)
            {

                ASSERTL0(nz <=2 ,"For single mode calculation can only have  nz <= 2");
                if(m_session->DefinesSolverInfo("BetaZero"))
                {
                    m_zeroMode = true;
                }
                int nz_loc = 2;
                m_locToGloMap[0] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure,nz_loc,m_zeroMode);
            }
            else 
            {
                // base mode
                int nz_loc = 1;
                m_locToGloMap[0] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure,nz_loc);

                if(nz > 1)
                {
                    nz_loc = 2;
                    // Assume all higher modes have the same boundary conditions and re-use mapping
                    m_locToGloMap[1] = MemoryManager<CoupledLocalToGlobalC0ContMap>::AllocateSharedPtr(m_session,m_graph,m_boundaryConditions,velocity,m_pressure->GetPlane(2),nz_loc,false);
                    // Note high order modes cannot be singular
                    for(i = 2; i < nz; ++i)
                    {
                        m_locToGloMap[i] = m_locToGloMap[1];
                    }
                }
            }
        }
        else if (expdim == 3)
        {
            //m_pressure = MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(pressure_exp);
            ASSERTL0(false,"Setup mapping aray");
        }
        else
        {
            ASSERTL0(false,"Exp dimension not recognised");
        }

        // creation of the extrapolation object
        if(m_equationType == eUnsteadyNavierStokes)
        {
            std::string vExtrapolation = "Standard";

            if (m_session->DefinesSolverInfo("Extrapolation"))
            {
                vExtrapolation = m_session->GetSolverInfo("Extrapolation");
            }

            m_extrapolation = GetExtrapolateFactory().CreateInstance(
                vExtrapolation,
                m_session,
                m_fields,
                m_pressure,
                m_velocity,
                m_advObject);
        }


        int nel  = m_fields[0]->GetNumElmts();
//        int n_vel  = m_fields.num_elements();


/*	std::stringstream sstm_bwdtrans;
	sstm_bwdtrans << "bwdtrans.txt";
	std::string result_bwdtrans = sstm_bwdtrans.str();
//      std::ofstream myfile_bwdtrans (result_bwdtrans);
        std::ofstream myfile_bwdtrans ("bwdtrans.txt");
	if (myfile_bwdtrans.is_open())
	{
		for(int n = 0; n < nel; ++n)
		{
			int eid = n;
//	                locExp = m_fields[m_velocity[0]]->GetExp(eid);
	                int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
	                int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
        		Array<OneD, NekDouble> coeffs(ncoeffs);
		        Array<OneD, NekDouble> phys  (nphys);

			for(int counter_ncoeffs = 0; counter_ncoeffs < ncoeffs; ++counter_ncoeffs)
			{

				Vmath::Zero(ncoeffs,coeffs,1);
		                coeffs[ counter_ncoeffs ] = 1.0;
		                m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);
	
				for(int counter_phys = 0; counter_phys < nphys; ++counter_phys)         
				{
					myfile_bwdtrans << std::setprecision(17) << phys[counter_phys] << " ";
				}
				myfile_bwdtrans << endl;
	
			}
			myfile_bwdtrans << endl;


        
		}
	        myfile_bwdtrans.close();
	}
	else std::cout << "Unable to open file";
*/
	std::stringstream sstm_cartmap0;
	sstm_cartmap0 << "cartmap0.txt";
//	std::string result_cartmap0 = sstm_cartmap0.str();
//	std::basic_string<char, std::char_traits<char>, std::allocator<char> > result_cartmap0 = sstm_cartmap0.str();
//        std::ofstream myfile_cartmap0 (result_cartmap0);
  /*      std::ofstream myfile_cartmap0 ("cartmap0.txt");
	if (myfile_cartmap0.is_open())
	{
		for(int n = 0; n < nel; ++n)
		{
			int eid = n;
	                StdRegions::StdExpansionSharedPtr locExp = m_fields[m_velocity[0]]->GetExp(eid);
	                int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
	                int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
		        int pqsize  = m_pressure->GetExp(eid)->GetTotPoints();
                        Array<OneD, NekDouble> deriv  (pqsize);
        		Array<OneD, NekDouble> coeffs(ncoeffs);
		        Array<OneD, NekDouble> phys  (nphys);

			for(int counter_phys = 0; counter_phys < nphys; ++counter_phys)
			{

				Vmath::Zero(nphys,phys,1);
		                phys[ counter_phys ] = 1.0;
		                locExp->PhysDeriv(MultiRegions::DirCartesianMap[0],phys, deriv);
	
				for(int counter_deriv = 0; counter_deriv < pqsize; ++counter_deriv)         
				{
					myfile_cartmap0 << std::setprecision(17) << deriv[counter_deriv] << " ";
				}
				myfile_cartmap0 << endl;
	
			}
			myfile_cartmap0 << endl;


        
		}
	        myfile_cartmap0.close();
	}
	else std::cout << "Unable to open file";

	std::stringstream sstm_cartmap1;
	sstm_cartmap1 << "cartmap1.txt";
	std::string result_cartmap1 = sstm_cartmap1.str();
//        std::ofstream myfile_cartmap1 (result_cartmap1);
        std::ofstream myfile_cartmap1 ("cartmap1.txt");
	if (myfile_cartmap1.is_open())
	{
		for(int n = 0; n < nel; ++n)
		{
			int eid = n;
	                StdRegions::StdExpansionSharedPtr locExp = m_fields[m_velocity[0]]->GetExp(eid);
	                int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
	                int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
		        int pqsize  = m_pressure->GetExp(eid)->GetTotPoints();
                        Array<OneD, NekDouble> deriv  (pqsize);
        		Array<OneD, NekDouble> coeffs(ncoeffs);
		        Array<OneD, NekDouble> phys  (nphys);

			for(int counter_phys = 0; counter_phys < nphys; ++counter_phys)
			{

				Vmath::Zero(nphys,phys,1);
		                phys[ counter_phys ] = 1.0;
		                locExp->PhysDeriv(MultiRegions::DirCartesianMap[1],phys, deriv);
	
				for(int counter_deriv = 0; counter_deriv < pqsize; ++counter_deriv)         
				{
					myfile_cartmap1 << std::setprecision(17) << deriv[counter_deriv] << " ";
				}
				myfile_cartmap1 << endl;
	
			}
			myfile_cartmap1 << endl;


        
		}
	        myfile_cartmap1.close();
	}
	else std::cout << "Unable to open file";
*/
        // locExp->IProductWRTBase(tmpphys,coeffs);                                for all tmpphys

/*	std::stringstream sstm_IP;
	sstm_IP << "IP.txt";
	std::string result_IP = sstm_IP.str();
//        std::ofstream myfile_IP (result_IP);
        std::ofstream myfile_IP ("IP.txt");
	if (myfile_IP.is_open())
	{
		for(int n = 0; n < nel; ++n)
		{
			int eid = n;
	                StdRegions::StdExpansionSharedPtr locExp = m_fields[m_velocity[0]]->GetExp(eid);
	                int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
	                int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
		        int pqsize  = m_pressure->GetExp(eid)->GetTotPoints();
                        Array<OneD, NekDouble> deriv  (pqsize);
        		Array<OneD, NekDouble> coeffs(ncoeffs);
		        Array<OneD, NekDouble> phys  (nphys);

			for(int counter_phys = 0; counter_phys < nphys; ++counter_phys)
			{

				Vmath::Zero(nphys,phys,1);
		                phys[ counter_phys ] = 1.0;
		                locExp->IProductWRTBase(phys,coeffs);
	
				for(int counter_ncoeffs = 0; counter_ncoeffs < ncoeffs; ++counter_ncoeffs)         
				{
					myfile_IP << std::setprecision(17) << coeffs[counter_ncoeffs] << " ";
				}
				myfile_IP << endl;
	
			}
			myfile_IP << endl;


        
		}
	        myfile_IP.close();
	}
	else std::cout << "Unable to open file";
*/
/*	std::stringstream sstm_IP_d0;
	sstm_IP_d0 << "IP_d0.txt";
	std::string result_IP_d0 = sstm_IP_d0.str();
//        std::ofstream myfile_IP (result_IP);
        std::ofstream myfile_IP_d0 ("IP_d0.txt");
	if (myfile_IP_d0.is_open())
	{
		for(int n = 0; n < nel; ++n)
		{
			int eid = n;
	                StdRegions::StdExpansionSharedPtr locExp = m_fields[m_velocity[0]]->GetExp(eid);
	                int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
	                int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
		        int pqsize  = m_pressure->GetExp(eid)->GetTotPoints();
                        Array<OneD, NekDouble> deriv  (pqsize);
        		Array<OneD, NekDouble> coeffs(ncoeffs);
		        Array<OneD, NekDouble> phys  (nphys);

			for(int counter_phys = 0; counter_phys < nphys; ++counter_phys)
			{

				Vmath::Zero(nphys,phys,1);
		                phys[ counter_phys ] = 1.0;
		                locExp->IProductWRTDerivBase(0, phys,coeffs);
	
				for(int counter_ncoeffs = 0; counter_ncoeffs < ncoeffs; ++counter_ncoeffs)         
				{
					myfile_IP_d0 << std::setprecision(17) << coeffs[counter_ncoeffs] << " ";
				}
				myfile_IP_d0 << endl;
	
			}
			myfile_IP_d0 << endl;


        
		}
	        myfile_IP_d0.close();
	}
	else std::cout << "Unable to open file";


	std::stringstream sstm_IP_d1;
	sstm_IP_d1 << "IP_d1.txt";
	std::string result_IP_d1 = sstm_IP_d1.str();
//        std::ofstream myfile_IP (result_IP);
        std::ofstream myfile_IP_d1 ("IP_d1.txt");
	if (myfile_IP_d1.is_open())
	{
		for(int n = 0; n < nel; ++n)
		{
			int eid = n;
	                StdRegions::StdExpansionSharedPtr locExp = m_fields[m_velocity[0]]->GetExp(eid);
	                int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
	                int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
		        int pqsize  = m_pressure->GetExp(eid)->GetTotPoints();
                        Array<OneD, NekDouble> deriv  (pqsize);
        		Array<OneD, NekDouble> coeffs(ncoeffs);
		        Array<OneD, NekDouble> phys  (nphys);

			for(int counter_phys = 0; counter_phys < nphys; ++counter_phys)
			{

				Vmath::Zero(nphys,phys,1);
		                phys[ counter_phys ] = 1.0;
		                locExp->IProductWRTDerivBase(1, phys,coeffs);
	
				for(int counter_ncoeffs = 0; counter_ncoeffs < ncoeffs; ++counter_ncoeffs)         
				{
					myfile_IP_d1 << std::setprecision(17) << coeffs[counter_ncoeffs] << " ";
				}
				myfile_IP_d1 << endl;
	
			}
			myfile_IP_d1 << endl;


        
		}
	        myfile_IP_d1.close();
	}
	else std::cout << "Unable to open file";
*/

	//  m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);

/*	std::stringstream sstm_IPp;
	sstm_IPp << "IPp.txt";
	std::string result_IPp = sstm_IPp.str();
//        std::ofstream myfile_IPp (result_IPp);
        std::ofstream myfile_IPp ("IPp.txt");
	if (myfile_IPp.is_open())
	{
		for(int n = 0; n < nel; ++n)
		{
			int eid = n;
	                StdRegions::StdExpansionSharedPtr locExp = m_fields[m_velocity[0]]->GetExp(eid);
	                int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
	                int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
		        int pqsize  = m_pressure->GetExp(eid)->GetTotPoints();
            		int psize   = m_pressure->GetExp(eid)->GetNcoeffs();
               		Array<OneD, NekDouble> pcoeffs(psize);
                        Array<OneD, NekDouble> deriv  (pqsize);
        		Array<OneD, NekDouble> coeffs(ncoeffs);
		        Array<OneD, NekDouble> phys  (nphys);

			for(int counter_phys = 0; counter_phys < pqsize; ++counter_phys)
			{

				Vmath::Zero(pqsize,deriv,1);
		                deriv[ counter_phys ] = 1.0;
		                m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
	
				for(int counter_ncoeffs = 0; counter_ncoeffs < psize; ++counter_ncoeffs)         
				{
					myfile_IPp << std::setprecision(17) << pcoeffs[counter_ncoeffs] << " ";
				}
				myfile_IPp << endl;
	
			}
			myfile_IPp << endl;

		        Array<OneD, unsigned int> bmap,imap;

		        locExp->GetBoundaryMap(bmap);
		        locExp->GetInteriorMap(imap);

	        	int nbmap = bmap.num_elements();
		        int nimap = imap.num_elements(); 
        */


//		std::stringstream sstmb;
//		sstmb << "bmap_" << n << ".txt";
//		std::string resultb = sstmb.str();


//          std::ofstream myfile_bmap (resultb);
/*          std::ofstream myfile_bmap ("bmap.txt");
	  if (myfile_bmap.is_open())
	  {
	    for (int i = 0; i < nbmap; i++)
	    {

		    {
			myfile_bmap << std::setprecision(17) << bmap[i] << " ";
		    }

	    }
                myfile_bmap.close();
	  }
	  else std::cout << "Unable to open file";

		std::stringstream sstmi;
		sstmi << "imap_" << n << ".txt";
		std::string resulti = sstmi.str();
*/
//	    std::cout << "saving current imap " << resulti << std::endl;

//          std::ofstream myfile_imap (resulti);
/*          std::ofstream myfile_imap ("imap.txt");
	  if (myfile_imap.is_open())
	  {
	    for (int i = 0; i < nimap; i++)
	    {

		    {
			myfile_imap << std::setprecision(17) << imap[i] << " ";
		    }

	    }
                myfile_imap.close();
	  }
	  else std::cout << "Unable to open file";
*/
      //      StdRegions::ConstFactorMap factors;
      //      factors[StdRegions::eFactorLambda] = lambda/m_kinvis;
      //      LocalRegions::MatrixKey helmkey(StdRegions::eHelmholtz,
      //                                      locExp->DetShapeType(),
      //                                      *locExp,
      //                                      factors);

      //          DNekScalMat &HelmMat = *(locExp->as<LocalRegions::Expansion>()
      //                                         ->GetLocMatrix(helmkey));
/*
		}
	        myfile_IPp.close();
	}
	else std::cout << "Unable to open file";
*/
    }

    void CoupledLinearNS_TT::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation)
    {

        
        int nz;
        if(m_singleMode)
        {
            
            NekDouble lambda_imag; 
            
            // load imaginary component of any potential shift
            // Probably should be called from DriverArpack but not yet
            // clear how to do this
            m_session->LoadParameter("imagShift",lambda_imag,NekConstants::kNekUnsetDouble);
            nz = 1;
            m_mat  = Array<OneD, CoupledSolverMatrices> (nz);
            
            ASSERTL1(m_npointsZ <=2,"Only expected a maxmimum of two planes in single mode linear NS solver");
            
            if(m_zeroMode)
            {
                SetUpCoupledMatrix(lambda,Advfield,IsLinearNSEquation,0,m_mat[0],m_locToGloMap[0],lambda_imag);
            }
            else
            {
                NekDouble beta =  2*M_PI/m_LhomZ; 
                NekDouble lam = lambda + m_kinvis*beta*beta;
                
                SetUpCoupledMatrix(lam,Advfield,IsLinearNSEquation,1,m_mat[0],m_locToGloMap[0],lambda_imag);
            }
        }
        else 
        {
            int n;
            if(m_npointsZ > 1)
            { 
                nz = m_npointsZ/2;
            }
            else
            {
                nz =  1;
            }
            
            m_mat  = Array<OneD, CoupledSolverMatrices> (nz);
            
            // mean mode or 2D mode.
            SetUpCoupledMatrix(lambda,Advfield,IsLinearNSEquation,0,m_mat[0],m_locToGloMap[0]);
            
            for(n = 1; n < nz; ++n)
            {
                NekDouble beta = 2*M_PI*n/m_LhomZ;
                
                NekDouble lam = lambda + m_kinvis*beta*beta;
                
                SetUpCoupledMatrix(lam,Advfield,IsLinearNSEquation,n,m_mat[n],m_locToGloMap[n]);
            }
        }
        


    }

    void CoupledLinearNS_TT::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation,const int HomogeneousMode, CoupledSolverMatrices &mat, CoupledLocalToGlobalC0ContMapSharedPtr &locToGloMap, const NekDouble lambda_imag)
    {

	if(use_Newton)        
	{
		IsLinearNSEquation = true;
	}

        int  n,i,j,k,eid;
        int  nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int  nvel   = m_velocity.num_elements();
        
        // if Advfield is defined can assume it is an Oseen or LinearNS equation
        bool AddAdvectionTerms = (Advfield ==  NullNekDoubleArrayofArray)? false: true; 
        //bool AddAdvectionTerms = true; // Temporary debugging trip
        
        // call velocity space Static condensation and fill block
        // matrices.  Need to set this up differently for Oseen and
        // Lin NS.  Ideally should make block diagonal for Stokes and
        // Oseen problems.
        DNekScalMatSharedPtr loc_mat;
        StdRegions::StdExpansionSharedPtr locExp;
        NekDouble one = 1.0;
        int nint,nbndry;
        int rows, cols;
        NekDouble zero = 0.0;
        Array<OneD, unsigned int> bmap,imap; 
        
        Array<OneD,unsigned int> nsize_bndry   (nel);
        Array<OneD,unsigned int> nsize_bndry_p1(nel);
        Array<OneD,unsigned int> nsize_int     (nel);
        Array<OneD,unsigned int> nsize_p       (nel);
        Array<OneD,unsigned int> nsize_p_m1    (nel);
        
	


        int nz_loc;
        
        if(HomogeneousMode) // Homogeneous mode flag
        {
            nz_loc = 2;
        }
        else
        {
            if(m_singleMode)
            {
                nz_loc = 2;
            }
            else
            {
                nz_loc = 1;
            }
        }
        
        // Set up block matrix sizes - 
        for(n = 0; n < nel; ++n)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(n);
            nsize_bndry[n] = nvel*m_fields[m_velocity[0]]->GetExp(eid)->NumBndryCoeffs()*nz_loc;
            nsize_bndry_p1[n] = nsize_bndry[n]+nz_loc;
            nsize_int[n] = (nvel*m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs()*nz_loc - nsize_bndry[n]);
            nsize_p[n] = m_pressure->GetExp(eid)->GetNcoeffs()*nz_loc;
            nsize_p_m1[n] = nsize_p[n]-nz_loc;
     /*       std::cout << "nsize_bndry[n] " << nsize_bndry[n] << std::endl;
            std::cout << "nsize_bndry_p1[n] " << nsize_bndry_p1[n] << std::endl;
            std::cout << "nsize_int[n] " << nsize_int[n] << std::endl;
            std::cout << "nsize_p[n] " << nsize_p[n] << std::endl;
            std::cout << "nsize_p_m1[n] " << nsize_p_m1[n] << std::endl; */
        }

	// need my mats for projection purposes: think about: have no splitting bnd/p/int could be beneficial for number of affine terms
/*	Eigen::MatrixXd D_all = Eigen::MatrixXd::Zero( nsize_int[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd Dbnd_all = Eigen::MatrixXd::Zero( nsize_p[0]*nel , nsize_bndry[0]*nel );
	Eigen::MatrixXd Dint_all = Eigen::MatrixXd::Zero( nsize_p[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd B_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd A_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_bndry[0]*nel );

	Eigen::MatrixXd D_adv_all = Eigen::MatrixXd::Zero( nsize_int[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd B_adv_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd C_adv_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd A_adv_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_bndry[0]*nel );

	Eigen::MatrixXd D_no_adv_all = Eigen::MatrixXd::Zero( nsize_int[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd B_no_adv_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd C_no_adv_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_int[0]*nel );
	Eigen::MatrixXd A_no_adv_all = Eigen::MatrixXd::Zero( nsize_bndry[0]*nel , nsize_bndry[0]*nel );
*/
//	Eigen::MatrixXd Dd_all = Eigen::MatrixXd::Zero(  nsize_int[0]*nel , nsize_int[0]*nel );
        
//	cout << typeid(nsize_int[0]).name() << endl;
//	cout << typeid(nsize_bndry[0]).name() << endl;
//	cout <<  Dd_all << endl;

        MatrixStorage blkmatStorage = eDIAGONAL;
        DNekScalBlkMatSharedPtr pAh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry_p1,nsize_bndry_p1,blkmatStorage);
        mat.m_BCinv = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
        mat.m_Btilde = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry,nsize_int,blkmatStorage);
        mat.m_Cinv = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_int,nsize_int,blkmatStorage);
        
        mat.m_D_bnd = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p,nsize_bndry,blkmatStorage);
        
        mat.m_D_int = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p,nsize_int,blkmatStorage);
        
        // Final level static condensation matrices. 
        DNekScalBlkMatSharedPtr pBh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_bndry_p1,nsize_p_m1,blkmatStorage);
        DNekScalBlkMatSharedPtr pCh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p_m1,nsize_bndry_p1,blkmatStorage);
        DNekScalBlkMatSharedPtr pDh = MemoryManager<DNekScalBlkMat>
        ::AllocateSharedPtr(nsize_p_m1,nsize_p_m1,blkmatStorage);
        
        
        Timer timer;
        timer.Start();
        for(n = 0; n < nel; ++n)
        {
            eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(n);
            nbndry = nsize_bndry[n];
            nint   = nsize_int[n];
              
	 //       cout << "nint: " << nint << endl;
	 //       cout << "nbndry: " << nbndry << endl;

            k = nsize_bndry_p1[n];
            DNekMatSharedPtr Ah = MemoryManager<DNekMat>::AllocateSharedPtr(k,k,zero);
            DNekMatSharedPtr A_adv = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry+1,nbndry+1,zero);  // the +1 should not be necessary but would have to change indexing for that
            DNekMatSharedPtr A_no_adv = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry+1,nbndry+1,zero);
            Array<OneD, NekDouble> Ah_data = Ah->GetPtr();
	    Array<OneD, NekDouble> A_adv_data = A_adv->GetPtr();
	    Array<OneD, NekDouble> A_no_adv_data = A_no_adv->GetPtr();
            int AhRows = k;
            DNekMatSharedPtr B  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            DNekMatSharedPtr B_adv  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            DNekMatSharedPtr B_no_adv  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            Array<OneD, NekDouble> B_data = B->GetPtr();
	    Array<OneD, NekDouble> B_adv_data = B_adv->GetPtr();
	    Array<OneD, NekDouble> B_no_adv_data = B_no_adv->GetPtr();
            DNekMatSharedPtr C  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            DNekMatSharedPtr C_adv  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            DNekMatSharedPtr C_no_adv  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            Array<OneD, NekDouble> C_data = C->GetPtr();
	    Array<OneD, NekDouble> C_adv_data = C_adv->GetPtr();
	    Array<OneD, NekDouble> C_no_adv_data = C_no_adv->GetPtr();
            DNekMatSharedPtr D  = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint,zero);
            DNekMatSharedPtr D_adv  = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint,zero);
            DNekMatSharedPtr D_no_adv  = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint,zero);
            Array<OneD, NekDouble> D_data = D->GetPtr();
            Array<OneD, NekDouble> D_adv_data = D_adv->GetPtr();
            Array<OneD, NekDouble> D_no_adv_data = D_no_adv->GetPtr();
            
            DNekMatSharedPtr Dbnd = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p[n],nsize_bndry[n],zero);
            DNekMatSharedPtr Dint = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p[n],nsize_int[n],zero);
            
            locExp = m_fields[m_velocity[0]]->GetExp(eid);
            locExp->GetBoundaryMap(bmap);
            locExp->GetInteriorMap(imap);
            StdRegions::ConstFactorMap factors;
            factors[StdRegions::eFactorLambda] = lambda/m_kinvis;
/*            LocalRegions::MatrixKey helmkey(StdRegions::eHelmholtz,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);
*/

            LocalRegions::MatrixKey helmkey(StdRegions::eHelmholtz,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);


/*            LocalRegions::MatrixKey helmkey(StdRegions::eLaplacian,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);

  */          
            
            int ncoeffs = m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs();
            int nphys   = m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints();
            int nbmap = bmap.num_elements();
            int nimap = imap.num_elements(); 

	 //   std::cout << "nimap " << nimap  << std::endl;
            
            Array<OneD, NekDouble> coeffs(ncoeffs);
            Array<OneD, NekDouble> phys  (nphys);
            int psize   = m_pressure->GetExp(eid)->GetNcoeffs();
            int pqsize  = m_pressure->GetExp(eid)->GetTotPoints();
            
            Array<OneD, NekDouble> deriv  (pqsize);
            Array<OneD, NekDouble> pcoeffs(psize);
            
            if(AddAdvectionTerms == false) // use static condensed managed matrices
            {
                // construct velocity matrices using statically
                // condensed elemental matrices and then construct
                // pressure matrix systems
                DNekScalBlkMatSharedPtr CondMat; 
                CondMat = locExp->GetLocStaticCondMatrix(helmkey);
                
                for(k = 0; k < nvel*nz_loc; ++k)
                {
                    DNekScalMat &SubBlock = *CondMat->GetBlock(0,0);
                    rows = SubBlock.GetRows();
                    cols = SubBlock.GetColumns(); 
                    for(i = 0; i < rows; ++i)
                    {
                        for(j = 0; j < cols; ++j)
                        {
                            (*Ah)(i+k*rows,j+k*cols) = m_kinvis*SubBlock(i,j);
                        }                        
                    }
                }
                
                for(k = 0; k < nvel*nz_loc; ++k)
                {
                    DNekScalMat &SubBlock  = *CondMat->GetBlock(0,1);
                    DNekScalMat &SubBlock1 = *CondMat->GetBlock(1,0);
                    rows = SubBlock.GetRows();
                    cols = SubBlock.GetColumns(); 
                    for(i = 0; i < rows; ++i)
                    {
                        for(j = 0; j < cols; ++j)
                        {
                            (*B)(i+k*rows,j+k*cols) = SubBlock(i,j);
                            (*C)(i+k*rows,j+k*cols) = m_kinvis*SubBlock1(j,i);
                        }
                    }
                }
                
                for(k = 0; k < nvel*nz_loc; ++k)
                {
                    DNekScalMat &SubBlock = *CondMat->GetBlock(1,1);
                    NekDouble inv_kinvis = 1.0/m_kinvis;
                    rows = SubBlock.GetRows();
                    cols = SubBlock.GetColumns(); 
                    for(i = 0; i < rows; ++i)
                    {
                        for(j = 0; j < cols; ++j)
                        {
                            (*D)(i+k*rows,j+k*cols) = inv_kinvis*SubBlock(i,j);
                        }
                    }
                }
                
                
                // Loop over pressure space and construct boundary block matrices.         
                for(i = 0; i < bmap.num_elements(); ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[bmap[i]] = 1.0;
                    m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);
                    
                    // Differentiation & Inner product wrt base. 
                    for(j = 0; j < nvel; ++j)
                    {
                        if( (nz_loc == 2)&&(j == 2)) // handle d/dz derivative
                        {
                            NekDouble beta =  -2*M_PI*HomogeneousMode/m_LhomZ;
                            
                            Vmath::Smul(m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints(), beta, phys,1,deriv,1);
                            
                            m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                            
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() + 
                                        ((nz_loc*j+1)*bmap.num_elements()+i)*nsize_p[n],1);
                            
                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() + 
                                        ((nz_loc*j)*bmap.num_elements()+i)*nsize_p[n]+psize,1);
                            
                        }
                        else
                        {
                            if(j < 2) // required for mean mode of homogeneous expansion 
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[j],phys,deriv);		
                                m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                                // copy into column major storage. 
                                for(k = 0; k < nz_loc; ++k)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dbnd->GetRawPtr() + 
                                                ((nz_loc*j+k)*bmap.num_elements()+i)*nsize_p[n]+ k*psize,1);
                                }
                            }
                        }
                    }
                }
                
                for(i = 0; i < imap.num_elements(); ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[imap[i]] = 1.0;
                    m_fields[m_velocity[0]]->GetExp(eid)->BwdTrans(coeffs,phys);
                    
                    // Differentiation & Inner product wrt base. 
                    for(j = 0; j < nvel; ++j)
                    {
                        if( (nz_loc == 2)&&(j == 2)) // handle d/dz derivative
                        {
                            NekDouble beta = -2*M_PI*HomogeneousMode/m_LhomZ;
                            
                            Vmath::Smul(m_fields[m_velocity[0]]->GetExp(eid)->GetTotPoints(), beta, phys,1,deriv,1); 
                            
                            m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                            
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() + 
                                        ((nz_loc*j+1)*imap.num_elements()+i)*nsize_p[n],1);
                            
                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() + 
                                        ((nz_loc*j)*imap.num_elements()+i)*nsize_p[n]+psize,1);
                            
                        }
                        else
                        {
                            if(j < 2) // required for mean mode of homogeneous expansion 
                            {
                                //m_fields[m_velocity[0]]->GetExp(eid)->PhysDeriv(j,phys, deriv);
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[j],phys,deriv);
                                
                                m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                                
                                // copy into column major storage. 
                                for(k = 0; k < nz_loc; ++k)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dint->GetRawPtr() +
                                                ((nz_loc*j+k)*imap.num_elements()+i)*nsize_p[n]+ k*psize,1);
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                // construct velocity matrices and pressure systems at
                // the same time resusing differential of velocity
                // space
                
	        LocalRegions::MatrixKey helmkey_l00(StdRegions::eLaplacian00,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);


		DNekScalMat &HelmMat = *(locExp->as<LocalRegions::Expansion>()
                                               ->GetLocMatrix(helmkey_l00));

                Array<OneD, const NekDouble> HelmMat_data = HelmMat.GetOwnedMatrix()->GetPtr();
              

/*		std::stringstream sstm_l00;
		sstm_l00 << "Lapl00_" << n << ".txt";
		std::string result_l00 = sstm_l00.str();
		const char* rr_l00 = result_l00.c_str();

        	  std::ofstream myfileHM_l00 (rr_l00);
		  if (myfileHM_l00.is_open())
		  {
		    for (int i = 0; i < HelmMat_data.num_elements(); i++)
		    {

			    {
				myfileHM_l00 << std::setprecision(17) << HelmMat_data[i] << " ";
			    }

		    }
        	        myfileHM_l00.close();
		  }
		  else std::cout << "Unable to open file";

*/



////////////////////////////////////////////////////////////////////////////////////////////////

/*
	        LocalRegions::MatrixKey helmkey_l11(StdRegions::eLaplacian11,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);


		HelmMat = *(locExp->as<LocalRegions::Expansion>()
                                               ->GetLocMatrix(helmkey_l11));

              
                HelmMat_data = HelmMat.GetOwnedMatrix()->GetPtr();
                

		std::stringstream sstm_l11;
		sstm_l11 << "Lapl11_" << n << ".txt";
		std::string result_l11 = sstm_l11.str();
		const char* rr_l11 = result_l11.c_str();

        	  std::ofstream myfileHM_l11 (rr_l11);
		  if (myfileHM_l11.is_open())
		  {
		    for (int i = 0; i < HelmMat_data.num_elements(); i++)
		    {

			    {
				myfileHM_l11 << std::setprecision(17) << HelmMat_data[i] << " ";
			    }

		    }
        	        myfileHM_l11.close();
		  }
		  else std::cout << "Unable to open file";
*/
////////////////////////////////////////////////////////////////////////////////////////////////

/*
	        LocalRegions::MatrixKey helmkey_l01(StdRegions::eLaplacian01,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);


		HelmMat = *(locExp->as<LocalRegions::Expansion>()
                                               ->GetLocMatrix(helmkey_l01));

              
                HelmMat_data = HelmMat.GetOwnedMatrix()->GetPtr();
                

		std::stringstream sstm_l01;
		sstm_l01 << "Lapl01_" << n << ".txt";
		std::string result_l01 = sstm_l01.str();
		const char* rr_l01 = result_l01.c_str();

        	  std::ofstream myfileHM_l01 (rr_l01);
		  if (myfileHM_l01.is_open())
		  {
		    for (int i = 0; i < HelmMat_data.num_elements(); i++)
		    {

			    {
				myfileHM_l01 << std::setprecision(17) << HelmMat_data[i] << " ";
			    }

		    }
        	        myfileHM_l01.close();
		  }
		  else std::cout << "Unable to open file";
*/
////////////////////////////////////////////////////////////////////////////////////////////////


/*	        LocalRegions::MatrixKey helmkey_w0(StdRegions::eWeakDeriv0,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);


		HelmMat = *(locExp->as<LocalRegions::Expansion>()
                                               ->GetLocMatrix(helmkey_w0));

              
                HelmMat_data = HelmMat.GetOwnedMatrix()->GetPtr();
                

		std::stringstream sstm_w0;
		sstm_w0 << "w0_" << n << ".txt";
		std::string result_w0 = sstm_w0.str();
		const char* rr_w0 = result_w0.c_str();

        	  std::ofstream myfileHM_w0 (rr_w0);
		  if (myfileHM_w0.is_open())
		  {
		    for (int i = 0; i < HelmMat_data.num_elements(); i++)
		    {

			    {
				myfileHM_w0 << std::setprecision(17) << HelmMat_data[i] << " ";
			    }

		    }
        	        myfileHM_w0.close();
		  }
		  else std::cout << "Unable to open file";
*/


////////////////////////////////////////////////////////////////////////////////////////////////


/*	        LocalRegions::MatrixKey helmkey_w1(StdRegions::eWeakDeriv1,
                                            locExp->DetShapeType(),
                                            *locExp,
                                            factors);


		HelmMat = *(locExp->as<LocalRegions::Expansion>()
                                               ->GetLocMatrix(helmkey_w1));

              
                HelmMat_data = HelmMat.GetOwnedMatrix()->GetPtr();
                

		std::stringstream sstm_w1;
		sstm_w1 << "w1_" << n << ".txt";
		std::string result_w1 = sstm_w1.str();
		const char* rr_w1 = result_w1.c_str();

        	  std::ofstream myfileHM_w1 (rr_w1);
		  if (myfileHM_w1.is_open())
		  {
		    for (int i = 0; i < HelmMat_data.num_elements(); i++)
		    {

			    {
				myfileHM_w1 << std::setprecision(17) << HelmMat_data[i] << " ";
			    }

		    }
        	        myfileHM_w1.close();
		  }
		  else std::cout << "Unable to open file";
*/


///////////////////////////////////////////////////////////////////////////////////////////////



                HelmMat = *(locExp->as<LocalRegions::Expansion>()
                                               ->GetLocMatrix(helmkey));

                DNekScalMatSharedPtr MassMat;
                
                HelmMat_data = HelmMat.GetOwnedMatrix()->GetPtr();
                NekDouble HelmMatScale = HelmMat.Scale();
                int HelmMatRows = HelmMat.GetRows();
                
                if((lambda_imag != NekConstants::kNekUnsetDouble)&&(nz_loc == 2))
                {
                    LocalRegions::MatrixKey masskey(StdRegions::eMass,
                                                    locExp->DetShapeType(),
                                                    *locExp);
                    MassMat = locExp->as<LocalRegions::Expansion>()
                                    ->GetLocMatrix(masskey);
                }
                
                Array<OneD, NekDouble> Advtmp;
                Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
                // Use ExpList phys array for temporaary storage
                Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
                int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(eid);
                int nv;
                int npoints = locExp->GetTotPoints();
                

	/*	std::stringstream sstm;
		sstm << "HelmMat_" << n << ".txt";
		std::string result = sstm.str();
		const char* rr = result.c_str();

        	  std::ofstream myfileHM (rr);
		  if (myfileHM.is_open())
		  {
		    for (int i = 0; i < HelmMat_data.num_elements(); i++)
		    {

			    {
				myfileHM << std::setprecision(17) << HelmMat_data[i] << " ";
			    }

		    }
        	        myfileHM.close();
		  }
		  else std::cout << "Unable to open file";
*/

                // Calculate derivative of base flow 
                if(IsLinearNSEquation)
                {
                    int nv1;
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(nv = 0; nv < nvel; ++nv)
                    {
                        for(nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }
                            
                            if((nv1 == 2)&&(m_HomogeneousType == eHomogeneous1D))
                            {
                                Vmath::Zero(npoints,AdvDeriv[nv*nvel+nv1],1); // dU/dz = 0
                            }
                            else
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],Advfield[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                            }
                        }
                    }
                }
                
                
//	        cout << "nbmap: " << nbmap << endl;
//	        cout << "nbndry: " << nbndry << endl;

                for(i = 0; i < nbmap; ++i)
                {
                    
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[bmap[i]] = 1.0;
                    locExp->BwdTrans(coeffs,phys);
                    
                    for(k = 0; k < nvel*nz_loc; ++k)
                    {
                        for(j = 0; j < nbmap; ++j)
                        {
                            //                            Ah_data[i+k*nbmap + (j+k*nbmap)*AhRows] += m_kinvis*HelmMat(bmap[i],bmap[j]);
                            Ah_data[i+k*nbmap + (j+k*nbmap)*AhRows] += m_kinvis*HelmMatScale*HelmMat_data[bmap[i] + HelmMatRows*bmap[j]];
//                            A_no_adv_data[i+k*nbmap + (j+k*nbmap)*AhRows] += m_kinvis*HelmMatScale*HelmMat_data[bmap[i] + HelmMatRows*bmap[j]];
                        }
                        
                        for(j = 0; j < nimap; ++j)
                        {
                            B_data[i+k*nbmap + (j+k*nimap)*nbndry] += m_kinvis*HelmMatScale*HelmMat_data[bmap[i]+HelmMatRows*imap[j]];
//                            B_no_adv_data[i+k*nbmap + (j+k*nimap)*nbndry] += m_kinvis*HelmMatScale*HelmMat_data[bmap[i]+HelmMatRows*imap[j]];
                        }
                    }
                    
                    if((lambda_imag != NekConstants::kNekUnsetDouble)&&(nz_loc == 2))
                    {
                        for(k = 0; k < nvel; ++k)
                        {
                            for(j = 0; j < nbmap; ++j)
                            {
                                Ah_data[i+2*k*nbmap + (j+(2*k+1)*nbmap)*AhRows] -= lambda_imag*(*MassMat)(bmap[i],bmap[j]);
                            }
                            
                            for(j = 0; j < nbmap; ++j)
                            {
                                Ah_data[i+(2*k+1)*nbmap + (j+2*k*nbmap)*AhRows] += lambda_imag*(*MassMat)(bmap[i],bmap[j]);
                            }
                            
                            for(j = 0; j < nimap; ++j)
                            {
                                B_data[i+2*k*nbmap + (j+(2*k+1)*nimap)*nbndry] -= lambda_imag*(*MassMat)(bmap[i],imap[j]);
                            }
                            
                            for(j = 0; j < nimap; ++j)
                            {
                                B_data[i+(2*k+1)*nbmap + (j+2*k*nimap)*nbndry] += lambda_imag*(*MassMat)(bmap[i],imap[j]);
                            }
                            
                        }
                    }
                    
                    
                    
                    for(k = 0; k < nvel; ++k)
                    {
                        if((nz_loc == 2)&&(k == 2)) // handle d/dz derivative
                        { 
                            NekDouble beta = -2*M_PI*HomogeneousMode/m_LhomZ;
                            
                            // Real Component
                            Vmath::Smul(npoints,beta,phys,1,deriv,1);
                            
                            m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() + 
                                        ((nz_loc*k+1)*bmap.num_elements()+i)*nsize_p[n],1);
                            
                            // Imaginary Component
                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dbnd->GetRawPtr() + 
                                        ((nz_loc*k)*bmap.num_elements()+i)*nsize_p[n]+psize,1);
                            
                            // now do advection terms
                            Vmath::Vmul(npoints, 
                                        Advtmp = Advfield[k] + phys_offset,
                                        1,deriv,1,tmpphys,1);
                            
                            locExp->IProductWRTBase(tmpphys,coeffs);
                            
                            
                            // real contribution
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    Ah_data[j+2*nv*nbmap + (i+(2*nv+1)*nbmap)*AhRows] +=
                                    coeffs[bmap[j]];
                                }
                                
                                for(j = 0; j < nimap; ++j)
                                {
                                    C_data[i+(2*nv+1)*nbmap + (j+2*nv*nimap)*nbndry] += 
                                    coeffs[imap[j]];
                                }
                            }
                            
                            Vmath::Neg(ncoeffs,coeffs,1);
                            // imaginary contribution
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    Ah_data[j+(2*nv+1)*nbmap + (i+2*nv*nbmap)*AhRows] +=
                                    coeffs[bmap[j]];
                                }
                                
                                for(j = 0; j < nimap; ++j)
                                {
                                    C_data[i+2*nv*nbmap + (j+(2*nv+1)*nimap)*nbndry] += 
                                    coeffs[imap[j]];
                                }
                            }
                        }
                        else
                        {
                            if(k < 2)
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[k],phys, deriv);
                                Vmath::Vmul(npoints, 
                                            Advtmp = Advfield[k] + phys_offset,
                                            1,deriv,1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                
                                for(nv = 0; nv < nvel*nz_loc; ++nv)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        Ah_data[j+nv*nbmap + (i+nv*nbmap)*AhRows] += coeffs[bmap[j]];
            //                            A_adv_data[j+nv*nbmap + (i+nv*nbmap)*AhRows] += coeffs[bmap[j]];
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        C_data[i+nv*nbmap + (j+nv*nimap)*nbndry] += coeffs[imap[j]];
          //                              C_adv_data[i+nv*nbmap + (j+nv*nimap)*nbndry] += coeffs[imap[j]];
                                    }
                                }
                                
                                // copy into column major storage. 
                                m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                                for(j = 0; j < nz_loc; ++j)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1, Dbnd->GetRawPtr() + ((nz_loc*k+j)*bmap.num_elements() + i)*nsize_p[n]+ j*psize,1);
                                }
                            }
                        }
                        
                        if(IsLinearNSEquation)
                        {
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints,phys,1,
                                            AdvDeriv[k*nvel+nv],
                                            1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        Ah_data[j+(k*nz_loc+n1)*nbmap + 
                                        (i+(nv*nz_loc+n1)*nbmap)*AhRows] +=
                                        coeffs[bmap[j]];
      //                                  A_adv_data[j+(k*nz_loc+n1)*nbmap + 
        //                                (i+(nv*nz_loc+n1)*nbmap)*AhRows] +=
          //                              coeffs[bmap[j]];
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        C_data[i+(nv*nz_loc+n1)*nbmap + 
                                        (j+(k*nz_loc+n1)*nimap)*nbndry] += 
                                        coeffs[imap[j]];
//                                        C_adv_data[i+(nv*nz_loc+n1)*nbmap + 
  //                                      (j+(k*nz_loc+n1)*nimap)*nbndry] += 
    //                                    coeffs[imap[j]];
                                    }
                                }
                            }                            
                        }
                    }
                }
                
                
                for(i = 0; i < nimap; ++i)
                {
                    // Fill element with mode
                    Vmath::Zero(ncoeffs,coeffs,1);
                    coeffs[imap[i]] = 1.0;
                    locExp->BwdTrans(coeffs,phys);
                    
                    for(k = 0; k < nvel*nz_loc; ++k)
                    {
                        for(j = 0; j < nbmap; ++j) // C set up as transpose
                        {
                            C_data[j+k*nbmap + (i+k*nimap)*nbndry] += m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*bmap[j]];
            //                C_no_adv_data[j+k*nbmap + (i+k*nimap)*nbndry] += m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*bmap[j]];
                        }
                        
                        for(j = 0; j < nimap; ++j)
                        {
                            D_data[i+k*nimap + (j+k*nimap)*nint] += m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*imap[j]];
          //                  D_no_adv_data[i+k*nimap + (j+k*nimap)*nint] += m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*imap[j]];
		//	    cout << "writing D_data at location " << i+k*nimap + (j+k*nimap)*nint << " the value " << m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*imap[j]] << endl;
                        }
                    }
                    
                    if((lambda_imag != NekConstants::kNekUnsetDouble)&&(nz_loc == 2))
                    {
                        for(k = 0; k < nvel; ++k)
                        {
                            for(j = 0; j < nbmap; ++j) // C set up as transpose
                            {
                                C_data[j+2*k*nbmap + (i+(2*k+1)*nimap)*nbndry] += lambda_imag*(*MassMat)(bmap[j],imap[i]);
                            }
                            
                            for(j = 0; j < nbmap; ++j) // C set up as transpose
                            {
                                C_data[j+(2*k+1)*nbmap + (i+2*k*nimap)*nbndry] -= lambda_imag*(*MassMat)(bmap[j],imap[i]);
                            }
                            
                            for(j = 0; j < nimap; ++j)
                            {
                                D_data[i+2*k*nimap + (j+(2*k+1)*nimap)*nint] -= lambda_imag*(*MassMat)(imap[i],imap[j]);
                            }
                            
                            for(j = 0; j < nimap; ++j)
                            {
                                D_data[i+(2*k+1)*nimap + (j+2*k*nimap)*nint] += lambda_imag*(*MassMat)(imap[i],imap[j]);
                            }
                        }
                    }
                    
                    
                    for(k = 0; k < nvel; ++k)
                    {
                        if((nz_loc == 2)&&(k == 2)) // handle d/dz derivative
                        { 
                            NekDouble beta = -2*M_PI*HomogeneousMode/m_LhomZ;
                            
                            // Real Component
                            Vmath::Smul(npoints,beta,phys,1,deriv,1);
                            m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() + 
                                        ((nz_loc*k+1)*imap.num_elements()+i)*nsize_p[n],1);
                            // Imaginary Component
                            Vmath::Neg(psize,pcoeffs,1);
                            Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                        Dint->GetRawPtr() + 
                                        ((nz_loc*k)*imap.num_elements()+i)*nsize_p[n]+psize,1);
                            
                            // Advfield[k] *d/dx_k to all velocity
                            // components on diagonal
                            Vmath::Vmul(npoints, Advtmp = Advfield[k] + phys_offset,1,deriv,1,tmpphys,1);
                            locExp->IProductWRTBase(tmpphys,coeffs);
                            
                            // Real Components
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    B_data[j+2*nv*nbmap + (i+(2*nv+1)*nimap)*nbndry] += 
                                    coeffs[bmap[j]];
                                }
                                
                                for(j = 0; j < nimap; ++j)
                                {
                                    D_data[j+2*nv*nimap + (i+(2*nv+1)*nimap)*nint] += 
                                    coeffs[imap[j]];
                                }
                            }
                            Vmath::Neg(ncoeffs,coeffs,1);
                            // Imaginary 
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                for(j = 0; j < nbmap; ++j)
                                {
                                    B_data[j+(2*nv+1)*nbmap + (i+2*nv*nimap)*nbndry] += 
                                    coeffs[bmap[j]];
                                }
                                
                                for(j = 0; j < nimap; ++j)
                                {
                                    D_data[j+(2*nv+1)*nimap + (i+2*nv*nimap)*nint] += 
                                    coeffs[imap[j]];
                                }
                            }
                            
                        }
                        else
                        {
                            if(k < 2)
                            {
                                // Differentiation & Inner product wrt base. 
                                //locExp->PhysDeriv(k,phys, deriv);
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[k],phys,deriv);
                                Vmath::Vmul(npoints, 
                                            Advtmp = Advfield[k] + phys_offset,
                                            1,deriv,1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                
                                for(nv = 0; nv < nvel*nz_loc; ++nv)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        B_data[j+nv*nbmap + (i+nv*nimap)*nbndry] += coeffs[bmap[j]];
	//				B_adv_data[j+nv*nbmap + (i+nv*nimap)*nbndry] += coeffs[bmap[j]];
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        D_data[j+nv*nimap + (i+nv*nimap)*nint] += coeffs[imap[j]];
      //                                  D_adv_data[j+nv*nimap + (i+nv*nimap)*nint] += coeffs[imap[j]];
	    //                            cout << "writing D_data advec at location " << j+nv*nimap + (i+nv*nimap)*nint << " the value " << coeffs[imap[j]] << endl;
                                    }
                                }
                                // copy into column major storage. 
                                m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                                for(j = 0; j < nz_loc; ++j)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dint->GetRawPtr() + 
                                                ((nz_loc*k+j)*imap.num_elements() + i)*nsize_p[n]+j*psize,1);
                                }
                            }
                        }
                        
                        if(IsLinearNSEquation)
                        {
                            int n1;
                            for(nv = 0; nv < nvel; ++nv)
                            {
                                // u'.Grad U terms 
                                Vmath::Vmul(npoints,phys,1, 
                                            AdvDeriv[k*nvel+nv],
                                            1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(j = 0; j < nbmap; ++j)
                                    {
                                        B_data[j+(k*nz_loc+n1)*nbmap + 
                                        (i+(nv*nz_loc+n1)*nimap)*nbndry] += 
                                        coeffs[bmap[j]];
//                                        B_adv_data[j+(k*nz_loc+n1)*nbmap + 
  //                                      (i+(nv*nz_loc+n1)*nimap)*nbndry] += 
    //                                    coeffs[bmap[j]];
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        D_data[j+(k*nz_loc+n1)*nimap +  
                                        (i+(nv*nz_loc+n1)*nimap)*nint] += 
                                        coeffs[imap[j]];
//                                        D_adv_data[j+(k*nz_loc+n1)*nimap +  
//                                        (i+(nv*nz_loc+n1)*nimap)*nint] += 
//                                        coeffs[imap[j]];
                                    }
                                }
                            }
                        }
                    }
                }
                

		// from here starts the multilevel static condensation
		// so collect data for ROM
		

/*                cout << "D->GetRows() " << D->GetRows() << std::endl ;
                cout << "D->GetColumns() " << D->GetColumns() << std::endl ;
                cout << "Dbnd->GetRows() " << Dbnd->GetRows() << std::endl ;
                cout << "Dbnd->GetColumns() " << Dbnd->GetColumns() << std::endl ;
                cout << "Dint->GetRows() " << Dint->GetRows() << std::endl ;
                cout << "Dint->GetColumns() " << Dint->GetColumns() << std::endl ;
                cout << "B->GetRows() " << B->GetRows() << std::endl ;
                cout << "B->GetColumns() " << B->GetColumns() << std::endl ;
                cout << "C->GetRows() " << C->GetRows() << std::endl ;
                cout << "C->GetColumns() " << C->GetColumns() << std::endl ;
                cout << "Ah->GetRows() " << Ah->GetRows() << std::endl ;
                cout << "Ah->GetColumns() " << Ah->GetColumns() << std::endl ; */

//		Array<OneD, Array<OneD, NekDouble> > A_save = Ah;           
//		DNekMatSharedPtr D_save = D;
//     		DNekMat D_save1 = *D_save;
//		Array<OneD, Array<OneD, NekDouble> > D_save2 = D_save1;

/*		Eigen::MatrixXd D_Matrix(D->GetRows(), D->GetColumns());
		Eigen::MatrixXd D_adv_Matrix(D_adv->GetRows(), D_adv->GetColumns());
		Eigen::MatrixXd D_no_adv_Matrix(D_no_adv->GetRows(), D_no_adv->GetColumns());
		Eigen::MatrixXd Dbnd_Matrix(Dbnd->GetRows(), Dbnd->GetColumns());
		Eigen::MatrixXd Dint_Matrix(Dint->GetRows(), Dint->GetColumns());
		Eigen::MatrixXd B_Matrix(B->GetRows(), B->GetColumns());
		Eigen::MatrixXd B_adv_Matrix(B_adv->GetRows(), B_adv->GetColumns());
		Eigen::MatrixXd B_no_adv_Matrix(B_no_adv->GetRows(), B_no_adv->GetColumns());
		Eigen::MatrixXd C_Matrix(C->GetRows(), C->GetColumns());
		Eigen::MatrixXd C_adv_Matrix(C_adv->GetRows(), C_adv->GetColumns());
		Eigen::MatrixXd C_no_adv_Matrix(C_no_adv->GetRows(), C_no_adv->GetColumns());
		Eigen::MatrixXd A_Matrix(Ah->GetRows() - 1, Ah->GetColumns() - 1);
		Eigen::MatrixXd A_adv_Matrix(A_adv->GetRows() - 1, A_adv->GetColumns() - 1);
		Eigen::MatrixXd A_no_adv_Matrix(A_no_adv->GetRows() - 1, A_no_adv->GetColumns() - 1);
*/
//		Eigen::VectorXd v(m_equ[0]->GetNpoints());

/*		for (int i_eig_dof = 0; i_eig_dof < D->GetRows(); i_eig_dof++)
		{
			for (int j_eig_dof = 0; j_eig_dof < D->GetColumns(); j_eig_dof++)
			{
				D_Matrix(i_eig_dof,j_eig_dof) = (*D)(i_eig_dof,j_eig_dof);
				D_adv_Matrix(i_eig_dof,j_eig_dof) = (*D_adv)(i_eig_dof,j_eig_dof);
				D_no_adv_Matrix(i_eig_dof,j_eig_dof) = (*D_no_adv)(i_eig_dof,j_eig_dof);    
			}
		}

		for (int i_eig_dof = 0; i_eig_dof < Dbnd->GetRows(); i_eig_dof++)
		{
			for (int j_eig_dof = 0; j_eig_dof < Dbnd->GetColumns(); j_eig_dof++)
			{
				Dbnd_Matrix(i_eig_dof,j_eig_dof) = (*Dbnd)(i_eig_dof,j_eig_dof); 
			}
		}

		for (int i_eig_dof = 0; i_eig_dof < Dint->GetRows(); i_eig_dof++)
		{
			for (int j_eig_dof = 0; j_eig_dof < Dint->GetColumns(); j_eig_dof++)
			{
				Dint_Matrix(i_eig_dof,j_eig_dof) = (*Dint)(i_eig_dof,j_eig_dof); 
			}
		}

		for (int i_eig_dof = 0; i_eig_dof < B->GetRows(); i_eig_dof++)
		{
			for (int j_eig_dof = 0; j_eig_dof < B->GetColumns(); j_eig_dof++)
			{
				B_Matrix(i_eig_dof,j_eig_dof) = (*B)(i_eig_dof,j_eig_dof);
				B_adv_Matrix(i_eig_dof,j_eig_dof) = (*B_adv)(i_eig_dof,j_eig_dof); 
				B_no_adv_Matrix(i_eig_dof,j_eig_dof) = (*B_no_adv)(i_eig_dof,j_eig_dof);  
			}
		}

		for (int i_eig_dof = 0; i_eig_dof < C->GetRows(); i_eig_dof++)
		{
			for (int j_eig_dof = 0; j_eig_dof < C->GetColumns(); j_eig_dof++)
			{
				C_Matrix(i_eig_dof,j_eig_dof) = (*C)(i_eig_dof,j_eig_dof);
				C_adv_Matrix(i_eig_dof,j_eig_dof) = (*C_adv)(i_eig_dof,j_eig_dof);
				C_no_adv_Matrix(i_eig_dof,j_eig_dof) = (*C_no_adv)(i_eig_dof,j_eig_dof);    
			}
		}

		for (int i_eig_dof = 0; i_eig_dof < Ah->GetRows() - 1; i_eig_dof++)
		{
			for (int j_eig_dof = 0; j_eig_dof < Ah->GetColumns() - 1; j_eig_dof++)
			{
				A_Matrix(i_eig_dof,j_eig_dof) = (*Ah)(i_eig_dof,j_eig_dof);
				A_adv_Matrix(i_eig_dof,j_eig_dof) = (*A_adv)(i_eig_dof,j_eig_dof);
				A_no_adv_Matrix(i_eig_dof,j_eig_dof) = (*A_no_adv)(i_eig_dof,j_eig_dof);    
			}
		}
*/




//		cout << D_no_adv_Matrix << endl;
		// would also need to keep the advection and convection parts seperate

//		A_all.block(n*nsize_bndry[0], n*nsize_bndry[0], nsize_bndry[0], nsize_bndry[0]) = *Ah;  // can I do that directly ? no, probably not
//		A_all.block(n*nsize_bndry[0], n*nsize_bndry[0], nsize_bndry[0], nsize_bndry[0]) = A_Matrix;
//		A_adv_all.block(n*nsize_bndry[0], n*nsize_bndry[0], nsize_bndry[0], nsize_bndry[0]) = A_adv_Matrix;
//		A_no_adv_all.block(n*nsize_bndry[0], n*nsize_bndry[0], nsize_bndry[0], nsize_bndry[0]) = A_no_adv_Matrix;

//		B_all.block(n*nsize_bndry[0], n*nsize_int[0], nsize_bndry[0], nsize_int[0]) = B_Matrix;
//		B_adv_all.block(n*nsize_bndry[0], n*nsize_int[0], nsize_bndry[0], nsize_int[0]) = B_adv_Matrix;
//		B_no_adv_all.block(n*nsize_bndry[0], n*nsize_int[0], nsize_bndry[0], nsize_int[0]) = B_no_adv_Matrix;

//		C_all.block(n*nsize_bndry[0], n*nsize_int[0], nsize_bndry[0], nsize_int[0]) = C_Matrix;
//		C_adv_all.block(n*nsize_bndry[0], n*nsize_int[0], nsize_bndry[0], nsize_int[0]) = C_adv_Matrix;
//		C_no_adv_all.block(n*nsize_bndry[0], n*nsize_int[0], nsize_bndry[0], nsize_int[0]) = C_no_adv_Matrix;

//		D_all.block(n*nsize_int[0], n*nsize_int[0], nsize_int[0], nsize_int[0]) = D_Matrix;
//		D_adv_all.block(n*nsize_int[0], n*nsize_int[0], nsize_int[0], nsize_int[0]) = D_adv_Matrix;
//		D_no_adv_all.block(n*nsize_int[0], n*nsize_int[0], nsize_int[0], nsize_int[0]) = D_no_adv_Matrix;

//		Dbnd_all.block(n*nsize_p[0], n*nsize_bndry[0], nsize_p[0], nsize_bndry[0]) = Dbnd_Matrix;
//		Dint_all.block(n*nsize_p[0], n*nsize_int[0], nsize_p[0], nsize_int[0]) = Dint_Matrix; 


/*		sing_B[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = B_elem[curr_elem,:,:]
		sing_Btilde[curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = C_elem[curr_elem,:,:]
		sing_C[curr_elem*nsize_int:curr_elem*nsize_int+nsize_int, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = D_elem[curr_elem,:,:]
		sing_Dbnd[curr_elem*nsize_p:curr_elem*nsize_p+nsize_p, curr_elem*nsize_bndry:curr_elem*nsize_bndry+nsize_bndry] = Dbnd_elem[curr_elem,:,:]
		sing_Dint[curr_elem*nsize_p:curr_elem*nsize_p+nsize_p, curr_elem*nsize_int:curr_elem*nsize_int+nsize_int] = Dint_elem[curr_elem,:,:] */


			// need my mats for projection purposes: think about: have no splitting bnd/p/int could be beneficial for number of affine terms
		/*	Eigen::MatrixXd D_all( nsize_int*nel , nsize_int*nel );
			Eigen::MatrixXd Dbnd_all( nsize_p*nel , nsize_bndry*nel );
			Eigen::MatrixXd Dint_all( nsize_p*nel , nsize_int*nel );
			Eigen::MatrixXd B_all( nsize_bndry*nel , nsize_int*nel );
			Eigen::MatrixXd C_all( nsize_bndry*nel , nsize_int*nel );
			Eigen::MatrixXd A_all( nsize_bndry*nel , nsize_bndry*nel );*/




            for(int n1 = 0; n1 < D->GetColumns(); ++n1)
            {
                for(int n2 = 0; n2 < D->GetRows(); ++n2)
                {
                    
//		    cout << "(*D)(n1,n2) " << (*D)(n1,n2) << std::endl ;

                }
            }

	




                D->Invert();
                (*B) = (*B)*(*D);
                
                
                // perform (*Ah) = (*Ah) - (*B)*(*C) but since size of
                // Ah is larger than (*B)*(*C) easier to call blas
                // directly
                Blas::Dgemm('N','T', B->GetRows(), C->GetRows(), 
                            B->GetColumns(), -1.0, B->GetRawPtr(),
                            B->GetRows(), C->GetRawPtr(), 
                            C->GetRows(), 1.0, 
                            Ah->GetRawPtr(), Ah->GetRows());
            }  



            mat.m_BCinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
            mat.m_Btilde->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,C));
            mat.m_Cinv->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,D));
            mat.m_D_bnd->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dbnd));
            mat.m_D_int->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, Dint));
            
            // Do matrix manipulations and get final set of block matries    
            // reset boundary to put mean mode into boundary system. 
            
            DNekMatSharedPtr Cinv,BCinv,Btilde; 
            DNekMat  DintCinvDTint, BCinvDTint_m_DTbnd, DintCinvBTtilde_m_Dbnd;
            
            Cinv   = D;
            BCinv  = B;  
            Btilde = C; 
            
            DintCinvDTint      = (*Dint)*(*Cinv)*Transpose(*Dint);
            BCinvDTint_m_DTbnd = (*BCinv)*Transpose(*Dint) - Transpose(*Dbnd);
            
            // This could be transpose of BCinvDint in some cases
            DintCinvBTtilde_m_Dbnd = (*Dint)*(*Cinv)*Transpose(*Btilde) - (*Dbnd); 
            
            // Set up final set of matrices. 
            DNekMatSharedPtr Bh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_bndry_p1[n],nsize_p_m1[n],zero);
            DNekMatSharedPtr Ch = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n],nsize_bndry_p1[n],zero);
            DNekMatSharedPtr Dh = MemoryManager<DNekMat>::AllocateSharedPtr(nsize_p_m1[n], nsize_p_m1[n],zero);
            Array<OneD, NekDouble> Dh_data = Dh->GetPtr();
            
            // Copy matrices into final structures. 
            int n1,n2;
            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(i = 0; i < psize-1; ++i)
                {
                    for(n2 = 0; n2 < nz_loc; ++n2)
                    {
                        for(j = 0; j < psize-1; ++j)
                        {
                            //(*Dh)(i+n1*(psize-1),j+n2*(psize-1)) = 
                            //-DintCinvDTint(i+1+n1*psize,j+1+n2*psize);
                            Dh_data[(i+n1*(psize-1)) + (j+n2*(psize-1))*Dh->GetRows()] = 
                            -DintCinvDTint(i+1+n1*psize,j+1+n2*psize);
                        }
                    }                    
                }
            }
            
            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(i = 0; i < nsize_bndry_p1[n]-nz_loc; ++i)
                {
                    (*Ah)(i,nsize_bndry_p1[n]-nz_loc+n1) = BCinvDTint_m_DTbnd(i,n1*psize);
                    (*Ah)(nsize_bndry_p1[n]-nz_loc+n1,i) = DintCinvBTtilde_m_Dbnd(n1*psize,i);
                }
            }
            
            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(n2 = 0; n2 < nz_loc; ++n2)
                {
                    (*Ah)(nsize_bndry_p1[n]-nz_loc+n1,nsize_bndry_p1[n]-nz_loc+n2) = 
                    -DintCinvDTint(n1*psize,n2*psize);
                }
            }
            
            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(j = 0; j < psize-1; ++j)
                {
                    for(i = 0; i < nsize_bndry_p1[n]-nz_loc; ++i)
                    {
                        (*Bh)(i,j+n1*(psize-1)) = BCinvDTint_m_DTbnd(i,j+1+n1*psize);
                        (*Ch)(j+n1*(psize-1),i) = DintCinvBTtilde_m_Dbnd(j+1+n1*psize,i);
                    }
                }
            }
            
            for(n1 = 0; n1 < nz_loc; ++n1)
            {
                for(n2 = 0; n2 < nz_loc; ++n2)
                {
                    for(j = 0; j < psize-1; ++j)
                    {
                        (*Bh)(nsize_bndry_p1[n]-nz_loc+n1,j+n2*(psize-1)) = -DintCinvDTint(n1*psize,j+1+n2*psize);
                        (*Ch)(j+n2*(psize-1),nsize_bndry_p1[n]-nz_loc+n1) = -DintCinvDTint(j+1+n2*psize,n1*psize);
                    }
                }
            }
            
            // Do static condensation
            Dh->Invert();
            (*Bh) = (*Bh)*(*Dh);
            //(*Ah) = (*Ah) - (*Bh)*(*Ch);
            Blas::Dgemm('N','N', Bh->GetRows(), Ch->GetColumns(), Bh->GetColumns(), -1.0,
                        Bh->GetRawPtr(), Bh->GetRows(), Ch->GetRawPtr(), Ch->GetRows(), 
                        1.0, Ah->GetRawPtr(), Ah->GetRows());
            
            
	    for(int n1 = 0; n1 < Ah->GetColumns(); ++n1)
            {
                for(int n2 = 0; n2 < Ah->GetRows(); ++n2)
                {
                    
  //                  cout << "D->GetRows() " << D->GetRows() << std::endl ;
  //                  cout << "D->GetColumns() " << D->GetColumns() << std::endl ;
//		    if (n == 5) { cout << "(*Ah)(n1,n2) " << (*Ah)(n1,n2) << std::endl ; }

                }
            }


            // Set matrices for later inversion. Probably do not need to be 
            // attached to class
            pAh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ah));
            pBh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Bh));
            pCh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Ch));
            pDh->SetBlock(n,n,loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Dh));    
        }
        timer.Stop();
//        cout << "Matrix Setup Costs: " << timer.TimePerTest(1) << endl;
        
	// end of the loop over the spectral elements
            
/*	    RB_A = A_all; // das ist ein bisschen doppelt, die RB mats sind im .h header definiert, die _all mats hier
	    RB_A_adv = A_adv_all;
	    RB_A_no_adv = A_no_adv_all;
	    RB_B = B_all;
	    RB_B_adv = B_adv_all;
	    RB_B_no_adv = B_no_adv_all;
	    RB_C = C_all;
	    RB_C_adv = C_adv_all;
	    RB_C_no_adv = C_no_adv_all;
	    RB_D = D_all;
	    RB_D_adv = D_adv_all;
	    RB_D_no_adv = D_no_adv_all;
	    RB_Dbnd = Dbnd_all;
	    RB_Dint = Dint_all; 
*/
//	cout << D_all << endl;


//        cout << RB_D << endl;
//        cout << RB_D_adv << endl;
//        cout << RB_D_no_adv << endl;
        
        timer.Start();
        // Set up global coupled boundary solver. 
        // This is a key to define the solution matrix type
        // currently we are giving it a argument of eLInearAdvectionReaction 
        // since this then makes the matrix storage of type eFull
        MultiRegions::GlobalLinSysKey key(StdRegions::eLinearAdvectionReaction,locToGloMap);
        mat.m_CoupledBndSys = MemoryManager<MultiRegions::GlobalLinSysDirectStaticCond>::AllocateSharedPtr(key,m_fields[0],pAh,pBh,pCh,pDh,locToGloMap);
        mat.m_CoupledBndSys->Initialise(locToGloMap);
        timer.Stop();
//        cout << "Multilevel condensation: " << timer.TimePerTest(1) << endl;
    }


    void CoupledLinearNS_TT::Solve(void)
    {

        const unsigned int ncmpt = m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > forcing_phys(ncmpt);
        Array<OneD, Array<OneD, NekDouble> > forcing     (ncmpt);

        for(int i = 0; i < ncmpt; ++i)
        {
            forcing_phys[i] = Array<OneD, NekDouble> (m_fields[m_velocity[0]]->GetNpoints(), 0.0);
            forcing[i]      = Array<OneD, NekDouble> (m_fields[m_velocity[0]]->GetNcoeffs(),0.0);
        }

        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            const NekDouble time = 0;
            (*x)->Apply(m_fields, forcing_phys, forcing_phys, time);
        }
        for (unsigned int i = 0; i < ncmpt; ++i)
        {
            bool waveSpace = m_fields[m_velocity[i]]->GetWaveSpace();
            m_fields[i]->SetWaveSpace(true);
            m_fields[i]->IProductWRTBase(forcing_phys[i], forcing[i]);
            m_fields[i]->SetWaveSpace(waveSpace);
        }

        SolveLinearNS(forcing);

    }

    void CoupledLinearNS_TT::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing)
    {
        int i,n;
        Array<OneD,  MultiRegions::ExpListSharedPtr> vel_fields(m_velocity.num_elements());
        Array<OneD, Array<OneD, NekDouble> > force(m_velocity.num_elements());
        
        if(m_HomogeneousType == eHomogeneous1D)
        {
            int ncoeffsplane = m_fields[m_velocity[0]]->GetPlane(0)->GetNcoeffs();
            for(n = 0; n < m_npointsZ/2; ++n)
            {
                // Get the a Fourier mode of velocity and forcing components. 
                for(i = 0; i < m_velocity.num_elements(); ++i)
                {
                    vel_fields[i] = m_fields[m_velocity[i]]->GetPlane(2*n);
                    // Note this needs to correlate with how we pass forcing
                    force[i] = forcing[i] + 2*n*ncoeffsplane;
                }
                
                SolveLinearNS(force,vel_fields,m_pressure->GetPlane(2*n),n);
            }
            for(i = 0; i < m_velocity.num_elements(); ++i)
            {
                m_fields[m_velocity[i]]->SetPhysState(false);
            }
            m_pressure->SetPhysState(false);
        }
        else
        {
            for(i = 0; i < m_velocity.num_elements(); ++i)
            {
                vel_fields[i] = m_fields[m_velocity[i]];
                // Note this needs to correlate with how we pass forcing
                force[i] = forcing[i];
            }
            SolveLinearNS(force,vel_fields,m_pressure);
        }
    }
    
    void CoupledLinearNS_TT::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing,  Array<OneD, MultiRegions::ExpListSharedPtr> &fields, MultiRegions::ExpListSharedPtr &pressure,  const int mode)
    {
        int i,j,k,n,eid,cnt,cnt1;
        int nbnd,nint,offset;
        int nvel = m_velocity.num_elements();
        int nel  = fields[0]->GetNumElmts();
        Array<OneD, unsigned int> bmap, imap; 
        
        Array<OneD, NekDouble > f_bnd(m_mat[mode].m_BCinv->GetRows());
        NekVector< NekDouble  > F_bnd(f_bnd.num_elements(), f_bnd, eWrapper);
        Array<OneD, NekDouble > f_int(m_mat[mode].m_BCinv->GetColumns());
        NekVector< NekDouble  > F_int(f_int.num_elements(),f_int, eWrapper);
        
        int nz_loc;
        int  nplanecoeffs = fields[m_velocity[0]]->GetNcoeffs();// this is fine since we pass the nplane coeff data. 
        
        if(mode) // Homogeneous mode flag
        {
            nz_loc = 2;
        }
        else
        {
            if(m_singleMode)
            {
                nz_loc = 2;
            }
            else
            {
                nz_loc = 1;
                if(m_HomogeneousType == eHomogeneous1D)
                {
                    // Zero fields to set complex mode to zero;
                    for(i = 0; i < fields.num_elements(); ++i)
                    {
                        Vmath::Zero(2*fields[i]->GetNcoeffs(),fields[i]->UpdateCoeffs(),1);
                    }
                    Vmath::Zero(2*pressure->GetNcoeffs(),pressure->UpdateCoeffs(),1);
                }
            }
        }

	// do need to write the forcing vectors
  /*        ofstream myfilef0 ("forcing0.txt");
	  if (myfilef0.is_open())
	  {
		for (int counter_row = 0; counter_row < forcing[0].num_elements(); ++counter_row)
		{
				myfilef0 << std::setprecision(17) << forcing[0][counter_row] << " ";
//				cout << "writing " << collected_trajectory_f1[counter_timesteps][counter_velocity] << " ";
			myfilef0 << "\n";
		}
                myfilef0.close();
	  }
	  else cout << "Unable to open file";

          ofstream myfilef1 ("forcing1.txt");
	  if (myfilef1.is_open())
	  {
		for (int counter_row = 0; counter_row < forcing[1].num_elements(); ++counter_row)
		{
				myfilef1 << std::setprecision(17) << forcing[1][counter_row] << " ";
//				cout << "writing " << collected_trajectory_f1[counter_timesteps][counter_velocity] << " ";
			myfilef1 << "\n";
		}
                myfilef1.close();
	  }
	  else cout << "Unable to open file";
*/        

        // Assemble f_bnd and f_int
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            eid = fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            fields[m_velocity[0]]->GetExp(eid)->GetBoundaryMap(bmap);
            fields[m_velocity[0]]->GetExp(eid)->GetInteriorMap(imap);
            nbnd   = bmap.num_elements();
            nint   = imap.num_elements();
            offset = fields[m_velocity[0]]->GetCoeff_Offset(eid);
            
            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(n = 0; n < nz_loc; ++n)
                {
                    for(k = 0; k < nbnd; ++k)
                    {
                        f_bnd[cnt+k] = forcing[j][n*nplanecoeffs + 
                        offset+bmap[k]];
                    }
                    for(k = 0; k < nint; ++k)
                    {
                        f_int[cnt1+k] = forcing[j][n*nplanecoeffs + 
                        offset+imap[k]];
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
        }
        
        Array<OneD, NekDouble > f_p(m_mat[mode].m_D_int->GetRows());
        NekVector<  NekDouble > F_p(f_p.num_elements(),f_p,eWrapper);
        NekVector<  NekDouble > F_p_tmp(m_mat[mode].m_Cinv->GetRows());
        
        // fbnd does not currently hold the pressure mean
        F_bnd = F_bnd - (*m_mat[mode].m_BCinv)*F_int;
        F_p_tmp = (*m_mat[mode].m_Cinv)*F_int;
        F_p = (*m_mat[mode].m_D_int) * F_p_tmp;
        
        // construct inner forcing 
        Array<OneD, NekDouble > bnd   (m_locToGloMap[mode]->GetNumGlobalCoeffs(),0.0);
        Array<OneD, NekDouble > fh_bnd(m_locToGloMap[mode]->GetNumGlobalCoeffs(),0.0);
        
        const Array<OneD,const int>& loctoglomap
        = m_locToGloMap[mode]->GetLocalToGlobalMap();
        const Array<OneD,const NekDouble>& loctoglosign
        = m_locToGloMap[mode]->GetLocalToGlobalSign();
        
/*          ofstream myfile1LocGloMap ("LocGloMap.txt");
	  if (myfile1LocGloMap.is_open())
	  {
		for (int counter = 0; counter < loctoglomap.num_elements(); ++counter)
		{
			myfile1LocGloMap << std::setprecision(17) << loctoglomap[counter] << " ";
			myfile1LocGloMap << "\n";
		}
                myfile1LocGloMap.close();
	  }
	  else cout << "Unable to open file";

          ofstream myfile2LocGloSign ("LocGloSign.txt");
	  if (myfile2LocGloSign.is_open())
	  {
		for (int counter = 0; counter < loctoglosign.num_elements(); ++counter)
		{
			myfile2LocGloSign << std::setprecision(17) << loctoglosign[counter] << " ";
			myfile2LocGloSign << "\n";
		}
                myfile2LocGloSign.close();
	  }
	  else cout << "Unable to open file";

        const Array<OneD,const int>& loctoglobndmap
        = m_locToGloMap[mode]->GetLocalToGlobalBndMap();
        const Array<OneD,const NekDouble>& loctoglobndsign
        = m_locToGloMap[mode]->GetLocalToGlobalBndSign();

          ofstream myfile1LocGloBndMap ("LocGloBndMap.txt");
	  if (myfile1LocGloBndMap.is_open())
	  {
		for (int counter = 0; counter < loctoglobndmap.num_elements(); ++counter)
		{
			myfile1LocGloBndMap << std::setprecision(17) << loctoglobndmap[counter] << " ";
			myfile1LocGloBndMap << "\n";
		}
                myfile1LocGloBndMap.close();
	  }
	  else cout << "Unable to open file";

          ofstream myfile2LocGloBndSign ("LocGloBndSign.txt");
	  if (myfile2LocGloBndSign.is_open())
	  {
		for (int counter = 0; counter < loctoglobndsign.num_elements(); ++counter)
		{
			myfile2LocGloBndSign << std::setprecision(17) << loctoglobndsign[counter] << " ";
			myfile2LocGloBndSign << "\n";
		}
                myfile2LocGloBndSign.close();
	  }
	  else cout << "Unable to open file";
*/
   /*         Array<OneD, Array<OneD, NekDouble> > glo_dof_to_phys_collector(number_of_global_dofs);
            Array<OneD, Array<OneD, NekDouble> > phys_to_glo_dof_collector(m_fields[0]->GetNpoints());
            Array<OneD, Array<OneD, NekDouble> > glo_dof_to_phys_collector_p(number_of_global_dofs);
            Array<OneD, Array<OneD, NekDouble> > phys_to_glo_dof_collector_p(pressure_field->GetNpoints());

	    // get a global_dof <-> phys map
	    for (int i_global_dofs = 0; i_global_dofs < number_of_global_dofs; i_global_dofs++)
	    {
              Array<OneD, NekDouble> init_global_test_field(number_of_global_dofs,0.0);
              Array<OneD, NekDouble> test_field(m_fields[0]->GetNpoints(),0.0);
              init_global_test_field[i_global_dofs] = 1.0;
              m_fields[0]->GlobalToLocal(init_global_test_field, local_test_field);
              m_fields[0]->BwdTrans(local_test_field, test_field);
	      glo_dof_to_phys_collector[i_global_dofs] = test_field;
	    } */

        offset = cnt = 0; 
        for(i = 0; i < nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nbnd = nz_loc*fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    fh_bnd[loctoglomap[offset+j*nbnd+k]] += 
                    loctoglosign[offset+j*nbnd+k]*f_bnd[cnt+k];
                }
                cnt += nbnd;
            }
            
            nint    = pressure->GetExp(eid)->GetNcoeffs();
            offset += nvel*nbnd + nint*nz_loc; 
        }
        
        offset = cnt1 = 0; 
        for(i = 0; i <  nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nbnd = nz_loc*fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            nint = pressure->GetExp(eid)->GetNcoeffs(); 
            
            for(n = 0; n < nz_loc; ++n)
            {
                for(j = 0; j < nint; ++j)
                {
                    fh_bnd[loctoglomap[offset + nvel*nbnd + n*nint+j]] = f_p[cnt1+j];
                }
                cnt1   += nint;
            }
            offset += nvel*nbnd + nz_loc*nint; 
        }
        
        //  Set Weak BC into f_bnd and Dirichlet Dofs in bnd
        const Array<OneD,const int>& bndmap
        = m_locToGloMap[mode]->GetBndCondCoeffsToGlobalCoeffsMap();
        
/*          ofstream myfilebndmap ("BndCondCoeffsToGlobalCoeffsMap.txt");
	  if (myfilebndmap.is_open())
	  {
		for (int counter = 0; counter < bndmap.num_elements(); ++counter)
		{
			{
				myfilebndmap << std::setprecision(17) << bndmap[counter] << " ";
			}
		}
                myfilebndmap.close();
	  }
	  else cout << "Unable to open file";   
*/
        // Forcing function with weak boundary conditions and
        // Dirichlet conditions
        int bndcnt=0;
        
        for(k = 0; k < nvel; ++k)
        {



            const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConds = fields[k]->GetBndConditions();
            Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp;
            if(m_HomogeneousType == eHomogeneous1D) 
            {
                bndCondExp = m_fields[k]->GetPlane(2*mode)->GetBndCondExpansions();
            }
            else
            {
                bndCondExp = m_fields[k]->GetBndCondExpansions();
            }
            
/*     		ofstream myfile_bnd_cond_elem ("bnd_cond_elem.txt");
		if (myfile_bnd_cond_elem.is_open())
	  	{
			myfile_bnd_cond_elem << std::setprecision(17) << bndCondExp.num_elements() << "\n";
                	myfile_bnd_cond_elem.close();
	  	}
	  	else cout << "Unable to open file"; 
*/
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                const Array<OneD, const NekDouble > bndCondCoeffs = bndCondExp[i]->GetCoeffs();

		std::stringstream sstm;
		sstm << "bndcond_k" << k << "_i_" << i << ".txt";
		std::string result = sstm.str();
		const char* rr = result.c_str();

//		std::string outname_txt = "bndcond_k" + std::to_string(k) + "_i_" + std::to_string(i) + ".txt";
//		const char* outname_t = outname_txt.c_str();
		const char* outname_t = rr;

       
/*     		ofstream myfile_t (outname_t);
		if (myfile_t.is_open())
	  	{
			for(int count_bnd = 0; count_bnd < bndCondCoeffs.num_elements(); count_bnd++)
			{
				myfile_t << std::setprecision(17) << bndCondCoeffs[count_bnd] << "\n";
			}
                	myfile_t.close();
	  	}
	  	else cout << "Unable to open file"; 
*/

//		cout << "i " << i << endl;
//		for(int count_bnd = 0; count_bnd < bndCondCoeffs.num_elements(); count_bnd++)
//		{
//			cout << "bndCondCoeffs[count_bnd] " << bndCondCoeffs[count_bnd] << endl;
//		}
 


               cnt = 0;
                for(n = 0; n < nz_loc; ++n)
                {
                    if(bndConds[i]->GetBoundaryConditionType() 
                        == SpatialDomains::eDirichlet)
                    {
                        for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                        {
                            if (m_equationType == eSteadyNavierStokes && m_initialStep == false)
                            {
                                //This condition set all the Dirichlet BC at 0 after
                                //the initial step of the Newton method
                                bnd[bndmap[bndcnt++]] = 0;
                            }
                            else
                            {
                                bnd[bndmap[bndcnt++]] = bndCondCoeffs[cnt++];
                            }
                        }
                    }
                    else
                    {                    
                        for(j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                        {
                            fh_bnd[bndmap[bndcnt++]]
                            += bndCondCoeffs[cnt++];
                        }
                    }
                }
            }
        }
        
//	cout <<	"m_locToGloMap[mode]->GetNumGlobalDirBndCoeffs() " << m_locToGloMap[mode]->GetNumGlobalDirBndCoeffs() << endl;
//	cout <<	"m_locToGloMap[mode]->GetNumGlobalCoeffs() " << m_locToGloMap[mode]->GetNumGlobalCoeffs() << endl;

/*          ofstream myfiledirbnd ("NumGlobalDirBndCoeffs.txt");
	  if (myfiledirbnd.is_open())
	  {
		myfiledirbnd << std::setprecision(17) << m_locToGloMap[mode]->GetNumGlobalDirBndCoeffs();
                myfiledirbnd.close();
	  }
	  else cout << "Unable to open file";   
*/

//	cout << "fh_bnd.num_elements() " << fh_bnd.num_elements() << std::endl;   
//	cout << "bnd.num_elements() " << bnd.num_elements() << std::endl;  
	
	for(i = 0; i <  fh_bnd.num_elements(); ++i)
	{
//		cout << "bnd[i] " << i << " " << bnd[i] << std::endl;   
	}
        m_mat[mode].m_CoupledBndSys->Solve(fh_bnd,bnd,m_locToGloMap[mode]);                 // actual multi-static-condensed solve here

/*	for(i = 0; i <  fh_bnd.num_elements(); ++i)
	{
		cout << "bnd[i] " << i << " " << bnd[i] << std::endl;   
	}
        
*/

        // unpack pressure and velocity boundary systems. 
        offset = cnt = 0; 
        int totpcoeffs = pressure->GetNcoeffs();
        Array<OneD, NekDouble> p_coeffs = pressure->UpdateCoeffs();
        for(i = 0; i <  nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nbnd = nz_loc*fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            nint = pressure->GetExp(eid)->GetNcoeffs(); 
            
            for(j = 0; j < nvel; ++j)
            {
                for(k = 0; k < nbnd; ++k)
                {
                    f_bnd[cnt+k] = loctoglosign[offset+j*nbnd+k]*bnd[loctoglomap[offset + j*nbnd + k]];
/*		    cout << "writing " << loctoglosign[offset+j*nbnd+k]*bnd[loctoglomap[offset + j*nbnd + k]] << " to " << cnt+k << std::endl;
		    cout << "loctoglosign[offset+j*nbnd+k] " << loctoglosign[offset+j*nbnd+k] << " loctoglomap[offset + j*nbnd + k] " << loctoglomap[offset + j*nbnd + k] << std::endl;
		    cout << "nvel " << nvel << std::endl;
		    cout << "nint " << nint << std::endl;
		    cout << "nbnd " << nbnd << std::endl;
		    cout << "nz_loc " << nz_loc << std::endl; */
                }
                cnt += nbnd;
            }
            offset += nvel*nbnd + nint*nz_loc;
        }
        
/*	for(i = 0; i <  f_bnd.num_elements(); ++i)
	{
		cout << "f_bnd[i] " << i << " " << f_bnd[i] << std::endl;   
	}
*/

        pressure->SetPhysState(false);
        
        offset = cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i)
        {
            eid  = fields[0]->GetOffset_Elmt_Id(i);
            nint = pressure->GetExp(eid)->GetNcoeffs(); 
            nbnd = fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            cnt1 = pressure->GetCoeff_Offset(eid);

//	    cout << "nint " << nint << std::endl;
//	    cout << "nbnd " << nbnd << std::endl;
//	    cout << "cnt1 " << cnt1 << std::endl;
            
            for(n = 0; n < nz_loc; ++n)
            {
                for(j = 0; j < nint; ++j)
                {
                    p_coeffs[n*totpcoeffs + cnt1+j] = 
                    f_p[cnt+j] = bnd[loctoglomap[offset + 
                    (nvel*nz_loc)*nbnd + 
                    n*nint + j]];
                }
                cnt += nint;
            }
            offset += (nvel*nbnd + nint)*nz_loc;
        }
        
/*	for(i = 0; i <  f_p.num_elements(); ++i)
	{
		cout << "f_p[i] " << i << " " << f_p[i] << std::endl;   
	}
*/

        // Back solve first level of static condensation for interior
        // velocity space and store in F_int
        F_int = F_int + Transpose(*m_mat[mode].m_D_int)*F_p
        - Transpose(*m_mat[mode].m_Btilde)*F_bnd;
        F_int = (*m_mat[mode].m_Cinv)*F_int;

	// --- writing the fields in loc format from here on, so what i want is also the f_bnd, f_p, f_int
//	curr_f_bnd = f_bnd;


	curr_f_bnd = Eigen::VectorXd::Zero(f_bnd.num_elements());
	for (int i_phys_dof = 0; i_phys_dof < f_bnd.num_elements(); i_phys_dof++)
	{
		curr_f_bnd(i_phys_dof) = f_bnd[i_phys_dof]; 
	}
	curr_f_p = Eigen::VectorXd::Zero(f_p.num_elements()); 
	for (int i_phys_dof = 0; i_phys_dof < f_p.num_elements(); i_phys_dof++)
	{
		curr_f_p(i_phys_dof) = f_p[i_phys_dof]; 
	}
	curr_f_int = Eigen::VectorXd::Zero(f_int.num_elements()); 
	for (int i_phys_dof = 0; i_phys_dof < f_int.num_elements(); i_phys_dof++)
	{
		curr_f_int(i_phys_dof) = f_int[i_phys_dof]; 
	}

        
        // Unpack solution from Bnd and F_int to v_coeffs 
        cnt = cnt1 = 0;
        for(i = 0; i < nel; ++i) // loop over elements
        {
            eid  = fields[m_velocity[0]]->GetOffset_Elmt_Id(i);
            fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
            fields[0]->GetExp(eid)->GetInteriorMap(imap);
            nbnd   = bmap.num_elements();
            nint   = imap.num_elements();
            offset = fields[0]->GetCoeff_Offset(eid);
            
	//    cout << "offset " << offset << std::endl;
	//    cout << "nint " << nint << std::endl;

            for(j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(n = 0; n < nz_loc; ++n)
                {
                    for(k = 0; k < nbnd; ++k)
                    {
                        fields[j]->SetCoeff(n*nplanecoeffs + 
                        offset+bmap[k],
                        f_bnd[cnt+k]);
                    }
                    
                    for(k = 0; k < nint; ++k)
                    {
                        fields[j]->SetCoeff(n*nplanecoeffs + 
                        offset+imap[k],
                        f_int[cnt1+k]);
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
        }
        
        for(j = 0; j < nvel; ++j) 
        {
            fields[j]->SetPhysState(false);
        }
    }

    NekDouble CoupledLinearNS_TT::Get_m_kinvis(void)
    {
	return m_kinvis;
    }

    void CoupledLinearNS_TT::Set_m_kinvis(NekDouble input)
    {
	m_kinvis = input;
    }

    Eigen::MatrixXd CoupledLinearNS_TT::Get_no_advection_matrix(void)
    {
/*	Eigen::MatrixXd no_adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	switch(globally_connected) {
		case 0:
			no_adv_matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = MtM * RB_A_no_adv;
			no_adv_matrix.block(0, RB_A.cols(), RB_Dbnd.cols(), RB_Dbnd.rows()) = -MtM * RB_Dbnd.transpose();
			no_adv_matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = MtM * RB_B_no_adv;
			no_adv_matrix.block(RB_A.rows(), 0, RB_Dbnd.rows(), RB_Dbnd.cols()) = -RB_Dbnd;
			no_adv_matrix.block(RB_A.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C_no_adv.transpose();
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols(), RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_no_adv;
			break;
		case 1:
			no_adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * RB_A_no_adv * Mtrafo;
			no_adv_matrix.block(0, nBndDofs, nBndDofs, RB_Dbnd.rows()) = -Mtrafo.transpose() * RB_Dbnd.transpose();
			no_adv_matrix.block(0, nBndDofs + RB_Dbnd.rows(), nBndDofs, RB_B.cols()) = Mtrafo.transpose() * RB_B_no_adv;
			no_adv_matrix.block(nBndDofs, 0, RB_Dbnd.rows(), nBndDofs) = -RB_Dbnd * Mtrafo;
			no_adv_matrix.block(nBndDofs, nBndDofs + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			no_adv_matrix.block(nBndDofs + RB_Dbnd.rows(), 0, RB_C.cols(), nBndDofs) = RB_C_no_adv.transpose() * Mtrafo;
			no_adv_matrix.block(nBndDofs + RB_Dbnd.rows(), nBndDofs, RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			no_adv_matrix.block(nBndDofs + RB_Dbnd.rows(), nBndDofs + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_no_adv;
			break;
		case 2:
			no_adv_matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = RB_A_no_adv;
			no_adv_matrix.block(0, RB_A.cols(), RB_Dbnd.cols(), RB_Dbnd.rows()) = -RB_Dbnd.transpose();
			no_adv_matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = RB_B_no_adv;
			no_adv_matrix.block(RB_A.rows(), 0, RB_Dbnd.rows(), RB_Dbnd.cols()) = -RB_Dbnd;
			no_adv_matrix.block(RB_A.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C_no_adv.transpose();
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols(), RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_no_adv;
			break;
	}
	return no_adv_matrix; */
    }

    Eigen::MatrixXd CoupledLinearNS_TT::Get_no_advection_matrix_pressure(void)
    {
/*	Eigen::MatrixXd no_adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	switch(globally_connected) {
		case 0:
			no_adv_matrix.block(0, RB_A.cols(), RB_Dbnd.cols(), RB_Dbnd.rows()) = -MtM * RB_Dbnd.transpose();
			no_adv_matrix.block(RB_A.rows(), 0, RB_Dbnd.rows(), RB_Dbnd.cols()) = -RB_Dbnd;
			no_adv_matrix.block(RB_A.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols(), RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			break;
		case 1:
			no_adv_matrix.block(0, nBndDofs, nBndDofs, RB_Dbnd.rows()) = -Mtrafo.transpose() * RB_Dbnd.transpose();
			no_adv_matrix.block(nBndDofs, 0, RB_Dbnd.rows(), nBndDofs) = -RB_Dbnd * Mtrafo;
			no_adv_matrix.block(nBndDofs, nBndDofs + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			no_adv_matrix.block(nBndDofs + RB_Dbnd.rows(), nBndDofs, RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			break;
		case 2:
			no_adv_matrix.block(0, RB_A.cols(), RB_Dbnd.cols(), RB_Dbnd.rows()) = -RB_Dbnd.transpose();
			no_adv_matrix.block(RB_A.rows(), 0, RB_Dbnd.rows(), RB_Dbnd.cols()) = -RB_Dbnd;
			no_adv_matrix.block(RB_A.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols(), RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			break;
	}			


	return no_adv_matrix; */
			
    }

    Eigen::MatrixXd CoupledLinearNS_TT::Get_no_advection_matrix_ABCD(void)
    {
/*	Eigen::MatrixXd no_adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	switch(globally_connected) {
		case 0:
			no_adv_matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = MtM * RB_A_no_adv;
			no_adv_matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = MtM * RB_B_no_adv;
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C_no_adv.transpose();
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_no_adv;
			break;
		case 1:
			no_adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * RB_A_no_adv * Mtrafo;
			no_adv_matrix.block(0, nBndDofs + RB_Dbnd.rows(), nBndDofs, RB_B.cols()) = Mtrafo.transpose() * RB_B_no_adv;
			no_adv_matrix.block(nBndDofs + RB_Dbnd.rows(), 0, RB_C.cols(), nBndDofs) = RB_C_no_adv.transpose() * Mtrafo;
			no_adv_matrix.block(nBndDofs + RB_Dbnd.rows(), nBndDofs + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_no_adv;
			break;
		case 2:
			no_adv_matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = RB_A_no_adv;
			no_adv_matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = RB_B_no_adv;
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C_no_adv.transpose();
			no_adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_no_adv;
			break;
	}

	return no_adv_matrix; */
    }
	
    Eigen::MatrixXd CoupledLinearNS_TT::Get_advection_matrix(void)
    {
/*	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	switch(globally_connected) {
		case 0:
			adv_matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = MtM * RB_A_adv;
			adv_matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = MtM * RB_B_adv;
			adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C_adv.transpose();
			adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_adv;
			break;
		case 1:
			adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * RB_A_adv * Mtrafo;
			adv_matrix.block(0, nBndDofs + RB_Dbnd.rows(), nBndDofs, RB_B.cols()) = Mtrafo.transpose() * RB_B_adv;
			adv_matrix.block(nBndDofs + RB_Dbnd.rows(), 0, RB_C.cols(), nBndDofs) = RB_C_adv.transpose() * Mtrafo;
			adv_matrix.block(nBndDofs + RB_Dbnd.rows(), nBndDofs + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_adv;
			break;
		case 2:
			adv_matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = RB_A_adv;
			adv_matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = RB_B_adv;
			adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C_adv.transpose();
			adv_matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D_adv;
			break;
	}

	return adv_matrix; */
    }	

    Eigen::MatrixXd CoupledLinearNS_TT::Get_complete_matrix(void)
    {
/*	Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	switch(globally_connected) {
		case 0:
			matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = MtM * RB_A;
			matrix.block(0, RB_A.cols(), RB_Dbnd.cols(), RB_Dbnd.rows()) = -MtM * RB_Dbnd.transpose();
			matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = MtM * RB_B;
			matrix.block(RB_A.rows(), 0, RB_Dbnd.rows(), RB_Dbnd.cols()) = -RB_Dbnd;
			matrix.block(RB_A.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C.transpose();
			matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols(), RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D;
			break;
		case 1:
			matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * RB_A * Mtrafo;
			matrix.block(0, nBndDofs, nBndDofs, RB_Dbnd.rows()) = -Mtrafo.transpose() * RB_Dbnd.transpose();
			matrix.block(0, nBndDofs + RB_Dbnd.rows(), nBndDofs, RB_B.cols()) = Mtrafo.transpose() * RB_B;
			matrix.block(nBndDofs, 0, RB_Dbnd.rows(), nBndDofs) = -RB_Dbnd * Mtrafo;
			matrix.block(nBndDofs, nBndDofs + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			matrix.block(nBndDofs + RB_Dbnd.rows(), 0, RB_C.cols(), nBndDofs) = RB_C.transpose() * Mtrafo;
			matrix.block(nBndDofs + RB_Dbnd.rows(), nBndDofs, RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			matrix.block(nBndDofs + RB_Dbnd.rows(), nBndDofs + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D;
			break;
		case 2:
			matrix.block(0, 0, RB_A.rows(), RB_A.cols()) = RB_A;
			matrix.block(0, RB_A.cols(), RB_Dbnd.cols(), RB_Dbnd.rows()) = -RB_Dbnd.transpose();
			matrix.block(0, RB_A.cols() + RB_Dbnd.rows(), RB_B.rows(), RB_B.cols()) = RB_B;
			matrix.block(RB_A.rows(), 0, RB_Dbnd.rows(), RB_Dbnd.cols()) = -RB_Dbnd;
			matrix.block(RB_A.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_Dint.rows(), RB_Dint.cols()) = -RB_Dint;	
			matrix.block(RB_A.rows() + RB_Dbnd.rows(), 0, RB_C.cols(), RB_C.rows()) = RB_C.transpose();
			matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols(), RB_Dint.cols(), RB_Dint.rows()) = -RB_Dint.transpose();
			matrix.block(RB_A.rows() + RB_Dbnd.rows(), RB_A.cols() + RB_Dbnd.rows(), RB_D.rows(), RB_D.cols()) = RB_D;
			break;
	}			
	return matrix; */
    }

    Eigen::MatrixXd CoupledLinearNS_TT::remove_cols_and_rows(Eigen::MatrixXd the_matrix, std::set<int> elements_to_be_removed)
    {

/*	Eigen::MatrixXd sm1;
	if (1)
	{
	Eigen::MatrixXd simplified_matrix = Eigen::MatrixXd::Zero(the_matrix.rows() - elements_to_be_removed.size(), the_matrix.cols() - elements_to_be_removed.size());
	int counter_row_simplified = 0;
	int counter_col_simplified = 0;
	for (int row_index=0; row_index < the_matrix.rows(); ++row_index)
	{
		if (!elements_to_be_removed.count(row_index))
		{
			for (int col_index=0; col_index < the_matrix.cols(); ++col_index)
			{
				if (!elements_to_be_removed.count(col_index))
				{
					simplified_matrix(counter_row_simplified, counter_col_simplified) = the_matrix(row_index, col_index);
					counter_col_simplified++;
				}			
			}
			counter_col_simplified = 0;
			counter_row_simplified++;
		}		
	}
	//return simplified_matrix;
	//	cout << "simplified_matrix.rows() " << simplified_matrix.rows() << " simplified_matrix.cols() " << simplified_matrix.cols() << endl;
	sm1 = simplified_matrix;
	} 
	
	cout << "simplified_matrix computation slow norm " << sm1.norm() << endl;
*/
	// or move whole blocks in-place ??
	// need to use iterators
	Eigen::MatrixXd simplified_matrix = Eigen::MatrixXd::Zero(the_matrix.rows() - elements_to_be_removed.size(), the_matrix.cols() - elements_to_be_removed.size());
	std::set<int>::iterator set_iterator_rows;
	std::set<int>::iterator set_iterator_cols; 
	int prev_iter_row = -1; // check if zero is to be removed
	int prev_iter_col = -1; // also need to move the last block??
	int passed_rows = 0;
	int passed_cols = 0;
//	Eigen::MatrixXd sm2;
	if (1)
	{
//		Eigen::MatrixXd simplified_matrix = Eigen::MatrixXd::Zero(the_matrix.rows() - elements_to_be_removed.size(), the_matrix.cols() - elements_to_be_removed.size());
		for (set_iterator_rows = elements_to_be_removed.begin(); set_iterator_rows != elements_to_be_removed.end(); set_iterator_rows++)
		{
			if ((prev_iter_row != *set_iterator_rows) && (0 != *set_iterator_rows))
			{
				prev_iter_col = -1;
				passed_cols = 0;
				for (set_iterator_cols = elements_to_be_removed.begin(); set_iterator_cols != elements_to_be_removed.end(); set_iterator_cols++)
				{
//					cout << " output set_iterator_rows " << *set_iterator_rows << endl; // might contain many consecutive elements
//					cout << " output set_iterator_cols " << *set_iterator_cols << endl; 
					if ((prev_iter_col != *set_iterator_cols) && (0 != *set_iterator_cols))
					{
						// do move a block
						simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, *set_iterator_rows - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, *set_iterator_rows - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1);
					}
					prev_iter_col = *set_iterator_cols;
					passed_cols++;
				}
				// last block
				simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, *set_iterator_rows - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, *set_iterator_rows - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1);
			}
			prev_iter_row = *set_iterator_rows;
			passed_rows++;
		}
		prev_iter_col = -1;
		passed_cols = 0;
		for (set_iterator_cols = elements_to_be_removed.begin(); set_iterator_cols != elements_to_be_removed.end(); set_iterator_cols++)
		{
			if ((prev_iter_col != *set_iterator_cols) && (0 != *set_iterator_cols))
			{
				// do move a block
				simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, the_matrix.rows() - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, the_matrix.rows() - prev_iter_row - 1, *set_iterator_cols - prev_iter_col - 1);
			}
			prev_iter_col = *set_iterator_cols;
			passed_cols++;
		}
		// last block
		simplified_matrix.block(prev_iter_row - passed_rows + 1, prev_iter_col - passed_cols + 1, the_matrix.rows() - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1) = the_matrix.block(prev_iter_row + 1, prev_iter_col + 1, the_matrix.rows() - prev_iter_row - 1, the_matrix.cols() - prev_iter_col - 1);


//		sm2 = simplified_matrix;
	}
//	cout << "simplified_matrix computation new norm " << sm2.norm() << endl;

/*	Eigen::MatrixXd simplified_matrix = Eigen::MatrixXd::Zero(the_matrix.rows() - elements_to_be_removed.size(), the_matrix.cols() - elements_to_be_removed.size());
	Eigen::MatrixXd intermediate_simplified_matrix = Eigen::MatrixXd::Zero(the_matrix.rows(), the_matrix.cols() - elements_to_be_removed.size());
	int counter_col_simplified = 0;
	for (int col_index = 0; col_index < the_matrix.cols(); ++col_index)
	{
		if (!elements_to_be_removed.count(col_index))
		{
			intermediate_simplified_matrix.col(counter_col_simplified) = the_matrix.col(col_index);
			counter_col_simplified++;
		}
	}
	int counter_row_simplified = 0;
	for (int row_index = 0; row_index < the_matrix.cols(); ++row_index)
	{
		if (!elements_to_be_removed.count(row_index))
		{
			simplified_matrix.row(counter_row_simplified) = intermediate_simplified_matrix.row(row_index);
			counter_row_simplified++;
		}
	}
*/
//	cout << "simplified_matrix.rows() " << simplified_matrix.rows() << " simplified_matrix.cols() " << simplified_matrix.cols() << endl;
//	cout << "simplified_matrix computation better norm " << simplified_matrix.norm() << endl;
//	cout << "simplified_matrix computation difference norm " << (simplified_matrix - sm2).norm() << endl; 

	return simplified_matrix; 



    }

    Eigen::VectorXd CoupledLinearNS_TT::remove_rows(Eigen::VectorXd the_vector, std::set<int> elements_to_be_removed)
    {
	Eigen::VectorXd simplified_vector = Eigen::VectorXd::Zero(the_vector.rows() - elements_to_be_removed.size());
	int counter_row_simplified = 0;
	for (int row_index=0; row_index < the_vector.rows(); ++row_index)
	{
		if (!elements_to_be_removed.count(row_index))
		{
			simplified_vector(counter_row_simplified) = the_vector(row_index);
			counter_row_simplified++;
		}		
	}
	return simplified_vector;
    }


	// since this function is virtual it actually can be instantiated here
    void CoupledLinearNS_TT::v_DoSolve(void)
    {


        switch(m_equationType)
        {
            case eUnsteadyStokes:
            case eUnsteadyNavierStokes:
                //AdvanceInTime(m_steps);
                UnsteadySystem::v_DoSolve();
                break;
            case eSteadyStokes:
            case eSteadyOseen:
            case eSteadyLinearisedNS:
            {
                Solve();				
                break;
            }
            case eSteadyNavierStokes:
            {	
                Timer Generaltimer;
                Generaltimer.Start();
                
                int Check(0);
                
                //Saving the init datas
                Checkpoint_Output(Check);
                Check++;
                
                cout<<"We execute INITIALLY SolveSteadyNavierStokes for m_kinvis = "<<m_kinvis<<" (<=> Re = "<<1/m_kinvis<<")"<<endl;
                SolveSteadyNavierStokes();
                
                while(m_kinvis > m_kinvisMin)
                {		
                    if (Check == 1)
                    {
                        cout<<"We execute SolveSteadyNavierStokes for m_kinvis = "<<m_kinvis<<" (<=> Re = "<<1/m_kinvis<<")"<<endl;
                        SolveSteadyNavierStokes();
                        Checkpoint_Output(Check);
                        Check++;
                    }
                    
                    Continuation();
                    
                    if (m_kinvis > m_kinvisMin)
                    {
                        cout<<"We execute SolveSteadyNavierStokes for m_kinvis = "<<m_kinvis<<" (<=> Re = "<<1/m_kinvis<<")"<<endl;
                        SolveSteadyNavierStokes();
                        Checkpoint_Output(Check);
                        Check++;
                    }
                }
                
                
                Generaltimer.Stop();
                cout<<"\nThe total calculation time is : " << Generaltimer.TimePerTest(1)/60 << " minute(s). \n\n";
                
                break;
            }
            case eNoEquationType:
            default:
                ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS");
        }

    }

    void CoupledLinearNS_TT::DoInitialiseAdv(Array<OneD, NekDouble> myAdvField_x, Array<OneD, NekDouble> myAdvField_y)
    {
	// only covers case eSteadyOseen

	Array<OneD, Array<OneD, NekDouble> > myAdvField(2);
	myAdvField[0] = myAdvField_x;
	myAdvField[1] = myAdvField_y;

        std::vector<std::string> fieldStr;
        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
             fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
        }
//	cout << "fieldStr[0] " << fieldStr[0] << endl;
//        EvaluateFunction(fieldStr,AdvField,"AdvectionVelocity"); // defined in EquationSystem

        SetUpCoupledMatrix(0.0, myAdvField, false);

    }

    void CoupledLinearNS_TT::setDBC(Eigen::MatrixXd collect_f_all)
    {
	no_dbc_in_loc = 0;
	no_not_dbc_in_loc = 0;
	for ( int index_c_f_bnd = 0; index_c_f_bnd < curr_f_bnd.size(); index_c_f_bnd++ )
	{
		if (collect_f_all(index_c_f_bnd,0) == collect_f_all(index_c_f_bnd,1))
		{
			no_dbc_in_loc++;
			elem_loc_dbc.insert(index_c_f_bnd);
		}
		else
		{
			no_not_dbc_in_loc++;
			elem_not_loc_dbc.insert(index_c_f_bnd);
		}
	}
    }

    void CoupledLinearNS_TT::setDBC_M(Eigen::MatrixXd collect_f_all)
    {
	// compute the dofs after Mtrafo multiplication to get global indices instead of local indices
	// Mtrafo = Eigen::MatrixXd (RB_A.rows(), nBndDofs);
	M_no_dbc_in_loc = 0;
	M_no_not_dbc_in_loc = 0;
	Eigen::VectorXd compare_vec1 = Mtrafo.transpose() * collect_f_all.block( 0, 0, curr_f_bnd.size(), 1);
	Eigen::VectorXd compare_vec2 = Mtrafo.transpose() * collect_f_all.block( 0, 1, curr_f_bnd.size(), 1);
	for ( int index_c_f_bnd = 0; index_c_f_bnd < compare_vec1.rows(); index_c_f_bnd++ )
	{
		if (compare_vec1(index_c_f_bnd) == compare_vec2(index_c_f_bnd))
		{
			M_no_dbc_in_loc++; // is actually no_dbc_in_global
			M_elem_loc_dbc.insert(index_c_f_bnd);
		}
		else
		{
			M_no_not_dbc_in_loc++;
			M_elem_not_loc_dbc.insert(index_c_f_bnd);
		}
	}
	M_truth_size = compare_vec1.rows() + curr_f_p.size() + curr_f_int.size();  // compare_vec1.rows() corresponds to nBndDofs
	M_truth_size_without_DBC = M_no_not_dbc_in_loc + curr_f_p.size() + curr_f_int.size();
	M_f_bnd_dbc = Eigen::VectorXd::Zero(M_no_dbc_in_loc);
	M_f_bnd_dbc_full_size = Eigen::VectorXd::Zero(M_truth_size);
	M_RB = Eigen::MatrixXd::Zero(M_truth_size_without_DBC, PODmodes.cols());
	// first multiply the bnd part in PODmodes with Mtrafo
	Eigen::MatrixXd M_PODmodes_bnd = Mtrafo.transpose() * PODmodes.block( 0, 0, curr_f_bnd.size(), RBsize );
	Eigen::MatrixXd M_collect_f_all_bnd = Mtrafo.transpose() * collect_f_all.block( 0, 0, curr_f_bnd.size(), Nmax );
	M_collect_f_all = Eigen::MatrixXd::Zero( M_truth_size , Nmax );
	int counter_all = 0;
	int counter_dbc = 0;
	for (int index=0; index < M_truth_size; ++index)  // take from the M_PODmodes_bnd if index is below compare_vec1.rows(), otherwise from PODmodes
	{
		if (!M_elem_loc_dbc.count(index))
		{
			if (index < compare_vec1.rows())
			{
				M_RB.row(counter_all) = M_PODmodes_bnd.row(index);
			}
			else
			{
				M_RB.row(counter_all) = PODmodes.row(index);
			}
			M_f_bnd_dbc_full_size(index) = 0;
			counter_all++;
		}
		else
		{
			M_f_bnd_dbc_full_size(index) = collect_f_all(index,0);
			M_f_bnd_dbc(counter_dbc) = collect_f_all(index,0);
			counter_dbc++;
		}
		if (index < compare_vec1.rows())
		{
			M_collect_f_all.row(index) = M_collect_f_all_bnd.row(index);
		}
		else
		{
			M_collect_f_all.row(index) = collect_f_all.row(index);
		}

	}
	elem_loc_dbc = M_elem_loc_dbc;
	f_bnd_dbc_full_size = M_f_bnd_dbc_full_size;
	RB = M_RB; // could be discussed if desired like that...
	// does not diminish size  collect_f_all = M_collect_f_all; // could be discussed if desired like that...
    }

    void CoupledLinearNS_TT::v_DoInitialise(void)
    {
        switch(m_equationType)
        {
            case eUnsteadyStokes:
            case eUnsteadyNavierStokes:
            {
                
//                LibUtilities::TimeIntegrationMethod intMethod;
//                std::string TimeIntStr = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");
//                int i;
//                for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
//                {
//                    if(boost::iequals(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr))
//                    {
//                        intMethod = (LibUtilities::TimeIntegrationMethod)i;
//                        break;
//                    }
//                }
//
//                ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
//
//                m_integrationScheme = LibUtilities::GetTimeIntegrationWrapperFactory().CreateInstance(LibUtilities::TimeIntegrationMethodMap[intMethod]);
                
                // Could defind this from IncNavierStokes class? 
                m_ode.DefineOdeRhs(&CoupledLinearNS::EvaluateAdvection, this);
                
                m_ode.DefineImplicitSolve(&CoupledLinearNS::SolveUnsteadyStokesSystem,this);
                
                // Set initial condition using time t=0
                
                SetInitialConditions(0.0);
                
            }
        case eSteadyStokes:
            SetUpCoupledMatrix(0.0);
            break;
        case eSteadyOseen:
            {

		Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.num_elements());
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    AdvField[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }
                
                ASSERTL0(m_session->DefinesFunction("AdvectionVelocity"),
                         "Advection Velocity section must be defined in "
                         "session file.");
                
                std::vector<std::string> fieldStr;
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                }

//		cout << "fieldStr[0] " << fieldStr[0] << endl;
//		cout << "fieldStr.size() " << fieldStr.size() << endl;

                EvaluateFunction(fieldStr,AdvField,"AdvectionVelocity"); // defined in EquationSystem

//		cout << "AdvField.num_elements() " << AdvField.num_elements() << endl;
//		cout << "AdvField[0].num_elements() " << AdvField[0].num_elements() << endl;                                

//		DefineRBspace();

                SetUpCoupledMatrix(0.0, AdvField, false);

            }
            break;
        case eSteadyNavierStokes:
            {				
                m_session->LoadParameter("KinvisMin", m_kinvisMin);
                m_session->LoadParameter("KinvisPercentage", m_KinvisPercentage);
                m_session->LoadParameter("Tolerence", m_tol);
                m_session->LoadParameter("MaxIteration", m_maxIt);
                m_session->LoadParameter("MatrixSetUpStep", m_MatrixSetUpStep);
                m_session->LoadParameter("Restart", m_Restart);
                
                
                DefineForcingTerm();
                
                if (m_Restart == 1)
                {
                    ASSERTL0(m_session->DefinesFunction("Restart"),
                             "Restart section must be defined in session file.");
                    
                    Array<OneD, Array<OneD, NekDouble> > Restart(m_velocity.num_elements());
                    for(int i = 0; i < m_velocity.num_elements(); ++i)
                    {
                        Restart[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                    }
                    std::vector<std::string> fieldStr;
                    for(int i = 0; i < m_velocity.num_elements(); ++i)
                    {
                        fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                    }
                    EvaluateFunction(fieldStr, Restart, "Restart");
                    
                    for(int i = 0; i < m_velocity.num_elements(); ++i)
                    {
                        m_fields[m_velocity[i]]->FwdTrans_IterPerExp(Restart[i], m_fields[m_velocity[i]]->UpdateCoeffs());
                    }
                    cout << "Saving the RESTART file for m_kinvis = "<< m_kinvis << " (<=> Re = " << 1/m_kinvis << ")" <<endl;
                }
                else //We solve the Stokes Problem
                {
                    
                    /*Array<OneD, Array<OneD, NekDouble> >ZERO(m_velocity.num_elements());
                     *					
                     *					for(int i = 0; i < m_velocity.num_elements(); ++i)
                     *					{				
                     *						ZERO[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                     *						m_fields[m_velocity[i]]->FwdTrans(ZERO[i], m_fields[m_velocity[i]]->UpdateCoeffs());
                     }*/
                    
                    SetUpCoupledMatrix(0.0);						
                    m_initialStep = true;
                    m_counter=1;
                    //SolveLinearNS(m_ForcingTerm_Coeffs);
                    Solve();
                    m_initialStep = false;
                    cout << "Saving the Stokes Flow for m_kinvis = "<< m_kinvis << " (<=> Re = " << 1/m_kinvis << ")" <<endl;
                }
            }
            break;
        case eSteadyLinearisedNS:
            {                
                SetInitialConditions(0.0);
                
                Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.num_elements());
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    AdvField[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0);
                }
                
                ASSERTL0(m_session->DefinesFunction("AdvectionVelocity"),
                         "Advection Velocity section must be defined in "
                         "session file.");
                
                std::vector<std::string> fieldStr;
                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    fieldStr.push_back(m_boundaryConditions->GetVariable(m_velocity[i]));
                }
                EvaluateFunction(fieldStr,AdvField,"AdvectionVelocity");
                
                SetUpCoupledMatrix(m_lambda,AdvField,true);
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for CoupledLinearNS");
        }

    }

    Eigen::VectorXd CoupledLinearNS_TT::trafoSnapshot(Eigen::VectorXd RB_via_POD, double kInvis)
    {

	


	return RB_via_POD;
    }

    void CoupledLinearNS_TT::set_MtM()
    {
	nBndDofs = m_locToGloMap[0]->GetNumGlobalBndCoeffs();  // number of global bnd dofs
        const Array<OneD,const int>& loctoglobndmap = m_locToGloMap[0]->GetLocalToGlobalBndMap();
        const Array<OneD,const NekDouble>& loctoglobndsign = m_locToGloMap[0]->GetLocalToGlobalBndSign();
//	Eigen::MatrixXd Mtrafo(RB_A.rows(), nBndDofs);
	M_truth_size = curr_f_bnd.size() + curr_f_p.size() + curr_f_int.size();  // compare_vec1.rows() corresponds to nBndDofs
	if (debug_mode)
	{
		cout << "Local dof size, also M_truth_size is " << curr_f_bnd.size() + curr_f_p.size() + curr_f_int.size() << endl;
	}
	M_truth_size_without_DBC = no_not_dbc_in_loc + curr_f_p.size() + curr_f_int.size();
	Mtrafo = Eigen::MatrixXd (f_bnd_size, nBndDofs);
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
	int nsize_bndry_p1 = loctoglobndmap.num_elements() / nel;
	int nsize_bndry = nsize_bndry_p1-1;
	for (int curr_elem = 0; curr_elem < nel; curr_elem++)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		int cnt_no_pp = curr_elem*nsize_bndry;
		for ( int index_ele = 0; index_ele < nsize_bndry_p1; index_ele++ )
		{
			int gid1 = loctoglobndmap[cnt+index_ele];
			int sign1 = loctoglobndsign[cnt+index_ele];
			if ((gid1 >= 0) && (index_ele < nsize_bndry))
			{
				Mtrafo(cnt_no_pp + index_ele, gid1) = sign1;
			}
		}
	}
	MtM = Mtrafo * Mtrafo.transpose();
	f_bnd_dbc = Eigen::VectorXd::Zero(no_dbc_in_loc);
	f_bnd_dbc_full_size = Eigen::VectorXd::Zero(PODmodes.rows());
	RB = Eigen::MatrixXd::Zero(PODmodes.rows() - no_dbc_in_loc, PODmodes.cols());
	int counter_all = 0;
	int counter_dbc = 0;
	for (int index=0; index < PODmodes.rows(); ++index)
	{
		if (!elem_loc_dbc.count(index))
		{
			RB.row(counter_all) = PODmodes.row(index);
			f_bnd_dbc_full_size(index) = 0;
			counter_all++;
		}
		else
		{
			f_bnd_dbc_full_size(index) = collect_f_all(index,0);
			f_bnd_dbc(counter_dbc) = collect_f_all(index,0);
			counter_dbc++;
		}
	}

    }

    void CoupledLinearNS_TT::gen_phys_base_vecs()
    {
	int RBsize = RB.cols();
	PhysBaseVec_x = Array<OneD, Array<OneD, double> > (RBsize); 
	PhysBaseVec_y = Array<OneD, Array<OneD, double> > (RBsize);
	
	if (debug_mode)
	{
		cout << " number of local dofs per velocity direction " << GetNcoeffs() << endl;
		cout << " number of quadrature dofs per velocity direction " << GetNpoints() << endl;
	}

	for (int curr_trafo_iter=0; curr_trafo_iter < RBsize; curr_trafo_iter++)
	{
		Eigen::VectorXd f_bnd = PODmodes.block(0, curr_trafo_iter, curr_f_bnd.size(), 1);
		Eigen::VectorXd f_int = PODmodes.block(curr_f_bnd.size()+curr_f_p.size(), curr_trafo_iter, curr_f_int.size(), 1);
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
	        Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(GetNcoeffs());
		Array<OneD, double> field_1(GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
	        int cnt = 0;
		int cnt1 = 0;
		int nvel = 2;
	        int nz_loc = 1;
	        int  nplanecoeffs = fields[0]->GetNcoeffs();
		int  nel  = m_fields[0]->GetNumElmts();
	        for(int i = 0; i < nel; ++i) 
	        {
	            int eid  = i;
	            fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
	            fields[0]->GetExp(eid)->GetInteriorMap(imap);
	            int nbnd   = bmap.num_elements();
	            int nint   = imap.num_elements();
	            int offset = fields[0]->GetCoeff_Offset(eid);
	            
	            for(int j = 0; j < nvel; ++j)
	            {
	                for(int n = 0; n < nz_loc; ++n)
	                {
	                    for(int k = 0; k < nbnd; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
	                    }
	                    
	                    for(int k = 0; k < nint; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
	                    }
	                    cnt  += nbnd;
	                    cnt1 += nint;
	                }
	            }
	        }
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);
		PhysBaseVec_x[curr_trafo_iter] = curr_PhysBaseVec_x;
		PhysBaseVec_y[curr_trafo_iter] = curr_PhysBaseVec_y;
		
	}

	eigen_phys_basis_x = Eigen::MatrixXd::Zero(GetNpoints(), RBsize);
	eigen_phys_basis_y = Eigen::MatrixXd::Zero(GetNpoints(), RBsize);
	for (int index_phys_base=0; index_phys_base<GetNpoints(); index_phys_base++)
	{
		for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
		{
			eigen_phys_basis_x(index_phys_base,index_RBsize) = PhysBaseVec_x[index_RBsize][index_phys_base];
			eigen_phys_basis_y(index_phys_base,index_RBsize) = PhysBaseVec_y[index_RBsize][index_phys_base];
		}
	}

	Eigen::VectorXd curr_col = eigen_phys_basis_x.col(0);
	double norm_curr_col = curr_col.norm();
	eigen_phys_basis_x.col(0) = curr_col / norm_curr_col;
	curr_col = eigen_phys_basis_y.col(0);
	norm_curr_col = curr_col.norm();
	eigen_phys_basis_y.col(0) = curr_col / norm_curr_col;
	Eigen::VectorXd orthogonal_complement;
	
	for (int orth_iter=1; orth_iter<RBsize; orth_iter++)
	{
		curr_col = eigen_phys_basis_x.col(orth_iter);
		Eigen::MatrixXd leftmostCols = eigen_phys_basis_x.leftCols(orth_iter);
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_x.col(orth_iter) = orthogonal_complement / norm_curr_col;
		curr_col = eigen_phys_basis_y.col(orth_iter);
		leftmostCols = eigen_phys_basis_y.leftCols(orth_iter);
		orthogonal_complement = curr_col - leftmostCols * leftmostCols.transpose() * curr_col;
		norm_curr_col = orthogonal_complement.norm();
		eigen_phys_basis_y.col(orth_iter) = orthogonal_complement / norm_curr_col;
	}

	orth_PhysBaseVec_x = Array<OneD, Array<OneD, double> > (RBsize); 
	orth_PhysBaseVec_y = Array<OneD, Array<OneD, double> > (RBsize); 
	for (int index_RBsize=0; index_RBsize<RBsize; index_RBsize++)
	{
		Array<OneD, double> curr_iter_x(GetNpoints());
		Array<OneD, double> curr_iter_y(GetNpoints());
		for (int index_phys_base=0; index_phys_base<GetNpoints(); index_phys_base++)	
		{
			curr_iter_x[index_phys_base] = eigen_phys_basis_x(index_phys_base,index_RBsize);
			curr_iter_y[index_phys_base] = eigen_phys_basis_y(index_phys_base,index_RBsize);			
		}
		orth_PhysBaseVec_x[index_RBsize] = curr_iter_x;
		orth_PhysBaseVec_y[index_RBsize] = curr_iter_y;			
	}
    }

    void CoupledLinearNS_TT::gen_proj_adv_terms()
    {
	RBsize = RB.cols();
	adv_mats_proj_x = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_mats_proj_y = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_x = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_y = Array<OneD, Eigen::VectorXd > (RBsize);
//	adv_vec_proj_x_newton = Array<OneD, Eigen::VectorXd > (RBsize);
//	adv_vec_proj_y_newton = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_x_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_y_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);


	Array<OneD, double> PhysBase_zero(GetNpoints(), 0.0);
	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];

		InitObject();

//		DoInitialiseAdv(curr_PhysBaseVec_x, PhysBase_zero); // call with parameter in phys state
		// needs to be replaced with a more gen. term		Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(RB_A.rows() + RB_Dbnd.rows() + RB_C.cols(), RB_A.cols() + RB_Dbnd.rows() + RB_B.cols() );
		Eigen::MatrixXd adv_matrix;
//		adv_matrix = Get_advection_matrix();      // <-- replace here with function -------------------------------------------------------------------------------------
		adv_matrix = gen_adv_mats_proj_x(curr_PhysBaseVec_x, use_Newton);
		Eigen::VectorXd add_to_rhs_adv(M_truth_size); // probably need this for adv and non-adv
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;
		Eigen::MatrixXd adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);

		adv_vec_proj_x_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);


		if (use_Newton)
		{
//			Eigen::VectorXd add_to_rhs_adv_newton(M_truth_size); 
//			add_to_rhs_adv_newton = adv_matrix * PODmodes * PODmodes.transpose() * collect_f_all.col(3);      
//			Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton, elem_loc_dbc);
			// alt: not working
//			adv_rhs_add_newton = adv_matrix_simplified * remove_rows(collect_f_all.col(3), elem_loc_dbc);
			// end alt
//			Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
//			adv_vec_proj_x_newton[trafo_iter] = adv_rhs_proj_newton;

			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
//				add_to_rhs_adv_newton_RB = adv_matrix_simplified * RB.col(RB_counter);
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;

				adv_vec_proj_x_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}


		Eigen::VectorXd adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		Eigen::MatrixXd adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		Eigen::VectorXd adv_rhs_proj = RB.transpose() * adv_rhs_add;

		adv_mats_proj_x[trafo_iter] = adv_mat_proj;
		adv_vec_proj_x[trafo_iter] = adv_rhs_proj;

		adv_vec_proj_y_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);
//		DoInitialiseAdv(PhysBase_zero , curr_PhysBaseVec_y ); // call with parameter in phys state
//		adv_matrix = Get_advection_matrix();      // <-- replace here with function -------------------------------------------------------------------------------------
		adv_matrix = gen_adv_mats_proj_y(curr_PhysBaseVec_y, use_Newton);
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;   
		adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);
		adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		adv_rhs_proj = RB.transpose() * adv_rhs_add;
		adv_mats_proj_y[trafo_iter] = adv_mat_proj;
		adv_vec_proj_y[trafo_iter] = adv_rhs_proj;

		if (use_Newton)
		{
//			Eigen::VectorXd add_to_rhs_adv_newton(M_truth_size); 
//			add_to_rhs_adv_newton = adv_matrix  * PODmodes * PODmodes.transpose() *  collect_f_all.col(3);      
//			Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton, elem_loc_dbc);
//			Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
//			adv_vec_proj_y_newton[trafo_iter] = adv_rhs_proj_newton;


			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
//				add_to_rhs_adv_newton_RB = adv_matrix_simplified * RB.col(RB_counter);
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;

				adv_vec_proj_y_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}
	}
    }

    Eigen::MatrixXd CoupledLinearNS_TT::gen_adv_mats_proj_x(Array<OneD, double> curr_PhysBaseVec_x, int use_Newton)
    {
	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 
	Array<OneD, Eigen::MatrixXd > A_elem(m_fields[0]->GetNumElmts()); 
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int

        // Calculate derivative of base flow 
        Array<OneD, Array<OneD, NekDouble> > Advfield(m_velocity.num_elements());
        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
               Advfield[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0); // use here the input vector
	}
	Advfield[0] = curr_PhysBaseVec_x;

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
		int npoints = locExp->GetTotPoints();
		int eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(curr_elem);
		int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(eid);
		Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
	        if(use_Newton) // formerly isLinearNSEquation
	        {
                    int nv1;
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(int nv = 0; nv < nvel; ++nv)
                    {
                        for(nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }
                            
//                            if((nv1 == 2)&&(m_HomogeneousType == eHomogeneous1D))
  //                          {
    //                            Vmath::Zero(npoints,AdvDeriv[nv*nvel+nv1],1); // dU/dz = 0
      //                      }
        //                    else
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],Advfield[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                            }
                        }
                    }
	        }


                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = curr_PhysBaseVec_x[curr_elem*nphys + i];
		}
		Array<OneD, double> Ah_ele_vec = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_x_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        Ah_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nbmap)*Ahrows] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        C_ele_vec[i+(nv*nz_loc+n1)*nbmap + (j+(k*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[imap[j]];
                                    }
                                }
                            } 
			 }
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_x_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            int n1;
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u'.Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys, coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        B_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        D_ele_vec[j+(k*nz_loc+n1)*nimap + (i+(nv*nz_loc+n1)*nimap)*nsize_int] += coeffs[imap[j]];
                                    }
                                }
                            }
                        }

			} // for (int k = 0; k < 2; ++k)

		} // for (int i = 0; i < nimap; ++i)

		A_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_bndry );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				A_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		}
		D_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_int , nsize_int ); 
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	Eigen::MatrixXd D_adv_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	Eigen::MatrixXd B_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd A_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		A_adv_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = A_elem[i];
		B_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = B_elem[i];
		C_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = C_elem[i];
		D_adv_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = D_elem[i];
	}

	switch(globally_connected) {
		case 0:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 1:
			adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_adv_all * Mtrafo;
			adv_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_adv_all;
			adv_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_adv_all.transpose() * Mtrafo;
			adv_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 2:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
	}


	return adv_matrix;

	
    }

    Eigen::MatrixXd CoupledLinearNS_TT::gen_adv_mats_proj_y(Array<OneD, double> curr_PhysBaseVec_y, int use_Newton)
    {
	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 
	Array<OneD, Eigen::MatrixXd > A_elem(m_fields[0]->GetNumElmts()); 
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int

        // Calculate derivative of base flow 
        Array<OneD, Array<OneD, NekDouble> > Advfield(m_velocity.num_elements());
        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
               Advfield[i] = Array<OneD, NekDouble> (m_fields[m_velocity[i]]->GetTotPoints(),0.0); // use here the input vector
	}
	Advfield[1] = curr_PhysBaseVec_y;

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
		int npoints = locExp->GetTotPoints();
		int eid = m_fields[m_velocity[0]]->GetOffset_Elmt_Id(curr_elem);
		int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(eid);
		Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
	        if(use_Newton) // formerly isLinearNSEquation
	        {
                    int nv1;
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(int nv = 0; nv < nvel; ++nv)
                    {
                        for(nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }
                            
//                            if((nv1 == 2)&&(m_HomogeneousType == eHomogeneous1D))
  //                          {
    //                            Vmath::Zero(npoints,AdvDeriv[nv*nvel+nv1],1); // dU/dz = 0
      //                      }
        //                    else
                            {
                                locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],Advfield[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                            }
                        }
                    }
	        }


                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_y_part[i] = curr_PhysBaseVec_y[curr_elem*nphys + i];
		}
		Array<OneD, double> Ah_ele_vec = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        Ah_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nbmap)*Ahrows] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        C_ele_vec[i+(nv*nz_loc+n1)*nbmap + (j+(k*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[imap[j]];
                                    }
                                }
                            } 
			 }
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_y_coeffs[int(imap[j])];
						}
					}
				}
                        if(use_Newton) // formerly isLinearNSEquation
                        {
                            int n1;
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u'.Grad U terms 
                                Vmath::Vmul(npoints, phys, 1, AdvDeriv[k*nvel+nv], 1, tmpphys, 1);
                                locExp->IProductWRTBase(tmpphys, coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1)
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        B_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nimap)*nsize_bndry] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        D_ele_vec[j+(k*nz_loc+n1)*nimap + (i+(nv*nz_loc+n1)*nimap)*nsize_int] += coeffs[imap[j]];
                                    }
                                }
                            }
                        }

			} // for (int k = 0; k < 2; ++k)

		} // for (int i = 0; i < nimap; ++i)

		A_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_bndry );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				A_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_bndry , nsize_int );
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		}
		D_elem[curr_elem] = Eigen::MatrixXd::Zero( nsize_int , nsize_int ); 
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	Eigen::MatrixXd D_adv_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	Eigen::MatrixXd B_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd A_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		A_adv_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = A_elem[i];
		B_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = B_elem[i];
		C_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = C_elem[i];
		D_adv_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = D_elem[i];
	}

	switch(globally_connected) {
		case 0:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 1:
			adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_adv_all * Mtrafo;
			adv_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_adv_all;
			adv_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_adv_all.transpose() * Mtrafo;
			adv_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 2:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
	}


	return adv_matrix;


    }


    Array<OneD, Array<OneD, Eigen::MatrixXd > > CoupledLinearNS_TT::gen_adv_mats_proj_y_2d(Array<OneD, double> curr_PhysBaseVec_y, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_y_2d)
    {
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  curr_adv_mats_proj_y_2d(number_elem_trafo);
	adv_vec_proj_y_2d = Array<OneD, Array<OneD, Eigen::VectorXd > >(number_elem_trafo);
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		curr_adv_mats_proj_y_2d[i] = Array<OneD, Eigen::MatrixXd > (2);
		curr_adv_mats_proj_y_2d[i][0] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		curr_adv_mats_proj_y_2d[i][1] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		adv_vec_proj_y_2d[i] = Array<OneD, Eigen::VectorXd > (2);
		adv_vec_proj_y_2d[i][0] = Eigen::VectorXd::Zero(RBsize);
		adv_vec_proj_y_2d[i][1] = Eigen::VectorXd::Zero(RBsize);
	}
	// not what I wanted:   cout << "sizeof curr_adv_mats_proj_x_2d " << sizeof(curr_adv_mats_proj_x_2d) << endl;
	
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc;
        nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 

	Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
//	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
//	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
//	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_y_part[i] = curr_PhysBaseVec_y[curr_elem*nphys + i];
		}
		Array<OneD, Array<OneD, double> > Ah_ele_vec(2);
		Ah_ele_vec[0] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[1] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, Array<OneD, double> > B_ele_vec(2);
		B_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > C_ele_vec(2); 
		C_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > D_ele_vec(2);
		D_ele_vec[0] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		D_ele_vec[1] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[0][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_x_coeffs[int(bmap[j])];
							Ah_ele_vec[1][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[0][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(imap[j])];
							C_ele_vec[1][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(imap[j])];
						}
					}

				}
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[0][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(bmap[j])];
							B_ele_vec[1][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[0][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_x_coeffs[int(imap[j])];
							D_ele_vec[1][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_y_coeffs[int(imap[j])];
						}
					}
				}
			}

		} // for (int i = 0; i < nimap; ++i)

		Ah_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		Ah_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		Ah_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem][0](i,j) = Ah_ele_vec[0][ i + j*Ahrows ];
				Ah_elem[curr_elem][1](i,j) = Ah_ele_vec[1][ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		B_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		B_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem][0](i,j) = B_ele_vec[0][ i + j*nsize_bndry ];
				B_elem[curr_elem][1](i,j) = B_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		C_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		C_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem][0](i,j) = C_ele_vec[0][ i + j*nsize_bndry ];
				C_elem[curr_elem][1](i,j) = C_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		D_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		D_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem][0](i,j) = D_ele_vec[0][ i + j*nsize_int];
				D_elem[curr_elem][1](i,j) = D_ele_vec[1][ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	// would need a function which (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		curr_adv_mats_proj_y_2d[i][0] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 0, adv_vec_proj_y_2d[i][0]);
		curr_adv_mats_proj_y_2d[i][1] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 1, adv_vec_proj_y_2d[i][1]);
	}


	return curr_adv_mats_proj_y_2d;

    }

    Array<OneD, Array<OneD, Eigen::MatrixXd > > CoupledLinearNS_TT::gen_adv_mats_proj_x_2d(Array<OneD, double> curr_PhysBaseVec_x, Array<OneD, Array<OneD, Eigen::VectorXd > > &adv_vec_proj_x_2d)
    {
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  curr_adv_mats_proj_x_2d(number_elem_trafo);
	adv_vec_proj_x_2d = Array<OneD, Array<OneD, Eigen::VectorXd > >(number_elem_trafo);
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		curr_adv_mats_proj_x_2d[i] = Array<OneD, Eigen::MatrixXd > (2);
		curr_adv_mats_proj_x_2d[i][0] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		curr_adv_mats_proj_x_2d[i][1] = Eigen::MatrixXd::Zero(RBsize, RBsize);
		adv_vec_proj_x_2d[i] = Array<OneD, Eigen::VectorXd > (2);
		adv_vec_proj_x_2d[i][0] = Eigen::VectorXd::Zero(RBsize);
		adv_vec_proj_x_2d[i][1] = Eigen::VectorXd::Zero(RBsize);
	}
	// not what I wanted:   cout << "sizeof curr_adv_mats_proj_x_2d " << sizeof(curr_adv_mats_proj_x_2d) << endl;
	
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap;
        int nz_loc;
        nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1; 

	Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
//	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
//	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
//	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = curr_PhysBaseVec_x[curr_elem*nphys + i];
		}
		Array<OneD, Array<OneD, double> > Ah_ele_vec(2);
		Ah_ele_vec[0] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[1] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, Array<OneD, double> > B_ele_vec(2);
		B_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > C_ele_vec(2); 
		C_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > D_ele_vec(2);
		D_ele_vec[0] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		D_ele_vec[1] = Array<OneD, double> (nsize_int*nsize_int, 0.0);
		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[0][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_x_coeffs[int(bmap[j])];
							Ah_ele_vec[1][j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[0][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(imap[j])];
							C_ele_vec[1][ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(imap[j])];
						}
					}

				}
			} // for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)

		for (int i = 0; i < nimap; ++i)
		{

			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			
			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[0][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_x_coeffs[int(bmap[j])];
							B_ele_vec[1][ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += adv_y_coeffs[int(bmap[j])];
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[0][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_x_coeffs[int(imap[j])];
							D_ele_vec[1][ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += adv_y_coeffs[int(imap[j])];
						}
					}
				}
			}

		} // for (int i = 0; i < nimap; ++i)

		Ah_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		Ah_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		Ah_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem][0](i,j) = Ah_ele_vec[0][ i + j*Ahrows ];
				Ah_elem[curr_elem][1](i,j) = Ah_ele_vec[1][ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		B_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		B_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem][0](i,j) = B_ele_vec[0][ i + j*nsize_bndry ];
				B_elem[curr_elem][1](i,j) = B_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		C_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		C_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem][0](i,j) = C_ele_vec[0][ i + j*nsize_bndry ];
				C_elem[curr_elem][1](i,j) = C_ele_vec[1][ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (2);
		D_elem[curr_elem][0] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		D_elem[curr_elem][1] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem][0](i,j) = D_ele_vec[0][ i + j*nsize_int];
				D_elem[curr_elem][1](i,j) = D_ele_vec[1][ i + j*nsize_int];
			}
		} 

	} // for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); curr_elem++)

	// would need a function which (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects
	for (int i = 0; i < number_elem_trafo; ++i)
	{

		time_t timer_1;
		time_t timer_2;
		time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

		curr_adv_mats_proj_x_2d[i][0] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 0, adv_vec_proj_x_2d[i][0]);

		time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
		double seconds = difftime(timer_1, timer_2);
//		cout << "time for a single adv_geo_mat_projector in seconds " << seconds << endl;


		curr_adv_mats_proj_x_2d[i][1] = adv_geo_mat_projector(Ah_elem, B_elem, C_elem, D_elem, i, 1, adv_vec_proj_x_2d[i][1]);
	}


	return curr_adv_mats_proj_x_2d;

    }

    Eigen::MatrixXd CoupledLinearNS_TT::adv_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > Ah_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &adv_vec_proj)
    {

	// eats up all the compute time...

	// this function (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects
	Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
        int nz_loc;
        nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
	Eigen::MatrixXd D_adv_all = Eigen::MatrixXd::Zero( nsize_int*nel, nsize_int*nel );
	Eigen::MatrixXd B_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel, nsize_int*nel );
	Eigen::MatrixXd C_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel, nsize_int*nel );
	Eigen::MatrixXd A_adv_all = Eigen::MatrixXd::Zero( nsize_bndry*nel, nsize_bndry*nel );

	time_t timer_1;
	time_t timer_2;
	time_t timer_3;
	time_t timer_4;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
//		cout << " i " << i << endl;
		int curr_elem_pos = get_curr_elem_pos(i);
		if (curr_elem_pos == curr_elem_trafo)
		{
			Eigen::MatrixXd curr_A_elem = Ah_elem[i][deriv_index].block(0, 0, Ah_elem[i][deriv_index].rows() - 1, Ah_elem[i][deriv_index].cols() - 1);
			Eigen::MatrixXd curr_B_elem = B_elem[i][deriv_index];
			Eigen::MatrixXd curr_C_elem = C_elem[i][deriv_index];
			Eigen::MatrixXd curr_D_elem = D_elem[i][deriv_index];

			A_adv_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = curr_A_elem;
			B_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_B_elem;
			C_adv_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_C_elem;
			D_adv_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = curr_D_elem;
	
		}
	}

	time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
	double seconds = difftime(timer_2, timer_1);
//	cout << "time within adv_geo_mat_projector for block write in seconds " << seconds << endl;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

	switch(globally_connected) {
		case 0:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 1:
			adv_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_adv_all * Mtrafo;
			adv_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_adv_all;
			adv_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_adv_all.transpose() * Mtrafo;
			adv_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
		case 2:
			adv_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_adv_all;
			adv_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_adv_all;
			adv_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_adv_all.transpose();
			adv_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_adv_all;
			break;
	}

	time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_2, timer_1);
//	cout << "time within adv_geo_mat_projector for adv_matrix write in seconds " << seconds << endl;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

	// also need to create appropriate vector right-hand-side contribution
	Eigen::VectorXd add_to_rhs_adv(M_truth_size); // probably need this for adv and non-adv
	add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;

	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
	Eigen::VectorXd adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc); // or is this eating up the time?
	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for remove_rows(add_to_rhs_adv, elem_loc_dbc) in seconds " << seconds << endl;

	adv_vec_proj = RB.transpose() * adv_rhs_add;

	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
	Eigen::MatrixXd adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);
	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for remove_cols_and_rows in seconds " << seconds << endl;

	// it is limiting the Eigen MatVec operation -- could compare to the Nektar++ MatVec operation -- I rather think one has to employ the sparsity structure of A

	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
	Eigen::MatrixXd adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB; // is this eating up the time?
	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for RB.transpose() * adv_matrix_simplified * RB in seconds " << seconds << endl;

//	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
//	Eigen::MatrixXd nn = (adv_matrix_simplified * RB.col(2)); // is this eating up the time?
//	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
//	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for (adv_matrix_simplified * RB.col(2)) " << seconds << endl;

//	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
//	Eigen::MatrixXd nn2 = RB.transpose(); // is this eating up the time?
//	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
//	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for RB.transpose() in seconds " << seconds << endl;

//	time(&timer_3);  /* get current time; same as: timer = time(NULL)  */
//	Eigen::MatrixXd nn3 = nn2 * nn; // is this eating up the time?
//	time(&timer_4);  /* get current time; same as: timer = time(NULL)  */
//	seconds = difftime(timer_4, timer_3);
//	cout << "time within adv_geo_mat_projector just for nn2 * nn in seconds " << seconds << endl; 

	time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
	seconds = difftime(timer_2, timer_1);
//	cout << "time within adv_geo_mat_projector for some MatVecs in seconds " << seconds << endl;

//	cout << " finished " << endl;

	return adv_mat_proj;

    }

    void CoupledLinearNS_TT::gen_proj_adv_terms_2d()
    {
	RBsize = RB.cols();

	adv_mats_proj_x = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_mats_proj_y = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_x = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_y = Array<OneD, Eigen::VectorXd > (RBsize);
	adv_vec_proj_x_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);
	adv_vec_proj_y_newton_RB = Array<OneD, Eigen::MatrixXd > (RBsize);

	// to be superseded by:
	adv_mats_proj_x_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > (RBsize); // should be RBsize x number_elem_trafo x 2 x RBsize x RBsize
	adv_mats_proj_y_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > (RBsize);
	adv_vec_proj_x_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > (RBsize);
	adv_vec_proj_y_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::VectorXd > > > (RBsize);

	Array<OneD, double> PhysBase_zero(GetNpoints(), 0.0);
	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		time_t timer_1;
		time_t timer_2;
		time(&timer_1);  /* get current time; same as: timer = time(NULL)  */

		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];

		adv_mats_proj_x_2d[trafo_iter] = gen_adv_mats_proj_x_2d(curr_PhysBaseVec_x, adv_vec_proj_x_2d[trafo_iter]);

		time(&timer_2);  /* get current time; same as: timer = time(NULL)  */
		double seconds = difftime(timer_1, timer_2);
		if (debug_mode)
		{
			cout << "time for a single gen_adv_mats_proj_x_2d in seconds " << seconds << endl;
			cout << "have 2 times RBsize of that " << endl;
		}

		adv_mats_proj_y_2d[trafo_iter] = gen_adv_mats_proj_y_2d(curr_PhysBaseVec_y, adv_vec_proj_y_2d[trafo_iter]);
	
	}
	

//	Array<OneD, double> PhysBase_zero(GetNpoints(), 0.0);
/*
	for(int trafo_iter = 0; trafo_iter < RBsize; trafo_iter++)
	{
		Array<OneD, double> curr_PhysBaseVec_x = orth_PhysBaseVec_x[trafo_iter];
		Array<OneD, double> curr_PhysBaseVec_y = orth_PhysBaseVec_y[trafo_iter];

		InitObject();

		DoInitialiseAdv(curr_PhysBaseVec_x, PhysBase_zero); // call with parameter in phys state
		// needs to be replaced with a more gen. term		Eigen::MatrixXd adv_matrix = Eigen::MatrixXd::Zero(RB_A.rows() + RB_Dbnd.rows() + RB_C.cols(), RB_A.cols() + RB_Dbnd.rows() + RB_B.cols() );
		Eigen::MatrixXd adv_matrix;
		adv_matrix = Get_advection_matrix();
		Eigen::VectorXd add_to_rhs_adv(M_truth_size); // probably need this for adv and non-adv
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;
		Eigen::MatrixXd adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);

		adv_vec_proj_x_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);

		if (use_Newton)
		{
			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
				adv_vec_proj_x_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}


		Eigen::VectorXd adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		Eigen::MatrixXd adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		Eigen::VectorXd adv_rhs_proj = RB.transpose() * adv_rhs_add;

		adv_mats_proj_x[trafo_iter] = adv_mat_proj;
		adv_vec_proj_x[trafo_iter] = adv_rhs_proj;

		adv_vec_proj_y_newton_RB[trafo_iter] = Eigen::MatrixXd::Zero(RBsize,RBsize);
		DoInitialiseAdv(PhysBase_zero , curr_PhysBaseVec_y ); // call with parameter in phys state
		adv_matrix = Get_advection_matrix();
		add_to_rhs_adv = adv_matrix * f_bnd_dbc_full_size;   
		adv_matrix_simplified = remove_cols_and_rows(adv_matrix, elem_loc_dbc);
		adv_rhs_add = remove_rows(add_to_rhs_adv, elem_loc_dbc);
		adv_mat_proj = RB.transpose() * adv_matrix_simplified * RB;
		adv_rhs_proj = RB.transpose() * adv_rhs_add;
		adv_mats_proj_y[trafo_iter] = adv_mat_proj;
		adv_vec_proj_y[trafo_iter] = adv_rhs_proj;

		if (use_Newton)
		{
			for(int RB_counter = 0; RB_counter < RBsize; RB_counter++)
			{			
				Eigen::VectorXd add_to_rhs_adv_newton_RB(M_truth_size); 
				add_to_rhs_adv_newton_RB = adv_matrix * PODmodes.col(RB_counter);      
				Eigen::VectorXd adv_rhs_add_newton = remove_rows(add_to_rhs_adv_newton_RB, elem_loc_dbc);
				Eigen::VectorXd adv_rhs_proj_newton = RB.transpose() * adv_rhs_add_newton;
				adv_vec_proj_y_newton_RB[trafo_iter].col(RB_counter) = adv_rhs_proj_newton;
			}
		}
	} */
    }


    Eigen::MatrixXd CoupledLinearNS_TT::DoTrafo(Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection, Array<OneD, Array<OneD, NekDouble> > snapshot_y_collection, Array<OneD, NekDouble> param_vector)
    {
	int Nmax = param_vector.num_elements();
	Eigen::MatrixXd collect_f_bnd( curr_f_bnd.size() , Nmax );
	Eigen::MatrixXd collect_f_p( curr_f_p.size() , Nmax );
	Eigen::MatrixXd collect_f_int( curr_f_int.size() , Nmax );
	for (int i=0; i<Nmax; i++)
	{
		Set_m_kinvis( param_vector[i] );	
//		cout << "CLNS_trafo.Get_m_kinvis " << CLNS_trafo.Get_m_kinvis() << endl;
	//	CLNS_trafo.DoInitialise();
		DoInitialiseAdv(snapshot_x_collection[i], snapshot_y_collection[i]); // replaces .DoInitialise();
		DoSolve();

		// compare the accuracy
		Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = UpdateFields();
		m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys());
		m_fields_t[1]->BwdTrans(m_fields_t[1]->GetCoeffs(), m_fields_t[1]->UpdatePhys());
		Array<OneD, NekDouble> out_field_trafo_x(GetNpoints(), 0.0);
		Array<OneD, NekDouble> out_field_trafo_y(GetNpoints(), 0.0);

		Eigen::VectorXd csx0_trafo(GetNpoints());
		Eigen::VectorXd csy0_trafo(GetNpoints());
		Eigen::VectorXd csx0(GetNpoints());
		Eigen::VectorXd csy0(GetNpoints());

		CopyFromPhysField(0, out_field_trafo_x); 
		CopyFromPhysField(1, out_field_trafo_y);
		for( int index_conv = 0; index_conv < GetNpoints(); ++index_conv)
		{
			csx0_trafo(index_conv) = out_field_trafo_x[index_conv];
			csy0_trafo(index_conv) = out_field_trafo_y[index_conv];
			csx0(index_conv) = snapshot_x_collection[i][index_conv];
			csy0(index_conv) = snapshot_y_collection[i][index_conv];
		}

		if (debug_mode)
		{
			cout << "csx0.norm() " << csx0.norm() << endl;
			cout << "csx0_trafo.norm() " << csx0_trafo.norm() << endl;
			cout << "csy0.norm() " << csy0.norm() << endl;
			cout << "csy0_trafo.norm() " << csy0_trafo.norm() << endl;
		}
		
		Eigen::VectorXd trafo_f_bnd = curr_f_bnd;
		Eigen::VectorXd trafo_f_p = curr_f_p;
		Eigen::VectorXd trafo_f_int = curr_f_int;

		collect_f_bnd.col(i) = trafo_f_bnd;
		collect_f_p.col(i) = trafo_f_p;
		collect_f_int.col(i) = trafo_f_int;

	}
	Eigen::MatrixXd collect_f_all( curr_f_bnd.size()+curr_f_p.size()+curr_f_int.size() , Nmax );
	collect_f_all.block(0,0,collect_f_bnd.rows(),collect_f_bnd.cols()) = collect_f_bnd;
	collect_f_all.block(collect_f_bnd.rows(),0,collect_f_p.rows(),collect_f_p.cols()) = collect_f_p;
	collect_f_all.block(collect_f_bnd.rows()+collect_f_p.rows(),0,collect_f_int.rows(),collect_f_int.cols()) = collect_f_int;

	return collect_f_all;
    }

    double CoupledLinearNS_TT::Geo_T(double w, int elemT, int index)
    {
/*
def Geo_T(w, elemT, index): # index 0: det, index 1,2,3,4: mat_entries
	if elemT == 0:
		T = np.matrix([[1, 0], [0, 1/w]])
	if elemT == 1:
		T = np.matrix([[1, 0], [0, 2/(3-w)]])
	if elemT == 2:
		T = np.matrix([[1, 0], [-(1-w)/2, 1]])
	if elemT == 3:
		T = np.matrix([[1, 0], [-(w-1)/2, 1]])
	if elemT == 4:
		T = np.matrix([[1, 0], [0, 1]])			
	if index == 0:
		return 1/np.linalg.det(T)
	if index == 1:
		return T[0,0]
	if index == 2:
		return T[0,1]
	if index == 3:
		return T[1,0]
	if index == 4:
		return T[1,1]


*/
	Eigen::Matrix2d T;
	if (elemT == 0) 
	{
		T << 1, 0, 0, 1/w;
	}
	else if (elemT == 1) 
	{
		T << 1, 0, 0, 2/(3-w);
	}
	else if (elemT == 2) 
	{
		T << 1, 0, -(1-w)/2, 1;
	}
	else if (elemT == 3) 
	{
		T << 1, 0, -(w-1)/2, 1;
	}
	else if (elemT == 4) 
	{
		T << 1, 0, 0, 1;
	}

//	cout << T << endl;

	if (index == 0) 
	{
		return 1/(T(0,0)*T(1,1) - T(0,1)*T(1,0)); // 1/det
	}
	else if (index == 1) 
	{
		return T(0,0);
	}
	else if (index == 2) 
	{
		return T(0,1);
	}
	else if (index == 3) 
	{
		return T(1,0);
	}
	else if (index == 4) 
	{
		return T(1,1);
	}

	return 0;


    }


    Eigen::MatrixXd CoupledLinearNS_TT::project_onto_basis( Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y)
    {
	Eigen::VectorXd c_snapshot_x(GetNpoints());
	Eigen::VectorXd c_snapshot_y(GetNpoints());
	for (int i = 0; i < GetNpoints(); ++i)
	{
		c_snapshot_x(i) = snapshot_x[i];
		c_snapshot_y(i) = snapshot_y[i];
	}
	Eigen::VectorXd curr_x_proj = eigen_phys_basis_x.transpose() * c_snapshot_x;
	Eigen::VectorXd curr_y_proj = eigen_phys_basis_y.transpose() * c_snapshot_y;
	curr_xy_projected = Eigen::MatrixXd::Zero(curr_x_proj.rows(), 2);
	curr_xy_projected.col(0) = curr_x_proj;
	curr_xy_projected.col(1) = curr_y_proj;

	return curr_xy_projected;

    }

    Eigen::MatrixXd CoupledLinearNS_TT::reproject_from_basis( Eigen::MatrixXd curr_xy_proj )
    {
		Eigen::VectorXd reproj_curr_x = eigen_phys_basis_x * curr_xy_proj.col(0);
		Eigen::VectorXd reproj_curr_y = eigen_phys_basis_y * curr_xy_proj.col(1);
		Eigen::MatrixXd curr_xy_reprojected = Eigen::MatrixXd::Zero(eigen_phys_basis_x.rows(), 2);
		curr_xy_reprojected.col(0) = reproj_curr_x;
		curr_xy_reprojected.col(1) = reproj_curr_y;

		return curr_xy_reprojected;

    }

    int CoupledLinearNS_TT::get_curr_elem_pos(int curr_elem)
    {
	// find within Array<OneD, std::set<int> > elements_trafo; the current entry
	for(int i = 0; i < elements_trafo.num_elements(); ++i)
	{
		if (elements_trafo[i].count(curr_elem))
		{
			return i;
		}
	}
    }

    Array<OneD, Array<OneD, NekDouble> > CoupledLinearNS_TT::trafo_current_para(Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y, Array<OneD, NekDouble> parameter_of_interest, Eigen::VectorXd & ref_f_bnd, Eigen::VectorXd & ref_f_p, Eigen::VectorXd & ref_f_int)
    {

	double w = parameter_of_interest[0];	
	double mKinvis = parameter_of_interest[1];
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 

	// verify and transform to bnd / p / int the snapshot data

        int nz_loc;
        nz_loc = 1;

        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nbndry = nsize_bndry;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
	int nint = nsize_int;
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
//	cout << "nsize_int " << nsize_int << endl;
//	cout << "nsize_p " << nsize_p << endl;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;

	Array<OneD, Eigen::MatrixXd > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1
	
	// set the newton forcing in terms of the current iterate
	myAdvField_Newton = Array<OneD, Array<OneD, NekDouble> > (2);
	myAdvField_Newton[0] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	myAdvField_Newton[1] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
	for (int i = 0; i<m_fields[0]->GetTotPoints(); ++i)
	{
		myAdvField_Newton[0][i] = snapshot_x[i];
		myAdvField_Newton[1][i] = snapshot_y[i];
	}


	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{
		int curr_elem_pos = get_curr_elem_pos(curr_elem);
		double detT = Geo_T(w, curr_elem_pos, 0);
		double Ta = Geo_T(w, curr_elem_pos, 1);
		double Tb = Geo_T(w, curr_elem_pos, 2);
		double Tc = Geo_T(w, curr_elem_pos, 3);
		double Td = Geo_T(w, curr_elem_pos, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td;
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int npoints = locExp->GetTotPoints();
                int phys_offset = m_fields[m_velocity[0]]->GetPhys_Offset(curr_elem);
                Array<OneD, Array<OneD, NekDouble> > AdvDeriv(nvel*nvel);
//                Array<OneD, NekDouble> tmpphys = m_fields[0]->UpdatePhys();
                Array<OneD, NekDouble> tmpphys = Array<OneD, NekDouble>(npoints,0.0);
//		cout << "ncoeffs " << ncoeffs << endl;
//		cout << "nphys " << nphys << endl;
//		cout << "pqsize " << pqsize << endl;   // pqsize == nphys and ncoeffs == nphys / 2 when?

                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = snapshot_x[curr_elem*nphys + i];
			curr_snap_y_part[i] = snapshot_y[curr_elem*nphys + i];
		}

                // Calculate derivative of base flow 
                if(use_Newton)
                {
                    int cnt = 0;
                    AdvDeriv[0] = Array<OneD, NekDouble>(nvel*nvel*npoints);
                    for(int nv = 0; nv < nvel; ++nv)
                    {
                        for(int nv1 = 0; nv1 < nvel; ++nv1)
                        {
                            if(cnt < nvel*nvel-1)
                            {
                                AdvDeriv[cnt+1] = AdvDeriv[cnt] + npoints;
                                ++cnt;
                            }
                            locExp->PhysDeriv(MultiRegions::DirCartesianMap[nv1],myAdvField_Newton[nv] + phys_offset, AdvDeriv[nv*nvel+nv1]);
                        }
                    }
                }
                

		Array<OneD, double> Ah_ele_vec(Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec(nsize_int*nsize_int, 0.0);
		Array<OneD, double> Dbnd_ele_vec(nsize_p*nsize_bndry, 0.0);
		Array<OneD, double> Dint_ele_vec(nsize_p*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] );
				}
				for (int j = 0; j < nimap; ++j)
				{
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] );
				}
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
/*					for (int il = 0; il < nphys; ++il)
					{
						cout << "curr_snap_x_part[il] " << curr_snap_x_part[il] << endl;
						cout << "deriv_0[il] " << deriv_0[il] << endl;
						cout << "tmpphys_x[il] " << tmpphys_x[il] << endl;
					} */
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Ta * adv_x_coeffs[int(bmap[j])] + Tc * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Ta * adv_x_coeffs[int(imap[j])] + Tc * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = detT * (Ta * pcoeffs_x[il] + Tc * pcoeffs_y[il]);
					}
				}
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
/*					for (int il = 0; il < nphys; ++il)
					{
						cout << "curr_snap_x_part[il] " << curr_snap_x_part[il] << endl;
						cout << "deriv_0[il] " << deriv_0[il] << endl;
						cout << "tmpphys_x[il] " << tmpphys_x[il] << endl;
					} */
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Tb * adv_x_coeffs[int(bmap[j])] + Td * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Tb * adv_x_coeffs[int(imap[j])] + Td * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = detT * (Tb * pcoeffs_x[il] + Td * pcoeffs_y[il]);
					}


				}


                        if(use_Newton)
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints,phys,1, AdvDeriv[k*nvel+nv], 1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1) // n1 is only zero
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        Ah_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nbmap)*Ahrows] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        C_ele_vec[i+(nv*nz_loc+n1)*nbmap + (j+(k*nz_loc+n1)*nimap)*nbndry] += coeffs[imap[j]];
                                    }
                                }
                            } 
			}

			} //for (int k = 0; k < 2; ++k)


		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					C_ele_vec[ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] );
				}
				for (int j = 0; j < nimap; ++j)
				{
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] );
				}
				if (k == 0)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_x_part.begin(), curr_snap_x_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Ta * adv_x_coeffs[int(bmap[j])] + Tc * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += detT * (Ta * adv_x_coeffs[int(imap[j])] + Tc * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = detT * (Ta * pcoeffs_x[il] + Tc * pcoeffs_y[il]);
					}
				}
				if (k == 1)
				{
					Array<OneD, double> tmpphys_x(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_0.begin(), tmpphys_x.begin(),  std::multiplies<double>() ); 
					Array<OneD, double> tmpphys_y(nphys, 0.0);
					std::transform( curr_snap_y_part.begin(), curr_snap_y_part.end(), deriv_1.begin(), tmpphys_y.begin(),  std::multiplies<double>() ); 
					locExp->IProductWRTBase(tmpphys_x, adv_x_coeffs);
					locExp->IProductWRTBase(tmpphys_y, adv_y_coeffs);
					for (int nv = 0; nv < 2; ++nv)
					{
						for (int j = 0; j < nbmap; ++j)
						{
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Tb * adv_x_coeffs[int(bmap[j])] + Td * adv_y_coeffs[int(bmap[j])]);
						}
						for (int j = 0; j < nimap; ++j)
						{
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ]  += detT * (Tb * adv_x_coeffs[int(imap[j])] + Td * adv_y_coeffs[int(imap[j])]);
						}
					}
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = detT * (Tb * pcoeffs_x[il] + Td * pcoeffs_y[il]);
					}
				}

                        if(use_Newton)
                        {
                            for(int nv = 0; nv < nvel; ++nv)
                            {
                                // u' . Grad U terms 
                                Vmath::Vmul(npoints,phys,1, AdvDeriv[k*nvel+nv], 1,tmpphys,1);
                                locExp->IProductWRTBase(tmpphys,coeffs);
                                
                                for(int n1 = 0; n1 < nz_loc; ++n1) // n1 is only zero
                                {
                                    for(int j = 0; j < nbmap; ++j)
                                    {
                                        B_ele_vec[j+(k*nz_loc+n1)*nbmap + (i+(nv*nz_loc+n1)*nimap)*nbndry] += coeffs[bmap[j]];
                                    }
                                    
                                    for(int j = 0; j < nimap; ++j)
                                    {
                                        D_ele_vec[j+(k*nz_loc+n1)*nimap + (i+(nv*nz_loc+n1)*nimap)*nint] += coeffs[imap[j]];
                                    }
                                }
                            } 
			}



			} //for (int k = 0; k < 2; ++k)


		}

 		// chosen choice: redo the nektar++ approach
		// possible alternatives:
		// build instead the sing_* matrices directly
		// or copy for test purposes from setupcoupledmats??

		Ah_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 
		Dbnd_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_bndry);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				Dbnd_elem[curr_elem](i,j) = Dbnd_ele_vec[ i + j*nsize_p];
			}
		} 
		Dint_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_int);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				Dint_elem[curr_elem](i,j) = Dint_ele_vec[ i + j*nsize_p];
			}
		} 

		
		D_elem[curr_elem] = D_elem[curr_elem].inverse();
		B_elem[curr_elem] = B_elem[curr_elem] * D_elem[curr_elem];
		Ah_elem[curr_elem].block(0,0,nsize_bndry,nsize_bndry) = Ah_elem[curr_elem].block(0,0,nsize_bndry,nsize_bndry) - B_elem[curr_elem] * C_elem[curr_elem].transpose();
		Eigen::MatrixXd Elem_Cinv = D_elem[curr_elem];
		Eigen::MatrixXd Elem_BCinv = B_elem[curr_elem];
		Eigen::MatrixXd Elem_Btilde = C_elem[curr_elem];
		Eigen::MatrixXd Elem_DintCinvDTint = Dint_elem[curr_elem] * Elem_Cinv * Dint_elem[curr_elem].transpose();
		Eigen::MatrixXd Elem_BCinvDTint_m_DTbnd = Elem_BCinv * Dint_elem[curr_elem].transpose() - Dbnd_elem[curr_elem].transpose();
		Eigen::MatrixXd Elem_DintCinvBTtilde_m_Dbnd = Dint_elem[curr_elem] * Elem_Cinv * Elem_Btilde.transpose() - Dbnd_elem[curr_elem];

		Eigen::MatrixXd Bh_curr_ele = Eigen::MatrixXd::Zero( nsize_bndry_p1, nsize_p_m1 );
		Eigen::MatrixXd Ch_curr_ele = Eigen::MatrixXd::Zero( nsize_p_m1, nsize_bndry_p1 );
		Eigen::MatrixXd Dh_curr_ele = Eigen::MatrixXd::Zero( nsize_p_m1, nsize_p_m1 );
		
		for (int i = 0; i < nsize_p_m1; ++i)
		{
			for (int j = 0; j < nsize_p_m1; ++j)
			{
				Dh_curr_ele(i,j) = -Elem_DintCinvDTint(i+1,j+1);
			}
		}

		for (int i = 0; i < nsize_bndry; ++i)
		{
			Ah_elem[curr_elem](i,nsize_bndry) = Elem_BCinvDTint_m_DTbnd(i,0);
			Ah_elem[curr_elem](nsize_bndry,i) = Elem_DintCinvBTtilde_m_Dbnd(0,i);
		}		
		Ah_elem[curr_elem](nsize_bndry,nsize_bndry) = -Elem_DintCinvDTint(0,0);

		for (int j = 0; j < nsize_p_m1; ++j)
		{
			for (int i = 0; i < nsize_bndry; ++i)
			{
				Bh_curr_ele(i,j) = Elem_BCinvDTint_m_DTbnd(i,j+1);
				Ch_curr_ele(j,i) = Elem_DintCinvBTtilde_m_Dbnd(j+1,i);
			}
		}

		for (int j = 0; j < nsize_p_m1; ++j)
		{
			Bh_curr_ele(nsize_bndry, j) = - Elem_DintCinvDTint(0, j+1);
			Ch_curr_ele(j, nsize_bndry) = - Elem_DintCinvDTint(j+1, 0);
		}

		Dh_curr_ele = Dh_curr_ele.inverse();
		Bh_curr_ele = Bh_curr_ele * Dh_curr_ele;
		Ah_elem[curr_elem] = Ah_elem[curr_elem] - Bh_curr_ele * Ch_curr_ele;
		Bh_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_p_m1  );
		Bh_elem[curr_elem] = Bh_curr_ele;
		Ch_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p_m1, nsize_bndry_p1 );
		Ch_elem[curr_elem] = Ch_curr_ele;
		Dh_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p_m1, nsize_p_m1 );
		Dh_elem[curr_elem] = Dh_curr_ele;

		///////////////////////////
		// temporary debugging
//		cout << "Dh_elem[curr_elem](0,0) " << Dh_elem[curr_elem](0,0) << endl;
//		cout << "Dh_elem[curr_elem](1,0) " << Dh_elem[curr_elem](1,0) << endl;
//		cout << "Dh_elem[curr_elem](0,1) " << Dh_elem[curr_elem](0,1) << endl;
//		cout << "Dh_elem[curr_elem](1,1) " << Dh_elem[curr_elem](1,1) << endl;
		/////////////////////////
		


//		sing_A[]



//		cout << "finished curr_elem " << curr_elem << endl;


	} // loop over curr_elem

	// nBndDofs already defined
	// NumDirBCs from GetNumGlobalDirBndCoeffs ??

	int NumDirBCs = m_locToGloMap[0]->GetNumGlobalDirBndCoeffs();
	nBndDofs = m_locToGloMap[0]->GetNumGlobalBndCoeffs();  // number of global bnd dofs, also in the .h
	int nGlobHomBndDofs = nBndDofs - NumDirBCs;
	int rows = nGlobHomBndDofs;
	if (debug_mode)
	{
	//	cout << "rows " << rows << endl;
	//	cout << "nBndDofs " << nBndDofs << endl;
	//	cout << "NumDirBCs " << NumDirBCs << endl;
	}
	Eigen::MatrixXd my_Gmat = Eigen::MatrixXd::Zero(rows, rows);
	int num_elem = m_fields[0]->GetNumElmts();
	Eigen::MatrixXd M_trafo = Eigen::MatrixXd::Zero(num_elem*nsize_bndry + num_elem, nGlobHomBndDofs);
	Eigen::MatrixXd M_trafo_no_pp = Eigen::MatrixXd::Zero(num_elem*nsize_bndry, nGlobHomBndDofs);

        const Array<OneD,const int>& loctoglobndmap = m_locToGloMap[0]->GetLocalToGlobalBndMap();
        const Array<OneD,const NekDouble>& loctoglobndsign = m_locToGloMap[0]->GetLocalToGlobalBndSign();

	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		int cnt_no_pp = curr_elem*nsize_bndry;
		Eigen::MatrixXd loc_Ah = Ah_elem[curr_elem];
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			int gid1 = loctoglobndmap[cnt + i] - NumDirBCs;
			int sign1 = loctoglobndsign[cnt + i];
			if (gid1 >= 0)
			{
				M_trafo(cnt+i, gid1) = sign1;
				if (i < nsize_bndry)
				{
					M_trafo_no_pp(cnt_no_pp + i, gid1) = sign1;
				}
				for (int j = 0; j < nsize_bndry_p1; ++j)
				{
					int gid2 = loctoglobndmap[cnt + j] - NumDirBCs;
					int sign2 = loctoglobndsign[cnt + j];
					if (gid2 >= 0)
					{
						my_Gmat(gid1,gid2) += sign1*sign2*loc_Ah(i,j);
					}
				}
			}
		}
	}

	int nGlobBndDofs = nBndDofs;

        Array<OneD, NekDouble > f_bnd(num_elem*nsize_bndry);
        NekVector< NekDouble  > F_bnd(f_bnd.num_elements(), f_bnd, eWrapper);
        Array<OneD, NekDouble > f_int(num_elem*nsize_int);
        NekVector< NekDouble  > F_int(f_int.num_elements(),f_int, eWrapper);

	Array<OneD, Array<OneD, NekDouble> > forcing(2); // local dofs
	forcing[0] = Array<OneD, NekDouble>(m_fields[m_velocity[0]]->GetNcoeffs(),0.0);
	forcing[1] = Array<OneD, NekDouble>(m_fields[m_velocity[0]]->GetNcoeffs(),0.0);

	// set the newton forcing in terms of the current iterate
//	myAdvField_Newton = Array<OneD, Array<OneD, NekDouble> > (2);
//	myAdvField_Newton[0] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
//	myAdvField_Newton[1] = Array<OneD, NekDouble> (m_fields[0]->GetTotPoints(),0.0);
        Array<OneD, Array<OneD, NekDouble> > AdvField(m_velocity.num_elements());
	Array<OneD, Array<OneD, NekDouble> > Eval_Adv(m_velocity.num_elements());
	Array<OneD, Array<OneD, NekDouble> > AdvTerm(m_velocity.num_elements());
        for(int il = 0; il < m_velocity.num_elements(); ++il)
        {
                AdvField[il] = Array<OneD, NekDouble> (m_fields[m_velocity[il]]->GetTotPoints(),0.0);
		Eval_Adv[il] = Array<OneD, NekDouble> (m_fields[m_velocity[il]]->GetTotPoints(),0.0);
		AdvTerm[il] = Array<OneD, NekDouble> (m_fields[m_velocity[il]]->GetNcoeffs(),0.0);
        }
/*	for (int i = 0; i<m_fields[0]->GetTotPoints(); ++i)
	{
		myAdvField_Newton[0][i] = snapshot_x[i];
		myAdvField_Newton[1][i] = snapshot_y[i];
	}
*/
	EvaluateAdvectionTerms(myAdvField_Newton, Eval_Adv);
        for (unsigned int i = 0; i < m_velocity.num_elements(); ++i)
        {
            bool waveSpace = m_fields[m_velocity[i]]->GetWaveSpace();
            m_fields[i]->SetWaveSpace(true);
//            m_fields[i]->IProductWRTBase(forcing_phys[i], forcing[i]);
            m_fields[m_velocity[i]]->IProductWRTBase(Eval_Adv[i], AdvTerm[i]); //(w, (u.grad)u)
	    if (use_Newton)
	    {
		for (unsigned int il = 0; il < forcing[i].num_elements(); ++il)
	        {
			forcing[i][il] = forcing[i][il] - AdvTerm[i][il];
//			cout << Eval_Adv[i][il] << endl;
		}
	    }
            m_fields[i]->SetWaveSpace(waveSpace);
        }


        // Assemble f_bnd and f_int

        int cnt = 0;
	int cnt1 = 0;
        for(int i = 0; i < num_elem; ++i) // loop over elements
        {
            int eid = i;
            m_fields[m_velocity[0]]->GetExp(eid)->GetBoundaryMap(bmap);
            m_fields[m_velocity[0]]->GetExp(eid)->GetInteriorMap(imap);
            int nbnd   = bmap.num_elements();
            int nint   = imap.num_elements();
            int offset = m_fields[m_velocity[0]]->GetCoeff_Offset(eid);
            
            for(int j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                    for(int k = 0; k < nbnd; ++k)
                    {
                        f_bnd[cnt+k] = forcing[j][ offset+bmap[k]];
                    }
                    for(int k = 0; k < nint; ++k)
                    {
                        f_int[cnt1+k] = forcing[j][ offset+imap[k]];
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
            }
        }

        Array<OneD, NekDouble > f_p(num_elem*nsize_p);
        NekVector<  NekDouble > F_p(f_p.num_elements(),f_p,eWrapper);
        NekVector<  NekDouble > F_p_tmp(num_elem*nsize_int);
        

	Eigen::VectorXd f_bnd_rhs_eigen = Eigen::VectorXd::Zero(f_bnd.num_elements());
	for (int i = 0; i < f_bnd.num_elements(); ++i)
	{
		f_bnd_rhs_eigen(i) = f_bnd[i];
	}
	Eigen::VectorXd f_p_rhs_eigen = Eigen::VectorXd::Zero(f_p.num_elements());
	for (int i = 0; i < f_p.num_elements(); ++i)
	{
		f_p_rhs_eigen(i) = f_p[i]; // should be undefined at this stage
	}
	Eigen::VectorXd f_int_rhs_eigen = Eigen::VectorXd::Zero(f_int.num_elements()); // num_elem*nsize_int
	Eigen::VectorXd f_p_tmp_eigen = Eigen::VectorXd::Zero(f_int.num_elements()); // num_elem*nsize_int
	for (int i = 0; i < f_int.num_elements(); ++i)
	{
		f_int_rhs_eigen(i) = f_int[i];
	}



        // fbnd does not currently hold the pressure mean
//        F_bnd = F_bnd - (*m_mat[mode].m_BCinv)*F_int;
//        F_p_tmp = (*m_mat[mode].m_Cinv)*F_int;
//        F_p = (*m_mat[mode].m_D_int) * F_p_tmp;

//	 should do these operations in Eigen
	for (int i = 0; i < num_elem; ++i)
	{
		f_bnd_rhs_eigen.segment(i*nsize_bndry, nsize_bndry) = f_bnd_rhs_eigen.segment(i*nsize_bndry, nsize_bndry) - B_elem[i] * f_int_rhs_eigen.segment(i*nsize_int, nsize_int);
		f_p_tmp_eigen.segment(i*nsize_int, nsize_int) = D_elem[i] * f_int_rhs_eigen.segment(i*nsize_int, nsize_int);
		f_p_rhs_eigen.segment(i*nsize_p, nsize_p) = Dint_elem[i] * f_p_tmp_eigen.segment(i*nsize_int, nsize_int);

	}


	////////////////////
	// temporary debugging
/*	cout << "f_bnd_rhs_eigen.size() " << f_bnd_rhs_eigen.size() << endl;
	cout << "f_bnd_rhs_eigen.norm() " << f_bnd_rhs_eigen.norm() << endl;
	cout << "f_p_rhs_eigen.size() " << f_p_rhs_eigen.size() << endl;
	cout << "f_p_rhs_eigen.norm() " << f_p_rhs_eigen.norm() << endl;

*/
	///////////////////



        Array<OneD, NekDouble > bnd   (m_locToGloMap[0]->GetNumGlobalCoeffs(),0.0);
        Array<OneD, NekDouble > fh_bnd(m_locToGloMap[0]->GetNumGlobalCoeffs(),0.0);

        const Array<OneD,const int>& loctoglomap = m_locToGloMap[0]->GetLocalToGlobalMap();
        const Array<OneD,const NekDouble>& loctoglosign = m_locToGloMap[0]->GetLocalToGlobalSign();
        
	int offset = 0;
	cnt = 0; 
        for(int i = 0; i < num_elem; ++i)
        {
            int eid  = i;
            int nbnd = nz_loc*m_fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            
            for(int j = 0; j < nvel; ++j)
            {
                for(int k = 0; k < nbnd; ++k)
                {
                    fh_bnd[loctoglomap[offset+j*nbnd+k]] += loctoglosign[offset+j*nbnd+k]*f_bnd_rhs_eigen(cnt+k);
                }
                cnt += nbnd;
            }
            
            int nint    = m_pressure->GetExp(eid)->GetNcoeffs();
            offset += nvel*nbnd + nint*nz_loc; 
        }
        
        offset = cnt1 = 0; 
        for(int i = 0; i <  num_elem; ++i)
        {
            int eid  = i;
            int nbnd = nz_loc*m_fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            int nint = m_pressure->GetExp(eid)->GetNcoeffs(); 
            
            for(int n = 0; n < nz_loc; ++n)
            {
                for(int j = 0; j < nint; ++j)
                {
                    fh_bnd[loctoglomap[offset + nvel*nbnd + n*nint+j]] = f_p_rhs_eigen(cnt1+j);
                }
                cnt1   += nint;
            }
            offset += nvel*nbnd + nz_loc*nint; 
        }

        //  Set Weak BC into f_bnd and Dirichlet Dofs in bnd
        const Array<OneD,const int>& bndmap = m_locToGloMap[0]->GetBndCondCoeffsToGlobalCoeffsMap();
	// Forcing function with weak boundary conditions and Dirichlet conditions
        int bndcnt=0;
        for(int k = 0; k < nvel; ++k)
        {
	    const Array<OneD, SpatialDomains::BoundaryConditionShPtr> bndConds = m_fields[k]->GetBndConditions();
            Array<OneD, const MultiRegions::ExpListSharedPtr> bndCondExp;
            bndCondExp = m_fields[k]->GetBndCondExpansions();
            for(int i = 0; i < bndCondExp.num_elements(); ++i)
            {
                const Array<OneD, const NekDouble > bndCondCoeffs = bndCondExp[i]->GetCoeffs();
		int cnt = 0;
                for(int n = 0; n < nz_loc; ++n)
                {
		    if (debug_mode)
		    {
//			    cout << "bndConds[i]->GetBoundaryConditionType() " << bndConds[i]->GetBoundaryConditionType() << endl;
//			    cout << "SpatialDomains::eDirichlet " << SpatialDomains::eDirichlet << endl;
		    }
                    if(bndConds[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        for(int j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                        {
                            if (m_equationType == eSteadyNavierStokes && m_initialStep == false)
                            {
                                //This condition set all the Dirichlet BC at 0 after
                                //the initial step of the Newton method
                                bnd[bndmap[bndcnt++]] = 0;
//				cout << "setting bnd[bndmap[bndcnt++]] = 0; " << endl;
                            }
                            else
                            {
//				cout << "setting bnd[bndmap[bndcnt++]] = bndCondCoeffs[cnt++] with " << bndCondCoeffs[cnt] << endl;
                                bnd[bndmap[bndcnt++]] = bndCondCoeffs[cnt++];

                            }
                        }
                    }
                    else
                    {                    
                        for(int j = 0; j < (bndCondExp[i])->GetNcoeffs(); j++)
                        {
 //                           cout << "setting fh_bnd[bndmap[bndcnt++]] += bndCondCoeffs[cnt++] with " << bndCondCoeffs[cnt] << endl;
                            fh_bnd[bndmap[bndcnt++]] += bndCondCoeffs[cnt++];

                        }
                    }
                }
            }
        }

	////////////////////////////
	// temporary debugging
//	cout << "fh_bnd.num_elements() " << fh_bnd.num_elements() << endl;
	double temp_norm = 0;
	for (int i = 0; i < fh_bnd.num_elements(); i++)
		temp_norm += fh_bnd[i];
//	temp_norm = sqrt(temp_norm);
//	cout << "fh_bnd norm " << temp_norm << endl;
	///////////////////////////

	////////////////////////////
	// temporary debugging
//	cout << "bnd.num_elements() " << bnd.num_elements() << endl;
	temp_norm = 0;
	for (int i = 0; i < bnd.num_elements(); i++)
	{
		temp_norm += bnd[i];
//		cout << bnd[i] << " ";
	}
//	temp_norm = sqrt(temp_norm);
//	cout << "bnd norm " << temp_norm << endl;
	///////////////////////////

	int nLocBndDofs = num_elem * nsize_bndry_p1;
	// ??	nGlobDofs = ??
	// nGlobBndDofs defined??, si, si, si come nBndDofs
	int nDirBndDofs = NumDirBCs;
	Eigen::VectorXd loc_dbc = Eigen::VectorXd::Zero(nLocBndDofs);
	Eigen::VectorXd loc_dbc_pt1 = Eigen::VectorXd::Zero(nLocBndDofs);
	Eigen::VectorXd loc_dbc_pt2 = Eigen::VectorXd::Zero(nLocBndDofs);
	Eigen::VectorXd V_GlobBnd = Eigen::VectorXd::Zero(nGlobBndDofs);

/*	cout << "bnd.num_elements() is the same as m_locToGloMap[0]->GetNumGlobalCoeffs() " << bnd.num_elements() << endl;
	cout << "nGlobBndDofs " << nGlobBndDofs << endl;
	cout << "m_locToGloMap[0]->AtLastLevel() " << m_locToGloMap[0]->AtLastLevel() << endl;	
	cout << "m_locToGloMap[0]->GetStaticCondLevel() " << m_locToGloMap[0]->GetStaticCondLevel() << endl;
	cout << "m_locToGloMap[0]->GetLowestStaticCondLevel() " << m_locToGloMap[0]->GetLowestStaticCondLevel() << endl;
*/
	for (int i = 0; i < nGlobBndDofs; ++i)
	{
		V_GlobBnd(i) = bnd[i];
	}
	
	for (int i = 0; i < nLocBndDofs; ++i)
	{
		loc_dbc(i) = loctoglobndsign[i] * V_GlobBnd(loctoglobndmap[i]);
	}

	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt = curr_elem*nsize_bndry_p1;
//		cout << "test " << endl;
		Eigen::MatrixXd loc_Ah = Ah_elem[curr_elem];
		Eigen::MatrixXd loc_Bh = Bh_elem[curr_elem];
//		cout << "test2 " << endl;
		Eigen::VectorXd lds = loc_dbc.segment(cnt,nsize_bndry_p1);
		Eigen::VectorXd loc_f_int_rhs_eigen = f_int_rhs_eigen.segment(curr_elem*nsize_int,nsize_int);
//		cout << "test3 " << endl;
//		cout << "loc_f_int_rhs_eigen " << loc_f_int_rhs_eigen << endl;
//		cout << "loc_f_int_rhs_eigen sizes " << loc_f_int_rhs_eigen.rows() << " " << loc_f_int_rhs_eigen.cols() << endl;
//		cout << "loc_Bh sizes " << loc_Bh.rows() << " " << loc_Bh.cols() << endl;	
//		loc_dbc.segment(cnt,nsize_bndry_p1) = loc_Bh * loc_f_int_rhs_eigen + loc_Ah * lds; // should not interfere with aliasing // something wrong with dims of loc_Bh * loc_f_int_rhs_eigen
		loc_dbc.segment(cnt,nsize_bndry_p1) = loc_Ah * lds;
/*		cout << "test4 " << endl;
		loc_dbc_pt1.segment(cnt,nsize_bndry_p1) = loc_Bh * loc_f_int_rhs_eigen;
		cout << "test5 " << endl;
		loc_dbc_pt2.segment(cnt,nsize_bndry_p1) = loc_Ah * lds;
		cout << "test6 " << endl; */
	}

//	cout << "loc_dbc_pt2.head(5) " << loc_dbc_pt2.head(5) << endl;

	// have here in Nektar:  V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;
	//  f_int(num_elem*nsize_int); is also same as f_int_rhs_eigen
//	cout << "f_int.num_elements() " << f_int.num_elements() << endl;
	// and using Bh_elem
/*	cout << "Bh_elem[0].cols() * num_elem " << Bh_elem[0].cols() * num_elem << endl;
	cout << "Bh_elem[0].rows() * num_elem " << Bh_elem[0].rows() * num_elem << endl;
	cout << "Ah_elem[0].cols() * num_elem " << Ah_elem[0].cols() * num_elem << endl;
	cout << "Ah_elem[0].rows() * num_elem " << Ah_elem[0].rows() * num_elem << endl;
*/
	Eigen::VectorXd V_GlobHomBndTmp = Eigen::VectorXd::Zero(nGlobHomBndDofs);
	Eigen::VectorXd tmp = Eigen::VectorXd::Zero(nGlobBndDofs);
	
	for (int i = 0; i < nLocBndDofs; ++i)
	{
		tmp(loctoglobndmap[i]) += loctoglobndsign[i] * loc_dbc(i);
	}

	offset = nDirBndDofs;
	V_GlobHomBndTmp = tmp.segment(offset, nGlobBndDofs-offset);
	// actually this, but is mixing Nektar and Eigen data structures V_GlobHomBndTmp = -V_GlobHomBndTmp + fh_bnd[NumDirBCs : nGlobHomBndDofs + NumDirBCs]
	Eigen::VectorXd fh_bnd_add1 = Eigen::VectorXd::Zero(nGlobHomBndDofs);
	for (int i = 0; i < nGlobHomBndDofs; ++i)
	{
		fh_bnd_add1(i) = fh_bnd[NumDirBCs + i] ;
	}
	V_GlobHomBndTmp = -V_GlobHomBndTmp + fh_bnd_add1;
	Eigen::VectorXd my_sys_in = V_GlobHomBndTmp;

	/////////////////////// actual solve here ////////////////////////////////
	Eigen::VectorXd my_Asolution = my_Gmat.colPivHouseholderQr().solve(my_sys_in);
	//////////////////////////////////////////////////////////////////////////
/*	cout << "my_Gmat.rows() " << my_Gmat.rows() << endl;
	cout << "my_Gmat.cols() " << my_Gmat.cols() << endl;
	cout << "my_sys_in.rows() " << my_sys_in.rows() << endl;
	cout << "my_sys_in.cols() " << my_sys_in.cols() << endl;
*/
	Eigen::VectorXd my_bnd_after = Eigen::VectorXd::Zero(m_locToGloMap[0]->GetNumGlobalCoeffs());
	for (int i = 0; i < m_locToGloMap[0]->GetNumGlobalCoeffs(); ++i)
	{
		my_bnd_after(i) = bnd[i];
	}
	my_bnd_after.segment(offset, nGlobBndDofs-offset) += my_Asolution;
	V_GlobBnd = my_bnd_after.segment(0, nGlobBndDofs);
	for (int i = 0; i < nLocBndDofs; ++i)
	{
		loc_dbc(i) = loctoglobndsign[i] * V_GlobBnd(loctoglobndmap[i]);
	}

	////////////////////////////
	// temporary debugging
/*	cout << "loc_dbc.num_elements() " << loc_dbc.num_elements() << endl;
	double temp_norm = 0;
	for (int i = 0; i < loc_dbc.num_elements(); i++)
		temp_norm += loc_dbc[i] * loc_dbc[i];
	temp_norm = sqrt(temp_norm);
	cout << "loc_dbc norm " << temp_norm << endl; */
	///////////////////////////

	Eigen::VectorXd Fint = Eigen::VectorXd::Zero(num_elem*nsize_p - num_elem);
	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt = curr_elem*nsize_bndry_p1;
		Eigen::MatrixXd loc_Ch = Ch_elem[curr_elem];
//		F_int[curr_elem*nsize_p_m1:curr_elem*nsize_p_m1+nsize_p_m1] = fh_bnd[nGlobHomBndDofs + NumDirBCs + curr_elem*nsize_p_m1: nGlobHomBndDofs + NumDirBCs + curr_elem*nsize_p_m1+nsize_p_m1] - np.dot(loc_Ch,  loc_dbc[cnt:cnt+nsize_bndry_p1])
		Eigen::VectorXd fh_bnd_add2 = Eigen::VectorXd::Zero(nsize_p_m1);
		for (int i = 0; i < nsize_p_m1; ++i)
		{
			fh_bnd_add2(i) = fh_bnd[NumDirBCs + nGlobHomBndDofs + curr_elem*nsize_p_m1 + i] ;
		}
		Fint.segment(curr_elem*nsize_p_m1,nsize_p_m1) = fh_bnd_add2 - loc_Ch * loc_dbc.segment(cnt, nsize_bndry_p1);  // actually an fh_bnd in case of body forcing as well
		Eigen::MatrixXd loc_Dh = Dh_elem[curr_elem];
		Fint.segment(curr_elem*nsize_p_m1,nsize_p_m1) = loc_Dh * Fint.segment(curr_elem*nsize_p_m1,nsize_p_m1).eval();
	}
	my_bnd_after.segment(nGlobBndDofs, m_locToGloMap[0]->GetNumGlobalCoeffs() - nGlobBndDofs) = Fint;

	for (int i = 0; i < m_locToGloMap[0]->GetNumGlobalCoeffs(); ++i)
	{
		 bnd[i] = my_bnd_after(i);
	}

	// can I now copy-paste from Nektar ?? -- try to do so:

//        const Array<OneD,const int>& loctoglomap = m_locToGloMap[0]->GetLocalToGlobalMap();
//        const Array<OneD,const NekDouble>& loctoglosign = m_locToGloMap[0]->GetLocalToGlobalSign();


	////////////////////////////
	// temporary debugging
//	cout << "bnd.num_elements() " << bnd.num_elements() << endl;
	temp_norm = 0;
	for (int i = 0; i < bnd.num_elements(); i++)
	{
		temp_norm += bnd[i] * bnd[i];
//		cout << bnd[i] << " ";
	}
	temp_norm = sqrt(temp_norm);
//	cout << "bnd norm " << temp_norm << endl;
	///////////////////////////


        // unpack pressure and velocity boundary systems. 
        offset = 0;
	cnt = 0; 
        int totpcoeffs = m_pressure->GetNcoeffs();
        Array<OneD, NekDouble> p_coeffs = m_pressure->UpdateCoeffs();
        for(int i = 0; i <  nel; ++i)
        {
            int eid  = i;
            int nbnd = nz_loc*m_fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            int nint = m_pressure->GetExp(eid)->GetNcoeffs(); 
            
            for(int j = 0; j < nvel; ++j)
            {
                for(int k = 0; k < nbnd; ++k)
                {
                    f_bnd[cnt+k] = loctoglosign[offset+j*nbnd+k]*bnd[loctoglomap[offset + j*nbnd + k]];
                }
                cnt += nbnd;
            }
            offset += nvel*nbnd + nint*nz_loc;
        }

        m_pressure->SetPhysState(false);
        
  //      Array<OneD, NekDouble > f_p(num_elem*nsize_p);
//        NekVector<  NekDouble > F_p(f_p.num_elements(),f_p,eWrapper);

        offset = 0; 
	cnt = 0; 
	cnt1 = 0;
        for(int i = 0; i < nel; ++i)
        {
            int eid  = i;
            int nint = m_pressure->GetExp(eid)->GetNcoeffs(); 
            int nbnd = m_fields[0]->GetExp(eid)->NumBndryCoeffs(); 
            cnt1 = m_pressure->GetCoeff_Offset(eid);
            
            for(int n = 0; n < nz_loc; ++n)
            {
                for(int j = 0; j < nint; ++j)
                {
                    p_coeffs[n*totpcoeffs + cnt1+j] = f_p[cnt+j] = bnd[loctoglomap[offset + (nvel*nz_loc)*nbnd + n*nint + j]];
                }
                cnt += nint;
            }
            offset += (nvel*nbnd + nint)*nz_loc;
        }

        // Back solve first level of static condensation for interior
        // velocity space and store in F_int
	// F_int should initially be empty without body forcing

//        F_int = F_int + Transpose(*m_mat[mode].m_D_int)*F_p - Transpose(*m_mat[mode].m_Btilde)*F_bnd;
//        F_int = (*m_mat[mode].m_Cinv)*F_int;
		// like to do the MatVec in Eigen and the transform into Nektar?? for f_p and bnd

	Eigen::VectorXd f_bnd_eigen = Eigen::VectorXd::Zero(f_bnd.num_elements());
	for (int i = 0; i < f_bnd.num_elements(); ++i)
	{
		f_bnd_eigen(i) = f_bnd[i];
	}
	Eigen::VectorXd f_p_eigen = Eigen::VectorXd::Zero(f_p.num_elements());
	for (int i = 0; i < f_p.num_elements(); ++i)
	{
		f_p_eigen(i) = f_p[i];
	}
	Eigen::VectorXd f_int_eigen = Eigen::VectorXd::Zero(f_int.num_elements()); // num_elem*nsize_int

	for (int curr_elem = 0; curr_elem < num_elem; ++curr_elem)
	{
		int cnt_Dint = curr_elem*nsize_p;
		int cnt_C = curr_elem*nsize_bndry;
		int cnt_D = curr_elem*nsize_int;
		Eigen::MatrixXd loc_Dint = Dint_elem[curr_elem];
		Eigen::MatrixXd loc_C = C_elem[curr_elem];
		Eigen::MatrixXd loc_D = D_elem[curr_elem];
		// f_int_rhs should come in as well...
//		cout << "f_p_eigen.rows, f_p_eigen.cols " << f_p_eigen.rows() << " " << f_p_eigen.cols() << endl;
//		cout << "f_bnd_eigen.rows, f_bnd_eigen.cols " << f_bnd_eigen.rows() << " " << f_bnd_eigen.cols() << endl;
//		Eigen::VectorXd nn =  f_p_eigen.segment(cnt_Dint, nsize_p);
//		cout << "nn.rows, nn.cols " << nn.rows() << " " << nn.cols() << endl;
		f_int_eigen.segment(curr_elem*nsize_int, nsize_int) = f_int_rhs_eigen.segment(curr_elem*nsize_int, nsize_int) + loc_Dint.transpose() * f_p_eigen.segment(cnt_Dint, nsize_p) - loc_C.transpose() * f_bnd_eigen.segment(cnt_C, nsize_bndry);  // also f_int_rhs here
		f_int_eigen.segment(curr_elem*nsize_int, nsize_int) = loc_D * f_int_eigen.segment(curr_elem*nsize_int, nsize_int);
	}

	for (int i = 0; i < f_int.num_elements(); ++i)
	{
		f_int[i] = f_int_eigen(i);
	}

	int nlc = m_locToGloMap[0]->GetNumLocalCoeffs();
	int nplanecoeffs = nlc;

        // Unpack solution from Bnd and F_int to v_coeffs 
        cnt = 0; 
	cnt1 = 0;
        for(int i = 0; i < nel; ++i) // loop over elements
        {
            int eid  = i;
            m_fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
            m_fields[0]->GetExp(eid)->GetInteriorMap(imap);
            int nbnd   = bmap.num_elements();
            int nint   = imap.num_elements();
            int offset = m_fields[0]->GetCoeff_Offset(eid);
            
            for(int j = 0; j < nvel; ++j) // loop over velocity fields 
            {
                for(int n = 0; n < nz_loc; ++n)
                {
                    for(int k = 0; k < nbnd; ++k)
                    {
                        m_fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd[cnt+k]);
                    }
                    
                    for(int k = 0; k < nint; ++k)
                    {
                        m_fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int[cnt1+k]);
                    }
                    cnt  += nbnd;
                    cnt1 += nint;
                }
            }
        }
        
        for(int j = 0; j < nvel; ++j) 
        {
            m_fields[j]->SetPhysState(false);
        }

	// the fields should be the same as the input:  Array<OneD, NekDouble> snapshot_x, Array<OneD, NekDouble> snapshot_y
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields_t = UpdateFields();
	m_fields_t[0]->BwdTrans(m_fields_t[0]->GetCoeffs(), m_fields_t[0]->UpdatePhys());
	m_fields_t[1]->BwdTrans(m_fields_t[1]->GetCoeffs(), m_fields_t[1]->UpdatePhys());
	Array<OneD, NekDouble> out_field_trafo_x(GetNpoints(), 0.0);
	Array<OneD, NekDouble> out_field_trafo_y(GetNpoints(), 0.0);

	Eigen::VectorXd csx0_trafo(GetNpoints());
	Eigen::VectorXd csy0_trafo(GetNpoints());
	Eigen::VectorXd csx0(GetNpoints());
	Eigen::VectorXd csy0(GetNpoints());

	CopyFromPhysField(0, out_field_trafo_x); 
	CopyFromPhysField(1, out_field_trafo_y);
	for( int index_conv = 0; index_conv < GetNpoints(); ++index_conv)
	{
		csx0_trafo(index_conv) = out_field_trafo_x[index_conv];
		csy0_trafo(index_conv) = out_field_trafo_y[index_conv];
		csx0(index_conv) = snapshot_x[index_conv];
		csy0(index_conv) = snapshot_y[index_conv];
	}

	if (debug_mode)
	{
		cout << "csx0.norm() " << csx0.norm() << endl;
		cout << "csx0_trafo.norm() " << csx0_trafo.norm() << endl;
		cout << "csy0.norm() " << csy0.norm() << endl;
		cout << "csy0_trafo.norm() " << csy0_trafo.norm() << endl;
	}

	Array<OneD, Array<OneD, NekDouble> > snapshot_result_phys_velocity_x_y(2);
	snapshot_result_phys_velocity_x_y[0] = Array<OneD, NekDouble>(GetNpoints());
	snapshot_result_phys_velocity_x_y[1] = Array<OneD, NekDouble>(GetNpoints());
	snapshot_result_phys_velocity_x_y[0] = out_field_trafo_x;
	snapshot_result_phys_velocity_x_y[1] = out_field_trafo_y;

	curr_f_bnd = Eigen::VectorXd::Zero(f_bnd.num_elements());
	for (int i_phys_dof = 0; i_phys_dof < f_bnd.num_elements(); i_phys_dof++)
	{
		curr_f_bnd(i_phys_dof) = f_bnd[i_phys_dof]; 
	}
	curr_f_p = Eigen::VectorXd::Zero(f_p.num_elements()); 
	for (int i_phys_dof = 0; i_phys_dof < f_p.num_elements(); i_phys_dof++)
	{
		curr_f_p(i_phys_dof) = f_p[i_phys_dof]; 
	}
	curr_f_int = Eigen::VectorXd::Zero(f_int.num_elements()); 
	for (int i_phys_dof = 0; i_phys_dof < f_int.num_elements(); i_phys_dof++)
	{
		curr_f_int(i_phys_dof) = f_int[i_phys_dof]; 
	}
	
	ref_f_bnd = curr_f_bnd;
	ref_f_p = curr_f_p;
	ref_f_int = curr_f_int;

	return snapshot_result_phys_velocity_x_y;



/*

		for i in range(0, nimap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(imap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			deriv_0 = np.dot(phys, (loc_cm0))
			deriv_1 = np.dot(phys, (loc_cm1))
			coeffs_0_0 = np.dot(deriv_0, loc_IP_d0)
			coeffs_0_1 = np.dot(deriv_0, loc_IP_d1)
			coeffs_1_0 = np.dot(deriv_1, loc_IP_d0)
			coeffs_1_1 = np.dot(deriv_1, loc_IP_d1)
			for k in range(0, 2):
				for j in range(0, nbmap):
					C_ele_vec[j+k*nbmap + (i+k*nimap)*nsize_bndry] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] )
				for j in range(0, nimap):
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] )
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					deriv_y = np.dot(phys, (loc_cm1))
					tmpphys = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					tmpphys_y = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv_y
					coeffs = np.dot(tmpphys, loc_IP)
					coeffs_y = np.dot(tmpphys_y, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Ta * coeffs[int(bmap[j])] + Tc * coeffs_y[int(bmap[j])])
						for j in range(0, nimap):
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += detT * (Ta * coeffs[int(imap[j])] + Tc * coeffs_y[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp)
					pcoeff_y = np.dot(deriv_y, loc_IPp)
					Dint_ele_vec[ (k*nimap + i)*nsize_p : ((k*nimap + i)*nsize_p + nsize_p) ] = detT * (Ta * pcoeff + Tc * pcoeff_y)
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							B_ele_vec[ j+nv*nbmap + (i+nv*nimap)*nsize_bndry ] += detT * (Tb * coeffs[int(bmap[j])] + Td * coeffs[int(bmap[j])])
						for j in range(0, nimap):
							D_ele_vec[ j+nv*nimap + (i+nv*nimap)*nsize_int ] += detT * (Tb * coeffs[int(imap[j])] + Td * coeffs[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp)
					Dint_ele_vec[ (k*nimap + i)*nsize_p : ((k*nimap + i)*nsize_p + nsize_p) ] = detT * (Tb * pcoeff + Td * pcoeff)

		for i in range(0, nbmap):
			coeffs = 0*np.arange(ncoeffs*1.0)
			coeffs[int(bmap[i])] = 1.0
			phys = np.dot(coeffs, loc_bwd_mat)
			deriv_0 = np.dot(phys, (loc_cm0))
			deriv_1 = np.dot(phys, (loc_cm1))
			coeffs_0_0 = np.dot(deriv_0, loc_IP_d0)
			coeffs_0_1 = np.dot(deriv_0, loc_IP_d1)
			coeffs_1_0 = np.dot(deriv_1, loc_IP_d0)
			coeffs_1_1 = np.dot(deriv_1, loc_IP_d1)
			for k in range(0, 2):
				for j in range(0, nbmap):
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += mKinvis * detT * ( c00*coeffs_0_0[int(bmap[j])] + c01*(coeffs_0_1[int(bmap[j])] + coeffs_1_0[int(bmap[j])]) + c11*coeffs_1_1[int(bmap[j])] )
				for j in range(0, nimap):
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += mKinvis * detT * ( c00*coeffs_0_0[int(imap[j])] + c01*(coeffs_0_1[int(imap[j])] + coeffs_1_0[int(imap[j])]) + c11*coeffs_1_1[int(imap[j])] )
			for k in range(0, 2):
				if (k == 0):
					deriv = np.dot(phys, (loc_cm0))
					deriv_y = np.dot(phys, (loc_cm1))
					tmpphys = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					tmpphys_y = u_x[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv_y
					coeffs = np.dot(tmpphys, loc_IP)
					coeffs_y = np.dot(tmpphys_y, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Ta * coeffs[int(bmap[j])] + Tc * coeffs_y[int(bmap[j])])
						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Ta * coeffs[int(imap[j])] + Tc * coeffs_y[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp)
					pcoeff_y = np.dot(deriv_y, loc_IPp)
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = detT * (Ta * pcoeff + Tc * pcoeff_y)
				if (k == 1):
					deriv = np.dot(phys, (loc_cm1))
					tmpphys = u_y[(curr_elem*nphys) : (curr_elem*nphys + nphys)] * deriv
					coeffs = np.dot(tmpphys, loc_IP)
					for nv in range(0, 2):
						for j in range(0, nbmap):
							Ah_ele_vec[ j+nv*nbmap + (i+nv*nbmap)*Ahrows ] += detT * (Tb * coeffs[int(bmap[j])] + Td * coeffs[int(bmap[j])])
						for j in range(0, nimap):
							C_ele_vec[ i+nv*nbmap + (j+nv*nimap)*nsize_bndry ] += detT * (Tb * coeffs[int(imap[j])] + Td * coeffs[int(imap[j])])
					pcoeff = np.dot(deriv, loc_IPp) 
					Dbnd_ele_vec[ (k*nbmap + i)*nsize_p : ((k*nbmap + i)*nsize_p + nsize_p) ] = detT * (Tb * pcoeff + Td * pcoeff)											


*/

/*	for curr_elem in range(0, num_elem):
		Ah_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		H1_bnd_ele_vec = 0*np.arange(Ahrows*Ahrows*1.0)
		B_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		C_ele_vec = 0*np.arange(nsize_bndry*nsize_int*1.0)
		D_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		DnoK_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		Dadv_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		H1_int_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		D1_ele_vec = 0*np.arange(nsize_int*nsize_int*1.0)
		Dbnd_ele_vec = 0*np.arange(nsize_bndry*nsize_p*1.0)
		Dint_ele_vec = 0*np.arange(nsize_int*nsize_p*1.0)
		loc_bwd_mat = bwdtrans[(curr_elem*ncoeffs) : (curr_elem*ncoeffs + ncoeffs) ,:]
		loc_cm0 = cartmap0[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_cm1 = cartmap1[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IP = IP[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IP_d0 = IP_d0[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IP_d1 = IP_d1[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		loc_IPp = IPp[(curr_elem*nphys) : (curr_elem*nphys + nphys) ,:]
		HelmMat = HelmMats[curr_elem,:]
		L00 = Lapl00[curr_elem,:]
		L11 = Lapl11[curr_elem,:]
		mimic_L00 = 0*Lapl00[curr_elem,:]
		curr_elem_pos = get_curr_elem_pos(curr_elem, geo_trafo_list)
		detT = Geo_T(w, curr_elem_pos, 0)
		Ta = Geo_T(w, curr_elem_pos, 1)
		Tb = Geo_T(w, curr_elem_pos, 2)
		Tc = Geo_T(w, curr_elem_pos, 3)
		Td = Geo_T(w, curr_elem_pos, 4)
		c00 = Ta*Ta + Tb*Tb
		c01 = Ta*Tc + Tb*Td
		c11 = Tc*Tc + Td*Td			
*/
    }


    void CoupledLinearNS_TT::do_geo_trafo()
    {
 
	// setting collect_f_all, making use of snapshot_x_collection, snapshot_y_collection

//	cout << "checking if all quantities set: " << endl;
//	cout << "Nmax " << Nmax << endl;
//	cout << "curr_f_bnd.size() " << curr_f_bnd.size() << endl;
//	cout << "curr_f_p.size() " << curr_f_p.size() << endl;
//	cout << "curr_f_int.size() " << curr_f_int.size() << endl;

	Eigen::MatrixXd collect_f_bnd( curr_f_bnd.size() , Nmax );
	Eigen::MatrixXd collect_f_p( curr_f_p.size() , Nmax );
	Eigen::MatrixXd collect_f_int( curr_f_int.size() , Nmax );

	Array<OneD, NekDouble> collected_qoi = Array<OneD, NekDouble> (Nmax);

        for(int i = 0; i < Nmax; ++i)
	{
//		cout << "general_param_vector[i][0] " << general_param_vector[i][0] << endl;
//		cout << "general_param_vector[i][1] " << general_param_vector[i][1] << endl;
//		cout << "curr_f_bnd.size() " << curr_f_bnd.size() << endl;
//		cout << "curr_f_p.size() " << curr_f_p.size() << endl;
//		cout << "curr_f_int.size() " << curr_f_int.size() << endl;
		// identify the geometry one - introduce in the .h some vars keeping the para_type
		// here now [0] is geometry 'w' and [1] is k_invis

		//for (int repeat_i = 0; repeat_i < 10; ++repeat_i)

		// setting collect_f_all, making use of snapshot_x_collection, snapshot_y_collection

		// take a timing for trafo_current_para

		Eigen::VectorXd ref_f_bnd;
		Eigen::VectorXd ref_f_p;
		Eigen::VectorXd ref_f_int;

		// if using Newton should set myAdvField

		Array<OneD, Array<OneD, NekDouble> > snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_x_collection[i], snapshot_y_collection[i], general_param_vector[i], ref_f_bnd, ref_f_p, ref_f_int); 

//		cout << "curr_f_bnd.norm() " << curr_f_bnd.norm() << endl;
//		cout << "ref_f_bnd.norm() " << ref_f_bnd.norm() << endl;

//		Eigen::VectorXd trafo_f_bnd = curr_f_bnd;
//		Eigen::VectorXd trafo_f_p = curr_f_p;
//		Eigen::VectorXd trafo_f_int = curr_f_int;

		if (do_trafo_check)
		{
			double L2error = 1;
			do
			{
				Array<OneD, Array<OneD, NekDouble> > prev_snapshot_result_phys_velocity_x_y = snapshot_result_phys_velocity_x_y;
				snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_result_phys_velocity_x_y[0], snapshot_result_phys_velocity_x_y[1], general_param_vector[i], ref_f_bnd, ref_f_p, ref_f_int); // setting collect_f_all, making use of snapshot_x_collection, snapshot_y_collection

//				cout << "curr_f_bnd.norm() " << curr_f_bnd.norm() << endl;
//				cout << "ref_f_bnd.norm() " << ref_f_bnd.norm() << endl;

//				trafo_f_bnd = curr_f_bnd;
//				trafo_f_p = curr_f_p;
//				trafo_f_int = curr_f_int;

				double L2error_x = L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]);
				double L2error_y = L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]);
				double L2error_x_ref = L2Error(0);
				double L2error_y_ref = L2Error(1);
				L2error = sqrt(L2error_x*L2error_x + L2error_y*L2error_y) / sqrt(L2error_x_ref*L2error_x_ref + L2error_y_ref*L2error_y_ref);
				cout << "relative L2error w.r.t. current iterate " << L2error << endl;
//				cout << " L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]) " << L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]) << endl;
//				cout << " L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]) " << L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]) << endl;
			}
			while ((L2error > 2e-5) && (!load_cO_snapshot_data_from_files));
		}

//		snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_result_phys_velocity_x_y[0], snapshot_result_phys_velocity_x_y[1], general_param_vector[i], ref_f_bnd, ref_f_p, ref_f_int);

//		cout << "curr_f_bnd.norm() " << curr_f_bnd.norm() << endl;
//		cout << "ref_f_bnd.norm() " << ref_f_bnd.norm() << endl;

		if (qoi_dof >= 0)
		{
			cout << "converged qoi dof " << snapshot_result_phys_velocity_x_y[1][qoi_dof] << endl;
			collected_qoi[i] = snapshot_result_phys_velocity_x_y[1][qoi_dof];
		}

		Eigen::VectorXd trafo_f_bnd = curr_f_bnd;
		Eigen::VectorXd trafo_f_p = curr_f_p;
		Eigen::VectorXd trafo_f_int = curr_f_int;

		collect_f_bnd.col(i) = trafo_f_bnd;
		collect_f_p.col(i) = trafo_f_p;
		collect_f_int.col(i) = trafo_f_int;

		// need to replace the snapshot data with the converged one for error computations
		// that means replace data in the snapshot_x_collection and snapshot_y_collection
		if (replace_snapshot_with_transformed)
		{
			snapshot_x_collection[i] = snapshot_result_phys_velocity_x_y[0];
			snapshot_y_collection[i] = snapshot_result_phys_velocity_x_y[1];
		}

//		cout << "collect_f_bnd.col(i).norm() " << collect_f_bnd.col(i).norm() << endl;

		// generate the correct string
		std::stringstream sstm;
		sstm << "Conv_Oseen_param" << i << ".fld";
		std::string filename = sstm.str();
		if (!load_cO_snapshot_data_from_files)
		{
			write_curr_field(filename);
		}

	}

	std::stringstream sstm;
	sstm << "FOM_qoi.txt";
	std::string LocROM_txt = sstm.str();
	const char* outname = LocROM_txt.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < Nmax; i0++)
		{
			myfile << collected_qoi[i0] << "\t";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 	


	collect_f_all = Eigen::MatrixXd::Zero( curr_f_bnd.size()+curr_f_p.size()+curr_f_int.size() , Nmax );
	collect_f_all.block(0,0,collect_f_bnd.rows(),collect_f_bnd.cols()) = collect_f_bnd;
	collect_f_all.block(collect_f_bnd.rows(),0,collect_f_p.rows(),collect_f_p.cols()) = collect_f_p;
	collect_f_all.block(collect_f_bnd.rows()+collect_f_p.rows(),0,collect_f_int.rows(),collect_f_int.cols()) = collect_f_int;

	// do the same for VV reference solutions
    for(int i = fine_grid_dir0*fine_grid_dir1; i < fine_grid_dir0*fine_grid_dir1; ++i)
	{

		cout << "\n attempting trafo for VV reference solutions \n \n";


		double err_threshold = 1e-9;
		int num_iter = 0;
		Eigen::VectorXd ref_f_bnd;
		Eigen::VectorXd ref_f_p;
		Eigen::VectorXd ref_f_int;
//		Array<OneD, Array<OneD, NekDouble> > snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_x_collection[0], snapshot_y_collection[0], fine_general_param_vector[i], ref_f_bnd, ref_f_p, ref_f_int);
		Array<OneD, Array<OneD, NekDouble> > snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_x_collection_VV[i], snapshot_y_collection_VV[i], fine_general_param_vector[i], ref_f_bnd, ref_f_p, ref_f_int);  

		if (1)
		{
			double L2error = 1;
			do
			{
				Array<OneD, Array<OneD, NekDouble> > prev_snapshot_result_phys_velocity_x_y = snapshot_result_phys_velocity_x_y;
				snapshot_result_phys_velocity_x_y = trafo_current_para(snapshot_result_phys_velocity_x_y[0], snapshot_result_phys_velocity_x_y[1], fine_general_param_vector[i], ref_f_bnd, ref_f_p, ref_f_int); // setting collect_f_all, making use of snapshot_x_collection, snapshot_y_collection


				double L2error_x = L2Error(0, prev_snapshot_result_phys_velocity_x_y[0]);
				double L2error_y = L2Error(1, prev_snapshot_result_phys_velocity_x_y[1]);
				double L2error_x_ref = L2Error(0);
				double L2error_y_ref = L2Error(1);
				L2error = sqrt(L2error_x*L2error_x + L2error_y*L2error_y) / sqrt(L2error_x_ref*L2error_x_ref + L2error_y_ref*L2error_y_ref);
				cout << "relative L2error w.r.t. current iterate " << L2error << endl;
				num_iter++;
				if (num_iter == 100)
				{
					num_iter = 0;
					err_threshold = err_threshold*10;
				}
			}
			while ((L2error > err_threshold));
		}

		std::stringstream sstm;
		sstm << "Conv_ref_Oseen_param" << i << ".fld";
		std::string filename = sstm.str();
		if (1)
		{
			write_curr_field(filename);
		}

	}

    }

    void CoupledLinearNS_TT::write_curr_field(std::string filename)
    {

        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.num_elements()+1);
        std::vector<std::string> variables(m_fields.num_elements()+1);
        int i;
        
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
	    if (debug_mode)
	    {
		    cout << "variables[i] " << variables[i] << endl;
	    }
        }

	if (debug_mode)
	{
		cout << "m_singleMode " << m_singleMode << endl;	
	}

	fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs(), 0.0);  
        variables[i] = "p"; 

	WriteFld(filename,m_fields[0],fieldcoeffs,variables);

    }


    void CoupledLinearNS_TT::compute_snapshots_geometry_params()
    {
	// generate matrices for affine form 
	// do full order solves using the affine form
	// can check that against simulations with appropriate .xml
	// cout << "starting compute geometry snapshots" << endl;
	// Nmax should be available as the total number of to be computed snapshots
	Array<OneD, NekDouble> zero_phys_init(GetNpoints(), 0.0);
	snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (Nmax);
	snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (Nmax);
        for(int i = 0; i < Nmax; ++i)
	{
		// what is the current parameter vector ?
		// cout << "general_param_vector[i][0] " << general_param_vector[i][0] << endl;
		// cout << "general_param_vector[i][1] " << general_param_vector[i][1] << endl;
		
		// assuming the loaded configuration is the reference config...
		// but then, how to solve this ?? -- do not have a solvable large-scale system ??
		cout << "Error: compute_snapshots_geometry_params() not implemented :) " << endl;
	}
    }

    void CoupledLinearNS_TT::run_local_ROM_online(std::set<int> current_cluster, int current_cluster_number)
    {
	// Question: how to init?
	// could use all-zero or the cluster-mean
	Array<OneD, NekDouble> cluster_mean_x(snapshot_x_collection[0].num_elements(), 0.0);
	Array<OneD, NekDouble> cluster_mean_y(snapshot_y_collection[0].num_elements(), 0.0);
	for (std::set<int>::iterator it=current_cluster.begin(); it!=current_cluster.end(); ++it)
	{
		for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
		{
			cluster_mean_x[i] += (1.0 / current_cluster.size()) * snapshot_x_collection[*it][i];
			cluster_mean_y[i] += (1.0 / current_cluster.size()) * snapshot_y_collection[*it][i];
		}
	}


	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(), 3);  // is of size M_truth_size
	// start sweeping 
	Eigen::MatrixXd collected_relative_L2errors_snaps = Eigen::VectorXd::Zero(Nmax);
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double w;
		if (parameter_space_dimension == 1)
		{
			current_nu = param_vector[current_index];
		}
		else if (parameter_space_dimension == 2)
		{
			Array<OneD, NekDouble> current_param = general_param_vector[current_index];
			w = current_param[0];	
			current_nu = current_param[1];
		}
		if (debug_mode)
		{
//			cout << " online phase current nu " << current_nu << endl;
//			cout << " online phase current w " << w << endl;
		}
		Set_m_kinvis( current_nu );
		if (use_Newton)
		{
//			DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
			DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
		}

		Eigen::MatrixXd curr_xy_proj = project_onto_basis(cluster_mean_x, cluster_mean_y);
//		Eigen::MatrixXd curr_xy_proj = project_onto_basis(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		Eigen::MatrixXd affine_mat_proj;
		Eigen::VectorXd affine_vec_proj;

		if (parameter_space_dimension == 1)
		{
			affine_mat_proj = gen_affine_mat_proj(current_nu);
			affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
//			cout << "aff mat 1d " << affine_mat_proj << endl;
//			cout << "aff vec 1d " << affine_vec_proj << endl;
		}
		else if (parameter_space_dimension == 2)
		{
			affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
//			cout << "aff mat 2d " << affine_mat_proj << endl;
			affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
//			cout << "aff vec 2d " << affine_vec_proj << endl;
//			Eigen::MatrixXd affine_mat_proj_1d = gen_affine_mat_proj(current_nu);
//			Eigen::VectorXd affine_vec_proj_1d = gen_affine_vec_proj(current_nu, current_index);
		}

		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
		double relative_change_error;
		int no_iter=0;
		// now start looping
		do
		{
			// for now only Oseen // otherwise need to do the DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
			Eigen::VectorXd prev_solve_affine = solve_affine;
			Eigen::VectorXd repro_solve_affine = RB * solve_affine;
			Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
			Array<OneD, double> field_x;
			Array<OneD, double> field_y;
			recover_snapshot_loop(reconstruct_solution, field_x, field_y);
			if (use_Newton)
			{
				DoInitialiseAdv(field_x, field_y);
			}
			curr_xy_proj = project_onto_basis(field_x, field_y);
			if (parameter_space_dimension == 1)
			{
				affine_mat_proj = gen_affine_mat_proj(current_nu);
				affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
			}
			else if (parameter_space_dimension == 2)
			{
				affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
				affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
			}
			solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
			relative_change_error = (solve_affine - prev_solve_affine).norm() / prev_solve_affine.norm();
//			cout << "relative_change_error " << relative_change_error << endl;
			no_iter++;
		} 
		while( ((relative_change_error > 1e-6) && (no_iter < 500)) );
		if (debug_mode)
			cout << " no_iterations " << no_iter << endl;
//		cout << "solve_affine " << solve_affine << endl;
		Eigen::VectorXd repro_solve_affine = RB * solve_affine;
		Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
		if (globally_connected == 1)
		{
			mat_compare.col(0) = M_collect_f_all.col(current_index);
		}
		else
		{
			mat_compare.col(0) = collect_f_all.col(current_index);
			if (debug_mode)
			{
				Eigen::VectorXd current_f_all = Eigen::VectorXd::Zero(collect_f_all.rows());
				current_f_all = collect_f_all.col(current_index);
				Eigen::VectorXd current_f_all_wo_dbc = remove_rows(current_f_all, elem_loc_dbc);
				Eigen::VectorXd proj_current_f_all_wo_dbc = RB.transpose() * current_f_all_wo_dbc;
//				cout << "proj_current_f_all_wo_dbc " << proj_current_f_all_wo_dbc << endl;
				Eigen::VectorXd correctRHS = affine_mat_proj * proj_current_f_all_wo_dbc;
				Eigen::VectorXd correction_RHS = correctRHS - affine_vec_proj;
//				cout << "correctRHS " << correctRHS << endl;
//				cout << "correction_RHS " << correction_RHS << endl;
			}
		}
		mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
		mat_compare.col(2) = mat_compare.col(1) - mat_compare.col(0);
//		cout << mat_compare << endl;
		cout << "relative error norm: " << mat_compare.col(2).norm() / mat_compare.col(0).norm() << " of snapshot number " << iter_index << endl;

		if (debug_mode)
		{
//			cout << "snapshot_x_collection.num_elements() " << snapshot_x_collection.num_elements() << " snapshot_x_collection[0].num_elements() " << snapshot_x_collection[0].num_elements() << endl;
		}

		if (write_ROM_field || (qoi_dof >= 0))
		{
			recover_snapshot_data(reconstruct_solution, current_index);  // this is setting the fields
		}

		Array<OneD, double> field_x;
		Array<OneD, double> field_y;
		recover_snapshot_loop(reconstruct_solution, field_x, field_y);
		double rel_ITHACA_L2error = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection[current_index], snapshot_y_collection[current_index]) / L2norm_ITHACA(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		cout << "relative L2 error : " << rel_ITHACA_L2error << " of snapshot number " << iter_index << endl;
		collected_relative_L2errors_snaps(current_index) = rel_ITHACA_L2error;



	} // for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	if (1)
	{
		std::stringstream sstm_VV;
		sstm_VV << "LocROM_cluster_snap" << current_cluster_number << ".txt";
		std::string LocROM_txt_VV = sstm_VV.str();
		const char* outname_VV = LocROM_txt_VV.c_str();
		ofstream myfile_VV (outname_VV);
		if (myfile_VV.is_open())
		{
			for (int iter_index = 0; iter_index < Nmax; ++iter_index)
			{
				myfile_VV << std::setprecision(17) << collected_relative_L2errors_snaps(iter_index) << "\t";
			}
			myfile_VV.close();
		}
		else cout << "Unable to open file"; 

	}
	if (use_fine_grid_VV)
	{
		// repeat the evaluation without the accuracy check
		// fine_general_param_vector is available already
		// start sweeping 
		Eigen::MatrixXd collected_qoi = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_L2errors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_Linferrors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_L2errors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_Linferrors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		int fine_grid_dir0_index = 0;
		int fine_grid_dir1_index = 0;
		double locROM_qoi;
		for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		{
//			cout << "fine_grid_iter_index " << iter_index << " of max " << fine_grid_dir0*fine_grid_dir1 << endl;
			int current_index = iter_index;
			double current_nu;
			double w;
			Array<OneD, NekDouble> current_param = fine_general_param_vector[current_index];
			w = current_param[0];	
			current_nu = current_param[1];
			if (debug_mode)
			{
//				cout << " VV online phase current nu " << current_nu << endl;
//				cout << " VV online phase current w " << w << endl;
			}
			Set_m_kinvis( current_nu );
			if (use_Newton)
			{
				DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
			}
			Eigen::MatrixXd curr_xy_proj = project_onto_basis(cluster_mean_x, cluster_mean_y);
			Eigen::MatrixXd affine_mat_proj;
			Eigen::VectorXd affine_vec_proj;
			affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
			affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
			Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
			double relative_change_error;
			int no_iter=0;
			Array<OneD, double> field_x;
			Array<OneD, double> field_y;
			// now start looping
			do
			{
				// for now only Oseen // otherwise need to do the DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
				Eigen::VectorXd prev_solve_affine = solve_affine;
				Eigen::VectorXd repro_solve_affine = RB * solve_affine;
				Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);

				recover_snapshot_loop(reconstruct_solution, field_x, field_y);
				if (use_Newton)
				{
					DoInitialiseAdv(field_x, field_y);
				}
				curr_xy_proj = project_onto_basis(field_x, field_y);
				if (parameter_space_dimension == 1)
				{
					affine_mat_proj = gen_affine_mat_proj(current_nu);
					affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
				}
				else if (parameter_space_dimension == 2)
				{
					affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
					affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
				}
				solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
				relative_change_error = (solve_affine - prev_solve_affine).norm() / prev_solve_affine.norm();
//				cout << "relative_change_error " << relative_change_error << endl;
				no_iter++;
			} 
			while( ((relative_change_error > 1e-5) && (no_iter < 100)) );
//			cout << "ROM solve no iters used " << no_iter << endl;
			Eigen::VectorXd repro_solve_affine = RB * solve_affine;
			Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
			if (write_ROM_field || (qoi_dof >= 0))
			{
				locROM_qoi = recover_snapshot_data(reconstruct_solution, 0);
			}
			collected_qoi(fine_grid_dir0_index, fine_grid_dir1_index) = locROM_qoi;
			if (use_fine_grid_VV_and_load_ref)
			{
				collected_relative_L2errors(fine_grid_dir0_index, fine_grid_dir1_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
				collected_relative_Linferrors(fine_grid_dir0_index, fine_grid_dir1_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
				if (use_non_unique_up_to_two)
				{
					collected_relative_L2errors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
					collected_relative_Linferrors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
				}
			}
			fine_grid_dir1_index++;
			if (fine_grid_dir1_index == fine_grid_dir1)
			{
				fine_grid_dir1_index = 0;
				fine_grid_dir0_index++;
			}			
		} // for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		std::stringstream sstm;
		sstm << "LocROM_cluster" << current_cluster_number << ".txt";
		std::string LocROM_txt = sstm.str();
		const char* outname = LocROM_txt.c_str();
		ofstream myfile (outname);
		if (myfile.is_open())
		{
			for (int i0 = 0; i0 < fine_grid_dir0; i0++)
			{
				for (int i1 = 0; i1 < fine_grid_dir1; i1++)
				{
					myfile << std::setprecision(17) << collected_qoi(i0,i1) << "\t";
				}
				myfile << "\n";
			}
			myfile.close();
		}
		else cout << "Unable to open file"; 

		if (use_fine_grid_VV_and_load_ref)
		{

			std::stringstream sstm_VV;
			sstm_VV << "LocROM_cluster_VV" << current_cluster_number << ".txt";
			std::string LocROM_txt_VV = sstm_VV.str();
			const char* outname_VV = LocROM_txt_VV.c_str();
			ofstream myfile_VV (outname_VV);
			if (myfile_VV.is_open())
			{
				for (int i0 = 0; i0 < fine_grid_dir0; i0++)
				{
					for (int i1 = 0; i1 < fine_grid_dir1; i1++)
					{
						myfile_VV << std::setprecision(17) << collected_relative_L2errors(i0,i1) << "\t";
					}
					myfile_VV << "\n";
				}
				myfile_VV.close();
			}
			else cout << "Unable to open file"; 

			std::stringstream sstm_VV_Linf;
			sstm_VV_Linf << "LocROM_cluster_VV_Linf" << current_cluster_number << ".txt";
			std::string LocROM_txt_VV_Linf = sstm_VV_Linf.str();
			const char* outname_VV_Linf = LocROM_txt_VV_Linf.c_str();
			ofstream myfile_VV_Linf (outname_VV_Linf);
			if (myfile_VV_Linf.is_open())
			{
				for (int i0 = 0; i0 < fine_grid_dir0; i0++)
				{
					for (int i1 = 0; i1 < fine_grid_dir1; i1++)
					{
						myfile_VV_Linf << std::setprecision(17) << collected_relative_Linferrors(i0,i1) << "\t";
					}
					myfile_VV_Linf << "\n";
				}
				myfile_VV_Linf.close();
			}
			else cout << "Unable to open file"; 
			
			if (use_non_unique_up_to_two)
			{
				std::stringstream sstm_VV;
				sstm_VV << "LocROM_cluster_VV_v2" << current_cluster_number << ".txt";
				std::string LocROM_txt_VV = sstm_VV.str();
				const char* outname_VV = LocROM_txt_VV.c_str();
				ofstream myfile_VV (outname_VV);
				if (myfile_VV.is_open())
				{
					for (int i0 = 0; i0 < fine_grid_dir0; i0++)
					{
						for (int i1 = 0; i1 < fine_grid_dir1; i1++)
						{
							myfile_VV << std::setprecision(17) << collected_relative_L2errors_v2(i0,i1) << "\t";
						}
						myfile_VV << "\n";
					}
					myfile_VV.close();
				}
				else cout << "Unable to open file"; 

				std::stringstream sstm_VV_Linf;
				sstm_VV_Linf << "LocROM_cluster_VV_Linf_v2" << current_cluster_number << ".txt";
				std::string LocROM_txt_VV_Linf = sstm_VV_Linf.str();
				const char* outname_VV_Linf = LocROM_txt_VV_Linf.c_str();
				ofstream myfile_VV_Linf (outname_VV_Linf);
				if (myfile_VV_Linf.is_open())
				{
					for (int i0 = 0; i0 < fine_grid_dir0; i0++)
					{
						for (int i1 = 0; i1 < fine_grid_dir1; i1++)
						{
							myfile_VV_Linf << std::setprecision(17) << collected_relative_Linferrors_v2(i0,i1) << "\t";
						}
						myfile_VV_Linf << "\n";
					}
					myfile_VV_Linf.close();
				}
				else cout << "Unable to open file"; 
			}

		}

	}
    }

    void CoupledLinearNS_TT::associate_VV_to_clusters(Array<OneD, std::set<int> > clusters)
    {
	if (1)
	{
	// for each VV point identify the next closest cluster snapshot
	Eigen::MatrixXd VV_cluster_association = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
	Eigen::MatrixXd VV_cluster_association_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1); // for if use_non_unique_up_to_two
	int fine_general_param_vector_index = 0;
	int no_clusters = clusters.num_elements();
	for (int i0 = 0; i0 < fine_grid_dir0; i0++)
	{
		for (int i1 = 0; i1 < fine_grid_dir1; i1++)
		{
			// find next snapshot location o_i
			int general_param_vector_index = find_closest_snapshot_location(fine_general_param_vector[fine_general_param_vector_index], general_param_vector);
			fine_general_param_vector_index++;
			// identify cluster in which o_i is
			for (int j = 0; j < no_clusters; ++j)
			{
				if (clusters[j].count(general_param_vector_index) == 1)
					VV_cluster_association(i0, i1) = j;
				if (use_non_unique_up_to_two)
				{
					if (clusters[j].count(general_param_vector_index + Nmax/2) == 1)
						VV_cluster_association_v2(i0, i1) = j;
				}
			}
		}
	}
	std::stringstream sstm;
	sstm << "VV_cluster_association.txt";
	std::string LocROM_txt = sstm.str();
	const char* outname = LocROM_txt.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < fine_grid_dir0; i0++)
		{
			for (int i1 = 0; i1 < fine_grid_dir1; i1++)
			{
				myfile << VV_cluster_association(i0,i1) << "\t";
			}
			myfile << "\n";
		}
		if (use_non_unique_up_to_two)
		{
			for (int i0 = 0; i0 < fine_grid_dir0; i0++)
			{
				for (int i1 = 0; i1 < fine_grid_dir1; i1++)
				{
					myfile << VV_cluster_association_v2(i0,i1) << "\t";
				}
				myfile << "\n";
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 	
	}
	if (1)
	{
	// for each VV point identify the next closest cluster snapshot
	Eigen::MatrixXd VV_cluster_association = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
	int fine_general_param_vector_index = 0;
	int no_clusters = clusters.num_elements();
	for (int i0 = 0; i0 < fine_grid_dir0; i0++)
	{
		for (int i1 = 0; i1 < fine_grid_dir1; i1++)
		{
			// find next snapshot location o_i
			int general_param_vector_index = find_closest_snapshot_location_l1(fine_general_param_vector[fine_general_param_vector_index], general_param_vector);
			fine_general_param_vector_index++;
			// identify cluster in which o_i is
			for (int j = 0; j < no_clusters; ++j)
			{
				if (clusters[j].count(general_param_vector_index) == 1)
					VV_cluster_association(i0, i1) = j;
			}
		}
	}
	std::stringstream sstm;
	sstm << "VV_cluster_association_l1.txt";
	std::string LocROM_txt = sstm.str();
	const char* outname = LocROM_txt.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < fine_grid_dir0; i0++)
		{
			for (int i1 = 0; i1 < fine_grid_dir1; i1++)
			{
				myfile << VV_cluster_association(i0,i1) << "\t";
			}
			myfile << "\n";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 	
	}
	if (1)
	{
	// for each VV point identify the next closest cluster snapshot
	Eigen::MatrixXd VV_cluster_association = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
	int fine_general_param_vector_index = 0;
	int no_clusters = clusters.num_elements();
	for (int i0 = 0; i0 < fine_grid_dir0; i0++)
	{
		for (int i1 = 0; i1 < fine_grid_dir1; i1++)
		{
			// find next snapshot location o_i
			int general_param_vector_index = find_closest_snapshot_location_linf(fine_general_param_vector[fine_general_param_vector_index], general_param_vector);
			fine_general_param_vector_index++;
			// identify cluster in which o_i is
			for (int j = 0; j < no_clusters; ++j)
			{
				if (clusters[j].count(general_param_vector_index) == 1)
					VV_cluster_association(i0, i1) = j;
			}
		}
	}
	std::stringstream sstm;
	sstm << "VV_cluster_association_linf.txt";
	std::string LocROM_txt = sstm.str();
	const char* outname = LocROM_txt.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int i0 = 0; i0 < fine_grid_dir0; i0++)
		{
			for (int i1 = 0; i1 < fine_grid_dir1; i1++)
			{
				myfile << VV_cluster_association(i0,i1) << "\t";
			}
			myfile << "\n";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 	
	}

    }

    int CoupledLinearNS_TT::find_closest_snapshot_location(Array<OneD, NekDouble> VV_point, Array<OneD, Array<OneD, NekDouble> > general_param_vector)
    {
	double min_distance;
	double min_distance_index;
	for (int i = 0; i < general_param_vector.num_elements(); ++i)
	{
		// elementary euclidean distances
		double distance = sqrt( (VV_point[0] - general_param_vector[i][0])*(VV_point[0] - general_param_vector[i][0]) + (VV_point[1] - general_param_vector[i][1])*(VV_point[1] - general_param_vector[i][1]) );
		if (i == 0)
		{
			min_distance = distance;
			min_distance_index = 0;
		}
		else if (distance < min_distance)
		{
			min_distance = distance;
			min_distance_index = i;
		}
	}
	return min_distance_index;
    }

    int CoupledLinearNS_TT::find_closest_snapshot_location_l1(Array<OneD, NekDouble> VV_point, Array<OneD, Array<OneD, NekDouble> > general_param_vector)
    {
	double min_distance;
	double min_distance_index;
	for (int i = 0; i < general_param_vector.num_elements(); ++i)
	{
		double distance = abs(VV_point[0] - general_param_vector[i][0]) + abs(VV_point[1] - general_param_vector[i][1]);
		if (i == 0)
		{
			min_distance = distance;
			min_distance_index = 0;
		}
		else if (distance < min_distance)
		{
			min_distance = distance;
			min_distance_index = i;
		}
	}
	return min_distance_index;
    }

    int CoupledLinearNS_TT::find_closest_snapshot_location_linf(Array<OneD, NekDouble> VV_point, Array<OneD, Array<OneD, NekDouble> > general_param_vector)
    {
	double min_distance;
	double min_distance_index;
	for (int i = 0; i < general_param_vector.num_elements(); ++i)
	{
		std::vector<double> dist_vector = std::vector<double> (2);
		dist_vector[0] = (VV_point[0] - general_param_vector[i][0]);
		dist_vector[1] = (VV_point[1] - general_param_vector[i][1]);
		double distance = max( dist_vector[0], dist_vector[1] );
		if (i == 0)
		{
			min_distance = distance;
			min_distance_index = 0;
		}
		else if (distance < min_distance)
		{
			min_distance = distance;
			min_distance_index = i;
		}
	}
	return min_distance_index;
    }

    double CoupledLinearNS_TT::lagrange_interp(double curr_param, int curr_index, int sparse_poly_approx_dimension)
	{
		double lagrange_value = 1;
		for (int i = 0; i < sparse_poly_approx_dimension; ++i)
		{
			if (general_param_vector[curr_index][0] != general_param_vector[i][0])
				lagrange_value *= (curr_param - general_param_vector[i][0]) / (general_param_vector[curr_index][0] - general_param_vector[i][0]);
		}
		return lagrange_value;
	}

	void CoupledLinearNS_TT::sparse_approx_VV(int sparse_poly_approx_dimension, double& max, double& mean)
	{
		Array<OneD, NekDouble> collect_rel_L2error(fine_grid_dir0*fine_grid_dir1);
		for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		{
			int current_index = iter_index;
			Array<OneD, NekDouble> current_param = fine_general_param_vector[current_index];
			double w = current_param[0];	
			double current_nu = current_param[1];
			Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].num_elements());
			Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].num_elements());
			for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
			{
				interpolant_x[i] = 0.0;
				interpolant_y[i] = 0.0;
			}
			for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
			{
				double lagrange_value = lagrange_interp(w, index_interpol_op, sparse_poly_approx_dimension);
				for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)	
				{
					interpolant_x[i] += snapshot_x_collection[index_interpol_op][i] * lagrange_value;
					interpolant_y[i] += snapshot_y_collection[index_interpol_op][i] * lagrange_value;
				}
			}
			double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
			collect_rel_L2error[iter_index] = rel_L2error;
		}
		mean = 0;
		max = 0;
		for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		{
			mean += collect_rel_L2error[iter_index] / (fine_grid_dir0*fine_grid_dir1);
			if (collect_rel_L2error[iter_index] > max)
				max = collect_rel_L2error[iter_index];
		}
	}


    void CoupledLinearNS_TT::compute_sparse_poly_approx()
	{
		int sparse_poly_approx_dimension = max_sparse_poly_approx_dimension;
		// L2 error works on the snapshot_x_collection and snapshot_y_collection
		// start sweeping 
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			int current_index = iter_index;
			double current_nu;
			double w;
			if (parameter_space_dimension == 1)
			{
				current_nu = param_vector[current_index];
			}
			else if (parameter_space_dimension == 2)
			{
				Array<OneD, NekDouble> current_param = general_param_vector[current_index];
				w = current_param[0];	
				current_nu = current_param[1];
			}		
			// only case of w parameter
			Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].num_elements());
			Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].num_elements());
			for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
			{
				interpolant_x[i] = 0.0;
				interpolant_y[i] = 0.0;
			}
			for (int index_interpol_op = 0; index_interpol_op < sparse_poly_approx_dimension; ++index_interpol_op)
			{
				double lagrange_value = lagrange_interp(w, index_interpol_op, sparse_poly_approx_dimension);
				for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)	
				{
					interpolant_x[i] += snapshot_x_collection[index_interpol_op][i] * lagrange_value;
					interpolant_y[i] += snapshot_y_collection[index_interpol_op][i] * lagrange_value;
				}
			}
			double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
			cout << "rel_L2error at parameter " << w << " is " << rel_L2error << endl;
		}
		Array<OneD, NekDouble> collect_max(sparse_poly_approx_dimension);
		Array<OneD, NekDouble> collect_mean(sparse_poly_approx_dimension);
		double max, mean;
		for (int approx_dim = 1; approx_dim <= sparse_poly_approx_dimension; ++approx_dim)
		{
			sparse_approx_VV(approx_dim, max, mean);
			cout << "max at dim " << approx_dim << " is " << max << endl;
			cout << "mean at dim " << approx_dim << " is " << mean << endl;
			collect_max[approx_dim-1] = max;
			collect_mean[approx_dim-1] = mean;
		}
		{
		std::stringstream sstm;
		sstm << "sparse_conv_mean.txt";
		std::string sparse_conv_mean = sstm.str();
		const char* outname = sparse_conv_mean.c_str();
		ofstream myfile (outname);
		if (myfile.is_open())
		{
			for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
			{
				myfile << std::setprecision(17) << collect_mean[i0] << "\t";
			}
			myfile.close();
		}
		else cout << "Unable to open file"; 
		}
		{
		std::stringstream sstm;
		sstm << "sparse_conv_max.txt";
		std::string sparse_conv_max = sstm.str();
		const char* outname = sparse_conv_max.c_str();
		ofstream myfile (outname);
		if (myfile.is_open())
		{
			for (int i0 = 0; i0 < sparse_poly_approx_dimension; i0++)
			{
				myfile << std::setprecision(17) << collect_max[i0] << "\t";
			}
			myfile.close();
		}
		else cout << "Unable to open file"; 
		}
	}


    void CoupledLinearNS_TT::compute_ANN_approx()
	{
		ANN_POD_coeffs = Eigen::MatrixXd::Zero(fine_grid_dir0*fine_grid_dir1, RBsize*2);
		std::string predANN_txt = "evaluate/pred_fsg.txt";
		const char* predANN_txt_t = predANN_txt.c_str();
		ifstream myfile_predANN_txt_t (predANN_txt_t);
		std::vector< std::vector<double> > all_double;
		if (myfile_predANN_txt_t.is_open())
		{
			std::string line;
			std::vector< std::vector<double> > all_double;
			int counter = 0;
			while ( getline( myfile_predANN_txt_t, line ) ) 
			{
			//	cout << line << endl;
				std::istringstream is( line );
				std::vector<double> nn = std::vector<double>( std::istream_iterator<double>(is), std::istream_iterator<double>() );
			//	cout << "nn.size() "  << nn.size() << endl;
				for (int i = 0; i < nn.size(); ++i)
				{
					//optimal_clusters[counter].insert(nn[i]);
					//cout << nn[i] << endl;
					//cout << "counter " << counter << endl;
					ANN_POD_coeffs(counter,i) = nn[i];
				}
				++counter;
//				      all_integers.push_back( std::vector<int>( std::istream_iterator<int>(is), std::istream_iterator<int>() ) );
			}

	/*			for (int i = 0; i < no_clusters; ++i)
				{
					for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
					{
						myfile_optimal_clustering_txt_t << *it << "\t";
					}
					myfile_optimal_clustering_txt_t << endl;
				} */
			myfile_predANN_txt_t.close(); 
		}
		else cout << "Unable to open file evaluate/pred_fsg.txt"; 
		// L2 error works on the snapshot_x_collection and snapshot_y_collection
		// start sweeping 
		double mean_L2 = 0;
		double max_L2 = 0;
		double mean_Linf = 0;
		double max_Linf = 0;
		for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		{
			int current_index = iter_index;
			double current_nu;
			double w;
			if (parameter_space_dimension == 1)
			{
				current_nu = param_vector[current_index];
			}
			else if (parameter_space_dimension == 2)
			{
				Array<OneD, NekDouble> current_param = fine_general_param_vector[current_index];
				w = current_param[0];	
				current_nu = current_param[1];
			}		
			Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].num_elements());
			Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].num_elements());
			for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
			{
				interpolant_x[i] = 0.0;
				interpolant_y[i] = 0.0;
			}
			for (int index_interpol_op = 0; index_interpol_op < RBsize; ++index_interpol_op)
			{
				for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)	
				{
					interpolant_x[i] += eigen_phys_basis_x(i, index_interpol_op) * ANN_POD_coeffs(iter_index, index_interpol_op);
					interpolant_y[i] += eigen_phys_basis_y(i, index_interpol_op) * ANN_POD_coeffs(iter_index, index_interpol_op + RBsize);
				}
			}
			double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
			double rel_Linferror = Linfnorm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
		//	cout << "rel_L2error at parameter " << w << " and " << current_nu << " is " << rel_L2error << endl;
			if (rel_L2error > max_L2) max_L2 = rel_L2error;
			mean_L2 += rel_L2error / (fine_grid_dir0*fine_grid_dir1);
			if (rel_Linferror > max_Linf) max_Linf = rel_Linferror;
			mean_Linf += rel_Linferror / (fine_grid_dir0*fine_grid_dir1);
		}
		cout << "mean_L2 ann " << mean_L2 << endl;
		cout << "max_L2 ann " << max_L2 << endl;
		cout << "mean_Linf ann " << mean_Linf << endl;
		cout << "max_Linf ann " << max_Linf << endl;
/*		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			cout << "test" << endl;
			int current_index = iter_index;
			double current_nu;
			double w;
			if (parameter_space_dimension == 1)
			{
				current_nu = param_vector[current_index];
			}
			else if (parameter_space_dimension == 2)
			{
				Array<OneD, NekDouble> current_param = general_param_vector[current_index];
				w = current_param[0];	
				current_nu = current_param[1];
			}		
			Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].num_elements());
			Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].num_elements());
			for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
			{
				interpolant_x[i] = 0.0;
				interpolant_y[i] = 0.0;
			}
			for (int index_interpol_op = 0; index_interpol_op < RBsize; ++index_interpol_op)
			{
				for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)	
				{
					interpolant_x[i] += eigen_phys_basis_x(i, index_interpol_op) * train_data_x(index_interpol_op, iter_index);
					interpolant_y[i] += eigen_phys_basis_y(i, index_interpol_op) * train_data_y(index_interpol_op, iter_index);
				}
			}
			double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]) / L2norm_ITHACA(snapshot_x_collection[iter_index], snapshot_y_collection[iter_index]);
			cout << "rel_L2error at parameter " << w << " and " << current_nu << " is " << rel_L2error << endl;
		}
		
		cout << "comparing the first coeffs, ANN:" << endl;
		for (int index_interpol_op = 0; index_interpol_op < RBsize; ++index_interpol_op)
		{
			cout << "ANN_POD_coeffs(0, index_interpol_op) " << ANN_POD_coeffs(0, index_interpol_op) << endl;
			cout << "ANN_POD_coeffs(0, index_interpol_op + RBsize) " << ANN_POD_coeffs(0, index_interpol_op + RBsize) << endl;
			cout << "train_data_x(index_interpol_op, 0) " << train_data_x(index_interpol_op, 0) << endl;
			cout << "train_data_y(index_interpol_op, 0) " << train_data_y(index_interpol_op, 0) << endl;
		}   */
		
	}


    void CoupledLinearNS_TT::compute_ANN_approx_cluster(int cluster_no)
	{
		cout << "start CoupledLinearNS_TT::compute_ANN_approx_cluster " << endl;
		Eigen::MatrixXd local_ANN_POD_coeffs = Eigen::MatrixXd::Zero(fine_grid_dir0*fine_grid_dir1, RBsize*2);
		cout << "size1 " << fine_grid_dir0*fine_grid_dir1 << endl;
		cout << "size2 " << RBsize*2 << endl;
		std::stringstream sstm;
		sstm << "pred_fsg_" << cluster_no << ".txt";
		std::string predANN_txt = sstm.str();
	//	std::string predANN_txt = "evaluate/pred_fsg.txt";
		const char* predANN_txt_t = predANN_txt.c_str();
		ifstream myfile_predANN_txt_t (predANN_txt_t);
		std::vector< std::vector<double> > all_double;
		if (myfile_predANN_txt_t.is_open())
		{
			std::string line;
			std::vector< std::vector<double> > all_double;
			int counter = 0;
			while ( getline( myfile_predANN_txt_t, line ) ) 
			{
			//	cout << line << endl;
				std::istringstream is( line );
				std::vector<double> nn = std::vector<double>( std::istream_iterator<double>(is), std::istream_iterator<double>() );
			//	cout << "nn.size() "  << nn.size() << endl;
				for (int i = 0; i < nn.size(); ++i)
				{
					//optimal_clusters[counter].insert(nn[i]);
					//cout << nn[i] << endl;
				//	cout << "counter " << counter << endl;
				//	cout << "i " << i << endl;
				//	cout << "size1 " << fine_grid_dir0*fine_grid_dir1 << endl;
				//	cout << "size2 " << RBsize*2 << endl;
					local_ANN_POD_coeffs(counter,i) = nn[i];
				}
				++counter;
//				      all_integers.push_back( std::vector<int>( std::istream_iterator<int>(is), std::istream_iterator<int>() ) );
			}

	/*			for (int i = 0; i < no_clusters; ++i)
				{
					for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
					{
						myfile_optimal_clustering_txt_t << *it << "\t";
					}
					myfile_optimal_clustering_txt_t << endl;
				} */
			myfile_predANN_txt_t.close(); 
		}
		else cout << "Unable to open file pred_fsg for current cluster"; 
		cout << "prediction loaded " << endl;
		// L2 error works on the snapshot_x_collection and snapshot_y_collection
		// start sweeping 
		double fine_grid_dir1_index = 0;
		double fine_grid_dir0_index = 0;
		Eigen::MatrixXd collected_relative_L2errors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_Linferrors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		{
			int current_index = iter_index;
			double current_nu;
			double w;
			if (parameter_space_dimension == 1)
			{
				current_nu = param_vector[current_index];
			}
			else if (parameter_space_dimension == 2)
			{
				Array<OneD, NekDouble> current_param = fine_general_param_vector[current_index];
				w = current_param[0];	
				current_nu = current_param[1];
			}		
			Array<OneD, NekDouble> interpolant_x(snapshot_x_collection[0].num_elements());
			Array<OneD, NekDouble> interpolant_y(snapshot_y_collection[0].num_elements());
			for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
			{
				interpolant_x[i] = 0.0;
				interpolant_y[i] = 0.0;
			}
			for (int index_interpol_op = 0; index_interpol_op < RBsize; ++index_interpol_op)
			{
				for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)	
				{
					interpolant_x[i] += eigen_phys_basis_x(i, index_interpol_op) * local_ANN_POD_coeffs(iter_index, index_interpol_op);
					interpolant_y[i] += eigen_phys_basis_y(i, index_interpol_op) * local_ANN_POD_coeffs(iter_index, index_interpol_op + RBsize);
				}
			}
			double rel_L2error = L2norm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
			double rel_Linferror = Linfnorm_abs_error_ITHACA(interpolant_x, interpolant_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
			collected_relative_L2errors(fine_grid_dir0_index, fine_grid_dir1_index) =  rel_L2error;
			collected_relative_Linferrors(fine_grid_dir0_index, fine_grid_dir1_index) =  rel_Linferror;
			fine_grid_dir1_index++;
			if (fine_grid_dir1_index == fine_grid_dir1)
			{
				fine_grid_dir1_index = 0;
				fine_grid_dir0_index++;
			}
		}
		std::stringstream sstm_L2;
		sstm_L2 << "POD_NN_local_L2_" << cluster_no << ".txt";
		std::string LocROM_sstm_L2 = sstm_L2.str();
		const char* outname_VV = LocROM_sstm_L2.c_str();
		ofstream myfile_VV (outname_VV);
		if (myfile_VV.is_open())
		{
			for (int i0 = 0; i0 < fine_grid_dir0; i0++)
			{
				for (int i1 = 0; i1 < fine_grid_dir1; i1++)
				{
					myfile_VV << std::setprecision(17) << collected_relative_L2errors(i0,i1) << "\t";
				}
				myfile_VV << "\n";
			}
			myfile_VV.close();
		}
		else cout << "Unable to open file"; 
		std::stringstream sstm_Linf;
		sstm_Linf << "POD_NN_local_Linf_" << cluster_no << ".txt";
		std::string LocROM_sstm_Linf = sstm_Linf.str();
		const char* outname_VV_Linf = LocROM_sstm_Linf.c_str();
		ofstream myfile_VV_Linf (outname_VV_Linf);
		if (myfile_VV_Linf.is_open())
		{
			for (int i0 = 0; i0 < fine_grid_dir0; i0++)
			{
				for (int i1 = 0; i1 < fine_grid_dir1; i1++)
				{
					myfile_VV_Linf << std::setprecision(17) << collected_relative_Linferrors(i0,i1) << "\t";
				}
				myfile_VV_Linf<< "\n";
			}
			myfile_VV_Linf.close();
		}
		else cout << "Unable to open file"; 
	}


    void CoupledLinearNS_TT::online_phase()
    {

	if (use_sparse_poly)
	{
		compute_sparse_poly_approx();
		return;
	}

	if (use_ANN)
	{
		compute_ANN_approx();
		return;
	}
	if (use_ANN_local)
	{
		//compute_ANN_approx();
		return;
	}

	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(), 3);  // is of size M_truth_size
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double w;
		if (parameter_space_dimension == 1)
		{
			current_nu = param_vector[current_index];
		}
		else if (parameter_space_dimension == 2)
		{
			Array<OneD, NekDouble> current_param = general_param_vector[current_index];
			w = current_param[0];	
			current_nu = current_param[1];
		}
		if (debug_mode)
		{
			cout << " online phase current nu " << current_nu << endl;
			cout << " online phase current w " << w << endl;
		}
		Set_m_kinvis( current_nu );
		if (use_Newton)
		{
			DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		}

		Eigen::MatrixXd curr_xy_proj = project_onto_basis(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		Eigen::MatrixXd affine_mat_proj;
		Eigen::VectorXd affine_vec_proj;

		if (parameter_space_dimension == 1)
		{
			affine_mat_proj = gen_affine_mat_proj(current_nu);
			affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
//			cout << "aff mat 1d " << affine_mat_proj << endl;
//			cout << "aff vec 1d " << affine_vec_proj << endl;
		}
		else if (parameter_space_dimension == 2)
		{
			affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
//			cout << "aff mat 2d " << affine_mat_proj << endl;
			affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
//			cout << "aff vec 2d " << affine_vec_proj << endl;
//			Eigen::MatrixXd affine_mat_proj_1d = gen_affine_mat_proj(current_nu);
//			Eigen::VectorXd affine_vec_proj_1d = gen_affine_vec_proj(current_nu, current_index);
		}

		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
//		cout << "solve_affine " << solve_affine << endl;
		Eigen::VectorXd repro_solve_affine = RB * solve_affine;
		Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
		if (globally_connected == 1)
		{
			mat_compare.col(0) = M_collect_f_all.col(current_index);
		}
		else
		{
			mat_compare.col(0) = collect_f_all.col(current_index);
			if (debug_mode)
			{
				Eigen::VectorXd current_f_all = Eigen::VectorXd::Zero(collect_f_all.rows());
				current_f_all = collect_f_all.col(current_index);
				Eigen::VectorXd current_f_all_wo_dbc = remove_rows(current_f_all, elem_loc_dbc);
				Eigen::VectorXd proj_current_f_all_wo_dbc = RB.transpose() * current_f_all_wo_dbc;
//				cout << "proj_current_f_all_wo_dbc " << proj_current_f_all_wo_dbc << endl;
				Eigen::VectorXd correctRHS = affine_mat_proj * proj_current_f_all_wo_dbc;
				Eigen::VectorXd correction_RHS = correctRHS - affine_vec_proj;
//				cout << "correctRHS " << correctRHS << endl;
//				cout << "correction_RHS " << correction_RHS << endl;
			}
		}
		mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
		mat_compare.col(2) = mat_compare.col(1) - mat_compare.col(0);
//		cout << mat_compare << endl;
//		cout << "relative euclidean error norm: " << mat_compare.col(2).norm() / mat_compare.col(0).norm() << " of snapshot number " << iter_index << endl;

		// compare this to the projection error onto RB, cannot use collect_f_all, is numerically unstable
/*		cout << "mat_compare.cols() " << mat_compare.cols() << endl;
		cout << "mat_compare.rows() " << mat_compare.rows() << endl;
		cout << "collect_f_all.cols() " << collect_f_all.cols() << endl;
		cout << "collect_f_all.rows() " << collect_f_all.rows() << endl;
		cout << "RB.cols() " << RB.cols() << endl;
		cout << "RB.rows() " << RB.rows() << endl; */

		
//		cout << "curr_xy_proj.cols() " << curr_xy_proj.cols() << endl;
	//	cout << "curr_xy_proj.rows() " << curr_xy_proj.rows() << endl; 


//		Eigen::VectorXd proj_solution = collect_f_all.transpose() * mat_compare.col(0);
//		Eigen::VectorXd reproj_solution = collect_f_all * proj_solution;
		Eigen::VectorXd FOM_solution = mat_compare.col(0);
		Eigen::VectorXd FOM_solution_wo_dbc = remove_rows(FOM_solution, elem_loc_dbc);
		Eigen::VectorXd proj_FOM_solution_wo_dbc = RB.transpose() * FOM_solution_wo_dbc;
		Eigen::VectorXd reproj_FOM_solution_wo_dbc = RB * proj_FOM_solution_wo_dbc;
		Eigen::VectorXd reconstruct_FOM_solution = reconstruct_solution_w_dbc(reproj_FOM_solution_wo_dbc);

//		cout << "reconstruct_FOM_solution.norm(): " << reconstruct_FOM_solution.norm() << endl;
//		cout << "reconstruct_solution.norm(): " << reconstruct_solution.norm() << endl;
		Eigen::VectorXd diff_projection = reconstruct_FOM_solution - mat_compare.col(0);

//		cout << "relative euclidean projection error norm: " << diff_projection.norm() / mat_compare.col(0).norm() << " of snapshot number " << iter_index << endl;


		// now only in RB:
		Eigen::VectorXd diff_projection_RB = reproj_FOM_solution_wo_dbc - remove_rows(mat_compare.col(0), elem_loc_dbc);
		Eigen::VectorXd diff_RB = repro_solve_affine - remove_rows(mat_compare.col(0), elem_loc_dbc);
//		cout << "relative euclidean RB projection error norm: " << diff_projection_RB.norm() / remove_rows(mat_compare.col(0), elem_loc_dbc).norm() << " of snapshot number " << iter_index << endl;
//		cout << "relative euclidean RB error norm: " << diff_RB.norm() / remove_rows(mat_compare.col(0), elem_loc_dbc).norm() << " of snapshot number " << iter_index << endl;

		// have to use curr_xy_proj for better approximations


		if (debug_mode)
		{
			cout << "snapshot_x_collection.num_elements() " << snapshot_x_collection.num_elements() << " snapshot_x_collection[0].num_elements() " << snapshot_x_collection[0].num_elements() << endl;
		}

		if (write_ROM_field || (qoi_dof >= 0))
		{
			recover_snapshot_data(reconstruct_solution, current_index); // this is setting the fields and fieldcoeffs
		}

	  	Eigen::VectorXd f_bnd = reconstruct_solution.head(curr_f_bnd.size());
		Eigen::VectorXd f_int = reconstruct_solution.tail(curr_f_int.size());
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
		Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(GetNcoeffs());
		Array<OneD, double> field_1(GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
		int cnt = 0;
		int cnt1 = 0;
		int nvel = 2;
		int nz_loc = 1;
		int  nplanecoeffs = fields[0]->GetNcoeffs();
		int  nel  = m_fields[0]->GetNumElmts();
		for(int i = 0; i < nel; ++i) 
		{
		      int eid  = i;
		      fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
		      fields[0]->GetExp(eid)->GetInteriorMap(imap);
		      int nbnd   = bmap.num_elements();
		      int nint   = imap.num_elements();
			      int offset = fields[0]->GetCoeff_Offset(eid);
		            
		      for(int j = 0; j < nvel; ++j)
		      {
		           for(int n = 0; n < nz_loc; ++n)
		           {
		                    for(int k = 0; k < nbnd; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
		                    }
		                    
		                    for(int k = 0; k < nint; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
		                    }
		                    cnt  += nbnd;
		                    cnt1 += nint;
		           }
		      }
		}
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);

        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.num_elements()+1);
        int i;
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
	    }
//		cout << "after recover_snapshot_data have fieldcoeffs[0].num_elements() " << fieldcoeffs[0].num_elements() << endl;

//		cout << "after recover_snapshot_data have curr_PhysBaseVec_x.num_elements() " << curr_PhysBaseVec_x.num_elements() << endl;
//		cout << "after recover_snapshot_data have snapshot_x_collection[current_index].num_elements() " << snapshot_x_collection[current_index].num_elements() << endl;

		// have the FOM_snapshot_solution projection available as curr_xy_proj
		// datastructure: 	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection;
		Eigen::VectorXd diff_x_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		Eigen::VectorXd snap_x = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		Eigen::VectorXd diff_y_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		Eigen::VectorXd snap_y = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		for (int index_recr = 0; index_recr < curr_PhysBaseVec_x.num_elements(); ++index_recr)
		{
			snap_x(index_recr) = snapshot_x_collection[current_index][index_recr];
			snap_y(index_recr) = snapshot_y_collection[current_index][index_recr];
			diff_x_RB_solve(index_recr) = curr_PhysBaseVec_x[index_recr] - snapshot_x_collection[current_index][index_recr];
			diff_y_RB_solve(index_recr) = curr_PhysBaseVec_y[index_recr] - snapshot_y_collection[current_index][index_recr];
		}

		Eigen::MatrixXd curr_xy_reproj = reproject_from_basis(curr_xy_proj);

		cout << "relative euclidean error norm in x coords: " << diff_x_RB_solve.norm() / snap_x.norm() << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean error norm in y coords: " << diff_y_RB_solve.norm() / snap_y.norm() << " of snapshot number " << iter_index << endl;

//		cout << "curr_xy_reproj.cols() " << curr_xy_reproj.cols() << endl;
//		cout << "curr_xy_reproj.rows() " << curr_xy_reproj.rows() << endl;
//		cout << "snap_x.rows() " << snap_x.rows() << endl;

		Eigen::VectorXd diff_x_proj = curr_xy_reproj.col(0) - snap_x;
		Eigen::VectorXd diff_y_proj = curr_xy_reproj.col(1) - snap_y;
		
		cout << "relative euclidean projection error norm in x coords: " << diff_x_proj.norm() / snap_x.norm() << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean projection error norm in y coords: " << diff_y_proj.norm() / snap_y.norm() << " of snapshot number " << iter_index << endl;

	}

	if (compute_smaller_model_errs)
	{
		// repeat the parameter sweep with decreasing RB sizes up to 1, but in a separate function for readability
		for (int i=1; i < RBsize; ++i)
		{
    		online_snapshot_check_with_smaller_basis(i);
		}

	}

	if (use_fine_grid_VV)
	{
		// repeat the evaluation without the accuracy check
		// fine_general_param_vector is available already
		// start sweeping 

		// Question: how to init?
		// could use all-zero or the cluster-mean
		Array<OneD, NekDouble> cluster_mean_x(snapshot_x_collection[0].num_elements(), 0.0);
		Array<OneD, NekDouble> cluster_mean_y(snapshot_y_collection[0].num_elements(), 0.0);
/*		for (std::set<int>::iterator it=current_cluster.begin(); it!=current_cluster.end(); ++it)
		{
			for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
			{
				cluster_mean_x[i] += (1.0 / current_cluster.size()) * snapshot_x_collection[*it][i];
				cluster_mean_y[i] += (1.0 / current_cluster.size()) * snapshot_y_collection[*it][i];
			}
		} */

		Eigen::MatrixXd collected_qoi = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_L2errors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_Linferrors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_L2errors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_Linferrors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		int fine_grid_dir0_index = 0;
		int fine_grid_dir1_index = 0;
		double locROM_qoi;
		for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		{
//			cout << "fine_grid_iter_index " << iter_index << " of max " << fine_grid_dir0*fine_grid_dir1 << endl;
			int current_index = iter_index;
			double current_nu;
			double w;
			Array<OneD, NekDouble> current_param = fine_general_param_vector[current_index];
			w = current_param[0];	
			current_nu = current_param[1];
			if (debug_mode)
			{
	//			cout << " VV online phase current nu " << current_nu << endl;
	//			cout << " VV online phase current w " << w << endl;
			}
			Set_m_kinvis( current_nu );
			if (use_Newton)
			{
				DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
			}

			// test for debug purpose
	/*		cluster_mean_x = snapshot_x_collection_VV[iter_index];
			cluster_mean_y = snapshot_y_collection_VV[iter_index];
			Eigen::VectorXd snap_x = Eigen::VectorXd::Zero(cluster_mean_x.num_elements());
			Eigen::VectorXd snap_y = Eigen::VectorXd::Zero(cluster_mean_x.num_elements());
			for (int index_recr = 0; index_recr < cluster_mean_x.num_elements(); ++index_recr)
			{
				snap_x(index_recr) = snapshot_x_collection_VV[iter_index][index_recr];
				snap_y(index_recr) = snapshot_y_collection_VV[iter_index][index_recr];
			}
			cout << "snap_x.norm() " << snap_x.norm() << endl;
			cout << "snap_y.norm() " << snap_y.norm() << endl;   */
/*			for (int index_recr = 0; index_recr < cluster_mean_x.num_elements(); ++index_recr)
			{
				snap_x(index_recr) = snapshot_x_collection[iter_index][index_recr];
				snap_y(index_recr) = snapshot_y_collection[iter_index][index_recr];
			}
			cout << "sample snap_x.norm() " << snap_x.norm() << endl;
			cout << "sample snap_y.norm() " << snap_y.norm() << endl;
*/

			Eigen::MatrixXd curr_xy_proj = project_onto_basis(cluster_mean_x, cluster_mean_y);
			Eigen::MatrixXd affine_mat_proj;
			Eigen::VectorXd affine_vec_proj;
			affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
			affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
			Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
			double relative_change_error;
			int no_iter=0;
			Array<OneD, double> field_x;
			Array<OneD, double> field_y;
			// now start looping
			do
			{
				// for now only Oseen // otherwise need to do the DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
				Eigen::VectorXd prev_solve_affine = solve_affine;
				Eigen::VectorXd repro_solve_affine = RB * solve_affine;
				Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);

				recover_snapshot_loop(reconstruct_solution, field_x, field_y);
				if (use_Newton)
				{
					DoInitialiseAdv(field_x, field_y);
				}
				curr_xy_proj = project_onto_basis(field_x, field_y);
				if (parameter_space_dimension == 1)
				{
					affine_mat_proj = gen_affine_mat_proj(current_nu);
					affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
				}
				else if (parameter_space_dimension == 2)
				{
					//cout << " VV online phase current nu " << current_nu << endl;
					//cout << " VV online phase current w " << w << endl;
					affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
					affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
				}
				solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
				relative_change_error = (solve_affine - prev_solve_affine).norm() / prev_solve_affine.norm();
				//cout << "relative_change_error " << relative_change_error << " no_iter " << no_iter << endl;
				no_iter++;
			} 
			while( ((relative_change_error > 1e-12) && (no_iter < 100)) );
//			cout << "ROM solve no iters used " << no_iter << endl;
			Eigen::VectorXd repro_solve_affine = RB * solve_affine;
			Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
			recover_snapshot_loop(reconstruct_solution, field_x, field_y);
			Eigen::VectorXd  snap_x = Eigen::VectorXd::Zero(cluster_mean_x.num_elements());
			Eigen::VectorXd  snap_y = Eigen::VectorXd::Zero(cluster_mean_x.num_elements());
			for (int index_recr = 0; index_recr < cluster_mean_x.num_elements(); ++index_recr)
			{
				snap_x(index_recr) = field_x[index_recr];
				snap_y(index_recr) = field_y[index_recr];
			}
			cout << "solved snap_x.norm() " << snap_x.norm() << endl;
			cout << "solved snap_y.norm() " << snap_y.norm() << endl;

			if (write_ROM_field || (qoi_dof >= 0))
			{
				locROM_qoi = recover_snapshot_data(reconstruct_solution, 0);
			}
			collected_qoi(fine_grid_dir0_index, fine_grid_dir1_index) = locROM_qoi;
			if (use_fine_grid_VV_and_load_ref)
			{
				collected_relative_L2errors(fine_grid_dir0_index, fine_grid_dir1_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
				collected_relative_Linferrors(fine_grid_dir0_index, fine_grid_dir1_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
				if (use_non_unique_up_to_two)
				{
					collected_relative_L2errors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
					collected_relative_Linferrors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
				}
			}
			fine_grid_dir1_index++;
			if (fine_grid_dir1_index == fine_grid_dir1)
			{
				fine_grid_dir1_index = 0;
				fine_grid_dir0_index++;
			}			
		} // for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		std::stringstream sstm;
		sstm << "VV_ROM_cluster.txt";
		std::string LocROM_txt = sstm.str();
		const char* outname = LocROM_txt.c_str();
		ofstream myfile (outname);
		if (myfile.is_open())
		{
			for (int i0 = 0; i0 < fine_grid_dir0; i0++)
			{
				for (int i1 = 0; i1 < fine_grid_dir1; i1++)
				{
					myfile << std::setprecision(17) << collected_qoi(i0,i1) << "\t";
				}
				myfile << "\n";
			}
			myfile.close();
		}
		else cout << "Unable to open file"; 

		if (use_fine_grid_VV_and_load_ref)
		{

			std::stringstream sstm_VV;
			sstm_VV << "VV_ROM_cluster_VV.txt";
			std::string LocROM_txt_VV = sstm_VV.str();
			const char* outname_VV = LocROM_txt_VV.c_str();
			ofstream myfile_VV (outname_VV);
			if (myfile_VV.is_open())
			{
				for (int i0 = 0; i0 < fine_grid_dir0; i0++)
				{
					for (int i1 = 0; i1 < fine_grid_dir1; i1++)
					{
						myfile_VV << std::setprecision(17) << collected_relative_L2errors(i0,i1) << "\t";
					}
					myfile_VV << "\n";
				}
				myfile_VV.close();
			}
			else cout << "Unable to open file"; 

			std::stringstream sstm_VV_Linf;
			sstm_VV_Linf << "VV_ROM_cluster_VV_Linf.txt";
			std::string LocROM_txt_VV_Linf = sstm_VV_Linf.str();
			const char* outname_VV_Linf = LocROM_txt_VV_Linf.c_str();
			ofstream myfile_VV_Linf (outname_VV_Linf);
			if (myfile_VV_Linf.is_open())
			{
				for (int i0 = 0; i0 < fine_grid_dir0; i0++)
				{
					for (int i1 = 0; i1 < fine_grid_dir1; i1++)
					{
						myfile_VV_Linf << std::setprecision(17) << collected_relative_Linferrors(i0,i1) << "\t";
					}
					myfile_VV_Linf << "\n";
				}
				myfile_VV_Linf.close();
			}
			else cout << "Unable to open file"; 
			
			if (use_non_unique_up_to_two)
			{
				std::stringstream sstm_VV;
				sstm_VV << "VV_ROM_cluster_VV_v2.txt";
				std::string LocROM_txt_VV = sstm_VV.str();
				const char* outname_VV = LocROM_txt_VV.c_str();
				ofstream myfile_VV (outname_VV);
				if (myfile_VV.is_open())
				{
					for (int i0 = 0; i0 < fine_grid_dir0; i0++)
					{
						for (int i1 = 0; i1 < fine_grid_dir1; i1++)
						{
							myfile_VV << std::setprecision(17) << collected_relative_L2errors_v2(i0,i1) << "\t";
						}
						myfile_VV << "\n";
					}
					myfile_VV.close();
				}
				else cout << "Unable to open file"; 

				std::stringstream sstm_VV_Linf;
				sstm_VV_Linf << "VV_ROM_cluster_VV_Linf_v2.txt";
				std::string LocROM_txt_VV_Linf = sstm_VV_Linf.str();
				const char* outname_VV_Linf = LocROM_txt_VV_Linf.c_str();
				ofstream myfile_VV_Linf (outname_VV_Linf);
				if (myfile_VV_Linf.is_open())
				{
					for (int i0 = 0; i0 < fine_grid_dir0; i0++)
					{
						for (int i1 = 0; i1 < fine_grid_dir1; i1++)
						{
							myfile_VV_Linf << std::setprecision(17) << collected_relative_Linferrors_v2(i0,i1) << "\t";
						}
						myfile_VV_Linf << "\n";
					}
					myfile_VV_Linf.close();
				}
				else cout << "Unable to open file"; 
			}

		}

		if (compute_smaller_model_errs)
		{
			// repeat the parameter sweep with decreasing RB sizes up to 1, but in a separate function for readability
			for (int i=0; i < RBsize; ++i)
			{
    			online_snapshot_check_with_smaller_basis_VV(i);
			}

		}

	}



    }

    void CoupledLinearNS_TT::online_snapshot_check_with_smaller_basis_VV(int reduction_int)
	{

		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
		cout << "initiate online_snapshot_check_with_smaller_basis_VV RB=" << RBsize - reduction_int <<  endl;
		cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;


		// Question: how to init?
		// could use all-zero or the cluster-mean
		Array<OneD, NekDouble> cluster_mean_x(snapshot_x_collection[0].num_elements(), 0.0);
		Array<OneD, NekDouble> cluster_mean_y(snapshot_y_collection[0].num_elements(), 0.0);
/*		for (std::set<int>::iterator it=current_cluster.begin(); it!=current_cluster.end(); ++it)
		{
			for (int i = 0; i < snapshot_x_collection[0].num_elements(); ++i)
			{
				cluster_mean_x[i] += (1.0 / current_cluster.size()) * snapshot_x_collection[*it][i];
				cluster_mean_y[i] += (1.0 / current_cluster.size()) * snapshot_y_collection[*it][i];
			}
		} */

		Eigen::MatrixXd collected_qoi = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_L2errors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_Linferrors = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_L2errors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		Eigen::MatrixXd collected_relative_Linferrors_v2 = Eigen::MatrixXd::Zero(fine_grid_dir0, fine_grid_dir1);
		int fine_grid_dir0_index = 0;
		int fine_grid_dir1_index = 0;
		double locROM_qoi;
		for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		{
//			cout << "fine_grid_iter_index " << iter_index << " of max " << fine_grid_dir0*fine_grid_dir1 << endl;
			int current_index = iter_index;
			double current_nu;
			double w;
			Array<OneD, NekDouble> current_param = fine_general_param_vector[current_index];
			w = current_param[0];	
			current_nu = current_param[1];
			if (debug_mode)
			{
//				cout << " VV online phase current nu " << current_nu << endl;
//				cout << " VV online phase current w " << w << endl;
			}
			Set_m_kinvis( current_nu );
			if (use_Newton)
			{
				DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
			}
			Eigen::MatrixXd curr_xy_proj = project_onto_basis(cluster_mean_x, cluster_mean_y);
			Eigen::MatrixXd affine_mat_proj;
			Eigen::VectorXd affine_vec_proj;
			affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
			affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
			Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);

			Eigen::MatrixXd affine_mat_proj_cut = affine_mat_proj.block(0,0,RBsize-reduction_int,RBsize-reduction_int);
			Eigen::VectorXd affine_vec_proj_cut = affine_vec_proj.head(RBsize-reduction_int);
			Eigen::VectorXd solve_affine_cut = affine_mat_proj_cut.colPivHouseholderQr().solve(affine_vec_proj_cut);
			solve_affine = solve_affine_cut;
	//		cout << "solve_affine " << solve_affine << endl;

	//		Eigen::VectorXd repro_solve_affine = RB * solve_affine;
			Eigen::VectorXd repro_solve_affine = RB.leftCols(RBsize-reduction_int) * solve_affine_cut;

			double relative_change_error;
			int no_iter=0;
			Array<OneD, double> field_x;
			Array<OneD, double> field_y;
			// now start looping
			do
			{
				// for now only Oseen // otherwise need to do the DoInitialiseAdv(cluster_mean_x, cluster_mean_y);
				Eigen::VectorXd prev_solve_affine = solve_affine;
				Eigen::VectorXd repro_solve_affine = RB.leftCols(RBsize-reduction_int) * solve_affine;
				Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);

				recover_snapshot_loop(reconstruct_solution, field_x, field_y);
				if (use_Newton)
				{
					DoInitialiseAdv(field_x, field_y);
				}
				curr_xy_proj = project_onto_basis(field_x, field_y);
				if (parameter_space_dimension == 1)
				{
					affine_mat_proj = gen_affine_mat_proj(current_nu);
					affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
				}
				else if (parameter_space_dimension == 2)
				{
					affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
					affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
					affine_mat_proj_cut = affine_mat_proj.block(0,0,RBsize-reduction_int,RBsize-reduction_int);
					affine_vec_proj_cut = affine_vec_proj.head(RBsize-reduction_int);
				}
//				solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
			    solve_affine_cut = affine_mat_proj_cut.colPivHouseholderQr().solve(affine_vec_proj_cut);
                solve_affine = solve_affine_cut;
				relative_change_error = (solve_affine - prev_solve_affine).norm() / prev_solve_affine.norm();
//				cout << "relative_change_error " << relative_change_error << endl;
				no_iter++;
			} 
			while( ((relative_change_error > 1e-12) && (no_iter < 100)) );
//			cout << "ROM solve no iters used " << no_iter << endl;
			repro_solve_affine = RB.leftCols(RBsize-reduction_int)  * solve_affine;
			Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
			if (write_ROM_field || (qoi_dof >= 0))
			{
				locROM_qoi = recover_snapshot_data(reconstruct_solution, 0);
			}
			collected_qoi(fine_grid_dir0_index, fine_grid_dir1_index) = locROM_qoi;
			if (use_fine_grid_VV_and_load_ref)
			{
				collected_relative_L2errors(fine_grid_dir0_index, fine_grid_dir1_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
				collected_relative_Linferrors(fine_grid_dir0_index, fine_grid_dir1_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index], snapshot_y_collection_VV[iter_index]);
				if (use_non_unique_up_to_two)
				{
					collected_relative_L2errors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = L2norm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / L2norm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
					collected_relative_Linferrors_v2(fine_grid_dir0_index, fine_grid_dir1_index) = Linfnorm_abs_error_ITHACA(field_x, field_y, snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]) / Linfnorm_ITHACA(snapshot_x_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1], snapshot_y_collection_VV[iter_index + fine_grid_dir0*fine_grid_dir1]);
				}
			}
			fine_grid_dir1_index++;
			if (fine_grid_dir1_index == fine_grid_dir1)
			{
				fine_grid_dir1_index = 0;
				fine_grid_dir0_index++;
			}			
		} // for (int iter_index = 0; iter_index < fine_grid_dir0*fine_grid_dir1; ++iter_index)
		std::stringstream sstm;
		sstm << "VV_ROM_cluster_reduc" << reduction_int << ".txt";
		std::string LocROM_txt = sstm.str();
		const char* outname = LocROM_txt.c_str();
		ofstream myfile (outname);
		if (myfile.is_open())
		{
			for (int i0 = 0; i0 < fine_grid_dir0; i0++)
			{
				for (int i1 = 0; i1 < fine_grid_dir1; i1++)
				{
					myfile << std::setprecision(17) << collected_qoi(i0,i1) << "\t";
				}
				myfile << "\n";
			}
			myfile.close();
		}
		else cout << "Unable to open file"; 

		if (use_fine_grid_VV_and_load_ref)
		{

			std::stringstream sstm_VV;
			sstm_VV << "VV_ROM_cluster_VV_reduc" << reduction_int << ".txt";
			std::string LocROM_txt_VV = sstm_VV.str();
			const char* outname_VV = LocROM_txt_VV.c_str();
			ofstream myfile_VV (outname_VV);
			if (myfile_VV.is_open())
			{
				for (int i0 = 0; i0 < fine_grid_dir0; i0++)
				{
					for (int i1 = 0; i1 < fine_grid_dir1; i1++)
					{
						myfile_VV << std::setprecision(17) << collected_relative_L2errors(i0,i1) << "\t";
					}
					myfile_VV << "\n";
				}
				myfile_VV.close();
			}
			else cout << "Unable to open file"; 

			std::stringstream sstm_VV_Linf;
			sstm_VV_Linf << "VV_ROM_cluster_VV_Linf_reduc" << reduction_int << ".txt";
			std::string LocROM_txt_VV_Linf = sstm_VV_Linf.str();
			const char* outname_VV_Linf = LocROM_txt_VV_Linf.c_str();
			ofstream myfile_VV_Linf (outname_VV_Linf);
			if (myfile_VV_Linf.is_open())
			{
				for (int i0 = 0; i0 < fine_grid_dir0; i0++)
				{
					for (int i1 = 0; i1 < fine_grid_dir1; i1++)
					{
						myfile_VV_Linf << std::setprecision(17) << collected_relative_Linferrors(i0,i1) << "\t";
					}
					myfile_VV_Linf << "\n";
				}
				myfile_VV_Linf.close();
			}
			else cout << "Unable to open file"; 
			
			if (use_non_unique_up_to_two)
			{
				std::stringstream sstm_VV;
				sstm_VV << "VV_ROM_cluster_VV_v2_reduc" << reduction_int << ".txt";
				std::string LocROM_txt_VV = sstm_VV.str();
				const char* outname_VV = LocROM_txt_VV.c_str();
				ofstream myfile_VV (outname_VV);
				if (myfile_VV.is_open())
				{
					for (int i0 = 0; i0 < fine_grid_dir0; i0++)
					{
						for (int i1 = 0; i1 < fine_grid_dir1; i1++)
						{
							myfile_VV << std::setprecision(17) << collected_relative_L2errors_v2(i0,i1) << "\t";
						}
						myfile_VV << "\n";
					}
					myfile_VV.close();
				}
				else cout << "Unable to open file"; 

				std::stringstream sstm_VV_Linf;
				sstm_VV_Linf << "VV_ROM_cluster_VV_Linf_v2_reduc" << reduction_int << ".txt";
				std::string LocROM_txt_VV_Linf = sstm_VV_Linf.str();
				const char* outname_VV_Linf = LocROM_txt_VV_Linf.c_str();
				ofstream myfile_VV_Linf (outname_VV_Linf);
				if (myfile_VV_Linf.is_open())
				{
					for (int i0 = 0; i0 < fine_grid_dir0; i0++)
					{
						for (int i1 = 0; i1 < fine_grid_dir1; i1++)
						{
							myfile_VV_Linf << std::setprecision(17) << collected_relative_Linferrors_v2(i0,i1) << "\t";
						}
						myfile_VV_Linf << "\n";
					}
					myfile_VV_Linf.close();
				}
				else cout << "Unable to open file"; 
			}

		}

	}

    void CoupledLinearNS_TT::online_snapshot_check_with_smaller_basis(int reduction_int)
	{

	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	cout << "initiate online_snapshot_check_with_smaller_basis RB=" << RBsize - reduction_int <<  endl;
	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;

	Eigen::MatrixXd mat_compare = Eigen::MatrixXd::Zero(f_bnd_dbc_full_size.rows(), 3);  // is of size M_truth_size
	Eigen::VectorXd collected_relative_euclidean_errors = Eigen::VectorXd::Zero(Nmax);
	// start sweeping 
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
		int current_index = iter_index;
		double current_nu;
		double w;
		if (parameter_space_dimension == 1)
		{
			current_nu = param_vector[current_index];
		}
		else if (parameter_space_dimension == 2)
		{
			Array<OneD, NekDouble> current_param = general_param_vector[current_index];
			w = current_param[0];	
			current_nu = current_param[1];
		}
		if (debug_mode)
		{
			cout << " online phase current nu " << current_nu << endl;
			cout << " online phase current w " << w << endl;
		}
		Set_m_kinvis( current_nu );
		if (use_Newton)
		{
			DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		}

		Eigen::MatrixXd curr_xy_proj = project_onto_basis(snapshot_x_collection[current_index], snapshot_y_collection[current_index]);
		Eigen::MatrixXd affine_mat_proj;
		Eigen::VectorXd affine_vec_proj;

		if (parameter_space_dimension == 1)
		{
			affine_mat_proj = gen_affine_mat_proj(current_nu);
			affine_vec_proj = gen_affine_vec_proj(current_nu, current_index);
//			cout << "aff mat 1d " << affine_mat_proj << endl;
//			cout << "aff vec 1d " << affine_vec_proj << endl;
		}
		else if (parameter_space_dimension == 2)
		{
			affine_mat_proj = gen_affine_mat_proj_2d(current_nu, w);
//			cout << "aff mat 2d " << affine_mat_proj << endl;
			affine_vec_proj = gen_affine_vec_proj_2d(current_nu, w, current_index);
//			cout << "aff vec 2d " << affine_vec_proj << endl;
//			Eigen::MatrixXd affine_mat_proj_1d = gen_affine_mat_proj(current_nu);
//			Eigen::VectorXd affine_vec_proj_1d = gen_affine_vec_proj(current_nu, current_index);
		}

		Eigen::MatrixXd affine_mat_proj_cut = affine_mat_proj.block(0,0,RBsize-reduction_int,RBsize-reduction_int);
		Eigen::VectorXd affine_vec_proj_cut = affine_vec_proj.head(RBsize-reduction_int);
		Eigen::VectorXd solve_affine_cut = affine_mat_proj_cut.colPivHouseholderQr().solve(affine_vec_proj_cut);
		Eigen::VectorXd solve_affine = affine_mat_proj.colPivHouseholderQr().solve(affine_vec_proj);
//		cout << "solve_affine " << solve_affine << endl;

//		Eigen::VectorXd repro_solve_affine = RB * solve_affine;
		Eigen::VectorXd repro_solve_affine = RB.leftCols(RBsize-reduction_int) * solve_affine_cut;

		Eigen::VectorXd reconstruct_solution = reconstruct_solution_w_dbc(repro_solve_affine);
		if (globally_connected == 1)
		{
			mat_compare.col(0) = M_collect_f_all.col(current_index);
		}
		else
		{
			mat_compare.col(0) = collect_f_all.col(current_index);
			if (debug_mode)
			{
				Eigen::VectorXd current_f_all = Eigen::VectorXd::Zero(collect_f_all.rows());
				current_f_all = collect_f_all.col(current_index);
				Eigen::VectorXd current_f_all_wo_dbc = remove_rows(current_f_all, elem_loc_dbc);
				Eigen::VectorXd proj_current_f_all_wo_dbc = RB.transpose() * current_f_all_wo_dbc;
//				cout << "proj_current_f_all_wo_dbc " << proj_current_f_all_wo_dbc << endl;
				Eigen::VectorXd correctRHS = affine_mat_proj * proj_current_f_all_wo_dbc;
				Eigen::VectorXd correction_RHS = correctRHS - affine_vec_proj;
//				cout << "correctRHS " << correctRHS << endl;
//				cout << "correction_RHS " << correction_RHS << endl;
			}
		}
		mat_compare.col(1) = reconstruct_solution; // sembra abbastanza bene
		mat_compare.col(2) = mat_compare.col(1) - mat_compare.col(0);
//		cout << mat_compare << endl;
//		cout << "relative euclidean error norm: " << mat_compare.col(2).norm() / mat_compare.col(0).norm() << " of snapshot number " << iter_index << endl;

		// compare this to the projection error onto RB, cannot use collect_f_all, is numerically unstable
/*		cout << "mat_compare.cols() " << mat_compare.cols() << endl;
		cout << "mat_compare.rows() " << mat_compare.rows() << endl;
		cout << "collect_f_all.cols() " << collect_f_all.cols() << endl;
		cout << "collect_f_all.rows() " << collect_f_all.rows() << endl;
		cout << "RB.cols() " << RB.cols() << endl;
		cout << "RB.rows() " << RB.rows() << endl; */

		
//		cout << "curr_xy_proj.cols() " << curr_xy_proj.cols() << endl;
	//	cout << "curr_xy_proj.rows() " << curr_xy_proj.rows() << endl; 


//		Eigen::VectorXd proj_solution = collect_f_all.transpose() * mat_compare.col(0);
//		Eigen::VectorXd reproj_solution = collect_f_all * proj_solution;
		Eigen::VectorXd FOM_solution = mat_compare.col(0);
		Eigen::VectorXd FOM_solution_wo_dbc = remove_rows(FOM_solution, elem_loc_dbc);
		Eigen::VectorXd proj_FOM_solution_wo_dbc = RB.transpose() * FOM_solution_wo_dbc;
		Eigen::VectorXd reproj_FOM_solution_wo_dbc = RB * proj_FOM_solution_wo_dbc;
		Eigen::VectorXd reconstruct_FOM_solution = reconstruct_solution_w_dbc(reproj_FOM_solution_wo_dbc);

//		cout << "reconstruct_FOM_solution.norm(): " << reconstruct_FOM_solution.norm() << endl;
//		cout << "reconstruct_solution.norm(): " << reconstruct_solution.norm() << endl;
		Eigen::VectorXd diff_projection = reconstruct_FOM_solution - mat_compare.col(0);

//		cout << "relative euclidean projection error norm: " << diff_projection.norm() / mat_compare.col(0).norm() << " of snapshot number " << iter_index << endl;


		// now only in RB:
		Eigen::VectorXd diff_projection_RB = reproj_FOM_solution_wo_dbc - remove_rows(mat_compare.col(0), elem_loc_dbc);
		Eigen::VectorXd diff_RB = repro_solve_affine - remove_rows(mat_compare.col(0), elem_loc_dbc);
//		cout << "relative euclidean RB projection error norm: " << diff_projection_RB.norm() / remove_rows(mat_compare.col(0), elem_loc_dbc).norm() << " of snapshot number " << iter_index << endl;
//		cout << "relative euclidean RB error norm: " << diff_RB.norm() / remove_rows(mat_compare.col(0), elem_loc_dbc).norm() << " of snapshot number " << iter_index << endl;

		// have to use curr_xy_proj for better approximations


		if (debug_mode)
		{
			cout << "snapshot_x_collection.num_elements() " << snapshot_x_collection.num_elements() << " snapshot_x_collection[0].num_elements() " << snapshot_x_collection[0].num_elements() << endl;
		}

		if (write_ROM_field || (qoi_dof >= 0))
		{
			recover_snapshot_data(reconstruct_solution, current_index); // this is setting the fields and fieldcoeffs
		}

	  	Eigen::VectorXd f_bnd = reconstruct_solution.head(curr_f_bnd.size());
		Eigen::VectorXd f_int = reconstruct_solution.tail(curr_f_int.size());
		Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
		Array<OneD, unsigned int> bmap, imap; 
		Array<OneD, double> field_0(GetNcoeffs());
		Array<OneD, double> field_1(GetNcoeffs());
		Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
		Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
		int cnt = 0;
		int cnt1 = 0;
		int nvel = 2;
		int nz_loc = 1;
		int  nplanecoeffs = fields[0]->GetNcoeffs();
		int  nel  = m_fields[0]->GetNumElmts();
		for(int i = 0; i < nel; ++i) 
		{
		      int eid  = i;
		      fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
		      fields[0]->GetExp(eid)->GetInteriorMap(imap);
		      int nbnd   = bmap.num_elements();
		      int nint   = imap.num_elements();
			      int offset = fields[0]->GetCoeff_Offset(eid);
		            
		      for(int j = 0; j < nvel; ++j)
		      {
		           for(int n = 0; n < nz_loc; ++n)
		           {
		                    for(int k = 0; k < nbnd; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
		                    }
		                    
		                    for(int k = 0; k < nint; ++k)
		                    {
		                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
		                    }
		                    cnt  += nbnd;
		                    cnt1 += nint;
		           }
		      }
		}
		Array<OneD, double> test_nn = fields[0]->GetCoeffs();
		fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
		fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);

        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.num_elements()+1);
        int i;
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
	    }
//		cout << "after recover_snapshot_data have fieldcoeffs[0].num_elements() " << fieldcoeffs[0].num_elements() << endl;

//		cout << "after recover_snapshot_data have curr_PhysBaseVec_x.num_elements() " << curr_PhysBaseVec_x.num_elements() << endl;
//		cout << "after recover_snapshot_data have snapshot_x_collection[current_index].num_elements() " << snapshot_x_collection[current_index].num_elements() << endl;

		// have the FOM_snapshot_solution projection available as curr_xy_proj
		// datastructure: 	Array<OneD, Array<OneD, NekDouble> > snapshot_x_collection;
		Eigen::VectorXd diff_x_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		Eigen::VectorXd snap_x = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		Eigen::VectorXd diff_y_RB_solve = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		Eigen::VectorXd snap_y = Eigen::VectorXd::Zero(curr_PhysBaseVec_x.num_elements());
		for (int index_recr = 0; index_recr < curr_PhysBaseVec_x.num_elements(); ++index_recr)
		{
			snap_x(index_recr) = snapshot_x_collection[current_index][index_recr];
			snap_y(index_recr) = snapshot_y_collection[current_index][index_recr];
			diff_x_RB_solve(index_recr) = curr_PhysBaseVec_x[index_recr] - snapshot_x_collection[current_index][index_recr];
			diff_y_RB_solve(index_recr) = curr_PhysBaseVec_y[index_recr] - snapshot_y_collection[current_index][index_recr];
		}

		Eigen::MatrixXd curr_xy_reproj = reproject_from_basis(curr_xy_proj);

		double rel_x = diff_x_RB_solve.norm() / snap_x.norm();
		double rel_y = diff_y_RB_solve.norm() / snap_y.norm();

		cout << "relative euclidean error norm in x coords: " << rel_x << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean error norm in y coords: " << rel_y << " of snapshot number " << iter_index << endl;

		collected_relative_euclidean_errors(iter_index) = sqrt( rel_x * rel_x + rel_y * rel_y );

//		cout << "curr_xy_reproj.cols() " << curr_xy_reproj.cols() << endl;
//		cout << "curr_xy_reproj.rows() " << curr_xy_reproj.rows() << endl;
//		cout << "snap_x.rows() " << snap_x.rows() << endl;

		Eigen::VectorXd diff_x_proj = curr_xy_reproj.col(0) - snap_x;
		Eigen::VectorXd diff_y_proj = curr_xy_reproj.col(1) - snap_y;
		
		cout << "relative euclidean projection error norm in x coords: " << diff_x_proj.norm() / snap_x.norm() << " of snapshot number " << iter_index << endl;
		cout << "relative euclidean projection error norm in y coords: " << diff_y_proj.norm() / snap_y.norm() << " of snapshot number " << iter_index << endl;

	} //for (int iter_index = 0; iter_index < Nmax; ++iter_index)

	std::stringstream sstm;
	sstm << "ROM_cluster_reduc" << reduction_int << ".txt";
	std::string ROM_txt = sstm.str();
	const char* outname = ROM_txt.c_str();
	ofstream myfile (outname);
	if (myfile.is_open())
	{
		for (int iter_index = 0; iter_index < Nmax; ++iter_index)
		{
			myfile << std::setprecision(17) << collected_relative_euclidean_errors(iter_index) << "\n";
		}
		myfile.close();
	}
	else cout << "Unable to open file"; 

	}

    void CoupledLinearNS_TT::recover_snapshot_loop(Eigen::VectorXd reconstruct_solution, Array<OneD, double> & field_x, Array<OneD, double> & field_y)
    {
	Eigen::VectorXd f_bnd = reconstruct_solution.head(curr_f_bnd.size());
	Eigen::VectorXd f_int = reconstruct_solution.tail(curr_f_int.size());
	Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
	Array<OneD, unsigned int> bmap, imap; 
	Array<OneD, double> field_0(GetNcoeffs());
	Array<OneD, double> field_1(GetNcoeffs());
	Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
	Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
	int cnt = 0;
	int cnt1 = 0;
	int nvel = 2;
	int nz_loc = 1;
	int  nplanecoeffs = fields[0]->GetNcoeffs();
	int  nel  = m_fields[0]->GetNumElmts();
	for(int i = 0; i < nel; ++i) 
	{
	      int eid  = i;
	      fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
	      fields[0]->GetExp(eid)->GetInteriorMap(imap);
	      int nbnd   = bmap.num_elements();
	      int nint   = imap.num_elements();
	      int offset = fields[0]->GetCoeff_Offset(eid);
	            
	      for(int j = 0; j < nvel; ++j)
	      {
	           for(int n = 0; n < nz_loc; ++n)
	           {
	                    for(int k = 0; k < nbnd; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
	                    }
	                    
	                    for(int k = 0; k < nint; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
	                    }
	                    cnt  += nbnd;
	                    cnt1 += nint;
	           }
	      }
	}
	Array<OneD, double> test_nn = fields[0]->GetCoeffs();
	fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
	fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);
	field_x = curr_PhysBaseVec_x;
	field_y = curr_PhysBaseVec_y;
    }

    double CoupledLinearNS_TT::recover_snapshot_data(Eigen::VectorXd reconstruct_solution, int current_index)
    {

	Eigen::VectorXd f_bnd = reconstruct_solution.head(curr_f_bnd.size());
	Eigen::VectorXd f_int = reconstruct_solution.tail(curr_f_int.size());
	Array<OneD, MultiRegions::ExpListSharedPtr> fields = UpdateFields(); 
	Array<OneD, unsigned int> bmap, imap; 
	Array<OneD, double> field_0(GetNcoeffs());
	Array<OneD, double> field_1(GetNcoeffs());
	Array<OneD, double> curr_PhysBaseVec_x(GetNpoints(), 0.0);
	Array<OneD, double> curr_PhysBaseVec_y(GetNpoints(), 0.0);
	int cnt = 0;
	int cnt1 = 0;
	int nvel = 2;
	int nz_loc = 1;
	int  nplanecoeffs = fields[0]->GetNcoeffs();
	int  nel  = m_fields[0]->GetNumElmts();
	for(int i = 0; i < nel; ++i) 
	{
	      int eid  = i;
	      fields[0]->GetExp(eid)->GetBoundaryMap(bmap);
	      fields[0]->GetExp(eid)->GetInteriorMap(imap);
	      int nbnd   = bmap.num_elements();
	      int nint   = imap.num_elements();
	      int offset = fields[0]->GetCoeff_Offset(eid);
	            
	      for(int j = 0; j < nvel; ++j)
	      {
	           for(int n = 0; n < nz_loc; ++n)
	           {
	                    for(int k = 0; k < nbnd; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+bmap[k], f_bnd(cnt+k));
	                    }
	                    
	                    for(int k = 0; k < nint; ++k)
	                    {
	                        fields[j]->SetCoeff(n*nplanecoeffs + offset+imap[k], f_int(cnt1+k));
	                    }
	                    cnt  += nbnd;
	                    cnt1 += nint;
	           }
	      }
	}
	Array<OneD, double> test_nn = fields[0]->GetCoeffs();
	fields[0]->BwdTrans_IterPerExp(fields[0]->GetCoeffs(), curr_PhysBaseVec_x);
	fields[1]->BwdTrans_IterPerExp(fields[1]->GetCoeffs(), curr_PhysBaseVec_y);

	Eigen::VectorXd eigen_phys_basis_x = Eigen::VectorXd::Zero(GetNpoints());
	Eigen::VectorXd eigen_phys_basis_y = Eigen::VectorXd::Zero(GetNpoints());
	Eigen::VectorXd eigen_phys_basis_x_snap = Eigen::VectorXd::Zero(GetNpoints());
	Eigen::VectorXd eigen_phys_basis_y_snap = Eigen::VectorXd::Zero(GetNpoints());

	for (int index_phys_base=0; index_phys_base<GetNpoints(); index_phys_base++)
	{
		eigen_phys_basis_x(index_phys_base) = curr_PhysBaseVec_x[index_phys_base];
		eigen_phys_basis_y(index_phys_base) = curr_PhysBaseVec_y[index_phys_base];
		
		eigen_phys_basis_x_snap(index_phys_base) = snapshot_x_collection[current_index][index_phys_base];
		eigen_phys_basis_y_snap(index_phys_base) = snapshot_y_collection[current_index][index_phys_base];

	}

	if (debug_mode)
	{
	/*	cout << "eigen_phys_basis_x.norm() " << eigen_phys_basis_x.norm() << endl;
		cout << "eigen_phys_basis_x_snap.norm() " << eigen_phys_basis_x_snap.norm() << endl;
		cout << "(eigen_phys_basis_x - eigen_phys_basis_x_snap).norm() " << (eigen_phys_basis_x - eigen_phys_basis_x_snap).norm() << endl;
		cout << "eigen_phys_basis_y.norm() " << eigen_phys_basis_y.norm() << endl;
		cout << "eigen_phys_basis_y_snap.norm() " << eigen_phys_basis_y_snap.norm() << endl;
		cout << "(eigen_phys_basis_y - eigen_phys_basis_y_snap).norm() " << (eigen_phys_basis_y - eigen_phys_basis_y_snap).norm() << endl; */
	}




        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.num_elements()+1);
        std::vector<std::string> variables(m_fields.num_elements()+1);
        int i;
        
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
	    if (debug_mode)
	    {
		//    cout << "variables[i] " << variables[i] << endl;
	    }
        }

	if (debug_mode)
	{
	//	cout << "m_singleMode " << m_singleMode << endl;	
	}

	fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs(), 0.0);  
        variables[i] = "p"; 

	if (write_ROM_field)
	{
		WriteFld("Test.fld",m_fields[0],fieldcoeffs,variables);
	}

	if (qoi_dof >= 0)
	{
		//cout << "converged qoi dof ROM " << curr_PhysBaseVec_y[qoi_dof] << endl;
	}
	return curr_PhysBaseVec_y[qoi_dof];

/*        fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());  
        // project pressure field to velocity space        
        if(m_singleMode==true)
        {
            Array<OneD, NekDouble > tmpfieldcoeffs (m_fields[0]->GetNcoeffs()/2);
            m_pressure->GetPlane(0)->BwdTrans_IterPerExp(m_pressure->GetPlane(0)->GetCoeffs(), m_pressure->GetPlane(0)->UpdatePhys());
            m_pressure->GetPlane(1)->BwdTrans_IterPerExp(m_pressure->GetPlane(1)->GetCoeffs(), m_pressure->GetPlane(1)->UpdatePhys()); 
            m_fields[0]->GetPlane(0)->FwdTrans_IterPerExp(m_pressure->GetPlane(0)->GetPhys(),fieldcoeffs[i]);
            m_fields[0]->GetPlane(1)->FwdTrans_IterPerExp(m_pressure->GetPlane(1)->GetPhys(),tmpfieldcoeffs);
            for(int e=0; e<m_fields[0]->GetNcoeffs()/2; e++)
            {
                fieldcoeffs[i][e+m_fields[0]->GetNcoeffs()/2] = tmpfieldcoeffs[e];
            }          
        }
        else
        {
            m_pressure->BwdTrans_IterPerExp(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());
            m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),fieldcoeffs[i]);
        }
        variables[i] = "p"; 
        
        std::string outname = m_sessionName + ".fld";
        
        WriteFld(outname,m_fields[0],fieldcoeffs,variables);
*/




    }
	
    void CoupledLinearNS_TT::offline_phase()
    {
	time_t timer_1;
	time_t timer_2;
//	  struct tm y2k = {0};
//	  double seconds;
//	  y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
//	  y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
//	  time(&timer);  /* get current time; same as: timer = time(NULL)  */
//	  seconds = difftime(timer,mktime(&y2k));
//	  printf ("%.f seconds since January 1, 2000 in the current timezone", seconds);

	int load_snapshot_data_from_files = m_session->GetParameter("load_snapshot_data_from_files");
	int number_of_snapshots = m_session->GetParameter("number_of_snapshots");
	if (m_session->DefinesParameter("parameter_space_dimension")) 
	{
		parameter_space_dimension = m_session->GetParameter("parameter_space_dimension");	
	}
	else
	{
		parameter_space_dimension = 1;
	}
	if (m_session->DefinesParameter("load_cO_snapshot_data_from_files")) 
	{
		load_cO_snapshot_data_from_files = m_session->GetParameter("load_cO_snapshot_data_from_files");	
	}
	else
	{
		load_cO_snapshot_data_from_files = 0;
	}
	if (m_session->DefinesParameter("do_trafo_check")) 
	{
		do_trafo_check = m_session->GetParameter("do_trafo_check");	
	}
	else
	{
		do_trafo_check = 1;
	}
	if (m_session->DefinesParameter("compute_smaller_model_errs")) 
	{
		compute_smaller_model_errs = m_session->GetParameter("compute_smaller_model_errs");	
	}
	else
	{
		compute_smaller_model_errs = 0;
	}
	if (m_session->DefinesParameter("qoi_dof")) 
	{
		qoi_dof = m_session->GetParameter("qoi_dof");	
	}
	else
	{
		qoi_dof = -1;
	}
	POD_tolerance = m_session->GetParameter("POD_tolerance");
	ref_param_index = m_session->GetParameter("ref_param_index");
	ref_param_nu = m_session->GetParameter("ref_param_nu");
	if (m_session->DefinesParameter("globally_connected")) // this sets how the truth system global coupling is enforced
	{
		globally_connected = m_session->GetParameter("globally_connected");
	}
	else
	{
		globally_connected = 0;
	}
	if (m_session->DefinesParameter("use_Newton")) // set if Newton or Oseen iteration
	{
		use_Newton = m_session->GetParameter("use_Newton");
	}
	else
	{
		use_Newton = 0;
	}
	if (m_session->DefinesParameter("use_ANN")) // after POD output data for ANN construction and then break
	{
		use_ANN = m_session->GetParameter("use_ANN");
	}
	else
	{
		use_ANN = 0;
	}
	if (m_session->DefinesParameter("use_ANN_local")) // after POD output data for ANN construction in local ROM and then break
	{
		use_ANN_local = m_session->GetParameter("use_ANN_local");
	}
	else
	{
		use_ANN_local = 0;
	}
	if (m_session->DefinesParameter("use_non_unique_up_to_two")) // set if Newton or Oseen iteration
	{
		use_non_unique_up_to_two = m_session->GetParameter("use_non_unique_up_to_two");
		number_of_snapshots = 2*number_of_snapshots;
	}
	else
	{
		use_non_unique_up_to_two = 0;
	}
	if (m_session->DefinesParameter("debug_mode")) // debug_mode with many extra information but very slow
	{
		debug_mode = m_session->GetParameter("debug_mode");
	}
	else
	{
		debug_mode = 0;
	} 
	if (m_session->DefinesParameter("write_ROM_field")) 
	{
		write_ROM_field = m_session->GetParameter("write_ROM_field");
	}
	else
	{
		write_ROM_field = 0;
	} 
	if (m_session->DefinesParameter("snapshot_computation_plot_rel_errors")) 
	{
		snapshot_computation_plot_rel_errors = m_session->GetParameter("snapshot_computation_plot_rel_errors");
	}
	else
	{
		snapshot_computation_plot_rel_errors = 0;
	} 
	Nmax = number_of_snapshots;
	if (parameter_space_dimension == 1)
	{
		param_vector = Array<OneD, NekDouble> (Nmax);
	        for(int i = 0; i < number_of_snapshots; ++i)
	        {
			// generate the correct string
			std::stringstream sstm;
			sstm << "param" << i;
			std::string result = sstm.str();
			// const char* rr = result.c_str();
		        param_vector[i] = m_session->GetParameter(result);
	        }
	}
	else if (parameter_space_dimension == 2)
	{
		general_param_vector = Array<OneD, Array<OneD, NekDouble> > (Nmax);
		number_of_snapshots_dir0 = m_session->GetParameter("number_of_snapshots_dir0");
		number_of_snapshots_dir1 = m_session->GetParameter("number_of_snapshots_dir1");
		if (m_session->DefinesParameter("use_fine_grid_VV")) 
		{
			use_fine_grid_VV = m_session->GetParameter("use_fine_grid_VV");
		}
		else
		{
			use_fine_grid_VV = 0;
		} 
		if (m_session->DefinesParameter("use_fine_grid_VV_and_load_ref")) 
		{
			use_fine_grid_VV_and_load_ref = m_session->GetParameter("use_fine_grid_VV_and_load_ref");
		}
		else
		{
			use_fine_grid_VV_and_load_ref = 0;
		} 
		if (m_session->DefinesParameter("use_fine_grid_VV_random")) 
		{
			use_fine_grid_VV_random = m_session->GetParameter("use_fine_grid_VV_random");
		}
		else
		{
			use_fine_grid_VV_random = 0;
		} 
		if (m_session->DefinesParameter("replace_snapshot_with_transformed")) 
		{
			replace_snapshot_with_transformed = m_session->GetParameter("replace_snapshot_with_transformed");
		}
		else
		{
			replace_snapshot_with_transformed = 1;
		} 
		if (m_session->DefinesParameter("fine_grid_dir0")) 
		{
			fine_grid_dir0 = m_session->GetParameter("fine_grid_dir0");
		}
		else
		{
			fine_grid_dir0 = 0;
		} 
		if (m_session->DefinesParameter("fine_grid_dir1")) 
		{
			fine_grid_dir1 = m_session->GetParameter("fine_grid_dir1");
		}
		else
		{
			fine_grid_dir1 = 0;
		} 
		if (m_session->DefinesParameter("use_sparse_poly")) 
		{
			use_sparse_poly = m_session->GetParameter("use_sparse_poly");
		}
		else
		{
			use_sparse_poly = 0;
		} 
		if (m_session->DefinesParameter("max_sparse_poly_approx_dimension")) 
		{
			max_sparse_poly_approx_dimension = m_session->GetParameter("max_sparse_poly_approx_dimension");
		}
		else
		{
			max_sparse_poly_approx_dimension = 1;
		} 
		Array<OneD, NekDouble> index_vector(parameter_space_dimension, 0.0);

//		for(int i = 0; i < parameter_space_dimension; ++i)
//		{
//			cout << "psdiv " << index_vector[i] << endl;
//		}
		int i_all = 0;
//		general_param_vector[i_all] = Array<OneD, NekDouble> (parameter_space_dimension);
		Array<OneD, NekDouble> parameter_point(parameter_space_dimension, 0.0);
		for(int i = 0; i < Nmax; ++i)
		{
			parameter_point = Array<OneD, NekDouble> (parameter_space_dimension, 0.0);
			general_param_vector[i] = parameter_point;
		}
		for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
		{
			// generate the correct string
			std::stringstream sstm;
			sstm << "param" << i0 << "_dir0";
			std::string result = sstm.str();
		        if (i0 == 0)
				start_param_dir0 = m_session->GetParameter(result);
		        if (i0 == number_of_snapshots_dir0-1)
				end_param_dir0 = m_session->GetParameter(result);
			for(int i1 = 0; i1 < number_of_snapshots_dir1; ++i1)
			{
				// generate the correct string
				std::stringstream sstm1;
				sstm1 << "param" << i1 << "_dir1";
				std::string result1 = sstm1.str();
			    if (i1 == 0)
					start_param_dir1 = m_session->GetParameter(result1);
			    if (i1 == number_of_snapshots_dir1-1)
					end_param_dir1 = m_session->GetParameter(result1);
				general_param_vector[i_all][0] = m_session->GetParameter(result);
			    general_param_vector[i_all][1] = m_session->GetParameter(result1);
				i_all = i_all + 1;
//				general_param_vector[i_all] = Array<OneD, NekDouble> (parameter_space_dimension);
			}
		}
		// if there is set use_non_unique_up_to_two then double the param vector
		if (use_non_unique_up_to_two)
		{
			// Nmax should already be correct
			// cout << "Nmax " << Nmax << endl;
			for (int i = 0; i < Nmax/2; ++i)
			{
				general_param_vector[i + Nmax/2][0] = general_param_vector[i][0];
				general_param_vector[i + Nmax/2][1] = general_param_vector[i][1];
			}

		}

		// output sample grid as *.txt
	        std::string sample_grid_txt = "sample_grid.txt";
		const char* outname_sample_grid_txt = sample_grid_txt.c_str();
		ofstream myfile_sample_grid_txt (outname_sample_grid_txt);
		i_all = 0;
		if (myfile_sample_grid_txt.is_open())
		{
			for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
			{
				for(int i1 = 0; i1 < number_of_snapshots_dir1; ++i1)
				{
					myfile_sample_grid_txt << std::setprecision(17) << general_param_vector[i_all][0] << "\t" << general_param_vector[i_all][1] << endl;
					i_all++;
				}
			}
			myfile_sample_grid_txt.close();
		}
		else cout << "Unable to open file"; 

		// output sample grids as *.txt
	        std::string sample_grid_d1_txt = "sample_grid_d1.txt";
		const char* outname_sample_grid_d1_txt = sample_grid_d1_txt.c_str();
		ofstream myfile_sample_grid_d1_txt (outname_sample_grid_d1_txt);
		if (myfile_sample_grid_d1_txt.is_open())
		{
			for(int i0 = 0; i0 < number_of_snapshots_dir0; ++i0)
			{
				myfile_sample_grid_d1_txt << std::setprecision(17) << general_param_vector[i0*number_of_snapshots_dir1][0] << endl;
			}
			myfile_sample_grid_d1_txt.close();
		}
		else cout << "Unable to open file"; 

	        std::string sample_grid_d2_txt = "sample_grid_d2.txt";
		const char* outname_sample_grid_d2_txt = sample_grid_d2_txt.c_str();
		ofstream myfile_sample_grid_d2_txt (outname_sample_grid_d2_txt);
		if (myfile_sample_grid_d2_txt.is_open())
		{
			for(int i1 = 0; i1 < number_of_snapshots_dir1; ++i1)
			{
				myfile_sample_grid_d2_txt << std::setprecision(17) << general_param_vector[i1][1] << endl;
			}
			myfile_sample_grid_d2_txt.close();
		}
		else cout << "Unable to open file"; 

		if (use_fine_grid_VV)
		{
		//	cout << "start_param_dir0 " << start_param_dir0 << endl;
		//	cout << "end_param_dir0 " << end_param_dir0 << endl;
		//	cout << "start_param_dir1 " << start_param_dir1 << endl;
		//	cout << "end_param_dir1 " << end_param_dir1 << endl;
			fine_general_param_vector = Array<OneD, Array<OneD, NekDouble> >(fine_grid_dir0*fine_grid_dir1);
			for(int i = 0; i < fine_grid_dir0*fine_grid_dir1; ++i)
			{
				parameter_point = Array<OneD, NekDouble> (parameter_space_dimension, 0.0);
				fine_general_param_vector[i] = parameter_point;
			}
			int current_index = 0;
			for (int i1 = 0; i1 < fine_grid_dir0; i1++)
			{
				for (int i2 = 0; i2 < fine_grid_dir1; i2++)
				{
					double p0 = start_param_dir0 + double(i1)/(fine_grid_dir0-1) * (end_param_dir0 - start_param_dir0);
					double p1 = start_param_dir1 + double(i2)/(fine_grid_dir1-1) * (end_param_dir1 - start_param_dir1);
					fine_general_param_vector[current_index][0] = p0;
					fine_general_param_vector[current_index][1] = p1;
					current_index++;
//					cout << "p0 " << p0 << " p1 " << p1 << endl;
				}
			}

			// or alternatively overwrite with given random vector
			if (use_fine_grid_VV_random)
			{
				int current_index = 0;
				for (int i1 = 0; i1 < fine_grid_dir0; i1++)
				{
					// generate the correct string
					std::stringstream sstm;
					sstm << "VV_param" << i1 << "_dir0";
					std::string result = sstm.str();
					double param_dir0 = m_session->GetParameter(result);
		
					for (int i2 = 0; i2 < fine_grid_dir1; i2++)
					{

						// generate the correct string
						std::stringstream sstm1;
						sstm1 << "VV_param" << i2 << "_dir1";
						std::string result1 = sstm1.str();
						double param_dir1 = m_session->GetParameter(result1);


						fine_general_param_vector[current_index][0] = param_dir0;
						fine_general_param_vector[current_index][1] = param_dir1;
						current_index++;
 					//	cout << "p0 " << param_dir0 << " p1 " << param_dir1 << endl;
					}
				}
				
			}

			// output fine sample grid as *.txt
		        std::string fine_sample_grid_txt = "fine_sample_grid.txt";
			const char* outname_fine_sample_grid_txt = fine_sample_grid_txt.c_str();
			ofstream myfile_fine_sample_grid_txt (outname_fine_sample_grid_txt);
			i_all = 0;
			if (myfile_fine_sample_grid_txt.is_open())
			{
				for(int i0 = 0; i0 < fine_grid_dir0; ++i0)
				{
					for(int i1 = 0; i1 < fine_grid_dir1; ++i1)
					{
						myfile_fine_sample_grid_txt << std::setprecision(17) << fine_general_param_vector[i_all][0] << "\t" << fine_general_param_vector[i_all][1] << endl;
						i_all++;
					}
				}
				myfile_fine_sample_grid_txt.close();
			}
			else cout << "Unable to open file"; 

			// write to file the VV grid
//				        std::string outname_txt = m_sessionName + ".txt";
		        std::string VV_grid_txt = "VV_grid_d1.txt";
			const char* outname_t = VV_grid_txt.c_str();
			ofstream myfile_t (outname_t);
			if (myfile_t.is_open())
			{
				for (int i1 = 0; i1 < fine_grid_dir0; i1++)
				{
					double p0 = start_param_dir0 + double(i1)/(fine_grid_dir0-1) * (end_param_dir0 - start_param_dir0);
					myfile_t << std::setprecision(17) << p0 << "\t";
				}
       				myfile_t.close();
			}
			else cout << "Unable to open file"; 
		        VV_grid_txt = "VV_grid_d2.txt";
			const char* outname_t2 = VV_grid_txt.c_str();
			ofstream myfile_t2 (outname_t2);
			if (myfile_t2.is_open())
			{
				for (int i2 = 0; i2 < fine_grid_dir1; i2++)
				{
					double p1 = start_param_dir1 + double(i2)/(fine_grid_dir1-1) * (end_param_dir1 - start_param_dir1);
					myfile_t2 << std::setprecision(17) << p1 << "\t";
				}
       				myfile_t2.close();
			}
			else cout << "Unable to open file"; 




		}
		int type_para1 = m_session->GetParameter("type_para1");
//		cout << "type para1 " << type_para1 << endl;
		int type_para2 = m_session->GetParameter("type_para2");
		//cout << "type para2 " << type_para2 << endl;
		number_elem_trafo = m_session->GetParameter("number_elem_trafo");

		elements_trafo = Array<OneD, std::set<int> > (number_elem_trafo);
//	    <P> elem_1 = {32,30,11,10,9,8,17,16,15,14,25,24,23,22}   	</P>   
//          <P> elem_2 = {12, 28, 34, 13,21,20,18,19,26,27}   		</P>   
//	    <P> elem_3 = {29, 31}   	</P>   
//	    <P> elem_4 = {33, 35}   	</P>   
//	    <P> elem_5 = {0,1,2,3,4,5,6,7}  
		elements_trafo[0].insert(32); // of course, this should work automatically
		elements_trafo[0].insert(30);
		elements_trafo[0].insert(11);
		elements_trafo[0].insert(10);
		elements_trafo[0].insert(9);
		elements_trafo[0].insert(8);
		elements_trafo[0].insert(17);
		elements_trafo[0].insert(16);
		elements_trafo[0].insert(15);
		elements_trafo[0].insert(14);
		elements_trafo[0].insert(25);
		elements_trafo[0].insert(24);
		elements_trafo[0].insert(23);
		elements_trafo[0].insert(22);

		elements_trafo[1].insert(12);
		elements_trafo[1].insert(28);
		elements_trafo[1].insert(34);
		elements_trafo[1].insert(13);
		elements_trafo[1].insert(21);
		elements_trafo[1].insert(20);
		elements_trafo[1].insert(18);
		elements_trafo[1].insert(19);
		elements_trafo[1].insert(26);
		elements_trafo[1].insert(27);

		elements_trafo[2].insert(29);
		elements_trafo[2].insert(31);

		elements_trafo[3].insert(33);
		elements_trafo[3].insert(35);

		elements_trafo[4].insert(0);
		elements_trafo[4].insert(1);
		elements_trafo[4].insert(2);
		elements_trafo[4].insert(3);
		elements_trafo[4].insert(4);
		elements_trafo[4].insert(5);
		elements_trafo[4].insert(6);
		elements_trafo[4].insert(7);

/*	def Geo_T(w, elemT, index): # index 0: det, index 1,2,3,4: mat_entries
	if elemT == 0:
		T = np.matrix([[1, 0], [0, 1/w]])
	if elemT == 1:
		T = np.matrix([[1, 0], [0, 2/(3-w)]])
	if elemT == 2:
		T = np.matrix([[1, 0], [-(1-w)/2, 1]])
	if elemT == 3:
		T = np.matrix([[1, 0], [-(w-1)/2, 1]])
	if elemT == 4:
		T = np.matrix([[1, 0], [0, 1]])			
	if index == 0:
		return 1/np.linalg.det(T)
	if index == 1:
		return T[0,0]
	if index == 2:
		return T[0,1]
	if index == 3:
		return T[1,0]
	if index == 4:
		return T[1,1]	*/

//		elements_trafo_matrix = Array<OneD, Eigen::Matrix2d > (number_elem_trafo); // but how to do this with symbolic computation?
		 

	}
//	cout << "parameter vector generated" << endl;
//	cout << "general_param_vector.num_elements() " << general_param_vector.num_elements() << endl;
//	cout << "general_param_vector[0].num_elements() " << general_param_vector[0].num_elements() << endl;
//	cout << "general_param_vector[1].num_elements() " << general_param_vector[1].num_elements() << endl;
//	cout << "general_param_vector[2].num_elements() " << general_param_vector[2].num_elements() << endl;
//	cout << "general_param_vector[19].num_elements() " << general_param_vector[19].num_elements() << endl;

/*	for(int i0 = 0; i0 < general_param_vector.num_elements(); ++i0)
	{
		for(int i1 = 0; i1 < general_param_vector[i0].num_elements(); ++i1)
		{
			cout << "general_param_vector[i0][i1] " << general_param_vector[i0][i1] << endl;
		}
	}*/

//	cout << Geo_T( 0.2 , 0, 0) << endl;;
	InitObject();
	if ( load_snapshot_data_from_files )
	{
		if (parameter_space_dimension == 1)
		{
			load_snapshots(number_of_snapshots);
		}
		else if (parameter_space_dimension == 2)
		{
			if (!load_cO_snapshot_data_from_files)
			{
				load_snapshots_geometry_params(number_of_snapshots); 
			}
			else 
			{
				load_snapshots_geometry_params_conv_Oseen(number_of_snapshots); 
			}
			if (use_fine_grid_VV && use_fine_grid_VV_and_load_ref)
			{
				// load the refs -- maybe 	InitObject();   has to happen first
				Array<OneD, NekDouble> collected_fine_grid_ref_qoi = Array<OneD, NekDouble> (fine_grid_dir0*fine_grid_dir1);
				if (use_non_unique_up_to_two)
				{
					collected_fine_grid_ref_qoi = Array<OneD, NekDouble> (2*fine_grid_dir0*fine_grid_dir1);
				}
				
				int nvelo = 2;
		        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot(nvelo); // for a 2D problem

		        for(int i = 0; i < nvelo; ++i)
		        {
		            test_load_snapshot[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);  // number of phys points
		        }

				cout << "no of phys points " << GetNpoints() << endl;

		        std::vector<std::string> fieldStr;
		        for(int i = 0; i < nvelo; ++i)
		        {
		           fieldStr.push_back(m_session->GetVariable(i));
		        }

				snapshot_x_collection_VV = Array<OneD, Array<OneD, NekDouble> > (fine_grid_dir0*fine_grid_dir1);
				snapshot_y_collection_VV = Array<OneD, Array<OneD, NekDouble> > (fine_grid_dir0*fine_grid_dir1);
				if (use_non_unique_up_to_two)
				{
					snapshot_x_collection_VV = Array<OneD, Array<OneD, NekDouble> > (2*fine_grid_dir0*fine_grid_dir1);
					snapshot_y_collection_VV = Array<OneD, Array<OneD, NekDouble> > (2*fine_grid_dir0*fine_grid_dir1);
				}

		        for(int i = 0; i < fine_grid_dir0*fine_grid_dir1; ++i)
		        {
					// generate the correct string
					std::stringstream sstm;
					sstm << "VV" << i+1;
					std::string result = sstm.str();
					const char* rr = result.c_str();

				        EvaluateFunction(fieldStr, test_load_snapshot, result);
					snapshot_x_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
					snapshot_y_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
					for (int j=0; j < GetNpoints(); ++j)
					{
						snapshot_x_collection_VV[i][j] = test_load_snapshot[0][j];
						snapshot_y_collection_VV[i][j] = test_load_snapshot[1][j];
					}
					collected_fine_grid_ref_qoi[i] = snapshot_y_collection_VV[i][qoi_dof];
		        }
				if (use_non_unique_up_to_two)
				{
			        for(int i = fine_grid_dir0*fine_grid_dir1; i < 2*fine_grid_dir0*fine_grid_dir1; ++i)
			        {
						// generate the correct string
						std::stringstream sstm;
						sstm << "VV" << i+1;
						std::string result = sstm.str();
						const char* rr = result.c_str();
	
				        EvaluateFunction(fieldStr, test_load_snapshot, result);
						snapshot_x_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
						snapshot_y_collection_VV[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
						for (int j=0; j < GetNpoints(); ++j)
						{
							snapshot_x_collection_VV[i][j] = test_load_snapshot[0][j];
							snapshot_y_collection_VV[i][j] = test_load_snapshot[1][j];
						}
						collected_fine_grid_ref_qoi[i] = snapshot_y_collection_VV[i][qoi_dof];
			        }
				}


				// write the FOM_qoi

				std::stringstream sstm;
				sstm << "ref_FOM_qoi.txt";
				std::string LocROM_txt = sstm.str();
				const char* outname = LocROM_txt.c_str();
				ofstream myfile (outname);
				if (myfile.is_open())
				{
					for (int i0 = 0; i0 < fine_grid_dir0*fine_grid_dir1; i0++)
					{
						myfile << collected_fine_grid_ref_qoi[i0] << "\t";
					}
					if (use_non_unique_up_to_two)
					{
						for (int i0 = fine_grid_dir0*fine_grid_dir1; i0 < 2*fine_grid_dir0*fine_grid_dir1; i0++)
						{
							myfile << collected_fine_grid_ref_qoi[i0] << "\t";
						}						
					}

					myfile.close();

				}
				else cout << "Unable to open file"; 
				

			}
		}
	}
	else
	{
		if (parameter_space_dimension == 1)
		{
			compute_snapshots(number_of_snapshots); // use something new for CMAME example
		}
		else if (parameter_space_dimension == 2)
		{
			compute_snapshots_geometry_params(); // currently not implemented
		}
	}

	cout << "snapshots available" << endl;


	DoInitialise(); 

//	sleep(10);
//	cout << "continue after sleep(10)" << endl;


	DoSolve();  // for setting dimensions
	f_bnd_size = curr_f_bnd.size();
	f_p_size = curr_f_p.size();
	f_int_size = curr_f_int.size();

	cout << "DoInitialise and DoSolve done" << endl;
	cout << "parameter_space_dimension " << parameter_space_dimension << endl;

	m_session->SetSolverInfo("SolverType", "CoupledLinearisedNS_trafoP");
	if (parameter_space_dimension == 1)
	{
		CoupledLinearNS_trafoP babyCLNS_trafo(m_session);
		babyCLNS_trafo.InitObject();
		babyCLNS_trafo.use_Newton = use_Newton;
		babyCLNS_trafo.debug_mode = debug_mode;
		collect_f_all = babyCLNS_trafo.DoTrafo(snapshot_x_collection, snapshot_y_collection, param_vector);
	}
	else if (parameter_space_dimension == 2)
	{
		do_geo_trafo(); // setting collect_f_all, making use of snapshot_x_collection, snapshot_y_collection
	}

	// insert here the route to LocalROMs
	if (m_session->DefinesParameter("use_LocROM")) // this sets how the truth system global coupling is enforced
	{
		use_LocROM = m_session->GetParameter("use_LocROM");
	}
	else
	{
		use_LocROM = 0;
	}	
	if (use_LocROM)
	{
		if (debug_mode)
		{
/*			cout << "testing a proper L2 norm availability for LocROM" << endl; // e.g. by calling L2norm_ITHACA
			for (int i=0; i<Nmax; ++i)
			{
				cout << "L2norm of snapshot " << i << " "  << L2norm_ITHACA(snapshot_x_collection[i], snapshot_y_collection[i]) << endl;
			} */
		}
		// in principle set the number of LocalClusters to use or determine a range of Clusters to consider
		// Open Q: normalize, and if so with which norm ?
		// Open Q: use phys / loc / loc_bnd_p_int / glo ?
		

		if (m_session->DefinesParameter("only_single_cluster")) 
		{
			only_single_cluster = m_session->GetParameter("only_single_cluster");	
		}
		else
		{
			only_single_cluster = 0;
		}

		if (m_session->DefinesParameter("which_single_cluster")) 
		{
			which_single_cluster = m_session->GetParameter("which_single_cluster");	
		}
		else
		{
			which_single_cluster = -1;
		}
		if (m_session->DefinesParameter("load_predef_cluster")) 
		{
			load_predef_cluster = m_session->GetParameter("load_predef_cluster");	
		}
		else
		{
			load_predef_cluster = 0;
		}


		// Test range:
//		int no_clusters = 3;
		int no_clusters = m_session->GetParameter("no_clusters");
		use_overlap_p_space = m_session->GetParameter("use_overlap_p_space");
		double optimal_CVT_energy = -1;
		Array<OneD, std::set<int> > optimal_clusters(no_clusters);
		std::srand ( unsigned ( std::time(0) ) );

//		cout << "ATTENTION: using pre-def clustering!" << endl;
	// 7er
/* 58 59 60 61 62 63 64 65 66 67 68 69 70 71
 35 38 39 40 41 42 43 44 46 47 48 49 50 51 52 53 54 55 56 57
 0 1 2 3 9 10 11 18
 7 8 14 15 16 17 22 23 24 25 26 30 31 32 33 34
 4 5 6 12 13 21
 19 20 27 28 36
 29 37 45
*/

// 8er
/* 21 29 37 45
 0 1 2 9
 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71
 23 24 25 26 30 31 32 33 34 35 38 39 40 41 42 43 44 46 47 48 49
 7 8 14 15 16 17 22
 10 11 18 19 20 27
 28 36
 3 4 5 6 12 13
*/
/*		optimal_clusters[0].insert(21);
		optimal_clusters[0].insert(29);
		optimal_clusters[0].insert(37);
		optimal_clusters[0].insert(45);

		optimal_clusters[1].insert(0);
		optimal_clusters[1].insert(1);
		optimal_clusters[1].insert(2);
		optimal_clusters[1].insert(9);

		for (int i = 50; i <= 71; ++i)
			optimal_clusters[2].insert(i);

		for (int i = 23; i <= 26; ++i)
			optimal_clusters[3].insert(i);
		for (int i = 30; i <= 35; ++i)
			optimal_clusters[3].insert(i);
		for (int i = 38; i <= 44; ++i)
			optimal_clusters[3].insert(i);
		for (int i = 46; i <= 49; ++i)
			optimal_clusters[3].insert(i);
 
		optimal_clusters[4].insert(7);
		optimal_clusters[4].insert(8);
		optimal_clusters[4].insert(14);
		optimal_clusters[4].insert(15);
		optimal_clusters[4].insert(16);
		optimal_clusters[4].insert(17);
		optimal_clusters[4].insert(22);

		optimal_clusters[5].insert(10);
		optimal_clusters[5].insert(11);
		optimal_clusters[5].insert(18);
		optimal_clusters[5].insert(19);
		optimal_clusters[5].insert(20);
		optimal_clusters[5].insert(27);

		optimal_clusters[6].insert(28);
		optimal_clusters[6].insert(36);

		optimal_clusters[7].insert(3);
		optimal_clusters[7].insert(4);
		optimal_clusters[7].insert(5);
		optimal_clusters[7].insert(6);
		optimal_clusters[7].insert(12);
		optimal_clusters[7].insert(13);
*/


		if (load_predef_cluster)
		{
		        std::string optimal_clustering_txt = "optimal_clustering.txt";
			const char* optimal_clustering_txt_t = optimal_clustering_txt.c_str();
			ifstream myfile_optimal_clustering_txt_t (optimal_clustering_txt_t);
			std::vector< std::vector<int> > all_integers;
			if (myfile_optimal_clustering_txt_t.is_open())
			{

				std::string line;
				std::vector< std::vector<int> > all_integers;
				int counter = 0;
				while ( getline( myfile_optimal_clustering_txt_t, line ) ) 
				{
				      cout << line << endl;
      				      std::istringstream is( line );
					std::vector<int> nn = std::vector<int>( std::istream_iterator<int>(is), std::istream_iterator<int>() );
					cout << "nn.size() "  << nn.size() << endl;
					for (int i = 0; i < nn.size(); ++i)
					{
						optimal_clusters[counter].insert(nn[i]);
					}
					++counter;
//				      all_integers.push_back( std::vector<int>( std::istream_iterator<int>(is), std::istream_iterator<int>() ) );
				}

	/*			for (int i = 0; i < no_clusters; ++i)
				{
					for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
					{
						myfile_optimal_clustering_txt_t << *it << "\t";
					}
					myfile_optimal_clustering_txt_t << endl;
				}
	      			myfile_optimal_clustering_txt_t.close(); */
			}
			else cout << "Unable to open file"; 

			cout << "the loaded clustering is: " << endl;
			cout << optimal_clusters.num_elements() << endl;
			for (int i1=0; i1 < optimal_clusters.num_elements(); ++i1)
			{
				for (std::set<int>::iterator it=optimal_clusters[i1].begin(); it!=optimal_clusters[i1].end(); ++it)
				{
					cout << *it << " ";
				}
				cout << endl;
			}

		}


	if (!load_predef_cluster)
	{
		for (int iter = 0; iter < 100; ++iter)
//		for (int iter = 0; iter < 0; ++iter)
		{
			Array<OneD, std::set<int> > clusters(no_clusters);
			double CVT_energy = 0;
			k_means_ITHACA(no_clusters, clusters, CVT_energy);
			if ((optimal_CVT_energy == -1) || (CVT_energy < optimal_CVT_energy))
			{
				optimal_CVT_energy = CVT_energy;
				optimal_clusters = clusters;
			}
		}
		if (debug_mode)
		{
			cout << "optimal final clustering" << endl;
			for (int i = 0; i < no_clusters; ++i)
			{
				for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
				{
					std::cout << ' ' << *it;
				}
				cout << endl;
			}
			cout << "optimal CVT_energy " << optimal_CVT_energy << endl;
		}

		// save the optimal final clustering
		cout << "optimal final clustering" << endl;
		for (int i = 0; i < no_clusters; ++i)
		{
			for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
			{
				std::cout << ' ' << *it;
			}
			cout << endl;
		}
		cout << "optimal CVT_energy " << optimal_CVT_energy << endl;

	        std::string optimal_clustering_txt = "optimal_clustering.txt";
		const char* optimal_clustering_txt_t = optimal_clustering_txt.c_str();
		ofstream myfile_optimal_clustering_txt_t (optimal_clustering_txt_t);
		if (myfile_optimal_clustering_txt_t.is_open())
		{
			for (int i = 0; i < no_clusters; ++i)
			{
				for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
				{
					myfile_optimal_clustering_txt_t << *it << "\t";
				}
				myfile_optimal_clustering_txt_t << endl;
			}
      			myfile_optimal_clustering_txt_t.close();
		}
		else cout << "Unable to open file"; 

	}


		Eigen::MatrixXd samples_cluster_association = Eigen::MatrixXd::Zero(number_of_snapshots_dir0, number_of_snapshots_dir1);
		int moving_index = 0;
		for (int i0 = 0; i0 < number_of_snapshots_dir0; i0++)
		{
			for (int i1 = 0; i1 < number_of_snapshots_dir1; i1++)
			{
				// identify cluster in which o_i is
				for (int j = 0; j < no_clusters; ++j)
				{
					if (optimal_clusters[j].count(moving_index) == 1)
						samples_cluster_association(i0, i1) = j;
				}
				moving_index++;
			}
		}

		Eigen::MatrixXd samples_cluster_association_v2 = Eigen::MatrixXd::Zero(number_of_snapshots_dir0, number_of_snapshots_dir1);		
		if (use_non_unique_up_to_two)
		{
			int moving_index_v2 = number_of_snapshots_dir0*number_of_snapshots_dir1;
			for (int i0 = 0; i0 < number_of_snapshots_dir0; i0++)
			{
				for (int i1 = 0; i1 < number_of_snapshots_dir1; i1++)
				{
					// identify cluster in which o_i is
					for (int j = 0; j < no_clusters; ++j)
					{
						if (optimal_clusters[j].count(moving_index_v2) == 1)
							samples_cluster_association_v2(i0, i1) = j;
					}
					moving_index_v2++;
				}
			}
			
		}

        std::string samples_optimal_clustering_txt = "samples_optimal_clustering.txt";
		const char* samples_optimal_clustering_txt_t = samples_optimal_clustering_txt.c_str();
		ofstream myfile_samples_optimal_clustering_txt_t (samples_optimal_clustering_txt_t);
		if (myfile_samples_optimal_clustering_txt_t.is_open())
		{
			for (int i0 = 0; i0 < number_of_snapshots_dir0; i0++)
			{
				for (int i1 = 0; i1 < number_of_snapshots_dir1; i1++)
				{
					myfile_samples_optimal_clustering_txt_t << samples_cluster_association(i0, i1)  << "\t";
				}
				myfile_samples_optimal_clustering_txt_t << endl;
			}
			if (use_non_unique_up_to_two)
			{
				for (int i0 = 0; i0 < number_of_snapshots_dir0; i0++)
				{
					for (int i1 = 0; i1 < number_of_snapshots_dir1; i1++)
					{
						myfile_samples_optimal_clustering_txt_t << samples_cluster_association_v2(i0, i1)  << "\t";
					}
					myfile_samples_optimal_clustering_txt_t << endl;
				}
			}
   			myfile_samples_optimal_clustering_txt_t.close();
		}
		else cout << "Unable to open file"; 


		// keep collect_f_all for later to continue regular execution, while a LocROM verification and validation module runs from here first
		cout << "use_fine_grid_VV " << use_fine_grid_VV << endl;
		evaluate_local_clusters(optimal_clusters);
		
	}   // 	if (use_LocROM)

	if (use_sparse_poly)
	{
		cout << "encountered use_sparse_poly, return from offline_phase" << endl;
		return ;
	}

	if (use_ANN_local)
	{
		cout << "encountered use_ANN_local, return from offline_phase" << endl;
		return ;
	}

	Eigen::BDCSVD<Eigen::MatrixXd> svd_collect_f_all(collect_f_all, Eigen::ComputeThinU);
//	cout << "svd_collect_f_all.singularValues() " << svd_collect_f_all.singularValues() << endl << endl;
	Eigen::VectorXd singular_values = svd_collect_f_all.singularValues();
	if (debug_mode)
	{
		cout << "sum singular values " << singular_values.sum() << endl << endl;
	}
	Eigen::VectorXd rel_singular_values = singular_values / singular_values.sum();
//	cout << "relative singular value percents: " << rel_singular_values << endl;
	// determine RBsize corresponding to the chosen POD_tolerance
	RBsize = 1; 
	Eigen::VectorXd cum_rel_singular_values = Eigen::VectorXd::Zero(singular_values.rows());
	for (int i = 0; i < singular_values.rows(); ++i)
	{
		cum_rel_singular_values(i) = singular_values.head(i+1).sum() / singular_values.sum();
		if (cum_rel_singular_values(i) < POD_tolerance)
		{
			RBsize = i+2;
		}		
	}
	if (debug_mode)
	{
		cout << "cumulative relative singular value percentages: " << cum_rel_singular_values << endl;
		cout << "RBsize: " << RBsize << endl;
	}
	


	Eigen::MatrixXd collect_f_all_PODmodes = svd_collect_f_all.matrixU(); // this is a local variable...
	// here probably limit to something like 99.99 percent of PODenergy, this will set RBsize
	setDBC(collect_f_all); // agnostic to RBsize
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
//	PODmodes = Eigen::MatrixXd::Zero(collect_f_all_PODmodes.rows(), collect_f_all_PODmodes.cols());
	PODmodes = Eigen::MatrixXd::Zero(collect_f_all_PODmodes.rows(), RBsize);  
	PODmodes = collect_f_all_PODmodes.leftCols(RBsize);
	set_MtM();
	cout << "RBsize: " << RBsize << endl;
	if (globally_connected == 1)
	{
		setDBC_M(collect_f_all);
	}
	if (debug_mode)
	{
		cout << "M_no_dbc_in_loc " << M_no_dbc_in_loc << endl;
		cout << "no_dbc_in_loc " << no_dbc_in_loc << endl;
		cout << "M_no_not_dbc_in_loc " << M_no_not_dbc_in_loc << endl;
		cout << "no_not_dbc_in_loc " <<	no_not_dbc_in_loc << endl;
	}
	//Eigen::VectorXd f_bnd_dbc_full_size = CLNS.f_bnd_dbc_full_size;
	// c_f_all_PODmodes_wo_dbc becomes CLNS.RB
	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc = RB;
	if (debug_mode)
	{
		cout << "c_f_all_PODmodes_wo_dbc.rows() " << c_f_all_PODmodes_wo_dbc.rows() << endl;
		cout << "c_f_all_PODmodes_wo_dbc.cols() " << c_f_all_PODmodes_wo_dbc.cols() << endl;
	}



	time(&timer_1);
	gen_phys_base_vecs();
	time(&timer_2);
	if (debug_mode)	
	{
		cout << "finished gen_phys_base_vecs in " << difftime(timer_2, timer_1) << " seconds" << endl;
	}

	if (use_ANN)
	{
		cout << "encountered use_ANN, writing train data and return from offline" << endl;
		train_data_x = Eigen::MatrixXd::Zero(RBsize, Nmax);
		train_data_y = Eigen::MatrixXd::Zero(RBsize, Nmax);
		for (int i = 0; i < Nmax; ++i) // for all original snapshots
		{
			// project original snapshot data onto the POD basis
			Eigen::MatrixXd snapshot_xy_proj = project_onto_basis(snapshot_x_collection[i], snapshot_y_collection[i]);
			//			cout << "snapshot_xy_proj.rows() " << snapshot_xy_proj.rows() << " snapshot_xy_proj.cols() " << snapshot_xy_proj.cols() << endl;
			train_data_x.col(i) = snapshot_xy_proj.col(0);
			train_data_y.col(i) = snapshot_xy_proj.col(1);
		}
		// write the training data to *.txt files


		ofstream myfileANN ("trainANN.txt");
		if (myfileANN.is_open())
		{
			for (int i = 0; i < RBsize; i++)
			{
				for (int j = 0; j < Nmax; j++)
				{
					myfileANN << std::setprecision(17) << train_data_x(i,j) << "\t";
				}
				myfileANN << "\n";
			}
			for (int i = 0; i < RBsize; i++)
			{
				for (int j = 0; j < Nmax; j++)
				{
					myfileANN << std::setprecision(17) << train_data_y(i,j) << "\t";
				}
				myfileANN << "\n";
			}
            myfileANN.close();
		}
		else cout << "Unable to open file"; 


		return;
	}

	if (parameter_space_dimension == 1)
	{
		gen_proj_adv_terms();
	}
	else if (parameter_space_dimension == 2)
	{
		gen_proj_adv_terms_2d();
	}
	if (debug_mode)	
	{
		cout << "finished gen_proj_adv_terms " << endl;
	}
	if (parameter_space_dimension == 1)
	{
		gen_reference_matrices();
	}
	else if (parameter_space_dimension == 2)
	{
		gen_reference_matrices_2d();
	}
	if (debug_mode)	
	{
		cout << "finished gen_reference_matrices " << endl;
	}
    }

    Eigen::MatrixXd CoupledLinearNS_TT::gen_affine_mat_proj(double current_nu)
    {
	Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
	for (int i = 0; i < RBsize; ++i)
	{
		recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_xy_projected(i,0) + adv_mats_proj_y[i] * curr_xy_projected(i,1);
	}
	Eigen::MatrixXd affine_mat_proj = the_const_one_proj + current_nu * the_ABCD_one_proj + recovered_affine_adv_mat_proj_xy;



/*	if (debug_mode)
	{
		Eigen::MatrixXd affine_mat = Get_complete_matrix();
		cout << "affine_mat.rows() " << affine_mat.rows() << " affine_mat.cols() " << affine_mat.cols() << endl;
		cout << "affine_mat_proj.rows() " << affine_mat_proj.rows() << " affine_mat_proj.cols() " << affine_mat_proj.cols() << endl;
		cout << "RB.rows() " << RB.rows() << " RB.cols() " << RB.cols() << endl;
		cout << "PODmodes.rows() " << PODmodes.rows() << " PODmodes.cols() " << PODmodes.cols() << endl;
		cout << "affine_mat.norm() "  << affine_mat.norm() << endl;
		cout << "affine_mat_proj.norm() " << affine_mat_proj.norm() << endl;
		Eigen::MatrixXd reproj_affine_mat = Eigen::MatrixXd::Zero(RB.rows(), RB.rows());
		reproj_affine_mat = RB * affine_mat_proj * RB.transpose();
		cout << "reproj_affine_mat.norm() " << reproj_affine_mat.norm() << endl;
		Eigen::MatrixXd affine_matrix_simplified = remove_cols_and_rows(affine_mat, elem_loc_dbc);
		cout << "affine_matrix_simplified.norm() "  << affine_matrix_simplified.norm() << endl;
		Eigen::MatrixXd reduced_affine_mat = Eigen::MatrixXd::Zero(RB.cols(), RB.cols());
		reduced_affine_mat = RB.transpose() * affine_matrix_simplified * RB;
		cout << "reduced_affine_mat.norm() "  << reduced_affine_mat.norm() << endl;

	} */

	return affine_mat_proj;
    }

    void CoupledLinearNS_TT::evaluate_local_clusters(Array<OneD, std::set<int> > optimal_clusters)
    {

	// determine the local cluster projection space and go through each offline phase
	int no_clusters = optimal_clusters.num_elements();
	cout << "no_clusters " << no_clusters << endl;
	if (use_fine_grid_VV)
		associate_VV_to_clusters(optimal_clusters);
	//int use_overlap_p_space = 1;
	Array<OneD, std::set<int> > optimal_clusters_mod = Array<OneD, std::set<int> >(no_clusters);
	Array<OneD, std::set<int> > optimal_clusters_orig = Array<OneD, std::set<int> >(no_clusters);
	optimal_clusters_orig = optimal_clusters;
	if (use_overlap_p_space)
	{
		// adjust local clusters accordingly
		for (int i = 0; i < no_clusters; ++i)
		{
			int index_all = 0;
			for (int j0 = 0; j0 < number_of_snapshots_dir0; ++j0)
			{
				for (int j1 = 0; j1 < number_of_snapshots_dir1; ++j1)
				{
					// is it a neighbour of cluster i?
					if (!optimal_clusters[i].count(index_all))
					{
						int is_neigh = 0;
						int neigh0;
						if (j1 > 0)
						{
							neigh0 = j1-1 + number_of_snapshots_dir1*j0;
							//cout << "neigh0 " << neigh0 << endl;
							if (optimal_clusters[i].count(neigh0))
								is_neigh = 1;
						}
						if (j1 < number_of_snapshots_dir1-1)
						{
							neigh0 = j1+1 + number_of_snapshots_dir1*j0;
							//cout << "neigh0 " << neigh0 << endl;
							if (optimal_clusters[i].count(neigh0))
								is_neigh = 1;
						}
						if (j0 > 0)
						{
							neigh0 = j1 + number_of_snapshots_dir1*(j0-1);
							//cout << "neigh0 " << neigh0 << endl;
							if (optimal_clusters[i].count(neigh0))
								is_neigh = 1;
						}
						if (j0 < number_of_snapshots_dir0-1)
						{
							neigh0 = j1 + number_of_snapshots_dir1*(j0+1);
							//cout << "neigh0 " << neigh0 << endl;
							if (optimal_clusters[i].count(neigh0))
								is_neigh = 1;
						}
						if (is_neigh)
						{
							optimal_clusters_mod[i].insert(index_all);
							//if (debug_mode)
							//	cout << " inserting to clust " << i << " the index " << index_all << endl;
						}
					}
					else
						optimal_clusters_mod[i].insert(index_all);

					if (use_non_unique_up_to_two)
					{
						// same again with Nmax/2 offset: is it a neighbour of cluster i?
						if (!optimal_clusters[i].count(index_all + Nmax/2))
						{
							int is_neigh = 0;
							int neigh0;
							if (j1 > 0)
							{
								neigh0 = j1-1 + number_of_snapshots_dir1*j0;
								//cout << "neigh0 " << neigh0 << endl;
								if (optimal_clusters[i].count(neigh0 + Nmax/2))
									is_neigh = 1;
							}
							if (j1 < number_of_snapshots_dir1-1)
							{
								neigh0 = j1+1 + number_of_snapshots_dir1*j0;
								//cout << "neigh0 " << neigh0 << endl;
								if (optimal_clusters[i].count(neigh0 + Nmax/2))
									is_neigh = 1;
							}
							if (j0 > 0)
							{
								neigh0 = j1 + number_of_snapshots_dir1*(j0-1);
								//cout << "neigh0 " << neigh0 << endl;
								if (optimal_clusters[i].count(neigh0 + Nmax/2))
									is_neigh = 1;
							}
							if (j0 < number_of_snapshots_dir0-1)
							{
								neigh0 = j1 + number_of_snapshots_dir1*(j0+1);
								//cout << "neigh0 " << neigh0 << endl;
								if (optimal_clusters[i].count(neigh0 + Nmax/2))
									is_neigh = 1;
							}
							if (is_neigh)
							{
								optimal_clusters_mod[i].insert(index_all + Nmax/2);
								//if (debug_mode)
								//	cout << " inserting to clust " << i << " the index " << index_all << endl;
							}
						}
						else
							optimal_clusters_mod[i].insert(index_all + Nmax/2);

					}
					index_all++;
				}
			}
		}
		optimal_clusters = optimal_clusters_mod;
		if (debug_mode)
		{
			cout << "optimal final clustering after possible transition regions " << endl;
			for (int i = 0; i < no_clusters; ++i)
			{
				for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
				{
					std::cout << ' ' << *it;
				}
				cout << endl;
			}
		}
	} // if (use_overlap_p_space)
//	for (int i = 0; i < no_clusters; ++i)
	if (only_single_cluster)
	{
	for (int i = which_single_cluster; i < which_single_cluster+1; ++i)
	{

		Eigen::MatrixXd local_collect_f_all_orig = Eigen::MatrixXd::Zero( collect_f_all.rows() , optimal_clusters_orig[i].size() );
		int j = 0;
		for (std::set<int>::iterator it=optimal_clusters_orig[i].begin(); it!=optimal_clusters_orig[i].end(); ++it)
		{
			local_collect_f_all_orig.col(j) = collect_f_all.col(*it);
			j++;
		}

//		Eigen::MatrixXd local_collect_f_all = Eigen::MatrixXd::Zero( collect_f_all.rows() , optimal_clusters_orig[i].size() );
		Eigen::MatrixXd local_collect_f_all_add = Eigen::MatrixXd::Zero( collect_f_all.rows() , optimal_clusters[i].size() - optimal_clusters_orig[i].size() );
		j = 0;
		for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
		{
			if (!optimal_clusters_orig[i].count(*it))
			{
				local_collect_f_all_add.col(j) = collect_f_all.col(*it); 
				j++;
			}
		}
		for (int j = 0; j < optimal_clusters[i].size(); ++j)
		{
//			local_collect_f_all.col(j) = collect_f_all.col(optimal_clusters[i][j]);
		}
		if (use_overlap_p_space)
			run_local_ROM_offline_add_transition(local_collect_f_all_orig,local_collect_f_all_add, which_single_cluster);
		else
			run_local_ROM_offline(local_collect_f_all_orig);
		run_local_ROM_online(optimal_clusters[i], i);
//		run_local_ROM_online(optimal_clusters_orig[i], i);
	}
	}
	else
	for (int i = 0; i < no_clusters; ++i)
	{

		Eigen::MatrixXd local_collect_f_all_orig = Eigen::MatrixXd::Zero( collect_f_all.rows() , optimal_clusters_orig[i].size() );
		int j = 0;
		for (std::set<int>::iterator it=optimal_clusters_orig[i].begin(); it!=optimal_clusters_orig[i].end(); ++it)
		{
			local_collect_f_all_orig.col(j) = collect_f_all.col(*it);
			j++;
		}

//		Eigen::MatrixXd local_collect_f_all = Eigen::MatrixXd::Zero( collect_f_all.rows() , optimal_clusters_orig[i].size() );
		Eigen::MatrixXd local_collect_f_all_add = Eigen::MatrixXd::Zero( collect_f_all.rows() , optimal_clusters[i].size() - optimal_clusters_orig[i].size() );
		j = 0;
		for (std::set<int>::iterator it=optimal_clusters[i].begin(); it!=optimal_clusters[i].end(); ++it)
		{
			if (!optimal_clusters_orig[i].count(*it))
			{
				local_collect_f_all_add.col(j) = collect_f_all.col(*it); 
				j++;
			}
		}
		for (int j = 0; j < optimal_clusters[i].size(); ++j)
		{
//			local_collect_f_all.col(j) = collect_f_all.col(optimal_clusters[i][j]);
		}
		if (use_overlap_p_space)
			run_local_ROM_offline_add_transition(local_collect_f_all_orig,local_collect_f_all_add, i);
		else
			run_local_ROM_offline(local_collect_f_all_orig);
		if (!use_ANN_local)
		{
			run_local_ROM_online(optimal_clusters[i], i);
		}
		else
		{
			compute_ANN_approx_cluster(i);
		}
//		run_local_ROM_online(optimal_clusters_orig[i], i);
	}
    }

    void CoupledLinearNS_TT::run_local_ROM_offline(Eigen::MatrixXd collect_f_all)
   {
	Eigen::BDCSVD<Eigen::MatrixXd> svd_collect_f_all(collect_f_all, Eigen::ComputeThinU);
	Eigen::VectorXd singular_values = svd_collect_f_all.singularValues();
	if (debug_mode)
	{
		cout << "sum singular values " << singular_values.sum() << endl << endl;
	}
	Eigen::VectorXd rel_singular_values = singular_values / singular_values.sum();
	RBsize = 1; 
	Eigen::VectorXd cum_rel_singular_values = Eigen::VectorXd::Zero(singular_values.rows());
	for (int i = 0; i < singular_values.rows(); ++i)
	{
		cum_rel_singular_values(i) = singular_values.head(i+1).sum() / singular_values.sum();
		if (cum_rel_singular_values(i) < POD_tolerance)
		{
			RBsize = i+2;
		}		
	}
	if (debug_mode)
	{
		cout << "cumulative relative singular value percentages: " << std::setprecision(17) << cum_rel_singular_values << endl;
		cout << "RBsize: " << RBsize << endl;
	}
	


	Eigen::MatrixXd collect_f_all_PODmodes = svd_collect_f_all.matrixU(); // this is a local variable...
	// here probably limit to something like 99.99 percent of PODenergy, this will set RBsize
	setDBC(collect_f_all); // agnostic to RBsize
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
//	PODmodes = Eigen::MatrixXd::Zero(collect_f_all_PODmodes.rows(), collect_f_all_PODmodes.cols());
	PODmodes = Eigen::MatrixXd::Zero(collect_f_all_PODmodes.rows(), RBsize);  
	PODmodes = collect_f_all_PODmodes.leftCols(RBsize);
	set_MtM();
	cout << "RBsize: " << RBsize << endl;
	if (globally_connected == 1)
	{
		setDBC_M(collect_f_all);
	}
	if (debug_mode)
	{
		cout << "M_no_dbc_in_loc " << M_no_dbc_in_loc << endl;
		cout << "no_dbc_in_loc " << no_dbc_in_loc << endl;
		cout << "M_no_not_dbc_in_loc " << M_no_not_dbc_in_loc << endl;
		cout << "no_not_dbc_in_loc " <<	no_not_dbc_in_loc << endl;
	}
	//Eigen::VectorXd f_bnd_dbc_full_size = CLNS.f_bnd_dbc_full_size;
	// c_f_all_PODmodes_wo_dbc becomes CLNS.RB
	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc = RB;
	if (debug_mode)
	{
		cout << "c_f_all_PODmodes_wo_dbc.rows() " << c_f_all_PODmodes_wo_dbc.rows() << endl;
		cout << "c_f_all_PODmodes_wo_dbc.cols() " << c_f_all_PODmodes_wo_dbc.cols() << endl;
	}
	gen_phys_base_vecs();
	if (debug_mode)	
	{
//		cout << "finished gen_phys_base_vecs in " << difftime(timer_2, timer_1) << " seconds" << endl;
	}
	if (parameter_space_dimension == 1)
	{
		gen_proj_adv_terms();
	}
	else if (parameter_space_dimension == 2)
	{
		gen_proj_adv_terms_2d();
	}
	if (debug_mode)	
	{
		cout << "finished gen_proj_adv_terms " << endl;
	}
	if (parameter_space_dimension == 1)
	{
		gen_reference_matrices();
	}
	else if (parameter_space_dimension == 2)
	{
		gen_reference_matrices_2d();
	}
	if (debug_mode)	
	{
		cout << "finished gen_reference_matrices " << endl;
	}

   }


    void CoupledLinearNS_TT::run_local_ROM_offline_add_transition(Eigen::MatrixXd collect_f_all_orig, Eigen::MatrixXd collect_f_all_add, int no_clust)
   {
//	Eigen::MatrixXd collect_f_all_orig_mod = Eigen::MatrixXd::Zero(collect_f_all_orig.rows(), collect_f_all_orig.cols() + collect_f_all_add.cols()); 
//	collect_f_all_orig_mod << collect_f_all_orig, collect_f_all_add;
	Eigen::BDCSVD<Eigen::MatrixXd> svd_collect_f_all(collect_f_all_orig, Eigen::ComputeThinU);
	Eigen::VectorXd singular_values = svd_collect_f_all.singularValues();
	if (debug_mode)
	{
		cout << "sum singular values " << singular_values.sum() << endl << endl;
	}
	Eigen::VectorXd rel_singular_values = singular_values / singular_values.sum();
	RBsize = 1; 
	Eigen::VectorXd cum_rel_singular_values = Eigen::VectorXd::Zero(singular_values.rows());
	for (int i = 0; i < singular_values.rows(); ++i)
	{
		cum_rel_singular_values(i) = singular_values.head(i+1).sum() / singular_values.sum();
		if (cum_rel_singular_values(i) < POD_tolerance)
		{
			RBsize = i+2;
		}		
	}
	if (debug_mode)
	{
		cout << "cumulative relative singular value percentages: " << std::setprecision(17) << cum_rel_singular_values << endl;
		cout << "RBsize: " << RBsize << endl;
	}

	Eigen::MatrixXd collect_f_all_PODmodes_int = svd_collect_f_all.matrixU(); // this is a local variable...
	Eigen::MatrixXd PODmodes_int = Eigen::MatrixXd::Zero(collect_f_all_PODmodes_int.rows(), RBsize); 
	PODmodes_int = collect_f_all_PODmodes_int.leftCols(RBsize);

	cout << "now adding" << endl;	
	collect_f_all_add = collect_f_all_add - PODmodes_int * PODmodes_int.transpose() * collect_f_all_add;

	
//	Eigen::MatrixXd cfac(collect_f_all_PODmodes.rows(), collect_f_all_PODmodes.cols()+collect_f_all_add.cols());
	Eigen::MatrixXd cfac(PODmodes_int.rows(), collect_f_all_add.cols());
	cfac << collect_f_all_add;
	Eigen::BDCSVD<Eigen::MatrixXd> svd_collect_f_allc(cfac, Eigen::ComputeThinU);
	singular_values = svd_collect_f_allc.singularValues();
	if (debug_mode)
	{
		cout << "sum singular values " << singular_values.sum() << endl << endl;
	}
	rel_singular_values = singular_values / singular_values.sum();
       	int RBsize_c = 1; 
	cum_rel_singular_values = Eigen::VectorXd::Zero(singular_values.rows());

	double overlap_POD_part = m_session->GetParameter("overlap_POD_part");	

	for (int i = 0; i < singular_values.rows(); ++i)
	{
		cum_rel_singular_values(i) = singular_values.head(i+1).sum() / singular_values.sum();
		if (cum_rel_singular_values(i) < POD_tolerance - POD_tolerance/overlap_POD_part)
//		if (cum_rel_singular_values(i) < POD_tolerance - POD_tolerance/2)
		{
			RBsize_c = i+2;
		}		
	}
	if (debug_mode)
	{
		cout << "cumulative relative singular value percentages: " << std::setprecision(17) << cum_rel_singular_values << endl;
		cout << "RBsize_c: " << RBsize_c << endl;
	}

//	RBsize_c = 2; 

	RBsize += RBsize_c;
	Eigen::MatrixXd collect_f_all_PODmodes(PODmodes_int.rows(), RBsize);
	collect_f_all_PODmodes << PODmodes_int, svd_collect_f_allc.matrixU().leftCols(RBsize_c); // this is a local variable...
	cout << "are total modes orth?? " << collect_f_all_PODmodes.transpose() * collect_f_all_PODmodes <<  endl;	

	// does yet another POD help?
/*	Eigen::MatrixXd cfac_ya(PODmodes_int.rows(), collect_f_all_PODmodes.cols());
	cfac_ya = collect_f_all_PODmodes;
	Eigen::BDCSVD<Eigen::MatrixXd> svd_cfac_ya(cfac_ya, Eigen::ComputeThinU);
	singular_values = svd_cfac_ya.singularValues();
	if (debug_mode)
	{
		cout << "sum singular values " << singular_values.sum() << endl << endl;
	}
	rel_singular_values = singular_values / singular_values.sum();
	cum_rel_singular_values = Eigen::VectorXd::Zero(singular_values.rows());
	for (int i = 0; i < singular_values.rows(); ++i)
	{
		cum_rel_singular_values(i) = singular_values.head(i+1).sum() / singular_values.sum();
	}
	cout << "cumulative relative singular value percentages: " << std::setprecision(17) << cum_rel_singular_values << endl;
	collect_f_all_PODmodes = svd_cfac_ya.matrixU();
*/

	cout << "dimension of proj. space " << collect_f_all_PODmodes.rows() << " by " << collect_f_all_PODmodes.cols() << endl;
	cout << "dimension of collected snapshots " << collect_f_all.rows() << " by " << collect_f_all.cols() << endl;
	// double check the projection error of all initial snapshots
	for (int iter_index = 0; iter_index < Nmax; ++iter_index)
	{
/*		cout << "projection error of snapshot number " << iter_index << endl;
		Eigen::VectorXd curr_f_all = collect_f_all.col(iter_index);
		Eigen::VectorXd proj_curr_f_all = collect_f_all_PODmodes.transpose() * curr_f_all;
		Eigen::VectorXd reproj_curr_f_all = collect_f_all_PODmodes * proj_curr_f_all;
		cout << "curr_f_all.norm() " << curr_f_all.norm() << endl;
		cout << "reproj_curr_f_all.norm() " << reproj_curr_f_all.norm() << endl;
		cout << "rel error " << (reproj_curr_f_all - curr_f_all).norm() / curr_f_all.norm() << endl; */
	}


	// here probably limit to something like 99.99 percent of PODenergy, this will set RBsize
	setDBC(collect_f_all); // agnostic to RBsize
	Array<OneD, MultiRegions::ExpListSharedPtr> m_fields = UpdateFields();
        int  nel  = m_fields[0]->GetNumElmts(); // number of spectral elements
//	PODmodes = Eigen::MatrixXd::Zero(collect_f_all_PODmodes.rows(), collect_f_all_PODmodes.cols());
//	PODmodes = Eigen::MatrixXd::Zero(collect_f_all_PODmodes.rows(), RBsize);  
//	PODmodes = collect_f_all_PODmodes.leftCols(RBsize);
	PODmodes = collect_f_all_PODmodes;
	set_MtM();
	cout << "RBsize: " << RBsize << endl;
	if (globally_connected == 1)
	{
		setDBC_M(collect_f_all);
	}
	if (debug_mode)
	{
		cout << "M_no_dbc_in_loc " << M_no_dbc_in_loc << endl;
		cout << "no_dbc_in_loc " << no_dbc_in_loc << endl;
		cout << "M_no_not_dbc_in_loc " << M_no_not_dbc_in_loc << endl;
		cout << "no_not_dbc_in_loc " <<	no_not_dbc_in_loc << endl;
	}
	//Eigen::VectorXd f_bnd_dbc_full_size = CLNS.f_bnd_dbc_full_size;
	// c_f_all_PODmodes_wo_dbc becomes CLNS.RB
	Eigen::MatrixXd c_f_all_PODmodes_wo_dbc = RB;
	if (debug_mode)
	{
		cout << "c_f_all_PODmodes_wo_dbc.rows() " << c_f_all_PODmodes_wo_dbc.rows() << endl;
		cout << "c_f_all_PODmodes_wo_dbc.cols() " << c_f_all_PODmodes_wo_dbc.cols() << endl;
	}
	gen_phys_base_vecs();
	if (debug_mode)	
	{
//		cout << "finished gen_phys_base_vecs in " << difftime(timer_2, timer_1) << " seconds" << endl;
	}
	
	if (use_ANN_local)
	{
		cout << "encountered use_ANN_local, writing train data and return from run_local_ROM_offline_add_transition for current cluster " << endl;
		train_data_x = Eigen::MatrixXd::Zero(RBsize, Nmax);
		train_data_y = Eigen::MatrixXd::Zero(RBsize, Nmax);
		for (int i = 0; i < Nmax; ++i) // for all original snapshots
		{
			// project original snapshot data onto the POD basis
			Eigen::MatrixXd snapshot_xy_proj = project_onto_basis(snapshot_x_collection[i], snapshot_y_collection[i]);
			//			cout << "snapshot_xy_proj.rows() " << snapshot_xy_proj.rows() << " snapshot_xy_proj.cols() " << snapshot_xy_proj.cols() << endl;
			train_data_x.col(i) = snapshot_xy_proj.col(0);
			train_data_y.col(i) = snapshot_xy_proj.col(1);
		}
		// write the training data to *.txt files
		std::stringstream sstm_l11;
		sstm_l11 << "trainANN_" << no_clust << ".txt";
		std::string result_l11 = sstm_l11.str();
		const char* rr_l11 = result_l11.c_str();
                 std::ofstream myfileANN (rr_l11);
		if (myfileANN.is_open())
		{
			for (int i = 0; i < RBsize; i++)
			{
				for (int j = 0; j < Nmax; j++)
				{
					myfileANN << std::setprecision(17) << train_data_x(i,j) << "\t";
				}
				myfileANN << "\n";
			}
			for (int i = 0; i < RBsize; i++)
			{
				for (int j = 0; j < Nmax; j++)
				{
					myfileANN << std::setprecision(17) << train_data_y(i,j) << "\t";
				}
				myfileANN << "\n";
			}
	         myfileANN.close();
		}
		else cout << "Unable to open file"; 
		cout << "finished writing train data" << endl;
		return;
	}
	
	if (parameter_space_dimension == 1)
	{
		gen_proj_adv_terms();
	}
	else if (parameter_space_dimension == 2)
	{
		gen_proj_adv_terms_2d();
	}
	if (debug_mode)	
	{
		cout << "finished gen_proj_adv_terms " << endl;
	}
	if (parameter_space_dimension == 1)
	{
		gen_reference_matrices();
	}
	else if (parameter_space_dimension == 2)
	{
		gen_reference_matrices_2d();
	}
	if (debug_mode)	
	{
		cout << "finished gen_reference_matrices " << endl;
	}

   }


    void CoupledLinearNS_TT::k_means_ITHACA(int no_clusters, Array<OneD, std::set<int> > &clusters, double &CVT_energy)
    {
	// should run many times with different (random) start values

	std::vector<int> myvector;
	for (int i=0; i<Nmax; ++i) myvector.push_back(i); 
	std::random_shuffle ( myvector.begin(), myvector.end() );	
	if (debug_mode)
	{
/*		std::cout << "myvector contains:";
		for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it)
			std::cout << ' ' << *it;
		std::cout << '\n'; */
	}

	// what to use as a cell equivalent - try the one by J. Burkardt, but it cannot change dynamically...
	// can use an array of sets, since 'no_clusters' is fixed
	// Array<OneD, std::set<int> > clusters(no_clusters); // the j-th myvector entry is the initial centroid of the j-th cluster
	for (int i = 0; i < Nmax; ++i)
	{
		Array<OneD, double> distances(no_clusters);
		for (int j = 0; j < no_clusters; ++j)
		{
			// compute distance of i-th snapshot to j-th centroid
			distances[j] = L2norm_abs_error_ITHACA(snapshot_x_collection[i], snapshot_y_collection[i], snapshot_x_collection[myvector[j]], snapshot_y_collection[myvector[j]]) * L2norm_abs_error_ITHACA(snapshot_x_collection[i], snapshot_y_collection[i], snapshot_x_collection[myvector[j]], snapshot_y_collection[myvector[j]]); 
		}
		// find the minimum distance and assign i-th snapshot to appropriate cluster
		int minElementIndex = std::min_element(distances.begin(),distances.end()) - distances.begin();
		double minElement = *std::min_element(distances.begin(), distances.end());
		if (debug_mode)
		{
		/*	std::cout << "minElementIndex:" << minElementIndex << ", minElement:" << minElement << '\n';
			std::cout << "distances contains:";
			for (Array<OneD, double>::iterator it=distances.begin(); it!=distances.end(); ++it)
				std::cout << ' ' << *it;
			std::cout << endl << "the correct initial cluster centroid index is thus: " << myvector[minElementIndex] << endl; */
		}
		clusters[minElementIndex].insert(i);
	}
	if (debug_mode)
	{
/*		for (int i = 0; i < no_clusters; ++i)
		{
			for (std::set<int>::iterator it=clusters[i].begin(); it!=clusters[i].end(); ++it)
			{
				std::cout << ' ' << *it;
			}
			cout << endl;
		} */
	}
	Array<OneD, Array< OneD, NekDouble > > mean_kj_x(no_clusters);
	Array<OneD, Array< OneD, NekDouble > > mean_kj_y(no_clusters);
	for (int iter = 0; iter < 10; iter++)
	{
		// find new centroids as mean of vectors in the clusters
		for (int i = 0; i < no_clusters; ++i)
		{
			mean_kj_x[i] = Array< OneD, NekDouble >(snapshot_x_collection[0].num_elements(),0.0);
			mean_kj_y[i] = Array< OneD, NekDouble >(snapshot_y_collection[0].num_elements(),0.0);
			int cluster_size = clusters[i].size();
			for (std::set<int>::iterator it=clusters[i].begin(); it!=clusters[i].end(); ++it)
			{
				// std::cout << ' ' << *it;
				for (int j = 0; j < snapshot_x_collection[0].num_elements(); ++j)
				{
					mean_kj_x[i][j] += (1.0/cluster_size) * snapshot_x_collection[*it][j];
					mean_kj_y[i][j] += (1.0/cluster_size) * snapshot_y_collection[*it][j];
				}		
			}
		}
		// re-determine corresponding clusters
		Array<OneD, std::set<int> > clusters(no_clusters); 
		for (int i = 0; i < Nmax; ++i)
		{
			Array<OneD, double> distances(no_clusters);
			for (int j = 0; j < no_clusters; ++j)
			{
				// compute distance of i-th snapshot to j-th centroid
				distances[j] = L2norm_abs_error_ITHACA(snapshot_x_collection[i], snapshot_y_collection[i], mean_kj_x[j], mean_kj_y[j]) * L2norm_abs_error_ITHACA(snapshot_x_collection[i], snapshot_y_collection[i], mean_kj_x[j], mean_kj_y[j]); 
			}
			// find the minimum distance and assign i-th snapshot to appropriate cluster
			int minElementIndex = std::min_element(distances.begin(),distances.end()) - distances.begin();
			double minElement = *std::min_element(distances.begin(), distances.end());
			if (debug_mode)
			{
		/*		std::cout << "minElementIndex:" << minElementIndex << ", minElement:" << minElement << '\n';
				std::cout << "distances contains:";
				for (Array<OneD, double>::iterator it=distances.begin(); it!=distances.end(); ++it)
					std::cout << ' ' << *it;
				std::cout << endl << "the correct cluster is thus: " << minElementIndex << endl; */
			}
			clusters[minElementIndex].insert(i);
		}
		if (debug_mode)
		{
/*			for (int i = 0; i < no_clusters; ++i)
			{
				for (std::set<int>::iterator it=clusters[i].begin(); it!=clusters[i].end(); ++it)
				{
					std::cout << ' ' << *it;
				}
				cout << endl;
			} */
		}
	}
	if (debug_mode)
	{
/*		cout << "final clustering" << endl;
		for (int i = 0; i < no_clusters; ++i)
		{
			for (std::set<int>::iterator it=clusters[i].begin(); it!=clusters[i].end(); ++it)
			{
				std::cout << ' ' << *it;
			}
			cout << endl;
		}*/
	}

	// CVT energy and silhouette score
	// double CVT_energy = 0.0;
	for (int i = 0; i < no_clusters; ++i)
	{
		for (std::set<int>::iterator it=clusters[i].begin(); it!=clusters[i].end(); ++it)
		{
			CVT_energy += L2norm_abs_error_ITHACA(snapshot_x_collection[*it], snapshot_y_collection[*it], mean_kj_x[i], mean_kj_y[i]) * L2norm_abs_error_ITHACA(snapshot_x_collection[*it], snapshot_y_collection[*it], mean_kj_x[i], mean_kj_y[i]); 
		}
	}
	if (debug_mode)
	{
	//	cout << "CVT energy " << CVT_energy << endl;
	}

    }

    double CoupledLinearNS_TT::L2norm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y )
    {
	Array< OneD, NekDouble > x_difference(component1_x.num_elements());
	Array< OneD, NekDouble > y_difference(component1_y.num_elements());
	for (int i = 0; i < component1_y.num_elements(); ++i)
	{
		x_difference[i] = component1_x[i] - component2_x[i];
		y_difference[i] = component1_y[i] - component2_y[i];
	}
	double result = L2norm_ITHACA(x_difference, y_difference);
	return result;
    }

    double CoupledLinearNS_TT::Linfnorm_abs_error_ITHACA( Array< OneD, NekDouble > component1_x, Array< OneD, NekDouble > component1_y, Array< OneD, NekDouble > component2_x, Array< OneD, NekDouble > component2_y )
    {
	Array< OneD, NekDouble > x_difference(component1_x.num_elements());
	Array< OneD, NekDouble > y_difference(component1_y.num_elements());
	for (int i = 0; i < component1_y.num_elements(); ++i)
	{
		x_difference[i] = component1_x[i] - component2_x[i];
		y_difference[i] = component1_y[i] - component2_y[i];
	}
	double result = Linfnorm_ITHACA(x_difference, y_difference);
	return result;
    }

    Eigen::MatrixXd CoupledLinearNS_TT::gen_affine_mat_proj_2d(double current_nu, double w)
    {
	// avail datastructure	Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > adv_mats_proj_x_2d;
//	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_const_one_proj_2d;
//	Array<OneD, Array<OneD, Eigen::MatrixXd > > the_ABCD_one_proj_2d;
	Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);
	Eigen::MatrixXd recovered_press_proj = Eigen::MatrixXd::Zero(RBsize, RBsize);
	Eigen::MatrixXd recovered_ABCD_proj = Eigen::MatrixXd::Zero(RBsize, RBsize);
	for (int index_elem = 0; index_elem < number_elem_trafo; ++index_elem)
	{
		double detT = Geo_T(w, index_elem, 0);
		double Ta = Geo_T(w, index_elem, 1);
		double Tb = Geo_T(w, index_elem, 2);
		double Tc = Geo_T(w, index_elem, 3);
		double Td = Geo_T(w, index_elem, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td;
		for (int i = 0; i < RBsize; ++i)
		{
//			recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_xy_projected(i,0) + adv_mats_proj_y[i] * curr_xy_projected(i,1);
			recovered_affine_adv_mat_proj_xy += detT * curr_xy_projected(i,0) * (Ta * adv_mats_proj_x_2d[i][index_elem][0] + Tc * adv_mats_proj_x_2d[i][index_elem][1]) + detT * curr_xy_projected(i,1) * (Tb * adv_mats_proj_y_2d[i][index_elem][0] + Td * adv_mats_proj_y_2d[i][index_elem][1]);

		}
		recovered_ABCD_proj += detT * (c00 * the_ABCD_one_proj_2d[index_elem][0] + c01*(the_ABCD_one_proj_2d[index_elem][1] + the_ABCD_one_proj_2d[index_elem][2]) + c11*the_ABCD_one_proj_2d[index_elem][3]);
		recovered_press_proj += detT * (Ta * the_const_one_proj_2d[index_elem][0] + Tc * the_const_one_proj_2d[index_elem][1] + Tb * the_const_one_proj_2d[index_elem][2] + Td * the_const_one_proj_2d[index_elem][3]);

	}

	Eigen::MatrixXd affine_mat_proj = recovered_press_proj + current_nu * recovered_ABCD_proj + recovered_affine_adv_mat_proj_xy;
//	Eigen::MatrixXd affine_mat_proj = the_const_one_proj + current_nu * the_ABCD_one_proj + recovered_affine_adv_mat_proj_xy;

/*	if (debug_mode)
	{
		Eigen::MatrixXd affine_mat = Get_complete_matrix();
		cout << "affine_mat.rows() " << affine_mat.rows() << " affine_mat.cols() " << affine_mat.cols() << endl;
		cout << "affine_mat_proj.rows() " << affine_mat_proj.rows() << " affine_mat_proj.cols() " << affine_mat_proj.cols() << endl;
		cout << "RB.rows() " << RB.rows() << " RB.cols() " << RB.cols() << endl;
		cout << "PODmodes.rows() " << PODmodes.rows() << " PODmodes.cols() " << PODmodes.cols() << endl;
		cout << "affine_mat.norm() "  << affine_mat.norm() << endl;
		cout << "affine_mat_proj.norm() " << affine_mat_proj.norm() << endl;
		Eigen::MatrixXd reproj_affine_mat = Eigen::MatrixXd::Zero(RB.rows(), RB.rows());
		reproj_affine_mat = RB * affine_mat_proj * RB.transpose();
		cout << "reproj_affine_mat.norm() " << reproj_affine_mat.norm() << endl;
		Eigen::MatrixXd affine_matrix_simplified = remove_cols_and_rows(affine_mat, elem_loc_dbc);
		cout << "affine_matrix_simplified.norm() "  << affine_matrix_simplified.norm() << endl;
		Eigen::MatrixXd reduced_affine_mat = Eigen::MatrixXd::Zero(RB.cols(), RB.cols());
		reduced_affine_mat = RB.transpose() * affine_matrix_simplified * RB;
		cout << "reduced_affine_mat.norm() "  << reduced_affine_mat.norm() << endl;

	} */

	return affine_mat_proj;
    }

    Eigen::VectorXd CoupledLinearNS_TT::gen_affine_vec_proj(double current_nu, int current_index)
    {
	Eigen::VectorXd recovered_affine_adv_rhs_proj_xy = Eigen::VectorXd::Zero(RBsize); 
	for (int i = 0; i < RBsize; ++i)
	{
		recovered_affine_adv_rhs_proj_xy -= adv_vec_proj_x[i] * curr_xy_projected(i,0) + adv_vec_proj_y[i] * curr_xy_projected(i,1);
	}	



	Eigen::VectorXd add_rhs_Newton = Eigen::VectorXd::Zero(RBsize); 
	Eigen::VectorXd recovered_affine_adv_rhs_proj_xy_newton = Eigen::VectorXd::Zero(RBsize);
	Eigen::MatrixXd recovered_affine_adv_rhs_proj_xy_newton_RB = Eigen::MatrixXd::Zero(RBsize,RBsize);  
	if (use_Newton)
	{
		// can I build the Newton-required term from recovered_affine_adv_mat_proj_xy and curr_xy_projected ?
		Eigen::MatrixXd recovered_affine_adv_mat_proj_xy = Eigen::MatrixXd::Zero(RBsize, RBsize);

		for (int i = 0; i < RBsize; ++i)
		{
			recovered_affine_adv_mat_proj_xy += adv_mats_proj_x[i] * curr_xy_projected(i,0) + adv_mats_proj_y[i] * curr_xy_projected(i,1);
			recovered_affine_adv_rhs_proj_xy_newton_RB -= adv_vec_proj_x_newton_RB[i] * curr_xy_projected(i,0) + adv_vec_proj_y_newton_RB[i] * curr_xy_projected(i,1);
		}
//		cout << "recovered_affine_adv_mat_proj_xy.rows() " << recovered_affine_adv_mat_proj_xy.rows() << " recovered_affine_adv_mat_proj_xy.cols() " << recovered_affine_adv_mat_proj_xy.cols() << endl;
//		cout << "collect_f_all.rows() " << collect_f_all.rows() << " collect_f_all.cols() " << collect_f_all.cols() << endl;
		Eigen::VectorXd current_f_all = Eigen::VectorXd::Zero(collect_f_all.rows());
		current_f_all = collect_f_all.col(current_index);
		Eigen::VectorXd current_f_all_wo_dbc = remove_rows(current_f_all, elem_loc_dbc);
		Eigen::VectorXd proj_current_f_all_wo_dbc = RB.transpose() * current_f_all_wo_dbc;
		proj_current_f_all_wo_dbc = PODmodes.transpose() * current_f_all;
//		cout << "proj_current_f_all_wo_dbc " << proj_current_f_all_wo_dbc << endl;

		if (use_Newton)
		{
			add_rhs_Newton = recovered_affine_adv_rhs_proj_xy_newton_RB.transpose() * proj_current_f_all_wo_dbc;
//			cout << "add_rhs_Newton " << add_rhs_Newton << endl;
			add_rhs_Newton = recovered_affine_adv_rhs_proj_xy_newton_RB * proj_current_f_all_wo_dbc;
//			cout << "add_rhs_Newton 2 " << add_rhs_Newton << endl;
		}


		for (int i = 0; i < RBsize; ++i)
		{
//			recovered_affine_adv_rhs_proj_xy_newton -= adv_vec_proj_x_newton[i] * curr_xy_projected(i,0) + adv_vec_proj_y_newton[i] * curr_xy_projected(i,1);
//			cout << "adv_vec_proj_x_newton_RB.rows() " << adv_vec_proj_x_newton_RB[i].rows() << " adv_vec_proj_x_newton_RB.cols() " << adv_vec_proj_x_newton_RB[i].cols() << endl;
//			cout << "adv_vec_proj_y_newton_RB.rows() " << adv_vec_proj_y_newton_RB[i].rows() << " adv_vec_proj_y_newton_RB.cols() " << adv_vec_proj_y_newton_RB[i].cols() << endl;
//			recovered_affine_adv_rhs_proj_xy_newton_RB -= adv_vec_proj_x_newton_RB[i] * curr_xy_projected(i,0) + adv_vec_proj_y_newton_RB[i] * curr_xy_projected(i,1);
		}	

//		cout << "recovered_affine_adv_rhs_proj_xy_newton " << -0.5*recovered_affine_adv_rhs_proj_xy_newton << endl;


	}
//	return -the_const_one_rhs_proj - current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy  -0.5*recovered_affine_adv_rhs_proj_xy_newton ;

//	cout << " the_const_one_rhs_proj " <<  the_const_one_rhs_proj  << endl;
//	cout << " the_ABCD_one_rhs_proj " <<  the_ABCD_one_rhs_proj  << endl;
//	cout << " recovered_affine_adv_rhs_proj_xy " <<  recovered_affine_adv_rhs_proj_xy  << endl;

	return -the_const_one_rhs_proj - current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy  -0.5*add_rhs_Newton ;  
    }

    Eigen::VectorXd CoupledLinearNS_TT::gen_affine_vec_proj_2d(double current_nu, double w, int current_index)
    {
	Eigen::VectorXd recovered_affine_adv_rhs_proj_xy = Eigen::VectorXd::Zero(RBsize); 
	Eigen::VectorXd recovered_press_proj = Eigen::VectorXd::Zero(RBsize);
	Eigen::VectorXd recovered_ABCD_proj = Eigen::VectorXd::Zero(RBsize);
	for (int index_elem = 0; index_elem < number_elem_trafo; ++index_elem)
	{
		double detT = Geo_T(w, index_elem, 0);
		double Ta = Geo_T(w, index_elem, 1);
		double Tb = Geo_T(w, index_elem, 2);
		double Tc = Geo_T(w, index_elem, 3);
		double Td = Geo_T(w, index_elem, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td;
		for (int i = 0; i < RBsize; ++i)
		{
	//		recovered_affine_adv_rhs_proj_xy -= adv_vec_proj_x[i] * curr_xy_projected(i,0) + adv_vec_proj_y[i] * curr_xy_projected(i,1);
			recovered_affine_adv_rhs_proj_xy -= detT * curr_xy_projected(i,0) * (Ta * adv_vec_proj_x_2d[i][index_elem][0] + Tc * adv_vec_proj_x_2d[i][index_elem][1]) + detT * curr_xy_projected(i,1) * (Tb * adv_vec_proj_y_2d[i][index_elem][0] + Td * adv_vec_proj_y_2d[i][index_elem][1]);
		}	
		recovered_ABCD_proj += detT * (c00 * the_ABCD_one_rhs_proj_2d[index_elem][0] + c01*(the_ABCD_one_rhs_proj_2d[index_elem][1] + the_ABCD_one_rhs_proj_2d[index_elem][2]) + c11*the_ABCD_one_rhs_proj_2d[index_elem][3]);
		recovered_press_proj += detT * (Ta * the_const_one_rhs_proj_2d[index_elem][0] + Tc * the_const_one_rhs_proj_2d[index_elem][1] + Tb * the_const_one_rhs_proj_2d[index_elem][2] + Td * the_const_one_rhs_proj_2d[index_elem][3]);

	}


	Eigen::VectorXd add_rhs_Newton = Eigen::VectorXd::Zero(RBsize); 

//	cout << " recovered_press_proj " <<  recovered_press_proj  << endl;
//	cout << " recovered_ABCD_proj " <<  recovered_ABCD_proj  << endl;
//	cout << " recovered_affine_adv_rhs_proj_xy " <<  recovered_affine_adv_rhs_proj_xy  << endl;


	return -recovered_press_proj - current_nu * recovered_ABCD_proj + recovered_affine_adv_rhs_proj_xy  -0.5*add_rhs_Newton ;  
//	return -the_const_one_rhs_proj - current_nu * the_ABCD_one_rhs_proj + recovered_affine_adv_rhs_proj_xy  -0.5*add_rhs_Newton ;  
    }


    Eigen::VectorXd CoupledLinearNS_TT::reconstruct_solution_w_dbc(Eigen::VectorXd reprojected_solve)
    {
	Eigen::VectorXd reconstruct_solution = Eigen::VectorXd::Zero(f_bnd_dbc_full_size.rows());  // is of size M_truth_size
	int counter_wo_dbc = 0;
	for (int row_index=0; row_index < f_bnd_dbc_full_size.rows(); ++row_index)
	{
		if (!elem_loc_dbc.count(row_index))
		{
			reconstruct_solution(row_index) = reprojected_solve(counter_wo_dbc);
			counter_wo_dbc++;
		}
		else
		{
			reconstruct_solution(row_index) = f_bnd_dbc_full_size(row_index);
		}
	}
	return reconstruct_solution;
    }

    void CoupledLinearNS_TT::gen_reference_matrices()
    {
	double current_nu = ref_param_nu;
	int current_index = ref_param_index;
	Set_m_kinvis( current_nu );
//	DoInitialiseAdv(snapshot_x_collection[current_index], snapshot_y_collection[current_index]); // why is this necessary? 
//	the_const_one = Get_no_advection_matrix_pressure();
	the_const_one = gen_no_advection_matrix_pressure();
//	cout << "co norm " << the_const_one.block(0, f_bnd_size, 10, 10) << endl;
//	cout << "co2 norm " << the_const_one2.block(0, f_bnd_size, 10, 10) << endl;
//	Eigen::MatrixXd diff = the_const_one - the_const_one2;
//	cout << "diff norm " << diff.block(0, f_bnd_size, 10, 10) << endl;
//	the_ABCD_one = Get_no_advection_matrix_ABCD();
	the_ABCD_one = gen_no_advection_matrix_ABCD();
	the_const_one_simplified = remove_cols_and_rows(the_const_one, elem_loc_dbc);
	the_ABCD_one_simplified = remove_cols_and_rows(the_ABCD_one, elem_loc_dbc);
	the_const_one_proj = RB.transpose() * the_const_one_simplified * RB;
	the_ABCD_one_proj = RB.transpose() * the_ABCD_one_simplified * RB;
	the_ABCD_one_rhs = the_ABCD_one * f_bnd_dbc_full_size;
	the_const_one_rhs = the_const_one * f_bnd_dbc_full_size;
	the_ABCD_one_rhs_simplified = remove_rows(the_ABCD_one_rhs, elem_loc_dbc);
	the_const_one_rhs_simplified = remove_rows(the_const_one_rhs, elem_loc_dbc);
	the_ABCD_one_rhs_proj = RB.transpose() * the_ABCD_one_rhs_simplified;
	the_const_one_rhs_proj = RB.transpose() * the_const_one_rhs_simplified;
//	cout << "the_const_one_rhs_proj " << the_const_one_rhs_proj << endl;
    }

    Eigen::MatrixXd CoupledLinearNS_TT::gen_no_advection_matrix_pressure()
    {
//	double mKinvis = 1;
//	double w = 1;
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 

	// verify and transform to bnd / p / int the snapshot data

        int nz_loc;
        nz_loc = 1;

        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;

	Array<OneD, Eigen::MatrixXd > Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
	Array<OneD, Eigen::MatrixXd > Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
	

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();

//		cout << "ncoeffs " << ncoeffs << endl;
//		cout << "nphys " << nphys << endl;
//		cout << "pqsize " << pqsize << endl;   // pqsize == nphys and ncoeffs == nphys / 2 when?

                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
//		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
//		Array<OneD, double> curr_snap_y_part(nphys, 0.0);

		Array<OneD, double> Dbnd_ele_vec(nsize_p*nsize_bndry, 0.0);
		Array<OneD, double> Dint_ele_vec(nsize_p*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
			            	int psize = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dbnd_ele_vec[ (k*nbmap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			for (int k = 0; k < 2; ++k)
			{
				if (k == 0)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
//						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il] + pcoeffs_y[il];
						Dint_ele_vec[ (k*nimap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		}

		Dbnd_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_bndry);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				Dbnd_elem[curr_elem](i,j) = Dbnd_ele_vec[ i + j*nsize_p];
			}
		} 
		Dint_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_p, nsize_int);
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				Dint_elem[curr_elem](i,j) = Dint_ele_vec[ i + j*nsize_p];
			}
		} 
	}
//	cout << " nsize_p " << nsize_p << endl;
//	cout << " nsize_bndry " << nsize_bndry << endl;
//	cout << "Dbnd_elem[0] 1..4 " << Dbnd_elem[0].block(0,0,3,3) << endl;
	Eigen::MatrixXd Dbnd_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_bndry*nel );
	Eigen::MatrixXd Dint_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_int*nel );
	Eigen::MatrixXd press_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		Dbnd_all.block(i*nsize_p, i*nsize_bndry, nsize_p, nsize_bndry) = Dbnd_elem[i];
		Dint_all.block(i*nsize_p, i*nsize_int, nsize_p, nsize_int) = Dint_elem[i];
	}

	switch(globally_connected) {
		case 0:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -MtM * Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 1:
			press_matrix.block(0, nBndDofs, nBndDofs, f_p_size) = -Mtrafo.transpose() * Dbnd_all.transpose();
			press_matrix.block(nBndDofs, 0, f_p_size, nBndDofs) = -Dbnd_all * Mtrafo;
			press_matrix.block(nBndDofs, nBndDofs + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(nBndDofs + f_p_size, nBndDofs, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 2:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
	}	

	// also need to create appropriate vector right-hand-side contribution

/*	Eigen::VectorXd add_to_rhs_press(M_truth_size); // probably need this for adv and non-adv
	add_to_rhs_press = press_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd press_rhs_add = remove_rows(add_to_rhs_press, elem_loc_dbc);
	press_vec_proj = RB.transpose() * press_rhs_add;
	Eigen::MatrixXd press_matrix_simplified = remove_cols_and_rows(press_matrix, elem_loc_dbc);
	Eigen::MatrixXd press_mat_proj = RB.transpose() * press_matrix_simplified * RB;
*/
	return press_matrix;

    }

    Eigen::MatrixXd CoupledLinearNS_TT::gen_no_advection_matrix_ABCD()
    {
//	double mKinvis = 1;
//	double w = 1;
	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;

	Array<OneD, Eigen::MatrixXd > Ah_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry_p1, nsize_bndry_p1
	Array<OneD, Eigen::MatrixXd > B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Eigen::MatrixXd > D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();

                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
		Array<OneD, double> Ah_ele_vec(Ahrows*Ahrows, 0.0);
		Array<OneD, double> B_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> C_ele_vec(nsize_bndry*nsize_int, 0.0);
		Array<OneD, double> D_ele_vec(nsize_int*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					Ah_ele_vec[ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_0_0[int(bmap[j])] + coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					B_ele_vec[i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_0_0[int(imap[j])] + coeffs_1_1[int(imap[j])];
				}
			} //for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					C_ele_vec[ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_0_0[int(bmap[j])] + coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					D_ele_vec[i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_0_0[int(imap[j])] + coeffs_1_1[int(imap[j])];
				}
			} //for (int k = 0; k < 2; ++k)
		}

 		// chosen choice: redo the nektar++ approach
		// possible alternatives:
		// build instead the sing_* matrices directly
		// or copy for test purposes from setupcoupledmats??

		Ah_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry_p1, nsize_bndry_p1);
		for (int i = 0; i < nsize_bndry_p1; ++i)
		{
			for (int j = 0; j < nsize_bndry_p1; ++j)
			{
				Ah_elem[curr_elem](i,j) = Ah_ele_vec[ i + j*Ahrows ];
			}
		} 
		B_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem](i,j) = B_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		C_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem](i,j) = C_ele_vec[ i + j*nsize_bndry ];
			}
		} 
		D_elem[curr_elem] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem](i,j) = D_ele_vec[ i + j*nsize_int];
			}
		} 



	}

	Eigen::MatrixXd A_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	Eigen::MatrixXd B_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd D_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	Eigen::MatrixXd ABCD_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		A_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = Ah_elem[i].block(0,0,nsize_bndry,nsize_bndry);
		B_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = B_elem[i];
		C_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = C_elem[i];
		D_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = D_elem[i];
	}
	switch(globally_connected) {
		case 0:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 1:
			ABCD_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_all * Mtrafo;
			ABCD_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_all;
			ABCD_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_all.transpose() * Mtrafo;
			ABCD_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 2:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
	}

/*	Eigen::VectorXd add_to_rhs_ABCD(M_truth_size);
	add_to_rhs_ABCD = ABCD_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd ABCD_rhs_add = remove_rows(add_to_rhs_ABCD, elem_loc_dbc);
	ABCD_vec_proj = RB.transpose() * ABCD_rhs_add;
	Eigen::MatrixXd ABCD_matrix_simplified = remove_cols_and_rows(ABCD_matrix, elem_loc_dbc);
	Eigen::MatrixXd ABCD_mat_proj = RB.transpose() * ABCD_matrix_simplified * RB;
*/

	return ABCD_matrix;



    }

    void CoupledLinearNS_TT::gen_reference_matrices_2d()
    {
	// should also loop through the structures, doing an elementwise assembly
	// in principle similar to the advection business		adv_mats_proj_x_2d = Array<OneD, Array<OneD, Array<OneD, Eigen::MatrixXd > > > (RBsize); // should be RBsize x number_elem_trafo x 2 x RBsize x RBsize
	the_const_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo); // should be number_elem_trafo x 4 x RBsize x RBsize
	the_ABCD_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo);
	the_ABCD_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);	
	the_const_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		the_const_one_proj_2d[i] = Array<OneD, Eigen::MatrixXd > (4);
		the_const_one_rhs_proj_2d[i] = Array<OneD, Eigen::VectorXd > (4);
		the_ABCD_one_proj_2d[i] = Array<OneD, Eigen::MatrixXd > (4);
		the_ABCD_one_rhs_proj_2d[i] = Array<OneD, Eigen::VectorXd > (4);
		for (int j = 0; j < 4; ++j)
		{
			the_const_one_proj_2d[i][j] = Eigen::MatrixXd::Zero(RBsize, RBsize);
			the_const_one_rhs_proj_2d[i][j] = Eigen::VectorXd::Zero(RBsize);
			the_ABCD_one_proj_2d[i][j] = Eigen::MatrixXd::Zero(RBsize, RBsize);
			the_ABCD_one_rhs_proj_2d[i][j] = Eigen::VectorXd::Zero(RBsize);
		}
	}	

	StdRegions::StdExpansionSharedPtr locExp;
        Array<OneD, unsigned int> bmap,imap; 
        int nz_loc = 1;
        int nel  = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel   = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
        int Ahrows = nsize_bndry_p1;
	
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  A_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_bndry
//	Array<OneD, Eigen::MatrixXd > Bh_elem(m_fields[0]->GetNumElmts()); // 
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  B_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  C_elem(m_fields[0]->GetNumElmts()); // , nsize_bndry, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  D_elem(m_fields[0]->GetNumElmts()); // , nsize_int, nsize_int
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  Dbnd_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_bndry
	Array<OneD, Array<OneD, Eigen::MatrixXd > >  Dint_elem(m_fields[0]->GetNumElmts()); // , nsize_p, nsize_int
//	Array<OneD, Eigen::MatrixXd > Ch_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_bndry_p1
//	Array<OneD, Eigen::MatrixXd > Dh_elem(m_fields[0]->GetNumElmts()); // nsize_p_m1, nsize_p_m1

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{

//		cout << "curr_elem d " << curr_elem << endl;

		A_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		B_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		C_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		D_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		Dbnd_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		Dint_elem[curr_elem] = Array<OneD, Eigen::MatrixXd > (4);
		for (int i = 0; i < 4; i++)
		{
			A_elem[curr_elem][i] = Eigen::MatrixXd::Zero(Ahrows-1, Ahrows-1);
			B_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
			C_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_bndry, nsize_int);
			D_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_int, nsize_int);
			Dbnd_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_p, nsize_bndry);
			Dint_elem[curr_elem][i] = Eigen::MatrixXd::Zero(nsize_p, nsize_int);
		}
	}

	for (int curr_elem = 0; curr_elem < m_fields[0]->GetNumElmts(); ++curr_elem)
	{

//		cout << "curr_elem " << curr_elem << endl;

/*		int curr_elem_pos = get_curr_elem_pos(curr_elem);
		double detT = Geo_T(w, curr_elem_pos, 0);
		double Ta = Geo_T(w, curr_elem_pos, 1);
		double Tb = Geo_T(w, curr_elem_pos, 2);
		double Tc = Geo_T(w, curr_elem_pos, 3);
		double Td = Geo_T(w, curr_elem_pos, 4);
		double c00 = Ta*Ta + Tb*Tb;
		double c01 = Ta*Tc + Tb*Td;
		double c11 = Tc*Tc + Td*Td; */
                locExp = m_fields[m_velocity[0]]->GetExp(curr_elem);
                locExp->GetBoundaryMap(bmap);
                locExp->GetInteriorMap(imap);
                int ncoeffs = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetNcoeffs();
                int nphys   = m_fields[m_velocity[0]]->GetExp(curr_elem)->GetTotPoints();
		int pqsize  = m_pressure->GetExp(curr_elem)->GetTotPoints();
                int nbmap = bmap.num_elements();
                int nimap = imap.num_elements();
/*		Array<OneD, double> curr_snap_x_part(nphys, 0.0);
		Array<OneD, double> curr_snap_y_part(nphys, 0.0);
		for (int i = 0; i < nphys; ++i)
		{
			curr_snap_x_part[i] = snapshot_x[curr_elem*nphys + i];
			curr_snap_y_part[i] = snapshot_y[curr_elem*nphys + i];
		} */

		Array<OneD, Array<OneD, double> > Ah_ele_vec(4);
		Ah_ele_vec[0] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[1] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[2] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Ah_ele_vec[3] = Array<OneD, double> (Ahrows*Ahrows, 0.0);
		Array<OneD, Array<OneD, double> > B_ele_vec(4);
		B_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[2] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		B_ele_vec[3] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > C_ele_vec(4);
		C_ele_vec[0] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[1] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[2] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		C_ele_vec[3] = Array<OneD, double> (nsize_bndry*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > D_ele_vec(4);
		D_ele_vec[0] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		D_ele_vec[1] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		D_ele_vec[2] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		D_ele_vec[3] = Array<OneD, double>(nsize_int*nsize_int, 0.0);
		Array<OneD, Array<OneD, double> > Dbnd_ele_vec(4);
		Dbnd_ele_vec[0] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Dbnd_ele_vec[1] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Dbnd_ele_vec[2] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Dbnd_ele_vec[3] = Array<OneD, double>(nsize_p*nsize_bndry, 0.0);
		Array<OneD, Array<OneD, double> > Dint_ele_vec(4);
		Dint_ele_vec[0] = Array<OneD, double> (nsize_p*nsize_int, 0.0);
		Dint_ele_vec[1] = Array<OneD, double> (nsize_p*nsize_int, 0.0);
		Dint_ele_vec[2] = Array<OneD, double> (nsize_p*nsize_int, 0.0);
		Dint_ele_vec[3] = Array<OneD, double> (nsize_p*nsize_int, 0.0);

		for (int i = 0; i < nbmap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);

			coeffs[bmap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);

			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);

			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					Ah_ele_vec[0][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_0_0[int(bmap[j])];
					Ah_ele_vec[1][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_0_1[int(bmap[j])];
					Ah_ele_vec[2][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_1_0[int(bmap[j])];
					Ah_ele_vec[3][ i+k*nbmap + (j+k*nbmap)*Ahrows ] += coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					B_ele_vec[0][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_0_0[int(imap[j])];
					B_ele_vec[1][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_0_1[int(imap[j])];
					B_ele_vec[2][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_1_0[int(imap[j])];
					B_ele_vec[3][ i+k*nbmap + (j+k*nimap)*nsize_bndry] += coeffs_1_1[int(imap[j])];
				}
				if (k == 0)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[0][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dbnd_ele_vec[1][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dbnd_ele_vec[2][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dbnd_ele_vec[3][ (k*nbmap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		} // for (int i = 0; i < nbmap; ++i)
		for (int i = 0; i < nimap; ++i)
		{
			Array<OneD, double> coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_x_coeffs(ncoeffs, 0.0);
			Array<OneD, double> adv_y_coeffs(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_0_1(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_0(ncoeffs, 0.0);
			Array<OneD, double> coeffs_1_1(ncoeffs, 0.0);
			Array<OneD, double> phys(nphys, 0.0);
			Array<OneD, double> deriv_0(pqsize, 0.0);
			Array<OneD, double> deriv_1(pqsize, 0.0);
			coeffs[imap[i]] = 1.0;
			m_fields[m_velocity[0]]->GetExp(curr_elem)->BwdTrans(coeffs,phys);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[0], phys, deriv_0);
			locExp->PhysDeriv(MultiRegions::DirCartesianMap[1], phys, deriv_1);
			locExp->IProductWRTDerivBase(0, deriv_0, coeffs_0_0);
			locExp->IProductWRTDerivBase(1, deriv_0, coeffs_0_1);
			locExp->IProductWRTDerivBase(0, deriv_1, coeffs_1_0);
			locExp->IProductWRTDerivBase(1, deriv_1, coeffs_1_1);
			for (int k = 0; k < 2; ++k)
			{
				for (int j = 0; j < nbmap; ++j)
				{
					C_ele_vec[0][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_0_0[int(bmap[j])];
					C_ele_vec[1][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_0_1[int(bmap[j])];
					C_ele_vec[2][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_1_0[int(bmap[j])];
					C_ele_vec[3][ j+k*nbmap + (i+k*nimap)*nsize_bndry ] += coeffs_1_1[int(bmap[j])];
				}
				for (int j = 0; j < nimap; ++j)
				{
					D_ele_vec[0][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_0_0[int(imap[j])];
					D_ele_vec[1][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_0_1[int(imap[j])];
					D_ele_vec[2][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_1_0[int(imap[j])];
					D_ele_vec[3][i+k*nimap + (j+k*nimap)*nsize_int] += coeffs_1_1[int(imap[j])];
				}
				if (k == 0)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[0][ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dint_ele_vec[1][ (k*nimap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
				if (k == 1)
				{
			            	int psize   = m_pressure->GetExp(curr_elem)->GetNcoeffs();
			               	Array<OneD, NekDouble> pcoeffs_x(psize);
			               	Array<OneD, NekDouble> pcoeffs_y(psize);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_0,pcoeffs_x);
					m_pressure->GetExp(curr_elem)->IProductWRTBase(deriv_1,pcoeffs_y);
					for (int il = 0; il < nsize_p; ++il)
					{
						Dint_ele_vec[2][ (k*nimap + i)*nsize_p + il ] = pcoeffs_x[il];
						Dint_ele_vec[3][ (k*nimap + i)*nsize_p + il ] = pcoeffs_y[il];
					}
				}
			} //for (int k = 0; k < 2; ++k)
		}

//		cout << "starting writing vec -> mat" << endl;

		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				A_elem[curr_elem][0](i,j) = Ah_ele_vec[0][ i + j*Ahrows ];
				A_elem[curr_elem][1](i,j) = Ah_ele_vec[1][ i + j*Ahrows ];
				A_elem[curr_elem][2](i,j) = Ah_ele_vec[2][ i + j*Ahrows ];
				A_elem[curr_elem][3](i,j) = Ah_ele_vec[3][ i + j*Ahrows ];
			}
		} 
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				B_elem[curr_elem][0](i,j) = B_ele_vec[0][ i + j*nsize_bndry ];
				B_elem[curr_elem][1](i,j) = B_ele_vec[1][ i + j*nsize_bndry ];
				B_elem[curr_elem][2](i,j) = B_ele_vec[2][ i + j*nsize_bndry ];
				B_elem[curr_elem][3](i,j) = B_ele_vec[3][ i + j*nsize_bndry ];
			}
		} 
		for (int i = 0; i < nsize_bndry; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				C_elem[curr_elem][0](i,j) = C_ele_vec[0][ i + j*nsize_bndry ];
				C_elem[curr_elem][1](i,j) = C_ele_vec[1][ i + j*nsize_bndry ];
				C_elem[curr_elem][2](i,j) = C_ele_vec[2][ i + j*nsize_bndry ];
				C_elem[curr_elem][3](i,j) = C_ele_vec[3][ i + j*nsize_bndry ];
			}
		} 
		for (int i = 0; i < nsize_int; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				D_elem[curr_elem][0](i,j) = D_ele_vec[0][ i + j*nsize_int];
				D_elem[curr_elem][1](i,j) = D_ele_vec[1][ i + j*nsize_int];
				D_elem[curr_elem][2](i,j) = D_ele_vec[2][ i + j*nsize_int];
				D_elem[curr_elem][3](i,j) = D_ele_vec[3][ i + j*nsize_int];
			}
		} 
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_bndry; ++j)
			{
				Dbnd_elem[curr_elem][0](i,j) = Dbnd_ele_vec[0][ i + j*nsize_p];
				Dbnd_elem[curr_elem][1](i,j) = Dbnd_ele_vec[1][ i + j*nsize_p];
				Dbnd_elem[curr_elem][2](i,j) = Dbnd_ele_vec[2][ i + j*nsize_p];
				Dbnd_elem[curr_elem][3](i,j) = Dbnd_ele_vec[3][ i + j*nsize_p];
			}
		} 
		for (int i = 0; i < nsize_p; ++i)
		{
			for (int j = 0; j < nsize_int; ++j)
			{
				Dint_elem[curr_elem][0](i,j) = Dint_ele_vec[0][ i + j*nsize_p];
				Dint_elem[curr_elem][1](i,j) = Dint_ele_vec[1][ i + j*nsize_p];
				Dint_elem[curr_elem][2](i,j) = Dint_ele_vec[2][ i + j*nsize_p];
				Dint_elem[curr_elem][3](i,j) = Dint_ele_vec[3][ i + j*nsize_p];
			}
		}

	} // loop over curr_elem

//	to be build in a seperate projector function:
//	the_const_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo); // should be number_elem_trafo x 4 x RBsize x RBsize
//	the_ABCD_one_proj_2d = Array<OneD, Array<OneD, Eigen::MatrixXd > > (number_elem_trafo);
//	the_ABCD_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);	
//	the_const_one_rhs_proj_2d = Array<OneD, Array<OneD, Eigen::VectorXd > > (number_elem_trafo);
	
	for (int i = 0; i < number_elem_trafo; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			the_const_one_proj_2d[i][j] = press_geo_mat_projector(Dbnd_elem, Dint_elem, i, j, the_const_one_rhs_proj_2d[i][j]);
//			cout << "finished the_const_one_proj_2d " << endl;
			the_ABCD_one_proj_2d[i][j] = ABCD_geo_mat_projector(A_elem, B_elem, C_elem, D_elem, i, j, the_ABCD_one_rhs_proj_2d[i][j]);
		}
	}
	


    }

    Eigen::MatrixXd CoupledLinearNS_TT::ABCD_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > A_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > B_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > C_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > D_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &ABCD_vec_proj)
    {
	Eigen::MatrixXd ABCD_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
        int nz_loc = 1;
        int nel = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
	Eigen::MatrixXd A_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_bndry*nel );
	Eigen::MatrixXd B_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero( nsize_bndry*nel , nsize_int*nel );
	Eigen::MatrixXd D_all = Eigen::MatrixXd::Zero( nsize_int*nel , nsize_int*nel );
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		int curr_elem_pos = get_curr_elem_pos(i);
		if (curr_elem_pos == curr_elem_trafo)
		{
			Eigen::MatrixXd curr_A_elem = A_elem[i][deriv_index];
			Eigen::MatrixXd curr_B_elem = B_elem[i][deriv_index];
			Eigen::MatrixXd curr_C_elem = C_elem[i][deriv_index];
			Eigen::MatrixXd curr_D_elem = D_elem[i][deriv_index];
			A_all.block(i*nsize_bndry, i*nsize_bndry, nsize_bndry, nsize_bndry) = curr_A_elem;
			B_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_B_elem;
			C_all.block(i*nsize_bndry, i*nsize_int, nsize_bndry, nsize_int) = curr_C_elem;
			D_all.block(i*nsize_int, i*nsize_int, nsize_int, nsize_int) = curr_D_elem;
		}
	}
/*	switch(globally_connected) {
		case 0:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -MtM * Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 1:
			press_matrix.block(0, nBndDofs, nBndDofs, f_p_size) = -Mtrafo.transpose() * Dbnd_all.transpose();
			press_matrix.block(nBndDofs, 0, f_p_size, nBndDofs) = -Dbnd_all * Mtrafo;
			press_matrix.block(nBndDofs, nBndDofs + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(nBndDofs + f_p_size, nBndDofs, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 2:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
	}	*/
	switch(globally_connected) {
		case 0:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = MtM * A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = MtM * B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 1:
			ABCD_matrix.block(0, 0, nBndDofs, nBndDofs) = Mtrafo.transpose() * A_all * Mtrafo;
			ABCD_matrix.block(0, nBndDofs + f_p_size, nBndDofs, f_int_size) = Mtrafo.transpose() * B_all;
			ABCD_matrix.block(nBndDofs + f_p_size, 0, f_int_size, nBndDofs) = C_all.transpose() * Mtrafo;
			ABCD_matrix.block(nBndDofs + f_p_size, nBndDofs + f_p_size, f_int_size, f_int_size) = D_all;
			break;
		case 2:
			ABCD_matrix.block(0, 0, f_bnd_size, f_bnd_size) = A_all;
			ABCD_matrix.block(0, f_bnd_size + f_p_size, f_bnd_size, f_int_size) = B_all;
			ABCD_matrix.block(f_bnd_size + f_p_size, 0, f_int_size, f_bnd_size) = C_all.transpose();
			ABCD_matrix.block(f_bnd_size + f_p_size, f_bnd_size + f_p_size, f_int_size, f_int_size) = D_all;
			break;
	}
	Eigen::VectorXd add_to_rhs_ABCD(M_truth_size);
	add_to_rhs_ABCD = ABCD_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd ABCD_rhs_add = remove_rows(add_to_rhs_ABCD, elem_loc_dbc);
	ABCD_vec_proj = RB.transpose() * ABCD_rhs_add;
	Eigen::MatrixXd ABCD_matrix_simplified = remove_cols_and_rows(ABCD_matrix, elem_loc_dbc);
	Eigen::MatrixXd ABCD_mat_proj = RB.transpose() * ABCD_matrix_simplified * RB;
	return ABCD_mat_proj;
    }
    
    Eigen::MatrixXd CoupledLinearNS_TT::press_geo_mat_projector(Array<OneD, Array<OneD, Eigen::MatrixXd > > Dbnd_elem, Array<OneD, Array<OneD, Eigen::MatrixXd > > Dint_elem, int curr_elem_trafo, int deriv_index, Eigen::VectorXd &press_vec_proj)
    {

	// this function (i) expands to full size, (ii) removes dbc cols and rows, (iii) projects

//	cout << "entering press_geo_mat_projector" << endl;

	Eigen::MatrixXd press_matrix = Eigen::MatrixXd::Zero(M_truth_size, M_truth_size);
        int nz_loc = 1;
        int nel = m_fields[m_velocity[0]]->GetNumElmts();
        int nvel = m_velocity.num_elements();
        int nsize_bndry = nvel*m_fields[m_velocity[0]]->GetExp(0)->NumBndryCoeffs()*nz_loc;
        int nsize_bndry_p1 = nsize_bndry+nz_loc;
        int nsize_int = (nvel*m_fields[m_velocity[0]]->GetExp(0)->GetNcoeffs()*nz_loc - nsize_bndry);
        int nsize_p = m_pressure->GetExp(0)->GetNcoeffs()*nz_loc;
        int nsize_p_m1 = nsize_p-nz_loc;
	Eigen::MatrixXd Dbnd_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_bndry*nel );
	Eigen::MatrixXd Dint_all = Eigen::MatrixXd::Zero( nsize_p*nel , nsize_int*nel );

	time_t timer_1;
	time_t timer_2;
	time_t timer_3;
	time_t timer_4;

	time(&timer_1);  /* get current time; same as: timer = time(NULL)  */
	for (int i = 0; i < m_fields[0]->GetNumElmts(); ++i)
	{
		int curr_elem_pos = get_curr_elem_pos(i);
		if (curr_elem_pos == curr_elem_trafo)
		{
			Eigen::MatrixXd curr_Dbnd_elem = Dbnd_elem[i][deriv_index];
			Eigen::MatrixXd curr_Dint_elem = Dint_elem[i][deriv_index];

			Dbnd_all.block(i*nsize_p, i*nsize_bndry, nsize_p, nsize_bndry) = curr_Dbnd_elem;
			Dint_all.block(i*nsize_p, i*nsize_int, nsize_p, nsize_int) = curr_Dint_elem;
	
		}
	}

	// better to use:
//	int f_bnd_size;
//	int f_p_size;
//	int f_int_size;

	switch(globally_connected) {
		case 0:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -MtM * Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 1:
			press_matrix.block(0, nBndDofs, nBndDofs, f_p_size) = -Mtrafo.transpose() * Dbnd_all.transpose();
			press_matrix.block(nBndDofs, 0, f_p_size, nBndDofs) = -Dbnd_all * Mtrafo;
			press_matrix.block(nBndDofs, nBndDofs + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(nBndDofs + f_p_size, nBndDofs, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
		case 2:
			press_matrix.block(0, f_bnd_size, f_bnd_size, f_p_size) = -Dbnd_all.transpose();
			press_matrix.block(f_bnd_size, 0, f_p_size, f_bnd_size) = -Dbnd_all;
			press_matrix.block(f_bnd_size, f_bnd_size + f_p_size, f_p_size, f_int_size) = -Dint_all;	
			press_matrix.block(f_bnd_size + f_p_size, f_bnd_size, f_int_size, f_p_size) = -Dint_all.transpose();
			break;
	}	

	// also need to create appropriate vector right-hand-side contribution

	Eigen::VectorXd add_to_rhs_press(M_truth_size); // probably need this for adv and non-adv
	add_to_rhs_press = press_matrix * f_bnd_dbc_full_size;
	Eigen::VectorXd press_rhs_add = remove_rows(add_to_rhs_press, elem_loc_dbc);
	press_vec_proj = RB.transpose() * press_rhs_add;
	Eigen::MatrixXd press_matrix_simplified = remove_cols_and_rows(press_matrix, elem_loc_dbc);
	Eigen::MatrixXd press_mat_proj = RB.transpose() * press_matrix_simplified * RB;

	return press_mat_proj;
    }

    void CoupledLinearNS_TT::compute_snapshots(int number_of_snapshots)
    {
	CoupledLinearNS_trafoP babyCLNS_trafo(m_session);
	babyCLNS_trafo.InitObject();
	babyCLNS_trafo.use_Newton = use_Newton;
	babyCLNS_trafo.snapshot_computation_plot_rel_errors = snapshot_computation_plot_rel_errors;
	Array<OneD, NekDouble> zero_phys_init(GetNpoints(), 0.0);
	snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
	snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
        for(int i = 0; i < number_of_snapshots; ++i)
	{

		Array<OneD, Array<OneD, NekDouble> > converged_solution = babyCLNS_trafo.DoSolve_at_param(zero_phys_init, zero_phys_init, param_vector[i]);
		snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		for (int j=0; j < GetNpoints(); ++j)
		{
			snapshot_x_collection[i][j] = converged_solution[0][j];
			snapshot_y_collection[i][j] = converged_solution[1][j];
		}
	}


    }

    void CoupledLinearNS_TT::load_snapshots_geometry_params(int number_of_snapshots)
    {

		// assuming it is prepared in the correct ordering - i.e. - outer loop dir0, inner loop dir1

		int nvelo = 2;
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot(nvelo); // for a 2D problem

        for(int i = 0; i < nvelo; ++i)
        {
            test_load_snapshot[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);  // number of phys points
        }

		cout << "a size of double  " << sizeof(double) << endl;	
		cout << "no of phys points " << GetNpoints() << endl;
		cout << "memory requirements for snapshots in bytes " << sizeof(double)*number_of_snapshots*2*GetNpoints() << endl;

        std::vector<std::string> fieldStr;
        for(int i = 0; i < nvelo; ++i)
        {
           fieldStr.push_back(m_session->GetVariable(i));
        }

		// better way - difficult:
        for(int i = 0; i < number_of_snapshots; ++i)
        {
		// generate the correct string
//		std::stringstream sstm;
//		sstm << "TestSnap" << i+1;
//		std::string result = sstm.str();
//		const char* rr = result.c_str();

//	        EvaluateFunction(fieldStr, test_load_snapshot, result);

//	        LibUtilities::FunctionType vType;
//        	vType = m_session->GetFunctionType(result, fieldStr[0], 0);
/*		cout << " checking if (vType == LibUtilities::eFunctionTypeExpression)  " << endl;
		if (vType == LibUtilities::eFunctionTypeExpression)
			cout << " check if (vType == LibUtilities::eFunctionTypeExpression) successful " << endl;
		cout << " checking if (vType == LibUtilities::eFunctionTypeFile)  " << endl;
		if (vType == LibUtilities::eFunctionTypeFile)
			cout << " check if (vType == LibUtilities::eFunctionTypeFile) successful " << endl;         <-- this usually
		cout << " checking if (vType == LibUtilities::eFunctionTypeTransientFile)  " << endl;
		if (vType == LibUtilities::eFunctionTypeTransientFile)
			cout << " check if (vType == LibUtilities::eFunctionTypeTransientFile) successful " << endl; */
		// std::string filename = m_session->GetFunctionFilename(result, fieldStr[0], 0);
		// cout << " determined filename " << filename << endl;
		// auto gen these
/*		std::stringstream sstm;
		sstm << "xml_channel_narrowROM_" << i << ".fld";
		std::string result = sstm.str();
		const char* rr = result.c_str();

		EvaluateFunctionFld(fieldStr[0], test_load_snapshot[0], pFunctionName, 0.0, 0);
*/
		}


		snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
		snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
        for(int i = 0; i < number_of_snapshots; ++i)
        {
			// generate the correct string
			std::stringstream sstm;
			sstm << "TestSnap" << i+1;
			std::string result = sstm.str();
			const char* rr = result.c_str();

	        EvaluateFunction(fieldStr, test_load_snapshot, result);
			snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
			snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
			for (int j=0; j < GetNpoints(); ++j)
			{
				snapshot_x_collection[i][j] = test_load_snapshot[0][j];
				snapshot_y_collection[i][j] = test_load_snapshot[1][j];
			}
        }

    }

    void CoupledLinearNS_TT::load_snapshots_geometry_params_conv_Oseen(int number_of_snapshots)
    {

		// assuming it is prepared in the correct ordering - i.e. - outer loop dir0, inner loop dir1

		int nvelo = 2;
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot(nvelo); // for a 2D problem

        for(int i = 0; i < nvelo; ++i)
        {
            test_load_snapshot[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);  // number of phys points
        }

		cout << "no of phys points " << GetNpoints() << endl;

        std::vector<std::string> fieldStr;
        for(int i = 0; i < nvelo; ++i)
        {
           fieldStr.push_back(m_session->GetVariable(i));
        }

		snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
		snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
		int unique_number_of_snapshots;
		if (use_non_unique_up_to_two)
		{
			unique_number_of_snapshots = number_of_snapshots / 2;
		}
		else
		{
			unique_number_of_snapshots = number_of_snapshots;
		}
        for(int i = 0; i < unique_number_of_snapshots; ++i)
        {
		// generate the correct string
		std::stringstream sstm;
		sstm << "TestSnap_cO_" << i+1;
		std::string result = sstm.str();
		const char* rr = result.c_str();

	        EvaluateFunction(fieldStr, test_load_snapshot, result);
		snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		for (int j=0; j < GetNpoints(); ++j)
		{
			snapshot_x_collection[i][j] = test_load_snapshot[0][j];
			snapshot_y_collection[i][j] = test_load_snapshot[1][j];
		}
        }
		if (use_non_unique_up_to_two)
		{
	        for(int i = unique_number_of_snapshots; i < 2*unique_number_of_snapshots; ++i)
	        {
				// generate the correct string
				std::stringstream sstm;
				sstm << "TestSnap_cO_" << i+1;
				std::string result = sstm.str();
				const char* rr = result.c_str();

		        EvaluateFunction(fieldStr, test_load_snapshot, result);
				snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
				snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
				for (int j=0; j < GetNpoints(); ++j)
				{
					snapshot_x_collection[i][j] = test_load_snapshot[0][j];
					snapshot_y_collection[i][j] = test_load_snapshot[1][j];
				}	
       		}
		}

    }

    void CoupledLinearNS_TT::load_snapshots(int number_of_snapshots)
    {
	// fill the fields snapshot_x_collection and snapshot_y_collection

	int nvelo = 2;
        Array<OneD, Array<OneD, NekDouble> > test_load_snapshot(nvelo); // for a 2D problem

        for(int i = 0; i < nvelo; ++i)
        {
            test_load_snapshot[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);  // number of phys points
        }
               
        std::vector<std::string> fieldStr;
        for(int i = 0; i < nvelo; ++i)
        {
           fieldStr.push_back(m_session->GetVariable(i));
        }

	snapshot_x_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);
	snapshot_y_collection = Array<OneD, Array<OneD, NekDouble> > (number_of_snapshots);

        for(int i = 0; i < number_of_snapshots; ++i)
        {
		// generate the correct string
		std::stringstream sstm;
		sstm << "TestSnap" << i+1;
		std::string result = sstm.str();
		const char* rr = result.c_str();

	        EvaluateFunction(fieldStr, test_load_snapshot, result);
		snapshot_x_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		snapshot_y_collection[i] = Array<OneD, NekDouble> (GetNpoints(), 0.0);
		for (int j=0; j < GetNpoints(); ++j)
		{
			snapshot_x_collection[i][j] = test_load_snapshot[0][j];
			snapshot_y_collection[i][j] = test_load_snapshot[1][j];
		}
        }
    }

    void CoupledLinearNS_TT::trafoSnapshot_simple(Eigen::MatrixXd RB_in)
    {
	// moved to CoupledLinearNS_trafoP
    }

    double CoupledLinearNS_TT::L2norm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y )
    {
	// the input comes in phys coords
        NekDouble L2norm = -1.0;
//	cout << "m_NumQuadPointsError " << m_NumQuadPointsError << endl; // should be 0
//	cout << "component_x.num_elements() " << component_x.num_elements() << endl; // should be nphys
        if (m_NumQuadPointsError == 0)
        {
		double L2norm_x = m_fields[0]->L2(component_x);
		double L2norm_y = m_fields[1]->L2(component_y);
		L2norm = sqrt( L2norm_x*L2norm_x + L2norm_y*L2norm_y );
        }
	return L2norm;
    }

    double CoupledLinearNS_TT::Linfnorm_ITHACA( Array< OneD, NekDouble > component_x, Array< OneD, NekDouble > component_y )
    {
	// the input comes in phys coords
        NekDouble Linfnorm = -1.0;
//	cout << "m_NumQuadPointsError " << m_NumQuadPointsError << endl; // should be 0
//	cout << "component_x.num_elements() " << component_x.num_elements() << endl; // should be nphys
        if (m_NumQuadPointsError == 0)
        {
		double Linfnorm_x = m_fields[0]->L2(component_x);
		double Linfnorm_y = m_fields[1]->L2(component_y);
		Linfnorm = max(Linfnorm_x, Linfnorm_y);
        }
	return Linfnorm;
    }

   void CoupledLinearNS_TT::DefineRBspace(Eigen::MatrixXd RB_in)
    {
	// unused
	
    }

    void CoupledLinearNS_TT::v_Output(void)
    {    
        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.num_elements()+1);
        std::vector<std::string> variables(m_fields.num_elements()+1);
        int i;
        
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
        }


        fieldcoeffs[i] = Array<OneD, NekDouble>(m_fields[0]->GetNcoeffs());  
        // project pressure field to velocity space        
        if(m_singleMode==true)
        {
            Array<OneD, NekDouble > tmpfieldcoeffs (m_fields[0]->GetNcoeffs()/2);
            m_pressure->GetPlane(0)->BwdTrans_IterPerExp(m_pressure->GetPlane(0)->GetCoeffs(), m_pressure->GetPlane(0)->UpdatePhys());
            m_pressure->GetPlane(1)->BwdTrans_IterPerExp(m_pressure->GetPlane(1)->GetCoeffs(), m_pressure->GetPlane(1)->UpdatePhys()); 
            m_fields[0]->GetPlane(0)->FwdTrans_IterPerExp(m_pressure->GetPlane(0)->GetPhys(),fieldcoeffs[i]);
            m_fields[0]->GetPlane(1)->FwdTrans_IterPerExp(m_pressure->GetPlane(1)->GetPhys(),tmpfieldcoeffs);
            for(int e=0; e<m_fields[0]->GetNcoeffs()/2; e++)
            {
                fieldcoeffs[i][e+m_fields[0]->GetNcoeffs()/2] = tmpfieldcoeffs[e];
            }          
        }
        else
        {
            m_pressure->BwdTrans_IterPerExp(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());
            m_fields[0]->FwdTrans_IterPerExp(m_pressure->GetPhys(),fieldcoeffs[i]);
        }
        variables[i] = "p"; 
        
        std::string outname = m_sessionName + ".fld";
        
        WriteFld(outname,m_fields[0],fieldcoeffs,variables);

/*
//        Array<OneD, NekDouble> glo_fieldcoeffs(m_fields[0]->GetNcoeffs(), 1234.5678);
        Array<OneD, NekDouble> glo_fieldcoeffs(fieldcoeffs[2].num_elements(), 1234.5678);
	m_fields[0]->LocalToGlobal(fieldcoeffs[2], glo_fieldcoeffs);
//	m_pressure->LocalToGlobal(fieldcoeffs[2], glo_fieldcoeffs); This method is not defined or valid for this class type


	int first_ngc_index = 0;

	while(glo_fieldcoeffs[first_ngc_index] > 1234.5679 || glo_fieldcoeffs[first_ngc_index] < 1234.5677)
	{
		first_ngc_index++;
	}

	// can construct a glodofphys here by going with FwdTrans_IterPerExp and LocalToGlobal

        Array<OneD, Array<OneD, NekDouble> > glo_coll_fieldcoeffs(m_fields[0]->GetNpoints());
	for(int counter_phys_index = 0; counter_phys_index < m_fields[0]->GetNpoints(); counter_phys_index++)
	{
	        Array<OneD, NekDouble> phys_fieldcoeffs(m_fields[0]->GetNpoints(), 0.0);
	        Array<OneD, NekDouble> loc_fieldcoeffs(m_fields[0]->GetNcoeffs(), 0.0);
	        Array<OneD, NekDouble> glo_fieldcoeffs(first_ngc_index, 0.0);
		phys_fieldcoeffs[counter_phys_index] = 1.0;
                m_fields[0]->FwdTrans_IterPerExp(phys_fieldcoeffs, loc_fieldcoeffs);
		m_fields[0]->LocalToGlobal(loc_fieldcoeffs, glo_fieldcoeffs);
//		cout << "cgi: " << counter_global_index << " value " << glo_fieldcoeffs[counter_global_index] << endl;
	}

          ofstream myfilegdp ("phys_glo_dof.txt");
	  if (myfilegdp.is_open())
	  {
		for (int i_global_dofs1 = 0; i_global_dofs1 < first_ngc_index; i_global_dofs1++)
		{
			for (int i_phys_dofs2 = 0; i_phys_dofs2 < m_fields[0]->GetNpoints(); i_phys_dofs2++)
			{
				myfilegdp << std::setprecision(17) << glo_coll_fieldcoeffs[i_phys_dofs2][i_global_dofs1] << "\t";
			}
			myfilegdp << "\n";
		}
                myfilegdp.close();
	  }
	  else cout << "Unable to open file"; 

        Array<OneD, Array<OneD, NekDouble> > glo_coll_fieldcoeffs_gdp(m_fields[0]->GetNpoints());
	for(int counter_glo_index = 0; counter_glo_index < first_ngc_index; counter_glo_index++)
	{
	        Array<OneD, NekDouble> phys_fieldcoeffs(m_fields[0]->GetNpoints(), 0.0);
	        Array<OneD, NekDouble> loc_fieldcoeffs(m_fields[0]->GetNcoeffs(), 0.0);
	        Array<OneD, NekDouble> glo_fieldcoeffs(first_ngc_index, 0.0);
		glo_fieldcoeffs[counter_glo_index] = 1.0;
                m_fields[0]->GlobalToLocal(glo_fieldcoeffs, loc_fieldcoeffs);
		m_fields[0]->BwdTrans_IterPerExp(loc_fieldcoeffs, phys_fieldcoeffs);
		
		glo_coll_fieldcoeffs_gdp[counter_glo_index] = phys_fieldcoeffs;
//		cout << "cgi: " << counter_global_index << " value " << glo_fieldcoeffs[counter_global_index] << endl;
	}

          ofstream myfilegdpr ("glo_dof_phys.txt");
	  if (myfilegdpr.is_open())
	  {
		for (int i_global_dofs1 = 0; i_global_dofs1 < first_ngc_index; i_global_dofs1++)
		{
			for (int i_phys_dofs2 = 0; i_phys_dofs2 < m_fields[0]->GetNpoints(); i_phys_dofs2++)
			{
				myfilegdpr << std::setprecision(17) << glo_coll_fieldcoeffs_gdp[i_global_dofs1][i_phys_dofs2] << "\t";
			}
			myfilegdpr << "\n";
		}
                myfilegdpr.close();
	  }
	  else cout << "Unable to open file"; 

	// compute the dof correction

        Array<OneD, Array<OneD, NekDouble> > glo_coll_fieldcoeffs_dof_corr(m_fields[0]->GetNpoints());
	for(int counter_glo_index = 0; counter_glo_index < first_ngc_index; counter_glo_index++)
	{
	        Array<OneD, NekDouble> phys_fieldcoeffs(m_fields[0]->GetNpoints(), 0.0);
	        Array<OneD, NekDouble> loc_fieldcoeffs(m_fields[0]->GetNcoeffs(), 0.0);
	        Array<OneD, NekDouble> glo_fieldcoeffs(first_ngc_index, 0.0);
		glo_fieldcoeffs[counter_glo_index] = 1.0;
                m_fields[0]->GlobalToLocal(glo_fieldcoeffs, loc_fieldcoeffs);
		m_fields[0]->BwdTrans(loc_fieldcoeffs, phys_fieldcoeffs);
                m_fields[0]->FwdTrans_IterPerExp(phys_fieldcoeffs, loc_fieldcoeffs);
		m_fields[0]->LocalToGlobal(loc_fieldcoeffs, glo_fieldcoeffs);
		glo_coll_fieldcoeffs_dof_corr[counter_glo_index] = glo_fieldcoeffs;
//		cout << "cgi: " << counter_global_index << " value " << glo_fieldcoeffs[counter_global_index] << endl;
	}

          ofstream myfilegdprdc ("glo_dof_corr.txt");
	  if (myfilegdprdc.is_open())
	  {
		for (int i_global_dofs1 = 0; i_global_dofs1 < first_ngc_index; i_global_dofs1++)
		{
			for (int i_phys_dofs2 = 0; i_phys_dofs2 < first_ngc_index; i_phys_dofs2++)
			{
				myfilegdprdc << std::setprecision(17) << glo_coll_fieldcoeffs_dof_corr[i_global_dofs1][i_phys_dofs2] << "\t";
			}
			myfilegdprdc << "\n";
		}
                myfilegdprdc.close();
	  }
	  else cout << "Unable to open file"; 





	// construct a LocalGloMapA 

        Array<OneD, Array<OneD, NekDouble> > glo_coll_fieldcoeffsA(m_fields[0]->GetNcoeffs());
	for(int counter_loc_index = 0; counter_loc_index < m_fields[0]->GetNcoeffs(); counter_loc_index++)
	{
	        Array<OneD, NekDouble> loc_fieldcoeffs(m_fields[0]->GetNcoeffs(), 0.0);
	        Array<OneD, NekDouble> glo_fieldcoeffs(first_ngc_index, 0.0);
		loc_fieldcoeffs[counter_loc_index] = 1.0;
		m_fields[0]->LocalToGlobal(loc_fieldcoeffs, glo_fieldcoeffs);
		glo_coll_fieldcoeffsA[counter_loc_index] = glo_fieldcoeffs;
//		cout << "cgi: " << counter_global_index << " value " << glo_fieldcoeffs[counter_global_index] << endl;
	}

          ofstream myfilegdpA ("LocGloMapMatA.txt");
	  if (myfilegdpA.is_open())
	  {
		for (int i_global_dofs1 = 0; i_global_dofs1 < first_ngc_index; i_global_dofs1++)
		{
			for (int i_phys_dofs2 = 0; i_phys_dofs2 < m_fields[0]->GetNcoeffs(); i_phys_dofs2++)
			{
				myfilegdpA << std::setprecision(17) << glo_coll_fieldcoeffsA[i_phys_dofs2][i_global_dofs1] << "\t";
			}
			myfilegdpA << "\n";
		}
                myfilegdpA.close();
	  }
	  else cout << "Unable to open file"; 

        std::string outname_txt = m_sessionName + ".txt";
	const char* outname_t = outname_txt.c_str();

        Array<OneD, NekDouble> glo_fieldcoeffs_x(first_ngc_index, 0.0);
        Array<OneD, NekDouble> glo_fieldcoeffs_y(first_ngc_index, 0.0);
	m_fields[0]->LocalToGlobal(m_fields[0]->GetCoeffs(), glo_fieldcoeffs_x);
	m_fields[1]->LocalToGlobal(m_fields[1]->GetCoeffs(), glo_fieldcoeffs_y);


          ofstream myfile_t (outname_t);
	  if (myfile_t.is_open())
	  {
		for (int i_glo_dofs = 0; i_glo_dofs < first_ngc_index; i_glo_dofs++)
		{
			myfile_t << std::setprecision(17) << glo_fieldcoeffs_x[i_glo_dofs] << "\n";
		}
		for (int i_glo_dofs = 0; i_glo_dofs < first_ngc_index; i_glo_dofs++)
		{
			myfile_t << std::setprecision(17) << glo_fieldcoeffs_y[i_glo_dofs] << "\n";
		}
                myfile_t.close();
	  }
	  else cout << "Unable to open file"; 
*/




//	const std::vector< std::string > filnam = m_session->GetFilenames();
//	cout << "loaded filename: " << filnam[0] << endl;  

    }

}

/**
 * $Log: CoupledLinearNS.cpp,v $
 **/
