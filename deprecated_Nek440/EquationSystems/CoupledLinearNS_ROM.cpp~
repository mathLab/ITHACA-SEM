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
#include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS_ROM.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

using namespace std;

namespace Nektar
{

    string CoupledLinearNS_ROM::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("CoupledLinearisedNS_ROM", CoupledLinearNS_ROM::create);


    CoupledLinearNS_ROM::CoupledLinearNS_ROM(const LibUtilities::SessionReaderSharedPtr &pSession):
        UnsteadySystem(pSession),
        CoupledLinearNS(pSession),
        m_zeroMode(false)
    {



    }

    void CoupledLinearNS_ROM::v_InitObject()
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


	std::stringstream sstm_bwdtrans;
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

	std::stringstream sstm_cartmap0;
	sstm_cartmap0 << "cartmap0.txt";
//	std::string result_cartmap0 = sstm_cartmap0.str();
//	std::basic_string<char, std::char_traits<char>, std::allocator<char> > result_cartmap0 = sstm_cartmap0.str();
//        std::ofstream myfile_cartmap0 (result_cartmap0);
        std::ofstream myfile_cartmap0 ("cartmap0.txt");
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

        // locExp->IProductWRTBase(tmpphys,coeffs);                                for all tmpphys

	std::stringstream sstm_IP;
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

	std::stringstream sstm_IP_d0;
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


	//  m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);

	std::stringstream sstm_IPp;
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
        


		std::stringstream sstmb;
		sstmb << "bmap_" << n << ".txt";
		std::string resultb = sstmb.str();


//          std::ofstream myfile_bmap (resultb);
          std::ofstream myfile_bmap ("bmap.txt");
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

//	    std::cout << "saving current imap " << resulti << std::endl;

//          std::ofstream myfile_imap (resulti);
          std::ofstream myfile_imap ("imap.txt");
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

      //      StdRegions::ConstFactorMap factors;
      //      factors[StdRegions::eFactorLambda] = lambda/m_kinvis;
      //      LocalRegions::MatrixKey helmkey(StdRegions::eHelmholtz,
      //                                      locExp->DetShapeType(),
      //                                      *locExp,
      //                                      factors);

      //          DNekScalMat &HelmMat = *(locExp->as<LocalRegions::Expansion>()
      //                                         ->GetLocMatrix(helmkey));

		}
	        myfile_IPp.close();
	}
	else std::cout << "Unable to open file";

    }

    void CoupledLinearNS_ROM::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation)
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

    void CoupledLinearNS_ROM::SetUpCoupledMatrix(const NekDouble lambda,  const Array< OneD, Array< OneD, NekDouble > > &Advfield, bool IsLinearNSEquation,const int HomogeneousMode, CoupledSolverMatrices &mat, CoupledLocalToGlobalC0ContMapSharedPtr &locToGloMap, const NekDouble lambda_imag)
    {
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
            nsize_int  [n] = (nvel*m_fields[m_velocity[0]]->GetExp(eid)->GetNcoeffs()*nz_loc - nsize_bndry[n]);
            nsize_p[n] = m_pressure->GetExp(eid)->GetNcoeffs()*nz_loc;
            nsize_p_m1[n] = nsize_p[n]-nz_loc;
      /*      std::cout << "nsize_bndry[n] " << nsize_bndry[n] << std::endl;
            std::cout << "nsize_bndry_p1[n] " << nsize_bndry_p1[n] << std::endl;
            std::cout << "nsize_int[n] " << nsize_int[n] << std::endl;
            std::cout << "nsize_p[n] " << nsize_p[n] << std::endl;
            std::cout << "nsize_p_m1[n] " << nsize_p_m1[n] << std::endl; */
        }
        
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
            Array<OneD, NekDouble> Ah_data = Ah->GetPtr();
            int AhRows = k;
            DNekMatSharedPtr B  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            Array<OneD, NekDouble> B_data = B->GetPtr();
            DNekMatSharedPtr C  = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry,nint,zero);
            Array<OneD, NekDouble> C_data = C->GetPtr();
            DNekMatSharedPtr D  = MemoryManager<DNekMat>::AllocateSharedPtr(nint, nint,zero);
            Array<OneD, NekDouble> D_data = D->GetPtr();
            
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
              

		std::stringstream sstm_l00;
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





////////////////////////////////////////////////////////////////////////////////////////////////


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

////////////////////////////////////////////////////////////////////////////////////////////////


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
                

		std::stringstream sstm;
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
                        }
                        
                        for(j = 0; j < nimap; ++j)
                        {
                            B_data[i+k*nbmap + (j+k*nimap)*nbndry] += m_kinvis*HelmMatScale*HelmMat_data[bmap[i]+HelmMatRows*imap[j]];
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
                                        Ah_data[j+nv*nbmap + (i+nv*nbmap)*AhRows] +=
                                        coeffs[bmap[j]];
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        C_data[i+nv*nbmap + (j+nv*nimap)*nbndry] += 
                                        coeffs[imap[j]];
                                    }
                                }
                                
                                // copy into column major storage. 
                                m_pressure->GetExp(eid)->IProductWRTBase(deriv,pcoeffs);
                                for(j = 0; j < nz_loc; ++j)
                                {
                                    Blas::Dcopy(psize,&(pcoeffs)[0],1,
                                                Dbnd->GetRawPtr() + 
                                                ((nz_loc*k+j)*bmap.num_elements() + i)*nsize_p[n]+ j*psize,1);
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
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        C_data[i+(nv*nz_loc+n1)*nbmap + 
                                        (j+(k*nz_loc+n1)*nimap)*nbndry] += 
                                        coeffs[imap[j]];
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
                        }
                        
                        for(j = 0; j < nimap; ++j)
                        {
                            D_data[i+k*nimap + (j+k*nimap)*nint] += m_kinvis*HelmMatScale*HelmMat_data[imap[i]+HelmMatRows*imap[j]];
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
                                        B_data[j+nv*nbmap + (i+nv*nimap)*nbndry] +=
                                        coeffs[bmap[j]];
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        D_data[j+nv*nimap + (i+nv*nimap)*nint] += 
                                        coeffs[imap[j]];
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
                                    }
                                    
                                    for(j = 0; j < nimap; ++j)
                                    {
                                        D_data[j+(k*nz_loc+n1)*nimap +  
                                        (i+(nv*nz_loc+n1)*nimap)*nint] += 
                                        coeffs[imap[j]];
                                    }
                                }
                            }
                        }
                    }
                }
                
                
            for(int n1 = 0; n1 < D->GetColumns(); ++n1)
            {
                for(int n2 = 0; n2 < D->GetRows(); ++n2)
                {
                    
  //                  cout << "D->GetRows() " << D->GetRows() << std::endl ;
  //                  cout << "D->GetColumns() " << D->GetColumns() << std::endl ;
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
        cout << "Matrix Setup Costs: " << timer.TimePerTest(1) << endl;
        
        
        timer.Start();
        // Set up global coupled boundary solver. 
        // This is a key to define the solution matrix type
        // currently we are giving it a argument of eLInearAdvectionReaction 
        // since this then makes the matrix storage of type eFull
        MultiRegions::GlobalLinSysKey key(StdRegions::eLinearAdvectionReaction,locToGloMap);
        mat.m_CoupledBndSys = MemoryManager<MultiRegions::GlobalLinSysDirectStaticCond>::AllocateSharedPtr(key,m_fields[0],pAh,pBh,pCh,pDh,locToGloMap);
        mat.m_CoupledBndSys->Initialise(locToGloMap);
        timer.Stop();
        cout << "Multilevel condensation: " << timer.TimePerTest(1) << endl;
    }


    void CoupledLinearNS_ROM::Solve(void)
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

    void CoupledLinearNS_ROM::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing)
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
    
    void CoupledLinearNS_ROM::SolveLinearNS(const Array<OneD, Array<OneD, NekDouble> > &forcing,  Array<OneD, MultiRegions::ExpListSharedPtr> &fields, MultiRegions::ExpListSharedPtr &pressure,  const int mode)
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
          ofstream myfilef0 ("forcing0.txt");
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
        
          ofstream myfile1LocGloMap ("LocGloMap.txt");
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
        
          ofstream myfilebndmap ("BndCondCoeffsToGlobalCoeffsMap.txt");
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
            
     		ofstream myfile_bnd_cond_elem ("bnd_cond_elem.txt");
		if (myfile_bnd_cond_elem.is_open())
	  	{
			myfile_bnd_cond_elem << std::setprecision(17) << bndCondExp.num_elements() << "\n";
                	myfile_bnd_cond_elem.close();
	  	}
	  	else cout << "Unable to open file"; 

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

       
     		ofstream myfile_t (outname_t);
		if (myfile_t.is_open())
	  	{
			for(int count_bnd = 0; count_bnd < bndCondCoeffs.num_elements(); count_bnd++)
			{
				myfile_t << std::setprecision(17) << bndCondCoeffs[count_bnd] << "\n";
			}
                	myfile_t.close();
	  	}
	  	else cout << "Unable to open file"; 


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

          ofstream myfiledirbnd ("NumGlobalDirBndCoeffs.txt");
	  if (myfiledirbnd.is_open())
	  {
		myfiledirbnd << std::setprecision(17) << m_locToGloMap[mode]->GetNumGlobalDirBndCoeffs();
                myfiledirbnd.close();
	  }
	  else cout << "Unable to open file";   

	cout << "fh_bnd.num_elements() " << fh_bnd.num_elements() << std::endl;   
	cout << "bnd.num_elements() " << bnd.num_elements() << std::endl;  
	
	for(i = 0; i <  fh_bnd.num_elements(); ++i)
	{
//		cout << "bnd[i] " << i << " " << bnd[i] << std::endl;   
	}
        m_mat[mode].m_CoupledBndSys->Solve(fh_bnd,bnd,m_locToGloMap[mode]);

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

	// since this function is virtual it actually can be instantiated here
    void CoupledLinearNS_ROM::v_DoSolve(void)
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

    void CoupledLinearNS_ROM::v_DoInitialise(void)
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
                EvaluateFunction(fieldStr,AdvField,"AdvectionVelocity");
                
                SetUpCoupledMatrix(0.0,AdvField,false);

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

    void CoupledLinearNS_ROM::v_Output(void)
    {    
        std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_fields.num_elements()+1);
        std::vector<std::string> variables(m_fields.num_elements()+1);
        int i;
        
        for(i = 0; i < m_fields.num_elements(); ++i)
        {
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i]   = m_boundaryConditions->GetVariable(i);
        }



//	cout << "m_pressure->GetPhys() " << m_pressure->GetPhys() << endl;
//	cout << "m_pressure->GetCoeffs() " << m_pressure->GetCoeffs() << endl;

	// use LocalToGlobal (const Array< OneD, const NekDouble > &inarray, Array< OneD, NekDouble > &outarray, bool useComm=true)
//        Array<OneD, NekDouble> glo_fieldcoeffs(m_fields[0]->GetNcoeffs(), 0.0);


/*	for(int counter_global_num = 6000; counter_global_num < fieldcoeffs[0].num_elements(); counter_global_num++)
	{
	        Array<OneD, NekDouble> glo_fieldcoeffs(counter_global_num, 0.0);
		cout << "cgi: " << counter_global_num << endl;
		m_fields[0]->LocalToGlobal(fieldcoeffs[0], glo_fieldcoeffs);
		cout << "cgi: " << counter_global_num << endl;

		try
		{
			cout << "cgi: " << counter_global_num << " value " << glo_fieldcoeffs[counter_global_num-2] << " ";
		}
	        catch (const std::runtime_error&)
	        {
        	//	return 1;
	        }
	} */

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
		glo_coll_fieldcoeffs[counter_phys_index] = glo_fieldcoeffs;
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





//	const std::vector< std::string > filnam = m_session->GetFilenames();
//	cout << "loaded filename: " << filnam[0] << endl;  

    }

}

/**
 * $Log: CoupledLinearNS.cpp,v $
 **/
