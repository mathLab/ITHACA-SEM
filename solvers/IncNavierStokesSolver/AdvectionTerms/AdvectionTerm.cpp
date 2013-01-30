///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionTerm.cpp
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
// Description: Base class for Navier-Stokes advection term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.h>
#include <string>


namespace Nektar
{

    AdvectionTermFactory& GetAdvectionTermFactory()
    {
        typedef Loki::SingletonHolder<AdvectionTermFactory,
            Loki::CreateUsingNew,
            Loki::NoDestroy > Type;
        return Type::Instance();
    }

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    AdvectionTerm::AdvectionTerm(
            const LibUtilities::SessionReaderSharedPtr&        pSession,
            const SpatialDomains::MeshGraphSharedPtr&          pGraph)
        : m_session(pSession),
          m_graph(pGraph)
	{
	}

    void AdvectionTerm::v_InitObject()
    {
        // Set space dimension for use in class
        m_spacedim = m_graph->GetSpaceDimension();
        m_expdim   = m_graph->GetMeshDimension();

        // Save the basename of input file name for output details.
        m_sessionName = m_session->GetFilename();
        m_sessionName = m_sessionName.substr(0, m_sessionName.find_last_of("."));
        
        if(m_session->DefinesSolverInfo("PROJECTION"))
        {
            std::string ProjectStr
                = m_session->GetSolverInfo("PROJECTION");
            
            if((ProjectStr == "Continuous")||(ProjectStr == "Galerkin")||
               (ProjectStr == "CONTINUOUS")||(ProjectStr == "GALERKIN"))
            {
                m_projectionType = MultiRegions::eGalerkin;
            }
            else if((ProjectStr == "MixedCGDG")||(ProjectStr == "Mixed_CG_Discontinuous"))
            {
                m_projectionType = MultiRegions::eMixed_CG_Discontinuous;
            }                        
            else if(ProjectStr == "DisContinuous")
            {
                m_projectionType = MultiRegions::eDiscontinuous;
            }
            else
            {
                ASSERTL0(false,"PROJECTION value not recognised");
            }
        }
        else
        {
            cerr << "Projection type not specified in SOLVERINFO,"
                "defaulting to continuous Galerkin" << endl;
            m_projectionType = MultiRegions::eGalerkin;
        }
		
        m_CoeffState = MultiRegions::eLocal;
        m_dealiasing = false;
        
        if(m_session->DefinesSolverInfo("DEALIASING"))
        {
            m_dealiasing = true;
        }

        m_session->MatchSolverInfo("SPECTRALHPDEALIASING","True",m_specHP_dealiasing,false);
        if(m_specHP_dealiasing == false)
        {
            m_session->MatchSolverInfo("SPECTRALHPDEALIASING","On",m_specHP_dealiasing,false);
        }

        m_session->MatchSolverInfo("ModeType","SingleMode",m_SingleMode,false);
        m_session->MatchSolverInfo("ModeType","HalfMode",m_HalfMode,false);
    }
    
    AdvectionTerm::~AdvectionTerm()
    {
    }
    
    void AdvectionTerm::DoAdvection(Array<OneD, MultiRegions::ExpListSharedPtr> &pFields, 
                                    const int nConvectiveFields, 
                                    const Array<OneD, int> &vel_loc, 
                                    const Array<OneD, const Array<OneD, NekDouble> > &pInarray, 
                                    Array<OneD, Array<OneD, NekDouble> > &pOutarray,
                                    NekDouble time,
                                    Array<OneD, NekDouble> &pWk)
    {
		int i,j;
        int VelDim           = vel_loc.num_elements();        
        int nqtot            = pFields[0]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
	
        ASSERTL1(nConvectiveFields == pInarray.num_elements(),"Number of convective fields and Inarray are not compatible");
        
        for(i = 0; i < VelDim; ++i)
        {
            if(pFields[i]->GetWaveSpace() && !m_SingleMode && !m_HalfMode)
            {
                j = vel_loc[i];
                velocity[i] = Array<OneD, NekDouble>(nqtot,0.0);
                pFields[i]->HomogeneousBwdTrans(pInarray[j],velocity[i]);
            }
            else 
            {
                velocity[i] = pInarray[vel_loc[i]];
            }
        }
        

        DoAdvection(pFields,velocity,pInarray,pOutarray,time,pWk);
    }

    void AdvectionTerm::DoAdvection(Array<OneD, MultiRegions::ExpListSharedPtr> &pFields, 
                                    const Array<OneD, const Array<OneD, NekDouble> > &velocity, 
                                    const Array<OneD, const Array<OneD, NekDouble> > &pInarray, 
                                    Array<OneD, Array<OneD, NekDouble> > &pOutarray,
                                    NekDouble time,
                                    Array<OneD, NekDouble> &pWk)
    {
        int i;
        int nqtot  = pFields[0]->GetTotPoints();
        int VelDim = velocity.num_elements();
        Array<OneD, NekDouble > Deriv;
	

        m_nConvectiveFields = pInarray.num_elements();

        // Set up Derivative work space;
        if(pWk.num_elements())
        {
            ASSERTL0(pWk.num_elements() >= nqtot*VelDim,"Workspace is not sufficient");
            Deriv = pWk;
        }
        else
        {
            Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
        }
	
	
        for(i=0; i< m_nConvectiveFields; ++i)
        {
            v_ComputeAdvectionTerm(pFields,velocity,pInarray[i],pOutarray[i],i,time,Deriv);
            Vmath::Neg(nqtot,pOutarray[i],1);
        }
    }
    
} //end of namespace
