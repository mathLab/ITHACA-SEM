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

#include <IncNavierStokesSolver/EquationSystems/AdvectionTerm.h>
#include <cstdio>
#include <cstdlib>

#include <cmath>

#include <string>
namespace Nektar
{
    /**
     * Basic construnctor
     */
    AdvectionTerm::AdvectionTerm(void)
    {     
		
    }
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    AdvectionTerm::AdvectionTerm(
            LibUtilities::CommSharedPtr                 pComm,
            LibUtilities::SessionReaderSharedPtr        pSession,
            SpatialDomains::MeshGraphSharedPtr          pGraph,
            SpatialDomains::BoundaryConditionsSharedPtr pBoundaryConditions)
        : m_comm(pComm),
          m_session(pSession),
          m_graph(pGraph),
          m_boundaryConditions(pBoundaryConditions)
	{
        // Set space dimension for use in class
        m_spacedim = pGraph->GetSpaceDimension();
        m_expdim   = pGraph->GetMeshDimension();

        // Save the basename of input file name for output details.
		m_sessionName = pSession->GetFilename();
		m_sessionName = m_sessionName.substr(0,
											 m_sessionName.find_last_of("."));

        if(m_boundaryConditions->SolverInfoExists("PROJECTION"))
        {
            std::string ProjectStr
            = m_boundaryConditions->GetSolverInfo("PROJECTION");

            if((ProjectStr == "Continuous")||(ProjectStr == "Galerkin")||
               (ProjectStr == "CONTINUOUS")||(ProjectStr == "GALERKIN"))
            {
                m_projectionType = EquationSystem::eGalerkin;
            }
            else if(ProjectStr == "DisContinuous")
            {
                m_projectionType = EquationSystem::eDiscontinuousGalerkin;
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
            m_projectionType = EquationSystem::eGalerkin;
        }
	}

	AdvectionTerm::~AdvectionTerm()
	{
	}
	
    void AdvectionTerm::v_DoAdvection(
                               Array<OneD, MultiRegions::ExpListSharedPtr > &m_fields,
                               const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                               Array<OneD, Array<OneD, NekDouble> > &outarray,
                               Array<OneD, NekDouble> &wk)
    {
        ASSERTL0(false,"This advection form is not defined in this class");
    }

	int AdvectionTerm::NoCaseStringCompare(const string & s1, const string& s2)
    {
        //if (s1.size() < s2.size()) return -1;
        //if (s1.size() > s2.size()) return 1;
		
        string::const_iterator it1=s1.begin();
        string::const_iterator it2=s2.begin();
		
        //stop when either string's end has been reached
        while ( (it1!=s1.end()) && (it2!=s2.end()) )
        {
            if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
            {
                // return -1 to indicate smaller than, 1 otherwise
                return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1;
            }
			
            //proceed to the next character in each string
            ++it1;
            ++it2;
        }
		
        size_t size1=s1.size();
        size_t size2=s2.size();// cache lengths
		
        //return -1,0 or 1 according to strings' lengths
        if (size1==size2)
        {
            return 0;
        }
		
        return (size1 < size2) ? -1 : 1;
    }
	
	
} //end of namespace
