///////////////////////////////////////////////////////////////////////////////
//
// File EulerCFE.h
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
// Description: Euler equations in conservative variables without artificial diffusion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_EULERCFE_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_EULERCFE_H

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>


namespace Nektar
{  

  enum ProblemType
  {           
    eGeneral,          ///< No problem defined - Default Inital data
    eIsentropicVortex, ///< Isentropic Vortex
    eRinglebFlow,      ///< Ringleb Flow
    SIZE_ProblemType   ///< Length of enum list
  };
  
  const char* const ProblemTypeMap[] =
    {
      "General",
      "IsentropicVortex",
      "RinglebFlow"
    };
  
  class EulerCFE : public CompressibleFlowSystem
  {
  public:
      friend class MemoryManager<EulerCFE>;

    /// Creates an instance of this class.
    static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
    {
      SolverUtils::EquationSystemSharedPtr p = MemoryManager<EulerCFE>::AllocateSharedPtr(pSession);
      p->InitObject();
      return p;
    }
    /// Name of class.
    static std::string className;
    
    virtual ~EulerCFE();

    ///< problem type selector
    ProblemType     m_problemType;   
    
  protected:

    EulerCFE(const LibUtilities::SessionReaderSharedPtr& pSession);

    virtual void v_InitObject();

    /// Print a summary of time stepping parameters.
    virtual void v_PrintSummary(std::ostream &out);

    void DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time);
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
			  Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time);
    virtual void v_SetInitialConditions(
        NekDouble               initialtime = 0.0,
        bool                    dumpInitialConditions = true);
    virtual void v_EvaluateExactSolution(
        unsigned int            field,
        Array<OneD, NekDouble> &outfield,
        const NekDouble         time = 0.0);
      
  private:

    void SetBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble> >            &physarray, 
        NekDouble                                        time);

    /// Isentropic Vortex Test Case.
    void GetExactIsentropicVortex(
        int                                              field, 
        Array<OneD, NekDouble>                          &outarray, 
        NekDouble                                        time);
    void SetInitialIsentropicVortex(
        NekDouble                                        initialtime);
    void SetBoundaryIsentropicVortex(
        int                                              bcRegion, 
        NekDouble                                        time, 
        int cnt, Array<OneD, Array<OneD, NekDouble> >   &physarray);

    /// Ringleb Flow Test Case.
    void GetExactRinglebFlow(
        int                                             field, 
        Array<OneD, NekDouble>                         &outarray);
    void SetInitialRinglebFlow(
        void);
    void SetBoundaryRinglebFlow(
        int                                              bcRegion, 
        NekDouble                                        time, 
        int                                              cnt, 
        Array<OneD, Array<OneD, NekDouble> >            &physarray);
  };
}
#endif
