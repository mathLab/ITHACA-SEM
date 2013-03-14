///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFE.h
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
// Description: NavierStokes equations in conservative variable
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_NAVIERSTOKESCFE_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_NAVIERSTOKESCFE_H

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

namespace Nektar
{  

  enum ProblemType
  {           
    eGeneral,          ///< No problem defined - Default Inital data
    SIZE_ProblemType   ///< Length of enum list
  };
  
  const char* const ProblemTypeMap[] =
    {
      "General"
    };
  
  /**
   * 
   * 
   **/
  class NavierStokesCFE : public CompressibleFlowSystem
  {
  public:
      friend class MemoryManager<NavierStokesCFE>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
    {
      SolverUtils::EquationSystemSharedPtr p = MemoryManager<NavierStokesCFE>::AllocateSharedPtr(pSession);
      p->InitObject();
      return p;
    }
    /// Name of class
    static std::string className;
    
    virtual ~NavierStokesCFE();

    ///< problem type selector
    ProblemType                                     m_problemType;   
    
  protected:
    NavierStokesCFE(const LibUtilities::SessionReaderSharedPtr& pSession);

    virtual void v_InitObject();
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
        NekDouble                               initialtime = 0.0,
        bool                                    dumpInitialConditions = true);
    
    private:
      void SetBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble> >             &physarray, 
        NekDouble                                         time);
  };
}
#endif
