///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFEAxisym.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_NSCFEAXISYM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_NSCFEAXISYM_H

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFE.h>

namespace Nektar
{
  /**
   *
   *
   **/
  class NavierStokesCFEAxisym : public NavierStokesCFE
  {
  public:
      friend class MemoryManager<NavierStokesCFEAxisym>;

    // Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
      SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<NavierStokesCFEAxisym>
                ::AllocateSharedPtr(pSession, pGraph);
      p->InitObject();
      return p;
    }
    // Name of class
    static std::string className;

    virtual ~NavierStokesCFEAxisym();

  protected:
    Array<OneD, Array<OneD, NekDouble> > m_viscousForcing;

    NavierStokesCFEAxisym(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph);

    virtual void v_InitObject();

    virtual void v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >   &pBwd);

    virtual void v_GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >         &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor);
    virtual void v_GetViscousFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >         &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        boost::ignore_unused(physfield, derivatives, viscousTensor);
        NEKERROR(ErrorUtil::efatal,
                 "Dealiased flux not implemented for axisymmetric case");
    }
  };
}
#endif
