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
  /**
   *
   *
   **/
  class NavierStokesCFE : public CompressibleFlowSystem
  {
  public:
      friend class MemoryManager<NavierStokesCFE>;

    // Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
      SolverUtils::EquationSystemSharedPtr p =
          MemoryManager<NavierStokesCFE>::AllocateSharedPtr(pSession, pGraph);
      p->InitObject();
      return p;
    }
    // Name of class
    static std::string className;

    virtual ~NavierStokesCFE();

  protected:
    std::string                         m_ViscosityType;
    NekDouble                           m_Cp;
    NekDouble                           m_Cv;
    NekDouble                           m_Prandtl;

    NekDouble                            m_Twall;
    NekDouble                           m_muRef;
    NekDouble                           m_thermalConductivityRef;
    Array<OneD, NekDouble>              m_mu;
    Array<OneD, NekDouble>              m_thermalConductivity;


    NavierStokesCFE(const LibUtilities::SessionReaderSharedPtr& pSession,
                    const SpatialDomains::MeshGraphSharedPtr& pGraph);
    
    void GetViscousFluxVectorConservVar(
        const int                                                nDim,
        const Array<OneD, Array<OneD, NekDouble> >               &inarray,
        const TensorOfArray3D<NekDouble>                         &qfields,
        TensorOfArray3D<NekDouble>                               &outarray,
        Array< OneD, int >                                       
            &nonZeroIndex       =   NullInt1DArray,    
        const Array<OneD, Array<OneD, NekDouble> >               
            &normal             =   NullNekDoubleArrayofArray,           
        const Array<OneD, NekDouble>                             
            &ArtifDiffFactor    =   NullNekDouble1DArray);
    void GetViscousSymmtrFluxConservVar(
            const int                                           nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >          &inaverg,
            const Array<OneD, Array<OneD, NekDouble > >         &inarray,
            TensorOfArray3D<NekDouble>                          &outarray,
            Array< OneD, int >                                  &nonZeroIndex,    
            const Array<OneD, Array<OneD, NekDouble> >          &normals);
    
    void SpecialBndTreat(
              Array<OneD,       Array<OneD, NekDouble> >    &consvar);

    void GetArtificialViscosity(
        const Array<OneD, Array<OneD, NekDouble> >  &inarray,
              Array<OneD,             NekDouble  >  &muav);

    void GetViscousFluxBilinearForm(
        const int                                           nSpaceDim,
        const int                                           FluxDirection,
        const int                                           DerivDirection,
        const Array<OneD, const Array<OneD, NekDouble> >    &inaverg,
        const Array<OneD, const Array<OneD, NekDouble> >    &injumpp,
        Array<OneD, Array<OneD, NekDouble> >                &outarray);
    
    virtual void v_InitObject();

    virtual void v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd);
    virtual void v_DoDiffusion_coeff(
        const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
        Array<OneD, Array<OneD, NekDouble> >                &outarray,
        const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >          &pBwd);

    virtual void v_GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >         &physfield,
        TensorOfArray3D<NekDouble>                         &derivatives,
        TensorOfArray3D<NekDouble>                         &viscousTensor);
    virtual void v_GetViscousFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >         &physfield,
        TensorOfArray3D<NekDouble>                         &derivatives,
        TensorOfArray3D<NekDouble>                         &viscousTensor);

    virtual void v_GetFluxPenalty(
        const Array<OneD, Array<OneD, NekDouble> > &uFwd,
        const Array<OneD, Array<OneD, NekDouble> > &uBwd,
              Array<OneD, Array<OneD, NekDouble> > &penaltyCoeff);

    void GetViscosityAndThermalCondFromTemp(
        const Array<OneD, NekDouble> &temperature,
              Array<OneD, NekDouble> &mu,
              Array<OneD, NekDouble> &thermalCond);

  };
}
#endif
