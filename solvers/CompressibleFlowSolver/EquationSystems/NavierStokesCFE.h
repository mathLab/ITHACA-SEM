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
#include <CompressibleFlowSolver/Misc/EquationOfState.h>
#include <LibUtilities/BasicUtils/Smath.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField1D.h>

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

    typedef std::function<void (
            const Array<OneD, NekDouble>    &,
            const NekDouble                 &,
            const Array<OneD, NekDouble>    &,
                DNekMatSharedPtr            &)> GetdFlux_dDeriv;

  protected:
    std::string                         m_ViscosityType;
    NekDouble                           m_mu;
    NekDouble                           m_thermalConductivity;
    NekDouble                           m_Cp;
    NekDouble                           m_Cv;
    NekDouble                           m_Prandtl;
    NekDouble                           m_Twall;
    NekDouble                           m_mu0;
    std::string                         m_physicalSensorType;
    std::string                         m_smoothing;
    MultiRegions::ContField2DSharedPtr  m_C0Project2DExp;
    MultiRegions::ContField3DSharedPtr  m_C0Project3DExp;

    /// Equation of system for computing temperature
    EquationOfStateSharedPtr            m_eos;

    Array<OneD, GetdFlux_dDeriv>        m_GetdFlux_dDeriv_Array;

    NavierStokesCFE(const LibUtilities::SessionReaderSharedPtr& pSession,
                    const SpatialDomains::MeshGraphSharedPtr& pGraph);

    void GetViscousFluxVectorConservVar(
        const int                                                       nConvectiveFields,
        const int                                                       nDim,
        const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &qfields,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &outarray,
              Array< OneD, int >                                        &nonZeroIndex       =   NullInt1DArray,
        const Array<OneD, Array<OneD, NekDouble> >                      &normal             =   NullNekDoubleArrayofArray,
        const Array<OneD, NekDouble>                                    &ArtifDiffFactor    =   NullNekDouble1DArray);

    void GetPrimDerivFromConsDeriv(
        const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &qfields,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &outarray);

    void SpecialBndTreat(
        const int                                           nConvectiveFields,
              Array<OneD,       Array<OneD, NekDouble> >    &consvar);
    void ApplyFluxBndConds(
        const int                                           nConvectiveFields,
              Array<OneD,       Array<OneD, NekDouble> >    &flux);

    void GetArtificialViscosity(
        const Array<OneD, Array<OneD, NekDouble> >  &inarray,
              Array<OneD,             NekDouble  >  &muav);

    void GetViscousFluxBilinearForm(
        const int                                                       nConvectiveFields,
        const int                                                       FluxDirection,
        const int                                                       DerivDirection,
        const Array<OneD, const Array<OneD, NekDouble> >                &inaverg,
        const Array<OneD, const Array<OneD, NekDouble> >                &injumpp,
        const Array<OneD, NekDouble>                                    &mu,
        const Array<OneD, const Array<OneD, NekDouble> >                &auxVars,
              Array<OneD, Array<OneD, NekDouble> >                      &outarray);

    void CalcAuxiVarForBilinearFom(
        const int                                                       nConvectiveFields,
        const Array<OneD, const Array<OneD, NekDouble> >                &inaverg,
        Array<OneD, NekDouble>                                          &mu,
        Array<OneD, Array<OneD, NekDouble> >                            &auxVars);

    void CalcViscosity(
        const Array<OneD, const Array<OneD, NekDouble> >                &inaverg,
        Array<OneD, NekDouble>                                          &mu);
    
    virtual void v_InitObject();

    virtual void v_ExtraFldOutput(
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables);

    virtual void v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd);
    virtual void v_DoDiffusion_coeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
               Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd);

    virtual void v_DoDiffusionFlux(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble>>              &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

    virtual void v_GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >         &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor);
    virtual void v_GetViscousFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >         &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor);


    virtual void v_GetViscousSymmtrFluxConservVar(
            const int                                                       nConvectiveFields,
            const int                                                       nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >                      &inaverg,
            const Array<OneD, Array<OneD, NekDouble > >                     &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &outarray,
            Array< OneD, int >                                              &nonZeroIndex,
            const Array<OneD, Array<OneD, NekDouble> >                      &normals);

    void GetPhysicalAV(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield);

    void GetTracePhysicalAV();
    void Ducros( Array<OneD, NekDouble> &field );
    void C0Smooth(Array<OneD, NekDouble> &field);

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
  virtual void v_MinusDiffusionFluxJacDirctn(
        const int                                                       nDirctn,
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray);
    virtual void v_MinusDiffusionFluxJacDirctnElmt(
            const int                                                       nConvectiveFields,
            const int                                                       nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
            const Array<OneD, Array<OneD,  Array<OneD, NekDouble> > >       &locDerv,
            const Array<OneD, NekDouble>                                    &locmu,
            const Array<OneD, NekDouble>                                    &locDmuDT,
            const Array<OneD, NekDouble>                                    &normals,
            DNekMatSharedPtr                                                &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray);

    virtual void v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr                            &explist,
        const Array<OneD, const Array<OneD, NekDouble> >                &normals,
        const int                                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,

        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray,
        const int                                                       nfluxDir);
    virtual void v_GetFluxDerivJacDirctnElmt(
        const int                                                       nConvectiveFields,
        const int                                                       nElmtPnt,
        const int                                                       nDervDir,
        const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
        const Array<OneD, NekDouble>                                    &locmu,
        const Array<OneD, Array<OneD, NekDouble> >                      &locnormal,
        DNekMatSharedPtr                                                &wspMat,
        Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray);
    
    virtual void v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr                            &explist,
        const Array<OneD, const Array<OneD, NekDouble> >                &normals,
        const int                                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
              Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac);


    virtual void v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                  DNekMatSharedPtr                              &fluxJac);

    virtual void v_CalphysDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >  &qfield);

    virtual void v_CalcMuDmuDT(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            Array<OneD, NekDouble>                                          &mu,
            Array<OneD, NekDouble>                                          &DmuDT);
      
    /**
     * @brief return part of viscous Jacobian:
     * \todo flux derived with Qx=[drho_dx,drhou_dx,drhov_dx,drhoE_dx]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void GetdFlux_dQx_2D(
        const Array<OneD, NekDouble> &normals,
        const NekDouble &mu,
        const Array<OneD, NekDouble> &U,
        DNekMatSharedPtr &OutputMatrix );

    /**
     * @brief return part of viscous Jacobian:
     * \todo flux derived with Qx=[drho_dy,drhou_dy,drhov_dy,drhoE_dy]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void GetdFlux_dQy_2D(
        const Array<OneD, NekDouble> &normals,
        const NekDouble &mu,
        const Array<OneD, NekDouble> &U,
        DNekMatSharedPtr &OutputMatrix );

    /**
     * @brief return part of viscous Jacobian derived with Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=0)= dF_dQx;
     */
    void GetdFlux_dQx_3D(
        const Array<OneD, NekDouble> &normals,
        const NekDouble &mu,
        const Array<OneD, NekDouble> &U,
        DNekMatSharedPtr &OutputMatrix );

    /**
     * @brief return part of viscous Jacobian derived with Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=1)= dF_dQy;
     */
    void GetdFlux_dQy_3D(
        const Array<OneD, NekDouble> &normals,
        const NekDouble &mu,
        const Array<OneD, NekDouble> &U,
        DNekMatSharedPtr &OutputMatrix );


    /**
     * @brief return part of viscous Jacobian derived with Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Output: 3D 4*5 Matrix (flux about rho is zero)
     * OutputMatrix(dir=2)= dF_dQz;
     */
    void GetdFlux_dQz_3D(
        const Array<OneD, NekDouble> &normals,
        const NekDouble &mu,
        const Array<OneD, NekDouble> &U,
        DNekMatSharedPtr &OutputMatrix );


    /**
     * @brief return part of viscous Jacobian
     * Input:
     * normals:Point normals
     * mu: dynamicviscosity
     * dmu_dT: mu's derivative with T using Sutherland's law
     * U=[rho,rhou,rhov,rhoE]
     * Output: 3*4 Matrix (the flux about rho is zero)
     * OutputMatrix dFLux_dU,  the matrix sign is consistent with SIPG
    */
    void GetdFlux_dU_2D(
        const Array<OneD, NekDouble>                        &normals,
        const NekDouble                                     mu,
        const NekDouble                                     dmu_dT,
        const Array<OneD, NekDouble>                        &U,
        const Array<OneD, const Array<OneD, NekDouble> >    &qfield,
              DNekMatSharedPtr                              &OutputMatrix);

    /**
     * @brief return part of viscous Jacobian
     * Input:
     * normals:Point normals
     * mu: dynamicviscosity
     * dmu_dT: mu's derivative with T using Sutherland's law
     * U=[rho,rhou,rhov,rhow,rhoE]
     * Output: 4*5 Matrix (the flux about rho is zero)
     * OutputMatrix dFLux_dU,  the matrix sign is consistent with SIPG
    */
    void GetdFlux_dU_3D(
        const Array<OneD, NekDouble>                        &normals,
        const NekDouble                                     mu,
        const NekDouble                                     dmu_dT,
        const Array<OneD, NekDouble>                        &U,
        const Array<OneD, const Array<OneD, NekDouble> >    &qfield,
              DNekMatSharedPtr                              &OutputMatrix);

#endif

  };
}
#endif
