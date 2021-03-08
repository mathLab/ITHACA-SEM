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
#include <CompressibleFlowSolver/Misc/EquationOfState.h>
#include <LibUtilities/BasicUtils/Smath.hpp>
#include <MultiRegions/ContField.h>

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
    /// flag to switch between constant viscosity and Sutherland
    /// an enum could be added for more options
    bool                                m_is_mu_variable{false};
    /// flag to switch between IP and LDG
    /// an enum could be added for more options
    bool                                m_is_diffIP{false};

    NekDouble                           m_Cp;
    NekDouble                           m_Cv;
    NekDouble                           m_Prandtl;
    NekDouble                           m_Twall;
    NekDouble                           m_mu0;
    std::string                         m_physicalSensorType;
    std::string                         m_smoothing;
    MultiRegions::ContFieldSharedPtr    m_C0ProjectExp;

    /// Equation of system for computing temperature
    EquationOfStateSharedPtr            m_eos;

    Array<OneD, GetdFlux_dDeriv>        m_GetdFlux_dDeriv_Array;
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

    virtual void v_InitObject();

    virtual void v_ExtraFldOutput(
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables);

    virtual void v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd);
    virtual void v_DoDiffusionCoeff(
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

      void GetPhysicalAV(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield);
      
      void GetTracePhysicalAV();
      void Ducros( Array<OneD, NekDouble> &field );
      void C0Smooth(Array<OneD, NekDouble> &field);
  
      virtual void v_MinusDiffusionFluxJacDirctnElmt(
            const int                                  nConvectiveFields,
            const int                                  nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >  &locVars,
            const TensorOfArray3D<NekDouble>            &locDerv,
            const Array<OneD, NekDouble>                &locmu,
            const Array<OneD, NekDouble>                &locDmuDT,
            const Array<OneD, NekDouble>                &normals,
            DNekMatSharedPtr                            &wspMat,
            Array<OneD, Array<OneD, NekDouble> >        &PntJacArray);

    virtual void v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr              &explist,
        const Array<OneD, const Array<OneD, NekDouble> >  &normals,
        const int                                         nDervDir,
        const Array<OneD, const Array<OneD, NekDouble> >  &inarray,
        TensorOfArray5D<NekDouble>                        &ElmtJacArray,
        const int                                         nfluxDir);
    virtual void v_GetFluxDerivJacDirctnElmt(
        const int                                         nConvectiveFields,
        const int                                         nElmtPnt,
        const int                                         nDervDir,
        const Array<OneD, Array<OneD, NekDouble> >        &locVars,
        const Array<OneD, NekDouble>                      &locmu,
        const Array<OneD, Array<OneD, NekDouble> >        &locnormal,
        DNekMatSharedPtr                                  &wspMat,
        Array<OneD, Array<OneD, NekDouble> >              &PntJacArray);
    
    virtual void v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr              &explist,
        const Array<OneD, const Array<OneD, NekDouble> >  &normals,
        const int                                         nDervDir,
        const Array<OneD, const Array<OneD, NekDouble> >  &inarray,
              Array<OneD, Array<OneD, DNekMatSharedPtr>>  &ElmtJac);

    virtual void v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                  DNekMatSharedPtr                              &fluxJac);

    virtual void v_CalcPhysDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield);

    virtual void v_CalcMuDmuDT(
            const Array<OneD, const Array<OneD, NekDouble> >     &inarray,
            Array<OneD, NekDouble>                               &mu,
            Array<OneD, NekDouble>                               &DmuDT);
      
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

    virtual void v_GetFluxPenalty(
        const Array<OneD, Array<OneD, NekDouble> > &uFwd,
        const Array<OneD, Array<OneD, NekDouble> > &uBwd,
              Array<OneD, Array<OneD, NekDouble> > &penaltyCoeff);

    void GetViscosityAndThermalCondFromTemp(
        const Array<OneD, NekDouble> &temperature,
              Array<OneD, NekDouble> &mu,
              Array<OneD, NekDouble> &thermalCond);

    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline void GetViscosityAndThermalCondFromTempKernel(
        const T& temperature, T& mu, T& thermalCond)
    {
        GetViscosityFromTempKernel(temperature, mu);
        NekDouble tRa = m_Cp / m_Prandtl;
        thermalCond = tRa * mu;
    }

    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline void GetViscosityFromTempKernel(
        const T& temperature, T& mu)
    {
        // Variable viscosity through the Sutherland's law
        if (m_is_mu_variable)
        {
            mu = m_varConv->GetDynamicViscosity(temperature);
        }
        else
        {
            mu = m_muRef;
        }
    }

    /**
     * @brief Calculate diffusion flux using the Jacobian form.
     *
     * @param in
     *
     * @param out
     * outarray[nvars] flux
     */
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline void GetViscousFluxBilinearFormKernel(
        const unsigned short nDim,
        const unsigned short FluxDirection,
        const unsigned short DerivDirection,
        // these need to be a pointers because of the custom allocator of vec_t
        const T*             inaverg,
        const T*             injumpp,
        const T&             mu,
        T*                   outarray)
    {
        // Constants
        unsigned short nDim_plus_one = nDim + 1;
        unsigned short FluxDirection_plus_one = FluxDirection + 1;
        unsigned short DerivDirection_plus_one = DerivDirection + 1;

        NekDouble gammaoPr = m_gamma / m_Prandtl;
        NekDouble one_minus_gammaoPr = 1.0 - gammaoPr;

        constexpr NekDouble OneThird = 1. / 3.;
        constexpr NekDouble TwoThird = 2. * OneThird;
        constexpr NekDouble FourThird = 4. * OneThird;


        if (DerivDirection == FluxDirection)
        {
            // rho flux always zero
            outarray[0] = 0.0; // store 1x

            // load 1/rho
            T oneOrho = 1.0 / inaverg[0]; // load 1x
            // get vel, vel^2, sum of vel^2
            std::array<T, 3> u, u2;
            T u2sum{};
            for (unsigned short d = 0; d < nDim; ++d)
            {
                u[d] = inaverg[d+1] * oneOrho; // load 1x
                u2[d] = u[d]*u[d];
                u2sum += u2[d];
            }

            // get E - sum v^2
            T E_minus_u2sum = inaverg[nDim_plus_one]; // load 1x
            E_minus_u2sum *= oneOrho;
            E_minus_u2sum -= u2sum;

            // get nu = mu/rho
            T nu = mu * oneOrho; // load 1x

            // ^^^^ above is almost the same for both loops

            T tmp1 = OneThird * u2[FluxDirection] + u2sum;
            tmp1 += gammaoPr * E_minus_u2sum;
            tmp1 *= injumpp[0]; //load 1x

            T tmp2 = gammaoPr * injumpp[nDim_plus_one] - tmp1; //load 1x

            // local var for energy output
            T outTmpE = 0.0;
            for (unsigned short d = 0; d < nDim; ++d)
            {
                unsigned short d_plus_one = d + 1;
                //flux[rhou, rhov, rhow]
                T outTmpD = injumpp[d_plus_one] - u[d] * injumpp[0];
                outTmpD *= nu;
                //flux rhoE
                T tmp3 = one_minus_gammaoPr * u[d];
                outTmpE +=  tmp3 * injumpp[d_plus_one];

                if (d == FluxDirection)
                {
                    outTmpD *= FourThird;
                    T tmp4 = OneThird * u[FluxDirection];
                    outTmpE += tmp4 * injumpp[FluxDirection_plus_one];
                }

                outarray[d_plus_one] = outTmpD; //store 1x
            }

            outTmpE += tmp2;
            outTmpE *= nu;
            outarray[nDim_plus_one] = outTmpE; //store 1x

        }
        else
        {
            // rho flux always zero
            outarray[0] = 0.0; // store 1x

            // load 1/rho
            T oneOrho = 1.0 / inaverg[0]; // load 1x
            // get vel, vel^2, sum of vel^2
            std::array<T, 3> u, u2;
            T u2sum{};
            for (unsigned short d = 0; d < nDim; ++d)
            {
                unsigned short d_plus_one = d + 1;
                u[d] = inaverg[d_plus_one] * oneOrho; // load 1x
                u2[d] = u[d]*u[d];
                u2sum += u2[d];
                // Not all directions are set
                // one could work out the one that is not set
                outarray[d_plus_one] = 0.0; // store 1x
            }

            // get E - sum v^2
            T E_minus_u2sum = inaverg[nDim_plus_one]; // load 1x
            E_minus_u2sum *= oneOrho;
            E_minus_u2sum -= u2sum;

            // get nu = mu/rho
            T nu = mu * oneOrho; // load 1x

            // ^^^^ above is almost the same for both loops

            T tmp1 = u[DerivDirection] * injumpp[0] -
                injumpp[DerivDirection_plus_one]; // load 2x
            tmp1 *= TwoThird;
            outarray[FluxDirection_plus_one] = nu * tmp1; // store 1x

            tmp1 = injumpp[FluxDirection_plus_one] - u[FluxDirection] * injumpp[0];
            outarray[DerivDirection_plus_one] = nu * tmp1; // store 1x

            tmp1 = OneThird * u[FluxDirection] * u[DerivDirection];
            tmp1 *= injumpp[0];

            T tmp2 = TwoThird * u[FluxDirection] *
                injumpp[DerivDirection_plus_one];

            tmp1 += tmp2;

            tmp1 = u[DerivDirection] * injumpp[FluxDirection_plus_one] -
                tmp1;
            outarray[nDim_plus_one] = nu * tmp1; // store 1x

        }
    }



    /**
     * @brief Return the flux vector for the IP diffusion problem, based on
     * conservative variables
     */
    template<bool IS_TRACE>
    void GetViscousFluxVectorConservVar(
        const int                                              nDim,
        const Array<OneD, Array<OneD, NekDouble> >             &inarray,
        const TensorOfArray3D<NekDouble>                       &qfields,
        TensorOfArray3D<NekDouble>                             &outarray,
        Array< OneD, int >                                     &nonZeroIndex,
        const Array<OneD, Array<OneD, NekDouble> >             &normal,
        const Array<OneD, NekDouble>                           &ArtifDiffFactor)
    {
        size_t nConvectiveFields = inarray.size();
        size_t nPts = inarray[0].size();
        size_t n_nonZero = nConvectiveFields - 1;

        // max outfield dimensions
        constexpr unsigned short nOutMax = 3 - 2 * IS_TRACE;
        constexpr unsigned short nVarMax = 5;
        constexpr unsigned short nDimMax = 3;

        // vector loop
        using namespace tinysimd;
        using vec_t = simd<NekDouble>;
        size_t sizeVec = (nPts / vec_t::width) * vec_t::width;
        size_t p = 0;

        for (; p < sizeVec; p += vec_t::width)
        {
            // there is a significant penalty to use std::vector
            alignas(vec_t::alignment) std::array<vec_t, nVarMax> inTmp,
                qfieldsTmp, outTmp;
            alignas(vec_t::alignment) std::array<vec_t, nDimMax> normalTmp;
            alignas(vec_t::alignment) std::array<vec_t, nVarMax * nOutMax>
                outArrTmp{{}};

            // rearrenge and load data
            for (size_t f = 0; f < nConvectiveFields; ++f)
            {
                inTmp[f].load(&(inarray[f][p]), is_not_aligned);
                // zero output vector
                if (IS_TRACE)
                {
                    outArrTmp[f] = 0.0;
                }
                else
                {
                    for (int d = 0; d < nDim; ++d)
                    {
                        outArrTmp[f + nConvectiveFields * d] = 0.0;
                    }
                }
            }
            if (IS_TRACE)
            {
                for (size_t d = 0; d < nDim; ++d)
                {
                    normalTmp[d].load(&(normal[d][p]), is_not_aligned);
                }
            }

            // get temp
            vec_t temperature = m_varConv->GetTemperature(inTmp.data());
            // get viscosity
            vec_t mu;
            GetViscosityFromTempKernel(temperature, mu);

            for (size_t nderiv = 0; nderiv < nDim; ++nderiv)
            {
                // rearrenge and load data
                for (size_t f = 0; f < nConvectiveFields; ++f)
                {
                    qfieldsTmp[f].load(&(qfields[nderiv][f][p]), is_not_aligned);
                }

                for (size_t d = 0; d < nDim; ++d)
                {
                    GetViscousFluxBilinearFormKernel(nDim, d, nderiv,
                        inTmp.data(), qfieldsTmp.data(), mu, outTmp.data());

                    if (IS_TRACE)
                    {
                        for (size_t f = 0; f < nConvectiveFields; ++f)
                        {
                            outArrTmp[f] += normalTmp[d] * outTmp[f];
                        }
                    }
                    else
                    {
                        for (size_t f = 0; f < nConvectiveFields; ++f)
                        {
                            outArrTmp[f + nConvectiveFields * d] += outTmp[f];
                        }
                    }
                }
            }

            // store data
            if (IS_TRACE)
            {
                for (int f = 0; f < nConvectiveFields; ++f)
                {
                    outArrTmp[f].store(&(outarray[0][f][p]), is_not_aligned);
                }
            }
            else
            {
                for (int d = 0; d < nDim; ++d)
                {
                    for (int f = 0; f < nConvectiveFields; ++f)
                    {
                        outArrTmp[f + nConvectiveFields * d].store(
                            &(outarray[d][f][p]), is_not_aligned);
                    }
                }
            }
        }

        // scalar loop
        for (; p < nPts; ++p)
        {
            std::array<NekDouble, nVarMax> inTmp, qfieldsTmp, outTmp;
            std::array<NekDouble, nDimMax> normalTmp;
            std::array<NekDouble, nVarMax * nOutMax> outArrTmp{{}};
            // rearrenge and load data
            for (int f = 0; f < nConvectiveFields; ++f)
            {
                inTmp[f] = inarray[f][p];
                // zero output vector
                if (IS_TRACE)
                {
                    outArrTmp[f] = 0.0;
                }
                else
                {
                    for (int d = 0; d < nDim; ++d)
                    {
                        outArrTmp[f + nConvectiveFields * d] = 0.0;
                    }
                }
            }
            

            if (IS_TRACE)
            {
                for (int d = 0; d < nDim; ++d)
                {
                    normalTmp[d] = normal[d][p];
                }
            }

            // get temp
            NekDouble temperature = m_varConv->GetTemperature(inTmp.data());
            // get viscosity
            NekDouble mu;
            GetViscosityFromTempKernel(temperature, mu);

            for (int nderiv = 0; nderiv < nDim; ++nderiv)
            {
                // rearrenge and load data
                for (int f = 0; f < nConvectiveFields; ++f)
                {
                    qfieldsTmp[f] = qfields[nderiv][f][p];
                }

                for (int d = 0; d < nDim; ++d)
                {
                    GetViscousFluxBilinearFormKernel(nDim, d, nderiv,
                        inTmp.data(), qfieldsTmp.data(), mu, outTmp.data());

                    if (IS_TRACE)
                    {
                        for (size_t f = 0; f < nConvectiveFields; ++f)
                        {
                            outArrTmp[f] += normalTmp[d] * outTmp[f];
                        }
                    }
                    else
                    {
                        for (size_t f = 0; f < nConvectiveFields; ++f)
                        {
                            outArrTmp[f + nConvectiveFields * d] += outTmp[f];
                        }
                    }

                }
            }

            // store data
            if (IS_TRACE)
            {
                for (int f = 0; f < nConvectiveFields; ++f)
                {
                    outarray[0][f][p] = outArrTmp[f];
                }
            }
            else
            {
                for (int d = 0; d < nDim; ++d)
                {
                    for (int f = 0; f < nConvectiveFields; ++f)
                    {
                        outarray[d][f][p] = outArrTmp[f + nConvectiveFields * d];
                    }
                }
            }
        }


        // this loop would need to be brought up into the main loop so that it
        // can be vectorized as well
        if (ArtifDiffFactor.size())
        {
            n_nonZero = nConvectiveFields;

            for (size_t p = 0; p < nPts; ++p)
            {
                for (int d = 0; d < nDim; ++d)
                {
                    if (IS_TRACE)
                    {
                        NekDouble tmp = ArtifDiffFactor[p] * normal[d][p];

                        for (int j = 0; j < nConvectiveFields; ++j)
                        {
                            outarray[0][j][p] += tmp * qfields[d][j][p];
                        }
                    }
                    else
                    {
                        for (int j = 0; j < nConvectiveFields; ++j)
                        {
                            outarray[d][j][p] += ArtifDiffFactor[p] * qfields[d][j][p];
                        }
                    }
                }
            }
        }

        // this is always the same, it should be initialized with the IP class
        nonZeroIndex = Array< OneD, int > {n_nonZero, 0};
        for (int i = 1; i < n_nonZero + 1; ++i)
        {
            nonZeroIndex[n_nonZero - i] =   nConvectiveFields - i;
        }
    }




  };
}
#endif
