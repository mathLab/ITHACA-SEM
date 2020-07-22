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
     * @param
     *
     * outarray[nvars] flux
     */
    template <class T, typename = typename std::enable_if
        <
            std::is_floating_point<T>::value ||
            tinysimd::is_vector_floating_point<T>::value
        >::type
    >
    inline void GetViscousFluxBilinearFormKernel(
        const int               nDim,
        const int               FluxDirection,
        const int               DerivDirection,
        const std::vector<T>    &inaverg,
        const std::vector<T>    &injumpp,
        const T                 &mu,
        std::vector<T>          &outarray)
    {
        // size_t nConvectiveFields   = inaverg.size();
        // size_t nPts = inaverg[nConvectiveFields - 1].size();
        // size_t nDim = nSpaceDim;

        // Constants
        size_t nDim_plus_one = nDim + 1;
        size_t FluxDirection_plus_one = FluxDirection + 1;
        size_t DerivDirection_plus_one = DerivDirection + 1;

        NekDouble gammaoPr = m_gamma / m_Prandtl;
        NekDouble one_minus_gammaoPr = 1.0 - gammaoPr;

        constexpr NekDouble OneThird = 1. / 3.;
        constexpr NekDouble TwoThird = 2. * OneThird;
        constexpr NekDouble FourThird = 4. * OneThird;


        if (DerivDirection == FluxDirection)
        {
            // necessary???
            outarray[0] = 0.0; // store 1x

            // load 1/rho
            T oneOrho = 1.0 / inaverg[0]; // load 1x
            // get vel, vel^2, sum of vel^2
            std::array<T, 3> u{}, u2{};
            T u2sum{};
            for (size_t d = 0; d < nDim; ++d)
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
            for (size_t d = 0; d < nDim; ++d)
            {
                size_t d_plus_one = d + 1;
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
            // Always zero
            outarray[0] = 0.0; // store 1x

            // load 1/rho
            T oneOrho = 1.0 / inaverg[0]; // load 1x
            // get vel, vel^2, sum of vel^2
            std::array<T, 3> u{}, u2{};
            T u2sum{};
            for (size_t d = 0; d < nDim; ++d)
            {
                size_t d_plus_one = d + 1;
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

            tmp1 = - u[FluxDirection] * injumpp[0] +
                injumpp[FluxDirection_plus_one];
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



  };
}
#endif
