//////////////////////////////////////////////////////////////////////////////
//
// File ContField3DHomogeneous1D.cpp
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
// Description: Field definition for 3D domain with boundary
// conditions and a 1D homogeneous direction
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField2D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        ContField3DHomogeneous1D::ContField3DHomogeneous1D(void):
            DisContField3DHomogeneous1D()
        {
        }

        ContField3DHomogeneous1D::ContField3DHomogeneous1D(
                                const ContField3DHomogeneous1D &In):
                                DisContField3DHomogeneous1D (In,false)
        {

            bool False = false;
            ContField2DSharedPtr zero_plane =
                    std::dynamic_pointer_cast<ContField2D> (In.m_planes[0]);

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n] =   MemoryManager<ContField2D>::
                                        AllocateSharedPtr(*zero_plane,False);
            }

            SetCoeffPhys();
        }

        ContField3DHomogeneous1D::ContField3DHomogeneous1D(
                            const ContField3DHomogeneous1D &In,
                            const SpatialDomains::MeshGraphSharedPtr &graph2D,
                            const std::string                        &variable):
            DisContField3DHomogeneous1D (In, false)
        {
            ContField2DSharedPtr zero_plane_old =
                    std::dynamic_pointer_cast<ContField2D> (In.m_planes[0]);

            ContField2DSharedPtr zero_plane =
                        MemoryManager<ContField2D>::
                                    AllocateSharedPtr(*zero_plane_old,graph2D,
                                                            variable);

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n] =   MemoryManager<ContField2D>::
                                        AllocateSharedPtr(*zero_plane,graph2D,
                                                            variable);
            }

            SetCoeffPhys();

            if(variable.compare("DefaultVar") != 0)
            {
                SpatialDomains::BoundaryConditions bcs(m_session, graph2D);
                SetupBoundaryConditions(m_homogeneousBasis->GetBasisKey(),
                                        m_lhom, bcs,variable);
            }
        }

        ContField3DHomogeneous1D::~ContField3DHomogeneous1D()
        {
        }

        ContField3DHomogeneous1D::ContField3DHomogeneous1D(
                          const LibUtilities::SessionReaderSharedPtr &pSession,
                          const LibUtilities::BasisKey &HomoBasis,
                          const NekDouble lhom,
                          const bool useFFT,
                          const bool dealiasing,
                          const SpatialDomains::MeshGraphSharedPtr &graph2D,
                          const std::string &variable,
                          const bool CheckIfSingularSystem,
                          const Collections::ImplementationType ImpType):
            DisContField3DHomogeneous1D(pSession,HomoBasis,lhom,
                                        useFFT,dealiasing)
        {
            int i,n,nel;
            ContField2DSharedPtr plane_zero;
            ContField2DSharedPtr plane_two;

            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);
            m_graph = graph2D;

            // Plane zero (k=0 - cos) - singularaty check required for Poisson
            // problems
            plane_zero = MemoryManager<ContField2D>::AllocateSharedPtr(
                                        pSession, graph2D, variable, false,
                                        CheckIfSingularSystem, ImpType);

            plane_two  = MemoryManager<ContField2D>::AllocateSharedPtr(
                                        pSession, graph2D, variable, false,
                                        false, ImpType);

            m_exp = MemoryManager<LocalRegions::ExpansionVector>
                                        ::AllocateSharedPtr();

            for(n = 0; n < m_planes.size(); ++n)
            {
                // Plane zero and one (k=0 - cos and sin) - singularaty check
                // required for Poisson problems
                if(m_transposition->GetK(n) == 0)
                {
                    m_planes[n] = MemoryManager<ContField2D>
                        ::AllocateSharedPtr(*plane_zero, graph2D, variable,
                                                false, CheckIfSingularSystem);
                }
                else
                {
                    // For k > 0 singularty check not required anymore -
                    // creating another ContField2D to avoid Assembly Map copy
                    // TODO: We may want to deal with it in a more efficient
                    // way in the future.
                    m_planes[n] = MemoryManager<ContField2D>
                            ::AllocateSharedPtr(*plane_two, graph2D, variable,
                                                false, false);
                }

                nel = m_planes[n]->GetExpSize();

                for(i = 0; i < nel; ++i)
                {
                    (*m_exp).push_back(m_planes[n]->GetExp(i));
                }
            }

            nel = GetExpSize();

            SetCoeffPhys();

            // Do not set up BCs if default variable
            if(variable.compare("DefaultVar") != 0)
            {
                SetupBoundaryConditions(HomoBasis,lhom,bcs,variable);
            }
        }


        void ContField3DHomogeneous1D::v_ImposeDirichletConditions(
                                                Array<OneD,NekDouble>& outarray)
        {
            Array<OneD, NekDouble> tmp;
            int ncoeffs = m_planes[0]->GetNcoeffs();

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->ImposeDirichletConditions(tmp = outarray +
                                                       n*ncoeffs);
            }
        }

        void ContField3DHomogeneous1D::v_FillBndCondFromField(void)
        {
            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->FillBndCondFromField();
            }
        }

        void ContField3DHomogeneous1D::v_FillBndCondFromField(const int nreg)
        {
            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->FillBndCondFromField(nreg);
            }
        }

        /**
         *
         */
        void  ContField3DHomogeneous1D::v_LocalToGlobal(bool useComm)
        {
            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->LocalToGlobal(useComm);
            }
        }


        /**
         *
         */
        void  ContField3DHomogeneous1D::v_GlobalToLocal(void)
        {
            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->GlobalToLocal();
            }
        }


        /**
         *
         */
        void ContField3DHomogeneous1D::v_SmoothField(
                Array<OneD,NekDouble> &field)
        {
            int cnt = 0;
            Array<OneD, NekDouble> tmp;

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->SmoothField(tmp = field + cnt);

                cnt  += m_planes[n]->GetTotPoints();
            }
        }


        void ContField3DHomogeneous1D::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const MultiRegions::VarFactorsMap &varfactors,
                const Array<OneD, const NekDouble> &dirForcing,
                const bool PhysSpaceForcing)
        {

            int n;
            int cnt = 0;
            int cnt1 = 0;
            NekDouble beta;
            StdRegions::ConstFactorMap new_factors;

            Array<OneD, NekDouble> e_out;
            Array<OneD, NekDouble> fce(inarray.size());
            Array<OneD, const NekDouble> wfce;

            // Fourier transform forcing function
            if(m_WaveSpace)
            {
                fce = inarray;
            }
            else
            {
                HomogeneousFwdTrans(inarray, fce);
            }

            bool smode = false;

            if (m_homogeneousBasis->GetBasisType() ==
                LibUtilities::eFourierHalfModeRe ||
                m_homogeneousBasis->GetBasisType() ==
                LibUtilities::eFourierHalfModeIm )
            {
                smode = true;
            }

            for(n = 0; n < m_planes.size(); ++n)
            {
                if(n != 1 || m_transposition->GetK(n) != 0 || smode)
                {

                    beta = 2*M_PI*(m_transposition->GetK(n))/m_lhom;
                    new_factors = factors;
                    // add in Homogeneous Fourier direction and SVV if turned on
                    new_factors[StdRegions::eFactorLambda] +=
                                                beta*beta*(1+GetSpecVanVisc(n));

                    wfce = (PhysSpaceForcing)? fce+cnt:fce+cnt1;
                    m_planes[n]->HelmSolve(wfce,
                                           e_out = outarray + cnt1,
                                           new_factors, varcoeff,
                                           varfactors, dirForcing,
                                           PhysSpaceForcing);
                }

                cnt  += m_planes[n]->GetTotPoints();
                cnt1 += m_planes[n]->GetNcoeffs();
            }
        }

        /**
         * Reset the GlobalLinSys Manager
         */
        void ContField3DHomogeneous1D::v_ClearGlobalLinSysManager(void)
        {
            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->ClearGlobalLinSysManager();
            }
        }

    } // end of namespace
} //end of namespace
