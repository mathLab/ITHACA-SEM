//////////////////////////////////////////////////////////////////////////////
//
// File DisContField3DHomogeneous1D.cpp
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
// conditions using LDG flux and a 1D homogeneous direction
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField2D.h>


namespace Nektar
{
    namespace MultiRegions
    {

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(void)
        : ExpList3DHomogeneous1D(),
          m_bndCondExpansions(),
          m_bndCondBndWeight(),
          m_bndConditions()
        {
        }

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::BasisKey               &HomoBasis,
            const NekDouble                             lhom,
            const bool                                  useFFT,
            const bool                                  dealiasing):
            ExpList3DHomogeneous1D(pSession,HomoBasis,lhom,useFFT,dealiasing),
              m_bndCondExpansions(),
              m_bndCondBndWeight(),
              m_bndConditions()
        {
        }

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(
            const DisContField3DHomogeneous1D &In,
            const bool                         DeclarePlanesSetCoeffPhys)
            : ExpList3DHomogeneous1D (In,false),
              m_bndCondExpansions    (In.m_bndCondExpansions),
              m_bndCondBndWeight     (In.m_bndCondBndWeight),
              m_bndConditions        (In.m_bndConditions)
        {
            if (DeclarePlanesSetCoeffPhys)
            {
                DisContField2DSharedPtr zero_plane =
                    std::dynamic_pointer_cast<DisContField2D> (In.m_planes[0]);

                for(int n = 0; n < m_planes.size(); ++n)
                {
                    m_planes[n] =
                        MemoryManager<DisContField2D>::
                          AllocateSharedPtr(*zero_plane, false);
                }

                SetCoeffPhys();
            }
        }

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::BasisKey               &HomoBasis,
            const NekDouble                             lhom,
            const bool                                  useFFT,
            const bool                                  dealiasing,
            const SpatialDomains::MeshGraphSharedPtr   &graph2D,
            const std::string                          &variable,
            const Collections::ImplementationType       ImpType):
            ExpList3DHomogeneous1D(pSession, HomoBasis, lhom, useFFT,
                                   dealiasing),
              m_bndCondExpansions(),
              m_bndCondBndWeight(),
              m_bndConditions()
        {
            int i, n, nel;
            DisContField2DSharedPtr plane_zero;
            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);

            // note that nzplanes can be larger than nzmodes
            m_planes[0] = plane_zero = MemoryManager<DisContField2D>::
                AllocateSharedPtr(pSession, graph2D, variable, true, false,
                                  ImpType);

            m_exp = MemoryManager<LocalRegions::ExpansionVector>
                ::AllocateSharedPtr();

            nel = m_planes[0]->GetExpSize();

            for (i = 0; i < nel; ++i)
            {
                (*m_exp).push_back(m_planes[0]->GetExp(i));
            }

            for (n = 1; n < m_planes.size(); ++n)
            {
                m_planes[n] = MemoryManager<DisContField2D>::
                    AllocateSharedPtr(*plane_zero, graph2D,
                                      variable, true, false);
                for(i = 0; i < nel; ++i)
                {
                    (*m_exp).push_back((*m_exp)[i]);
                }
            }

            // Set up trace object.
            Array<OneD, ExpListSharedPtr> trace(m_planes.size());
            for (n = 0; n < m_planes.size(); ++n)
            {
                trace[n] = m_planes[n]->GetTrace();
            }

            m_trace = MemoryManager<ExpList2DHomogeneous1D>::AllocateSharedPtr(
                               pSession, HomoBasis, lhom, useFFT,
                               dealiasing, trace);

            // Setup default optimisation information
            nel = GetExpSize();

            SetCoeffPhys();

            // Do not set up BCs if default variable
            if(variable.compare("DefaultVar") != 0)
            {
                SetupBoundaryConditions(HomoBasis, lhom, bcs, variable);
            }

            SetUpDG();
        }

        /**
         * @brief Default destructor
         */
        DisContField3DHomogeneous1D::~DisContField3DHomogeneous1D()
        {
        }

        void DisContField3DHomogeneous1D::SetupBoundaryConditions(
            const LibUtilities::BasisKey       &HomoBasis,
            const NekDouble                     lhom,
            SpatialDomains::BoundaryConditions &bcs,
            const std::string                   variable)
        {
            int n, cnt = 0;

            // Setup an ExpList2DHomogeneous1D expansion for boundary
            // conditions and link to class declared in m_planes
            const SpatialDomains::BoundaryRegionCollection  &bregions =
                bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions =
                bcs.GetBoundaryConditions();

            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpListSharedPtr>(
                bregions.size());
            m_bndConditions = m_planes[0]->UpdateBndConditions();

            m_bndCondBndWeight = Array<OneD, NekDouble> {bregions.size(),0.0};

            int nplanes = m_planes.size();
            Array<OneD, MultiRegions::ExpListSharedPtr>
                PlanesBndCondExp(nplanes);

            for (auto &it : bregions)
            {
                SpatialDomains::BoundaryConditionShPtr boundaryCondition =
                    GetBoundaryCondition(bconditions, it.first, variable);
                for (n = 0; n < nplanes; ++n)
                {
                    PlanesBndCondExp[n] = m_planes[n]->
                        UpdateBndCondExpansion(cnt);
                }

                // Initialise comm for the boundary regions
                auto comm = boundaryCondition->GetComm();
                int size = boundaryCondition->GetComm()->GetSize();

		if(size > 1)
		{
                    // It seems to work either way
                    // comm->SplitComm(1,size);
                    comm->SplitComm(m_StripZcomm->GetSize(),
                       size/m_StripZcomm->GetSize());
		}

                m_bndCondExpansions[cnt++] =
                    MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::
                    AllocateSharedPtr(m_session, HomoBasis, lhom,
                                      m_useFFT, false,
                                      PlanesBndCondExp, comm);
            }
            v_EvaluateBoundaryConditions(0.0, variable);
        }

        void DisContField3DHomogeneous1D::v_HelmSolve(
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

            // Transform forcing function in half-physical space if required
            if (m_WaveSpace)
            {
                fce = inarray;
            }
            else
            {
                HomogeneousFwdTrans(inarray,fce);
            }

            for (n = 0; n < m_planes.size(); ++n)
            {
                if (n != 1 || m_transposition->GetK(n) != 0)
                {
                    beta = 2*M_PI*(m_transposition->GetK(n))/m_lhom;
                    new_factors = factors;
                    // add in Homogeneous Fourier direction and SVV if turned on
                    new_factors[StdRegions::eFactorLambda] +=
                        beta*beta*(1+GetSpecVanVisc(n));

                    wfce = (PhysSpaceForcing)? fce+cnt:fce+cnt1;
                    m_planes[n]->HelmSolve(
                        wfce,
                        e_out = outarray + cnt1,
                        new_factors, varcoeff, varfactors,
                        dirForcing,
                        PhysSpaceForcing);
                }

                cnt  += m_planes[n]->GetTotPoints();
                cnt1 += m_planes[n]->GetNcoeffs();
            }
        }

        void DisContField3DHomogeneous1D::v_EvaluateBoundaryConditions(
            const NekDouble   time,
            const std::string varName,
            const NekDouble   x2_in,
            const NekDouble   x3_in)
        {
	    boost::ignore_unused(x2_in,x3_in);
            int i;
            int npoints;
            int nbnd = m_bndCondExpansions.size();
            MultiRegions::ExpListSharedPtr locExpList;

            for (i = 0; i < nbnd; ++i)
            {
                if (time == 0.0 || m_bndConditions[i]->IsTimeDependent())
                {
                    locExpList = m_bndCondExpansions[i];
                    npoints    = locExpList->GetNpoints();

                    Array<OneD, NekDouble> x0(npoints, 0.0);
                    Array<OneD, NekDouble> x1(npoints, 0.0);
                    Array<OneD, NekDouble> x2(npoints, 0.0);
                    Array<OneD, NekDouble> valuesFile(npoints, 1.0), valuesExp(npoints, 1.0);

                    locExpList->GetCoords(x0, x1, x2);

                    if (m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eDirichlet)
                    {
                        SpatialDomains::DirichletBCShPtr bcPtr = std::static_pointer_cast<
                            SpatialDomains::DirichletBoundaryCondition>(
                                m_bndConditions[i]);
			std::string filebcs = bcPtr->m_filename;
			std::string exprbcs = bcPtr->m_expr;

                        if (filebcs != "")
                        {
                            ExtractFileBCs(filebcs, bcPtr->GetComm(), varName, locExpList);
                            valuesFile = locExpList->GetPhys();
                        }

                        if (exprbcs != "")
                        {
                            LibUtilities::Equation  condition =
                                std::static_pointer_cast<SpatialDomains::
                                    DirichletBoundaryCondition >(
                                    m_bndConditions[i])->m_dirichletCondition;

                            condition.Evaluate(x0, x1, x2, time, valuesExp);
                        }

                        Vmath::Vmul(npoints, valuesExp, 1, valuesFile, 1,
                                    locExpList->UpdatePhys(), 1);

                        // set wave space to false since have set up phys values
                        locExpList->SetWaveSpace(false);

                        locExpList->FwdTrans_BndConstrained(
                            locExpList->GetPhys(),
                            locExpList->UpdateCoeffs());
                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                             == SpatialDomains::eNeumann)
                    {
                        SpatialDomains::NeumannBCShPtr bcPtr = std::static_pointer_cast<
                            SpatialDomains::NeumannBoundaryCondition>(
                                m_bndConditions[i]);

			std::string filebcs = bcPtr->m_filename;

                        if (filebcs != "")
                        {
                            ExtractFileBCs(filebcs, bcPtr->GetComm(), varName, locExpList);
                        }
                        else
                        {
                            LibUtilities::Equation condition = std::
                                static_pointer_cast<SpatialDomains::
                                                    NeumannBoundaryCondition>(
                                    m_bndConditions[i])->m_neumannCondition;

                            condition.Evaluate(x0, x1, x2, time,
                                               locExpList->UpdatePhys());

                            locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                        locExpList->UpdateCoeffs());
                        }
                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                             == SpatialDomains::eRobin)
                    {
                        SpatialDomains::RobinBCShPtr bcPtr = std::static_pointer_cast<
                            SpatialDomains::RobinBoundaryCondition>(
                                m_bndConditions[i]);
			std::string filebcs = bcPtr->m_filename;

                        if (filebcs != "")
                        {
                            ExtractFileBCs(filebcs, bcPtr->GetComm(), varName, locExpList);
                        }
                        else
                        {
                            LibUtilities::Equation condition = std::
                                static_pointer_cast<SpatialDomains::
                                                    RobinBoundaryCondition>(
                                    m_bndConditions[i])->m_robinFunction;

                            condition.Evaluate(x0, x1, x2, time,
                                               locExpList->UpdatePhys());

                        }

                        locExpList->IProductWRTBase(locExpList->GetPhys(),
                                                    locExpList->UpdateCoeffs());

                    }
                    else if (m_bndConditions[i]->GetBoundaryConditionType()
                             == SpatialDomains::ePeriodic)
                    {
                        continue;
                    }
                    else
                    {
                        ASSERTL0(false, "This type of BC not implemented yet");
                    }
                }
            }
        }

        std::shared_ptr<ExpList> &DisContField3DHomogeneous1D::
            v_UpdateBndCondExpansion(int i)
        {
            return UpdateBndCondExpansion(i);
        }

        Array<OneD, SpatialDomains::BoundaryConditionShPtr>
            &DisContField3DHomogeneous1D::v_UpdateBndConditions()
        {
            return UpdateBndConditions();
        }

        void DisContField3DHomogeneous1D::GetBoundaryToElmtMap(
            Array<OneD, int> &ElmtID,
            Array<OneD,int> &EdgeID)
        {

            if(m_BCtoElmMap.size() == 0)
            {
                Array<OneD, int> ElmtID_tmp;
                Array<OneD, int> EdgeID_tmp;

                m_planes[0]->GetBoundaryToElmtMap(ElmtID_tmp, EdgeID_tmp);
                int nel_per_plane = m_planes[0]->GetExpSize();
                int nplanes = m_planes.size();

                int MapSize = ElmtID_tmp.size();

                m_BCtoElmMap = Array<OneD, int>(nplanes*MapSize);
                m_BCtoEdgMap = Array<OneD, int>(nplanes*MapSize);

                // If this mesh (or partition) has no BCs, skip this step
                if (MapSize > 0)
                {
                    int i ,j, n, cnt;
                    int cntPlane = 0;
                    for (cnt=n=0; n < m_bndCondExpansions.size(); ++n)
                    {
                        int planeExpSize = m_planes[0]
                                                ->GetBndCondExpansions()[n]
                                                ->GetExpSize();
                        for (i = 0; i < planeExpSize ; ++i, ++cntPlane)
                        {
                            for(j = 0; j < nplanes; j++)
                            {
                                m_BCtoElmMap[cnt+i+j*planeExpSize] =
                                        ElmtID_tmp[cntPlane]+j*nel_per_plane;
                                m_BCtoEdgMap[cnt+i+j*planeExpSize] =
                                        EdgeID_tmp[cntPlane];
                            }
                        }
                        cnt += m_bndCondExpansions[n]->GetExpSize();
                    }
                }
            }
            ElmtID = m_BCtoElmMap;
            EdgeID = m_BCtoEdgMap;
        }

        void DisContField3DHomogeneous1D::v_GetBndElmtExpansion(int i,
                            std::shared_ptr<ExpList> &result,
                            const bool DeclareCoeffPhysArrays)
        {
            int n, cnt, nq;
            int offsetOld, offsetNew;

            std::vector<unsigned int> eIDs;
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);

            // Skip other boundary regions
            for (cnt = n = 0; n < i; ++n)
            {
                cnt += m_bndCondExpansions[n]->GetExpSize();
            }

            // Populate eIDs with information from BoundaryToElmtMap
            for (n = 0; n < m_bndCondExpansions[i]->GetExpSize(); ++n)
            {
                eIDs.push_back(ElmtID[cnt+n]);
            }

            // Create expansion list
            result =
                MemoryManager<ExpList3DHomogeneous1D>::AllocateSharedPtr
                    (*this, eIDs);

            // Copy phys and coeffs to new explist
            if ( DeclareCoeffPhysArrays)
            {
                Array<OneD, NekDouble> tmp1, tmp2;
                for (n = 0; n < result->GetExpSize(); ++n)
                {
                    nq = GetExp(ElmtID[cnt+n])->GetTotPoints();
                    offsetOld = GetPhys_Offset(ElmtID[cnt+n]);
                    offsetNew = result->GetPhys_Offset(n);
                    Vmath::Vcopy(nq, tmp1 = GetPhys()+ offsetOld, 1,
                                tmp2 = result->UpdatePhys()+ offsetNew, 1);

                    nq = GetExp(ElmtID[cnt+n])->GetNcoeffs();
                    offsetOld = GetCoeff_Offset(ElmtID[cnt+n]);
                    offsetNew = result->GetCoeff_Offset(n);
                    Vmath::Vcopy(nq, tmp1 = GetCoeffs()+ offsetOld, 1,
                                tmp2 = result->UpdateCoeffs()+ offsetNew, 1);
                }
            }

            // Set wavespace value
            result->SetWaveSpace(GetWaveSpace());
        }

        void DisContField3DHomogeneous1D::GetBCValues(       Array<OneD, NekDouble> & BndVals,
                                                       const Array<OneD, NekDouble> & TotField,
                                                             int                      BndID )
        {
            LocalRegions::ExpansionSharedPtr   elmt;
            LocalRegions::Expansion1DSharedPtr temp_BC_exp;

            Array<OneD, const NekDouble> tmp_Tot;
            Array<OneD, NekDouble>       tmp_BC;

            int cnt = 0;
            int pos = 0;
            int exp_size, exp_size_per_plane, elmtID, boundaryID;
            int offset, exp_dim;

            for (int k = 0; k < m_planes.size(); k++)
            {
                for (int n = 0; n < m_bndConditions.size(); ++n)
                {
                    exp_size = m_bndCondExpansions[n]->GetExpSize();
                    exp_size_per_plane = exp_size/m_planes.size();

                    for (int i = 0; i < exp_size_per_plane; i++)
                    {
                        if(n == BndID)
                        {
                            elmtID      = m_BCtoElmMap[cnt];
                            boundaryID  = m_BCtoEdgMap[cnt];
                            exp_dim     = m_bndCondExpansions[n]->
                                GetExp(i+k*exp_size_per_plane)->GetTotPoints();
                            offset      = GetPhys_Offset(elmtID);
                            elmt        = GetExp(elmtID);
                            temp_BC_exp = std::dynamic_pointer_cast<
                                LocalRegions::Expansion1D>(
                                    m_bndCondExpansions[n]->GetExp(
                                        i + k * exp_size_per_plane )
                                );

                            elmt->GetEdgePhysVals(boundaryID, temp_BC_exp,
                                                  tmp_Tot = TotField + offset,
                                                  tmp_BC = BndVals + pos);
                            pos        += exp_dim;
                        }
                        cnt++;
                    }
                }
            }
        }

        void DisContField3DHomogeneous1D::NormVectorIProductWRTBase(
            Array<OneD, const NekDouble> &V1,
            Array<OneD, const NekDouble> &V2,
            Array<OneD, NekDouble>       &outarray,
            int                           BndID)
        {
            LocalRegions::ExpansionSharedPtr   elmt;
            LocalRegions::Expansion1DSharedPtr temp_BC_exp;

            Array<OneD, NekDouble> tmp_V1;
            Array<OneD, NekDouble> tmp_V2;
            Array<OneD, NekDouble> tmp_outarray;

            int cnt = 0;
            int exp_size, exp_size_per_plane, elmtID, Phys_offset, Coef_offset;

            for(int k = 0; k < m_planes.size(); k++)
            {
                for(int n = 0; n < m_bndConditions.size(); ++n)
                {
                    exp_size = m_bndCondExpansions[n]->GetExpSize();
                    exp_size_per_plane = exp_size/m_planes.size();

                    for(int i = 0; i < exp_size_per_plane; i++)
                    {
                        if(n == BndID)
                        {
                            elmtID = m_BCtoElmMap[cnt];

                            Phys_offset = m_bndCondExpansions[n]->
                                GetPhys_Offset(i+k*exp_size_per_plane);
                            Coef_offset = m_bndCondExpansions[n]->
                                GetCoeff_Offset(i+k*exp_size_per_plane);

                            elmt = GetExp(elmtID);
                            temp_BC_exp = std::dynamic_pointer_cast<
                                LocalRegions::Expansion1D>(
                                    m_bndCondExpansions[n]->GetExp(
                                        i + k * exp_size_per_plane )
                                );

                            temp_BC_exp->NormVectorIProductWRTBase(
                                tmp_V1 = V1 + Phys_offset,
                                tmp_V2 = V2 + Phys_offset,
                                tmp_outarray = outarray + Coef_offset);
                        }
                        cnt++;
                    }
                }
            }
        }

        void DisContField3DHomogeneous1D::v_ExtractTracePhys(
                        Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(m_physState == true,
                     "Field must be in physical state to extract trace space.");

            v_ExtractTracePhys(m_phys, outarray);
        }

        /**
         * @brief This method extracts the trace (edges in 2D) for each plane
         * from the field @a inarray and puts the values in @a outarray.
         *
         * It assumes the field is C0 continuous so that it can overwrite the
         * edge data when visited by the two adjacent elements.
         *
         * @param inarray   An array containing the 2D data from which we wish
         *                  to extract the edge data.
         * @param outarray  The resulting edge information.
         */
        void DisContField3DHomogeneous1D::v_ExtractTracePhys(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray)
        {
            int nPoints_plane = m_planes[0]->GetTotPoints();
            int nTracePts = m_planes[0]->GetTrace()->GetTotPoints();

            for (int i = 0; i < m_planes.size(); ++i)
            {
                Array<OneD, NekDouble> inarray_plane(nPoints_plane, 0.0);
                Array<OneD, NekDouble> outarray_plane(nPoints_plane, 0.0);

                Vmath::Vcopy(nPoints_plane,
                             &inarray[i*nPoints_plane], 1,
                             &inarray_plane[0], 1);

                m_planes[i]->ExtractTracePhys(inarray_plane, outarray_plane);

                Vmath::Vcopy(nTracePts,
                             &outarray_plane[0], 1,
                             &outarray[i*nTracePts], 1);
            }
        }

        /**
         */
        void DisContField3DHomogeneous1D::v_GetBoundaryNormals( int                                    i,
                                                                Array<OneD, Array<OneD, NekDouble> > & normals )
        {
            int                              expdim = GetCoordim(0);
            int                              coordim = 3;
            Array<OneD, NekDouble>           tmp;
            LocalRegions::ExpansionSharedPtr elmt;

            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);

            // Initialise result
            normals = Array<OneD, Array<OneD, NekDouble> > (coordim);
            for (int j = 0; j < coordim; ++j)
            {
                normals[j] = Array<OneD, NekDouble> ( GetBndCondExpansions()[i]->GetTotPoints(), 0.0 );
            }

            // Skip other boundary regions
            int cnt = 0;
            for( int n = 0; n < i; ++n )
            {
                cnt += GetBndCondExpansions()[n]->GetExpSize();
            }

            int offset;
            for( int n = 0; n < GetBndCondExpansions()[i]->GetExpSize(); ++n )
            {
                offset = GetBndCondExpansions()[i]->GetPhys_Offset(n);
                int nq = GetBndCondExpansions()[i]->GetExp(n)->GetTotPoints();

                elmt   = GetExp(ElmtID[cnt+n]);
                const Array<OneD, const Array<OneD, NekDouble> > normalsElmt
                            = elmt->GetSurfaceNormal(EdgeID[cnt+n]);
                // Copy to result
                for (int j = 0; j < expdim; ++j)
                {
                    Vmath::Vcopy(nq, normalsElmt[j], 1,
                                     tmp = normals[j] + offset, 1);
                }
            }
        }

        /**
        * @brief Set up all DG member variables and maps.
        */
        void DisContField3DHomogeneous1D::SetUpDG()
        {
            const int nPlanes     = m_planes.size();
            const int nTracePlane = m_planes[0]->GetTrace()->GetExpSize();

            // Get trace map from first plane.
            AssemblyMapDGSharedPtr traceMap = m_planes[0]->GetTraceMap();
            const Array<OneD, const int> &traceBndMap
                = traceMap->GetBndCondIDToGlobalTraceID();
            int mapSize = traceBndMap.size();

            // Set up trace boundary map
            m_traceBndMap = Array<OneD, int>(nPlanes * mapSize);

            int i, n, e, cnt = 0, cnt1 = 0;

            for (i = 0; i < m_bndCondExpansions.size(); ++i)
            {
                int nExp      = m_bndCondExpansions[i]->GetExpSize();
                int nPlaneExp = nExp / nPlanes;

                if (m_bndConditions[i]->GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }

                for (n = 0; n < nPlanes; ++n)
                {
                    const int offset = n * nTracePlane;
                    for (e = 0; e < nPlaneExp; ++e)
                    {
                        m_traceBndMap[cnt++] = offset + traceBndMap[cnt1+e];
                    }
                }

                cnt1 += nPlaneExp;
            }
        }
    } // end of namespace
} //end of namespace
