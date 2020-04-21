///////////////////////////////////////////////////////////////////////////////
//
// File Expansion2D.cpp
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
// Description: File for Expansion2D routines
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion3D.h>
#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/Geometry2D.h>
#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
    namespace LocalRegions
    {
        Expansion2D::Expansion2D(SpatialDomains::Geometry2DSharedPtr pGeom) :
            StdExpansion(), Expansion(pGeom), StdExpansion2D()
        {
            m_elementFaceLeft  = -1;
            m_elementFaceRight = -1;
        }

        void Expansion2D::v_AddEdgeNormBoundaryInt(
            const int                           edge,
            const ExpansionSharedPtr           &EdgeExp,
            const Array<OneD, const NekDouble> &Fx,
            const Array<OneD, const NekDouble> &Fy,
                  Array<OneD,       NekDouble> &outarray)
        {
            ASSERTL1(GetCoordim() == 2,
                     "Routine only set up for two-dimensions");

            const Array<OneD, const Array<OneD, NekDouble> > normals
                                    = GetEdgeNormal(edge);

            if (m_requireNeg.size() == 0)
            {
                m_requireNeg.resize(GetNedges());

                for (int i = 0; i < GetNedges(); ++i)
                {
                    m_requireNeg[i] = false;

                    Expansion1DSharedPtr edgeExp =
                                    m_edgeExp[i].lock()->as<Expansion1D>();

                    if (edgeExp->GetRightAdjacentElementExp())
                    {
                        if (edgeExp->GetRightAdjacentElementExp()->GetGeom2D()
                            ->GetGlobalID() == GetGeom2D()->GetGlobalID())
                        {
                            m_requireNeg[i] = true;
                        }
                    }
                }
            }

            // We allow the case of mixed polynomial order by supporting only
            // those modes on the edge common to both adjoining elements. This
            // is enforced here by taking the minimum size and padding with
            // zeros.
            int nquad_e = min(EdgeExp->GetNumPoints(0),
                              int(normals[0].size()));

            int nEdgePts = EdgeExp->GetTotPoints();
            Array<OneD, NekDouble> edgePhys(nEdgePts);
            Vmath::Vmul (nquad_e, normals[0], 1, Fx, 1, edgePhys, 1);
            Vmath::Vvtvp(nquad_e, normals[1], 1, Fy, 1, edgePhys, 1,
                                  edgePhys,   1);

            Expansion1DSharedPtr locExp = EdgeExp->as<Expansion1D>();

            if (m_requireNeg[edge])
            {
                if (locExp->GetRightAdjacentElementExp()->GetGeom()->GetGlobalID()
                    == m_geom->GetGlobalID())
                {
                    Vmath::Neg(nquad_e, edgePhys, 1);
                }
            }

            AddEdgeNormBoundaryInt(edge, EdgeExp, edgePhys, outarray);
        }

        void Expansion2D::v_AddEdgeNormBoundaryInt(
            const int                           edge,
            const ExpansionSharedPtr           &EdgeExp,
            const Array<OneD, const NekDouble> &Fn,
                  Array<OneD,       NekDouble> &outarray)
        {
            int i;

            if (m_requireNeg.size() == 0)
            {
                m_requireNeg.resize(GetNedges());

                for (i = 0; i < GetNedges(); ++i)
                {
                    m_requireNeg[i] = false;

                    Expansion1DSharedPtr edgeExp =
                                m_edgeExp[i].lock()->as<Expansion1D>();

                    if (edgeExp->GetRightAdjacentElementExp())
                    {
                        if (edgeExp->GetRightAdjacentElementExp()->GetGeom2D()
                            ->GetGlobalID() == GetGeom2D()->GetGlobalID())
                        {
                            m_requireNeg[i] = true;
                        }
                    }
                }
            }

            StdRegions::IndexMapKey ikey(
                StdRegions::eEdgeToElement, DetShapeType(),
                GetBasisNumModes(0), GetBasisNumModes(1), 0,
                edge, GetEorient(edge));
            StdRegions::IndexMapValuesSharedPtr map =
                StdExpansion::GetIndexMap(ikey);

            // Order of the element
            int order_e = map->size();
            // Order of the trace
            int n_coeffs = EdgeExp->GetNcoeffs();

            Array<OneD, NekDouble> edgeCoeffs(n_coeffs);
            if(n_coeffs!=order_e) // Going to orthogonal space
            {
                EdgeExp->FwdTrans(Fn, edgeCoeffs);
                Expansion1DSharedPtr locExp = EdgeExp->as<Expansion1D>();

                if (m_requireNeg[edge])
                {
                    Vmath::Neg(n_coeffs, edgeCoeffs, 1);
                }

                Array<OneD, NekDouble> coeff(n_coeffs,0.0);
                LibUtilities::BasisType btype = ((LibUtilities::BasisType) 1); //1-->Ortho_A
                LibUtilities::BasisKey bkey_ortho(btype,EdgeExp->GetBasis(0)->GetNumModes(),EdgeExp->GetBasis(0)->GetPointsKey());
                LibUtilities::BasisKey bkey(EdgeExp->GetBasis(0)->GetBasisType(),EdgeExp->GetBasis(0)->GetNumModes(),EdgeExp->GetBasis(0)->GetPointsKey());
                LibUtilities::InterpCoeff1D(bkey,edgeCoeffs,bkey_ortho,coeff);

                // Cutting high frequencies
                for(i = order_e; i < n_coeffs; i++)
                {
                    coeff[i] = 0.0;
                }

                LibUtilities::InterpCoeff1D(bkey_ortho,coeff,bkey,edgeCoeffs);

                StdRegions::StdMatrixKey masskey(StdRegions::eMass,LibUtilities::eSegment,*EdgeExp);
                EdgeExp->MassMatrixOp(edgeCoeffs, edgeCoeffs, masskey);
            }
            else
            {
                EdgeExp->IProductWRTBase(Fn, edgeCoeffs);

                Expansion1DSharedPtr locExp = EdgeExp->as<Expansion1D>();

                if (m_requireNeg[edge])
                {
                    Vmath::Neg(n_coeffs, edgeCoeffs, 1);
                }
            }

            // Implementation for all the basis except Gauss points
            if(EdgeExp->GetBasis(0)->GetBasisType() !=
                   LibUtilities::eGauss_Lagrange)
            {
                // add data to outarray if forward edge normal is outwards
                for(i = 0; i < order_e; ++i)
                {
                    outarray[(*map)[i].index] +=
                        (*map)[i].sign * edgeCoeffs[i];
                }
            }
            else
            {
                int nCoeffs0, nCoeffs1;
                int j;

                StdRegions::ConstFactorMap factors;
                factors[StdRegions::eFactorGaussEdge] = edge;
                StdRegions::StdMatrixKey key(StdRegions::eGaussDG,
                                             DetShapeType(),*this,factors);

                DNekMatSharedPtr mat_gauss = m_stdMatrixManager[key];

                switch(edge)
                {
                    case 0:
                    {
                        nCoeffs1 = m_base[1]->GetNumModes();

                        for(i = 0; i < order_e; ++i)
                        {
                            for(j = 0; j < nCoeffs1; j++)
                            {
                                outarray[(*map)[i].index + j*order_e] +=
                                    mat_gauss->GetPtr()[j]*
                                    (*map)[i].sign*edgeCoeffs[i];
                            }
                        }
                        break;
                    }
                    case 1:
                    {
                        nCoeffs0 = m_base[0]->GetNumModes();

                        for(i = 0; i < order_e; ++i)
                        {
                            for(j = 0; j < nCoeffs0; j++)
                            {
                                outarray[(*map)[i].index - j] +=
                                    mat_gauss->GetPtr()[order_e - 1 -j]*
                                    (*map)[i].sign*edgeCoeffs[i];
                            }
                        }
                        break;
                    }
                    case 2:
                    {
                        nCoeffs1 = m_base[1]->GetNumModes();

                        for(i = 0; i < order_e; ++i)
                        {
                            for(j = 0; j < nCoeffs1; j++)
                            {
                                outarray[(*map)[i].index - j*order_e] +=
                                    mat_gauss->GetPtr()[order_e - 1 - j]*
                                    (*map)[i].sign*edgeCoeffs[i];
                            }
                        }
                        break;
                    }
                    case 3:
                    {
                        nCoeffs0 = m_base[0]->GetNumModes();

                        for(i = 0; i < order_e; ++i)
                        {
                            for(j = 0; j < nCoeffs0; j++)
                            {
                                outarray[(*map)[i].index + j] +=
                                    mat_gauss->GetPtr()[j]*
                                    (*map)[i].sign*edgeCoeffs[i];
                            }
                        }
                        break;
                    }
                    default:
                        ASSERTL0(false,"edge value (< 3) is out of range");
                        break;
                }
            }
        }

        void Expansion2D::SetTraceToGeomOrientation(
            Array<OneD, ExpansionSharedPtr> &EdgeExp,
            Array<OneD, NekDouble> &inout)
        {
            int i, cnt = 0;
            int nedges = GetNedges();
            Array<OneD, NekDouble> e_tmp;

            for(i = 0; i < nedges; ++i)
            {
                EdgeExp[i]->SetCoeffsToOrientation(GetEorient(i),
                                                   e_tmp = inout + cnt,
                                                   e_tmp = inout + cnt);
                cnt += GetEdgeNcoeffs(i);
            }
        }

        /**
         * Computes the C matrix entries due to the presence of the identity
         * matrix in Eqn. 32.
         */
        void Expansion2D::AddNormTraceInt(
            const int                        dir,
            Array<OneD, const NekDouble>    &inarray,
            Array<OneD, ExpansionSharedPtr> &EdgeExp,
            Array<OneD, NekDouble>          &outarray,
            const StdRegions::VarCoeffMap   &varcoeffs)
        {
            int i,e,cnt;
            int order_e,nquad_e;
            int nedges = GetNedges();

            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                nquad_e = EdgeExp[e]->GetNumPoints(0);

                const Array<OneD, const Array<OneD, NekDouble> > &normals
                    = GetEdgeNormal(e);
                Array<OneD, NekDouble> edgeCoeffs(order_e);
                Array<OneD, NekDouble> edgePhys  (nquad_e);

                for(i = 0; i < order_e; ++i)
                {
                    edgeCoeffs[i] = inarray[i+cnt];
                }
                cnt += order_e;

                EdgeExp[e]->BwdTrans(edgeCoeffs, edgePhys);

                // Multiply by variable coefficient
                /// @TODO: Document this
                // StdRegions::VarCoeffType VarCoeff[3] = {StdRegions::eVarCoeffD00,
                //                                         StdRegions::eVarCoeffD11,
                //                                         StdRegions::eVarCoeffD22};
                // StdRegions::VarCoeffMap::const_iterator x;
                // Array<OneD, NekDouble> varcoeff_work(nquad_e);

                // if ((x = varcoeffs.find(VarCoeff[dir])) != varcoeffs.end())
                // {
                //     GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                //     Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
                // }

                if (varcoeffs.find(StdRegions::eVarCoeffMF1x) !=
                        varcoeffs.end())
                {
                    // MMF case
                    Array<OneD, NekDouble> ncdotMF_e =
                            v_GetnEdgecdotMF(dir, e, EdgeExp[e], normals,
                                             varcoeffs);

                    Vmath::Vmul(nquad_e, ncdotMF_e, 1,
                                         edgePhys,  1,
                                         edgePhys,  1);
                }
                else
                {
                    Vmath::Vmul(nquad_e, normals[dir], 1,
                                         edgePhys,     1,
                                         edgePhys,     1);
                }

                AddEdgeBoundaryInt(e, EdgeExp[e], edgePhys, outarray, varcoeffs);
            }
        }

        void Expansion2D::AddNormTraceInt(
            const int                             dir,
            Array<OneD, ExpansionSharedPtr>      &EdgeExp,
            Array<OneD, Array<OneD, NekDouble> > &edgeCoeffs,
            Array<OneD, NekDouble>               &outarray)
        {
            int e;
            int nquad_e;
            int nedges = GetNedges();

            for(e = 0; e < nedges; ++e)
            {
                nquad_e = EdgeExp[e]->GetNumPoints(0);

                Array<OneD, NekDouble> edgePhys(nquad_e);
                const Array<OneD, const Array<OneD, NekDouble> > &normals
                    = GetEdgeNormal(e);

                EdgeExp[e]->BwdTrans(edgeCoeffs[e], edgePhys);

                Vmath::Vmul(nquad_e, normals[dir], 1, edgePhys, 1, edgePhys, 1);

                AddEdgeBoundaryInt(e, EdgeExp[e], edgePhys, outarray);
            }
        }


        /**
         * For a given edge add the \tilde{F}_1j contributions
         */
        void Expansion2D::AddEdgeBoundaryInt(
            const int                          edge,
            ExpansionSharedPtr                &EdgeExp,
            Array<OneD, NekDouble>            &edgePhys,
            Array<OneD, NekDouble>            &outarray,
            const StdRegions::VarCoeffMap     &varcoeffs)
        {
            int i;
            int order_e = EdgeExp->GetNcoeffs();
            int nquad_e = EdgeExp->GetNumPoints(0);
            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;
            Array<OneD, NekDouble> coeff(order_e);

            GetEdgeToElementMap(edge, v_GetEorient(edge), map, sign);

            StdRegions::VarCoeffType VarCoeff[3] = {StdRegions::eVarCoeffD00,
                                                    StdRegions::eVarCoeffD11,
                                                    StdRegions::eVarCoeffD22};
            StdRegions::VarCoeffMap::const_iterator x;

            /// @TODO Variable coeffs
            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
            {
                Array<OneD, NekDouble> work(nquad_e);
                GetPhysEdgeVarCoeffsFromElement(
                    edge, EdgeExp, x->second, work);
                Vmath::Vmul(nquad_e, work, 1, edgePhys, 1, edgePhys, 1);
            }

            EdgeExp->IProductWRTBase(edgePhys, coeff);

            // add data to out array
            for(i = 0; i < order_e; ++i)
            {
                outarray[map[i]] += sign[i]*coeff[i];
            }
        }

        // This method assumes that data in EdgeExp is orientated according to
        // elemental counter clockwise format AddHDGHelmholtzTraceTerms with
        // directions
        void Expansion2D::AddHDGHelmholtzTraceTerms(
            const NekDouble                     tau,
            const Array<OneD, const NekDouble> &inarray,
            Array<OneD, ExpansionSharedPtr>    &EdgeExp,
            const StdRegions::VarCoeffMap      &dirForcing,
            Array<OneD,NekDouble>              &outarray)
        {
            ASSERTL0(&inarray[0] != &outarray[0],
                     "Input and output arrays use the same memory");

            int e, cnt, order_e, nedges = GetNedges();
            Array<OneD, const NekDouble> tmp;

            cnt = 0;

            for(e = 0; e < nedges; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                Array<OneD, NekDouble> edgeCoeffs(order_e);
                Array<OneD, NekDouble> edgePhys  (EdgeExp[e]->GetTotPoints());

                Vmath::Vcopy(order_e, tmp = inarray + cnt, 1, edgeCoeffs, 1);
                EdgeExp[e]->BwdTrans(edgeCoeffs, edgePhys);
                AddHDGHelmholtzEdgeTerms(
                    tau, e, EdgeExp, edgePhys, dirForcing, outarray);

                cnt += order_e;
            }
        }

        // evaluate additional terms in HDG edges. Not that this assumes that
        // edges are unpacked into local cartesian order.
        void Expansion2D::AddHDGHelmholtzEdgeTerms(
            const NekDouble                  tau,
            const int                        edge,
            Array<OneD, ExpansionSharedPtr> &EdgeExp,
            Array<OneD, NekDouble>          &edgePhys,
            const StdRegions::VarCoeffMap   &varcoeffs,
            Array<OneD, NekDouble>          &outarray)
        {
            bool mmf = (varcoeffs.find(StdRegions::eVarCoeffMF1x) !=
                        varcoeffs.end());
            int i, j, n;
            int nquad_e = EdgeExp[edge]->GetNumPoints(0);
            int order_e = EdgeExp[edge]->GetNcoeffs();
            int coordim = mmf ? 2 : GetCoordim();
            int ncoeffs = GetNcoeffs();

            Array<OneD, NekDouble> inval   (nquad_e);
            Array<OneD, NekDouble> outcoeff(order_e);
            Array<OneD, NekDouble> tmpcoeff(ncoeffs);

            const Array<OneD, const Array<OneD, NekDouble> > &normals
                                = GetEdgeNormal(edge);

            Array<OneD,unsigned int> emap;
            Array<OneD,int> sign;

            StdRegions::Orientation edgedir = GetEorient(edge);

            DNekVec                Coeffs  (ncoeffs,outarray,eWrapper);
            DNekVec                Tmpcoeff(ncoeffs,tmpcoeff,eWrapper);

            GetEdgeToElementMap(edge,edgedir,emap,sign);

            StdRegions::MatrixType DerivType[3] =
            {
                StdRegions::eWeakDeriv0,
                StdRegions::eWeakDeriv1,
                StdRegions::eWeakDeriv2
            };

            StdRegions::VarCoeffType VarCoeff[3] =
            {
                StdRegions::eVarCoeffD00,
                StdRegions::eVarCoeffD11,
                StdRegions::eVarCoeffD22
            };

            StdRegions::VarCoeffMap::const_iterator x;
            /// @TODO: What direction to use here??
            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
            {
                Array<OneD, NekDouble> work(nquad_e);
                GetPhysEdgeVarCoeffsFromElement(
                    edge, EdgeExp[edge], x->second, work);
                Vmath::Vmul(nquad_e, work, 1, edgePhys, 1, edgePhys, 1);
            }

            //================================================================
            // Add F = \tau <phi_i,in_phys>
            // Fill edge and take inner product
            EdgeExp[edge]->IProductWRTBase(edgePhys, outcoeff);
            // add data to out array
            for(i = 0; i < order_e; ++i)
            {
                outarray[emap[i]] += sign[i] * tau * outcoeff[i];
            }
            //================================================================

            //===============================================================
            // Add -\sum_i D_i^T M^{-1} G_i + E_i M^{-1} G_i =
            //                         \sum_i D_i M^{-1} G_i term

            // Two independent direction
            DNekScalMatSharedPtr invMass;
            for(n = 0; n < coordim; ++n)
            {
                if (mmf)
                {
                    StdRegions::VarCoeffMap Weight;
                    Weight[StdRegions::eVarCoeffMass] = v_GetMFMag(n,varcoeffs);

                    MatrixKey invMasskey( StdRegions::eInvMass,
                                          DetShapeType(), *this,
                                          StdRegions::NullConstFactorMap,
                                          Weight);

                    invMass = GetLocMatrix(invMasskey);

                    Array<OneD, NekDouble> ncdotMF_e =
                            v_GetnEdgecdotMF(n, edge, EdgeExp[edge], normals,
                                             varcoeffs);

                    Vmath::Vmul(nquad_e, ncdotMF_e, 1, edgePhys, 1, inval, 1);
                }
                else
                {
                    Vmath::Vmul(nquad_e, normals[n], 1, edgePhys, 1, inval, 1);
                    invMass = GetLocMatrix(StdRegions::eInvMass);
                }

                // Multiply by variable coefficient
                /// @TODO: Document this (probably not needed)
//                StdRegions::VarCoeffMap::const_iterator x;
//                if ((x = varcoeffs.find(VarCoeff[n])) != varcoeffs.end())
//                {
//                    GetPhysEdgeVarCoeffsFromElement(edge,EdgeExp[edge],x->second,varcoeff_work);
//                    Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[edge]->GetPhys(),1,EdgeExp[edge]->UpdatePhys(),1);
//                }

                EdgeExp[edge]->IProductWRTBase(inval, outcoeff);

                // M^{-1} G
                for(i = 0; i < ncoeffs; ++i)
                {
                    tmpcoeff[i] = 0;
                    for(j = 0; j < order_e; ++j)
                    {
                        tmpcoeff[i] += (*invMass)(i,emap[j]) * sign[j]
                                                             * outcoeff[j];
                    }
                }

                if (mmf)
                {
                    StdRegions::VarCoeffMap VarCoeffDirDeriv;
                    VarCoeffDirDeriv[StdRegions::eVarCoeffMF] =
                            v_GetMF(n,coordim,varcoeffs);
                    VarCoeffDirDeriv[StdRegions::eVarCoeffMFDiv] =
                            v_GetMFDiv(n,varcoeffs);

                    MatrixKey Dmatkey( StdRegions::eWeakDirectionalDeriv,
                                       DetShapeType(), *this,
                                       StdRegions::NullConstFactorMap,
                                       VarCoeffDirDeriv);

                    DNekScalMat &Dmat = *GetLocMatrix(Dmatkey);

                    Coeffs = Coeffs  + Dmat*Tmpcoeff;
                }
                else
                {
                    if(varcoeffs.find(VarCoeff[n]) != varcoeffs.end())
                    {
                        MatrixKey mkey(DerivType[n], DetShapeType(), *this,
                                StdRegions::NullConstFactorMap, varcoeffs);

                        DNekScalMat &Dmat = *GetLocMatrix(mkey);
                        Coeffs = Coeffs  + Dmat*Tmpcoeff;
                    }
                    else
                    {
                        DNekScalMat &Dmat = *GetLocMatrix(DerivType[n]);
                        Coeffs = Coeffs  + Dmat*Tmpcoeff;
                    }
                }
            }
        }


        /**
         * Extracts the variable coefficients along an edge
         */
        void Expansion2D::GetPhysEdgeVarCoeffsFromElement(
            const int                           edge,
            ExpansionSharedPtr                 &EdgeExp,
            const Array<OneD, const NekDouble> &varcoeff,
            Array<OneD,NekDouble>              &outarray)
        {
            Array<OneD, NekDouble> tmp(GetNcoeffs());
            Array<OneD, NekDouble> edgetmp(EdgeExp->GetNcoeffs());

            // FwdTrans varcoeffs
            FwdTrans(varcoeff, tmp);

            // Map to edge
            Array<OneD,unsigned int>    emap;
            Array<OneD, int>            sign;
            StdRegions::Orientation edgedir = GetEorient(edge);
            GetEdgeToElementMap(edge,edgedir,emap,sign);

            for (unsigned int i = 0; i < EdgeExp->GetNcoeffs(); ++i)
            {
                edgetmp[i] = tmp[emap[i]];
            }

            // BwdTrans
            EdgeExp->BwdTrans(edgetmp, outarray);
        }


        /**
         * Computes matrices needed for the HDG formulation. References to
         * equations relate to the following paper:
         *   R. M. Kirby, S. J. Sherwin, B. Cockburn, To CG or to HDG: A
         *   Comparative Study, J. Sci. Comp P1-30
         *   DOI 10.1007/s10915-011-9501-7
         */
        DNekMatSharedPtr Expansion2D::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            // (Z^e)^{-1} (Eqn. 33, P22)
            case StdRegions::eHybridDGHelmholtz:
                {
                    ASSERTL1(IsBoundaryInteriorExpansion(),
                             "HybridDGHelmholtz matrix not set up "
                             "for non boundary-interior expansions");

                    int       i,j,k;
                    NekDouble lambdaval = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    NekDouble tau       = mkey.GetConstFactor(StdRegions::eFactorTau);
                    int       ncoeffs   = GetNcoeffs();
                    int       nedges    = GetNedges();
                    int       shapedim  = 2;
                    const StdRegions::VarCoeffMap &varcoeffs
                                        = mkey.GetVarCoeffs();
                    bool      mmf       =
                        (varcoeffs.find(StdRegions::eVarCoeffMF1x) !=
                         varcoeffs.end());

                    Array<OneD,unsigned int> emap;
                    Array<OneD,int> sign;
                    StdRegions::Orientation edgedir = StdRegions::eForwards;
                    ExpansionSharedPtr EdgeExp;
                    ExpansionSharedPtr EdgeExp2;

                    int order_e, coordim = GetCoordim();
                    DNekScalMat  &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                           StdRegions::eWeakDeriv1,
                                                           StdRegions::eWeakDeriv2};
                    DNekMat LocMat(ncoeffs,ncoeffs);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,ncoeffs);
                    DNekMat &Mat = *returnval;
                    Vmath::Zero(ncoeffs*ncoeffs,Mat.GetPtr(),1);

                    StdRegions::VarCoeffType Coeffs[3] = {StdRegions::eVarCoeffD00,
                                                          StdRegions::eVarCoeffD11,
                                                          StdRegions::eVarCoeffD22};

                    StdRegions::VarCoeffMap::const_iterator x;

                    for(i=0;  i < coordim; ++i)
                    {
                        if (mmf)
                        {
                            if(i < shapedim)
                            {
                                StdRegions::VarCoeffMap VarCoeffDirDeriv;
                                VarCoeffDirDeriv[StdRegions::eVarCoeffMF] =
                                        v_GetMF(i,shapedim,varcoeffs);
                                VarCoeffDirDeriv[StdRegions::eVarCoeffMFDiv] =
                                        v_GetMFDiv(i,varcoeffs);

                                MatrixKey Dmatkey(StdRegions::
                                                        eWeakDirectionalDeriv,
                                                  DetShapeType(), *this,
                                                  StdRegions::
                                                        NullConstFactorMap,
                                                  VarCoeffDirDeriv);

                                DNekScalMat &Dmat = *GetLocMatrix(Dmatkey);

                                StdRegions::VarCoeffMap Weight;
                                Weight[StdRegions::eVarCoeffMass] =
                                        v_GetMFMag(i,mkey.GetVarCoeffs());

                                MatrixKey invMasskey( StdRegions::eInvMass,
                                                      DetShapeType(), *this,
                                                      StdRegions::
                                                            NullConstFactorMap,
                                                      Weight);

                                DNekScalMat &invMass =
                                        *GetLocMatrix(invMasskey);

                                Mat = Mat + Dmat*invMass*Transpose(Dmat);
                            }
                        }
                        else if(mkey.HasVarCoeff(Coeffs[i]))
                        {
                            MatrixKey DmatkeyL( DerivType[i],
                                                DetShapeType(), *this,
                                                StdRegions::NullConstFactorMap,
                                                mkey.GetVarCoeffAsMap(
                                                        Coeffs[i]));

                            MatrixKey DmatkeyR( DerivType[i],
                                                DetShapeType(), *this);

                            DNekScalMat &DmatL = *GetLocMatrix(DmatkeyL);
                            DNekScalMat &DmatR = *GetLocMatrix(DmatkeyR);
                            Mat = Mat + DmatL*invMass*Transpose(DmatR);
                        }
                        else
                        {
                            DNekScalMat &Dmat = *GetLocMatrix(DerivType[i]);
                            Mat = Mat + Dmat*invMass*Transpose(Dmat);
                        }

                    }

                    // Add Mass Matrix Contribution for Helmholtz problem
                    DNekScalMat  &Mass = *GetLocMatrix(StdRegions::eMass);
                    Mat = Mat + lambdaval*Mass;

                    // Add tau*E_l using elemental mass matrices on each edge
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp = GetEdgeExp(i);
                        EdgeExp2 = GetEdgeExp(i);
                        order_e = EdgeExp->GetNcoeffs();
                        int nq = EdgeExp->GetNumPoints(0);
                        GetEdgeToElementMap(i,edgedir,emap,sign);

                        // @TODO: Document
                        StdRegions::VarCoeffMap edgeVarCoeffs;
                        if (mkey.HasVarCoeff(StdRegions::eVarCoeffD00))
                        {
                            Array<OneD, NekDouble> mu(nq);
                            GetPhysEdgeVarCoeffsFromElement(i, EdgeExp2, mkey.GetVarCoeff(StdRegions::eVarCoeffD00), mu);
                            edgeVarCoeffs[StdRegions::eVarCoeffMass] = mu;
                        }
                        DNekScalMat &eMass = *EdgeExp->GetLocMatrix(StdRegions::eMass, StdRegions::NullConstFactorMap, edgeVarCoeffs);
                        //DNekScalMat &eMass = *EdgeExp->GetLocMatrix(StdRegions::eMass);

                        for(j = 0; j < order_e; ++j)
                        {
                            for(k = 0; k < order_e; ++k)
                            {
                                Mat(emap[j],emap[k]) = Mat(emap[j],emap[k]) + tau*sign[j]*sign[k]*eMass(j,k);
                            }
                        }
                    }
                }
                break;
            // U^e (P22)
            case StdRegions::eHybridDGLamToU:
                {
                    int       i,j,k;
                    int       nbndry  = NumDGBndryCoeffs();
                    int       ncoeffs = GetNcoeffs();
                    int       nedges  = GetNedges();
                    NekDouble tau     = mkey.GetConstFactor(StdRegions::eFactorTau);

                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);

                    Array<OneD, ExpansionSharedPtr>  EdgeExp(nedges);
                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry);
                    DNekMat &Umat = *returnval;

                    // Z^e matrix
                    MatrixKey newkey(StdRegions::eInvHybridDGHelmholtz, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat  &invHmat = *GetLocMatrix(newkey);

                    Array<OneD,unsigned int> emap;
                    Array<OneD,int> sign;

                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = GetEdgeExp(i);
                    }

                    // for each degree of freedom of the lambda space
                    // calculate Umat entry
                    // Generate Lambda to U_lambda matrix
                    for(j = 0; j < nbndry; ++j)
                    {
                        // standard basis vectors e_j
                        Vmath::Zero(nbndry,&lambda[0],1);
                        Vmath::Zero(ncoeffs,&f[0],1);
                        lambda[j] = 1.0;

                        SetTraceToGeomOrientation(EdgeExp,lambda);

                        // Compute F = [I   D_1 M^{-1}   D_2 M^{-1}] C e_j
                        AddHDGHelmholtzTraceTerms(
                            tau, lambda, EdgeExp, mkey.GetVarCoeffs(), f);

                        // Compute U^e_j
                        Ulam = invHmat*F; // generate Ulam from lambda

                        // fill column of matrix
                        for(k = 0; k < ncoeffs; ++k)
                        {
                            Umat(k,j) = Ulam[k];
                        }
                    }
                }
                break;
            // Q_0, Q_1, Q_2 matrices (P23)
            // Each are a product of a row of Eqn 32 with the C matrix.
            // Rather than explicitly computing all of Eqn 32, we note each
            // row is almost a multiple of U^e, so use that as our starting
            // point.
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
                {
                    int i        = 0;
                    int j        = 0;
                    int k        = 0;
                    int dir      = 0;
                    int nbndry   = NumDGBndryCoeffs();
                    int ncoeffs  = GetNcoeffs();
                    int nedges   = GetNedges();
                    int shapedim = 2;

                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);
                    Array<OneD, ExpansionSharedPtr>  EdgeExp(nedges);

                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);

                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry);
                    DNekMat &Qmat = *returnval;

                    // Lambda to U matrix
                    MatrixKey lamToUkey(StdRegions::eHybridDGLamToU, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &lamToU = *GetLocMatrix(lamToUkey);

                    // Inverse mass matrix
                    DNekScalMat &invMass = *GetLocMatrix(StdRegions::eInvMass);


                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = GetEdgeExp(i);
                    }

                    //Weak Derivative matrix
                    DNekScalMatSharedPtr Dmat;
                    switch(mkey.GetMatrixType())
                    {
                    case StdRegions::eHybridDGLamToQ0:
                        dir = 0;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv0);
                        break;
                    case StdRegions::eHybridDGLamToQ1:
                        dir = 1;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv1);
                        break;
                    case StdRegions::eHybridDGLamToQ2:
                        dir = 2;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv2);
                        break;
                    default:
                        ASSERTL0(false,"Direction not known");
                        break;
                    }

                    const StdRegions::VarCoeffMap &varcoeffs =
                            mkey.GetVarCoeffs();
                    if (varcoeffs.find(StdRegions::eVarCoeffMF1x) !=
                            varcoeffs.end())
                    {
                        StdRegions::VarCoeffMap VarCoeffDirDeriv;
                        VarCoeffDirDeriv[StdRegions::eVarCoeffMF] =
                                v_GetMF(dir,shapedim,varcoeffs);
                        VarCoeffDirDeriv[StdRegions::eVarCoeffMFDiv] =
                                v_GetMFDiv(dir,varcoeffs);

                        MatrixKey Dmatkey( StdRegions::eWeakDirectionalDeriv,
                                           DetShapeType(), *this,
                                           StdRegions::NullConstFactorMap,
                                           VarCoeffDirDeriv);

                        Dmat = GetLocMatrix(Dmatkey);

                        StdRegions::VarCoeffMap Weight;
                        Weight[StdRegions::eVarCoeffMass] =
                                v_GetMFMag(dir,mkey.GetVarCoeffs());

                        MatrixKey invMasskey( StdRegions::eInvMass,
                                              DetShapeType(), *this,
                                              StdRegions::NullConstFactorMap,
                                              Weight);

                        invMass = *GetLocMatrix(invMasskey);
                    }
                    else
                    {
                        StdRegions::MatrixType DerivType[3] = {
                            StdRegions::eWeakDeriv0,
                            StdRegions::eWeakDeriv1,
                            StdRegions::eWeakDeriv2};

                        Dmat = GetLocMatrix(DerivType[dir]);

                        MatrixKey invMasskey( StdRegions::eInvMass,
                                              DetShapeType(), *this);
                        invMass = *GetLocMatrix(invMasskey);
                    }

                    // for each degree of freedom of the lambda space
                    // calculate Qmat entry
                    // Generate Lambda to Q_lambda matrix
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        lambda[j] = 1.0;

                        // for lambda[j] = 1 this is the solution to ulam
                        for(k = 0; k < ncoeffs; ++k)
                        {
                            Ulam[k] = lamToU(k,j);
                        }

                        // -D^T ulam
                        Vmath::Neg(ncoeffs,&ulam[0],1);
                        F = Transpose(*Dmat)*Ulam;

                        SetTraceToGeomOrientation(EdgeExp,lambda);

                        // Add the C terms resulting from the I's on the
                        // diagonals of Eqn 32
                        AddNormTraceInt(dir,lambda,EdgeExp,f,mkey.GetVarCoeffs());

                        // finally multiply by inverse mass matrix
                        Ulam = invMass*F;

                        // fill column of matrix (Qmat is in column major format)
                        Vmath::Vcopy(ncoeffs,&ulam[0],1,&(Qmat.GetPtr())[0]+j*ncoeffs,1);
                    }
                }
                break;
                        // Matrix K (P23)
            case StdRegions::eHybridDGHelmBndLam:
                {
                    int i,j,e,cnt;
                    int order_e, nquad_e;
                    int nbndry  = NumDGBndryCoeffs();
                    int coordim = GetCoordim();
                    int nedges  = GetNedges();
                    NekDouble tau = mkey.GetConstFactor(StdRegions::eFactorTau);
                    StdRegions::VarCoeffMap::const_iterator x;
                    const StdRegions::VarCoeffMap &varcoeffs
                                = mkey.GetVarCoeffs();
                    bool mmf    =
                        (varcoeffs.find(StdRegions::eVarCoeffMF1x) !=
                         varcoeffs.end());

                    Array<OneD,NekDouble>       work, varcoeff_work;
                    Array<OneD,const Array<OneD, NekDouble> > normals;
                    Array<OneD, ExpansionSharedPtr>  EdgeExp(nedges);
                    Array<OneD, NekDouble> lam(nbndry);

                    Array<OneD,unsigned int>    emap;
                    Array<OneD, int>            sign;
                    StdRegions::Orientation edgedir;

                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry, nbndry);
                    DNekMat &BndMat = *returnval;

                    DNekScalMatSharedPtr LamToQ[3];

                    // Matrix to map Lambda to U
                    MatrixKey LamToUkey(StdRegions::eHybridDGLamToU, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LamToU = *GetLocMatrix(LamToUkey);

                    // Matrix to map Lambda to Q0
                    MatrixKey LamToQ0key(StdRegions::eHybridDGLamToQ0, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    LamToQ[0] = GetLocMatrix(LamToQ0key);

                    // Matrix to map Lambda to Q1
                    MatrixKey LamToQ1key(StdRegions::eHybridDGLamToQ1, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    LamToQ[1] = GetLocMatrix(LamToQ1key);

                    // Matrix to map Lambda to Q2 for 3D coordinates
                    if (coordim == 3)
                    {
                        MatrixKey LamToQ2key(StdRegions::eHybridDGLamToQ2, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                        LamToQ[2] = GetLocMatrix(LamToQ2key);
                    }

                    // Set up edge segment expansions from local geom info
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = GetEdgeExp(i);
                    }

                    // Set up matrix derived from <mu, Q_lam.n - \tau (U_lam - Lam) >
                    for(i = 0; i < nbndry; ++i)
                    {
                        cnt = 0;

                        Vmath::Zero(nbndry,lam,1);
                        lam[i] = 1.0;
                        SetTraceToGeomOrientation(EdgeExp,lam);

                        for(e = 0; e < nedges; ++e)
                        {
                            order_e = EdgeExp[e]->GetNcoeffs();
                            nquad_e = EdgeExp[e]->GetNumPoints(0);

                            normals = GetEdgeNormal(e);
                            edgedir = GetEorient(e);

                            work = Array<OneD,NekDouble>(nquad_e);
                            varcoeff_work = Array<OneD, NekDouble>(nquad_e);

                            GetEdgeToElementMap(e,edgedir,emap,sign);


                            StdRegions::VarCoeffType VarCoeff[3] = {StdRegions::eVarCoeffD00,
                                                                    StdRegions::eVarCoeffD11,
                                                                    StdRegions::eVarCoeffD22};

                            // Q0 * n0 (BQ_0 terms)
                            Array<OneD, NekDouble> edgeCoeffs(order_e);
                            Array<OneD, NekDouble> edgePhys  (nquad_e);
                            for(j = 0; j < order_e; ++j)
                            {
                                edgeCoeffs[j] = sign[j]*(*LamToQ[0])(emap[j],i);
                            }

                            EdgeExp[e]->BwdTrans(edgeCoeffs, edgePhys);
// @TODO Var coeffs
                            // Multiply by variable coefficient
//                            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
//                            {
//                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
//                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
//                            }
                            if (mmf)
                            {
                                Array<OneD, NekDouble> ncdotMF =
                                        v_GetnEdgecdotMF(0, e, EdgeExp[e],
                                                         normals, varcoeffs);
                                Vmath::Vmul(nquad_e, ncdotMF,  1,
                                                     edgePhys, 1,
                                                     work,     1);
                            }
                            else
                            {
                                Vmath::Vmul(nquad_e, normals[0], 1,
                                                     edgePhys,   1,
                                                     work,       1);
                            }

                            // Q1 * n1 (BQ_1 terms)
                            for(j = 0; j < order_e; ++j)
                            {
                                edgeCoeffs[j] = sign[j]*(*LamToQ[1])(emap[j],i);
                            }

                            EdgeExp[e]->BwdTrans(edgeCoeffs, edgePhys);

// @TODO var coeffs
                            // Multiply by variable coefficients
//                            if ((x = varcoeffs.find(VarCoeff[1])) != varcoeffs.end())
//                            {
//                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
//                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
//                            }

                            if (mmf)
                            {
                                Array<OneD, NekDouble> ncdotMF =
                                        v_GetnEdgecdotMF(1, e, EdgeExp[e],
                                                         normals, varcoeffs);
                                Vmath::Vvtvp(nquad_e, ncdotMF,  1,
                                                      edgePhys, 1,
                                                      work,     1,
                                                      work,     1);
                            }
                            else
                            {
                                Vmath::Vvtvp(nquad_e, normals[1], 1,
                                                      edgePhys,   1,
                                                      work,       1,
                                                      work,       1);
                            }

                            // Q2 * n2 (BQ_2 terms)
                            if (coordim == 3)
                            {
                                for(j = 0; j < order_e; ++j)
                                {
                                    edgeCoeffs[j] = sign[j]*(*LamToQ[2])(emap[j],i);
                                }

                                EdgeExp[e]->BwdTrans(edgeCoeffs, edgePhys);
// @TODO var coeffs
                                // Multiply by variable coefficients
//                                if ((x = varcoeffs.find(VarCoeff[2])) != varcoeffs.end())
//                                {
//                                    GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
//                                    Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
//                                }

                                Vmath::Vvtvp(nquad_e, normals[2], 1, edgePhys, 1,
                                                      work,       1, work,     1);
                            }

                            // - tau (ulam - lam)
                            // Corresponds to the G and BU terms.
                            for(j = 0; j < order_e; ++j)
                            {
                                edgeCoeffs[j] = sign[j]*LamToU(emap[j],i) - lam[cnt+j];
                            }

                            EdgeExp[e]->BwdTrans(edgeCoeffs, edgePhys);

                            // Multiply by variable coefficients
                            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_e,varcoeff_work,1,edgePhys,1,edgePhys,1);
                            }

                            Vmath::Svtvp(nquad_e,-tau,edgePhys,1,
                                         work,1,work,1);
/// TODO: Add variable coeffs
                            EdgeExp[e]->IProductWRTBase(work, edgeCoeffs);

                            EdgeExp[e]->SetCoeffsToOrientation(edgeCoeffs, edgedir);

                            for(j = 0; j < order_e; ++j)
                            {
                                BndMat(cnt+j,i) = edgeCoeffs[j];
                            }

                            cnt += order_e;
                        }
                    }
                }
                break;
            //HDG postprocessing
            case StdRegions::eInvLaplacianWithUnityMean:
                {
                    MatrixKey lapkey(StdRegions::eLaplacian, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat  &LapMat = *GetLocMatrix(lapkey);

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(LapMat.GetRows(),LapMat.GetColumns());
                    DNekMatSharedPtr lmat = returnval;

                    (*lmat) = LapMat;

                    // replace first column with inner product wrt 1
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmp(nq);
                    Array<OneD, NekDouble> outarray(m_ncoeffs);
                    Vmath::Fill(nq,1.0,tmp,1);
                    IProductWRTBase(tmp, outarray);

                    Vmath::Vcopy(m_ncoeffs,&outarray[0],1,
                                 &(lmat->GetPtr())[0],1);

                    lmat->Invert();
                }
                break;
            default:
                ASSERTL0(false,"This matrix type cannot be generated from this class");
                break;
            }

            return returnval;
        }

        //Evaluate Coefficients of weak deriviative in the direction dir
        //given the input coefficicents incoeffs and the imposed
        //boundary values in EdgeExp (which will have its phys space updated);
        void Expansion2D::v_DGDeriv(
            int                                   dir,
            const Array<OneD, const NekDouble>   &incoeffs,
            Array<OneD, ExpansionSharedPtr>      &EdgeExp,
            Array<OneD, Array<OneD, NekDouble> > &edgeCoeffs,
            Array<OneD, NekDouble>               &out_d)
        {
            StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                   StdRegions::eWeakDeriv1,
                                                   StdRegions::eWeakDeriv2};

            int ncoeffs = GetNcoeffs();

            DNekScalMat &InvMass = *GetLocMatrix(StdRegions::eInvMass);
            DNekScalMat &Dmat    = *GetLocMatrix(DerivType[dir]);

            Array<OneD, NekDouble> coeffs = incoeffs;
            DNekVec     Coeffs  (ncoeffs,coeffs, eWrapper);

            Coeffs = Transpose(Dmat)*Coeffs;
            Vmath::Neg(ncoeffs, coeffs,1);

            // Add the boundary integral including the relevant part of
            // the normal
            AddNormTraceInt(dir, EdgeExp, edgeCoeffs, coeffs);

            DNekVec Out_d (ncoeffs,out_d,eWrapper);

            Out_d  = InvMass*Coeffs;
        }

        enum BndToLocMatrixMapType
        {
            eBndToFullMatrixCG,
            eBndToBndMatrixCG,
            eBndToTraceMatrixDG
        };

        void Expansion2D::v_AddRobinMassMatrix(const int edge,
                                        const Array<OneD, const NekDouble > &primCoeffs,
                                               DNekMatSharedPtr &inoutmat)
        {
            ASSERTL1(IsBoundaryInteriorExpansion(),
                     "Not set up for non boundary-interior expansions");
            ASSERTL1(inoutmat->GetRows() == inoutmat->GetColumns(),
                     "Assuming that input matrix was square");
            int i,j;
            int id1,id2;
            ExpansionSharedPtr edgeExp = m_edgeExp[edge].lock();
            int order_e = edgeExp->GetNcoeffs();

            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;

            StdRegions::VarCoeffMap varcoeffs;
            varcoeffs[StdRegions::eVarCoeffMass] = primCoeffs;

            LocalRegions::MatrixKey mkey(StdRegions::eMass,LibUtilities::eSegment, *edgeExp, StdRegions::NullConstFactorMap, varcoeffs);
            DNekScalMat &edgemat = *edgeExp->GetLocMatrix(mkey);

            // Now need to identify a map which takes the local edge
            // mass matrix to the matrix stored in inoutmat;
            // This can currently be deduced from the size of the matrix

            // - if inoutmat.m_rows() == v_NCoeffs() it is a full
            //   matrix system

            // - if inoutmat.m_rows() == v_NumBndCoeffs() it is a
            //  boundary CG system

            // - if inoutmat.m_rows() == v_NumDGBndCoeffs() it is a
            //  trace DG system
            int rows = inoutmat->GetRows();

            if (rows == GetNcoeffs())
            {
                GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);
            }
            else if(rows == NumBndryCoeffs())
            {
                int nbndry = NumBndryCoeffs();
                Array<OneD,unsigned int> bmap(nbndry);

                GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);

                GetBoundaryMap(bmap);

                for(i = 0; i < order_e; ++i)
                {
                    for(j = 0; j < nbndry; ++j)
                    {
                        if(map[i] == bmap[j])
                        {
                            map[i] = j;
                            break;
                        }
                    }
                    ASSERTL1(j != nbndry,"Did not find number in map");
                }
            }
            else if (rows == NumDGBndryCoeffs())
            {
                // possibly this should be a separate method
                int cnt = 0;
                map  = Array<OneD, unsigned int> (order_e);
                sign = Array<OneD, int> (order_e,1);

                for(i = 0; i < edge; ++i)
                {
                    cnt += GetEdgeNcoeffs(i);
                }

                for(i = 0; i < order_e; ++i)
                {
                    map[i] = cnt++;
                }
                // check for mapping reversal
                if(GetEorient(edge) == StdRegions::eBackwards)
                {
                    switch(edgeExp->GetBasis(0)->GetBasisType())
                    {
                    case LibUtilities::eGauss_Lagrange:
                        reverse( map.get() , map.get()+order_e);
                        break;
                    case LibUtilities::eGLL_Lagrange:
                        reverse( map.get() , map.get()+order_e);
                        break;
                    case LibUtilities::eModified_A:
                        {
                            swap(map[0],map[1]);
                            for(i = 3; i < order_e; i+=2)
                            {
                                sign[i] = -1;
                            }
                        }
                        break;
                    default:
                        ASSERTL0(false,"Edge boundary type not valid for this method");
                    }
                }
            }
            else
            {
                ASSERTL0(false,"Could not identify matrix type from dimension");
            }

            for(i = 0; i < order_e; ++i)
            {
                id1 = map[i];
                for(j = 0; j < order_e; ++j)
                {
                    id2 = map[j];
                    (*inoutmat)(id1,id2) +=  edgemat(i,j)*sign[i]*sign[j];
                }
            }
        }

        /**
         * Given an edge and vector of element coefficients:
         * - maps those elemental coefficients corresponding to the edge into
         *   an edge-vector.
         * - resets the element coefficients
         * - multiplies the edge vector by the edge mass matrix
         * - maps the edge coefficients back onto the elemental coefficients
         */
        void Expansion2D::v_AddRobinEdgeContribution(const int edgeid,
                                                     const Array<OneD,
                                                     const NekDouble> &primCoeffs,
                                                     const Array<OneD, NekDouble> &incoeffs,
                                                     Array<OneD, NekDouble> &coeffs)
        {
            ASSERTL1(IsBoundaryInteriorExpansion(),
                     "Not set up for non boundary-interior expansions");
            int i;
            ExpansionSharedPtr edgeExp = m_edgeExp[edgeid].lock();
            int order_e = edgeExp->GetNcoeffs();

            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;

            StdRegions::VarCoeffMap varcoeffs;
            varcoeffs[StdRegions::eVarCoeffMass] = primCoeffs;

            LocalRegions::MatrixKey mkey(StdRegions::eMass,LibUtilities::eSegment,
                                         *edgeExp, StdRegions::NullConstFactorMap,
                                         varcoeffs);
            DNekScalMat &edgemat = *edgeExp->GetLocMatrix(mkey);

            NekVector<NekDouble> vEdgeCoeffs (order_e);

            GetEdgeToElementMap(edgeid,v_GetEorient(edgeid),map,sign);

            for (i = 0; i < order_e; ++i)
            {
                vEdgeCoeffs[i] = incoeffs[map[i]]*sign[i];
            }

            vEdgeCoeffs = edgemat * vEdgeCoeffs;

            for (i = 0; i < order_e; ++i)
            {
                coeffs[map[i]] += vEdgeCoeffs[i]*sign[i];
            }
        }

        DNekMatSharedPtr Expansion2D::v_BuildVertexMatrix(
            const DNekScalMatSharedPtr &r_bnd)
        {
            MatrixStorage storage = eFULL;
            DNekMatSharedPtr m_vertexmatrix;

            int nVerts, vid1, vid2, vMap1, vMap2;
            NekDouble VertexValue;

            nVerts = GetNverts();

            m_vertexmatrix =
                MemoryManager<DNekMat>::AllocateSharedPtr(
                    nVerts, nVerts, 0.0, storage);
            DNekMat &VertexMat = (*m_vertexmatrix);

            for (vid1 = 0; vid1 < nVerts; ++vid1)
            {
                vMap1 = GetVertexMap(vid1);

                for (vid2 = 0; vid2 < nVerts; ++vid2)
                {
                    vMap2 = GetVertexMap(vid2);
                    VertexValue = (*r_bnd)(vMap1, vMap2);
                    VertexMat.SetValue(vid1, vid2, VertexValue);
                }
            }

            return m_vertexmatrix;
        }

        Array<OneD, unsigned int> Expansion2D::v_GetEdgeInverseBoundaryMap(
            int eid)
        {
            int n, j;
            int nEdgeCoeffs;
            int nBndCoeffs = NumBndryCoeffs();

            Array<OneD, unsigned int> bmap(nBndCoeffs);
            GetBoundaryMap(bmap);

            // Map from full system to statically condensed system (i.e reverse
            // GetBoundaryMap)
            map<int, int> invmap;
            for (j = 0; j < nBndCoeffs; ++j)
            {
                invmap[bmap[j]] = j;
            }

            // Number of interior edge coefficients
            nEdgeCoeffs = GetEdgeNcoeffs(eid) - 2;

            const SpatialDomains::Geometry2DSharedPtr &geom = GetGeom2D();

            Array<OneD, unsigned int> edgemaparray(nEdgeCoeffs);
            Array<OneD, unsigned int> maparray    (nEdgeCoeffs);
            Array<OneD, int>          signarray   (nEdgeCoeffs, 1);
            StdRegions::Orientation eOrient      = geom->GetEorient(eid);

            // maparray is the location of the edge within the matrix
            GetEdgeInteriorMap(eid, eOrient, maparray, signarray);

            for (n = 0; n < nEdgeCoeffs; ++n)
            {
                edgemaparray[n] = invmap[maparray[n]];
            }

            return edgemaparray;
        }

        void Expansion2D::v_SetUpPhysNormals(const int edge)
        {
            ComputeEdgeNormal(edge);
        }

        const StdRegions::NormalVector &Expansion2D::v_GetEdgeNormal(const int edge) const
        {
            std::map<int, StdRegions::NormalVector>::const_iterator x;
            x = m_edgeNormals.find(edge);
            ASSERTL0 (x != m_edgeNormals.end(),
                        "Edge normal not computed.");
            return x->second;
        }

        const StdRegions::NormalVector &Expansion2D::v_GetSurfaceNormal(
                const int id) const
        {
            return v_GetEdgeNormal(id);
        }

        void Expansion2D::ReOrientEdgePhysMap(
            const int                                   nvert,
            const StdRegions::Orientation               orient,
            const int                                   nq0,
            Array<OneD, int>                            &idmap)
        {
            boost::ignore_unused(nvert);

            if (idmap.size() != nq0)
            {
                idmap = Array<OneD, int>(nq0);
            }
            switch (orient)
            {
                case StdRegions::eForwards:
                    // Fwd
                    for (int i = 0; i < nq0; ++i)
                    {
                        idmap[i] = i;
                    }
                    break;
                case StdRegions::eBackwards:
                {
                    // Bwd
                    for (int i = 0; i < nq0; ++i)
                    {
                        idmap[i] = nq0-1-i;
                    }
                }
                    break;
                default:
                    ASSERTL0(false, "Unknown orientation");
                    break;
            }
        }

        Array<OneD, NekDouble> Expansion2D::v_GetMF(
            const int dir,
            const int shapedim,
            const StdRegions::VarCoeffMap   &varcoeffs)
        {
            return Expansion::v_GetMF(dir,shapedim,varcoeffs);
        }

        Array<OneD, NekDouble> Expansion2D::v_GetMFDiv(
            const int dir,
            const StdRegions::VarCoeffMap   &varcoeffs)
        {
            return Expansion::v_GetMFDiv(dir,varcoeffs);
        }

        Array<OneD, NekDouble> Expansion2D::v_GetMFMag(
            const int dir,
            const StdRegions::VarCoeffMap   &varcoeffs)
        {
            return Expansion::v_GetMFMag(dir,varcoeffs);
        }

        // Compute edgenormal \cdot vector
        Array<OneD, NekDouble> Expansion2D::v_GetnEdgecdotMF(
            const int dir,
            const int edge,
            ExpansionSharedPtr &EdgeExp_e,
            const Array<OneD, const Array<OneD, NekDouble> > &normals,
            const StdRegions::VarCoeffMap   &varcoeffs)
        {
            int nquad_e = EdgeExp_e->GetNumPoints(0);
            int coordim = GetCoordim();
            int nquad0  = m_base[0]->GetNumPoints();
            int nquad1  = m_base[1]->GetNumPoints();
            int nqtot   = nquad0*nquad1;

            StdRegions::VarCoeffType MMFCoeffs[15] =
            {
                StdRegions::eVarCoeffMF1x,
                StdRegions::eVarCoeffMF1y,
                StdRegions::eVarCoeffMF1z,
                StdRegions::eVarCoeffMF1Div,
                StdRegions::eVarCoeffMF1Mag,
                StdRegions::eVarCoeffMF2x,
                StdRegions::eVarCoeffMF2y,
                StdRegions::eVarCoeffMF2z,
                StdRegions::eVarCoeffMF2Div,
                StdRegions::eVarCoeffMF2Mag,
                StdRegions::eVarCoeffMF3x,
                StdRegions::eVarCoeffMF3y,
                StdRegions::eVarCoeffMF3z,
                StdRegions::eVarCoeffMF3Div,
                StdRegions::eVarCoeffMF3Mag
            };

            StdRegions::VarCoeffMap::const_iterator MFdir;

            Array<OneD, NekDouble> ncdotMF(nqtot,0.0);
            Array<OneD, NekDouble> tmp(nqtot);
            Array<OneD, NekDouble> tmp_e(nquad_e);
            for (int k=0; k<coordim; k++)
            {
                MFdir = varcoeffs.find(MMFCoeffs[dir*5+k]);
                tmp = MFdir->second;

                GetPhysEdgeVarCoeffsFromElement(edge, EdgeExp_e, tmp, tmp_e);

                Vmath::Vvtvp(nquad_e,
                             &tmp_e[0],      1,
                             &normals[k][0], 1,
                             &ncdotMF[0],    1,
                             &ncdotMF[0],    1);
            }
            return ncdotMF;
        }

        NekDouble Expansion2D::v_VectorFlux(
            const Array<OneD, Array<OneD, NekDouble> > &vec)
        {
            const Array<OneD, const Array<OneD, NekDouble> >
                &normals = GetLeftAdjacentElementExp()->
                GetFaceNormal(GetLeftAdjacentElementFace());

            int nq = GetTotPoints();
            Array<OneD, NekDouble > Fn(nq);
            Vmath::Vmul (nq, &vec[0][0], 1, &normals[0][0], 1, &Fn[0], 1);
            Vmath::Vvtvp(nq, &vec[1][0], 1, &normals[1][0], 1, &Fn[0], 1, &Fn[0], 1);
            Vmath::Vvtvp(nq, &vec[2][0], 1, &normals[2][0], 1, &Fn[0], 1, &Fn[0], 1);

            return StdExpansion::Integral(Fn);
        }
    }
}
