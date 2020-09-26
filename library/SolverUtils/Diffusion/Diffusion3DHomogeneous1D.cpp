///////////////////////////////////////////////////////////////////////////////
//
// File: Diffusion3DHomogeneous1D.cpp
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
// Description: LDG diffusion 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Diffusion/Diffusion3DHomogeneous1D.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string Diffusion3DHomogeneous1D::type[] = {
            GetDiffusionFactory().RegisterCreatorFunction(
                "LDG3DHomogeneous1D",       Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRDG3DHomogeneous1D",     Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRSD3DHomogeneous1D",     Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRHU3DHomogeneous1D",     Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRcmin3DHomogeneous1D",   Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRcinf3DHomogeneous1D",   Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LDGNS3DHomogeneous1D",     Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRDGNS3DHomogeneous1D",   Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRSDNS3DHomogeneous1D",   Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRHUNS3DHomogeneous1D",   Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRcminNS3DHomogeneous1D", Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRcinfNS3DHomogeneous1D", Diffusion3DHomogeneous1D::create)
        };

        /**
         * @brief Diffusion3DHomogeneous1D  uses the 2D WeakDG approach
         * to compute the diffusion term looping on the planes in the z
         * direction and adding the flux in z direction at the end.
         */
        Diffusion3DHomogeneous1D::Diffusion3DHomogeneous1D(std::string diffType)
        {
            // Strip trailing string "3DHomogeneous1D" to determine 2D diffusion
            // type, and create a diffusion object for the plane.
            m_diffType = diffType.substr(0, diffType.length()-15);
            m_planeDiff = GetDiffusionFactory().CreateInstance(m_diffType, m_diffType);
        }

        /**
         * @brief Initiliase Diffusion3DHomogeneous1D objects and store
         * them before starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void Diffusion3DHomogeneous1D::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            int nConvectiveFields = pFields.size();

            Array<OneD, MultiRegions::ExpListSharedPtr> pFields_plane0(
                nConvectiveFields);

            // Initialise the plane advection object.
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                pFields_plane0[i] = pFields[i]->GetPlane(0);
            }
            m_planeDiff->InitObject(pSession, pFields_plane0);

            m_numPoints      = pFields[0]->GetTotPoints();
            m_planes         = pFields[0]->GetZIDs();
            m_numPlanes      = m_planes.size();
            m_numPointsPlane = m_numPoints/m_numPlanes;
            m_homoLen        = pFields[0]->GetHomoLen();
            m_trans          = pFields[0]->GetTransposition();
            m_planeCounter = 0;

            if (m_diffType == "LDG")
            {
                // Set viscous flux for LDG
                m_planeDiff->SetFluxVector(m_fluxVector);
            }
            else if (m_diffType == "LDGNS")
            {
                // Set viscous flux for LDGNS
                m_planeDiff->SetFluxVectorNS(m_fluxVectorNS);
                // Set penalty flux
                m_planeDiff->SetFluxPenaltyNS(m_fluxPenaltyNS);
            }
            else if (m_diffType == "LFRDGNS" ||
                     m_diffType == "LFRHUNS" ||
                     m_diffType == "LFRSDNS" )
            {
                // Set viscous flux for FR cases
                m_planeDiff->SetFluxVectorNS(m_fluxVectorNS);
            }

            m_fieldsPlane   = Array<OneD, MultiRegions::ExpListSharedPtr>
                                                            (nConvectiveFields);

            if (m_fluxVectorNS)
            {
                m_inarrayPlane  = Array<OneD, Array<OneD, NekDouble> >
                                                        (nConvectiveFields - 1);
            }
            else
            {
                m_inarrayPlane  = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
            }
            m_outarrayPlane = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
            m_planePos      = Array<OneD, unsigned int>     (m_numPlanes);

            for (int i = 0; i < m_numPlanes; ++i)
            {
                m_planePos[i] = i * m_numPointsPlane;
            }

            if (m_fluxVectorNS)
            {
                m_homoDerivStore = Array<OneD, Array<OneD, NekDouble> >(
                    nConvectiveFields);
                m_homoDerivPlane = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(
                    m_numPlanes);

                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    m_homoDerivStore[i] = Array<OneD, NekDouble>(m_numPoints);
                }

                for (int i = 0; i < m_numPlanes; ++i)
                {
                    m_homoDerivPlane[i] = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);

                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        m_homoDerivPlane[i][j] = Array<OneD, NekDouble>(
                            m_numPointsPlane,
                            m_homoDerivStore[j] + m_planePos[i]);
                    }
                }
            }
        }

        /**
         * @brief Calculate WeakDG Diffusion for the linear problems
         * using an LDG interface flux and the the flux in the third direction.
         */
        void Diffusion3DHomogeneous1D::v_Diffuse(
            const std::size_t                                 nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            boost::ignore_unused(pFwd, pBwd);

            Array<OneD, NekDouble> tmp(m_numPoints), tmp2;
            Array<OneD, Array<OneD, NekDouble> > viscHComp;
            const int nPointsTot = fields[0]->GetNpoints();
            NekDouble beta;


            if (m_fluxVectorNS)
            {
                viscHComp = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (int i = 0; i < nConvectiveFields - 1; ++i)
                {
                    fields[0]->PhysDeriv(2, inarray[i], m_homoDerivStore[i]);
                    viscHComp[i] = Array<OneD, NekDouble>(m_numPoints);
                }
            }


            for (int i = 0; i < m_numPlanes; ++i)
            {
                // Set up memory references for fields, inarray and outarray for
                // this plane.
                for (int j = 0; j < inarray.size(); ++j)
                {
                    m_inarrayPlane [j] = Array<OneD, NekDouble>(
                        m_numPointsPlane, tmp2 = inarray [j] + m_planePos[i]);
                }

                for (int j = 0; j < nConvectiveFields; ++j)
                {
                    m_fieldsPlane  [j] = fields[j]->GetPlane(i);
                    m_outarrayPlane[j] = Array<OneD, NekDouble>(
                        m_numPointsPlane, tmp2 = outarray[j] + m_planePos[i]);
                }


                if (m_fluxVectorNS)
                {
                    m_planeDiff->SetHomoDerivs(m_homoDerivPlane[i]);
                }


                if (m_diffType == "LDGNS")
                {
                    // Store plane Fwd/Bwd traces
                    std::size_t nTracePts = m_fieldsPlane[0]->GetTrace()
                        ->GetTotPoints();
                    std::size_t nScalar = m_inarrayPlane.size();
                    Array<OneD, Array<OneD, NekDouble> > Fwd(nScalar);
                    Array<OneD, Array<OneD, NekDouble> > Bwd(nScalar);
                    {
                        for(std::size_t k = 0; k < nScalar; ++k)
                        {
                            Fwd[k] = Array<OneD, NekDouble>(nTracePts, 0.0);
                            Bwd[k] = Array<OneD, NekDouble>(nTracePts, 0.0);
                            m_fieldsPlane[k]->GetFwdBwdTracePhys(
                                m_inarrayPlane[k], Fwd[k], Bwd[k]);
                        }
                    }

                    m_planeDiff->Diffuse(nConvectiveFields,
                                         m_fieldsPlane,
                                         m_inarrayPlane,
                                         m_outarrayPlane,
                                         Fwd,
                                         Bwd);
                }
                else
                {
                    m_planeDiff->Diffuse(nConvectiveFields,
                                         m_fieldsPlane,
                                         m_inarrayPlane,
                                         m_outarrayPlane);
                }

                if (m_fluxVectorNS)
                {
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                        &viscTensor = m_planeDiff->GetFluxTensor();

                    // Extract H (viscTensor[2])
                    for (int j = 0; j < nConvectiveFields - 1; ++j)
                    {
                        Vmath::Vcopy(m_numPointsPlane,
                                     viscTensor[2][j+1], 1,
                                     tmp2 = viscHComp[j] + m_planePos[i], 1);
                    }
                }
            }



            if (m_fluxVectorNS)
            {
                for (int j = 0; j < nConvectiveFields - 1; ++j)
                {
                    fields[j+1]->PhysDeriv(2, viscHComp[j], tmp);
                    Vmath::Vadd(nPointsTot, outarray[j+1], 1, tmp, 1,
                        outarray[j+1], 1);
                }
            }
            else
            {
                for (int j = 0; j < nConvectiveFields; ++j)
                {
                    fields[j]->HomogeneousFwdTrans(inarray[j], tmp);

                    for (int i = 0; i < m_numPlanes; ++i)
                    {
                        beta  = 2*M_PI*m_trans->GetK(i)/m_homoLen;
                        beta *= beta;

                        Vmath::Smul(m_numPointsPlane,
                                    beta,
                                    &tmp[0] + i*m_numPointsPlane, 1,
                                    &tmp[0] + i*m_numPointsPlane, 1);
                    }

                    fields[0]->HomogeneousBwdTrans(tmp, tmp);

                    Vmath::Vsub(nPointsTot, outarray[j], 1, tmp, 1,
                                outarray[j], 1);
                }
            }
        }
    }// close namespace SolverUtils
}// close namespace nektar++
