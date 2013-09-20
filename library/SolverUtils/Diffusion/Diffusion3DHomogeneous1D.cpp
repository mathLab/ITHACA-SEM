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
// Description: LDG diffusion 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/Diffusion3DHomogeneous1D.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string Diffusion3DHomogeneous1D::type[] = {
            GetDiffusionFactory().RegisterCreatorFunction(
                "LDG3DHomogeneous1D",   Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFR3DHomogeneous1D",   Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LDGNS3DHomogeneous1D", Diffusion3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRNS3DHomogeneous1D", Diffusion3DHomogeneous1D::create)
        };

        /**
         * @brief Diffusion3DHomogeneous1D  uses the 2D WeakDG approach 
         * to compute the diffusion term looping on the planes in the z
         * direction and adding the flux in z direction at the end.
         */
        Diffusion3DHomogeneous1D::Diffusion3DHomogeneous1D(std::string diffType)
        {
            // Strip trailing string "3DHomogeneous1D" to determine 2D advection
            // type, and create an advection object for the plane.
            string name = diffType.substr(0, diffType.length()-15);
            m_planeDiff = GetDiffusionFactory().CreateInstance(name, name);
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
            int nConvectiveFields = pFields.num_elements();

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
            m_numPlanes      = m_planes.num_elements();
            m_numPointsPlane = m_numPoints/m_numPlanes;
            m_homoLen        = pFields[0]->GetHomoLen();

            // Set Riemann solver and flux vector callback for this plane.
            m_planeDiff->SetRiemannSolver(m_riemann);
            //m_planeDiff->SetFluxVector   (
            //    &Diffusion3DHomogeneous1D::ModifiedFluxVector, this);
            m_planeCounter = 0;

            // Override Riemann solver scalar and vector callbacks.
            map<string, RSScalarFuncType>::iterator it1;
            map<string, RSScalarFuncType> scalars = m_riemann->GetScalars();

            for (it1 = scalars.begin(); it1 != scalars.end(); ++it1)
            {
                boost::shared_ptr<HomoRSScalar> tmp = MemoryManager<HomoRSScalar>
                   ::AllocateSharedPtr(it1->second, m_numPlanes);
                m_riemann->SetScalar(it1->first, &HomoRSScalar::Exec, tmp);
            }

            m_fluxVecStore = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(
                nConvectiveFields);

            // Set up storage for flux vector.
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                m_fluxVecStore[i] = Array<OneD, Array<OneD, NekDouble> >(3);
                for (int j = 0; j < 3; ++j)
                {
                    m_fluxVecStore[i][j] = Array<OneD, NekDouble>(m_numPoints);
                }
            }

            m_fluxVecPlane = Array<OneD, Array<OneD,
                          Array<OneD, Array<OneD, NekDouble> > > >(m_numPlanes);
            m_fieldsPlane   = Array<OneD, MultiRegions::ExpListSharedPtr>
                                                            (nConvectiveFields);
            m_inarrayPlane  = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
            m_outarrayPlane = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
            m_planePos      = Array<OneD, unsigned int>     (m_numPlanes);

            // Set up memory reference which links fluxVecPlane to fluxVecStore.
            for (int i = 0; i < m_numPlanes; ++i)
            {
                m_planePos[i] = i * m_numPointsPlane;
                m_fluxVecPlane[i] =
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(
                        nConvectiveFields);

                for (int j = 0; j < nConvectiveFields; ++j)
                {
                    m_fluxVecPlane[i][j] =
                        Array<OneD, Array<OneD, NekDouble> >(3);
                    for (int k = 0; k < 3; ++k)
                    {
                        m_fluxVecPlane[i][j][k] = Array<OneD, NekDouble>(
                            m_numPointsPlane,
                            m_fluxVecStore[j][k] + m_planePos[i]);
                    }
                }
            }
        }
        
        /**
         * @brief Calculate WeakDG Diffusion for the linear problems
         * using an LDG interface flux and the the flux in the third direction.
         */
        void Diffusion3DHomogeneous1D::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            Array<OneD, NekDouble> tmp(m_numPoints), tmp2;
            const int nPointsTot = fields[0]->GetNpoints();
            int i, j;
            NekDouble beta;

            for (i = 0; i < m_numPlanes; ++i)
            {
                // Set up memory references for fields, inarray and outarray for
                // this plane.
                for (int j = 0; j < nConvectiveFields; ++j)
                {
                    m_fieldsPlane  [j] = fields[j]->GetPlane(i);
                    m_inarrayPlane [j] = Array<OneD, NekDouble>(
                        m_numPointsPlane, tmp2 = inarray [j] + m_planePos[i]);
                    m_outarrayPlane[j] = Array<OneD, NekDouble>(
                        m_numPointsPlane, tmp2 = outarray[j] + m_planePos[i]);
                }

                m_planeDiff->Diffuse(nConvectiveFields,
                                     m_fieldsPlane,
                                     m_inarrayPlane,
                                     m_outarrayPlane);
            }

            for (j = 0; j < nConvectiveFields; ++j)
            {
                fields[j]->HomogeneousFwdTrans(inarray[j], tmp);

                for (i = 0; i < m_numPlanes; ++i)
                {
                    beta  = 2*M_PI*m_trans->GetK(i)/m_homoLen;
                    beta *= beta;

                    Vmath::Smul(m_numPointsPlane,
                                beta,
                                &tmp[0] + i*m_numPointsPlane, 1,
                                &tmp[0] + i*m_numPointsPlane, 1);
                }

                fields[0]->HomogeneousBwdTrans(tmp, tmp);

                Vmath::Vsub(nPointsTot, outarray[j], 1, tmp, 1, outarray[j], 1);
            }
        }

        void Diffusion3DHomogeneous1D::ModifiedFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
        {
            // Return section of flux vector for this plane.
            outarray = m_fluxVecPlane[m_planeCounter];

            // Increment the plane counter.
            m_planeCounter = (m_planeCounter + 1) % m_numPlanes;
        }
    }// close namespace SolverUtils
}// close namespace nektar++
