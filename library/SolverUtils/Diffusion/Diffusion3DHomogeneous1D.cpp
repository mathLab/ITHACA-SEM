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

using namespace std;

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
            m_trans          = pFields[0]->GetTransposition();
            m_planeCounter = 0;
            m_planeDiff->SetFluxVectorNS(m_fluxVectorNS);

            if (m_riemann)
            {
                // Set Riemann solver and flux vector callback for this plane.
                m_planeDiff->SetRiemannSolver(m_riemann);

                // Override Riemann solver scalar and vector callbacks.
                map<string, RSScalarFuncType>::iterator it1;
                map<string, RSScalarFuncType> scalars = m_riemann->GetScalars();

                for (it1 = scalars.begin(); it1 != scalars.end(); ++it1)
                {
                    boost::shared_ptr<HomoRSScalar> tmp =
                        MemoryManager<HomoRSScalar>
                            ::AllocateSharedPtr(it1->second, m_numPlanes);
                    m_riemann->SetScalar(it1->first, &HomoRSScalar::Exec, tmp);
                }
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
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            
            Array<OneD, NekDouble> tmp(m_numPoints), tmp2;
            Array<OneD, Array<OneD, NekDouble> > viscHComp;
            const int nPointsTot = fields[0]->GetNpoints();
            int i, j;
            NekDouble beta;
            

            if (m_fluxVectorNS)
            {
                viscHComp = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (i = 0; i < nConvectiveFields - 1; ++i)
                {
                    fields[0]->PhysDeriv(2, inarray[i], m_homoDerivStore[i]);
                    viscHComp[i] = Array<OneD, NekDouble>(m_numPoints);
                }
            }
            

            for (i = 0; i < m_numPlanes; ++i)
            {
                // Set up memory references for fields, inarray and outarray for
                // this plane.
                for (int j = 0; j < inarray.num_elements(); ++j)
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


                
                m_planeDiff->Diffuse(nConvectiveFields,
                                     m_fieldsPlane,
                                     m_inarrayPlane,
                                     m_outarrayPlane);
                
                if (m_fluxVectorNS)
                {
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscTensor = m_planeDiff->GetFluxTensor();

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
                for (j = 0; j < nConvectiveFields - 1; ++j)
                {
                    fields[j+1]->PhysDeriv(2, viscHComp[j], tmp);
                    Vmath::Vadd(nPointsTot, outarray[j+1], 1, tmp, 1, outarray[j+1], 1);
                }
            }
            else
            {
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

                    Vmath::Vsub(nPointsTot, outarray[j], 1, tmp, 1,
                                outarray[j], 1);
                }
            }
        }
    }// close namespace SolverUtils
}// close namespace nektar++
