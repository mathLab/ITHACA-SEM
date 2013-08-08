///////////////////////////////////////////////////////////////////////////////
//
// File: Advection3DHomogeneous1D.cpp
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
// Description: FR advection 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/Advection3DHomogeneous1D.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <StdRegions/StdSegExp.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <iomanip>


namespace Nektar
{
    namespace SolverUtils
    {
        std::string Advection3DHomogeneous1D::type[] = {
            GetAdvectionFactory().RegisterCreatorFunction(
                "WeakDG3DHomogeneous1D",   Advection3DHomogeneous1D::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                "FRDG3DHomogeneous1D",   Advection3DHomogeneous1D::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                "FRDG3DHomogeneous1D",   Advection3DHomogeneous1D::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                "FRSD3DHomogeneous1D",   Advection3DHomogeneous1D::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                "FRHU3DHomogeneous1D",   Advection3DHomogeneous1D::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                "FRcmin3DHomogeneous1D", Advection3DHomogeneous1D::create),
            GetAdvectionFactory().RegisterCreatorFunction(
                "FRcinf3DHomogeneous1D", Advection3DHomogeneous1D::create)
        };

        /**
         * @brief AdvectionFR uses the Flux Reconstruction (FR) approach to
         * compute the advection term. The implementation is only for segments,
         * quadrilaterals and hexahedra at the moment.
         *
         * \todo Extension to triangles, tetrahedra and other shapes.
         * (Long term objective)
         */
        Advection3DHomogeneous1D::Advection3DHomogeneous1D(std::string advType)
          : m_advType(advType)
        {
            string advName = advType.substr(
                0, advType.length()-15);//advType.find_first_of("3DHomogeneous1D"));
            m_planeAdv = GetAdvectionFactory().CreateInstance(advName, advName);
        }

        /**
         * @brief Initiliase Advection3DHomogeneous1D objects and store them
         * before starting the time-stepping.
         *
         * This routine calls the virtual functions #v_SetupMetrics,
         * #v_SetupCFunctions and #v_SetupInterpolationMatrices to
         * initialise the objects needed by AdvectionFR.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void Advection3DHomogeneous1D::v_InitObject(
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
            m_planeAdv->InitObject(pSession, pFields_plane0);

            nPointsTot       = pFields[0]->GetTotPoints();
            nCoeffs          = pFields[0]->GetNcoeffs();
            m_planes         = pFields[0]->GetZIDs();
            num_planes       = m_planes.num_elements();
            nPointsTot_plane = nPointsTot/num_planes;
            nCoeffs_plane    = nCoeffs/num_planes;

            // Set Riemann solver and flux vector callback for this plane.
            m_planeAdv->SetRiemannSolver(m_riemann);
            m_planeAdv->SetFluxVector   (
                &Advection3DHomogeneous1D::ModifiedFluxVector, this);
            m_planeCounter = 0;

            m_fluxVecStore = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(
                nConvectiveFields);

            // Set up storage for flux vector.
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                m_fluxVecStore[i] = Array<OneD, Array<OneD, NekDouble> >(3);
                for (int j = 0; j < 3; ++j)
                {
                    m_fluxVecStore[i][j] = Array<OneD, NekDouble>(nPointsTot);
                }
            }

            m_fluxVecPlane = Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >(
                num_planes);
            m_fieldsPlane   = Array<OneD, MultiRegions::ExpListSharedPtr>(
                nConvectiveFields);
            m_inarrayPlane  = Array<OneD, Array<OneD, NekDouble> > (
                nConvectiveFields);
            m_outarrayPlane = Array<OneD, Array<OneD, NekDouble> > (
                nConvectiveFields);
            m_advVelPlane   = Array<OneD, Array<OneD, NekDouble> > (3);

            // Set up memory reference which links fluxVecPlane to fluxVecStore.
            for (int i = 0; i < num_planes; ++i)
            {
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
                            nPointsTot_plane, m_fluxVecStore[j][k] + i*nPointsTot_plane);
                    }
                }
            }
        }

        /**
         * @brief Compute the advection term at each time-step using the Flux
         * Reconstruction approach (FR) looping on the planes.
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Solution at the previous time-step.
         * @param outarray            Advection term to be passed at the
         *                            time integration class.
         */
        void Advection3DHomogeneous1D::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            Array<OneD, NekDouble> tmp(nPointsTot);
            int nVel = advVel.num_elements();

            // Call flux vector function for entire domain.
            m_fluxVector(inarray, m_fluxVecStore);

            for (int i = 0; i < num_planes; ++i)
            {
                // Set up memory references.
                for (int j = 0; j < nConvectiveFields; ++j)
                {
                    Array<OneD, NekDouble> tmp2;
                    m_fieldsPlane[j]   = fields[j]->GetPlane(i);
                    m_inarrayPlane[j]  = Array<OneD, NekDouble>(nPointsTot_plane, tmp2 = inarray[j] + i*nPointsTot_plane);
                    m_outarrayPlane[j] = Array<OneD, NekDouble>(nPointsTot_plane, tmp2 = outarray[j] + i*nPointsTot_plane);
                }

                /*
                Array<OneD, NekDouble> tmp2;
                for (int j = 0; j < nVel; ++j)
                {
                    cout << "advVel num elements = " << advVel[j].num_elements() << endl;
                    cout << "offset = " << (i*nPointsTot_plane) << endl;
                    if (advVel[j].num_elements() != 0)
                    {
                        m_advVelPlane[j] = Array<OneD, NekDouble>(
                            nPointsTot_plane, tmp2 = advVel[j] + i*nPointsTot_plane);
                    }
                }
                */

                m_planeAdv->Advect(nConvectiveFields, m_fieldsPlane,
                                   m_advVelPlane, m_inarrayPlane,
                                   m_outarrayPlane);
            }

            // Calculate Fourier derivative
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                fields[0]->PhysDeriv(2, m_fluxVecStore[i][2], tmp);

                Vmath::Vadd(nPointsTot, outarray[i], 1, tmp, 1,
                                        outarray[i], 1);
            }
        }

        void Advection3DHomogeneous1D::ModifiedFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
        {
            // Increment the plane counter.
            m_planeCounter = (m_planeCounter + 1) % num_planes;
            outarray = m_fluxVecPlane[m_planeCounter];
        }
    }
}
