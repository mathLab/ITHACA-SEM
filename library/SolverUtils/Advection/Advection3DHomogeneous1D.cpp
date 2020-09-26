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

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Advection/Advection3DHomogeneous1D.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <iostream>
#include <iomanip>

using namespace std;

namespace Nektar
{
    namespace SolverUtils
    {
        std::string Advection3DHomogeneous1D::type[] = {
            GetAdvectionFactory().RegisterCreatorFunction(
                "WeakDG3DHomogeneous1D", Advection3DHomogeneous1D::create),
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
            // Strip trailing string "3DHomogeneous1D" to determine 2D advection
            // type, and create an advection object for the plane.
            string advName = advType.substr(0, advType.length()-15);
            m_planeAdv = GetAdvectionFactory().CreateInstance(advName, advName);
        }

        /**
         * @brief Initiliase Advection3DHomogeneous1D objects and store them
         * before starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void Advection3DHomogeneous1D::v_InitObject(
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
            m_planeAdv->InitObject(pSession, pFields_plane0);

            m_numPoints      = pFields[0]->GetTotPoints();
            m_planes         = pFields[0]->GetZIDs();
            m_numPlanes      = m_planes.size();
            m_numPointsPlane = m_numPoints/m_numPlanes;

            // Set Riemann solver and flux vector callback for this plane.
            m_planeAdv->SetRiemannSolver(m_riemann);
            m_planeAdv->SetFluxVector   (
                &Advection3DHomogeneous1D::ModifiedFluxVector, this);
            m_planeCounter = 0;

            // Override Riemann solver scalar and vector callbacks.
            map<string, RSScalarFuncType> scalars = m_riemann->GetScalars();
            map<string, RSVecFuncType> vectors = m_riemann->GetVectors();

            for (auto &it1 : scalars)
            {
                std::shared_ptr<HomoRSScalar> tmp = MemoryManager<HomoRSScalar>
                    ::AllocateSharedPtr(it1.second, m_numPlanes);
                m_riemann->SetScalar(it1.first, &HomoRSScalar::Exec, tmp);
            }

            for (auto &it2 : vectors)
            {
                std::shared_ptr<HomoRSVector> tmp = MemoryManager<HomoRSVector>
                    ::AllocateSharedPtr(it2.second, m_numPlanes, it2.first);
                m_riemann->SetVector(it2.first, &HomoRSVector::Exec, tmp);
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
            m_advVelPlane   = Array<OneD, Array<OneD, NekDouble> > (3);

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
         * @brief Compute the advection operator for a given input @a inarray
         * and put the result in @a outarray.
         *
         * @param nConvectiveFields   Number of fields to advect.
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Input which will be advected.
         * @param outarray            Computed advection.
         */
        void Advection3DHomogeneous1D::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            boost::ignore_unused(pFwd, pBwd);

            Array<OneD, NekDouble> tmp(m_numPoints), tmp2;
            int nVel = advVel.size();

            // Call solver's flux vector function to compute the flux vector on
            // the entire domain.
            m_fluxVector(inarray, m_fluxVecStore);

            // Loop over each plane.
            for (int i = 0; i < m_numPlanes; ++i)
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

                for (int j = 0; j < nVel; ++j)
                {
                    if (advVel[j].size() != 0)
                    {
                        m_advVelPlane[j] = Array<OneD, NekDouble>(
                            m_numPointsPlane, tmp2 = advVel[j] + m_planePos[i]);
                    }
                }

                // Compute advection term for this plane.
                m_planeAdv->Advect(nConvectiveFields, m_fieldsPlane,
                                   m_advVelPlane, m_inarrayPlane,
                                   m_outarrayPlane, time);
            }

            // Calculate Fourier derivative and add to final result.
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                fields[0]->PhysDeriv(2, m_fluxVecStore[i][2], tmp);

                Vmath::Vadd(m_numPoints, outarray[i], 1, tmp, 1,
                                         outarray[i], 1);
            }
        }

        void Advection3DHomogeneous1D::ModifiedFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &outarray)
        {
            boost::ignore_unused(inarray);

            // Return section of flux vector for this plane.
            outarray = m_fluxVecPlane[m_planeCounter];

            // Increment the plane counter.
            m_planeCounter = (m_planeCounter + 1) % m_numPlanes;
        }
    }
}
