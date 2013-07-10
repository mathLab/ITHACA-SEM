///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG3DHomogeneous1D.cpp
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
// Description: Weak DG advection 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/AdvectionWeakDG3DHomogeneous1D.h>

#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string AdvectionWeakDG3DHomogeneous1D::type = 
            GetAdvectionFactory().RegisterCreatorFunction(
                "WeakDG3DHomogeneous1D", 
                AdvectionWeakDG3DHomogeneous1D::create);
        
        AdvectionWeakDG3DHomogeneous1D::AdvectionWeakDG3DHomogeneous1D()
        {
            string advName = "WeakDG";
            m_planeAdv = GetAdvectionFactory().CreateInstance(advName, advName);
        }
        
        /**
         * @brief Initiliase AdvectionWeakDG3DHomogeneous1D objects and store 
         * them before starting the time-stepping.
         *
         * This routine calls the virtual functions #v_SetupMetrics,
         * #v_SetupCFunctions and #v_SetupInterpolationMatrices to
         * initialise the objects needed by AdvectionFR.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void AdvectionWeakDG3DHomogeneous1D::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            int nConvectiveFields = pFields.num_elements();
            
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields_plane0
                                                            (nConvectiveFields);
            
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                pFields_plane0[i] = pFields[i]->GetPlane(0);
            }
            
            m_planeAdv->InitObject(pSession, pFields_plane0);
            
            nPointsTot      = pFields[0]->GetTotPoints();
            nCoeffs         = pFields[0]->GetNcoeffs();
            
            planes = pFields[0]->GetZIDs();
            num_planes = planes.num_elements();
            nPointsTot_plane = nPointsTot/num_planes;
            
            m_planeAdv->SetRiemannSolver(m_riemann);
            m_planeAdv->SetFluxVectorVec(m_fluxVector);
            
            fluxvector =
                Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
            
            for (j = 0; j < nConvectiveFields; j ++)
            {
                fluxvector[j] = Array<OneD, NekDouble>(nPointsTot, 0.0);
            }
            
            fluxvector_homo =
                Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                            (num_planes);
            outarray_homo = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
            fields_plane =
                Array <OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                                                            (num_planes);
            inarray_plane = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                            (num_planes);
            outarray_plane = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                            (num_planes);
            advVel_plane = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                            (num_planes);
            
            for (i = 0; i < num_planes; ++i)
            {
                fields_plane[i] = Array<OneD, MultiRegions::ExpListSharedPtr>
                                                            (nConvectiveFields);
                inarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
                outarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
                fluxvector_homo[i] =
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                            (nConvectiveFields);
                advVel_plane[i] = Array<OneD, Array<OneD, NekDouble> >(3);
                
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    inarray_plane[i][j] = Array<OneD, NekDouble>
                                                        (nPointsTot_plane, 0.0);
                    outarray_plane[i][j] = Array<OneD, NekDouble>
                                                        (nPointsTot_plane, 0.0);
                    outarray_homo[j] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                    
                    fluxvector_homo[i][j] = Array<OneD, Array<OneD, NekDouble> >
                                                                            (3);
                    
                    for (int k = 0; k < 3 ; ++k)
                    {
                        fluxvector_homo[i][j][k] =
                            Array<OneD, NekDouble>(nPointsTot_plane, 0.0);
                        
                        advVel_plane[i][k] = Array<OneD, NekDouble>
                                                        (nPointsTot_plane, 0.0);
                    }
                }
            }
        }
        
        /**
        * @brief Compute the advection term at each time-step using the 
        * Discontinuous Galerkin approach (DG) looping on the planes.
        *
        * @param nConvectiveFields   Number of fields.
        * @param fields              Pointer to fields.
        * @param advVel              Advection velocities.
        * @param inarray             Solution at the previous time-step.
        * @param outarray            Advection term to be passed at the
        *                            time integration class.
        *
        */
        void AdvectionWeakDG3DHomogeneous1D::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int nVel = advVel.num_elements();
            
            for (i = 0; i < num_planes; ++i)
            {
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    fields_plane[i][j]= fields[j]->GetPlane(i);
                    
                    Vmath::Vcopy(nPointsTot_plane,
                                 &inarray[j][i * nPointsTot_plane], 1,
                                 &inarray_plane[i][j][0], 1);
                }
                
                for (j = 0; j < nVel; j ++)
                {
                    if (advVel[j].num_elements() != 0 )
                    {   Vmath::Vcopy(nPointsTot_plane,
                                     &advVel[j][i * nPointsTot_plane], 1,
                                     &advVel_plane[i][j][0], 1);
                    }
                }
             
                m_planeAdv->Advect(nConvectiveFields, fields_plane[i],
                                   advVel_plane[i], inarray_plane[i],
                                   outarray_plane[i]);
            
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &outarray_plane[i][j][0], 1,
                                 &outarray[j][i * nPointsTot_plane], 1);
                }
                
                m_planeAdv->FluxVec(fluxvector_homo[i]);
                
                for ( j = 0; j < nConvectiveFields; ++j)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &fluxvector_homo[i][j][2][0], 1,
                                 &fluxvector[j][i * nPointsTot_plane], 1);
                }
            }
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                fields[0]->PhysDeriv(2, fluxvector[i], outarray_homo[i]);
               
                Vmath::Vadd(nPointsTot,
                            outarray[i], 1,
                            outarray_homo[i], 1,
                            outarray[i], 1);
            }
        }
    }
}
