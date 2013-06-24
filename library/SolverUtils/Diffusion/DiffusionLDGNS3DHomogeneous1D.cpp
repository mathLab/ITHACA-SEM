///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDGNS.cpp
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
// Description: LDGNS diffusion 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionLDGNS3DHomogeneous1D.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLDGNS3DHomogeneous1D ::type = GetDiffusionFactory().
            RegisterCreatorFunction("LDGNS3DHomogeneous1D",
                                    DiffusionLDGNS3DHomogeneous1D ::create);
        
        /**
         * @brief DiffusionLDGNS3DHomogeneous1D uses the 2D WeakDG approach 
         * to compute the diffusion term looping on the planes in the z
         * direction and adding the flux in z direction at the end.
         */
        DiffusionLDGNS3DHomogeneous1D ::DiffusionLDGNS3DHomogeneous1D ()
        {
            string diffName = "LDGNS";
            m_planeDiff = GetDiffusionFactory().CreateInstance(
                                                            diffName, diffName);
        }
        
        /**
         * @brief Initiliase DiffusionLDGNS3DHomogeneous1D objects and store them
         * before starting the time-stepping.
         *
         * This routine calls the virtual functions #v_SetupMetrics,
         * #v_SetupCFunctions and #v_SetupInterpolationMatrices to
         * initialise the objects needed by DiffusionLFR.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLDGNS3DHomogeneous1D ::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            int nConvectiveFields = pFields.num_elements();
            
            Array<OneD, MultiRegions::ExpListSharedPtr>
            pFields_plane0(nConvectiveFields);
            
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                pFields_plane0[i] = pFields[i]->GetPlane(0);
            }
            
            m_planeDiff->InitObject(pSession, pFields_plane0);
            
            spaceDim   = 3;
            nPointsTot = pFields[0]->GetTotPoints();
            nCoeffs    = pFields[0]->GetNcoeffs();
            
            planes     = pFields[0]->GetZIDs();
            num_planes = planes.num_elements();
            
            nPointsTot_plane = nPointsTot/num_planes;
            nCoeffs_plane    = nCoeffs/num_planes;
        }
        
        /**
         * @brief Calculate WeakDG Diffusion for the Navier-Stokes (NS) equations
         * using an LDG interface flux and the the flux in the third direction.
         *
         * The equations that need a diffusion operator are those related
         * with the velocities and with the energy.
         */
        void DiffusionLDGNS3DHomogeneous1D ::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {            
            int i, j, k;
            int nVariables = inarray.num_elements();
               
            DiffusionLDGNSSharedPtr diffLDGNS = boost::dynamic_pointer_cast<
                DiffusionLDGNS>(m_planeDiff);
            
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nConvectiveFields);
            for (j = 0; j < nConvectiveFields; j ++)
            {
                fluxvector[j] = Array<OneD, NekDouble>(nPointsTot);
            }
            
            Array<OneD, Array<OneD, NekDouble> > derivatives01_homo(nVariables);
            for (i = 0; i < nVariables; ++i)
            {
                derivatives01_homo[i] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                fields[0]->PhysDeriv(2, inarray[i], derivatives01_homo[i]);
            }
            
            Array <OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                fields_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                inarray_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                derivatives01_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                outarray_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                fluxvector_homo(num_planes);
            
            m_planeDiff->SetRiemannSolver(m_riemann);
            m_planeDiff->SetFluxVectorVecNS(m_fluxVectorNS);
            
            for (i = 0; i < num_planes; ++i)
            {
                fields_plane[i] = Array<OneD, MultiRegions::ExpListSharedPtr>
                    (nConvectiveFields);
                inarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                    (nVariables);
                derivatives01_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                    (nVariables);
                outarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                    (nConvectiveFields);
                fluxvector_homo[i] =
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(spaceDim);
                
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    fields_plane[i][j]= fields[j]->GetPlane(i);
                        outarray_plane[i][j] = Array<OneD, NekDouble>
                    (nPointsTot_plane, 0.0);
                    
                }
                
                for (j = 0; j < nVariables; j ++)
                {
                    inarray_plane[i][j] = Array<OneD, NekDouble>
                        (nPointsTot_plane, 0.0);
                    Vmath::Vcopy(nPointsTot_plane,
                                 &inarray[j][i * nPointsTot_plane], 1,
                                 &inarray_plane[i][j][0], 1);
                    
                    derivatives01_plane[i][j] = Array<OneD, NekDouble>
                        (nPointsTot_plane, 0.0);
                    Vmath::Vcopy(nPointsTot_plane,
                                 &derivatives01_homo[j][i * nPointsTot_plane], 1,
                                 &derivatives01_plane[i][j][0], 1);
                }
                
                for (j = 0; j < spaceDim ; ++j)
                {
                    fluxvector_homo[i][j] =
                        Array<OneD, Array<OneD, NekDouble> > 
                            (nConvectiveFields);
                    
                    for (int k = 0; k < nConvectiveFields ; ++k)
                    {
                        fluxvector_homo[i][j][k] =
                            Array<OneD, NekDouble>(nPointsTot_plane, 0.0);
                    }
                }
                
                // Set the first order derivatives in 3rd direction 
                diffLDGNS->SetHomoDerivs(derivatives01_plane[i]);
                
                m_planeDiff->Diffuse(nConvectiveFields,
                                     fields_plane[i],
                                     inarray_plane[i],
                                     outarray_plane[i]);
                
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &outarray_plane[i][j][0], 1,
                                 &outarray[j][i * nPointsTot_plane], 1);
                }
                
                m_planeDiff->FluxVec(fluxvector_homo[i]);
                
                for ( j = 0; j < nConvectiveFields; ++j)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &fluxvector_homo[i][2][j][0], 1,
                                 &fluxvector[j][i * nPointsTot_plane], 1);
                }
            }
            
            Array<OneD, Array<OneD, NekDouble> > outarray_homo(nConvectiveFields);
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                outarray_homo[i] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                
                fields[0]->PhysDeriv(2, fluxvector[i], outarray_homo[i]);
                
                Vmath::Vadd(nPointsTot,
                            outarray[i], 1,
                            outarray_homo[i], 1,
                            outarray[i], 1);
            }    
        }
    }//end of namespace SolverUtils
}//end of namespace Nektar}
