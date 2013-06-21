///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDG3DHomogeneous1D.cpp
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

#include <SolverUtils/Diffusion/DiffusionLDG3DHomogeneous1D.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLDG3DHomogeneous1D ::type = GetDiffusionFactory().
            RegisterCreatorFunction(
                "LDG3DHomogeneous1D", DiffusionLDG3DHomogeneous1D ::create);
        
        /**
         * @brief DiffusionLDG3DHomogeneous1D  uses the 2D WeakDG approach 
         * to compute the diffusion term looping on the planes in the z
         * direction and adding the flux in z direction at the end.
         */
        DiffusionLDG3DHomogeneous1D ::DiffusionLDG3DHomogeneous1D ()
        {
            string diffName = "LDG";
            m_planeDiff = GetDiffusionFactory().CreateInstance(diffName, diffName);
        }
        
        /**
         * @brief Initiliase DiffusionLDG3DHomogeneous1D objects and store 
         * them before starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLDG3DHomogeneous1D ::v_InitObject(
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
            
            spaceDim = 3;
            nPointsTot      = pFields[0]->GetTotPoints();
            nCoeffs         = pFields[0]->GetNcoeffs();
            
            planes = pFields[0]->GetZIDs();
            num_planes = planes.num_elements();
            
            nPointsTot_plane = nPointsTot/num_planes;
            nCoeffs_plane = nCoeffs/num_planes;

        }
        
        /**
         * @brief Calculate WeakDG Diffusion for the linear problems
         * using an LDG interface flux and the the flux in the third direction.
         *
         */
        void DiffusionLDG3DHomogeneous1D ::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int i, j, k;
                        
            Array <OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                fields_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                inarray_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                outarray_plane(num_planes);
            
            for (i = 0; i < num_planes; ++i)
            {
                fields_plane[i] = Array<OneD, MultiRegions::ExpListSharedPtr>
                    (nConvectiveFields);
                inarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                    (nConvectiveFields);
                outarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                    (nConvectiveFields);
                
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    fields_plane[i][j]= fields[j]->GetPlane(i);
                    inarray_plane[i][j] = Array<OneD, NekDouble>
                        (nPointsTot_plane, 0.0);
                    outarray_plane[i][j] = Array<OneD, NekDouble>
                        (nPointsTot_plane, 0.0);
                    
                    Vmath::Vcopy(nPointsTot_plane,
                                 &inarray[j][i * nPointsTot_plane], 1,
                                 &inarray_plane[i][j][0], 1);
                }
                
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
            }
            
            Array<OneD, Array<OneD, NekDouble> > flux(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > flux_homo(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > outarray_z(nConvectiveFields);

            NekDouble beta;
            int Homolen = fields[0]->GetHomoLen();
            
            for (j = 0; j < nConvectiveFields; j++)
            {
                flux[j] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                flux_homo[j] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                outarray_z[j] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                
                // Transform flux in Fourier space
                fields[0]->HomogeneousFwdTrans(inarray[j], flux[j]);
            }
            
            m_transpositionLDG = fields[0]->GetTransposition();

            for (i = 0; i < num_planes; ++i)
            {
                beta = 2*M_PI*(m_transpositionLDG->GetK(i))/Homolen;
                
                for (j = 0; j < nConvectiveFields; j++)
                {
                    // Derivative in Fourier space
                    Vmath::Smul(nPointsTot_plane,
                                beta*beta ,
                                &flux[j][0] + i*nPointsTot_plane, 1,
                                &flux_homo[j][0] + i*nPointsTot_plane, 1);
                }
            }
            
            for  (j = 0; j < nConvectiveFields; ++j)
            {
                // Transform back in physical space
                fields[0]->HomogeneousBwdTrans(flux_homo[j], outarray_z[j]);
                
                Vmath::Vsub(nPointsTot,
                            outarray[j], 1,
                            outarray_z[j], 1,
                            outarray[j], 1);
            }
        }
    }// close namespace SolverUtils
}// close namespace nektar++
