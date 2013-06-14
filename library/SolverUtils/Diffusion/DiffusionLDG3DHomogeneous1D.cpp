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
            RegisterCreatorFunction("LDG3DHomogeneous1D", DiffusionLDG3DHomogeneous1D ::create);
        
        DiffusionLDG3DHomogeneous1D ::DiffusionLDG3DHomogeneous1D ()
        {
            string diffName = "LDG";
            m_planeDiff = GetDiffusionFactory().CreateInstance(diffName, diffName);
        }
        
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
        }
        
        void DiffusionLDG3DHomogeneous1D ::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int num;
            int i, j, k;
            int spaceDim = 3;
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            
            Array<OneD, unsigned int> planes;
            planes = fields[0]->GetZIDs();
            int num_planes = planes.num_elements();
            
            int nPointsTot_plane = nPointsTot/num_planes;
            int nCoeffs_plane = nCoeffs/num_planes;
            
            Array<OneD, Array<OneD, NekDouble> >
                outarray_homo(nConvectiveFields);
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
            
            for  (i = 0; i < nConvectiveFields; ++i)
            {
                outarray_homo[i] = Array<OneD, NekDouble>(nPointsTot, 0.0);
                
                // to change
                fields[0]->PhysDeriv(2, inarray[i], outarray_homo[i]);
                fields[0]->PhysDeriv(2, outarray_homo[i], outarray_homo[i]);
                
                Vmath::Vsub(nPointsTot,
                            outarray[i], 1,
                            outarray_homo[i], 1,
                            outarray[i], 1);
            }
        }
    }
}
