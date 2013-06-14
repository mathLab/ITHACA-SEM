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
        
        void AdvectionWeakDG3DHomogeneous1D::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            int num;
            int i, j, k;
            int nVel = advVel.num_elements();
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            
            Array<OneD, unsigned int> planes;
            planes = fields[0]->GetZIDs();
            int num_planes = planes.num_elements();
            
            int nPointsTot_plane = nPointsTot/num_planes;
            int nCoeffs_plane = nCoeffs/num_planes;
             
            m_planeAdv->SetRiemannSolver(m_riemann);
            m_planeAdv->SetFluxVectorVec(m_fluxVector);
            
            
            Array<OneD, Array<OneD, NekDouble> > fluxvector(nConvectiveFields);
            for (j = 0; j < nConvectiveFields; j ++)
            {
                fluxvector[j] = Array<OneD, NekDouble>(nPointsTot);
            }
            
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                fluxvector_homo(num_planes);
            Array<OneD, Array<OneD, NekDouble> >
                outarray_homo(nConvectiveFields);
            Array <OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                fields_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                inarray_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                outarray_plane(num_planes);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                advVel_plane(num_planes);
            
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
                advVel_plane[i] = Array<OneD, Array<OneD, NekDouble> >(nVel);
                
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
                
                for (j = 0; j < nVel; j ++)
                {
                    advVel_plane[i][j] = Array<OneD, NekDouble>
                                                        (nPointsTot_plane, 0.0);
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
                    fluxvector_homo[i][j] =
                        Array<OneD, Array<OneD, NekDouble> >(nVel);
                    
                    for (int k = 0; k < nVel ; ++k)
                    {
                        fluxvector_homo[i][j][k] =
                            Array<OneD, NekDouble>(nPointsTot_plane, 0.0);
                    }
                    
                    Vmath::Vcopy(nPointsTot_plane,
                                 &outarray_plane[i][j][0], 1,
                                 &outarray[j][i * nPointsTot_plane], 1);
                }
                
                m_fluxVector(inarray_plane[i], fluxvector_homo[i]);
                
                for ( j = 0; j < nConvectiveFields; ++j)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &fluxvector_homo[i][j][2][0], 1,
                                 &fluxvector[j][i * nPointsTot_plane], 1);
                }
            }
            
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
    }
}
