///////////////////////////////////////////////////////////////////////////////
//
// File Interp.cpp
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
// Description: Definition of Interpolation methods 
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Basis.h>

namespace Nektar
{    
    namespace LibUtilities
    {     

        // Physical Space Interpolation methods


        // 1D Interpolation
        void Interp1D(const BasisKey &fbasis0, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0, 
                      Array<OneD, NekDouble> &to)
        {
            Interp1D(fbasis0.GetPointsKey(),from,tbasis0.GetPointsKey(),to);
        }
        
        void Interp1D(const PointsKey &fpoints0, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0, 
                      Array<OneD, NekDouble> &to)
        {
            if(fpoints0 == tpoints0) //check to see if the same
            {
                Vmath::Vcopy(fpoints0.GetNumPoints(),from,1,to,1);
            }
            else // interpolate
            {
                DNekMatSharedPtr I0;
                
                I0 = PointsManager()[fpoints0]->GetI(tpoints0);
                
                NekVector<NekDouble> in(fpoints0.GetNumPoints(),from,eWrapper);
                NekVector<NekDouble> out(tpoints0.GetNumPoints(),to,eWrapper);
                
                out  = (*I0)*in;
            }
        }

        void Interp1D(const BasisKey &fbasis0, 
                      const NekDouble *from,  
                      const BasisKey &tbasis0, 
                      NekDouble *to)
        {
            Interp1D(fbasis0.GetPointsKey(),from,tbasis0.GetPointsKey(),to);
        }
        
        void Interp1D(const PointsKey &fpoints0, 
                      const NekDouble *from,  
                      const PointsKey &tpoints0, 
                      NekDouble *to)
        {
            if(fpoints0 == tpoints0) //check to see if the same
            {
                Vmath::Vcopy(fpoints0.GetNumPoints(),from,1,to,1);
            }
            else // interpolate
            {
                
                DNekMatSharedPtr I0;
                
                I0 = PointsManager()[fpoints0]
                    ->GetI(tpoints0);
                
                Blas::Dgemv('N', tpoints0.GetNumPoints(), fpoints0.GetNumPoints(), 
                            1.0, I0->GetPtr().get(), tpoints0.GetNumPoints(), 
                            from, 1, 0.0, to, 1);
            }
        }

        // 2D Interpolation
        void Interp2D(const BasisKey &fbasis0, 
                      const BasisKey &fbasis1, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0,
                      const BasisKey &tbasis1,
                      Array<OneD, NekDouble> &to)
        {
            Interp2D(fbasis0.GetPointsKey(),fbasis1.GetPointsKey(),from.data(),
                     tbasis0.GetPointsKey(),tbasis1.GetPointsKey(),to.data());
        }

        void Interp2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      Array<OneD, NekDouble> &to)
        {
            Interp2D(fpoints0,fpoints1,from.data(),tpoints0,tpoints1,to.data());
        }

        void Interp2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const NekDouble *from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      NekDouble *to)
        {
            // default interpolation
            if((fpoints0 == tpoints0)&&(fpoints1 == tpoints1))
            {
                Vmath::Vcopy(tpoints0.GetNumPoints()*tpoints1.GetNumPoints(),
                             from,1,to,1);
                return;
            }

            DNekMatSharedPtr I0,I1;
            Array<OneD, NekDouble> wsp(tpoints1.GetNumPoints()*fpoints0.GetNumPoints()); // fnp0*tnp1

            int fnp0 = fpoints0.GetNumPoints();
            int fnp1 = fpoints1.GetNumPoints();
            int tnp0 = tpoints0.GetNumPoints();
            int tnp1 = tpoints1.GetNumPoints();

            if(fpoints1 == tpoints1)
            {
                Vmath::Vcopy(fnp0*tnp1,from,1,wsp.get(),1);
            }
            else
            {
                I1 = PointsManager()[fpoints1]->GetI(tpoints1);
                Blas::Dgemm('N', 'T', fnp0, tnp1, fnp1, 1.0, from, fnp0,
                             I1->GetPtr().get(), tnp1, 0.0,  wsp.get(), fnp0);     
            }

            if(fpoints0 == tpoints0)
            {
                Vmath::Vcopy(tnp0*tnp1,wsp.get(),1,to,1);
            }
            else
            {
                I0 = PointsManager()[fpoints0]->GetI(tpoints0);
                Blas::Dgemm('N', 'N', tnp0, tnp1, fnp0, 1.0, I0->GetPtr().get(),
                            tnp0, wsp.get(), fnp0, 0.0, to, tnp0);     
            }
        }



        // 3D interpolation
        void Interp3D(const BasisKey &fbasis0, 
                      const BasisKey &fbasis1, 
                      const BasisKey &fbasis2, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0,
                      const BasisKey &tbasis1,
                      const BasisKey &tbasis2,
                      Array<OneD, NekDouble> &to)
        {
            Interp3D(fbasis0.GetPointsKey(),fbasis1.GetPointsKey(),
                     fbasis2.GetPointsKey(),from.data(),
                     tbasis0.GetPointsKey(),tbasis1.GetPointsKey(),
                     tbasis2.GetPointsKey(),to.data());
        }

        
        void Interp3D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const PointsKey &fpoints2, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      const PointsKey &tpoints2,
                      Array<OneD, NekDouble> &to)
        {
            Interp3D(fpoints0,fpoints1,fpoints2,from.data(),
                     tpoints0,tpoints1,tpoints2,to.data());
        }

        void Interp3D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1,
                      const PointsKey &fpoints2,  
                      const NekDouble *from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      const PointsKey &tpoints2,
                      NekDouble *to)
        {
            int i;
            DNekMatSharedPtr I0, I1, I2;
            
            int fnp0 = fpoints0.GetNumPoints();
            int fnp1 = fpoints1.GetNumPoints();
            int fnp2 = fpoints2.GetNumPoints();
            int tnp0 = tpoints0.GetNumPoints();
            int tnp1 = tpoints1.GetNumPoints();
            int tnp2 = tpoints2.GetNumPoints();
            
            Array<OneD, NekDouble> wsp1(tnp0*tnp1*fnp2);
            Array<OneD, NekDouble> wsp2(tnp0*fnp1*fnp2);
            
            I0 = PointsManager()[fpoints0]->GetI(tpoints0);
            I1 = PointsManager()[fpoints1]->GetI(tpoints1);
            I2 = PointsManager()[fpoints2]->GetI(tpoints2);   

            Blas::Dgemm('N', 'N', tnp0, fnp1*fnp2, fnp0, 1.0, I0->GetPtr().get(),
                        tnp0, from, fnp0, 0.0, wsp2.get(), tnp0);

            for(i = 0; i < fnp2; i++)
            {
                Blas::Dgemm('N', 'T', tnp0,  tnp1,  fnp1, 1.0, wsp2.get()+i*tnp0*fnp1, 
                            tnp0, I1->GetPtr().get(), tnp1, 0.0, wsp1.get()+i*tnp0*tnp1, tnp0);
            }
            
            Blas::Dgemm('N', 'T', tnp0*tnp1,  tnp2,  fnp2, 1.0, wsp1.get(), 
                        tnp0*tnp1, I2->GetPtr().get(), tnp2, 0.0, to, tnp0*tnp1);         
        }



    } // end of namespace
} // end of namespace

