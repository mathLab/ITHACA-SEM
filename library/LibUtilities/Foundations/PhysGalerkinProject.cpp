///////////////////////////////////////////////////////////////////////////////
//
// File PhysGalerkinProject.cpp
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
// Description: Definition of Physical Space Galerkin Projection methods 
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/Foundations/PhysGalerkinProject.h>
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
        void PhysGalerkinProject1D(const BasisKey &fbasis0, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0, 
                      Array<OneD, NekDouble> &to)
        {
            PhysGalerkinProject1D(fbasis0.GetPointsKey(),from,tbasis0.GetPointsKey(),to);
        }
        
        void PhysGalerkinProject1D(const PointsKey &fpoints0, 
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
                DNekMatSharedPtr GP0;
                
                GP0 = PointsManager()[tpoints0]->GetGalerkinProjection(fpoints0);
                
                NekVector<NekDouble> in(fpoints0.GetNumPoints(),from,eWrapper);
                NekVector<NekDouble> out(tpoints0.GetNumPoints(),to,eWrapper);
                
                GP0->Transpose();
                out  = (*GP0)*in;
            }
        }

        void PhysGalerkinProject1D(const BasisKey &fbasis0, 
                      const NekDouble *from,  
                      const BasisKey &tbasis0, 
                      NekDouble *to)
        {
            PhysGalerkinProject1D(fbasis0.GetPointsKey(),from,tbasis0.GetPointsKey(),to);
        }
        
        void PhysGalerkinProject1D(const PointsKey &fpoints0, 
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
                
                DNekMatSharedPtr GP0;
                
                GP0 = PointsManager()[tpoints0]
                    ->GetGalerkinProjection(fpoints0);
                
                Blas::Dgemv('T', tpoints0.GetNumPoints(), fpoints0.GetNumPoints(), 
                            1.0, GP0->GetPtr().get(), tpoints0.GetNumPoints(), 
                            from, 1, 0.0, to, 1);
            }
        }

        // 2D Interpolation
        void PhysGalerkinProject2D(const BasisKey &fbasis0, 
                      const BasisKey &fbasis1, 
                      const Array<OneD, const NekDouble>& from,  
                      const BasisKey &tbasis0,
                      const BasisKey &tbasis1,
                      Array<OneD, NekDouble> &to)
        {
            PhysGalerkinProject2D(fbasis0.GetPointsKey(),fbasis1.GetPointsKey(),from.data(),
                     tbasis0.GetPointsKey(),tbasis1.GetPointsKey(),to.data());
        }

        void PhysGalerkinProject2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const Array<OneD, const NekDouble>& from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      Array<OneD, NekDouble> &to)
        {
            PhysGalerkinProject2D(fpoints0,fpoints1,from.data(),tpoints0,tpoints1,to.data());
        }

        void PhysGalerkinProject2D(const PointsKey &fpoints0, 
                      const PointsKey &fpoints1, 
                      const NekDouble *from,  
                      const PointsKey &tpoints0,
                      const PointsKey &tpoints1,
                      NekDouble *to)
        {
            DNekMatSharedPtr GP0,GP1;
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
                GP1 = PointsManager()[tpoints1]->GetGalerkinProjection(fpoints1);
                Blas::Dgemm('N', 'T', fnp0, tnp1, fnp1, 1.0, from, fnp0,
                             GP1->GetPtr().get(), tnp1, 0.0,  wsp.get(), fnp0);     
            }

            if(fpoints0 == tpoints0)
            {
                Vmath::Vcopy(tnp0*tnp1,wsp.get(),1,to,1);
            }
            else
            {
                GP0 = PointsManager()[tpoints0]->GetGalerkinProjection(fpoints0);
                Blas::Dgemm('N', 'N', tnp0, tnp1, fnp0, 1.0, 
                            GP0->GetPtr().get(),
                            tnp0, wsp.get(), fnp0, 0.0, to, tnp0);     
            }
        }

    } // end of namespace
} // end of namespace

