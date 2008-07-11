///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionContField2D.cpp
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
// Description: Advection Field definition for 2D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AdvectionContField2D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        AdvectionContField2D::AdvectionContField2D(SpatialDomains::MeshGraph2D &graph2D,
                                                   SpatialDomains::BoundaryConditions &bcs):
            m_fields(bcs.GetNumVariables())
        {
            int i;
            for( i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i] = MemoryManager<ContField2D>::AllocateSharedPtr(graph2D,bcs,i);
            }
        }


        void AdvectionContField2D::SetInitialConditions(SpatialDomains::BoundaryConditions &bcs, int initialtime)
        {
            int i,j;
            int nq = m_fields[0]->GetPointsTot();

            Array<OneD,NekDouble> x0(nq);
            Array<OneD,NekDouble> x1(nq);
            Array<OneD,NekDouble> x2(nq);

            m_fields[0]->GetCoords(x0,x1,x2);

            for(i = 0 ; i < m_fields.num_elements(); i++)
            {
                SpatialDomains::ConstInitialConditionShPtr ifunc = bcs.GetInitialCondition(i);
                for(j = 0; j < nq; j++)
                {
                    ((m_fields[i])->UpdatePhys())[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
                    (m_fields[i])->SetPhysState(true);
                }
            }
        }

        void AdvectionContField2D::LinearAdvectionOperation(const Array<OneD,const NekDouble>& a, const Array<OneD,const NekDouble>& b)
        {
            int i;
            int nq = m_fields[0]->GetPointsTot();
            
            for(i = 0 ; i < m_fields.num_elements(); i++)
            {
                Array<OneD, NekDouble> fluxvector0(nq);
                Array<OneD, NekDouble> fluxvector1(nq);
            
                Vmath::Smul(nq,a[i],m_fields[i]->GetPhys(),1,fluxvector0,1);
                Vmath::Smul(nq,b[i],m_fields[i]->GetPhys(),1,fluxvector1,1);

                m_fields[i]->PhysDeriv(fluxvector0,fluxvector0);
                m_fields[i]->PhysDeriv(fluxvector1,NullNekDouble1DArray,fluxvector1);

                Vmath::Vadd(nq,fluxvector0,1,fluxvector1,1,fluxvector0,1);
                Vmath::Smul(nq,-1.0,fluxvector0,1,m_fields[i]->UpdatePhys(),1);
                m_fields[i]->SetPhysState(true);
            }
            
        }

        void AdvectionContField2D::FwdTrans(const AdvectionContField2D &In)
        {  
            int i;
            
            for(i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i]->FwdTrans(*m_fields[i]);
            }

        }

        void AdvectionContField2D::BwdTrans(const AdvectionContField2D &In)
        {  
            int i;
            
            for(i = 0 ; i < m_fields.num_elements(); i++)
            {
                m_fields[i]->BwdTrans(*m_fields[i]);
            }

        }
        


    } // end of namespace
} //end of namespace
