///////////////////////////////////////////////////////////////////////////////
//
// File ExpList.h
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
// Description: Expansion list top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H
#define NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H

#include <MultiRegions/MultiRegions.hpp>
#include <StdRegions/StdExpansion.h>

#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
  namespace MultiRegions
  {

    class ExpList
    {
    public:
      ExpList();

      ExpList(const ExpList &in);
      virtual ~ExpList();
      
      inline int GetNcoeffs(void) const
      {
    return m_ncoeffs;
      }
      
      inline int GetPointsTot(void) const
      {
    return m_npoints;
      }
      
      inline void SetTransState(const TransState transState)
      {
    m_transState = transState;
      }

      inline TransState GetTransState(void) const 
      {
          return m_transState; 
      }
      
      inline void SetPhys(ConstArray<OneD, NekDouble> &inarray)
      {
          Vmath::Vcopy(m_npoints,&inarray[0],1,&m_phys[0],1);
      }

      inline void SetPhysState(const bool physState)
      {
    m_physState = physState;
      }

      inline bool GetPhysState(void) const
      {
          return m_physState;
      }
      
      NekDouble PhysIntegral (void);
      void   IProductWRTBase (const ExpList &Sin);
      void   FwdTrans        (const ExpList &Sin);
      void   BwdTrans        (const ExpList &Sin); 
      void   PhysDeriv       (ExpList &S0, ExpList &S1, ExpList &S2); 
      
      void   GeneralMatrixOp(const StdRegions::MatrixType     mtype,
                             const ConstArray<OneD,NekDouble> &inarray,
                             Array<OneD, NekDouble>          &outarray,
                             NekDouble lambda = 1.0);

      void   GetCoords(Array<OneD, NekDouble> &coord_0,
               Array<OneD, NekDouble> &coord_1 = NullNekDouble1DArray,
               Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);

      void   WriteToFile(std::ofstream &out);
    
      DNekScalBlkMatSharedPtr  SetupBlockMatrix(StdRegions::MatrixType mtype, NekDouble scalar = 0.0);

      inline int GetCoordim(int eid)
      {
          ASSERTL2(eid <= (*m_exp).size(),"eid is larger than number of elements");
          
          return (*m_exp)[eid]->GetCoordim();
      }
      

      inline const ConstArray<OneD, NekDouble> &GetCoeffs() const 
      {
          return m_coeffs;
      }

      inline const ConstArray<OneD, NekDouble> &GetPhys()  const
      {
          return m_phys;
      }


      NekDouble Linf (const ExpList &Sol);
      NekDouble L2   (const ExpList &Sol);


      inline int GetExpSize(void)
      {
      return (*m_exp).size();
      }
      

      inline StdRegions::StdExpansionSharedPtr& GetExp(int n)
      {
      return (*m_exp)[n];
      }
      
    protected:
      int m_ncoeffs; 
      int m_npoints;
      Array<OneD, NekDouble> m_coeffs;
      Array<OneD, NekDouble> m_phys;

      TransState m_transState;
      bool       m_physState;
     
      boost::shared_ptr<StdRegions::StdExpansionVector> m_exp;
      
      inline Array<OneD, NekDouble> &UpdateCoeffs()
      {
          m_transState = eLocal;
          return m_coeffs;
      }

      inline Array<OneD, NekDouble> &UpdatePhys()
      {
          m_physState = true;
          return m_phys;
      }

      NekDouble PhysIntegral(const ConstArray<OneD, NekDouble> &inarray);
      void   IProductWRTBase(const ConstArray<OneD, NekDouble> &inarray, 
                 Array<OneD, NekDouble> &outarray);
      
      void   PhysDeriv(const ConstArray<OneD, NekDouble> &inarray,
               Array<OneD, NekDouble> &out_d0, 
               Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
               Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);
      void   FwdTrans (const ConstArray<OneD, NekDouble> &inarray,
               Array<OneD, NekDouble> &outarray);

      void   BwdTrans (const ConstArray<OneD, NekDouble> &inarray, 
               Array<OneD, NekDouble> &outarray); 
      
    private:

    };

    static ExpList NullExpList;
    
  } //end of namespace
} //end of namespace
  
#endif // EXPLIST_H

/**
* $Log: ExpList.h,v $
* Revision 1.23  2007/07/26 08:40:49  sherwin
* Update to use generalised i/o hooks in Helmholtz1D
*
* Revision 1.22  2007/07/22 23:04:20  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.21  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.20  2007/07/17 07:11:05  sherwin
* Chaned definition of NullExpList
*
* Revision 1.19  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.18  2007/07/13 09:02:24  sherwin
* Mods for Helmholtz solver
*
* Revision 1.17  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.16  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
