///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LocalRegions/TetExp.cpp,v $ 
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/TetExp.h>

namespace Nektar
{
  namespace LocalRegions 
  {
    
    TetExp::TetExp(const StdRegions::BasisKey &Ba, 
           const StdRegions::BasisKey &Bb, 
           const StdRegions::BasisKey &Bc):
      StdRegions::StdTetExp(Ba,Bb,Bc)
    {
    }
  
  
    TetExp::TetExp(const StdRegions::BasisKey &Ba, 
           const StdRegions::BasisKey &Bb, 
           const StdRegions::BasisKey &Bc, double *coeffs, 
           double *phys):
      StdRegions::StdTetExp(Ba,Bb,Bc,coeffs,phys)
    {
    }
  
  
    TetExp::TetExp(TetExp &T):StdRegions::StdTetExp(T)
    {  
    }
    
  
    TetExp::~TetExp()
    {  
    }
  
  }//end of namespace
}//end of namespace

/** 
 *    $Log: TetExp.cpp,v $
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.9  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/



