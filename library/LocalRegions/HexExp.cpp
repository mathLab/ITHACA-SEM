///////////////////////////////////////////////////////////////////////////////
//
// File HexExp.cpp
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
// Description: Methods for Hex expansion in local regoins
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/HexExp.h>

namespace Nektar
{
  namespace LocalRegions 
  {

 
    HexExp::HexExp(const StdRegions::BasisKey &Ba, 
		   const StdRegions::BasisKey &Bb, 
		   const StdRegions::BasisKey &Bc):
      StdRegions::StdHexExp(Ba,Bb,Bc)
    {
    }
  
 
    HexExp::HexExp(const StdRegions::BasisKey &Ba,
		   const StdRegions::BasisKey &Bb, 
		   const StdRegions::BasisKey &Bc, double *coeffs, 
		   double *phys):
      StdRegions::StdHexExp(Ba,Bb,Bc,coeffs,phys)
    {
    }
  
 
    HexExp::HexExp(HexExp &T):StdRegions::StdHexExp(T)
    {
    }
  
    HexExp::~HexExp()
    {
    }
        
  }//end of namespace
}//end of namespace

/** 
 *    $Log: HexExp.cpp,v $
 *    Revision 1.8  2006/03/12 07:43:31  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
