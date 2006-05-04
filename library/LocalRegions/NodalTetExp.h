///////////////////////////////////////////////////////////////////////////////
//
// File NodalTetExp.h
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
// Description: Header for NodalTetExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTETEXP_H
#define NODALTETEXP_H

#include <StdRegions/StdBasis.h>
#include <StdRegions/NodalBasisManager.h>

#include <LocalRegions/TetExp.h>

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
  namespace LocalRegions 
  {
  
    class NodalTetExp: public TetExp
    {
    public:
      NodalTetExp(const StdRegions::Basis &Ba, const StdRegions::Basis &Bb, 
		  const StdRegions::Basis &Bc):
	TetExp(Ba,Bb,Bc)
      {
      }
    
    protected:

    private:
            
    };

  } //end of namespace
} //end of namespace

#endif // NODALTETEXP_H


/** 
 *    $Log: NodalTetExp.h,v $
 *    Revision 1.7  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.6  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
