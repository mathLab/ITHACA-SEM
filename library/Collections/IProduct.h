///////////////////////////////////////////////////////////////////////////////
//
// File: IProduct.h
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
// Description: IProduct operators in 2D for multiple calls in different operators
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_COLLECTIONS_IPRODUCT_H
#define NEKTAR_LIBRARY_COLLECTIONS_IPRODUCT_H

namespace Nektar
{
namespace Collections
{

void QuadIProduct(bool colldir0, bool colldir1, int numElmt,
                  int  nquad0,   int  nquad1,
                  int  nmodes0,  int  nmodes1,
                  const Array<OneD, const NekDouble> &base0,
                  const Array<OneD, const NekDouble> &base1,
                  const Array<OneD, const NekDouble> &jac,
                  const Array<OneD, const NekDouble> &input,
                  Array<OneD, NekDouble> &output,
                  Array<OneD, NekDouble> &wsp);

void TriIProduct(bool sortTopVertex, int numElmt, int  nquad0,
                 int nquad1, int nmodes0, int  nmodes1,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp);

void HexIProduct(bool colldir0, bool colldir1, bool colldir2, int numElmt,
                 int  nquad0,   int  nquad1,   int  nquad2,
                 int  nmodes0,  int  nmodes1,  int  nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp);

void PrismIProduct(bool sortTopVert, int numElmt,
                 int nquad0,  int  nquad1,  int nquad2,
                 int nmodes0, int  nmodes1, int nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                   Array<OneD, NekDouble> &wsp);


void PyrIProduct(bool sortTopVert, int numElmt,
                 int nquad0,  int  nquad1,  int nquad2,
                 int nmodes0, int  nmodes1, int nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                   Array<OneD, NekDouble> &wsp);


void TetIProduct(bool sortTopEdge, int numElmt,
                 int  nquad0,   int  nquad1,   int  nquad2,
                 int  nmodes0,  int  nmodes1,  int  nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp);

}
}
#endif
