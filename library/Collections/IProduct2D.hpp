///////////////////////////////////////////////////////////////////////////////
//
// File: IProduct2D.hpp
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

#ifndef NEKTAR_LIBRARY_COLLECTIONS_IPRODUCT2D_HPP
#define NEKTAR_LIBRARY_COLLECTIONS_IPRODUCT2D_HPP

namespace Nektar 
{
    namespace Collections 
    {
        void QuadIProduct(bool colldir0, bool colldir1, int numElmt, 
                          int  nquad0,   int  nquad1, 
                          int  nmodes0,  int  nmodes1, 
                          const Array<OneD, const NekDouble> &base0,
                          const Array<OneD, const NekDouble> &base1,
                          const Array<OneD, const NekDouble> &input, 
                          Array<OneD, NekDouble> &output,
                          Array<OneD, NekDouble> &wsp)
        {
            int totpoints = nquad0*nquad1;
            int totmodes  = nmodes0*nmodes1;
            
            if(colldir0 && colldir1)
            {
                Vmath::Vcopy(numElmt*totmodes,wsp.get(),1,output.get(),1);
            }
            else
            {                 
                Array<OneD, NekDouble> wsp1 = wsp  + max(totpoints,totmodes)*numElmt; 
                if(colldir0)
                {
                    for(int i = 0; i < nquad0; ++i)
                    {
                        Vmath::Vcopy(nquad1*numElmt,&wsp[i],nquad0,
                                     &wsp1[i*nquad1*numElmt],1);
                    }
                }
                else
                {
                    Blas::Dgemm('T','N', nquad1*numElmt,nmodes0,nquad0,1.0,
                                &wsp[0],nquad0, base0.get(), nquad0, 
                                0.0,&wsp1[0], nquad1*numElmt);
                }
                
                
                if(numElmt > 1)
                {
                    
                    if(colldir1)
                    {
                        for(int i = 0; i < nquad1; ++i)
                        {
                            Vmath::Vcopy(numElmt*nmodes0,&wsp1[i],nquad1,
                                             &wsp[i*numElmt*nmodes0],1);
                        }
                    }
                    else
                    {
                        
                        Blas::Dgemm('T','N', numElmt*nmodes0,  nmodes1, nquad1,
                                    1.0, &wsp1[0], nquad1, base1.get(),   nquad1,
                                    0.0, &wsp[0], numElmt*nmodes0);
                    }
                    
                    for(int i = 0; i < totmodes; ++i)
                    {
                        Vmath::Vcopy(numElmt,&wsp[i*numElmt],1,&output[i],totmodes);
                    }
                }
                else
                {
                    if(colldir1)
                    {
                        for(int i = 0; i < nquad1; ++i)
                        {
                            Vmath::Vcopy(numElmt*nmodes0,&wsp1[i],nquad1,
                                         &output[i*numElmt*nmodes0],1);
                        }
                    }
                    else
                    {
                        Blas::Dgemm('T','N', nmodes0,  nmodes1, nquad1,
                                    1.0, &wsp1[0], nquad1, base1.get(), nquad1,
                                    0.0, &output[0], nmodes0);
                    }
                }
            }
        }
    }
}
#endif
