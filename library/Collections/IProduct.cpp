///////////////////////////////////////////////////////////////////////////////
//
// File: IProduct.cpp
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
// Description: IProduct operators for multiple calls in different operators
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Collection.h>
#include <Collections/IProduct.h>

using namespace std;

namespace Nektar
{
namespace Collections
{

/**
 *
 */
void QuadIProduct(bool colldir0, bool colldir1, int numElmt,
                  int  nquad0,   int  nquad1,
                  int  nmodes0,  int  nmodes1,
                  const Array<OneD, const NekDouble> &base0,
                  const Array<OneD, const NekDouble> &base1,
                  const Array<OneD, const NekDouble> &jac,
                  const Array<OneD, const NekDouble> &input,
                  Array<OneD, NekDouble> &output,
                  Array<OneD, NekDouble> &wsp)
{
    int totpoints = nquad0*nquad1;
    int totmodes  = nmodes0*nmodes1;

    Vmath::Vmul(numElmt*totpoints,jac,1,input,1,wsp,1);

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


/**
 *
 */
void TriIProduct(bool sortTopVertex, int numElmt, int  nquad0,
                 int nquad1, int nmodes0, int  nmodes1,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp)
{
    int totmodes  = LibUtilities::StdTriData::getNumberOfCoefficients(
                                                    nmodes0, nmodes1);
    int totpoints = nquad0*nquad1;

    Vmath::Vmul(numElmt*totpoints,jac,1,input,1,wsp,1);

    Array<OneD, NekDouble> wsp1 = wsp + max(totpoints,totmodes)*numElmt;

    Blas::Dgemm('T','N', nquad1*numElmt,nmodes0,nquad0,1.0,&wsp[0],nquad0,
                base0.get(), nquad0, 0.0, &wsp1[0], nquad1*numElmt);

    int i, mode;
    // Inner product with respect to 'b' direction
    for (mode=i=0; i < nmodes0; ++i)
    {
        Blas::Dgemm('T', 'N', nmodes1-i, numElmt,  nquad1,
                    1.0, base1.get()+mode*nquad1,  nquad1,
                    wsp1.get() + i*nquad1*numElmt, nquad1, 0.0,
                    &output[mode], totmodes);

        mode += nmodes1 - i;
    }

    // fix for modified basis by splitting top vertex mode
    if (sortTopVertex)
    {
        Blas::Dgemv('T', nquad1,numElmt,1.0,wsp1.get()+nquad1*numElmt,nquad1,
                    base1.get()+nquad1,1,1.0, &output[1],totmodes);
    }
}


/**
 *
 */
void HexIProduct(bool colldir0, bool colldir1, bool colldir2, int numElmt,
                 int  nquad0,   int  nquad1,  int nquad2,
                 int  nmodes0,  int  nmodes1, int nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp)
{
    int totmodes  = nmodes0*nmodes1*nmodes2;
    int totpoints = nquad0 *nquad1 *nquad2;


    if(colldir0 && colldir1 && colldir2)
    {

        Vmath::Vmul(numElmt*totpoints,jac,1,input,1,output,1);
    }
    else
    {
        Vmath::Vmul(numElmt*totpoints,jac,1,input,1,wsp,1);

        // Assign second half of workspace for 2nd DGEMM operation.
        Array<OneD, NekDouble> wsp1 = wsp  + totpoints*numElmt;

        // note sure what criterion we should use to swap around these
        // strategies
        if(numElmt < nmodes0 || 1)
        {
            Array<OneD, NekDouble> wsp2 = wsp1 + nmodes0*nquad1*nquad2;

            //loop over elements
            for(int n = 0; n < numElmt; ++n)
            {
                if(colldir0)
                {

                    for(int i = 0; i < nmodes0; ++i)
                    {
                        Vmath::Vcopy(nquad1*nquad2,&wsp[n*totpoints] + i,nquad0,
                                     wsp1.get()+nquad1*nquad2*i,1);
                    }
                }
                else
                {
                    Blas::Dgemm('T', 'N', nquad1*nquad2, nmodes0, nquad0,
                                1.0, &wsp[n*totpoints],  nquad0,
                                base0.get(), nquad0,
                                0.0, wsp1.get(),  nquad1*nquad2);
                }


                if(colldir1)
                {
                    // reshuffle data for next operation.
                    for(int i = 0; i < nmodes1; ++i)
                    {
                        Vmath::Vcopy(nquad2*nmodes0,wsp1.get()+i,nquad1,
                                     wsp2.get()+nquad2*nmodes0*i,1);
                    }
                }
                else
                {
                    Blas::Dgemm('T', 'N', nquad2*nmodes0, nmodes1, nquad1,
                                1.0, wsp1.get(),  nquad1,
                                base1.get(), nquad1,
                                0.0, wsp2.get(),  nquad2*nmodes0);
                }

                if(colldir2)
                {
                    // reshuffle data for next operation.
                    for(int i = 0; i < nmodes2; ++i)
                    {
                        Vmath::Vcopy(nmodes0*nmodes1,wsp2.get()+i,nquad2,
                                     &output[n*totmodes]+nmodes0*nmodes1*i,1);
                    }
                }
                else
                {
                    Blas::Dgemm('T', 'N', nmodes0*nmodes1, nmodes2, nquad2,
                                1.0, wsp2.get(),  nquad2,
                                base2.get(), nquad2,
                                0.0, &output[n*totmodes], nmodes0*nmodes1);
                }
            }
        }
        else
        {
            Array<OneD, NekDouble> wsp2 = wsp1 + numElmt*(max(totpoints,
                                                              totmodes));

            if(colldir0)
            {
                for(int i = 0; i < nquad0; ++i)
                {
                    Vmath::Vcopy(nquad1*nquad2*numElmt,&wsp[i],nquad0,
                                 &wsp1[i*nquad1*nquad2*numElmt],1);
                }
            }
            else
            {
                // large degmm but copy at end.
                Blas::Dgemm('T','N', nquad1*nquad2*numElmt, nmodes0, nquad0,
                            1.0, &wsp[0],  nquad0,  base0.get(),   nquad0,
                            0.0, &wsp1[0], nquad1*nquad2*numElmt);
            }

            if(colldir1)
            {
                for(int i = 0; i < nquad1; ++i)
                {
                    Vmath::Vcopy(nquad2*numElmt*nmodes0,&wsp1[i],nquad1,
                                 &wsp2[i*nquad2*numElmt*nmodes0],1);
                }
            }
            else
            {
                Blas::Dgemm('T','N', nquad2*numElmt*nmodes0,  nmodes1, nquad1,
                            1.0, &wsp1[0],   nquad1, base1.get(),   nquad1,
                            0.0, &wsp2[0],  nquad2*numElmt*nmodes0);
            }


            if(numElmt > 1)
            {
                if(colldir2)
                {
                    for(int i = 0; i < nquad2; ++i)
                    {
                        Vmath::Vcopy(nmodes0*nmodes1,&wsp2[i],nquad2,
                                     &output[i*nmodes0*nmodes1],1);
                    }
                }
                else
                {
                    Blas::Dgemm('T', 'N', numElmt*nmodes0*nmodes1, nmodes2,
                                nquad2, 1.0, &wsp2[0],  nquad2,
                                base2.get(),   nquad2, 0.0,
                                &wsp1[0],  numElmt*nmodes0*nmodes1);
                }

                for(int i = 0; i < totmodes; ++i)
                {
                    Vmath::Vcopy(numElmt, &wsp1[i*numElmt], 1,
                                          &output[i],       totmodes);
                }

            }
            else
            {
                if(colldir2)
                {
                    for(int i = 0; i < nquad2; ++i)
                    {
                        Vmath::Vcopy(nmodes0*nmodes1,&wsp2[i],nquad2,
                                     &output[i*nmodes0*nmodes1],1);
                    }
                }
                else
                {
                    Blas::Dgemm('T','N', numElmt*nmodes0*nmodes1, nmodes2,
                                nquad2, 1.0, &wsp2[0],  nquad2,
                                base2.get(),   nquad2, 0.0,
                                &output[0],  numElmt*nmodes0*nmodes1);
                }
            }
        }
    }
}


/**
 *
 */
void PrismIProduct(bool sortTopVertex, int numElmt,
                 int nquad0,  int  nquad1,  int nquad2,
                 int nmodes0, int  nmodes1, int nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp)
{
    int totmodes  = LibUtilities::StdPrismData::getNumberOfCoefficients(
                                                    nmodes0,nmodes1,nmodes2);
    int totpoints = nquad0*nquad1*nquad2;
    int cnt;
    int mode, mode1;

    Vmath::Vmul(numElmt*totpoints,jac,1,input,1,wsp,1);

    Array<OneD, NekDouble> wsp1 = wsp + numElmt * nquad2
                                                * (max(nquad0*nquad1,
                                                       nmodes0*nmodes1));

    // Perform iproduct  with respect to the  '0' direction
    Blas::Dgemm('T', 'N', nquad1*nquad2*numElmt, nmodes0, nquad0,
                1.0, wsp.get(), nquad0, base0.get(),
                nquad0, 0.0, wsp1.get(), nquad1*nquad2*numElmt);


    // Perform iproduct  with respect to the  '1' direction
    Blas::Dgemm('T', 'N', nquad2*numElmt*nmodes0, nmodes1, nquad1,
                1.0, wsp1.get(), nquad1, base1.get(),
                nquad1, 0.0, wsp.get(), nquad2*numElmt*nmodes0);


    // Inner product with respect to the '2' direction (not sure if it would
    // be better to swap loops?)
    mode = mode1 = cnt = 0;
    for(int i = 0; i < nmodes0; ++i)
    {
        cnt = i*nquad2*numElmt;
        for(int j = 0; j < nmodes1; ++j)
        {
            Blas::Dgemm('T', 'N', nmodes2-i, numElmt, nquad2,
                        1.0, base2.get()+mode*nquad2, nquad2,
                             wsp.get()+j*nquad2*numElmt*nmodes0 + cnt,   nquad2,
                        0.0, output.get()+mode1,    totmodes);
            mode1 += nmodes2-i;
        }
        mode  += nmodes2-i;
    }

    // fix for modified basis by splitting top vertex mode
    if (sortTopVertex)
    {
        // top singular vertex
        // ((1+a)/2 components entry into (1+c)/2)
        // Could be made into an mxv if we have specialised base1[1]
        for(int j =0; j < nmodes1; ++j)
        {
            Blas::Dgemv('T', nquad2,numElmt,1.0,
                        wsp.get()+j*nquad2*numElmt*nmodes0+nquad2*numElmt,
                        nquad2, base2.get()+nquad2,1,1.0,
                        &output[j*nmodes2+1], totmodes);
        }
    }
}


/**
 *
 */
void PyrIProduct(bool sortTopVertex, int numElmt,
                 int nquad0,  int  nquad1,  int nquad2,
                 int nmodes0, int  nmodes1, int nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp)
{
    int totmodes  = LibUtilities::StdPyrData::getNumberOfCoefficients(
                                                    nmodes0,nmodes1,nmodes2);
    int totpoints = nquad0*nquad1*nquad2;
    int cnt;
    int mode, mode1;

    ASSERTL1(wsp.size() >= numElmt*(nquad1*nquad2*nmodes0 +
                                            nquad2*max(nquad0*nquad1,nmodes0*nmodes1)),
             "Insufficient workspace size");

    Vmath::Vmul(numElmt*totpoints,jac,1,input,1,wsp,1);

    Array<OneD, NekDouble> wsp1 = wsp + numElmt * nquad2
                                                * (max(nquad0*nquad1,
                                                       nmodes0*nmodes1));

    // Perform iproduct  with respect to the  '0' direction
    Blas::Dgemm('T', 'N', nquad1*nquad2*numElmt, nmodes0, nquad0,
                1.0, wsp.get(), nquad0, base0.get(),
                nquad0, 0.0, wsp1.get(), nquad1*nquad2*numElmt);

    // Inner product with respect to the '1' direction
    mode = 0;
    for(int i=0; i < nmodes0; ++i)
    {
        Blas::Dgemm('T', 'N', nquad2*numElmt, nmodes1, nquad1,
                    1.0, wsp1.get()+ i*nquad1*nquad2*numElmt, nquad1,
                    base1.get(),           nquad1,
                    0.0, wsp.get() + mode*nquad2*numElmt,nquad2*numElmt);
        mode  += nmodes1;
    }

    // Inner product with respect to the '2' direction
    mode = mode1 = cnt = 0;
    for(int i = 0; i < nmodes0; ++i)
    {
        for(int j = 0; j < nmodes1; ++j, ++cnt)
        {
            int ijmax = max(i,j);
            Blas::Dgemm('T', 'N', nmodes2-ijmax, numElmt, nquad2,
                        1.0, base2.get()+mode*nquad2, nquad2,
                             wsp.get()+cnt*nquad2*numElmt,   nquad2,
                        0.0, output.get()+mode1,    totmodes);
            mode  += nmodes2-ijmax;
            mode1 += nmodes2-ijmax;
        }

        //increment mode in case order1!=order2
        for(int j = nmodes1; j < nmodes2; ++j)
        {
            int ijmax = max(i,j);
            mode += nmodes2-ijmax;
        }
    }

    // fix for modified basis for top singular vertex component
    // Already have evaluated (1+c)/2 (1-b)/2 (1-a)/2
    if(sortTopVertex)
    {
        for(int n = 0; n < numElmt; ++n)
        {
            // add in (1+c)/2 (1+b)/2 component
            output[1+n*totmodes] += Blas::Ddot(nquad2,
                             base2.get()+nquad2,1,
                             &wsp[nquad2*numElmt + n*nquad2],1);

            // add in (1+c)/2 (1-b)/2 (1+a)/2 component
            output[1+n*totmodes] += Blas::Ddot(nquad2,
                             base2.get()+nquad2,1,
                             &wsp[nquad2*nmodes1*numElmt+n*nquad2],1);

            // add in (1+c)/2 (1+b)/2 (1+a)/2 component
            output[1+n*totmodes] += Blas::Ddot(nquad2,
                             base2.get()+nquad2,1,
                             &wsp[nquad2*(nmodes1+1)*numElmt+n*nquad2],1);
        }
    }
}



/**
 *
 */
void TetIProduct(bool sortTopEdge, int numElmt,
                 int nquad0,  int  nquad1,  int nquad2,
                 int nmodes0, int  nmodes1, int nmodes2,
                 const Array<OneD, const NekDouble> &base0,
                 const Array<OneD, const NekDouble> &base1,
                 const Array<OneD, const NekDouble> &base2,
                 const Array<OneD, const NekDouble> &jac,
                 const Array<OneD, const NekDouble> &input,
                 Array<OneD, NekDouble> &output,
                 Array<OneD, NekDouble> &wsp)
{
    int totmodes  = LibUtilities::StdTetData::getNumberOfCoefficients(
                                                    nmodes0,nmodes1,nmodes2);
    int totpoints = nquad0*nquad1*nquad2;
    int cnt;
    int mode, mode1;

    Vmath::Vmul(numElmt*totpoints,jac,1,input,1,wsp,1);

    Array<OneD, NekDouble> wsp1 = wsp +
        nquad2*numElmt*(max(nquad0*nquad1,nmodes0*(2*nmodes1-nmodes0+1)/2));


    // Perform iproduct  with respect to the  '0' direction
    Blas::Dgemm('T', 'N', nquad1*nquad2*numElmt, nmodes0, nquad0,
                1.0, wsp.get(), nquad0, base0.get(),
                nquad0, 0.0, wsp1.get(), nquad1*nquad2*numElmt);

    // Inner product with respect to the '1' direction
    mode = 0;
    for(int i=0; i < nmodes0; ++i)
    {
        Blas::Dgemm('T', 'N', nquad2*numElmt, nmodes1-i, nquad1,
                    1.0, wsp1.get()+ i*nquad1*nquad2*numElmt, nquad1,
                    base1.get() + mode*nquad1,           nquad1,
                    0.0, wsp.get() + mode*nquad2*numElmt,nquad2*numElmt);
        mode  += nmodes1-i;
    }


    // fix for modified basis by splitting top vertex mode
    if (sortTopEdge)
    {
        // base singular vertex and singular edge (1+b)/2
        // ((1+a)/2 components entry into (1+b)/2)
        // Could be made into an mxm if we have specialised base1[1]
        for(int n = 0; n < numElmt; ++n)
        {
            Blas::Dgemv('T', nquad1, nquad2,
                        1.0, wsp1.get()+numElmt*nquad1*nquad2 +
                        n*nquad1*nquad2, nquad1,
                        base1.get()+nquad1,  1, 1.0,
                        wsp.get()+nquad2*numElmt + n*nquad2, 1);
        }
    }

    // Inner product with respect to the '2' direction
    mode = mode1 = cnt = 0;
    for(int i = 0; i < nmodes0; ++i)
    {
        for(int j = 0; j < nmodes1-i; ++j, ++cnt)
        {
            Blas::Dgemm('T', 'N', nmodes2-i-j, numElmt, nquad2,
                        1.0, base2.get()+mode*nquad2, nquad2,
                             wsp.get()+cnt*nquad2*numElmt,   nquad2,
                        0.0, output.get()+mode1,    totmodes);
            mode  += nmodes2-i-j;
            mode1 += nmodes2-i-j;
        }

        //increment mode in case order1!=order2
        mode +=  (nmodes2-nmodes1)*(nmodes2-nmodes1+1)/2;
    }

    // fix for modified basis for top singular vertex component
    // Already have evaluated (1+c)/2 (1-b)/2 (1-a)/2
    if(sortTopEdge)
    {
        for(int n = 0; n < numElmt; ++n)
        {
            // add in (1+c)/2 (1+b)/2 component
            output[1+n*totmodes] += Blas::Ddot(nquad2,
                             base2.get()+nquad2,1,
                             &wsp[nquad2*numElmt + n*nquad2],1);

            // add in (1+c)/2 (1-b)/2 (1+a)/2 component
            output[1+n*totmodes] += Blas::Ddot(nquad2,
                             base2.get()+nquad2,1,
                             &wsp[nquad2*nmodes1*numElmt+n*nquad2],1);
        }
    }

}

}
}
