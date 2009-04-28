///////////////////////////////////////////////////////////////////////////////
//
// File Expansion1D.h
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
// Description: Header file for Expansion1D routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION1D_H

#include <LocalRegions/Expansion.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        
        class Expansion1D: public Expansion
        {
        public:
            
#if 0             
            void AddHDGHelmholtzMatrixBoundaryTerms(const NekDouble tau, 
                                              const Array<OneD, const NekDouble> &inarray,
                                              Array<OneD,NekDouble> &outarray);
#endif
            
            void AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                           const Array<OneD, const NekDouble> &inarray,                                           Array<OneD,NekDouble> &outarray);

            void AddNormTraceInt(const int dir, 
                                 Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,NekDouble> &outarray);

            protected:
            DNekMatSharedPtr GenMatrix(const StdRegions::StdMatrixKey &mkey);
            
            private:
            // Do not add members here since it may lead to conflicts.
            // Only use this class for member functions
            
            
            virtual void v_AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                     const Array<OneD,const NekDouble> &inarray,  Array<OneD,NekDouble> &outarray)
            {
                AddHDGHelmholtzTraceTerms(tau,inarray,outarray);
            }
 
            virtual void v_GetBoundaryMap(Array<OneD, unsigned int> &maparray)
            {
                NEKERROR(ErrorUtil::efatal,"Method does not exist for this shape" );
            }
            

            virtual int v_GetCoordim(void)
            {
                NEKERROR(ErrorUtil::efatal,  "Methods not valid in this class");        
                return -1;
            }

            virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }
        };
        
        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<Expansion1D> Expansion1DSharedPtr;
        typedef std::vector< Expansion1DSharedPtr > Expansion1DVector;
        typedef std::vector< Expansion1DSharedPtr >::iterator Expansion1DVectorIter;
        
    } //end of namespace
} //end of namespace

#define EXPANSION1D_H
#endif

/** 
 *    $Log: Expansion1D.h,v $
 *    Revision 1.2  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
