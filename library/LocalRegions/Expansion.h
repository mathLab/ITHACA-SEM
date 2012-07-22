///////////////////////////////////////////////////////////////////////////////
//
// File Expansion.h
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
// Description: Header file for Expansion routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPANSION_H

#include <StdRegions/StdExpansion1D.h>
#include <StdRegions/StdExpansion0D.h>
#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        
        class Expansion : virtual public StdRegions::StdExpansion
        {
            public:
                LOCAL_REGIONS_EXPORT Expansion(); // default constructor. 
                LOCAL_REGIONS_EXPORT virtual ~Expansion() {}

                DNekScalMatSharedPtr GetLocMatrix(const LocalRegions::MatrixKey &mkey)
                {
                    return v_GetLocMatrix(mkey);
                }

                DNekScalMatSharedPtr GetLocMatrix(const StdRegions::MatrixType mtype,
                            const StdRegions::ConstFactorMap &factors = StdRegions::NullConstFactorMap,
                            const StdRegions::VarCoeffMap &varcoeffs = StdRegions::NullVarCoeffMap)
                {
                    MatrixKey mkey(mtype, DetExpansionType(), *this, factors, varcoeffs);
                    return GetLocMatrix(mkey);
                }


            protected:
                virtual DNekScalMatSharedPtr v_GetLocMatrix(const LocalRegions::MatrixKey &mkey)
                {
                    NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
                    return NullDNekScalMatSharedPtr;
                }

            private:

        };
        
        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<Expansion> ExpansionSharedPtr;
        typedef boost::weak_ptr<Expansion> ExpansionWeakPtr;
        typedef std::vector< ExpansionSharedPtr > ExpansionVector;
        typedef std::vector< ExpansionSharedPtr >::iterator ExpansionVectorIter;
        
    } //end of namespace
} //end of namespace

#define EXPANSION_H
#endif

/** 
 *    $Log: Expansion.h,v $
 *    Revision 1.5  2009/11/09 15:43:51  sehunchun
 *    HDG2DManifold Solver with Variable coefficients
 *
 *    Revision 1.4  2009/10/06 09:27:17  cbiotto
 *    Adding GetNedges virtual function
 *
 *    Revision 1.3  2008/12/18 14:08:24  pvos
 *    NekConstants update
 *
 *    Revision 1.2  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
