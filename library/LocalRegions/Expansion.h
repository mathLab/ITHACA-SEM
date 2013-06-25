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
#define EXPANSION_H

#include <StdRegions/StdExpansion.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <vector>

namespace Nektar
{
    namespace LocalRegions 
    {

        class Expansion;
        class MatrixKey;

        // type defines for use of PrismExp in a boost vector
        typedef boost::shared_ptr<Expansion> ExpansionSharedPtr;
        typedef boost::weak_ptr<Expansion> ExpansionWeakPtr;
        typedef std::vector< ExpansionSharedPtr > ExpansionVector;
        typedef std::vector< ExpansionSharedPtr >::iterator ExpansionVectorIter;

        class Expansion : virtual public StdRegions::StdExpansion
        {
            public:
                LOCAL_REGIONS_EXPORT Expansion(); // default constructor. 
                LOCAL_REGIONS_EXPORT virtual ~Expansion();

                LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr GetLocMatrix(const LocalRegions::MatrixKey &mkey);

                LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr GetLocMatrix(const StdRegions::MatrixType mtype,
                            const StdRegions::ConstFactorMap &factors = StdRegions::NullConstFactorMap,
                            const StdRegions::VarCoeffMap &varcoeffs = StdRegions::NullVarCoeffMap);

                LOCAL_REGIONS_EXPORT DNekMatSharedPtr BuildTransformationMatrix(
                    const DNekScalMatSharedPtr &r_bnd, 
                    const StdRegions::MatrixType matrixType);

                LOCAL_REGIONS_EXPORT DNekMatSharedPtr BuildVertexMatrix(
                    const DNekScalMatSharedPtr &r_bnd);

            protected:
                virtual DNekScalMatSharedPtr v_GetLocMatrix(const LocalRegions::MatrixKey &mkey);

                virtual DNekMatSharedPtr v_BuildTransformationMatrix(
                    const DNekScalMatSharedPtr &r_bnd, 
                    const StdRegions::MatrixType matrixType);

                virtual DNekMatSharedPtr v_BuildVertexMatrix(
                    const DNekScalMatSharedPtr &r_bnd);  

            private:

        };

    } //end of namespace
} //end of namespace

#endif
