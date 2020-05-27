///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion0D.h
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 0d expansion. Typically this inolves physical
// space operations.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef STDEXP0D_H
#define STDEXP0D_H

#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion0D: virtual public StdExpansion
        {

        public:

            STD_REGIONS_EXPORT StdExpansion0D();
            STD_REGIONS_EXPORT StdExpansion0D(int numcoeffs, const LibUtilities::BasisKey &Ba);
            STD_REGIONS_EXPORT StdExpansion0D(const StdExpansion0D &T);
            STD_REGIONS_EXPORT virtual ~StdExpansion0D();
			
            STD_REGIONS_EXPORT void PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                                                    Array<OneD,       NekDouble>& outarray);
			
            
        protected:
            std::map<int, NormalVector> m_vertexNormals;
			
                        STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                                         const Array<OneD, const NekDouble>& coords,
                                         const Array<OneD, const NekDouble>& physvals);
        private:

            // Virtual Functions ----------------------------------------

            virtual int v_GetCoordim(void)
            {
                return 1;
            }

            virtual int v_GetShapeDimension() const
            {
                return 1;
            }

            virtual int v_GetNedges() const
            {
                return 0;
            }

            virtual int v_GetNfaces() const
            {
                return 0;
            }

        };

        typedef std::shared_ptr<StdExpansion0D> StdExpansion0DSharedPtr;

    } //end of namespace
} //end of namespace

#endif //STDEXP0D_H
