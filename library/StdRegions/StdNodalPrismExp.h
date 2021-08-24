///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTetExp.h
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
// Description: Header field for nodal tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDNODALPRISMEXP_H
#define STDNODALPRISMEXP_H

#include <StdRegions/StdPrismExp.h>

namespace Nektar
{
    namespace StdRegions
    {
        class StdNodalPrismExp: public StdPrismExp
        {
        public:
            STD_REGIONS_EXPORT StdNodalPrismExp();
            STD_REGIONS_EXPORT StdNodalPrismExp(
                const LibUtilities::BasisKey &Ba, 
                const LibUtilities::BasisKey &Bb,
                const LibUtilities::BasisKey &Bc,
                const LibUtilities::PointsType Ntype);
            STD_REGIONS_EXPORT StdNodalPrismExp(const StdNodalPrismExp &T);
            STD_REGIONS_EXPORT ~StdNodalPrismExp();

            //-------------------------------
            // Nodal basis specific routines
            //-------------------------------
            STD_REGIONS_EXPORT void NodalToModal(
                const Array<OneD, const NekDouble>& inarray, 
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT void NodalToModalTranspose(
                const Array<OneD, const NekDouble>& inarray, 
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT void ModalToNodal(
                const Array<OneD, const NekDouble>& inarray, 
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT void GetNodalPoints(
                Array<OneD, const NekDouble> &x,
                Array<OneD, const NekDouble> &y,
                Array<OneD, const NekDouble> &z);
            STD_REGIONS_EXPORT DNekMatSharedPtr GenNBasisTransMatrix();

        protected:
            LibUtilities::PointsKey m_nodalPointsKey;

            STD_REGIONS_EXPORT virtual const LibUtilities::PointsKey v_GetNodalPointsKey() const
            {
                return m_nodalPointsKey;
            };


            STD_REGIONS_EXPORT virtual bool  v_IsNodalNonTensorialExp();
            
            //---------------------------------------
            // Transforms
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_BwdTrans(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFac(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_FwdTrans(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);

            
            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray,
                bool                                mult = true);
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(
                const int                           dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase_SumFac(
                const int                           dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            
            
            //---------------------------------------
            // Evaluation functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_FillMode(
                const int               mode,
                Array<OneD, NekDouble> &outarray);

            
            //---------------------------------------
            // Mapping functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual int  v_GetVertexMap(const int localVertexId,
                                                           bool useCoeffPacking = false);
            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(
                Array<OneD, unsigned int>& outarray);
            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(
                Array<OneD, unsigned int>& outarray);


            //---------------------------------------
            // Wrapper functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                const StdMatrixKey &mkey);
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                const StdMatrixKey &mkey);
        };
        
        typedef std::shared_ptr<StdNodalPrismExp> StdNodalPrismExpSharedPtr;
    } //end of namespace
} //end of namespace
#endif //STDNODALTETEXP_H
