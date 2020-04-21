///////////////////////////////////////////////////////////////////////////////
//
// File ExpListHomogeneous1D.h
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
// Description: Base class for expansions which are homogeneous in 1
// direction
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLISTHOMO1D_H
#define EXPLISTHOMO1D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <LibUtilities/FFT/NektarFFT.h>
#include <LibUtilities/Communication/Transposition.h>

namespace Nektar
{
    namespace MultiRegions
    {

        enum Homogeneous1DMatType
        {
            eForwardsCoeffSpace1D,
            eBackwardsCoeffSpace1D,
            eForwardsPhysSpace1D,
            eBackwardsPhysSpace1D
        };

        /// A map between homo matrix keys and their associated block
        /// matrices.
        typedef std::map< Homogeneous1DMatType, DNekBlkMatSharedPtr> Homo1DBlockMatrixMap;
        /// A shared pointer to a BlockMatrixMap.
        typedef std::shared_ptr<Homo1DBlockMatrixMap> Homo1DBlockMatrixMapShPtr;

        // Forward declaration for typedefs
        class ExpListHomogeneous1D;

        /// Shared pointer to an ExpList3DHomogeneous1D object.
        typedef std::shared_ptr<ExpListHomogeneous1D>   ExpListHomogeneous1DSharedPtr;
        /// Vector of pointers to ExpList3DHomogeneous1D objects.
        typedef std::vector< ExpListHomogeneous1DSharedPtr > ExpListHomogeneous1DVector;

        /// Abstraction of a two-dimensional multi-elemental expansion which
        /// is merely a collection of local expansions.
        class ExpListHomogeneous1D: public ExpList
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT ExpListHomogeneous1D();

            MULTI_REGIONS_EXPORT ExpListHomogeneous1D(const LibUtilities::SessionReaderSharedPtr &pSession,const LibUtilities::BasisKey &HomoBasis, const NekDouble lz, const bool useFFT, const bool dealiasing);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ExpListHomogeneous1D(const ExpListHomogeneous1D &In);

            MULTI_REGIONS_EXPORT ExpListHomogeneous1D(const ExpListHomogeneous1D &In,
                                                      const std::vector<unsigned int> &eIDs);

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpListHomogeneous1D();

            MULTI_REGIONS_EXPORT void Homogeneous1DTrans(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD, NekDouble> &outarray, 
                                                         bool IsForwards, 
                                                         bool Shuff = true,
                                                         bool UnShuff = true);

            inline void HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                            Array<OneD, NekDouble> &outarray, 
                                            bool Shuff = true,
                                            bool UnShuff = true);

            inline void HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray,
                                            bool Shuff = true,
                                            bool UnShuff = true);

            LibUtilities::BasisSharedPtr  GetHomogeneousBasis(void)
            {
                return m_homogeneousBasis;
            }

            MULTI_REGIONS_EXPORT void PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                                Array<OneD, NekDouble> &out_d0,
                                                Array<OneD, NekDouble> &out_d1,
                                                Array<OneD, NekDouble> &out_d2);

            MULTI_REGIONS_EXPORT void PhysDeriv(Direction edir,
                                                const Array<OneD, const NekDouble> &inarray,
                                                Array<OneD, NekDouble> &out_d);

            ExpListSharedPtr &GetPlane(int n)
            {
                return m_planes[n];
            }

            LibUtilities::TranspositionSharedPtr      m_transposition;
            LibUtilities::CommSharedPtr               m_StripZcomm;
        protected:

            /// FFT variables
            bool                                    m_useFFT;
            LibUtilities::NektarFFTSharedPtr        m_FFT;

            LibUtilities::NektarFFTSharedPtr        m_FFT_deal;

            Array<OneD,NekDouble>                   m_tmpIN;
            Array<OneD,NekDouble>                   m_tmpOUT;

            /// Definition of the total number of degrees of freedom and
            /// quadrature points. Sets up the storage for \a m_coeff and \a
            ///  m_phys.
            LibUtilities::BasisSharedPtr    m_homogeneousBasis;
            NekDouble                       m_lhom;  ///< Width of homogeneous direction
            Homo1DBlockMatrixMapShPtr       m_homogeneous1DBlockMat;
            Array<OneD, ExpListSharedPtr>   m_planes;

            std::unordered_map<int, int>  m_zIdToPlane;
            
            DNekBlkMatSharedPtr GenHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype) const;
            
            DNekBlkMatSharedPtr GetHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype) const;

            NekDouble GetSpecVanVisc(const int k)
            {
                NekDouble returnval = 0.0;

                if(m_specVanVisc.size())
                {
                    returnval = m_specVanVisc[k];
                }

                return returnval;
            }

            //  virtual functions
            virtual void v_SetHomo1DSpecVanVisc(Array<OneD, NekDouble> visc)
            {
                m_specVanVisc = visc;
            }

            virtual int v_GetNumElmts(void)
            {
                return m_planes[0]->GetExpSize();
            }

            virtual LibUtilities::BasisSharedPtr  v_GetHomogeneousBasis(void)
            {
                return GetHomogeneousBasis();
            }

            virtual void v_FwdTrans(const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray);

            virtual void v_FwdTrans_IterPerExp(const Array<OneD,const NekDouble> &inarray,
                                               Array<OneD,      NekDouble> &outarray);

            virtual void v_FwdTrans_BndConstrained(const Array<OneD,const NekDouble> &inarray,
                                               Array<OneD,      NekDouble> &outarray);

            virtual void v_BwdTrans(const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray);
            
            virtual void v_BwdTrans_IterPerExp(const Array<OneD,const NekDouble> &inarray,
                                               Array<OneD,      NekDouble> &outarray);
            
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD, NekDouble> &outarray);
            
            virtual void v_IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray, 
                                                      Array<OneD, NekDouble> &outarray);

            virtual std::vector<LibUtilities::FieldDefinitionsSharedPtr> v_GetFieldDefinitions(void);

            virtual void v_GetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef);

            virtual void v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata);

            virtual void v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractDataToCoeffs(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field, Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractCoeffsToCoeffs(
                                                 const std::shared_ptr<ExpList> &fromExpList, const Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs);

            virtual void v_WriteVtkPieceData(std::ostream &outfile, int expansion,
                                             std::string var);

            virtual void v_PhysInterp1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray);

            virtual void v_PhysGalerkinProjection1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray);


            virtual void v_HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                               Array<OneD, NekDouble> &outarray, 
                                               bool Shuff = true,
                                               bool UnShuff = true);
            
            virtual void v_HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                               Array<OneD, NekDouble> &outarray, 
                                               bool Shuff = true,
                                               bool UnShuff = true);

            virtual void v_DealiasedProd(const Array<OneD, NekDouble> &inarray1,
                                         const Array<OneD, NekDouble> &inarray2,
                                         Array<OneD, NekDouble> &outarray);

            virtual void v_DealiasedDotProd(
                const Array<OneD, Array<OneD, NekDouble> > &inarray1,
                const Array<OneD, Array<OneD, NekDouble> > &inarray2,
                Array<OneD, Array<OneD, NekDouble> > &outarray);

            virtual void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2);

            virtual void v_PhysDeriv(Direction edir,
                                     const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD, NekDouble> &out_d);

            virtual LibUtilities::TranspositionSharedPtr v_GetTransposition(void);

            virtual Array<OneD, const unsigned int> v_GetZIDs(void);

            virtual ExpListSharedPtr &v_GetPlane(int n)
            {
                return GetPlane(n);
            }

            virtual NekDouble v_GetHomoLen(void);

            virtual void v_SetHomoLen(const NekDouble lhom);

        private:

            //Padding operations variables
            bool m_dealiasing;
            int m_padsize;

            /// Spectral vanishing Viscosity coefficient for stabilisation
            Array<OneD, NekDouble> m_specVanVisc;
        };

        inline void ExpListHomogeneous1D::HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                              Array<OneD, NekDouble> &outarray, 
                                                              bool Shuff,
                                                              bool UnShuff)
        {
            v_HomogeneousFwdTrans(inarray,outarray,Shuff,UnShuff);
        }
    
        inline void ExpListHomogeneous1D::HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                              Array<OneD, NekDouble> &outarray, 
                                                              bool Shuff,
                                                              bool UnShuff)
        {
            v_HomogeneousBwdTrans(inarray,outarray,Shuff,UnShuff);
        }

    } //end of namespace
} //end of namespace

#endif//EXPLISTHOMO1D_H
