///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalC0ContMap.h
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
// Description: C0-continuous Local to Global mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_LOCALTOGLOBALC0CONTMAP_H
#define MULTIREGIONS_LOCALTOGLOBALC0CONTMAP_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalBaseMap.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>


#include <LocalRegions/PointExp.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace MultiRegions
    {
        const static map<int,int> NullIntIntMap;
        
    class LocalToGlobalC0ContMap: public LocalToGlobalBaseMap
        {
        public:            

            LocalToGlobalC0ContMap(); 
            
            LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                   const StdRegions::StdExpansionVector &locExpVector);

            // Constructor for the 1D expansion mappings
            LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                   const StdRegions::StdExpansionVector &locExpVector, 
                                   const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const map<int,int>& periodicVerticesId);

            // Constructor for the 2D expansion mappings
            LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                   const StdRegions::StdExpansionVector &locExpVector, 
                                   const Array<OneD, const ExpList1DSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const map<int,int>& periodicVerticesId,
                                   const map<int,int>& periodicEdgesId);

            // Constructor for the 3D expansion mappings
            LocalToGlobalC0ContMap(const int numLocalCoeffs, 
                                   const StdRegions::StdExpansionVector &locExpVector, 
                                   const Array<OneD, const ExpList2DSharedPtr> &bndCondExp,
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions,
                                   const map<int,int>& periodicVerticesId,
                                   const map<int,int>& periodicEdgesId,
                                   const map<int,int>& periodicFacesId);
            
            ~LocalToGlobalC0ContMap(); 
                
            int GetLocalToGlobalMap(const int i) const
            {
                return m_localToGlobalMap[i];
            }
            
            inline const Array<OneD,const int>&  GetLocalToGlobalMap(void)
            {
                return m_localToGlobalMap;
            }

            NekDouble GetLocalToGlobalSign(const int i) const
            {
                if(m_signChange)
                {
                    return m_localToGlobalSign[i];
                }
                else
                {
                    return 1.0;
                }
            }

            inline int GetNumLocalCoeffs() const
            {
                return m_numLocalCoeffs;
            }

            inline int GetNumGlobalCoeffs() const
            {
                return m_numGlobalCoeffs;
            }
            
            inline const void LocalToGlobal(const Array<OneD, const NekDouble>& loc, 
                                      Array<OneD, NekDouble>& global) const
            {
                if(m_signChange)
                {
                    Vmath::Scatr(m_numLocalCoeffs, m_localToGlobalSign.get(), loc.get(), m_localToGlobalMap.get(), global.get());
                }
                else
                {                
                    Vmath::Scatr(m_numLocalCoeffs, loc.get(), m_localToGlobalMap.get(), global.get());
                }
            }
            
            inline const void GlobalToLocal(const Array<OneD, const NekDouble>& global, 
                                      Array<OneD, NekDouble>& loc) const
            {
                if(m_signChange)
                {
                    Vmath::Gathr(m_numLocalCoeffs, m_localToGlobalSign.get(), global.get(), m_localToGlobalMap.get(), loc.get());
                }
                else
                {
                    Vmath::Gathr(m_numLocalCoeffs, global.get(), m_localToGlobalMap.get(), loc.get());
                }
            }
            
            inline const void Assemble(const Array<OneD, const NekDouble> &loc, 
                                 Array<OneD, NekDouble> &global) const
            {
                ASSERTL1(loc.get() != global.get(),"Local and Global Arrays cannot be the same");

                Vmath::Zero(m_numGlobalCoeffs, global.get(), 1);
                
                if(m_signChange)
                {
                    Vmath::Assmb(m_numLocalCoeffs, m_localToGlobalSign.get(), loc.get(), m_localToGlobalMap.get(), global.get());
                }
                else
                {
                    Vmath::Assmb(m_numLocalCoeffs, loc.get(), m_localToGlobalMap.get(), global.get());
                }
            } 

            inline int GetFullSystemBandWidth() const
            {
                return m_fullSystemBandWidth;
            }
            
        protected:
            int m_numLocalCoeffs;      //< number of local coefficients
            int m_numGlobalCoeffs;     // Total number of global coefficients
            Array<OneD,int> m_localToGlobalMap;  //< integer map of local coeffs to global space 
            Array<OneD,NekDouble> m_localToGlobalSign; //< integer sign of local coeffs to global space 

            int m_fullSystemBandWidth;
 
        private:
            void SetUp1DExpansionC0ContMap(const int numLocalCoeffs, 
                                           const StdRegions::StdExpansionVector &locExpVector, 
                                           const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp = 
                                               LocalRegions::NullPointExpSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions = 
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const map<int,int>& periodicVerticesId = NullIntIntMap);
            
            void SetUp2DExpansionC0ContMap(const int numLocalCoeffs, 
                                           const StdRegions::StdExpansionVector &locExpVector, 
                                           const Array<OneD, const MultiRegions::ExpList1DSharedPtr> &bndCondExp = 
                                               NullExpList1DSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions = 
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const map<int,int>& periodicVerticesId = NullIntIntMap,
                                           const map<int,int>& periodicEdgesId = NullIntIntMap);

            void SetUp3DExpansionC0ContMap(const int numLocalCoeffs, 
                                           const StdRegions::StdExpansionVector &locExpVector, 
                                           const Array<OneD, const ExpList2DSharedPtr> &bndCondExp = 
                                               NullExpList2DSharedPtrArray,
                                           const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndConditions = 
                                               SpatialDomains::NullBoundaryConditionShPtrArray,
                                           const map<int,int>& periodicVerticesId = NullIntIntMap,
                                           const map<int,int>& periodicEdgesId = NullIntIntMap,
                                           const map<int,int>& periodicFacesId = NullIntIntMap);

            void CalculateBndSystemBandWidth(const StdRegions::StdExpansionVector &locExpVector);


            void CalculateFullSystemBandWidth(const StdRegions::StdExpansionVector &locExpVector);
        };
        typedef boost::shared_ptr<LocalToGlobalC0ContMap>  LocalToGlobalC0ContMapSharedPtr;
        
    } // end of namespace
} // end of namespace

#endif //MULTIREGIONS_LOCALTOGLOBALC0CONTMAP_H

/**
* $Log: LocalToGlobalC0ContMap.h,v $
* Revision 1.6  2009/01/06 21:04:42  sherwin
* Constified GlobalToLocal, LocalToGlobal and Assemble calls that take input and output argument
*
* Revision 1.5  2008/12/17 17:08:53  pvos
* Performance updates
*
* Revision 1.4  2008/11/05 16:15:24  pvos
* Added bandwith calculation routines
*
* Revision 1.3  2008/10/04 19:53:04  sherwin
* Added check that input and output arrays are different
*
* Revision 1.2  2008/09/17 13:46:40  pvos
* Added LocalToGlobalC0ContMap for 3D expansions
*
* Revision 1.1  2008/09/16 13:36:06  pvos
* Restructured the LocalToGlobalMap classes
*
**/

