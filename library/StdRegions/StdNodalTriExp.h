///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTriExp.h
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
// Description: Header field for Nodal triangle routines built upon
// StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDNODALTRIEXP_H
#define STDNODALTRIEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/LocalRegionsDeclarations.hpp>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {
    class StdNodalTriExp: public StdTriExp
        {            
        public:            
            //Constructors
            STD_REGIONS_EXPORT StdNodalTriExp(void);
            
            STD_REGIONS_EXPORT StdNodalTriExp(const LibUtilities::BasisKey &Ba, 
                           const LibUtilities::BasisKey &Bb,
                           const LibUtilities::PointsType Ntype);
            
            //Copy Constructor
            STD_REGIONS_EXPORT StdNodalTriExp(const StdNodalTriExp &T);

            //Destructor
            STD_REGIONS_EXPORT ~StdNodalTriExp();

            ExpansionType DetExpansionType() const
            {
                return eTriangle;
            }

            /////////////////////////////
            /// Nodal basis specific routines
            ///////////////////////////

            STD_REGIONS_EXPORT void NodalToModal();
            STD_REGIONS_EXPORT void NodalToModal(const Array<OneD, const NekDouble>& inarray, 
                              Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT void NodalToModalTranspose();
            STD_REGIONS_EXPORT void NodalToModalTranspose(const Array<OneD, const NekDouble>& inarray, 
                                       Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT void ModalToNodal();
            STD_REGIONS_EXPORT void ModalToNodal(const Array<OneD, const NekDouble>& inarray, 
                              Array<OneD, NekDouble> &outarray);

            void GetNodalPoints(Array<OneD, const NekDouble> &x, 
                                Array<OneD, const NekDouble> &y)
            {
                LibUtilities::PointsManager()[*m_nodalPointsKey]->GetPoints(x,y);
            }
            
            //////////////////////////////
            // Integration Methods
            //////////////////////////////
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                StdNodalTriExp::IProductWRTBase_SumFac(inarray,outarray);
            }

            void IProductWRTDerivBase(const int dir, 
                                      const Array<OneD, const NekDouble>& inarray, 
                                      Array<OneD, NekDouble> & outarray)
            {
                StdNodalTriExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
            }
            
            /** \brief Fill outarray with nodal mode \a mode of expansion
             *   and put in m_phys
             */
            STD_REGIONS_EXPORT void FillMode(const int mode, Array<OneD, NekDouble> &outarray);
            
            STD_REGIONS_EXPORT void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");
            
            //-----------------------------
            // Evaluations Methods
            //-----------------------------
            
            void BwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray)
            {
                StdNodalTriExp::BwdTrans_SumFac(inarray,outarray);
            }

            STD_REGIONS_EXPORT void FwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT void GetBoundaryMap(Array<OneD, unsigned int>& outarray);

            STD_REGIONS_EXPORT void GetInteriorMap(Array<OneD, unsigned int>& outarray);

            int GetVertexMap(const int localVertexId)
            {
                ASSERTL0((localVertexId>=0)&&(localVertexId<=2),
                         "Local Vertex ID must be between 0 and 2");                
                return localVertexId;
            }
 
            STD_REGIONS_EXPORT void GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                    Array<OneD, unsigned int> &maparray,
                                    Array<OneD, int> &signarray);

            STD_REGIONS_EXPORT void GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                     Array<OneD, unsigned int> &maparray,
                                     Array<OneD, int> &signarray);

            void MassMatrixOp(const Array<OneD, const NekDouble> &inarray, 
                              Array<OneD,NekDouble> &outarray,
                              const StdMatrixKey &mkey)
            {              
                StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
            }

            void LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {                
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }

            void LaplacianMatrixOp(const int k1, const int k2, 
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {           
                StdExpansion::LaplacianMatrixOp_MatFree(k1,k2,inarray,outarray,mkey);
            }

            void WeakDerivMatrixOp(const int i,
                                   const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {
                StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
            }
            
            void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const StdMatrixKey &mkey)
            {
                StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }       
            
        protected:            
            boost::shared_ptr<LibUtilities::PointsKey> m_nodalPointsKey;

            void BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
            {
                Array<OneD, NekDouble> tmp(m_ncoeffs);
                NodalToModal(inarray,tmp);
                StdTriExp::BwdTrans_SumFac(tmp,outarray);
            }

            void IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                        Array<OneD, NekDouble> &outarray)
            {
                StdTriExp::IProductWRTBase_SumFac(inarray,outarray);
                NodalToModalTranspose(outarray,outarray);    
            }          
/*             void IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& base0,  */
/*                                         const Array<OneD, const NekDouble>& base1, */
/*                                         const Array<OneD, const NekDouble>& inarray,  */
/*                                         Array<OneD, NekDouble> &outarray) */
/*             { */
/*                 StdTriExp::IProductWRTBase_SumFac(base0,base1,inarray,outarray); */
/*                 NodalToModalTranspose(outarray,outarray);      */
/*             }  */

            void IProductWRTDerivBase_SumFac(const int dir, 
                                             const Array<OneD, const NekDouble>& inarray, 
                                             Array<OneD, NekDouble> & outarray)
            {
                StdTriExp::IProductWRTDerivBase_SumFac(dir,inarray,outarray);
                NodalToModalTranspose(outarray,outarray);
            }
                        
            STD_REGIONS_EXPORT DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT DNekMatSharedPtr GenNBasisTransMatrix();
            
        private:
            
            virtual int v_NumBndryCoeffs() const
            {
                return 3 + (GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
            } 

            virtual ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                return GetEdgeBasisType(i);
            }

            virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey) 
            {
                return GenMatrix(mkey);
            }

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }

            //////////////////////////////
            // Integration Methods
            //////////////////////////////
            
            virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray)
            {
                return Integral(inarray);
            }
            
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray, outarray);
            }

            virtual void v_IProductWRTDerivBase (const int dir, 
                                                 const Array<OneD, const NekDouble> &inarray, 
                                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTDerivBase(dir,inarray,outarray);
            }
            
            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
            {
                FillMode(mode, outarray);
            }
            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            virtual void v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray, out_d0, out_d1);
            }
            
            //-----------------------------
            // Evaluations Methods
            //-----------------------------
            
            virtual void v_BwdTrans(const Array<OneD, const NekDouble>& inarray,
                                    Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray,outarray);
            }
            
            /** \brief Virtual call to StdTriExp::FwdTrans */
            virtual void v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                    Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray,outarray);
            }
            
            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords)
            {
                return StdTriExp::PhysEvaluate(coords);
            }
            
            
            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v")
            {
                WriteToFile(outfile,format,dumpVar,var);
            }

            virtual void v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
            {
                GetBoundaryMap(outarray);
            }

            virtual void v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
            {
                GetInteriorMap(outarray);
            }

            virtual int v_GetVertexMap(int localVertexId)
            {
                return GetVertexMap(localVertexId);
            }
 
            virtual void v_GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                              Array<OneD, unsigned int> &maparray,
                                               Array<OneD, int> &signarray)
            {
                GetEdgeInteriorMap(eid,edgeOrient,maparray,signarray);
            }

            virtual void v_GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                               Array<OneD, unsigned int> &maparray,
                                               Array<OneD, int> &signarray)
            {
                GetEdgeToElementMap(eid,edgeOrient,maparray,signarray);
            }
                       
            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {
                WriteCoeffsToFile(outfile);
            }

            virtual void v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                BwdTrans_SumFac(inarray,outarray);
            }    
            
            virtual void v_IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                                  Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase_SumFac(inarray,outarray);
            }            
            
            virtual void v_IProductWRTDerivBase_SumFac(const int dir,
                                                       const Array<OneD, const NekDouble>& inarray, 
                                                       Array<OneD, NekDouble> &outarray)
            {
                IProductWRTDerivBase_SumFac(dir,inarray,outarray);
            }    


            virtual void v_MassMatrixOp(const Array<OneD, const NekDouble> &inarray, 
                                        Array<OneD,NekDouble> &outarray,
                                        const StdMatrixKey &mkey)
            {
                MassMatrixOp(inarray,outarray,mkey);
            }  
            
            virtual void v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey)
            {
                LaplacianMatrixOp(inarray,outarray,mkey);
            }

            virtual void v_LaplacianMatrixOp(const int k1, const int k2, 
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey)
            {
                LaplacianMatrixOp(k1,k2,inarray,outarray,mkey);
            }

            virtual void v_WeakDerivMatrixOp(const int i,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey)
            {
                WeakDerivMatrixOp(i,inarray,outarray,mkey);
            }
            
            virtual void v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const StdMatrixKey &mkey)
            {
                HelmholtzMatrixOp(inarray,outarray,mkey);
            }  
            
        };        
        typedef boost::shared_ptr<StdNodalTriExp> StdNodalTriExpSharedPtr;
    } //end of namespace
} //end of namespace

#endif //STDNODALTRIEXP_H

/**
 * $Log: StdNodalTriExp.h,v $
 * Revision 1.30  2009/01/21 16:58:39  pvos
 * Added additional geometric factors to improve efficiency
 *
 * Revision 1.29  2008/11/24 10:31:14  pvos
 * Changed name from _PartitionedOp to _MatFree
 *
 * Revision 1.28  2008/11/19 16:02:47  pvos
 * Added functionality for variable Laplacian coeffcients
 *
 * Revision 1.27  2008/11/05 16:08:15  pvos
 * Added elemental optimisation functionality
 *
 * Revision 1.26  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.25  2008/07/19 21:12:54  sherwin
 * Removed MapTo function and made orientation convention anticlockwise in UDG routines
 *
 * Revision 1.24  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.23  2008/06/05 15:06:06  pvos
 * Added documentation
 *
 * Revision 1.22  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.21  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.20  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.19  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.18  2008/04/02 22:18:10  pvos
 * Update for 2D local to global mapping
 *
 * Revision 1.17  2008/03/18 14:15:45  pvos
 * Update for nodal triangular helmholtz solver
 *
 * Revision 1.16  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.15  2007/10/03 11:37:51  sherwin
 * Updates relating to static condensation implementation
 *
 * Revision 1.14  2007/07/22 23:04:27  bnelson
 * Backed out Nektar::ptr.
 *
 * Revision 1.13  2007/07/20 02:16:54  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.12  2007/07/12 12:55:16  sherwin
 * Simplified Matrix Generation
 *
 * Revision 1.11  2007/07/11 06:35:24  sherwin
 * Update after manager reshuffle
 *
 * Revision 1.10  2007/05/31 19:13:12  pvos
 * Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 * Revision 1.9  2007/05/15 05:18:24  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.8  2007/04/18 09:44:01  sherwin
 * Working version for StdNodalTri. Removed lapack.cpp from compile.
 *
 * Revision 1.7  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.6  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.5  2007/01/15 21:13:46  sherwin
 * Nodal stuff correction and added Field Classes
 *
 * Revision 1.4  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.3  2006/06/01 14:46:16  kirby
 * *** empty log message ***
 *
 * Revision 1.2  2006/05/23 15:08:19  jfrazier
 * Minor cleanup to correct compile warnings.
 *
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.16  2006/04/25 20:23:33  jfrazier
 * Various fixes to correct bugs, calls to ASSERT, etc.
 *
 * Revision 1.15  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.14  2006/03/06 17:12:46  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.13  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.12  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 **/

