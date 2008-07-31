///////////////////////////////////////////////////////////////////////////////
//
// File StdQuadExp.h
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
// Description: Header field for Quadrilateral routines built upon
// StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STDQUADEXP_H
#define NEKTAR_LIB_STDREGIONS_STDQUADEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdExpMap.h>

namespace Nektar
{
    namespace StdRegions
    {

    class StdQuadExp: public StdExpansion2D
        {

        public:

            StdQuadExp();

            /** \brief Constructor using BasisKey class for quadrature
             *  points and order definition 
             */
            StdQuadExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb);

            /** \brief Copy Constructor */
            StdQuadExp(const StdQuadExp &T);

            /** \brief Destructor */
            ~StdQuadExp();

            /** \brief Return Shape of region, using ShapeType enum list. 
             *  i.e. Quadrilateral
             */
            ExpansionType DetExpansionType() const
            {
                return eQuadrilateral;
            };


            virtual bool v_IsBoundaryInteriorExpansion()
            {
                bool returnval = false;
                
                if((m_base[0]->GetBasisType() == LibUtilities::eModified_A)
                   ||(m_base[0]->GetBasisType() == LibUtilities::eGLL_Lagrange))
                {
                    if((m_base[1]->GetBasisType() == LibUtilities::eModified_A)
                       ||(m_base[1]->GetBasisType() == LibUtilities::eGLL_Lagrange))
                    {
                        returnval = true;
                    }
                }
                
                return returnval;
            }


            /** \brief Fill outarray with mode \a mode of expansion
             *
             *    Note for quadrilateral expansions _base[0] (i.e. p)  modes run 
             *  fastest
             */
            void FillMode(const int mode, Array<OneD, NekDouble> &outarray);

            //////////////////////////////
            // Integration Methods
            //////////////////////////////

            NekDouble Integral(const Array<OneD, const NekDouble>& inarray);

            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                                inarray,outarray,1);
            }

            void IProductWRTDerivBase(const int dir, 
                                      const Array<OneD, const NekDouble>& inarray, 
                                      Array<OneD, NekDouble> & outarray);
            
            /** \brief Calculate the inner product of inarray with respect to
             *  the basis B=base0*base1 and put into outarray
             *
             *  \f$ 
             *  \begin{array}{rcl}
             *  I_{pq} = (\phi_q \phi_q, u) & = & \sum_{i=0}^{nq_0}
             *  \sum_{j=0}^{nq_1}
             *  \phi_p(\xi_{0,i}) \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} 
             *  \xi_{1,j}) \\
             *  & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
             *  \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} 
             *  \end{array}
             *  \f$ 
             *
             *  where
             *
             *  \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
             *
             *  which can be implemented as
             *
             *  \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) 
             *  \tilde{u}_{i,j} = {\bf B_1 U}  \f$
             *  \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} = 
             *  {\bf B_0 F}  \f$
             */
            void IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                                 const Array<OneD, const NekDouble>& base1,
                                 const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray,
                                 int coll_check);

            //----------------------------------
            // Generate Matrix Routine
            //----------------------------------

            DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey);

            //----------------------------
            // Differentiation Methods
            //----------------------------

            /** \brief Calculate the derivative of the physical points 
             *
             *  For quadrilateral region can use the Tensor_Deriv function
             *  defined under StdExpansion.
             */
            void PhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                           Array<OneD, NekDouble> &out_d0, 
                           Array<OneD, NekDouble> &out_d1,
                           Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);

            void PhysDeriv(const int dir, 
                           const Array<OneD, const NekDouble>& inarray,
                           Array<OneD, NekDouble> &outarray);

            void StdPhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                              Array<OneD, NekDouble> &out_d0,
                              Array<OneD, NekDouble> &out_d1,
                              Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray, out_d0, out_d1);
            }

            void StdPhysDeriv (const int dir, 
                               const Array<OneD, const NekDouble>& inarray, 
                               Array<OneD, NekDouble> &outarray)
            {
                PhysDeriv(dir,inarray,outarray);
            }

            //----------------------------
            // Evaluations Methods
            //-----------------------------

            void BwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);
            void FwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            NekDouble PhysEvaluate(Array<OneD, const NekDouble>& coords)
            {
                return  StdExpansion2D::PhysEvaluate(coords); 
            }

            void GetBoundaryMap(Array<OneD, unsigned int>& outarray);

            void GetInteriorMap(Array<OneD, unsigned int>& outarray);

            int GetVertexMap(const int localVertexId);
 
            void GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                    Array<OneD, unsigned int> &maparray);

            void GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                     Array<OneD, unsigned int> &maparray,
                                     Array<OneD, int> &signarray);

            void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true);
            void WriteCoeffsToFile(std::ofstream &outfile);

            int GetEdgeNcoeffs(const int i) const
            {
                ASSERTL2((i > 0)&&(i < 3),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasisNumModes(0);
                }
                else
                {
                    return  GetBasisNumModes(1); 
                }
            }

            int GetEdgeNumPoints(const int i) const
            {
                ASSERTL2((i > 0)&&(i < 3),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetNumPoints(0);
                }
                else
                {
                    return  GetNumPoints(1); 
                }
            }

            const LibUtilities::BasisType GetEdgeBasisType(const int i) const
            {
                ASSERTL2((i >= 0)&&(i <= 3),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasisType(0);
                }
                else
                {
                    return  GetBasisType(1);
                }

            }


            const LibUtilities::BasisKey DetEdgeBasisKey(const int i) const
            {
                ASSERTL2((i >= 0)&&(i <= 3),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasis(0)->GetBasisKey();
                }
                else
                {
                    return  GetBasis(1)->GetBasisKey();
                }

            }


            void GetCoords(Array<OneD, NekDouble> &coords_0, 
                           Array<OneD, NekDouble> &coords_1);

            void LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray);

            void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const double lambda);            

        protected:


        private:

            virtual int v_GetNverts() const
            {
                return 4;
            }

            virtual int v_GetNedges() const
            {
                return 4;
            }

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }

            virtual int v_GetEdgeNumPoints(const int i) const
            {
                return GetEdgeNumPoints(i);
            }

            virtual int v_NumBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                         GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                         GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");

                return 4 + 2*(GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
            } 

            virtual int v_NumDGBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                         GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                         GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");

                return  2*GetBasisNumModes(0) + 2*GetBasisNumModes(1);
            } 

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                return GetEdgeBasisType(i);
            }

            
            virtual const LibUtilities::BasisKey v_DetEdgeBasisKey(const int i) const
            {
                return DetEdgeBasisKey(i);
            }

            virtual ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }

            virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                     Array<OneD, NekDouble> &coords_1,
                                     Array<OneD, NekDouble> &coords_2)
            {
                GetCoords(coords_0,coords_1);
            }


            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &array)
            {
                FillMode(mode,array);
            }

            virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray )
            {
                return Integral(inarray);
            }

            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            virtual void v_IProductWRTDerivBase (const int dir, 
                                                 const Array<OneD, const NekDouble> &inarray, 
                                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTDerivBase(dir,inarray,outarray);
            }

            virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }

            virtual void v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray,out_d0, out_d1);
            }

            virtual void v_PhysDeriv(const int dir, 
                                     const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &outarray)
            {
                PhysDeriv(dir,inarray,outarray);
            }

            virtual void v_StdPhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                                        Array<OneD, NekDouble> &out_d0,
                                        Array<OneD, NekDouble> &out_d1,
                                        Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                StdPhysDeriv(inarray, out_d0, out_d1);
            }

            virtual void   v_StdPhysDeriv (const int dir, 
                                           const Array<OneD, const NekDouble>& inarray, 
                                           Array<OneD, NekDouble> &outarray)
            {
                StdPhysDeriv(dir,inarray, outarray);                
            }

            virtual void v_BwdTrans(const Array<OneD, const NekDouble>& inarray, 
                                    Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray, outarray);
            }

            virtual void v_FwdTrans(const Array<OneD, const NekDouble>& inarray, 
                                    Array<OneD, NekDouble> &outarray)
            {
                FwdTrans(inarray, outarray);
            }

            virtual NekDouble v_PhysEvaluate(Array<OneD, const NekDouble>& Lcoords)
            {
                return PhysEvaluate(Lcoords);
            }

            virtual void v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
            {
                GetBoundaryMap(outarray);
            }

            virtual void v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
            {
                GetInteriorMap(outarray);
            }

            virtual int v_GetVertexMap(const int localVertexId)
            {
                return GetVertexMap(localVertexId);
            }
 
            virtual void v_GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                              Array<OneD, unsigned int> &maparray)
            {
                GetEdgeInteriorMap(eid,edgeOrient,maparray);
            }

            virtual void v_GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                               Array<OneD, unsigned int> &maparray,
                                               Array<OneD, int> &signarray)
            {
                GetEdgeToElementMap(eid,edgeOrient,maparray,signarray);
            }

            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true)
            {
                WriteToFile(outfile,format,dumpVar);
            }

            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {
                WriteCoeffsToFile(outfile);
            }

            virtual void v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
            {
                LaplacianMatrixOp(inarray,outarray);
            }
            
            virtual void v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const double lambda)
            {
                HelmholtzMatrixOp(inarray,outarray,lambda);
            }  
        };
        typedef boost::shared_ptr<StdQuadExp> StdQuadExpSharedPtr;

    } //end of namespace
} //end of namespace

#endif //STDQUADEXP_H

/**
 * $Log: StdQuadExp.h,v $
 * Revision 1.37  2008/07/31 11:10:15  sherwin
 * Updates for handling EdgeBasisKey for use with DG advection. Depracated GetEdgeBasis and added DetEdgeBasisKey
 *
 * Revision 1.36  2008/07/29 22:21:15  sherwin
 * A bunch of mods for DG advection and separaring the GetGeom calls into GetGeom1D ...
 *
 * Revision 1.35  2008/07/19 21:12:54  sherwin
 * Removed MapTo function and made orientation convention anticlockwise in UDG routines
 *
 * Revision 1.34  2008/07/12 16:42:28  sherwin
 * Added GetEdgeNumPoints Routiner
 *
 * Revision 1.33  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.32  2008/07/02 14:08:56  pvos
 * Implementation of HelmholtzMatOp and LapMatOp on shape level
 *
 * Revision 1.31  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.30  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.29  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.28  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.27  2008/04/02 22:18:10  pvos
 * Update for 2D local to global mapping
 *
 * Revision 1.26  2008/03/18 14:15:45  pvos
 * Update for nodal triangular helmholtz solver
 *
 * Revision 1.25  2008/03/12 15:25:09  pvos
 * Clean up of the code
 *
 * Revision 1.24  2008/02/29 19:15:19  sherwin
 * Update for UDG stuff
 *
 * Revision 1.23  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.22  2007/12/06 22:44:47  pvos
 * 2D Helmholtz solver updates
 *
 * Revision 1.21  2007/11/08 16:55:14  pvos
 * Updates towards 2D helmholtz solver
 *
 * Revision 1.20  2007/10/03 11:37:51  sherwin
 * Updates relating to static condensation implementation
 *
 * Revision 1.19  2007/07/20 02:16:54  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.18  2007/07/12 12:55:16  sherwin
 * Simplified Matrix Generation
 *
 * Revision 1.17  2007/07/10 20:41:51  kirby
 * more fixes
 *
 * Revision 1.16  2007/07/10 19:27:57  kirby
 * Update for new matrix structures
 *
 * Revision 1.15  2007/05/31 19:13:12  pvos
 * Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 * Revision 1.14  2007/05/15 05:18:24  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.13  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.12  2007/04/08 03:36:58  jfrazier
 * Updated to use SharedArray consistently and minor reformatting.
 *
 * Revision 1.11  2007/04/06 08:44:43  sherwin
 * Update to make 2D regions work at StdRegions level
 *
 * Revision 1.10  2007/04/05 15:20:11  sherwin
 * Updated 2D stuff to comply with SharedArray philosophy
 *
 * Revision 1.9  2007/04/05 11:40:21  pvincent
 * *** empty log message ***
 *
 * Revision 1.8  2007/01/17 16:36:58  pvos
 * updating doxygen documentation
 *
 * Revision 1.7  2007/01/17 16:05:41  pvos
 * updated doxygen documentation
 *
 * Revision 1.6  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.5  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
 *
 * Revision 1.4  2006/07/02 17:16:18  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.3  2006/06/01 14:13:36  kirby
 * *** empty log message ***
 *
 * Revision 1.2  2006/05/23 15:08:19  jfrazier
 * Minor cleanup to correct compile warnings.
 *
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.34  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.33  2006/03/05 23:17:53  sherwin
 *
 * Corrected to allow MMatrix1D and MMatrix2D to execute properly
 *
 * Revision 1.32  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.31  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.30  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/




