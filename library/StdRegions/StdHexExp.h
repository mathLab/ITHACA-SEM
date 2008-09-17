///////////////////////////////////////////////////////////////////////////////
//
// File StdHexExp.h
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
// Description: Hex routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGSION_STDHEXEXP_H
#define NEKTAR_LIBS_STDREGSION_STDHEXEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>


namespace Nektar
{
    namespace StdRegions
    {

    class StdHexExp: public StdExpansion3D
        {

        public:
        
            StdHexExp();

            /** \brief Constructor using BasisKey class for quadrature
             *  points and order definition 
             */
            StdHexExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc);

            /** \brief Constructor using BasisKey class for quadrature
             *  points and order definition where m_coeffs and m_phys 
             *  are all set. 
             */
            StdHexExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc,
                      double *coeffs, double *phys);

            /** \brief Copy Constructor */
            StdHexExp(const StdHexExp &T);

            /** \brief Destructor */
            ~StdHexExp();
            

            /** \brief Return Shape of region, using  ShapeType enum list. 
             *  i.e. Hexahedron
             */
      
            ExpansionType DetExpansionType() const
            {
                return eHexahedron;
            }

            void GetBoundaryMap(Array<OneD, unsigned int>& outarray);
            
            void GetInteriorMap(Array<OneD, unsigned int>& outarray);
            
            int GetVertexMap(const int localVertexId);
 
            void GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                    Array<OneD, unsigned int> &maparray,
                                    Array<OneD, int> &signarray);

            void GetFaceInteriorMap(const int fid, const FaceOrientation faceOrient,
                                    Array<OneD, unsigned int> &maparray,
                                    Array<OneD, int> &signarray);
            
            void GetFaceToElementMap(const int fid, const FaceOrientation faceOrient,
                                     Array<OneD, unsigned int> &maparray,
                                     Array<OneD, int> &signarray);

            int NumBndryCoeffs() const
            {
                ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                         GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                         GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");
                ASSERTL1(GetBasisType(2) == LibUtilities::eModified_A ||
                         GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                         "BasisType is not a boundary interior form");

                int nmodes0 = m_base[0]->GetNumModes();
                int nmodes1 = m_base[1]->GetNumModes();
                int nmodes2 = m_base[2]->GetNumModes();

                return ( 2*( nmodes0*nmodes1 + nmodes0*nmodes2 + nmodes1*nmodes2) -
                         4*( nmodes0 + nmodes1 + nmodes2 ) + 8 );
            }

            int GetFaceNcoeffs(const int i) const
            {
                ASSERTL2((i >= 0) && (i <= 5), "face id is out of range");
                if((i == 0) || (i == 5))
                {
                    return GetBasisNumModes(0)*GetBasisNumModes(1);
                }
                else if((i == 1) || (i == 3))
                {
                    return GetBasisNumModes(0)*GetBasisNumModes(2);
                }
                else
                {
                    return GetBasisNumModes(1)*GetBasisNumModes(2);
                }
            }

            int GetFaceIntNcoeffs(const int i) const
            {
                ASSERTL2((i >= 0) && (i <= 5), "face id is out of range");
                if((i == 0) || (i == 5))
                {
                    return (GetBasisNumModes(0)-2)*(GetBasisNumModes(1)-2);
                }
                else if((i == 1) || (i == 3))
                {
                    return (GetBasisNumModes(0)-2)*(GetBasisNumModes(2)-2);
                }
                else
                {
                    return (GetBasisNumModes(1)-2)*(GetBasisNumModes(2)-2);
                }
            }

            const int GetEdgeNcoeffs(const int i) const
            {
                ASSERTL2((i >= 0)&&(i <= 11),"edge id is out of range");

                if((i == 0)||(i == 2)||(i == 8)||(i == 10))
                {
                    return  GetBasisNumModes(0);
                }
                else if((i == 1)||(i == 3)||(i == 9)||(i == 11))
                {
                    return  GetBasisNumModes(1); 
                }
                else
                {
                    return GetBasisNumModes(2);
                }

            }

            const LibUtilities::BasisType GetEdgeBasisType(const int i) const
            {
                ASSERTL2((i >= 0)&&(i <= 11),"edge id is out of range");

                if((i == 0)||(i == 2)||(i==8)||(i==10))
                {
                    return  GetBasisType(0);
                }
                else if((i == 1)||(i == 3)||(i == 9)||(i == 11))
                {
                    return  GetBasisType(1);
                }
                else
                {
                    return GetBasisType(2);
                }

            }

            /** \brief Fill outarray with mode \a mode of expansion
             *
             *    Note for hexahedral expansions _base[0] (i.e. p)  modes run 
             *  fastest
             */

            void FillMode(const int mode, Array<OneD, NekDouble> &outarray);

            //////////////////////////////
            // Integration Methods
            //////////////////////////////

            NekDouble Integral3D(const Array<OneD, const NekDouble>& inarray, 
                                 const Array<OneD, const NekDouble>& wx,
                                 const Array<OneD, const NekDouble>& wy, 
                                 const Array<OneD, const NekDouble>& wz);
            NekDouble Integral(const Array<OneD, const NekDouble>& inarray);

            /** \brief  Inner product of \a inarray over region with respect to the 
		expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata() and return in \a outarray 
	
		Wrapper call to StdHexExp::IProductWRTBase
	
		Input:\n
	
		- \a inarray: array of function evaluated at the physical collocation points
	
		Output:\n
	
		- \a outarray: array of inner product with respect to each basis over region

            */
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata(),inarray,outarray,1);
            }

            //----------------------------
            // Differentiation Methods
            //----------------------------

            /** \brief Calculate the deritive of the physical points 
             *
             *  For quadrilateral region can use the Tensor_Deriv function
             *  defined under StdExpansion.
             */
            void PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                           Array<OneD, NekDouble> &out_d0,
                           Array<OneD, NekDouble> &out_d1,
                           Array<OneD, NekDouble> &out_d2);

            void StdPhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                              Array<OneD, NekDouble> &out_d0,
                              Array<OneD, NekDouble> &out_d1,
                              Array<OneD, NekDouble> &out_d2)
            {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }
                                   
            void GetCoords(Array<OneD, NekDouble> &coords_0, 
                           Array<OneD, NekDouble> &coords_1, Array<OneD, NekDouble> &coords_2);

            //----------------------------
            // Evaluations Methods
            //---------------------------

            void BwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            void FwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);
            NekDouble PhysEvaluate(Array<OneD, const NekDouble>& coords)
            {
                return  StdExpansion3D::PhysEvaluate(coords);  
            }
            void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true);
            void WriteCoeffsToFile(std::ofstream &outfile);
            
            //----------------------------------
            // Local Matrix Routines
            //----------------------------------
            DNekMatSharedPtr GenMassMatrix();

            DNekMatSharedPtr GenLaplacianMatrix();

            DNekMatSharedPtr GenLaplacianMatrix(const int i, const int j);

            DNekMatSharedPtr GenWeakDerivMatrix(const int i);

            DNekMatSharedPtr GenNBasisTransMatrix();

            DNekMatSharedPtr GenBwdTransMatrix();


            DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
            {
                return StdExpansion::CreateGeneralMatrix(mkey);
            }

        protected:

            /** 
		\brief Calculate the inner product of inarray with respect to
		the basis B=base0*base1*base2 and put into outarray:
		
		\f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
		\sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
		\psi_{p}^{a} (\xi_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{r}^{a} (\xi_{3k})
		w_i w_j w_k u(\xi_{1,i} \xi_{2,j} \xi_{3,k})	     
		J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\xi_{1,i})
		\sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
		J_{i,j,k} \end{array} \f$ \n
		
		where
		
		\f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a ( \xi_1) \psi_{q}^a (\xi_2) \psi_{r}^a (\xi_3) \f$ \n
		
		which can be implemented as \n
		\f$f_{r} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{r}^a u(\xi_{1i},\xi_{2j},\xi_{3k})
		J_{i,j,k} = {\bf B_3 U}   \f$ \n
		\f$ g_{q} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{r} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
		\f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k})  = {\bf B_1 G} \f$

            **/
            void IProductWRTBase(const Array<OneD, const NekDouble>& bx, 
                                 const Array<OneD, const NekDouble>& by, 
                                 const Array<OneD, const NekDouble>& bz, 
                                 const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> & outarray, 
                                 int coll_check);

        private:

            virtual int v_GetNverts() const
            {
                return 8;
            }

            virtual int v_GetNedges() const
            {
                return 12;
            }

            virtual int v_GetNfaces() const
            {
                return 6;
            }

            virtual ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            };

            virtual int v_NumBndryCoeffs() const
            {
                return NumBndryCoeffs();
            } 

            virtual int v_GetFaceNcoeffs(const int i) const
            {
                return GetFaceNcoeffs(i);
            } 

            virtual int v_GetFaceIntNcoeffs(const int i) const
            {
                return GetFaceIntNcoeffs(i);
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
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray)
            {
                GetEdgeInteriorMap(eid,edgeOrient,maparray,signarray);
            }    

            virtual void v_GetFaceInteriorMap(const int fid, const FaceOrientation faceOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int>& signarray)
            {
                GetFaceInteriorMap(fid,faceOrient,maparray,signarray);
            }

            virtual void v_GetFaceToElementMap(const int fid, const FaceOrientation faceOrient,
                                               Array<OneD, unsigned int> &maparray,
                                               Array<OneD, int>& signarray)
            {
                GetFaceToElementMap(fid,faceOrient,maparray,signarray);
            }
            
            virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey) 
            {
                return GenMatrix(mkey);
            }

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                return GetEdgeBasisType(i);
            }

            virtual void v_GetCoords(Array<OneD, NekDouble> &coords_x,
                                     Array<OneD, NekDouble> &coords_y,
                                     Array<OneD, NekDouble> &coords_z)
            {
                GetCoords(coords_x, coords_y, coords_z);
            }
                                              
            virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray )
            {
                return Integral(inarray);
            }
                        
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(inarray, outarray);
            }

            virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
            {
                return FillMode(mode, outarray);
            }

            virtual void v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &out_d0,
                                     Array<OneD, NekDouble> &out_d1,
                                     Array<OneD, NekDouble> &out_d2)
            {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }

            virtual void v_StdPhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                                        Array<OneD, NekDouble> &out_d0,
                                        Array<OneD, NekDouble> &out_d1,
                                        Array<OneD, NekDouble> &out_d2)
            {
                StdPhysDeriv(inarray, out_d0, out_d1, out_d2);
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

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }
            
            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true)            
            {
                WriteToFile(outfile,format,dumpVar);
            }

            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {
                WriteCoeffsToFile(outfile);
            }

        };
        typedef boost::shared_ptr<StdHexExp> StdHexExpSharedPtr;


    } //end of namespace
} //end of namespace

#endif //STDHEXEXP_H

/**
 * $Log: StdHexExp.h,v $
 * Revision 1.26  2008/09/15 13:18:08  pvos
 * Added more hexahedron mapping routines
 *
 * Revision 1.25  2008/09/12 11:26:39  pvos
 * Updates for mappings in 3D
 *
 * Revision 1.24  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.23  2008/06/16 22:45:51  ehan
 * Populated the function GetFaceToElementMap(..)
 *
 * Revision 1.22  2008/06/05 15:06:06  pvos
 * Added documentation
 *
 * Revision 1.21  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.20  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.19  2008/05/15 22:40:05  ehan
 * Clean up the codes
 *
 * Revision 1.18  2008/05/15 04:14:37  ehan
 * Added virtual function v_CreatStdMatrix()
 *
 * Revision 1.17  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.16  2008/03/25 08:39:45  ehan
 * Added GetEdgeNcoeffs() and GetEdgeBasisType().
 *
 * Revision 1.15  2008/03/17 10:37:32  pvos
 * Clean up of the code
 *
 * Revision 1.14  2008/01/20 06:09:37  bnelson
 * Fixed visual c++ compile errors.
 *
 * Revision 1.13  2008/01/08 22:30:43  ehan
 * Clean up the codes.
 *
 * Revision 1.12  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.11  2007/12/01 00:52:32  ehan
 * Completed implementing and testing following functions:
 * Integral, IProductWRTBase, PhysDeriv. BwdTrans, FwdTrans, and PhysEvaluate.
 *
 * Revision 1.10  2007/07/20 02:16:54  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.9  2007/07/10 21:05:16  kirby
 * even more fixes
 *
 * Revision 1.7  2007/01/17 16:36:58  pvos
 * updating doxygen documentation
 *
 * Revision 1.6  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.5  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
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
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.30  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.29  2006/03/06 17:12:45  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.28  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.27  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.26  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 *
 **/



