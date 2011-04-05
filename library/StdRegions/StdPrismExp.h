//////////////////////////////////////////////////////////////////////////////
//
// File StdPrismExp.h
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
// Description: Header field for prismatic routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STDPRISMEXP_H
#define NEKTAR_LIB_STDREGIONS_STDPRISMEXP_H

#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {
        namespace StdPrismData
        {
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                int nCoef = 0;
                for( int a = 0; a < Na; ++a )
                {
                    for( int b = 0; b < Nb; ++b )
                    {
                        for( int c = 0; c < Nc - a; ++c )
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }
        
        class StdPrismExp: public StdExpansion3D
        {
        
        public:
        
            STD_REGIONS_EXPORT StdPrismExp();
        
            /** \brief Constructor using BasisKey class for quadrature
             *    points and order definition 
             */
            STD_REGIONS_EXPORT StdPrismExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc);
        
            /** \brief Constructor using BasisKey class for quadrature
             *  points and order definition where m_coeffs and m_phys are all
             *    set. 
             */
            STD_REGIONS_EXPORT StdPrismExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc,
                        double *coeffs, double *phys);
        
            /** \brief Copy Constructor */
            STD_REGIONS_EXPORT StdPrismExp(const StdPrismExp &T);
        
            /** \brief Destructor */
            STD_REGIONS_EXPORT ~StdPrismExp();
        
            /** \brief Return Shape of region, using  ShapeType enum list.
             *  i.e. Prism 
             */
            ExpansionType DetExpansionType() const
            {
                return ePrism;
            }
                                                
            STD_REGIONS_EXPORT void GetFaceToElementMap(const int fid, const FaceOrientation faceOrient,
                                     Array<OneD, unsigned int> &maparray,
                                     Array<OneD, int> &signarray);
                                    
            int GetFaceNcoeffs(const int i) const
            {
                ASSERTL2((i >= 0) && (i <= 4), "face id is out of range");
                if((i == 0))
                {
                    return GetBasisNumModes(0)*GetBasisNumModes(1);
                }
                else if((i == 1) || (i == 3))
                {
                    int P = GetBasisNumModes(0)-1, Q = GetBasisNumModes(2)-1;
                    return Q+1 + (P*(1 + 2*Q - P))/2;
                }
                else
                {
                    return GetBasisNumModes(1)*GetBasisNumModes(2);
                }
            }
            
            const int GetEdgeNcoeffs(const int i) const
            {
                ASSERTL2((i >= 0)&&(i <= 8),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasisNumModes(0);
                }
                else if((i == 1)||(i == 3)||(i==8))
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
                ASSERTL2((i >= 0)&&(i <= 8),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasisType(0);
                }
                else if((i == 1)||(i == 3)||(i==8))
                {
                    return  GetBasisType(1);
                }
                else
                {
                    return GetBasisType(2);
                }

            }

            STD_REGIONS_EXPORT NekDouble Integral3D(const Array<OneD, const NekDouble>& inarray, 
                                 const Array<OneD, const NekDouble>& wx,
                                 const Array<OneD, const NekDouble>& wy,
                                 const Array<OneD, const NekDouble>& wz);
            STD_REGIONS_EXPORT NekDouble Integral(const Array<OneD, const NekDouble>& inarray);
        
            /** \brief  Inner product of \a inarray over region with respect to the 
                expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata() and return in \a outarray 
	
                Wrapper call to StdPrismExp::IProductWRTBase
	
                Input:\n
	
                - \a inarray: array of function evaluated at the physical collocation points
	
                Output:\n
	
                - \a outarray: array of inner product with respect to each basis over region

            */
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata(),inarray,outarray);
            }

            STD_REGIONS_EXPORT void PhysDeriv(const Array<OneD, const NekDouble>& u_physical, 
                           Array<OneD, NekDouble> &out_dxi1, 
                           Array<OneD, NekDouble> &out_dxi2,
                           Array<OneD, NekDouble> &out_dxi3 );

            void StdPhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                              Array<OneD, NekDouble> &out_d0,
                              Array<OneD, NekDouble> &out_d1,
                              Array<OneD, NekDouble> &out_d2)
            {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
            }
                       
            /** \brief Backward tranform for triangular elements
             *
             *  \b Note: That 'r' (base[2]) runs fastest in this element
             */
            STD_REGIONS_EXPORT void BwdTrans(const Array<OneD, const NekDouble>& inarray, 
                          Array<OneD, NekDouble> &outarray);
            STD_REGIONS_EXPORT void FwdTrans( const Array<OneD, const NekDouble>& inarray,  Array<OneD, NekDouble> &outarray);

            /** \brief Single Point Evaluation */
            STD_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& xi);
         
            STD_REGIONS_EXPORT void GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z);
            STD_REGIONS_EXPORT void FillMode(const int mode, Array<OneD, NekDouble> &outarray);        
            STD_REGIONS_EXPORT void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");
            STD_REGIONS_EXPORT void WriteCoeffsToFile(std::ofstream &outfile);
                       
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
                \psi_{p}^{a} (\bar \eta_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{pr}^{b} (\xi_{3k})
                w_i w_j w_k u(\bar \eta_{1,i} \xi_{2,j} \xi_{3,k})	     
                J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\bar \eta_{1,i})
                \sum_{j=0}^{nq_1} \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{pr}^b u(\bar \eta_{1i},\xi_{2j},\xi_{3k})
                J_{i,j,k} \end{array} \f$ \n
            
                where
            
                \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\bar \eta_1) \psi_{q}^a (\xi_2) \psi_{pr}^b (\xi_3) \f$ \n
            
                which can be implemented as \n
                \f$f_{pr} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{pr}^b u(\bar \eta_{1i},\xi_{2j},\xi_{3k})
                J_{i,j,k} = {\bf B_3 U}   \f$ \n
                \f$ g_{q} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{pr} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
                \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k})  = {\bf B_1 G} \f$

            **/
            // Interior prism implementation based on Spen's book page 119. and 608.  
            STD_REGIONS_EXPORT void IProductWRTBase(const Array<OneD, const NekDouble>& bx, 
                                 const Array<OneD, const NekDouble>& by, 
                                 const Array<OneD, const NekDouble>& bz, 
                                 const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> & outarray); 
               
        private:
        
            virtual int v_GetNverts() const
            {
                return 6;
            }
        
            virtual int v_GetNedges() const
            {
                return 9;
            }
        
            virtual int v_GetNfaces() const
            {
                return 5;
            }

            virtual ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }

            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }

            virtual int v_GetFaceNcoeffs(const int i) const
            {
                return GetFaceNcoeffs(i);
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

            virtual int  v_CalcNumberOfCoefficients(const std::vector<unsigned int> &nummodes, int &modes_offset)
            {
                int nmodes = StdRegions::StdPrismData::getNumberOfCoefficients(nummodes[modes_offset],nummodes[modes_offset+1],nummodes[modes_offset+2]);
                modes_offset += 3;
                
                return nmodes;
            }

            virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
            {
                return GetEdgeBasisType(i);
            }

            virtual void v_GetCoords( Array<OneD, NekDouble> &coords_x,
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

           virtual void v_PhysDirectionalDeriv(const Array<OneD, const NekDouble>& inarray,
                                                const Array<OneD, const NekDouble>& direction,
                                                Array<OneD, NekDouble> &outarray)
            {
                ASSERTL0(false,"This method is not defined or valid for this class type");
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
      
            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& Lcoords)
            {
                return PhysEvaluate(Lcoords);
            }

            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v")
            {
                WriteToFile(outfile,format,dumpVar,var);
            }

            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {
                WriteCoeffsToFile(outfile);
            }

        };
        typedef boost::shared_ptr<StdPrismExp> StdPrismExpSharedPtr;
    
    } //end of namespace
} //end of namespace

#endif //STDPRISMEXP_H

/**
 * $Log: StdPrismExp.h,v $
 * Revision 1.24  2009/04/27 21:32:45  sherwin
 * Updated WriteToField method
 *
 * Revision 1.23  2009/04/20 16:11:47  sherwin
 * Mods to handle output and optimise DG work
 *
 * Revision 1.22  2009/01/01 02:38:05  ehan
 * cleaned up the code
 *
 * Revision 1.21  2008/11/24 21:06:36  ehan
 * Added virtual functions for necessary mapping routines (boundary map, vertex map, interior map) of Prism.
 *
 * Revision 1.20  2008/11/17 09:02:38  ehan
 * Added necessary mapping routines
 *
 * Revision 1.19  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.18  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.17  2008/06/16 22:46:19  ehan
 * Populated the function GetFaceToElementMap(..)
 *
 * Revision 1.16  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.15  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.14  2008/05/15 22:41:37  ehan
 * Added WriteToFile() function and its virtual function
 *
 * Revision 1.13  2008/05/15 04:14:48  ehan
 * Added virtual function v_CreatStdMatrix()
 *
 * Revision 1.12  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.11  2008/03/25 08:39:55  ehan
 * Added GetEdgeNcoeffs() and GetEdgeBasisType().
 *
 * Revision 1.10  2008/03/17 10:37:12  pvos
 * Clean up of the code
 *
 * Revision 1.9  2008/01/20 06:09:38  bnelson
 * Fixed visual c++ compile errors.
 *
 * Revision 1.8  2008/01/08 22:48:20  ehan
 * Fixed the call signature of a shadowed virtual function: Added a const qualifier to the passed parameter StdMatrixKey in the virtual function v_GenMatrix().  This enables Nektar to generate the correct standard mass matrix at initialization time.
 *
 * Revision 1.7  2008/01/03 12:32:44  ehan
 * Fixed errors from StdMatrix to StdMatrixKey.
 *
 * Revision 1.6  2007/12/28 23:20:20  ehan
 * Completed implementing and testing following functions:
 * Integral, IProductWRTBase, PhysDeriv. BwdTrans, FwdTrans, and PhysEvaluate.
 *
 * Revision 1.5  2007/07/20 02:16:54  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.4  2007/07/10 21:05:17  kirby
 * even more fixes
 *
 * Revision 1.3  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.2  2006/07/02 17:16:18  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.23  2006/03/06 17:12:46  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.22  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.21  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 **/

