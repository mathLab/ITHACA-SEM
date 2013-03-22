///////////////////////////////////////////////////////////////////////////////
//
// File StdPyrExp.h
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
// Description: Header field for pyramidic routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGIONS_STDPYREXP_H
#define NEKTAR_LIBS_STDREGIONS_STDPYREXP_H

#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdPyrExp : virtual public StdExpansion3D
        {
        
        public:
    
            STD_REGIONS_EXPORT StdPyrExp();

            STD_REGIONS_EXPORT StdPyrExp(const LibUtilities::BasisKey &Ba, 
                                         const LibUtilities::BasisKey &Bb, 
                                         const LibUtilities::BasisKey &Bc);

            STD_REGIONS_EXPORT StdPyrExp(const LibUtilities::BasisKey &Ba, 
                                         const LibUtilities::BasisKey &Bb, 
                                         const LibUtilities::BasisKey &Bc,
                                         NekDouble *coeffs, 
                                         NekDouble *phys);

            STD_REGIONS_EXPORT StdPyrExp(const StdPyrExp &T);

            STD_REGIONS_EXPORT ~StdPyrExp();

            //---------------------------------------
            // Integration/public 3D methods
            //---------------------------------------
            STD_REGIONS_EXPORT void TripleTensorProduct(
                const Array<OneD, const NekDouble>& fx, 
                const Array<OneD, const NekDouble>& gy, 
                const Array<OneD, const NekDouble>& hz, 
                const Array<OneD, const NekDouble>& inarray, 
                Array<OneD, NekDouble> & outarray);

            STD_REGIONS_EXPORT NekDouble TripleInnerProduct(
                const Array<OneD, const NekDouble>& fxyz,
                const Array<OneD, const NekDouble>& wx,
                const Array<OneD, const NekDouble>& wy,
                const Array<OneD, const NekDouble>& wz);

            STD_REGIONS_EXPORT NekDouble Integral3D(
                const Array<OneD, const NekDouble>& inarray,
                const Array<OneD, const NekDouble>& wx,
                const Array<OneD, const NekDouble>& wy,
                const Array<OneD, const NekDouble>& wz);

            STD_REGIONS_EXPORT void WriteCoeffsToFile(std::ofstream &outfile);

        protected:
            //---------------------------------------
            // Differentiation/integration Methods
            //---------------------------------------

            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2);

            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble>& outarray);

            STD_REGIONS_EXPORT virtual void v_StdPhysDeriv(
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &out_d0,
                Array<OneD, NekDouble> &out_d1,
                Array<OneD, NekDouble> &out_d2);

            STD_REGIONS_EXPORT virtual NekDouble v_Integral(
                const Array<OneD, const NekDouble>& inarray);


            //---------------------------------------
            // Transforms
            //---------------------------------------

            STD_REGIONS_EXPORT virtual void v_BwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_FwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);


            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            /*
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,       NekDouble> &outarray); 
            */

            //---------------------------------------
            // Evaluation functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble>& xi);

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble>& xi,
                const Array<OneD, const NekDouble>& physvals);

            STD_REGIONS_EXPORT virtual void v_GetCoords(
                Array<OneD, NekDouble> & xi_x, 
                Array<OneD, NekDouble> & xi_y, 
                Array<OneD, NekDouble> & xi_z);

            STD_REGIONS_EXPORT virtual void v_FillMode(
                const int mode, 
                Array<OneD, NekDouble> &outarray);  

            //---------------------------------------
            // Helper functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual int v_GetNverts() const;
            STD_REGIONS_EXPORT virtual int v_GetNedges() const;
            STD_REGIONS_EXPORT virtual int v_GetNfaces() const;
            STD_REGIONS_EXPORT virtual LibUtilities::ShapeType v_DetShapeType() const;
            STD_REGIONS_EXPORT virtual int v_NumBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual int v_GetEdgeNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_GetFaceNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_GetFaceIntNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_CalcNumberOfCoefficients(
                const std::vector<unsigned int> &nummodes, 
                int &modes_offset);
            STD_REGIONS_EXPORT virtual LibUtilities::BasisType v_GetEdgeBasisType(
                const int i) const;
            /*
            STD_REGIONS_EXPORT virtual void v_WriteToFile(
                std::ofstream &outfile,
                OutputFormat format,
                const bool dumpVar = true,
                std::string var = "v");
            */

            //---------------------------------------
            // Mappings
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_GetFaceToElementMap(
                const int                  fid,
                const Orientation      faceOrient,
                Array<OneD, unsigned int> &maparray,
                Array<OneD,          int> &signarray,
                int                        nummodesA=-1,
                int                        nummodesB=-1);
            STD_REGIONS_EXPORT virtual int  v_GetVertexMap(int localVertexId);
            /*
            STD_REGIONS_EXPORT virtual void v_GetEdgeInteriorMap(
                const int eid,
                const Orientation edgeOrient,
                Array<OneD, unsigned int> &maparray,
                Array<OneD, int> &signarray);
            STD_REGIONS_EXPORT virtual void v_GetFaceInteriorMap(
                const int fid,
                const Orientation faceOrient,
                Array<OneD, unsigned int> &maparray,
                Array<OneD, int>& signarray);
            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(
                Array<OneD, unsigned int> &outarray);
            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(
                Array<OneD, unsigned int>& outarray);
            */

            //---------------------------------------
            // Wrapper functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey);
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey);

        private:
            //---------------------------------------
            // Private helper functions
            //---------------------------------------
            STD_REGIONS_EXPORT int GetMode(int I, int J, int K);
            STD_REGIONS_EXPORT int GetTetMode(int I, int J, int K);
            STD_REGIONS_EXPORT void MultiplyByQuadratureMetric(
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray);


            /*
            // int GetBasisNumModes   (const int dir)

            int GetFaceNcoeffs(const int i) const
            {
            }

            const int GetEdgeNcoeffs(const int i) const
            {

            }

            
            const LibUtilities::BasisType GetEdgeBasisType(const int i) const
            {
                ASSERTL2((i >= 0)&&(i <= 7),"edge id is out of range");

                if((i == 0)||(i == 2))
                {
                    return  GetBasisType(0);
                }
                else if((i == 1)||(i == 3))
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
            */

            /*
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

            STD_REGIONS_EXPORT void FillMode(const int mode, Array<OneD, NekDouble> &outarray);

            */
            /** \brief Backward tranform for triangular elements
             *
             *  \b Note: That 'r' (base[2]) runs fastest in this element
             */
            /*
            STD_REGIONS_EXPORT void BwdTrans(const Array<OneD, const NekDouble>& inarray, Array<OneD, NekDouble> &outarray);   
            STD_REGIONS_EXPORT void FwdTrans(const Array<OneD, const NekDouble>& inarray,  Array<OneD, NekDouble> &outarray);
            */
            /** \brief Single Point Evaluation */
            /*
            STD_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& xi);
            STD_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& xi, const Array<OneD, const NekDouble> & physvals);
       
            STD_REGIONS_EXPORT void GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z);
            STD_REGIONS_EXPORT void WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");
            STD_REGIONS_EXPORT void WriteCoeffsToFile(std::ofstream &outfile);
                    
            DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
            {
                return StdExpansion::CreateGeneralMatrix(mkey);
            }               
                                           
        
        protected:
            */

            /** 
                \brief Calculate the inner product of inarray with respect to
                the basis B=base0*base1*base2 and put into outarray:
              
                \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
                \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2}
                \psi_{p}^{a} (\bar \eta_{1i}) \psi_{q}^{a} (\eta_{2j}) \psi_{pqr}^{c} (\eta_{3k})
                w_i w_j w_k u(\bar \eta_{1,i} \eta_{2,j} \eta_{3,k})	     
                J_{i,j,k}\\ & = & \sum_{i=0}^{nq_0} \psi_p^a(\bar \eta_{1,i})
                \sum_{j=0}^{nq_1} \psi_{q}^a(\eta_{2,j}) \sum_{k=0}^{nq_2} \psi_{pqr}^c u(\bar \eta_{1i},\eta_{2j},\eta_{3k})
                J_{i,j,k} \end{array} \f$ \n
            
                where
            
                \f$\phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\bar \eta_1) \psi_{q}^a (\eta_2) \psi_{pqr}^c (\eta_3) \f$ \n
            
                which can be implemented as \n
                \f$f_{pqr} (\xi_{3k}) = \sum_{k=0}^{nq_3} \psi_{pqr}^c u(\bar \eta_{1i},\eta_{2j},\eta_{3k})
                J_{i,j,k} = {\bf B_3 U}  \f$ \n
                \f$ g_{pq} (\xi_{3k}) = \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{pqr} (\xi_{3k})  = {\bf B_2 F}  \f$ \n
                \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0} \psi_{p}^a (\xi_{3k}) g_{pq} (\xi_{3k})  = {\bf B_1 G} \f$

            **/
            // Interior pyramid implementation based on Spen's book page 108. 113. and 609.
            /*
            STD_REGIONS_EXPORT void IProductWRTBase(const Array<OneD, const NekDouble>& bx, 
                                 const Array<OneD, const NekDouble>& by, 
                                 const Array<OneD, const NekDouble>& bz, 
                                 const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> & outarray);      
            */
            /*
        private:
        
            virtual int v_GetNverts() const
            {
                return 5;
            }
        
            virtual int v_GetNedges() const
            {
                return 8;
            }
        
            virtual int v_GetNfaces() const
            {
                return 5;
            }

            virtual LibUtilities::ShapeType v_DetShapeType() const
            {
                return DetShapeType();
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

            virtual void v_GetEdgeInteriorMap(const int eid, const Orientation edgeOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray)
            {
                GetEdgeInteriorMap(eid,edgeOrient,maparray,signarray);
            } 
                      
            virtual void v_GetFaceToElementMap(const int fid, const Orientation faceOrient,
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

            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& Lcoords)
            {
                return PhysEvaluate(Lcoords);
            }

            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& Lcoords,  const Array<OneD, const NekDouble> & physvals)
            {
                return PhysEvaluate(Lcoords, physvals);
            }
        
            virtual int v_GetEdgeNcoeffs(const int i) const
            {
                return GetEdgeNcoeffs(i);
            }
                
            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v")
            {
                WriteToFile(outfile,format,dumpVar,var);
            }

            virtual void v_WriteCoeffsToFile(std::ofstream &outfile)
            {
                WriteCoeffsToFile(outfile);
            }
            */
        };    
        typedef boost::shared_ptr<StdPyrExp> StdPyrExpSharedPtr;
    } //end of namespace
} //end of namespace

#endif //STDPYREXP_H
