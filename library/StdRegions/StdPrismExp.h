///////////////////////////////////////////////////////////////////////////////
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


namespace Nektar
{
    namespace StdRegions
    {
    
    class StdPrismExp: public StdExpansion3D
    {
        
    public:
        
        StdPrismExp();
        
        /** \brief Constructor using BasisKey class for quadrature
         *    points and order definition 
         */
        StdPrismExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc);
        
        /** \brief Constructor using BasisKey class for quadrature
         *  points and order definition where m_coeffs and m_phys are all
         *    set. 
         */
        StdPrismExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc,
            double *coeffs, double *phys);
        
        /** \brief Copy Constructor */
        StdPrismExp(const StdPrismExp &T);
        
        /** \brief Destructor */
        ~StdPrismExp();
        
        /** \brief Return Shape of region, using  ShapeType enum list.
         *  i.e. Prism 
         */
        ShapeType DetShapeType() const
        {
        return ePrism;
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

        NekDouble Integral3D(const Array<OneD, const NekDouble>& inarray, 
                    const Array<OneD, const NekDouble>& wx,
                    const Array<OneD, const NekDouble>& wy, 
                    const Array<OneD, const NekDouble>& wz);
        NekDouble Integral(const Array<OneD, const NekDouble>& inarray);
        void IProductWRTBase(const Array<OneD, const NekDouble>& inarray, Array<OneD, NekDouble> &outarray);
        void IProductWRTBase(  const Array<OneD, const NekDouble>& bx, 
                               const Array<OneD, const NekDouble>& by, 
                               const Array<OneD, const NekDouble>& bz, 
                               const Array<OneD, const NekDouble>& inarray, 
                               Array<OneD, NekDouble> & outarray );
        void PhysDeriv( Array<OneD, NekDouble> &out_d0,
                        Array<OneD, NekDouble> &out_d1,
                        Array<OneD, NekDouble> &out_d2);
        void PhysDeriv(const Array<OneD, const NekDouble>& u_physical, 
                       Array<OneD, NekDouble> &out_dxi1, 
                       Array<OneD, NekDouble> &out_dxi2,
                       Array<OneD, NekDouble> &out_dxi3 );
                       
         void BwdTrans(const Array<OneD, const NekDouble>& inarray, 
                       Array<OneD, NekDouble> &outarray);
         void FwdTrans( const Array<OneD, const NekDouble>& inarray,  Array<OneD, NekDouble> &outarray);
         NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& xi);
         void GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z);
         void FillMode(const int mode, Array<OneD, NekDouble> &outarray);        
         void GenLapMatrix(double * outarray);
                       
         DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
         {
             return StdExpansion::CreateGeneralMatrix(mkey);
         }

        

        
        
    protected:
        

        
        
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

        virtual ShapeType v_DetShapeType() const
        {
            return DetShapeType();
        }
        
        virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey) 
        {
            return GenMatrix(mkey);
        }
        
        virtual void v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            return FillMode(mode, outarray);
        }

        virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray )
        {
            return Integral(inarray);
        }
            
        virtual void v_GetCoords(
            Array<OneD, NekDouble> &coords_x,
            Array<OneD, NekDouble> &coords_y,
            Array<OneD, NekDouble> &coords_z)
        {
            GetCoords(coords_x, coords_y, coords_z);
        }
        
        virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray, outarray);
        }

        /** \brief Virtual call to GenMassMatrix */

        virtual void v_PhysDeriv( Array<OneD, NekDouble> &out_d0,
                                Array<OneD, NekDouble> &out_d1,
                                Array<OneD, NekDouble> &out_d2)
        {
                PhysDeriv(out_d0, out_d1, out_d2);
        }
        virtual void v_StdPhysDeriv( Array<OneD, NekDouble> &out_d0,
                                Array<OneD, NekDouble> &out_d1,
                                Array<OneD, NekDouble> &out_d2)
        {
                                
                PhysDeriv(out_d0, out_d1, out_d2);                
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
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
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
        
        virtual int v_GetEdgeNcoeffs(const int i) const
        {
            return GetEdgeNcoeffs(i);
        }

        virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const
        {
            return GetEdgeBasisType(i);
        }
        
        virtual void v_GenMassMatrix(Array<OneD, NekDouble> & outarray)
        {
                std::cout << "Implement me" << std::endl;
                return;
        } 
        
        virtual void v_GenLapMatrix (Array<OneD, NekDouble> & outarray)
        {
            std::cout << "Implement me" << std::endl;
            return;
        }
        
        virtual DNekMatSharedPtr v_GetMassMatrix() 
        {
            std::cout << "Implement me" << std::endl;
            int foo = 0;
            return DNekMatSharedPtr();
        } 
        
        virtual DNekMatSharedPtr v_GetLapMatrix()
        {
            std::cout << "Implement me" << std::endl;
            int foo = 0;
            return DNekMatSharedPtr();
        }  

        
        
    };
        typedef boost::shared_ptr<StdPrismExp> StdPrismExpSharedPtr;
    
    } //end of namespace
} //end of namespace

#endif //STDPRISMEXP_H

/**
 * $Log: StdPrismExp.h,v $
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

