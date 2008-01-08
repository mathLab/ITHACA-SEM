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


namespace Nektar
{
    namespace StdRegions
    {
    
    class StdPyrExp: public StdExpansion3D
    {
        
    public:
    
        StdPyrExp();
        
        /** \brief Constructor using BasisKey class for quadrature
         *  points and order definition 
         */
        StdPyrExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc);
        
        /** \brief Constructor using BasisKey class for quadrature
         *    points and order definition where m_coeffs and m_phys are all
         *  set
         */
        StdPyrExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc,
              double *coeffs, double *phys);
        
        /** \brief Copy Constructor */
        StdPyrExp(const StdPyrExp &T);
        
        /** \brief Destructor */
        ~StdPyrExp();
        
        /** \brief Return Shape of region, using  ShapeType enum list. 
         *  i.e. Pyramid
         */
        ShapeType DetShapeType() const
        {
        return ePyramid;
        };
        
        NekDouble Integral3D(const ConstArray<OneD, NekDouble>& inarray, 
                             const ConstArray<OneD, NekDouble>& wx,
                             const ConstArray<OneD, NekDouble>& wy, 
                             const ConstArray<OneD, NekDouble>& wz);
       NekDouble Integral(const ConstArray<OneD, NekDouble>& inarray);
       void IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> &outarray);
       void IProductWRTBase(const ConstArray<OneD, NekDouble>& bx, 
                            const ConstArray<OneD, NekDouble>& by, 
                            const ConstArray<OneD, NekDouble>& bz, 
                            const ConstArray<OneD, NekDouble>& inarray, 
                            Array<OneD, NekDouble> & outarray );
       void PhysDeriv( Array<OneD, NekDouble> &out_d0,
                       Array<OneD, NekDouble> &out_d1,
                       Array<OneD, NekDouble> &out_d2);
       void PhysDeriv(const ConstArray<OneD, NekDouble>& u_physical, 
                            Array<OneD, NekDouble> &out_dxi1, 
                            Array<OneD, NekDouble> &out_dxi2,
                            Array<OneD, NekDouble> &out_dxi3 );  
       void FillMode(const int mode, Array<OneD, NekDouble> &outarray);
       void BwdTrans(const ConstArray<OneD, NekDouble>& inarray, Array<OneD, NekDouble> &outarray);   
       void FwdTrans(const ConstArray<OneD, NekDouble>& inarray,  Array<OneD, NekDouble> &outarray); 
       NekDouble PhysEvaluate(const ConstArray<OneD, NekDouble>& xi);                 
       void GetCoords( Array<OneD, NekDouble> & xi_x, Array<OneD, NekDouble> & xi_y, Array<OneD, NekDouble> & xi_z);
       void GenLapMatrix(double * outarray);
                    
       DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
       {
           return StdExpansion::CreateGeneralMatrix(mkey);
       }               
        
        
                     
        
        
    protected:
        
        
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

        virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble>& inarray )
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
        
        virtual void v_IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray,
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

        virtual void v_PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                                Array<OneD, NekDouble> &out_d0,
                                Array<OneD, NekDouble> &out_d1,
                                Array<OneD, NekDouble> &out_d2)
        {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }                                  
        virtual void v_StdPhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
                                Array<OneD, NekDouble> &out_d0,
                                Array<OneD, NekDouble> &out_d1,
                                Array<OneD, NekDouble> &out_d2)
        {
                PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }

        
        virtual void v_BwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
        {
            BwdTrans(inarray, outarray);
        }

        virtual void v_FwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
        {
            FwdTrans(inarray, outarray);
        }

        virtual NekDouble v_PhysEvaluate(const ConstArray<OneD, NekDouble>& Lcoords)
        {
            return PhysEvaluate(Lcoords);
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
        } 
        
        virtual DNekMatSharedPtr v_GetLapMatrix()
        {
            std::cout << "Implement me" << std::endl;
            int foo = 0;
        } 
        
        
    };
    
    } //end of namespace
} //end of namespace

#endif //STDPYREXP_H

/**
 * $Log: StdPyrExp.h,v $
 * Revision 1.7  2008/01/03 12:33:14  ehan
 * Fixed errors from StdMatrix to StdMatrixKey.
 *
 * Revision 1.6  2008/01/03 10:40:52  ehan
 * Added basis, differentiation, backward transform, forward transform, interpolation, integration, and physEval functions.
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

