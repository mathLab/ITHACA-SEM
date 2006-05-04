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

#include <StdRegions/StdExpansion2D.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdMatrix.h>

#include <StdRegions/NodalBasisManager.h>

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/LocalRegionsDeclarations.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        class StdNodalTriExp: public StdTriExp
        {

        public:

            //Constructors
            StdNodalTriExp(const BasisKey &Ba, const BasisKey &Bb,
                NodalBasisType Ntype);

            //Copy Constructor
            StdNodalTriExp(StdNodalTriExp &T);

            //Destructor
            ~StdNodalTriExp();

            ShapeType DetShapeType()
            {
                return eTriangle;
            }

            /////////////////////////////
            /// Nodal basis specific routines
            ///////////////////////////

            void GenNBasisTransMatrix(double * outarray);


            StdMatContainer * GetNBasisTransMatrix();

            void NodalToModal();
            void ModalToNodal();

            int  GetNodalPoints(const double * &x, const double * &y)
            {
                const double *z;
                return NBasisManagerSingleton::Instance().GetNodePoints(m_nbtype,
                    m_base[0]->GetBasisOrder(),x,y,z);
            }

            //////////////////////////////
            /// Integration Methods
            //////////////////////////////

            void IProductWRTBase(const double * inarray, double * outarray);

            void FillMode(const int mode, double *outarray);

            void GenMassMatrix(double * outarray);

            void GenLapMatrix(double * outarray);

            StdMatContainer * GetMassMatrix();

            StdMatContainer * GetLapMatrix();

            //-----------------------------
            // Differentiation Methods
            //-----------------------------

            void Deriv(const double *inarray, double * outarray_d0,
                double *outarray_d1);

            void Deriv(double * outarray_d0, double *outarray_d1);

            //-----------------------------
            // Evaluations Methods
            //-----------------------------

            void BwdTrans(double * outarray);

            void FwdTrans(const double * inarray);

            double Evaluate(const double * coords);


        protected:

            static StdMatrix s_elmtmats;

            /// All Expansions share the same NodalBasisManager
            typedef Loki::SingletonHolder<NodalBasisManager> NBasisManagerSingleton;

            inline void IProductWRTBase(const double *base0, const double *base1,
                const double *inarray, double *outarray);

        private:

            NodalBasisType m_nbtype;

            virtual ShapeType v_DetShapeType()
            {
                DetShapeType();
            }

            virtual void v_GenNBasisTransMatrix(double * outarray)
            {
                GenNBasisTransMatrix(outarray);
            }

            virtual StdMatContainer * v_NBasisTransMatrix()
            {
                return GetNBasisTransMatrix();
            }


            //////////////////////////////
            /// Integration Methods
            //////////////////////////////

            virtual double v_Integral(const double *inarray )
            {
                return Integral(inarray);
            }

            virtual void v_IProductWRTBase(const double * inarray, double * outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            virtual void v_FillMode(const int mode, double *outarray)
            {
                FillMode(mode,outarray);
            }

            virtual void v_GenMassMatrix(double * outarray)
            {
                GenMassMatrix(outarray);
            }

            virtual void v_GenLapMatrix(double * outarray)
            {
                GenLapMatrix(outarray);
            }

            virtual StdMatContainer * v_GetMassMatrix()
            {
                return GetMassMatrix();
            }

            virtual StdMatContainer * v_GetLapMatrix()
            {
                return GetLapMatrix();
            }

            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            virtual void v_Deriv(const double *inarray, double * outarray_d0,
                double *outarray_d1)
            {
                Deriv(inarray,outarray_d0,outarray_d1);
            }

            virtual void v_Deriv(double * outarray_d0, double *outarray_d1)
            {
                Deriv(this->m_phys,outarray_d0,outarray_d1);
            }

            //-----------------------------
            // Evaluations Methods
            //-----------------------------

            virtual void v_BwdTrans(double * outarray)
            {
                return BwdTrans(outarray);
            }

            virtual void v_FwdTrans(const double * inarray)
            {
                return FwdTrans(inarray);
            }

            virtual double v_Evaluate(const double * coords)
            {
                return Evaluate(coords);
            }

            virtual int v_GetNodalPoints(const double * &x, const double* &y)
            {
                return GetNodalPoints(x,y);
            }

        };

    } //end of namespace
} //end of namespace

#endif //STDNODALTRIEXP_H

/**
* $Log: StdNodalTriExp.h,v $
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

