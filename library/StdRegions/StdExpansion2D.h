///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion2D.h
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 2D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXP2D_H
#define STDEXP2D_H

#include <StdRegions/StdExpansion.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion2D: public StdExpansion
        {
        public:
            StdExpansion2D();
            StdExpansion2D(const BasisKey &Ba, const BasisKey &Bb, int numcoeffs,
                double *coeffs, double *phys, bool spaceowner);
            StdExpansion2D(const StdExpansion2D &T);
            ~StdExpansion2D();

            // Generic operations in different elements

            /** \brief Calculate the 2D derivative in the local collapsed
            coordinate at the physical points  */
            void  TensorDeriv(double *outarray_d0,double *outarray_d1);

            /** \brief Calculate the 2D derivative in the local collapsed
            coordinate at the physical points  */
            void  TensorDeriv(const double *inarray, double *outarray_d0,
                double *outarray_d1);

            void  Deriv(double *outarray_d1, double *outarray_d2)
            {
                v_Deriv(outarray_d1, outarray_d2);
            }

            void  StdDeriv(double *outarray_d1, double *outarray_d2)
            {
                v_StdDeriv(outarray_d1, outarray_d2);
            }

            void Deriv(const double *inarray, double *outarray_d1,
                double *outarray_d2)
            {
                v_Deriv(inarray, outarray_d1, outarray_d2);
            }

            void StdDeriv(const double *inarray, double *outarray_d1,
                double *outarray_d2)
            {
                v_StdDeriv(inarray, outarray_d1, outarray_d2);
            }

            /** \brief Evaluate a function at points coords which is assumed
            to be in local collapsed coordinate format. The function is
            assumed to be in physical space */
            double PhysEvaluate(const double *coords);

            double Integral(const double *inarray, const double *w0, const double* w1);

            int GetNodalPoints(const double * &x, const double * &y)
            {
                return v_GetNodalPoints(x,y);
            }

        protected:

        private:

            // Virtual Functions ----------------------------------------

	    virtual int v_GetNverts() = 0;
	    virtual int v_GetNedges() = 0;
	    virtual int v_GetNfaces()
	    {
		ASSERTL0(false,"This function is only valid for 3D expansions");
		return 0;
	    }

            virtual ShapeType v_DetShapeType()                = 0;

            virtual StdMatContainer *v_GetMassMatrix()        = 0;
            virtual StdMatContainer *v_GetLapMatrix()         = 0;

            virtual int v_GetNodalPoints(const double * &x, const double * &y)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
                return 0;
            }

            virtual void v_GenNBasisTransMatrix(double * outarray)
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
            }

            virtual StdMatContainer *v_GetNBasisTransMatrix()
            {
                ASSERTL0(false, "This function is only valid for nodal expansions");
                return NULL;
            }


            virtual void   v_BwdTrans (double *outarray)      = 0;
            virtual void   v_FwdTrans (const double *inarray) = 0;

            virtual double v_Integral(const double *inarray ) = 0;
            virtual double v_Evaluate(const double * coords) = 0;

            virtual void   v_Deriv(double *outarray_d1, double *outarray_d2)    = 0;
            virtual void   v_StdDeriv(double *outarray_d1, double *outarray_d2) = 0;
            virtual void   v_Deriv(const double *inarray, double *outarray_d1,
                double *outarray_d2) = 0;
            virtual void   v_StdDeriv(const double *inarray, double *outarray_d1,
                double *outarray_d2) = 0;


	    virtual int v_GetCoordim(void)
	    {
                return 2; 
	    }


        };

    } //end of namespace
} //end of namespace

#endif //STDEXP2D_H

/**
* $Log: StdExpansion2D.h,v $
* Revision 1.3  2006/07/02 17:16:18  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.2  2006/06/13 18:05:02  sherwin
* Modifications to make MultiRegions demo ProjectLoc2D execute properly.
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.19  2006/05/02 21:21:12  sherwin
* Corrected libraries to compile new version of spatialdomains and demo Graph1D
*
* Revision 1.18  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.17  2006/03/06 17:12:45  sherwin
*
* Updated to properly execute all current StdRegions Demos.
*
* Revision 1.16  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.15  2006/03/04 20:26:54  bnelson
* Added comments after #endif.
*
* Revision 1.14  2006/03/01 22:59:12  sherwin
*
* First working version of Project1D
*
* Revision 1.13  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.12  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/



