///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTetExp.h
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
// Description: Header field for nodal tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDNODALTETEXP_H
#define STDNODALTETEXP_H

#include <StdRegions/StdExpansion.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdMatrix.h>

#include <StdRegions/NodalBasisManager.h>

namespace Nektar
{
  namespace StdRegions
  {

  class StdNodalTetExp: public StdTetExp
  {

  public:

    //Constructors
    StdNodalTetExp(const BasisKey &Ba, const BasisKey &Bb, const Basis *Bc);

    StdNodalTetExp(const BasisKey &Ba, const BasisKey &Bb, const Basis *Bc,
           double *coeffs, double *phys);

    /// Copy Constructor
    StdNodalTetExp(StdNodalTetExp &T);

    /// Destructor
    ~StdNodalTetExp();

    ShapeType DetShapeType()
    {
      return eTetrahedron;
    }

    //////////////////////////////
    /// Integration Methods
    //////////////////////////////

    double Integral(const double *inarray);
    void IProductWRTBase(const double * inarray, double * outarray);
    void FillMode(const int mode, double *outarray);

    void GenMassMatrix(double * outarray);
    void GenLapMatrix(double * outarray);

    StdMatContainer * GetMassMatrix();
    StdMatContainer * GetLapMatrix();

    //-----------------------------
    // Differentiation Methods
    //-----------------------------

    //-----------------------------
    // Evaluations Methods
    //-----------------------------

    void BwdTrans(double * outarray);

    void FwdTrans(const double * inarray);

    double Evaluate(double x, double y, double z);

  protected:

    static StdMatrix s_elmtmats;

    inline void IProductWRTBase(const double *base0, const double *base1,
                const double *inarray, double *outarray);

  private:
    virtual ShapeType V_DetShapeType()
    {
		return DetShapeType();
    }

    virtual double V_Integral(const double *inarray )
    {
		return Integral(inarray);
    }

    virtual void V_IProduct_WRT_B(const double * inarray, double * outarray)
    {
		IProductWRTBase(inarray,outarray);
    }

    virtual void V_FillMode(const int mode, double *outarray)
    {
		FillMode(mode,outarray);
    }

    virtual void V_GenMassMatrix(double * outarray)
    {
		GenMassMatrix(outarray);
    }

    virtual void V_GenLapMatrix(double * outarray)
    {
		GenLapMatrix(outarray);
    }

    virtual StdMatContainer * V_GetMassMatrix()
    {
		return GetMassMatrix();
    }

    virtual StdMatContainer * V_GetLapMatrix()
    {
		return GetLapMatrix();
    }

    virtual void V_BwdTrans(double * outarray)
    {
		BwdTrans(outarray);
    }

    virtual void V_FwdTrans(const double * inarray)
    {
		FwdTrans(inarray);
    }

    virtual double V_Evaluate(double x, double y, double z)
    {
		return Evaluate(x,y,z);
    }

    virtual double V_Linf(const double *sol)
    {
		return Linf(sol);
    }

    virtual double V_Linf()
    {
		return Linf();
    }

    virtual double V_L2(const double *sol)
    {
		return L2(sol);
    }

    virtual double V_L2()
    {
		return L2();
    }
  };

  } //end of namespace
} //end of namespace
#endif //STDNODALTETEXP_H

/**
 * $Log: StdNodalTetExp.h,v $
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.6  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.5  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 **/

