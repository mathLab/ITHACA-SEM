///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Helmholtz operator implementations
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MatrixFreeOps/Operator.hpp>
#include <MatrixFreeOps/Util.hpp>

#include <Collections/Operator.h>
#include <Collections/Collection.h>
#include <Collections/IProduct.h>

using namespace std;

namespace Nektar {
namespace Collections {

using LibUtilities::eSegment;
using LibUtilities::eQuadrilateral;
using LibUtilities::eTriangle;
using LibUtilities::eHexahedron;
using LibUtilities::eTetrahedron;
using LibUtilities::ePrism;
using LibUtilities::ePyramid;

/**
 * @brief Helmholtz operator using LocalRegions implementation.
 */
class Helmholtz_NoCollection : public Operator
{
    public:
        OPERATOR_CREATE(Helmholtz_NoCollection)

        virtual ~Helmholtz_NoCollection()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &entry0,
                      Array<OneD, NekDouble> &entry1,
                      Array<OneD, NekDouble> &entry2,
                      Array<OneD, NekDouble> &entry3,
                      Array<OneD, NekDouble> &wsp,
                const StdRegions::ConstFactorMap  &factors) final
        {
            boost::ignore_unused(entry2,entry3,wsp);

            unsigned int nmodes = m_expList[0]->GetNcoeffs();
            Array<OneD, NekDouble> tmp;

            for(int n = 0; n < m_numElmt; ++n)
            {
                StdRegions::StdMatrixKey mkey(StdRegions::eHelmholtz,
                                              (m_expList)[n]->DetShapeType(),
                                              *(m_expList)[n], factors); 
                m_expList[n]->GeneralMatrixOp(entry0 + n *nmodes,
                                              tmp = entry1 + n * nmodes,
                                              mkey);
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
        }

    protected:
        int                                         m_dim;
        int                                         m_coordim;
        vector<StdRegions::StdExpansionSharedPtr>   m_expList;

    private:
        Helmholtz_NoCollection(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData)
        {
            m_expList = pCollExp;
            m_dim     = pCollExp[0]->GetNumBases();
            m_coordim = pCollExp[0]->GetCoordim();
        }
};

/// Factory initialisation for the Helmholtz_NoCollection operators
OperatorKey Helmholtz_NoCollection::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eNoCollection,true),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eNoCollection,true),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eNoCollection,true),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    eHelmholtz, eNoCollection,false),
        Helmholtz_NoCollection::create,
        "Helmholtz_NoCollection_Hex")
};


/**
 * @brief Helmholtz operator using LocalRegions implementation.
 */
class Helmholtz_IterPerExp : public Operator
{
    public:
        OPERATOR_CREATE(Helmholtz_IterPerExp)

        virtual ~Helmholtz_IterPerExp()
        {
        }

        virtual void operator()(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD, NekDouble> &output,
                      Array<OneD, NekDouble> &output1,
                      Array<OneD, NekDouble> &output2,
                      Array<OneD, NekDouble> &wsp,
                const StdRegions::ConstFactorMap   &factors) final
        {
            boost::ignore_unused(output1,output2);

            const int nCoeffs = m_stdExp->GetNcoeffs();
            const int nPhys   = m_stdExp->GetTotPoints();

            ASSERTL1(input.size() >= m_numElmt*nCoeffs,
                     "input array size is insufficient");
            ASSERTL1(output.size() >= m_numElmt*nCoeffs,
                     "output array size is insufficient");
            
            Array<OneD, NekDouble> tmpphys, t1; 
            Array<OneD, Array<OneD, NekDouble> > dtmp(3);
            Array<OneD, Array<OneD, NekDouble> > tmp(3);

            auto x = factors.find(StdRegions::eFactorLambda);
            ASSERTL1(x != factors.end(),
                     "Constant factor not defined: "
                     + std::string(StdRegions::ConstFactorTypeMap
                                   [StdRegions::eFactorLambda]));
            NekDouble lambda = x->second; 
            
            tmpphys = wsp; 
            for(int i = 1; i < m_coordim+1; ++i)
            {
                dtmp[i-1] = wsp + i*nPhys;
                tmp [i-1] = wsp + (i+m_coordim)*nPhys; 
            }

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->BwdTrans(input + i*nCoeffs, tmpphys);

                // local derivative 
                m_stdExp->PhysDeriv(tmpphys, dtmp[0], dtmp[1], dtmp[2]);

                // determine mass matrix term 
                Vmath::Vmul(nPhys,m_jac+i*nPhys,1,tmpphys,1,tmpphys,1);
                m_stdExp->IProductWRTBase(tmpphys,t1 = output + i*nCoeffs); 
                Vmath::Smul(nCoeffs,lambda,output + i*nCoeffs,1,
                            t1 = output+i*nCoeffs,1);
                
                // calculate full derivative
                for(int j = 0; j < m_coordim; ++j)
                {
                    Vmath::Vmul(nPhys,m_derivFac[j*m_dim].origin() + i*nPhys,1,
                                &dtmp[0][0],1,&tmp[j][0],1);

                    for(int k = 1; k < m_dim; ++k)
                    {
                        Vmath::Vvtvp (nPhys, m_derivFac[j*m_dim+k].origin()
                                      + i*nPhys, 1, &dtmp[k][0], 1,
                                      &tmp[j][0],   1,  &tmp[j][0],   1);
                    }
                }

                // calculate dx/dxi tmp[0] + dy/dxi tmp[2] + dz/dxi tmp[3]
                for(int j = 0; j < m_dim; ++j)
                {
                    Vmath::Vmul (nPhys,m_derivFac[j].origin() + i*nPhys,1,
                                 &tmp[0][0], 1, &dtmp[j][0],1);

                    for(int k = 1; k < m_coordim; ++k)
                    {
                        Vmath::Vvtvp (nPhys, m_derivFac[j +k*m_dim].origin()
                                      + i*nPhys, 1, &tmp[k][0], 1,
                                      &dtmp[j][0], 1, &dtmp[j][0], 1);
                    }
                }

                // calculate Iproduct WRT Std Deriv
                for(int j = 0; j < m_dim; ++j)
                {
                    // multiply by Jacobian
                    Vmath::Vmul(nPhys,m_jac+i*nPhys,1,dtmp[j],1,dtmp[j],1);
                    m_stdExp->IProductWRTDerivBase(j,dtmp[j],tmp[0]);
                    Vmath::Vadd(nCoeffs,tmp[0],1,output+i*nCoeffs,1,
                                t1 = output+i*nCoeffs,1);
                }
            }
        }

        virtual void operator()(
                      int                           dir,
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output,
                      Array<OneD,       NekDouble> &wsp)
        {
            boost::ignore_unused(dir, input, output, wsp);
            NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
        }

    protected:
        Array<TwoD, const NekDouble>    m_derivFac;
        Array<OneD, const NekDouble>    m_jac;
        int                             m_dim;
        int                             m_coordim;
    
    private:
        Helmholtz_IterPerExp(
                vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                CoalescedGeomDataSharedPtr                pGeomData)
            : Operator(pCollExp, pGeomData)
        {
            LibUtilities::PointsKeyVector PtsKey = m_stdExp->GetPointsKeys();
            m_dim      = PtsKey.size();
            m_coordim  = pCollExp[0]->GetCoordim();
            int nqtot  = m_stdExp->GetTotPoints();

            m_derivFac = pGeomData->GetDerivFactors(pCollExp);
            m_jac      = pGeomData->GetJac(pCollExp);
            m_wspSize = (2*m_coordim+1)*nqtot;
        }
};

/// Factory initialisation for the Helmholtz_IterPerExp operators
OperatorKey Helmholtz_IterPerExp::m_typeArr[] = {
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eSegment,       eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Seg"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle,      eHelmholtz, eIterPerExp,true),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_NodalTri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Tet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron,   eHelmholtz, eIterPerExp,true),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_NodalTet"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid,       eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism,         eHelmholtz, eIterPerExp,true),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_NodalPrism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron,    eHelmholtz, eIterPerExp,false),
        Helmholtz_IterPerExp::create,
        "Helmholtz_IterPerExp_Hex")
};
    
/**
 * @brief Helmholtz operator using matrix free operators.
 */
class Helmholtz_MatrixFree : public Operator
{
public:
    OPERATOR_CREATE(Helmholtz_MatrixFree)
    
    virtual ~Helmholtz_MatrixFree()
    {
    }
    
    virtual void operator()(
            const Array<OneD, const NekDouble> &input,
                  Array<OneD,       NekDouble> &output0,
                  Array<OneD,       NekDouble> &output1,
                  Array<OneD,       NekDouble> &output2,
                  Array<OneD,       NekDouble> &wsp,
                const StdRegions::ConstFactorMap   &factors) final
    {
        boost::ignore_unused(output1,output2,wsp);

        // Set lambda for this call
        auto x = factors.find(StdRegions::eFactorLambda);
        ASSERTL1(x != factors.end(),
                 "Constant factor not defined: "
                 + std::string(StdRegions::ConstFactorTypeMap[StdRegions::eFactorLambda]));
        m_oper->SetLambda(x->second);
        
        if (m_isPadded)
        {
            // copy into padded vector
            Vmath::Vcopy(m_nmtot, input, 1, m_input, 1);
            (*m_oper)(m_input, m_output);
            Vmath::Vcopy(m_nmtot, m_output, 1, output0, 1);
        }
        else
        {
            (*m_oper)(input, output0);
        }

    }
    
    virtual void operator()(int                           dir,
                            const Array<OneD, const NekDouble> &input,
                            Array<OneD,       NekDouble> &output,
                            Array<OneD,       NekDouble> &wsp)
    {
        boost::ignore_unused(dir,input,output,wsp);
        NEKERROR(ErrorUtil::efatal, "Not valid for this operator.");
    }
    
private:
    std::shared_ptr<MatrixFree::Helmholtz> m_oper;
    /// flag for padding
    bool m_isPadded{false};
    /// padded or unpadded input/output vectors
    Array<OneD, NekDouble> m_input;
    Array<OneD, NekDouble> m_output;
    unsigned int m_nmtot; 
    
    Helmholtz_MatrixFree(vector<StdRegions::StdExpansionSharedPtr> pCollExp,
                         CoalescedGeomDataSharedPtr                pGeomData)
        : Operator(pCollExp, pGeomData)
    {
        const auto nmElmt = pCollExp[0]->GetStdExp()->GetNcoeffs();
        const auto nqElmt = pCollExp[0]->GetStdExp()->GetTotPoints();

        int coordim = pCollExp[0]->GetCoordim();
        
        // Padding if needed
        using vec_t = tinysimd::simd<NekDouble>;
        const auto nElmtNoPad = pCollExp.size();
        auto nElmtPad = nElmtNoPad;
        
        m_nmtot = nElmtNoPad*nmElmt; 

        if (nElmtNoPad % vec_t::width != 0)
        {
            m_isPadded = true;
            nElmtPad = nElmtNoPad + vec_t::width -
                (nElmtNoPad % vec_t::width);
            m_input  = Array<OneD, NekDouble>{nmElmt * nElmtPad, 0.0};
            m_output = Array<OneD, NekDouble>{nmElmt * nElmtPad, 0.0};
        }
        else
        {
            m_output = Array<OneD, NekDouble>{nmElmt*nElmtNoPad, 0.0};
        }
       
        // Check if deformed
        bool deformed{pGeomData->IsDeformed(pCollExp)};

        // Size of jacobian
        int jacSizeNoPad{nElmtNoPad};
        int jacSizePad{nElmtPad};

        if (deformed)
        {
            jacSizeNoPad = nElmtNoPad * nqElmt;
            jacSizePad = nElmtPad * nqElmt;
        }
        
        // Get Jacobian
        Array<OneD, NekDouble> jac{jacSizePad, 0.0};
        Vmath::Vcopy(jacSizeNoPad, pGeomData->GetJac(pCollExp), 1, jac, 1);
        
        // Get derivative factors
        const auto dim = pCollExp[0]->GetStdExp()->GetShapeDimension();
        Array<TwoD, NekDouble> df(dim * coordim, jacSizePad, 0.0);

        if (deformed)
        {
            for (unsigned int j = 0; j < dim * coordim; ++j)
            {
                Vmath::Vcopy(jacSizeNoPad,
                             &(pGeomData->GetDerivFactors(pCollExp))[j][0], 1,
                             &df[j][0], 1);
            }
        }
        else
        {
            for (unsigned int e = 0; e < nElmtNoPad; ++e)
            {
                for (unsigned int j = 0; j < dim * coordim; ++j)
                {
                    df[j][e] =
                        (pGeomData->GetDerivFactors(pCollExp))[j][e*nqElmt];
                }
            }
        }

        // Basis vector.
        std::vector<LibUtilities::BasisSharedPtr> basis(dim);
        for (auto i = 0; i < dim; ++i)
        {
            basis[i] = pCollExp[0]->GetBasis(i);
        }
        
        // Get shape type
        auto shapeType = pCollExp[0]->GetStdExp()->DetShapeType();
        
        // Generate operator string and create operator.
        std::string op_string = "Helmholtz";
        op_string += MatrixFree::GetOpstring(shapeType, deformed);
        auto oper = MatrixFree::GetOperatorFactory().
            CreateInstance(op_string, basis, nElmtPad);
        
        // Set Jacobian
        oper->SetJac(jac);
        
        // Store derivative factor
        oper->SetDF(df);
        
        m_oper = std::dynamic_pointer_cast<MatrixFree::Helmholtz>(oper);
        ASSERTL0(m_oper, "Failed to cast pointer.");
        
    }
};

/// Factory initialisation for the Helmholtz_MatrixFree operators
OperatorKey Helmholtz_MatrixFree::m_typeArr[] =
{
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eQuadrilateral, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Quad"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTriangle, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Tri"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eHexahedron, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Hex"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePrism, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Prism"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(ePyramid, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Pyr"),
    GetOperatorFactory().RegisterCreatorFunction(
        OperatorKey(eTetrahedron, eHelmholtz, eMatrixFree, false),
        Helmholtz_MatrixFree::create, "Helmholtz_MatrixFree_Tet"),
};

}
}
