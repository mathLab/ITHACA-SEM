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
                const StdRegions::ConstFactorMap   &factors) final
        {
            boost::ignore_unused(entry2,entry3,wsp,factors);

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
            Vmath::Vcopy(m_nqtot, input, 1, m_input, 1);
            (*m_oper)(m_input, m_output);
            Vmath::Vcopy(m_nqtot, m_output, 1, output0, 1);
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
    unsigned int m_nqtot; 
    
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
        
        m_nqtot = nElmtNoPad*nqElmt; 

        if (nElmtNoPad % vec_t::width != 0)
        {
            m_isPadded = true;
            nElmtPad = nElmtNoPad + vec_t::width -
                (nElmtNoPad % vec_t::width);
            m_input = Array<OneD, NekDouble>{nmElmt * nElmtPad, 0.0};
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
};

}
}
