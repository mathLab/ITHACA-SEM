#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include <iostream>

#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <Collections/Operator.h>
#include "VecData.hpp"
// #include "AVXAssembly.h"

namespace Nektar {
namespace AVX {


/// Operator base class
class Operator
{
public:
    virtual ~Operator()
    {
    }

    /// Returns the vector width that's used in the operator, e.g. AVX = 4,
    /// AVX-512 = 8.
    virtual size_t VectorWidth()
    {
        return 1;
    }

    /// Number of gigflops required to compute the operator.
    virtual NekDouble GFlops()
    {
        return 0.0;
    }

    /// Number of degrees of freedom that this operator will process.
    virtual NekDouble Ndof()
    {
        return 0.0;
    }

    /// Number of stores to memory
    virtual NekDouble NStores()
    {
        return 0.0;
    }

    /// Number of loads from memory
    virtual NekDouble NLoads()
    {
        return 0.0;
    }

    /// This operator requires derivative factors.
    virtual bool NeedsDF()
    {
        return false;
    }

    /// This operator requires Jacobian.
    virtual bool NeedsJac()
    {
        return false;
    }

    virtual void SetJac(const Array<OneD, const NekDouble> &jac) = 0;

    virtual void SetDF(const Array<TwoD, const NekDouble> &df) = 0;

    // virtual void set_asmMap(MultiRegions::AssemblyMapSharedPtr asmMap)
    // {
    // }

    // /// Provides a reference function for operators which can be used to
    // /// validate correctness compared to Nektar++.
    // ///
    // /// In 2D this is taken as \f$ \sin x \cos y \f$, and in 3D \f$ \sin x \cos
    // /// y \sin z \f$.
    // static void RefFn(MultiRegions::ExpListSharedPtr expList,
    //                   Array<OneD, NekDouble> &in)
    // {
    //     const int nq = expList->GetNpoints();
    //     const int dim = expList->GetExp(0)->GetShapeDimension();

    //     if(dim == 2){
    //         Array<OneD, NekDouble> xc(nq), yc(nq);

    //         expList->GetCoords(xc,yc);
    //         for(int i = 0; i < nq; i++){
    //             in[i] = sin(xc[i]) * cos(yc[i]);
    //         }
    //     }
    //     else if(dim == 3){
    //         Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq);

    //         expList->GetCoords(xc,yc,zc);
    //         for(int i = 0; i < nq; i++){
    //             in[i] = sin(xc[i]) * cos(yc[i]) * sin(zc[i]);
    //         }
    //     }
    // }

};

typedef std::shared_ptr<Operator> OperatorSharedPtr;

/// Base class for backwards transform operator.
class BwdTrans : virtual public Operator
{
public:
    BwdTrans(std::vector<LibUtilities::BasisSharedPtr> basis,
             int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~BwdTrans()
    {
    }

    virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0; //Abstract Method

    // void Ref(MultiRegions::ExpListSharedPtr expList,
    //          Array<OneD, NekDouble> &ref_exp,
    //          Array<OneD, NekDouble> &ref_bwd)
    // {
    //     Array<OneD, NekDouble> ref_fn(expList->GetNpoints());
    //     this->RefFn(expList, ref_fn);
    //     expList->FwdTrans_IterPerExp(ref_fn, ref_exp);
    //     expList->BwdTrans(ref_exp, ref_bwd);
    // }

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class IProduct : virtual public Operator
{
public:
    IProduct(std::vector<LibUtilities::BasisSharedPtr> basis,
             int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~IProduct()
    {
    }

    virtual bool NeedsJac() override
    {
        return true;
    }

    virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0;

    // void Ref(MultiRegions::ExpListSharedPtr expList,
    //             Array<OneD, NekDouble> &ref_fn,
    //             Array<OneD, NekDouble> &ref_iprod)
    // {
    //     this->RefFn(expList, ref_fn);
    //     expList->IProductWRTBase(ref_fn, ref_iprod);
    // }

protected:
    /// Vector of tensor product basis directions
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class PhysDeriv : virtual public Operator
{
public:
    PhysDeriv(std::vector<LibUtilities::BasisSharedPtr> basis,
             int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~PhysDeriv()
    {
    }

    virtual bool NeedsDF() override
    {
        return true;
    }

    virtual void operator()(const Array<OneD, const NekDouble> &in,
                                Array<OneD,       NekDouble> &out_d0,
                                Array<OneD,       NekDouble> &out_d1) = 0;

    virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output0,
        Array<OneD, NekDouble> &output1,
        Array<OneD, NekDouble> &output2) = 0;

    // void Ref(MultiRegions::ExpListSharedPtr expList,
    //             Array<OneD, NekDouble> &ref_fn,
    //             Array<OneD, NekDouble> &d0,
    //             Array<OneD, NekDouble> &d1,
    //             Array<OneD, NekDouble> &d2 = NullNekDouble1DArray)
    // {
    //     const int dim = expList->GetExp(0)->GetShapeDimension();

    //     this->RefFn(expList, ref_fn);

    //     if(dim == 2){
    //         expList->PhysDeriv(ref_fn, d0, d1);
    //     }
    //     else if(dim == 3){
    //         expList->PhysDeriv(ref_fn, d0, d1, d2);
    //     }
    // }

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class IProductWRTDerivBase : virtual public Operator
{
public:
    IProductWRTDerivBase(std::vector<LibUtilities::BasisSharedPtr> basis,
              int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~IProductWRTDerivBase()
    {
    }

    virtual bool NeedsJac() override
    {
        return true;
    }

    virtual bool NeedsDF() override
    {
        return true;
    }

    virtual void operator()(
        const Array<OneD, Array<OneD, NekDouble>> &in,
        Array<OneD, NekDouble> &out) = 0;

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

class Helmholtz : virtual public Operator
{
public:
    Helmholtz(std::vector<LibUtilities::BasisSharedPtr> basis,
              int nElmt) :
        m_basis(basis), m_nElmt(nElmt)
    {
    }

    virtual ~Helmholtz()
    {
    }

    virtual bool NeedsJac() override
    {
        return true;
    }

    virtual bool NeedsDF() override
    {
        return true;
    }

    virtual void operator()(
        const Array<OneD, const NekDouble> &input,
        Array<OneD, NekDouble> &output) = 0;

    // void Ref(MultiRegions::ExpListSharedPtr expList,
    //          Array<OneD, NekDouble> &ref_exp,
    //          Array<OneD, NekDouble> &ref_helm)
    // {
    //     Array<OneD, NekDouble> ref_fn(expList->GetNpoints());
    //     this->RefFn(expList, ref_fn);
    //     expList->FwdTrans_IterPerExp(ref_fn, ref_exp);

    //     StdRegions::ConstFactorMap factors;
    //     factors[StdRegions::eFactorLambda] = 1.0;

    //     MultiRegions::GlobalMatrixKey mkey(
    //         StdRegions::eHelmholtz,
    //         MultiRegions::NullAssemblyMapSharedPtr,
    //         factors);

    //     expList->GeneralMatrixOp_IterPerExp(mkey, ref_exp, ref_helm);
    // }

protected:
    std::vector<LibUtilities::BasisSharedPtr> m_basis;
    int m_nElmt;
};

// template<int VW>
// class HelmholtzGlobal : virtual public Operator
// {
//     public:
//     HelmholtzGlobal(std::vector<LibUtilities::BasisSharedPtr> basis,
//             int nElmt) :
//             m_basis(basis), m_nElmt(nElmt)
//     {
//     }

//     virtual ~HelmholtzGlobal()
//     {
//     }

//     virtual bool NeedsJac() override
//     {
//         return true;
//     }

//     virtual bool NeedsDF() override
//     {
//         return true;
//     }

//     virtual void operator()(
//         const Array<OneD, const NekDouble> &globalIn,
//         Array<OneD,       NekDouble> &globalOut,
//         AVXAssembly<VW>              &l2g) = 0;


//     // void Ref(MultiRegions::ExpListSharedPtr expList,
//     //          Array<OneD, NekDouble> &ref_exp,
//     //          Array<OneD, NekDouble> &ref_helm)
//     // {
//     //     Array<OneD, NekDouble> ref_fn(expList->GetNpoints());
//     //     this->RefFn(expList, ref_fn);
//     //     expList->FwdTrans_IterPerExp(ref_fn, ref_exp);

//     //     StdRegions::ConstFactorMap factors;
//     //     factors[StdRegions::eFactorLambda] = 1.0;

//     //     MultiRegions::GlobalMatrixKey mkey(
//     //         StdRegions::eHelmholtz,
//     //         MultiRegions::NullAssemblyMapSharedPtr,
//     //         factors);

//     //     expList->GeneralMatrixOp_IterPerExp(mkey, ref_exp, ref_helm);
//     // }

//     // void Ref2(MultiRegions::ExpListSharedPtr expList,
//     //         Array<OneD, NekDouble> &ref_exp,
//     //         Array<OneD, NekDouble> &ref_helm,
//     //         const MultiRegions::AssemblyMapSharedPtr asmMap)
//     // {
//     //     Array<OneD, NekDouble> ref_fn(expList->GetNpoints());
//     //     Array<OneD, NekDouble> ref_exp_local (expList->GetNcoeffs());
//     //     this->RefFn(expList, ref_fn);
//     //     expList->FwdTrans_IterPerExp(ref_fn, ref_exp_local);

//     //     StdRegions::ConstFactorMap factors;
//     //     factors[StdRegions::eFactorLambda] = 1.0;

//     //     MultiRegions::GlobalMatrixKey mkey(
//     //         StdRegions::eHelmholtz,
//     //         MultiRegions::NullAssemblyMapSharedPtr,
//     //         factors);

//     //     Array<OneD, NekDouble> ref_helm_local(expList->GetNcoeffs());
//     //     expList->GeneralMatrixOp_IterPerExp(mkey, ref_exp_local, ref_helm_local);

//     //     asmMap->Assemble(ref_helm_local, ref_helm);
//     // }

// protected:
//     std::vector<LibUtilities::BasisSharedPtr> m_basis;
//     int m_nElmt;

// };

template <int VW, int DIM, bool DEFORMED = false>
class AVXHelper : virtual public Operator
{
protected:
    AVXHelper(std::vector<LibUtilities::BasisSharedPtr> basis,
              int nElmt)
        : Operator()
    {
        // Sanity check: no padding yet!
        ASSERTL0(nElmt % VW == 0, "Number of elements not divisible by vector "
                 "width, padding not yet implemented.");

        // Calculate number of 'blocks', i.e. meta-elements
        m_nBlocks = nElmt / VW;

        // Depending on element dimension, set up basis information, quadrature,
        // etc, inside vectorised environment.
        for (int i = 0; i < DIM; ++i)
        {
            const Array<OneD, const NekDouble> bdata = basis[i]->GetBdata();
            const Array<OneD, const NekDouble> dbdata = basis[i]->GetDbdata();
            const Array<OneD, const NekDouble> w = basis[i]->GetW();

            m_nm[i] = basis[i]->GetNumModes();
            m_nq[i] = basis[i]->GetNumPoints();

            m_bdata[i].resize(bdata.num_elements());
            m_bdata_orig[i].resize(bdata.num_elements());
            for (auto j = 0; j < bdata.num_elements(); ++j)
            {
                m_bdata[i][j] = bdata[j];
                m_bdata_orig[i][j] = bdata[j];
            }

            m_dbdata[i].resize(dbdata.num_elements());
            m_dbdata_orig[i].resize(dbdata.num_elements());
            for (auto j = 0; j < dbdata.num_elements(); ++j)
            {
                m_dbdata[i][j] = dbdata[j];
                m_dbdata_orig[i][j] = dbdata[j];
            }

            NekDouble fac = 1.0;
            if (basis[i]->GetPointsType() == LibUtilities::eGaussRadauMAlpha1Beta0)
            {
                fac = 0.5;
            }
            else if (basis[i]->GetPointsType() == LibUtilities::eGaussRadauMAlpha2Beta0)
            {
                fac = 0.25;
            }

            m_w[i].resize(w.num_elements());
            m_w_orig[i].resize(w.num_elements());
            for (auto j = 0; j < w.num_elements(); ++j)
            {
                m_w[i][j] = fac * w[j];
                m_w_orig[i][j] = fac * w[j];
            }

            auto D = basis[i]->GetD()->GetPtr();
            m_D[i].resize(D.num_elements());
            m_D_orig[i].resize(D.num_elements());
            for(int j = 0; j < D.num_elements(); j++){
                m_D[i][j] = D[j];
                m_D_orig[i][j] = D[j];
            }

            auto Z = basis[i]->GetZ();
            m_Z[i].resize(Z.num_elements());
            m_Z_orig[i].resize(Z.num_elements());
            for(int j = 0; j < Z.num_elements(); j++){
                m_Z[i][j] = Z[j];
                m_Z_orig[i][j] = Z[j];
            }

        }
    }

    /// Set up Jacobian array for those operators that require geometric
    /// information.
    void SetJac(const Array<OneD, const NekDouble> &jac) override
    {
        if (DEFORMED)
        {
            int nq = m_nq[0];
            for (int i = 1; i < DIM; i++)
            {
                nq *= m_nq[i];
            }

            m_jac.resize(m_nBlocks*nq);

            AlignedVector<NekDouble> tmp(VW);
            for (size_t block = 0; block < m_nBlocks; ++block)
            {
                for(size_t q = 0; q < nq; q++)
                {
                    for (int j = 0; j < VW; ++j)
                    {
                        tmp[j] = jac[block*nq*VW + nq*j + q]; //Unvalidated until I can get an actual deformed mesh.
                    }

                    //Order is [block][quadpt]
                    m_jac[block*nq + q] = &tmp[0];
                }
            }
        }
        else{
            m_jac.resize(m_nBlocks);

            AlignedVector<NekDouble> tmp(VW);
            for (size_t i = 0; i < m_nBlocks; ++i)
            {
                for (int j = 0; j < VW; ++j)
                {
                    tmp[j] = jac[VW*i+j];
                }
                m_jac[i] = &tmp[0];
            }
        }
    }

    virtual void SetDF(const Array<TwoD, const NekDouble> &df) override
    {
        constexpr unsigned int n_df = DIM * DIM;
        AlignedVector<NekDouble> vec(VW);

        if (DEFORMED)
        {
            int nq = m_nq[0];
            for (int i = 1; i < DIM; ++i)
            {
                nq *= m_nq[i];
            }

            m_df.resize(m_nBlocks * n_df*nq);
            auto *df_ptr = &m_df[0];

            for (int e = 0; e < m_nBlocks; ++e)
            {
                for (int q = 0; q < nq; q++)
                {
                    for (int dir = 0; dir < n_df; ++dir, ++df_ptr)
                    {
                        for (int j = 0; j < VW; ++j)
                        {
                            vec[j] = df[dir][(VW*e + j)*nq + q];
                        }
                        *df_ptr = &vec[0];
                    }
                }
            }
        }
        else
        {
            m_df.resize(m_nBlocks * n_df);
            for (int e = 0; e < m_nBlocks; ++e)
            {
                for (int dir = 0; dir < n_df; ++dir)
                {
                    for(int j = 0; j < VW; ++j)
                    {
                        vec[j] = df[dir][VW*e + j];
                    }
                    // Must have all VW elemnts aligned to do a load.
                    m_df[e*n_df + dir] = &vec[0];
                }
            }
        }

    }

    virtual size_t VectorWidth() override
    {
        return VW;
    }

    int m_nBlocks;
    std::array<int, DIM> m_nm, m_nq;
    std::array<AlignedVector<VecData<NekDouble, VW>>, DIM> m_bdata;
    std::array<AlignedVector<VecData<NekDouble, VW>>, DIM> m_dbdata;
    std::array<AlignedVector<VecData<NekDouble, VW>>, DIM> m_D; //Derivatives
    std::array<AlignedVector<VecData<NekDouble, VW>>, DIM> m_Z; //Zeroes
    std::array<AlignedVector<VecData<NekDouble, VW>>, DIM> m_w; //Weights
    AlignedVector<VecData<NekDouble, VW>> m_df; //Chain rule function deriviatives for each element (00, 10, 20, 30...)
    AlignedVector<VecData<NekDouble, VW>> m_jac;

    // Original basis
    std::array<AlignedVector<NekDouble>, DIM> m_bdata_orig;
    std::array<AlignedVector<NekDouble>, DIM> m_dbdata_orig;
    std::array<AlignedVector<NekDouble>, DIM> m_D_orig; //Derivatives
    std::array<AlignedVector<NekDouble>, DIM> m_Z_orig; //Zeroes
    std::array<AlignedVector<NekDouble>, DIM> m_w_orig; //Weights
};

using OperatorFactory = LibUtilities::NekFactory<
    std::string,
    Operator,
    std::vector<LibUtilities::BasisSharedPtr>,
    int
    >;

OperatorFactory &GetOperatorFactory();

/// Helper function, get operator string
std::string GetOpstring(LibUtilities::ShapeType shape, bool deformed=false);

} // namespace AVX
} // namespace Nektar

#endif