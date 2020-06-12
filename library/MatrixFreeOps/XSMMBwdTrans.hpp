#ifndef XSMMBWDTRANS_HPP
#define XSMMBWDTRANS_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>
#include <boost/align/aligned_allocator.hpp>
#include <boost/align/is_aligned.hpp>

#include <libxsmm.h>

using namespace Nektar;

struct XSMMBwdTransQuad
{
    XSMMBwdTransQuad(LibUtilities::BasisSharedPtr basis, int nElmt)
        : m_basis(basis), m_bdata(basis->GetBdata().num_elements()), m_nElmt(nElmt),
          m_bdataTrans(basis->GetBdata().num_elements())
    {
        const int nq = basis->GetNumPoints(), nm = basis->GetNumModes();
        m_nm = nm;
        m_nq = nq;

        // Query/JIT functions
        NekDouble alpha = 1.0, beta = 0.0;
        m_func0 = libxsmm_dmmdispatch(
            m_nq, m_nm, m_nm, NULL, NULL, NULL, &alpha, &beta, NULL, NULL);
        m_func1 = libxsmm_dmmdispatch(
            m_nq, m_nq, m_nm, NULL, NULL, NULL, &alpha, &beta, NULL, NULL);
        ASSERTL0(m_func0 && m_func1, "Didn't compile");

        m_wsp.resize(nq * nm);

        for (int j = 0; j < nm; ++j)
        {
            for (int i = 0; i < nq; ++i)
            {
                m_bdata[i+j*nq] = basis->GetBdata()[i + j * nq];
                m_bdataTrans[j + i * nm] = m_bdata[i + j * nq];
            }
        }
    }

    void operator()(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        const int nqTot  = m_nq * m_nq;
        const int nmTot  = m_nm * m_nm;

        const NekDouble *in = &input[0];
        NekDouble *out = &output[0];
        NekDouble alpha = 1.0, beta = 0.0;

        for (int i = 0; i < m_nElmt; ++i)
        {
            // Direction 1
            m_func0(&m_bdata[0], in, &m_wsp[0]);

            // Direction 2
            m_func1(&m_wsp[0], &m_bdataTrans[0], out);

            in  += nmTot;
            out += nqTot;
        }
    }

private:
    LibUtilities::BasisSharedPtr m_basis;
    std::vector<NekDouble, boost::alignment::aligned_allocator<NekDouble, 64>> m_bdata, m_bdataTrans, m_wsp;
    libxsmm_dmmfunction m_func0, m_func1;
    //libxsmm_mmfunction<NekDouble> m_func0, m_func1;
    int m_nElmt, m_nq, m_nm;
};

struct XSMMBwdTransHex
{
    XSMMBwdTransHex(LibUtilities::BasisSharedPtr basis0,
                    LibUtilities::BasisSharedPtr basis1,
                    LibUtilities::BasisSharedPtr basis2,
                    int nElmt)
        : m_basis0(basis0), m_basis1(basis1), m_basis2(basis2),
          m_nq0(basis0->GetNumPoints()), m_nm0(basis0->GetNumModes()),
          m_nq1(basis1->GetNumPoints()), m_nm1(basis1->GetNumModes()),
          m_nq2(basis2->GetNumPoints()), m_nm2(basis2->GetNumModes()),
          m_bdata0   (basis0->GetBdata()),
          //m_bdata0(basis0->GetBdata().num_elements()),
          //m_bdata1(basis1->GetBdata().num_elements()),
          //m_bdata2(basis2->GetBdata().num_elements()),
          m_bdata1T(basis1->GetBdata().num_elements()),
          m_bdata2T(basis2->GetBdata().num_elements()),
          m_nElmt(nElmt)
    {
        const int nm12   = m_nm1 * m_nm2;
        const int nq0nm1 = m_nq0 * m_nm1;
        const int nq01   = m_nq0 * m_nq1;

        m_wsp .resize(m_nq0 * m_nm1 * m_nm2);
        m_wsp2.resize(m_nq0 * m_nq1 * m_nm2);

        // Query/JIT functions
        NekDouble alpha = 1.0, beta = 0.0;
        m_func0 = libxsmm_dmmdispatch(
            m_nq0, nm12, m_nm0, NULL, NULL, NULL, &alpha, &beta, NULL, NULL);
        m_func1 = libxsmm_dmmdispatch(
            m_nq0, m_nq1, m_nm1, NULL, NULL, NULL, &alpha, &beta, NULL, NULL);
        m_func2 = libxsmm_dmmdispatch(
            nq01, m_nq2, m_nm2, NULL, NULL, NULL, &alpha, &beta, NULL, NULL);

        Array<OneD, const NekDouble> base1 = basis1->GetBdata();
        Array<OneD, const NekDouble> base2 = basis2->GetBdata();

        for (int j = 0; j < m_nm1; ++j)
        {
            for (int i = 0; i < m_nq1; ++i)
            {
                m_bdata1T[j + i * m_nm1] = base1[i + j * m_nq1];
            }
        }

        for (int j = 0; j < m_nm2; ++j)
        {
            for (int i = 0; i < m_nq2; ++i)
            {
                m_bdata2T[j + i * m_nm2] = base2[i + j * m_nq2];
            }
        }
    }

    void operator()(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        const int nqTot  = m_nq0 * m_nq1 * m_nq2;
        const int nmTot  = m_nm0 * m_nm1 * m_nm2;
        const int nm12   = m_nm1 * m_nm2;
        const int nq0nm1 = m_nq0 * m_nm1;
        const int nq01   = m_nq0 * m_nq1;

        const NekDouble *in = &input[0];
        NekDouble *out = &output[0], *tmp1, *tmp2;

        for (int i = 0; i < m_nElmt; ++i)
        {
            tmp1 = &m_wsp[0];
            tmp2 = &m_wsp2[0];

            // Direction 1
            m_func0(&m_bdata0[0], in, tmp1);

            // Loop for direction 2
            for (int r = 0; r < m_nm2; ++r)
            {
                m_func1(tmp1, &m_bdata1T[0], tmp2);
                tmp1 += nq0nm1;
                tmp2 += nq01;
            }

            // Third direction
            m_func2(&m_wsp2[0], &m_bdata2T[0], out);

            in  += nmTot;
            out += nqTot;
        }
    }

private:
    LibUtilities::BasisSharedPtr m_basis0, m_basis1, m_basis2;
    //std::vector<NekDouble, boost::alignment::aligned_allocator<NekDouble, 64>>
    //m_bdata0, m_bdata1, m_bdata2, m_bdata1T, m_bdata2T, m_wsp, m_wsp2;
    Array<OneD, const NekDouble>    m_bdata0;
    Array<OneD, NekDouble>          m_bdata1T;
    Array<OneD, NekDouble>          m_bdata2T;
    std::vector<NekDouble>          m_wsp;
    std::vector<NekDouble>          m_wsp2;
    libxsmm_dmmfunction m_func0, m_func1, m_func2;
    int m_nElmt, m_nq0, m_nq1, m_nq2, m_nm0, m_nm1, m_nm2;
};

struct XSMMBwdTransTet
{
    XSMMBwdTransTet(LibUtilities::BasisSharedPtr basis0,
                    LibUtilities::BasisSharedPtr basis1,
                    LibUtilities::BasisSharedPtr basis2,
                    int nElmt)
        : m_basis0(basis0), m_basis1(basis1), m_basis2(basis2),
          m_nm0(basis0->GetNumModes()),
          m_nm1(basis1->GetNumModes()),
          m_nm2(basis2->GetNumModes()),
          m_nq0(basis0->GetNumPoints()),
          m_nq1(basis1->GetNumPoints()),
          m_nq2(basis2->GetNumPoints()),
          m_bdata0(basis0->GetBdata().get(),
                   basis0->GetBdata().get()+basis0->GetBdata().num_elements()),
          m_bdata1(basis1->GetBdata().get(),
                   basis1->GetBdata().get()+basis1->GetBdata().num_elements()),
          m_bdata2(basis2->GetBdata().get(),
                   basis2->GetBdata().get()+basis2->GetBdata().num_elements()),
          m_nElmt(nElmt)
    {
        // Query/JIT functions
        NekDouble alpha = 1.0, beta = 0.0;

        for(int i = 0, cnt = 0; i < m_nm0; ++i)
        {
            for(int j = 0; j < m_nm1-i; ++j, ++cnt)
            {
                m_func0.push_back(
                    libxsmm_dmmdispatch(
                        m_nq2, 1, m_nm2-i-j, NULL, NULL, NULL, &alpha, &beta, NULL, NULL));
            }

            m_func1.push_back(
                libxsmm_dmmdispatch(
                    m_nq1, m_nq2, m_nm1-i, NULL, NULL, NULL, &alpha, &beta, NULL, NULL));
        }

        m_func2 = libxsmm_dmmdispatch(
            m_nq0, m_nq1*m_nq2, m_nm0, NULL, NULL, NULL, &alpha, &beta, NULL, NULL);

        m_wsp.resize(m_nq2*m_nm0*(2*m_nm1-m_nm0+1)/2+m_nq2*m_nq1*m_nm0);
    }

    void operator()(
        const Array<OneD, const NekDouble> &input,
              Array<OneD,       NekDouble> &output)
    {
        NekDouble *tmp = &m_wsp[0], *tmp1 = &m_wsp[m_nq2*m_nm0*(2*m_nm1-m_nm0+1)/2];
        int i, j, mode, mode1, cnt;

        std::vector<NekDouble, boost::alignment::aligned_allocator<NekDouble, 64>>
            tmp1_trans(m_nq1 * m_nq2 * m_nm0);

        NekDouble const *inptr = &input[0];
        NekDouble *outptr = &output[0];

        const int nqTot  = m_nq0 * m_nq1 * m_nq2;
        const int nmTot  = LibUtilities::StdTetData::getNumberOfCoefficients(
            m_nm0, m_nm1, m_nm2);

        for (int e = 0; e < m_nElmt; ++e)
        {
            // Perform summation over '2' direction
            mode = mode1 = cnt = 0;
            for(i = 0; i < m_nm0; ++i)
            {
                for(j = 0; j < m_nm1-i; ++j, ++cnt)
                {
                    m_func0[cnt](&m_bdata2[mode*m_nq2], inptr + mode1, tmp + cnt*m_nq2);
                    //Blas::Dgemv('N', m_nq2, m_nm2-i-j, 1.0,
                    //            &m_bdata2[mode*m_nq2], m_nq2, inptr + mode1,
                    //            1, 0.0, tmp+cnt*m_nq2, 1);
                    mode  += m_nm2-i-j;
                    mode1 += m_nm2-i-j;
                }

                //increment mode in case order1!=order2
                for(j = m_nm1-i; j < m_nm2-i; ++j)
                {
                    mode += m_nm2-i-j;
                }
            }

            // fix for modified basis by adding split of top singular vertex
            // mode - currently (1+c)/2 x (1-b)/2 x (1-a)/2 component is
            // evaluated top singular vertex - (1+c)/2 x (1+b)/2 x (1-a)/2
            // component
            Blas::Daxpy(m_nq2, inptr[1], &m_bdata2[m_nq2], 1, &tmp[m_nq2], 1);

            // top singular vertex - (1+c)/2 x (1-b)/2 x (1+a)/2 component
            Blas::Daxpy(m_nq2, inptr[1], &m_bdata2[m_nq2], 1, &tmp[m_nm1*m_nq2], 1);

            // Perform summation over '1' direction
            mode = 0;
            for(i = 0; i < m_nm0; ++i)
            {
                //libxsmm_otrans(&tmp1_trans[0], tmp + mode*m_nq2,
                //               sizeof(NekDouble), m_nq2, m_nm1-i, m_nq2, m_nm1-i);

                m_func1[i](&m_bdata1[mode*m_nq1], &tmp1_trans[0], tmp1+i*m_nq1*m_nq2);

                // Blas::Dgemm('N', 'N', m_nq1, m_nq2, m_nm1-i, 1.0,
                //             &m_bdata1[mode*m_nq1], m_nq1, &tmp1_trans[0],
                //             m_nm1-i, 0.0, tmp1+i*m_nq1*m_nq2, m_nq1);
                mode  += m_nm1-i;
            }

            for(i = 0; i < m_nq2; ++i)
            {
                Blas::Daxpy(m_nq1, tmp[m_nq2+i], &m_bdata1[m_nq1], 1,
                            &tmp1[m_nq1*m_nq2]+i*m_nq1, 1);
            }

            //libxsmm_otrans(&tmp1_trans[0], tmp1, sizeof(NekDouble),
            //               m_nq1*m_nq2, m_nm0, m_nq1*m_nq2, m_nm0);

            // Perform summation over '0' direction
            m_func2(&m_bdata0[0], &tmp1_trans[0], outptr);
            // Blas::Dgemm('N', 'N', m_nq0, m_nq1*m_nq2, m_nm0,
            //             1.0, &m_bdata0[0], m_nq0,
            //             &tmp1_trans[0], m_nm0, 0.0, outptr, m_nq0);

            inptr += nmTot;
            outptr += nqTot;
        }
    }

private:
    LibUtilities::BasisSharedPtr m_basis0, m_basis1, m_basis2;
    std::vector<NekDouble, boost::alignment::aligned_allocator<NekDouble, 64>>
        m_bdata0, m_bdata1, m_bdata2, m_wsp;
    std::vector<libxsmm_dmmfunction> m_func0, m_func1;
    libxsmm_dmmfunction m_func2;
    int m_nElmt, m_nq0, m_nq1, m_nq2, m_nm0, m_nm1, m_nm2;
};

#endif
