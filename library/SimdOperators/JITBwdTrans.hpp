#ifndef JITBWDTRANS_HPP
#define JITBWDTRANS_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Basis.h>

#include <boost/align/aligned_allocator.hpp>

#include <immintrin.h>
#include <asmjit/asmjit.h>

using namespace Nektar;
using namespace asmjit;

struct JITBwdTransSeg
{
    typedef void (*JITFunc)(const double *, double *, double *);

    JITBwdTransSeg(LibUtilities::BasisSharedPtr basis, int nElmt)
        : m_basis(basis), m_bdata(basis->GetBdata().num_elements()),
          m_nElmt(nElmt), m_nElmtBlocks(nElmt / 4)
    {
        // Create vector width basis
        const Array<OneD, const NekDouble> bdata = basis->GetBdata();
        m_bdata.resize(bdata.num_elements());

        for (auto i = 0; i < m_bdata.size(); ++i)
        {
            m_bdata[i] = _mm256_set1_pd(bdata[i]);
        }
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out)
    {
        const int nm = m_basis->GetNumModes();
        const int nq = m_basis->GetNumPoints();
        const int nmBlocks = nm*nm*4, nqBlocks = nq*nq*4;
        auto *inptr = &in[0];
        auto *outptr = &out[0];

        // Build JIT kernel for single operation.
        JitRuntime runtime;
        CodeHolder code;

        code.init(runtime.getCodeInfo());
        X86Assembler c(&code);

        auto bs = sizeof(__m256d);

        for (int i = 0; i < nq; ++i)
        {
            c.vmovups(x86::ymm0, x86::ptr(x86::rdi));
            c.vmovups(x86::ymm1, x86::ptr(x86::rsi, i*bs));
            c.vmulpd(x86::ymm0, x86::ymm0, x86::ymm1);

            for (int p = 1; p < nm; ++p)
            {
                c.vmovups(x86::ymm1, x86::ptr(x86::rdi, p*bs));
                c.vmovups(x86::ymm2, x86::ptr(x86::rsi, (p*nq+i)*bs));
                c.vfmadd231pd(x86::ymm0, x86::ymm1, x86::ymm2);
            }

            c.vmovupd(x86::ptr(x86::rdx, i*bs), x86::ymm0);
        }

        c.ret();
        Error err = runtime.add(&m_func, &code);
        std::cout << "err = " << err << std::endl;

        for (int e = 0; e < m_nElmtBlocks; ++e)
        {
            m_func(inptr, (double *)&m_bdata[0], outptr);

            inptr += nmBlocks;
            outptr += nqBlocks;
        }
    }

private:
    LibUtilities::BasisSharedPtr m_basis;
    /// Padded basis
    std::vector<__m256d, boost::alignment::aligned_allocator<__m256d, 64>> m_bdata;
    int m_nElmt;
    int m_nElmtBlocks;
    JITFunc m_func;
};

inline void dot(X86Assembler &c, int n, X86Mem a, int aInc, X86Mem b, int bInc, X86Mem out)
{
    auto bs = sizeof(__m256);
    X86Ymm regs[16] = {
        x86::ymm0, x86::ymm1, x86::ymm2, x86::ymm3, x86::ymm4, x86::ymm5,
        x86::ymm6, x86::ymm7, x86::ymm8, x86::ymm9, x86::ymm10, x86::ymm11,
        x86::ymm12, x86::ymm13, x86::ymm14, x86::ymm15
    };

    int cnt = n;

    while (cnt > 0)
    {
        int nreg = 16, rhlf = nreg / 2, toproc = std::min(cnt, rhlf);
        //std::cout << "START: cnt = " << cnt << " rhlf = " << rhlf << " toproc = " << toproc << endl; 

        // accumulate into ymm0. need two other registers for everything else.

        // first do loads of a and b
        for (int regcnt = 0; regcnt < toproc; ++regcnt)
        {
            //std::cout << "ymm" << regcnt << " <- a" << std::endl;
            c.vmovupd(regs[regcnt], a);
            a.addOffset(aInc * bs);
        }

        for (int regcnt = rhlf; regcnt < rhlf+toproc; ++regcnt)
        {
            //std::cout << "ymm" << regcnt << " <- b" << std::endl;
            c.vmovupd(regs[regcnt], b);
            b.addOffset(bInc * bs);
        }

        // First time around, vmul into ymm0
        if (cnt == n)
        {
            c.vmulpd(regs[0], regs[0], regs[rhlf]);
            //std::cout << "vmulpd ymm0, ymm1, ymm" << rhlf << std::endl;
            for (int p = 1; p < toproc; ++p)
            {
                c.vfmadd231pd(regs[0], regs[p], regs[rhlf+p]);
                //                std::cout << "vfmadd231pd ymm0, ymm" << p+1 << ", ymm" << rhlf+p << std::endl;
            }
        }
        else
        {
            std::cout <<" LOLL" << std::endl;
            // do somethung lol
            for (int p = 0; p < toproc; ++p)
            {
                c.vfmadd231pd(regs[0], regs[p+1], regs[rhlf+p+1]);
                //std::cout << "vfmadd231pd ymm0, ymm" << p+1 << ", ymm" << rhlf+p << std::endl;
            }
        }

        cnt -= toproc;
        //std::cout << "END: cnt = " << cnt << std::endl;
    }

    // c.vmovupd(x86::ymm0, a);
    // c.vmovupd(x86::ymm1, b);

    // for (int p = 1; p < n; ++p)
    // {
    //     a.addOffset(aInc * bs);
    //     b.addOffset(bInc * bs);
    //     c.vmovupd(x86::ymm1, a);
    //     c.vmovupd(x86::ymm2, b);
    //     c.vfmadd231pd(x86::ymm0, x86::ymm1, x86::ymm2);
    // }

    c.vmovupd(out, x86::ymm0);
}

struct JITBwdTransQuad
{
    typedef void (*JITFunc)(const double *, double *, double *);

    JITBwdTransQuad(LibUtilities::BasisSharedPtr basis, int nElmt)
        : m_basis(basis), m_bdata(basis->GetBdata().num_elements()),
          m_nElmt(nElmt), m_nElmtBlocks(nElmt / 4)
    {
        // Create vector width basis
        const Array<OneD, const NekDouble> bdata = basis->GetBdata();
        m_bdata.resize(bdata.num_elements());

        for (auto i = 0; i < m_bdata.size(); ++i)
        {
            m_bdata[i] = _mm256_set1_pd(bdata[i]);
        }

        // Build JIT kernel for single operation.
        CodeHolder code;
        code.init(m_runtime.getCodeInfo());
        X86Assembler c(&code);

        const int nm = m_basis->GetNumModes();
        const int nq = m_basis->GetNumPoints();
        const int nmBlocks = nm*nm, nqBlocks = nq*nq;

        auto bs = sizeof(__m256d);

        Label loop = c.newLabel();

        // Allocate wsp on the stack.
        c.push(x86::rbp);
        c.mov(x86::rbp, x86::rsp);
        c.sub(x86::rsp, nm * nq * bs);

        int cnt_iq = 0, cnt_pq, cnt_ij = 0;

        // registers: rdi = in
        //            rsi = bdata
        //            rdx = out
        // stack:     rsp = wsp

        c.mov(x86::rax, m_nElmtBlocks);
        c.bind(loop);

        for (int i = 0; i < nq; ++i)
        {
            cnt_pq = 0;
            for (int q = 0; q < nm; ++q, ++cnt_iq)
            {
                dot(c, nm, x86::ptr(x86::rdi, cnt_pq*bs), 1, x86::ptr(x86::rsi, i*bs), nq, x86::ptr(x86::rsp, cnt_iq*bs));
                cnt_pq += nm;
            }
        }

        cnt_ij = 0;
        for (int j = 0; j < nq; ++j)
        {
            cnt_iq = 0;
            for (int i = 0; i < nq; ++i, ++cnt_ij)
            {
                dot(c, nm, x86::ptr(x86::rsp, cnt_iq*bs), 1, x86::ptr(x86::rsi, j*bs), nq, x86::ptr(x86::rdx, cnt_ij*bs));
                cnt_iq += nm;
            }
        }

        // update registers to next block
        c.add(x86::rdi, nmBlocks * bs);
        c.add(x86::rdx, nqBlocks * bs);

        // jump if we're not done
        c.dec(x86::rax);
        c.jnz(loop);

        // restore stack
        c.mov(x86::rsp, x86::rbp);
        c.pop(x86::rbp);
        c.ret();

        Error err = m_runtime.add(&m_func, &code);
    }

    void operator()(const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out)
    {
        m_func(&in[0], (double *)&m_bdata[0], &out[0]);
    }

private:
    LibUtilities::BasisSharedPtr m_basis;
    /// Padded basis
    std::vector<__m256d, boost::alignment::aligned_allocator<__m256d, 64>> m_bdata;
    int m_nElmt;
    int m_nElmtBlocks;
    JITFunc m_func;
    JitRuntime m_runtime;
 };

#endif
