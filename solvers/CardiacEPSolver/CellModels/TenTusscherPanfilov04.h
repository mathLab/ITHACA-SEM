///////////////////////////////////////////////////////////////////////////////
//
// File TenTusscherPanfilov04.h
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
// Description: TenTusscherPanfilov04 cell model
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_TENTUSSCHERPANFILOV04_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_TENTUSSCHERPANFILOV04_H

#include <CardiacEPSolver/CellModels/CellModel.h>

namespace Nektar
{
    class TenTusscherPanfilov04: public CellModel
    {

    public:
        /// Creates an instance of this class
        static CellModelSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const int nq)
        {
            return MemoryManager<TenTusscherPanfilov04>
                                        ::AllocateSharedPtr(pSession, nq);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        TenTusscherPanfilov04(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const int nq);

        /// Destructor
        virtual ~TenTusscherPanfilov04();

    protected:
        /// Computes the reaction terms $f(u,v)$ and $g(u,v)$.
        virtual void v_Update(
                const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        /// Prints a summary of the model parameters.
        virtual void v_PrintSummary(std::ostream &out);

    private:
        int cell_type;
        NekDouble R;
        NekDouble T;
        NekDouble F;
        NekDouble C_m;
        NekDouble S;
        NekDouble rho;
        NekDouble V_C;
        NekDouble V_SR;
        NekDouble K_o;
        NekDouble Na_o;
        NekDouble Ca_o;
        NekDouble g_Na;
        NekDouble g_K1;
        NekDouble g_to;
        NekDouble g_Kr;
        NekDouble g_Ks;
        NekDouble p_KNa;
        NekDouble g_Ca_L;
        NekDouble k_NaCa;
        NekDouble gamma;
        NekDouble K_m_Ca;
        NekDouble K_m_Na_i;
        NekDouble K_sat;
        NekDouble alpha;
        NekDouble P_NaK;
        NekDouble K_m_K;
        NekDouble K_m_Na;
        NekDouble g_pk;
        NekDouble g_pCa;
        NekDouble K_pCa;
        NekDouble g_b_Na;
        NekDouble g_b_Ca;
        NekDouble V_maxup;
        NekDouble K_up;
        NekDouble a_rel;
        NekDouble b_rel;
        NekDouble c_rel;
        NekDouble V_leak;
        NekDouble Buf_c;
        NekDouble K_bufc;
        NekDouble Buf_sr;
        NekDouble K_bufsr;
    };

}

#endif
