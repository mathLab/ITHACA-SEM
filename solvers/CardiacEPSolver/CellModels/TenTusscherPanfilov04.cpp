///////////////////////////////////////////////////////////////////////////////
//
// File TenTusscherPanfilov04.cpp
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
// Description: TenTusscherPanfilov04 ionic atrial cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <CardiacEPSolver/CellModels/TenTusscherPanfilov04.h>

using namespace std;

namespace Nektar
{
    std::string TenTusscherPanfilov04::className
              = GetCellModelFactory().RegisterCreatorFunction(
                        "TenTusscherPanfilov04",
                        TenTusscherPanfilov04::create,
                         "Ionic model of human atrial cell electrophysiology.");
    
    
    /**
    *
    */
    TenTusscherPanfilov04::TenTusscherPanfilov04(
                const LibUtilities::SessionReaderSharedPtr& pSession, const int nq)
            : CellModel(pSession, nq)
    {
        ASSERTL0(pSession->GetVariables().size() == 18,
                 "TenTusscherPanfilov04 cell model requires 18  variables.");

        m_nq = nq;
        cell_type = 0; // 0 = epi, 1 = M, 2 = endo

        R = 8.3143;
        T = 310.0;
        F = 96.4867;
        pSession->LoadParameter("Cm_ion", C_m);
        S = 0.2;
        rho = 162;
        V_C = 16404;
        V_SR = 1094;
        K_o = 5.4;      // millimolar
        Na_o = 140.0;   // millimolar
        Ca_o = 2.0;
        g_Na = 14.838;     // nanoS_per_picoF
        g_K1 = 5.405;    // nanoS_per_picoF
        g_to = (cell_type < 2) ? 0.294 : 0.073;  // nanoS_per_picoF
        g_Kr = 0.096;
        g_Ks = (cell_type % 2) ? 0.062 : 0.245;
        p_KNa = 0.03;
        g_Ca_L = 0.000175;
        k_NaCa = 1000;
        gamma = 0.35;
        K_m_Ca = 1.38;
        K_m_Na_i = 87.5;
        K_sat = 0.1;
        alpha = 2.5;
        P_NaK = 1.362;
        K_m_K = 1;
        K_m_Na = 40.0;
        g_pk = 0.0146;
        g_pCa = 0.025;
        K_pCa = 0.0005;
        g_b_Na = 0.00029;
        g_b_Ca = 0.000592;
        V_maxup = 0.000425;
        K_up = 0.00025;
        a_rel = 16.464;
        b_rel = 0.25;
        c_rel = 8.232;
        V_leak = 0.00008;
        Buf_c = 0.15;
        K_bufc = 0.001;
        Buf_sr = 10.0;
        K_bufsr = 0.3;
    }
    
    
    
    /**
    *
    */
    TenTusscherPanfilov04::~TenTusscherPanfilov04()
    {
        
    }
    
    
    
    void TenTusscherPanfilov04::v_Update(
                     const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                           Array<OneD,        Array<OneD, NekDouble> >&outarray,
                     const NekDouble time)
    {
        ASSERTL0(inarray.get() != outarray.get(),
                "Must have different arrays for input and output.");

        // Variables
        //  0   V    membrane potential
        //  1   -    unused (for bidomain only)
        //  2   m    fast sodium current m gate
        //  3   h    fast sodium current h gate
        //  4   j    fast sodium current j gate
        //  5   x_r1 rapid delayed rectifier K current gate
        //  6   x_r2 rapid delayed rectifier K current gate
        //  7   x_s  slow delayed rectifier K current gate
        //  8   d    L_type calcium gate
        //  9   f    L-type calcium gate
        //  10  f_Ca L-type calcium gate
        //  11  s    Ca release u gate
        //  12  r    Ca release v gate
        //  13  g    Ca release w gate
        //  14  Na_i Sodium
        //  15  Ca_i Calcium
        //  16  K_i  Potassium
        //  17  Ca_SR  Calcium up
        int nvariables = inarray.num_elements();
        int n          = m_nq;
        int i          = 0;
        NekDouble alpha, beta;

        Array<OneD, NekDouble> &tmp = outarray[12];
        Array<OneD, NekDouble> &tmp2 = outarray[13];

        // E_Na
        Array<OneD, NekDouble> &tmp_E_na = outarray[15];
        Vmath::Sdiv(n, Na_o, inarray[17], 1, tmp_E_na, 1);
        Vmath::Vlog(n, tmp_E_na, 1, tmp_E_na, 1);
        Vmath::Smul(n, R*T/F, tmp_E_na, 1, tmp_E_na, 1);

        // Sodium I_Na
        Array<OneD, NekDouble> &tmp_I_Na = outarray[16];
        Vmath::Vsub(n, inarray[0], 1, tmp_E_na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[2], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[2], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[2], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[3], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[4], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Smul(n, C_m*g_Na, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Na, 1, outarray[0], 1);
        Vmath::Smul(n, -1.0, tmp_I_Na, 1, outarray[17], 1);

        // Background current, sodium
        Array<OneD, NekDouble> &tmp_I_b_Na = outarray[13];
        Vmath::Vsub(n, inarray[0], 1, tmp_E_na, 1, tmp_I_b_Na, 1);
        Vmath::Smul(n, C_m*g_b_Na, tmp_I_b_Na, 1, tmp_I_b_Na, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_b_Na, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[17], 1, tmp_I_b_Na, 1, outarray[17], 1);

        // V - E_K
        Array<OneD, NekDouble> &tmp_V_E_k = outarray[12];
        Vmath::Sdiv(n, K_o, inarray[19], 1, tmp_V_E_k, 1);
        Vmath::Vlog(n, tmp_V_E_k, 1, tmp_V_E_k, 1);
        Vmath::Smul(n, R*T/F, tmp_V_E_k, 1, tmp_V_E_k, 1);
        Vmath::Vsub(n, inarray[0], 1, tmp_V_E_k, 1, tmp_V_E_k, 1);

        // Potassium I_K1
        /// @todo update I_K1
        Array<OneD, NekDouble> &tmp_I_K1 = outarray[13];
        Vmath::Sadd(n, 80.0, inarray[0], 1, tmp_I_K1, 1);
        Vmath::Smul(n, 0.07, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Vexp(n, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Sadd(n, 1.0, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Vdiv(n, tmp_V_E_k, 1, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Smul(n, C_m*g_K1, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_K1, 1, outarray[0], 1);
        Vmath::Smul(n, -1.0, tmp_I_K1, 1, outarray[19], 1);

        // Transient Outward K+ current
        Array<OneD, NekDouble> &tmp_I_to = outarray[13];
        Vmath::Vmul(n, inarray[12], 1, tmp_V_E_k, 1, tmp_I_to, 1);
        Vmath::Vmul(n, inarray[11], 1, tmp_I_to, 1, tmp_I_to, 1);
        Vmath::Smul(n, C_m*g_to, tmp_I_to, 1, tmp_I_to, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_to, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[19], 1, tmp_I_to, 1, outarray[19], 1);

        // Ultrarapid Delayed rectifier K+ current
//        Array<OneD, NekDouble> &tmp_I_kur = outarray[16];
//        Vmath::Sadd(n, -15.0, inarray[0], 1, tmp_I_kur, 1);
//        Vmath::Smul(n, -1.0/13.0, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Vexp(n, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Sadd(n, 1.0, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Sdiv(n, 0.05, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Sadd(n, 0.005, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Vmul(n, tmp_V_E_k,  1, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Vmul(n, inarray[7], 1, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Vmul(n, inarray[7], 1, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Vmul(n, inarray[7], 1, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Vmul(n, inarray[8], 1, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Smul(n, C_m, tmp_I_kur, 1, tmp_I_kur, 1);
//        Vmath::Vsub(n, outarray[0], 1, tmp_I_kur, 1, outarray[0], 1);
//        Vmath::Vsub(n, outarray[19], 1, tmp_I_kur, 1, outarray[19], 1);

        // Rapid delayed outward rectifier K+ current
        Array<OneD, NekDouble> &tmp_I_Kr = outarray[13];
        Vmath::Sadd(n, 15.0, inarray[0], 1, tmp_I_Kr, 1);
        Vmath::Smul(n, 1.0/22.4, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vexp(n, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Sadd(n, 1.0, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vdiv(n, tmp_V_E_k, 1, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vmul(n, inarray[9], 1, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Smul(n, C_m*g_Kr, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Kr, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[19], 1, tmp_I_Kr, 1, outarray[19], 1);

        // Slow delayed outward rectifier K+ Current
        Array<OneD, NekDouble> &tmp_I_Ks = outarray[13];
        Vmath::Vmul(n, inarray[10], 1, tmp_V_E_k, 1, tmp_I_Ks, 1);
        Vmath::Vmul(n, inarray[10], 1, tmp_I_Ks, 1, tmp_I_Ks, 1);
        Vmath::Smul(n, C_m*g_Ks, tmp_I_Ks, 1, tmp_I_Ks, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Ks, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[19], 1, tmp_I_Ks, 1, outarray[19], 1);

        // Background current, calcium
        Array<OneD, NekDouble> &tmp_I_b_Ca = outarray[2];
        Vmath::Sdiv(n, Ca_o, inarray[18], 1, tmp_I_b_Ca, 1);
        Vmath::Vlog(n, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Smul(n, 0.5*R*T/F, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Vsub(n, inarray[0], 1, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Smul(n, C_m*g_b_Ca, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_b_Ca, 1, outarray[0], 1);

        // L-Type Ca2+ current
        Array<OneD, NekDouble> &tmp_I_Ca_L = outarray[3];
        Vmath::Sadd(n, -65.0, inarray[0], 1, tmp_I_Ca_L, 1);
        Vmath::Vmul(n, inarray[11], 1, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Vmul(n, inarray[12], 1, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Vmul(n, inarray[13], 1, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Smul(n, C_m*g_Ca_L, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Ca_L, 1, outarray[0], 1);

        // Na-K Pump Current
        Array<OneD, NekDouble> &tmp_f_Na_k = outarray[12];
        Vmath::Smul(n, -F/R/T, inarray[0], 1, tmp_f_Na_k, 1);
        Vmath::Vexp(n, tmp_f_Na_k, 1, tmp, 1);
        Vmath::Smul(n, 0.0365*sigma, tmp, 1, tmp, 1);
        Vmath::Smul(n, -0.1*F/R/T, inarray[0], 1, tmp_f_Na_k, 1);
        Vmath::Vexp(n, tmp_f_Na_k, 1, tmp_f_Na_k, 1);
        Vmath::Smul(n, 0.1245, tmp_f_Na_k, 1, tmp_f_Na_k, 1);
        Vmath::Vadd(n, tmp_f_Na_k, 1, tmp, 1, tmp_f_Na_k, 1);
        Vmath::Sadd(n, 1.0, tmp_f_Na_k, 1, tmp_f_Na_k, 1);

        Array<OneD, NekDouble> &tmp_I_Na_K = outarray[13];
        Vmath::Sdiv(n, K_m_Na_i, inarray[17], 1, tmp_I_Na_K, 1);
        Vmath::Vpow(n, tmp_I_Na_K, 1, 1.5, tmp_I_Na_K, 1);
        Vmath::Sadd(n, 1.0, tmp_I_Na_K, 1, tmp_I_Na_K, 1);
        Vmath::Vmul(n, tmp_f_Na_k, 1, tmp_I_Na_K, 1, tmp_I_Na_K, 1);
        Vmath::Sdiv(n, C_m*I_Na_K_max*K_o/(K_o+K_i), tmp_I_Na_K, 1, tmp_I_Na_K, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Na_K, 1, outarray[0], 1);
        Vmath::Svtvp(n, -3.0, tmp_I_Na_K, 1, outarray[17], 1, outarray[17], 1);
        Vmath::Svtvp(n, 2.0, tmp_I_Na_K, 1, outarray[19], 1, outarray[19], 1);

        // Na-Ca exchanger current
        Array<OneD, NekDouble> &tmp_I_Na_Ca = outarray[4];
        Vmath::Smul(n, (gamma-1)*F/R/T, inarray[0], 1, tmp, 1);
        Vmath::Vexp(n, tmp, 1, tmp, 1);
        Vmath::Smul(n, K_sat, tmp, 1, tmp_I_Na_Ca, 1);
        Vmath::Sadd(n, 1.0, tmp_I_Na_Ca, 1, tmp_I_Na_Ca, 1);
        Vmath::Smul(n, (K_m_Na*K_m_Na*K_m_Na + Na_o*Na_o*Na_o)*(K_m_Ca + Ca_o), tmp_I_Na_Ca, 1, tmp_I_Na_Ca, 1);

        Vmath::Smul(n, Na_o*Na_o*Na_o, tmp, 1, tmp2, 1);
        Vmath::Vmul(n, tmp2, 1, inarray[18], 1, tmp2, 1);
        Vmath::Smul(n, gamma*F/R/T, inarray[0], 1, tmp, 1);
        Vmath::Vexp(n, tmp, 1, tmp, 1);
        Vmath::Vmul(n, inarray[17], 1, tmp, 1, tmp, 1);
        Vmath::Vmul(n, inarray[17], 1, tmp, 1, tmp, 1);
        Vmath::Vmul(n, inarray[17], 1, tmp, 1, tmp, 1);
        Vmath::Svtvm(n, Ca_o, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Smul(n, C_m*I_NaCa_max, tmp, 1, tmp, 1);
        Vmath::Vdiv(n, tmp, 1, tmp_I_Na_Ca, 1, tmp_I_Na_Ca, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Na_Ca, 1, outarray[0], 1);
        Vmath::Svtvp(n, -3.0, tmp_I_Na_Ca, 1, outarray[17], 1, outarray[17], 1);

        // Calcium Pump current
        Array<OneD, NekDouble> &tmp_I_p_Ca = outarray[5];
        Vmath::Sadd(n, 0.0005, inarray[18], 1, tmp_I_p_Ca, 1);
        Vmath::Vdiv(n, inarray[18], 1, tmp_I_p_Ca, 1, tmp_I_p_Ca, 1);
        Vmath::Smul(n, C_m*I_p_Ca_max, tmp_I_p_Ca, 1, tmp_I_p_Ca, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_p_Ca, 1, outarray[0], 1);

        ///@todo I_pK current
        ///@todo add I_stim to K_i current
        ///@todo I_ax

        // Scale currents by capacitance
        Vmath::Smul(n, 1.0/C_m, outarray[0], 1, outarray[0], 1);

        // Scale sodium and potassium by FV_i
        Vmath::Smul(n, 1.0/F/V_C, outarray[14], 1, outarray[14], 1);
        Vmath::Smul(n, 1.0/F/V_C, outarray[16], 1, outarray[16], 1);

        // I_tr
        Array<OneD, NekDouble> &tmp_I_tr = outarray[6];
        Vmath::Vsub(n, inarray[21], 1, inarray[20], 1, tmp_I_tr, 1);
        Vmath::Smul(n, 1.0/tau_tr, tmp_I_tr, 1, tmp_I_tr, 1);

        // I_leak (86)
        Array<OneD, NekDouble> &tmp_I_leak = outarray[5];
        Vmath::Vsub(n, inarray[17], 1, inarray[15], 1, tmp_I_leak, 1);
        Vmath::Smul(n, V_leak, tmp_I_leak, 1, tmp_I_leak, 1);

        // I_up
        Array<OneD, NekDouble> &tmp_I_up = outarray[6];
        Vmath::Vmul(n, inarray[15], 1, inarray[15], 1, tmp_I_up, 1);
/// @TODO HERE
        Vmath::Sdiv(n, NSR_K_up, inarray[18], 1, tmp_I_up, 1);
        Vmath::Sadd(n, 1.0, tmp_I_up, 1, tmp_I_up, 1);
        Vmath::Sdiv(n, NSR_I_up_max, tmp_I_up, 1, tmp_I_up, 1);

        // I_rel
        Array<OneD, NekDouble> &tmp_I_rel = outarray[7];
        Vmath::Vsub(n, inarray[20], 1, inarray[18], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[14], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[14], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[15], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[16], 1, tmp_I_rel, 1);
        Vmath::Smul(n, JSR_K_rel, tmp_I_rel, 1, tmp_I_rel, 1);

        // Ca_i current (94)
        Vmath::Svtvm(n, 2.0, tmp_I_Na_Ca, 1, tmp_I_p_Ca, 1, outarray[15], 1);
        Vmath::Vsub(n, outarray[15], 1, tmp_I_Ca_L, 1, outarray[15], 1);
        Vmath::Vsub(n, outarray[15], 1, tmp_I_b_Ca, 1, outarray[15], 1);
        Vmath::Smul(n, 0.5/F/V_C, outarray[15], 1, outarray[15], 1);
        Vmath::Vadd(n, tmp_I_leak, 1, outarray[15], 1, outarray[15], 1);
        Vmath::Svtvp(n, tmp_I_up, 1, outarray[15], 1, outarray[15], 1);
        Vmath::Svtvp(n, tmp_I_rel, 1, outarray[15], 1, outarray[15], 1);

        // Calcium concentration (18)
        Vmath::Vdiv(n, tmp_B1, 1, tmp_B2, 1, outarray[18], 1);

        // Calcium up (21)
        Vmath::Vsub(n, tmp_I_up, 1, tmp_I_up_leak, 1, outarray[21], 1);
        Vmath::Svtvp(n, -JSR_V_rel/JSR_V_up, tmp_I_tr, 1, outarray[21], 1, outarray[21], 1);

        // Calcium rel (20)
        Vmath::Vsub(n, tmp_I_tr, 1, tmp_I_rel, 1, tmp, 1);
        Vmath::Sadd(n, Km_Csqn, inarray[20], 1, outarray[20], 1);
        Vmath::Vmul(n, outarray[20], 1, outarray[20], 1, outarray[20], 1);
        Vmath::Sdiv(n, Csqn_max*Km_Csqn, outarray[20], 1, outarray[20], 1);
        Vmath::Sadd(n, 1.0, outarray[20], 1, outarray[20], 1);
        Vmath::Vdiv(n, tmp, 1, outarray[20], 1, outarray[20], 1);

        // Process gating variables
        const NekDouble * v;
        const NekDouble * x;
        NekDouble * x_new;
        // m
        for (i = 0, v = &inarray[0][0], x = &inarray[2][0], x_new = &outarray[2][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = (*v == (-47.13)) ? 3.2 : (0.32*(*v+47.13))/(1.0-exp((-0.1)*(*v + 47.13)));
            beta  = 0.08*exp(-(*v)/11.0);
            *x_new = alpha - (*x)*(alpha + beta);
        }
        // h
        for (i = 0, v = &inarray[0][0], x = &inarray[3][0], x_new = &outarray[3][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = (*v >= -40.0) ? 0.0 : 0.135*exp(-((*v)+80.0)/6.8);
            beta  = (*v >= -40.0) ? 1.0/(0.13*(1.0+exp(-(*v + 10.66)/11.1)))
                    : 3.56*exp(0.079*(*v))+310000.0*exp(0.35*(*v));
            *x_new = alpha - (*x)*(alpha + beta);
        }
        // j
        for (i = 0, v = &inarray[0][0], x = &inarray[4][0], x_new = &outarray[4][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = (*v >= -40.0) ? 0.0
                    : (-127140.0*exp(0.2444*(*v))-3.474e-05*exp(-0.04391*(*v)))*(((*v)+37.78)/(1.0+exp(0.311*((*v)+79.23))));
            beta  = (*v >= -40.0) ? (0.3*exp(-2.535e-07*(*v))/(1.0+exp(-0.1*(*v+32.0))))
                    : 0.1212*exp(-0.01052*(*v))/(1.0+exp(-0.1378*(*v+40.14)));
            *x_new = alpha - (*x)*(alpha + beta);
        }
        // oa
        for (i = 0, v = &inarray[0][0], x = &inarray[5][0], x_new = &outarray[5][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 0.65/(exp(-(*v+10.0)/8.5) + exp(-(*v-30.0)/59.0));
            beta  = 0.65/(2.5 + exp((*v+82.0)/17.0));
            *x_new = K_Q10*(alpha + beta)*(1.0/(1.0+exp(-(*v+20.47)/17.54)) - *x);
        }
        // oi
        for (i = 0, v = &inarray[0][0], x = &inarray[6][0], x_new = &outarray[6][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 1.0/(18.53 + exp((*v+113.7)/10.95));
            beta  = 1.0/(35.56 + exp(-(*v+1.26)/7.44));
            *x_new = K_Q10*(alpha + beta)*(1/(1+exp((*v+43.1)/5.3)) - *x);
        }
        // ua
        for (i = 0, v = &inarray[0][0], x = &inarray[7][0], x_new = &outarray[7][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 0.65/(exp(-(*v+10.0)/8.5)+exp(-(*v-30.0)/59.0));
            beta  = 0.65/(2.5+exp((*v+82.0)/17.0));
            *x_new = K_Q10*(alpha + beta)*(1.0/(1+exp(-(*v+30.3)/9.6)) - *x);
        }
        // ui
        for (i = 0, v = &inarray[0][0], x = &inarray[8][0], x_new = &outarray[8][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 1.0/(21.0 + exp(-(*v-185.0)/28.0));
            beta  = exp((*v-158.0)/16.0);
            *x_new = K_Q10*(alpha + beta)*(1.0/(1+exp((*v-99.45)/27.48)) - *x);
        }
        // xr
        for (i = 0, v = &inarray[0][0], x = &inarray[9][0], x_new = &outarray[9][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 0.0003*(*v+14.1)/(1-exp(-(*v+14.1)/5.0));
            beta  = 7.3898e-5*(*v-3.3328)/(exp((*v-3.3328)/5.1237)-1.0);
            *x_new = (alpha + beta)*(1.0/(1+exp(-(*v+14.1)/6.5)) - *x);
        }
        // xs
        for (i = 0, v = &inarray[0][0], x = &inarray[10][0], x_new = &outarray[10][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 4e-5*(*v-19.9)/(1.0-exp(-(*v-19.9)/17.0));
            beta  = 3.5e-5*(*v-19.9)/(exp((*v-19.9)/9.0)-1.0);
            *x_new = 2.0*(alpha + beta)*(1.0/sqrt(1.0+exp(-(*v-19.9)/12.7)) - *x);
        }
        // d
        for (i = 0, v = &inarray[0][0], x = &inarray[11][0], x_new = &outarray[11][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 1.0/(1.0 + exp(-(*v+10)/8.0));
            beta  = (1-exp(-(*v+10.0)/6.24))/(0.035*(*v+10.0)*(1+exp(-(*v+10.0)/6.24)));
            *x_new = (alpha - *x)/beta;
        }
        // f
        for (i = 0, v = &inarray[0][0], x = &inarray[12][0], x_new = &outarray[12][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            //alpha = 1.0/(1.0 + exp((*v+28.0)/6.9));
            alpha = exp((-(*v + 28.0)) / 6.9) / (1.0 + exp((-(*v + 28.0)) / 6.9));
            beta  = 9.0/(0.0197*exp(-0.0337*0.0337*(*v+10.0)*(*v+10.0))+0.02);
            *x_new = (alpha - *x)/beta;
        }
        // f_Ca
        for (i = 0, v = &inarray[0][0], x = &inarray[13][0], x_new = &outarray[13][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 1.0/(1.0+inarray[18][i]/0.00035);
            beta  = 2.0;
            *x_new = (alpha - *x)/beta;
        }

        Array<OneD, NekDouble> &tmp_Fn = outarray[16];
        Vmath::Svtsvtp(n, 0.5*5e-13/F, tmp_I_Ca_L, 1, -0.2*5e-13/F, tmp_I_Na_Ca, 1, tmp_Fn, 1);
        Vmath::Svtvm(n, 1e-12*JSR_V_rel, tmp_I_rel, 1, tmp_Fn, 1, tmp_Fn, 1);

        // u
        for (i = 0, v = &tmp_Fn[0], x = &inarray[14][0], x_new = &outarray[14][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 1.0/(1.0 + exp(-(*v - 3.4175e-13)/1.367e-15));
            beta  = 8.0;
            *x_new = (alpha - *x)/beta;
        }
        // v
        for (i = 0, v = &tmp_Fn[0], x = &inarray[15][0], x_new = &outarray[15][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 1.0 - 1.0/(1.0 + exp(-(*v - 6.835e-14)/13.67e-16));
            beta  = 1.91 + 2.09*(1.0+exp(-(*v - 3.4175e-13)/13.67e-16));
            *x_new = (alpha - *x)/beta;
        }
        // w
        for (i = 0, v = &inarray[0][0], x = &inarray[16][0], x_new = &outarray[16][0];
                i < n; ++i, ++v, ++x, ++x_new)
        {
            alpha = 1.0 - 1.0/(1.0 + exp(-(*v - 40.0)/17.0));
            beta  = 6.0*(1.0-exp(-(*v-7.9)/5.0))/(1.0+0.3*exp(-(*v-7.9)/5.0))/(*v-7.9);
            *x_new = (alpha - *x)/beta;
        }

    }


    /**
    *
    */
    void TenTusscherPanfilov04::v_PrintSummary(std::ostream &out)
    {
        out << "	Cell model      : TenTusscherPanfilov04" << std::endl;
    }

}
