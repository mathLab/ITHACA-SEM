///////////////////////////////////////////////////////////////////////////////
//
// File CourtemancheRamirezNattel.cpp
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
// Description: Courtemanche-Ramirez-Nattel ionic atrial cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <CardiacEPSolver/CellModels/CourtemancheRamirezNattel98.h>

using namespace std;

namespace Nektar
{
    std::string CourtemancheRamirezNattel98::className
              = GetCellModelFactory().RegisterCreatorFunction(
                        "CourtemancheRamirezNattel98",
                        CourtemancheRamirezNattel98::create,
                         "Ionic model of human atrial cell electrophysiology.");
    
    // Register cell model variants
    std::string CourtemancheRamirezNattel98::lookupIds[2] = {
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "Original", CourtemancheRamirezNattel98::eOriginal),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "AF", CourtemancheRamirezNattel98::eAF)
    };
    
    // Register default variant
    std::string CourtemancheRamirezNattel98::def =
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                    "CellModelVariant", "Original");

    /**
    *
    */
    CourtemancheRamirezNattel98::CourtemancheRamirezNattel98(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField)
            : CellModel(pSession, pField)
    {
        model_variant = pSession->GetSolverInfoAsEnum<
                CourtemancheRamirezNattel98::Variants>("CellModelVariant");

        C_m = 100;      // picoF
        g_Na = 7.8;     // nanoS_per_picoF
        g_K1 = 0.09;    // nanoS_per_picoF
        g_Kr = 0.029411765;
        g_Ks = 0.12941176;
        g_b_Na = 0.0006744375;
        g_b_Ca = 0.001131;
        R = 8.3143;
        T = 310.0;
        F = 96.4867;
        Na_o = 140.0;   // millimolar
        K_o = 5.4;      // millimolar
        sigma = 1.0/7.0*(exp(Na_o/67.3)-1);
        K_i = 1.5;
        K_m_Na_i = 10.0;
        I_Na_K_max = 0.59933874;
        I_NaCa_max = 1600.0;
        gamma = 0.35;
        Ca_o = 1.8;
        K_m_Na = 87.5;
        K_m_Ca = 1.38;
        K_sat = 0.1;
        I_p_Ca_max = 0.275;
        Trpn_max = 0.07;
        Km_Trpn = 0.0005;
        Cmdn_max = 0.05;
        Csqn_max = 10.0;
        Km_Cmdn = 0.00238;
        Km_Csqn = 0.8;
        NSR_I_up_max = 0.005;
        NSR_I_Ca_max = 15.0;
        NSR_K_up = 0.00092;
        JSR_K_rel = 30.0;
        JSR_V_cell = 20100.0;
        JSR_V_rel = 0.0048 * JSR_V_cell;
        JSR_V_up = 0.0552 * JSR_V_cell;
        tau_tr = 180.0;
        K_Q10 = 3.0;
        V_i = 0.68*JSR_V_cell;

        switch (model_variant) {
            case eOriginal:
                g_to = 0.1652;  // nanoS_per_picoF
                g_Kur_scaling = 1.0;
                g_Ca_L = 0.12375;
                break;
            case eAF:
                g_to = 0.0826;  // nanoS_per_picoF
                g_Kur_scaling = 0.5;
                g_Ca_L = 0.037125;
                break;
        }

        m_nvar = 21;

        // List gates and concentrations
        m_gates.push_back(1);
        m_gates.push_back(2);
        m_gates.push_back(3);
        m_gates.push_back(4);
        m_gates.push_back(5);
        m_gates.push_back(6);
        m_gates.push_back(7);
        m_gates.push_back(8);
        m_gates.push_back(9);
        m_gates.push_back(10);
        m_gates.push_back(11);
        m_gates.push_back(12);
        m_gates.push_back(13);
        m_gates.push_back(14);
        m_gates.push_back(15);
        m_concentrations.push_back(16);
        m_concentrations.push_back(17);
        m_concentrations.push_back(18);
        m_concentrations.push_back(19);
        m_concentrations.push_back(20);
    }
    
    
    
    /**
    *
    */
    CourtemancheRamirezNattel98::~CourtemancheRamirezNattel98()
    {
        
    }
    
    
    
    void CourtemancheRamirezNattel98::v_Update(
                     const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                           Array<OneD,        Array<OneD, NekDouble> >&outarray,
                     const NekDouble time)
    {
        ASSERTL0(inarray.get() != outarray.get(),
                "Must have different arrays for input and output.");

        // Variables
        //  0   V    membrane potential
        //  2   m    fast sodium current m gate
        //  3   h    fast sodium current h gate
        //  4   j    fast sodium current j gate
        //  5   o_a  transient outward potassium o_a gate
        //  6   o_i  transient outward potassium o_i gate
        //  7   u_a  ultra-rapid delayed rectifier K current gate
        //  8   u_i  ultra-rapid delayed rectifier K current gate
        //  9   x_r  rapid delayed rectifier K current gate
        //  10  x_s  slow delayed rectifier K current gate
        //  11  d    L_type calcium gate
        //  12  f    L-type calcium gate
        //  13  f_Ca L-type calcium gate
        //  14  u    Ca release u gate
        //  15  v    Ca release v gate
        //  16  w    Ca release w gate
        //  17  Na_i Sodium
        //  18  Ca_i Calcium
        //  19  K_i  Potassium
        //  20  Ca_rel Calcium Rel
        //  21  Ca_up  Calcium up
        int n = m_nq;
        int i = 0;
        NekDouble alpha, beta;
        Vmath::Zero(n, outarray[0], 1);

        Array<OneD, NekDouble> &tmp = outarray[11];
        Array<OneD, NekDouble> &tmp2 = outarray[12];

        // E_Na
        Array<OneD, NekDouble> &tmp_E_na = outarray[14];
        Vmath::Sdiv(n, Na_o, inarray[16], 1, tmp_E_na, 1);
        Vmath::Vlog(n, tmp_E_na, 1, tmp_E_na, 1);
        Vmath::Smul(n, R*T/F, tmp_E_na, 1, tmp_E_na, 1);

        // Sodium I_Na
        Array<OneD, NekDouble> &tmp_I_Na = outarray[15];
        Vmath::Vsub(n, inarray[0], 1, tmp_E_na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[1], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[1], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[1], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[2], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vmul(n, inarray[3], 1, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Smul(n, C_m*g_Na, tmp_I_Na, 1, tmp_I_Na, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Na, 1, outarray[0], 1);
        Vmath::Smul(n, -1.0, tmp_I_Na, 1, outarray[16], 1);

        // Background current, sodium
        Array<OneD, NekDouble> &tmp_I_b_Na = outarray[15];
        Vmath::Vsub(n, inarray[0], 1, tmp_E_na, 1, tmp_I_b_Na, 1);
        Vmath::Smul(n, C_m*g_b_Na, tmp_I_b_Na, 1, tmp_I_b_Na, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_b_Na, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[16], 1, tmp_I_b_Na, 1, outarray[16], 1);

        // V - E_K
        Array<OneD, NekDouble> &tmp_V_E_k = outarray[14];
        Vmath::Sdiv(n, K_o, inarray[18], 1, tmp_V_E_k, 1);
        Vmath::Vlog(n, tmp_V_E_k, 1, tmp_V_E_k, 1);
        Vmath::Smul(n, R*T/F, tmp_V_E_k, 1, tmp_V_E_k, 1);
        Vmath::Vsub(n, inarray[0], 1, tmp_V_E_k, 1, tmp_V_E_k, 1);

        // Potassium I_K1
        Array<OneD, NekDouble> &tmp_I_K1 = outarray[15];
        Vmath::Sadd(n, 80.0, inarray[0], 1, tmp_I_K1, 1);
        Vmath::Smul(n, 0.07, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Vexp(n, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Sadd(n, 1.0, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Vdiv(n, tmp_V_E_k, 1, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Smul(n, C_m*g_K1, tmp_I_K1, 1, tmp_I_K1, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_K1, 1, outarray[0], 1);
        Vmath::Smul(n, -1.0, tmp_I_K1, 1, outarray[18], 1);

        // Transient Outward K+ current
        Array<OneD, NekDouble> &tmp_I_to = outarray[15];
        Vmath::Vmul(n, inarray[5], 1, tmp_V_E_k, 1, tmp_I_to, 1);
        Vmath::Vmul(n, inarray[4], 1, tmp_I_to, 1, tmp_I_to, 1);
        Vmath::Vmul(n, inarray[4], 1, tmp_I_to, 1, tmp_I_to, 1);
        Vmath::Vmul(n, inarray[4], 1, tmp_I_to, 1, tmp_I_to, 1);
        Vmath::Smul(n, C_m*g_to, tmp_I_to, 1, tmp_I_to, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_to, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[18], 1, tmp_I_to, 1, outarray[18], 1);

        // Ultrarapid Delayed rectifier K+ current
        Array<OneD, NekDouble> &tmp_I_kur = outarray[15];
        Vmath::Sadd(n, -15.0, inarray[0], 1, tmp_I_kur, 1);
        Vmath::Smul(n, -1.0/13.0, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Vexp(n, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Sadd(n, 1.0, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Sdiv(n, 0.05, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Sadd(n, 0.005, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Vmul(n, tmp_V_E_k,  1, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Vmul(n, inarray[6], 1, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Vmul(n, inarray[6], 1, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Vmul(n, inarray[6], 1, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Vmul(n, inarray[7], 1, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Smul(n, C_m*g_Kur_scaling, tmp_I_kur, 1, tmp_I_kur, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_kur, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[18], 1, tmp_I_kur, 1, outarray[18], 1);

        // Rapid delayed outward rectifier K+ current
        Array<OneD, NekDouble> &tmp_I_Kr = outarray[15];
        Vmath::Sadd(n, 15.0, inarray[0], 1, tmp_I_Kr, 1);
        Vmath::Smul(n, 1.0/22.4, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vexp(n, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Sadd(n, 1.0, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vdiv(n, tmp_V_E_k, 1, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vmul(n, inarray[8], 1, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Smul(n, C_m*g_Kr, tmp_I_Kr, 1, tmp_I_Kr, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Kr, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[18], 1, tmp_I_Kr, 1, outarray[18], 1);

        // Slow delayed outward rectifier K+ Current
        Array<OneD, NekDouble> &tmp_I_Ks = outarray[15];
        Vmath::Vmul(n, inarray[9], 1, tmp_V_E_k, 1, tmp_I_Ks, 1);
        Vmath::Vmul(n, inarray[9], 1, tmp_I_Ks, 1, tmp_I_Ks, 1);
        Vmath::Smul(n, C_m*g_Ks, tmp_I_Ks, 1, tmp_I_Ks, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Ks, 1, outarray[0], 1);
        Vmath::Vsub(n, outarray[18], 1, tmp_I_Ks, 1, outarray[18], 1);

        // Background current, calcium
        Array<OneD, NekDouble> &tmp_I_b_Ca = outarray[1];
        Vmath::Sdiv(n, Ca_o, inarray[17], 1, tmp_I_b_Ca, 1);
        Vmath::Vlog(n, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Smul(n, 0.5*R*T/F, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Vsub(n, inarray[0], 1, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Smul(n, C_m*g_b_Ca, tmp_I_b_Ca, 1, tmp_I_b_Ca, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_b_Ca, 1, outarray[0], 1);

        // L-Type Ca2+ current
        Array<OneD, NekDouble> &tmp_I_Ca_L = outarray[2];
        Vmath::Sadd(n, -65.0, inarray[0], 1, tmp_I_Ca_L, 1);
        Vmath::Vmul(n, inarray[10], 1, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Vmul(n, inarray[11], 1, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Vmul(n, inarray[12], 1, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Smul(n, C_m*g_Ca_L, tmp_I_Ca_L, 1, tmp_I_Ca_L, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Ca_L, 1, outarray[0], 1);

        // Na-K Pump Current
        Array<OneD, NekDouble> &tmp_f_Na_k = outarray[14];
        Vmath::Smul(n, -F/R/T, inarray[0], 1, tmp_f_Na_k, 1);
        Vmath::Vexp(n, tmp_f_Na_k, 1, tmp, 1);
        Vmath::Smul(n, 0.0365*sigma, tmp, 1, tmp, 1);
        Vmath::Smul(n, -0.1*F/R/T, inarray[0], 1, tmp_f_Na_k, 1);
        Vmath::Vexp(n, tmp_f_Na_k, 1, tmp_f_Na_k, 1);
        Vmath::Smul(n, 0.1245, tmp_f_Na_k, 1, tmp_f_Na_k, 1);
        Vmath::Vadd(n, tmp_f_Na_k, 1, tmp, 1, tmp_f_Na_k, 1);
        Vmath::Sadd(n, 1.0, tmp_f_Na_k, 1, tmp_f_Na_k, 1);

        Array<OneD, NekDouble> &tmp_I_Na_K = outarray[15];
        Vmath::Sdiv(n, K_m_Na_i, inarray[16], 1, tmp_I_Na_K, 1);
        Vmath::Vpow(n, tmp_I_Na_K, 1, 1.5, tmp_I_Na_K, 1);
        Vmath::Sadd(n, 1.0, tmp_I_Na_K, 1, tmp_I_Na_K, 1);
        Vmath::Vmul(n, tmp_f_Na_k, 1, tmp_I_Na_K, 1, tmp_I_Na_K, 1);
        Vmath::Sdiv(n, C_m*I_Na_K_max*K_o/(K_o+K_i), tmp_I_Na_K, 1, tmp_I_Na_K, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Na_K, 1, outarray[0], 1);
        Vmath::Svtvp(n, -3.0, tmp_I_Na_K, 1, outarray[16], 1, outarray[16], 1);
        Vmath::Svtvp(n, 2.0, tmp_I_Na_K, 1, outarray[18], 1, outarray[18], 1);

        // Na-Ca exchanger current
        Array<OneD, NekDouble> &tmp_I_Na_Ca = outarray[3];
        Vmath::Smul(n, (gamma-1)*F/R/T, inarray[0], 1, tmp, 1);
        Vmath::Vexp(n, tmp, 1, tmp, 1);
        Vmath::Smul(n, K_sat, tmp, 1, tmp_I_Na_Ca, 1);
        Vmath::Sadd(n, 1.0, tmp_I_Na_Ca, 1, tmp_I_Na_Ca, 1);
        Vmath::Smul(n, (K_m_Na*K_m_Na*K_m_Na + Na_o*Na_o*Na_o)*(K_m_Ca + Ca_o), tmp_I_Na_Ca, 1, tmp_I_Na_Ca, 1);

        Vmath::Smul(n, Na_o*Na_o*Na_o, tmp, 1, tmp2, 1);
        Vmath::Vmul(n, tmp2, 1, inarray[17], 1, tmp2, 1);
        Vmath::Smul(n, gamma*F/R/T, inarray[0], 1, tmp, 1);
        Vmath::Vexp(n, tmp, 1, tmp, 1);
        Vmath::Vmul(n, inarray[16], 1, tmp, 1, tmp, 1);
        Vmath::Vmul(n, inarray[16], 1, tmp, 1, tmp, 1);
        Vmath::Vmul(n, inarray[16], 1, tmp, 1, tmp, 1);
        Vmath::Svtvm(n, Ca_o, tmp, 1, tmp2, 1, tmp, 1);
        Vmath::Smul(n, C_m*I_NaCa_max, tmp, 1, tmp, 1);
        Vmath::Vdiv(n, tmp, 1, tmp_I_Na_Ca, 1, tmp_I_Na_Ca, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_Na_Ca, 1, outarray[0], 1);
        Vmath::Svtvp(n, -3.0, tmp_I_Na_Ca, 1, outarray[16], 1, outarray[16], 1);

        // Calcium Pump current
        Array<OneD, NekDouble> &tmp_I_p_Ca = outarray[4];
        Vmath::Sadd(n, 0.0005, inarray[17], 1, tmp_I_p_Ca, 1);
        Vmath::Vdiv(n, inarray[17], 1, tmp_I_p_Ca, 1, tmp_I_p_Ca, 1);
        Vmath::Smul(n, C_m*I_p_Ca_max, tmp_I_p_Ca, 1, tmp_I_p_Ca, 1);
        Vmath::Vsub(n, outarray[0], 1, tmp_I_p_Ca, 1, outarray[0], 1);

        // Scale currents by capacitance
        Vmath::Smul(n, 1.0/C_m, outarray[0], 1, outarray[0], 1);

        // Scale sodium and potassium by FV_i
        Vmath::Smul(n, 1.0/F/V_i, outarray[16], 1, outarray[16], 1);
        Vmath::Smul(n, 1.0/F/V_i, outarray[18], 1, outarray[18], 1);

        // I_tr
        Array<OneD, NekDouble> &tmp_I_tr = outarray[5];
        Vmath::Vsub(n, inarray[20], 1, inarray[19], 1, tmp_I_tr, 1);
        Vmath::Smul(n, 1.0/tau_tr, tmp_I_tr, 1, tmp_I_tr, 1);

        // I_up_leak
        Array<OneD, NekDouble> &tmp_I_up_leak = outarray[6];
        Vmath::Smul(n, NSR_I_up_max/NSR_I_Ca_max, inarray[20], 1, tmp_I_up_leak, 1);

        // I_up
        Array<OneD, NekDouble> &tmp_I_up = outarray[7];
        Vmath::Sdiv(n, NSR_K_up, inarray[17], 1, tmp_I_up, 1);
        Vmath::Sadd(n, 1.0, tmp_I_up, 1, tmp_I_up, 1);
        Vmath::Sdiv(n, NSR_I_up_max, tmp_I_up, 1, tmp_I_up, 1);

        // I_rel
        Array<OneD, NekDouble> &tmp_I_rel = outarray[8];
        Vmath::Vsub(n, inarray[19], 1, inarray[17], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[13], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[13], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[14], 1, tmp_I_rel, 1);
        Vmath::Vmul(n, tmp_I_rel, 1, inarray[15], 1, tmp_I_rel, 1);
        Vmath::Smul(n, JSR_K_rel, tmp_I_rel, 1, tmp_I_rel, 1);

        // B1
        Array<OneD, NekDouble> &tmp_B1 = outarray[9];
        Vmath::Svtvm(n, 2.0, tmp_I_Na_Ca, 1, tmp_I_p_Ca, 1, tmp_B1, 1);
        Vmath::Vsub(n, tmp_B1, 1, tmp_I_Ca_L, 1, tmp_B1, 1);
        Vmath::Vsub(n, tmp_B1, 1, tmp_I_b_Ca, 1, tmp_B1, 1);
        Vmath::Smul(n, 0.5/F, tmp_B1, 1, tmp_B1, 1);
        Vmath::Svtvp(n, JSR_V_up, tmp_I_up_leak, 1, tmp_B1, 1, tmp_B1, 1);
        Vmath::Svtvp(n, -JSR_V_up, tmp_I_up, 1, tmp_B1, 1, tmp_B1, 1);
        Vmath::Svtvp(n, JSR_V_rel, tmp_I_rel, 1, tmp_B1, 1, tmp_B1, 1);
        Vmath::Smul(n, 1.0/V_i, tmp_B1, 1, tmp_B1, 1);

        // B2
        Array<OneD, NekDouble> &tmp_B2 = outarray[10];
        Vmath::Sadd(n, Km_Cmdn, inarray[17], 1, tmp_B2, 1);
        Vmath::Vmul(n, tmp_B2, 1, tmp_B2, 1, tmp_B2, 1);
        Vmath::Sdiv(n, Cmdn_max*Km_Cmdn, tmp_B2, 1, tmp_B2, 1);
        Vmath::Sadd(n, Km_Trpn, inarray[17], 1, tmp, 1);
        Vmath::Vmul(n, tmp, 1, tmp, 1, tmp, 1);
        Vmath::Sdiv(n, Trpn_max*Km_Trpn, tmp, 1, tmp, 1);
        Vmath::Vadd(n, tmp, 1, tmp_B2, 1, tmp_B2, 1);
        Vmath::Sadd(n, 1.0, tmp_B2, 1, tmp_B2, 1);

        // Calcium concentration (18)
        Vmath::Vdiv(n, tmp_B1, 1, tmp_B2, 1, outarray[17], 1);

        // Calcium up (21)
        Vmath::Vsub(n, tmp_I_up, 1, tmp_I_up_leak, 1, outarray[20], 1);
        Vmath::Svtvp(n, -JSR_V_rel/JSR_V_up, tmp_I_tr, 1, outarray[20], 1, outarray[20], 1);

        // Calcium rel (20)
        Vmath::Vsub(n, tmp_I_tr, 1, tmp_I_rel, 1, tmp, 1);
        Vmath::Sadd(n, Km_Csqn, inarray[19], 1, outarray[19], 1);
        Vmath::Vmul(n, outarray[19], 1, outarray[19], 1, outarray[19], 1);
        Vmath::Sdiv(n, Csqn_max*Km_Csqn, outarray[19], 1, outarray[19], 1);
        Vmath::Sadd(n, 1.0, outarray[19], 1, outarray[19], 1);
        Vmath::Vdiv(n, tmp, 1, outarray[19], 1, outarray[19], 1);

        // Process gating variables
        const NekDouble * v;
        const NekDouble * x;
        NekDouble * x_tau;
        NekDouble * x_new;
        // m
        for (i = 0, v = &inarray[0][0], x = &inarray[1][0], x_new = &outarray[1][0], x_tau = &m_gates_tau[0][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = (*v == (-47.13)) ? 3.2 : (0.32*(*v+47.13))/(1.0-exp((-0.1)*(*v + 47.13)));
            beta  = 0.08*exp(-(*v)/11.0);
            *x_tau = 1.0/(alpha + beta);
            *x_new = alpha*(*x_tau);
        }
        // h
        for (i = 0, v = &inarray[0][0], x = &inarray[2][0], x_new = &outarray[2][0], x_tau = &m_gates_tau[1][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = (*v >= -40.0) ? 0.0 : 0.135*exp(-((*v)+80.0)/6.8);
            beta  = (*v >= -40.0) ? 1.0/(0.13*(1.0+exp(-(*v + 10.66)/11.1)))
                    : 3.56*exp(0.079*(*v))+310000.0*exp(0.35*(*v));
            *x_tau = 1.0/(alpha + beta);
            *x_new = alpha*(*x_tau);
        }
        // j
        for (i = 0, v = &inarray[0][0], x = &inarray[3][0], x_new = &outarray[3][0], x_tau = &m_gates_tau[2][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = (*v >= -40.0) ? 0.0
                    : (-127140.0*exp(0.2444*(*v))-3.474e-05*exp(-0.04391*(*v)))*(((*v)+37.78)/(1.0+exp(0.311*((*v)+79.23))));
            beta  = (*v >= -40.0) ? (0.3*exp(-2.535e-07*(*v))/(1.0+exp(-0.1*(*v+32.0))))
                    : 0.1212*exp(-0.01052*(*v))/(1.0+exp(-0.1378*(*v+40.14)));
            *x_tau = 1.0/(alpha + beta);
            *x_new = alpha*(*x_tau);
        }
        // oa
        for (i = 0, v = &inarray[0][0], x = &inarray[4][0], x_new = &outarray[4][0], x_tau = &m_gates_tau[3][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = 0.65/(exp(-(*v+10.0)/8.5) + exp(-(*v-30.0)/59.0));
            beta  = 0.65/(2.5 + exp((*v+82.0)/17.0));
            *x_tau = 1.0/K_Q10/(alpha + beta);
            *x_new = (1.0/(1.0+exp(-(*v+20.47)/17.54)));
        }
        // oi
        for (i = 0, v = &inarray[0][0], x = &inarray[5][0], x_new = &outarray[5][0], x_tau = &m_gates_tau[4][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = 1.0/(18.53 + exp((*v+113.7)/10.95));
            beta  = 1.0/(35.56 + exp(-(*v+1.26)/7.44));
            *x_tau = 1.0/K_Q10/(alpha + beta);
            *x_new = (1.0/(1.0+exp((*v+43.1)/5.3)));
        }
        // ua
        for (i = 0, v = &inarray[0][0], x = &inarray[6][0], x_new = &outarray[6][0], x_tau = &m_gates_tau[5][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = 0.65/(exp(-(*v+10.0)/8.5)+exp(-(*v-30.0)/59.0));
            beta  = 0.65/(2.5+exp((*v+82.0)/17.0));
            *x_tau = 1.0/K_Q10/(alpha + beta);
            *x_new = 1.0/(1.0+exp(-(*v+30.3)/9.6));
        }
        // ui
        for (i = 0, v = &inarray[0][0], x = &inarray[7][0], x_new = &outarray[7][0], x_tau = &m_gates_tau[6][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = 1.0/(21.0 + exp(-(*v-185.0)/28.0));
            beta  = exp((*v-158.0)/16.0);
            *x_tau = 1.0/K_Q10/(alpha + beta);
            *x_new = 1.0/(1.0+exp((*v-99.45)/27.48));
        }
        // xr
        for (i = 0, v = &inarray[0][0], x = &inarray[8][0], x_new = &outarray[8][0], x_tau = &m_gates_tau[7][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = 0.0003*(*v+14.1)/(1-exp(-(*v+14.1)/5.0));
            beta  = 7.3898e-5*(*v-3.3328)/(exp((*v-3.3328)/5.1237)-1.0);
            *x_tau = 1.0/(alpha + beta);
            *x_new = 1.0/(1+exp(-(*v+14.1)/6.5));
        }
        // xs
        for (i = 0, v = &inarray[0][0], x = &inarray[9][0], x_new = &outarray[9][0], x_tau = &m_gates_tau[8][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            alpha = 4e-5*(*v-19.9)/(1.0-exp(-(*v-19.9)/17.0));
            beta  = 3.5e-5*(*v-19.9)/(exp((*v-19.9)/9.0)-1.0);
            *x_tau = 0.5/(alpha + beta);
            *x_new = 1.0/sqrt(1.0+exp(-(*v-19.9)/12.7));
        }
        // d
        for (i = 0, v = &inarray[0][0], x = &inarray[10][0], x_new = &outarray[10][0], x_tau = &m_gates_tau[9][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            *x_tau  = (1-exp(-(*v+10.0)/6.24))/(0.035*(*v+10.0)*(1+exp(-(*v+10.0)/6.24)));
            *x_new = 1.0/(1.0 + exp(-(*v+10)/8.0));
        }
        // f
        for (i = 0, v = &inarray[0][0], x = &inarray[11][0], x_new = &outarray[11][0], x_tau = &m_gates_tau[10][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            //alpha = 1.0/(1.0 + exp((*v+28.0)/6.9));
            *x_tau  = 9.0/(0.0197*exp(-0.0337*0.0337*(*v+10.0)*(*v+10.0))+0.02);
            *x_new = exp((-(*v + 28.0)) / 6.9) / (1.0 + exp((-(*v + 28.0)) / 6.9));
        }
        // f_Ca
        for (i = 0, v = &inarray[0][0], x = &inarray[12][0], x_new = &outarray[12][0], x_tau = &m_gates_tau[11][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            *x_tau  = 2.0;
            *x_new = 1.0/(1.0+inarray[17][i]/0.00035);
        }

        Array<OneD, NekDouble> &tmp_Fn = outarray[15];
        Vmath::Svtsvtp(n, 0.5*5e-13/F, tmp_I_Ca_L, 1, -0.2*5e-13/F, tmp_I_Na_Ca, 1, tmp_Fn, 1);
        Vmath::Svtvm(n, 1e-12*JSR_V_rel, tmp_I_rel, 1, tmp_Fn, 1, tmp_Fn, 1);

        // u
        for (i = 0, v = &tmp_Fn[0], x = &inarray[13][0], x_new = &outarray[13][0], x_tau = &m_gates_tau[12][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            *x_tau  = 8.0;
            *x_new = 1.0/(1.0 + exp(-(*v - 3.4175e-13)/1.367e-15));
        }
        // v
        for (i = 0, v = &tmp_Fn[0], x = &inarray[14][0], x_new = &outarray[14][0], x_tau = &m_gates_tau[13][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            *x_tau  = 1.91 + 2.09/(1.0+exp(-(*v - 3.4175e-13)/13.67e-16));
            *x_new = 1.0 - 1.0/(1.0 + exp(-(*v - 6.835e-14)/13.67e-16));
        }
        // w
        for (i = 0, v = &inarray[0][0], x = &inarray[15][0], x_new = &outarray[15][0], x_tau = &m_gates_tau[14][0];
                i < n; ++i, ++v, ++x, ++x_new, ++x_tau)
        {
            *x_tau  = 6.0*(1.0-exp(-(*v-7.9)/5.0))/(1.0+0.3*exp(-(*v-7.9)/5.0))/(*v-7.9);
            *x_new = 1.0 - 1.0/(1.0 + exp(-(*v - 40.0)/17.0));
        }

    }


    /**
    *
    */
    void CourtemancheRamirezNattel98::v_GenerateSummary(SummaryList& s)
    {
        SolverUtils::AddSummaryItem(s, "Cell model","CourtemancheRamirezNattel98");
        SolverUtils::AddSummaryItem(s, "Cell model var.", lookupIds[model_variant]);
    }


    void CourtemancheRamirezNattel98::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, -81.0,      m_cellSol[0],  1);
        Vmath::Fill(m_nq, 2.908e-03,  m_cellSol[1],  1);
        Vmath::Fill(m_nq, 9.649e-01,  m_cellSol[2],  1);
        Vmath::Fill(m_nq, 9.775e-01,  m_cellSol[3],  1);
        Vmath::Fill(m_nq, 3.043e-02,  m_cellSol[4],  1);
        Vmath::Fill(m_nq, 9.992e-01,  m_cellSol[5],  1);
        Vmath::Fill(m_nq, 4.966e-03,  m_cellSol[6],  1);
        Vmath::Fill(m_nq, 9.986e-01,  m_cellSol[7],  1);
        Vmath::Fill(m_nq, 3.296e-05,  m_cellSol[8],  1);
        Vmath::Fill(m_nq, 1.869e-02,  m_cellSol[9],  1);
        Vmath::Fill(m_nq, 1.367e-04,  m_cellSol[10], 1);
        Vmath::Fill(m_nq, 9.996e-01,  m_cellSol[11], 1);
        Vmath::Fill(m_nq, 7.755e-01,  m_cellSol[12], 1);
        Vmath::Fill(m_nq, 2.35e-112,  m_cellSol[13], 1);
        Vmath::Fill(m_nq, 1.0,        m_cellSol[14], 1);
        Vmath::Fill(m_nq, 0.9992,     m_cellSol[15], 1);
        Vmath::Fill(m_nq, 1.117e+01,  m_cellSol[16], 1);
        Vmath::Fill(m_nq, 1.013e-04,  m_cellSol[17], 1);
        Vmath::Fill(m_nq, 1.39e+02,   m_cellSol[18], 1);
        Vmath::Fill(m_nq, 1.488,      m_cellSol[19], 1);
        Vmath::Fill(m_nq, 1.488,      m_cellSol[20], 1);
    }

    std::string CourtemancheRamirezNattel98::v_GetCellVarName(unsigned int idx)
    {
        switch (idx)
        {
            case 0:  return "u";
            case 1:  return "m";
            case 2:  return "h";
            case 3:  return "j";
            case 4:  return "o_a";
            case 5:  return "o_i";
            case 6:  return "u_a";
            case 7:  return "u_i";
            case 8:  return "x_r";
            case 9: return "x_s";
            case 10: return "d";
            case 11: return "f";
            case 12: return "f_Ca";
            case 13: return "U";
            case 14: return "V";
            case 15: return "W";
            case 16: return "Na_i";
            case 17: return "Ca_i";
            case 18: return "K_i";
            case 19: return "Ca_rel";
            case 20: return "Ca_up";
            default: return "unknown";
        }
    }

}
