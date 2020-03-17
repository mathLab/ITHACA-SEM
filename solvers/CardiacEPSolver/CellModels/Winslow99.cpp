///////////////////////////////////////////////////////////////////////////////
//
// File Winslow99.cpp
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
// Description: Winslow 1999 cell model
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
//#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <CardiacEPSolver/CellModels/Winslow99.h>

namespace Nektar
{
    std::string Winslow99::className
              = GetCellModelFactory().RegisterCreatorFunction(
                        "Winslow99",
                        Winslow99::create,
                         "Winslow 1999 cell model.");
    
    
    /**
    *
    */
    Winslow99::Winslow99(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const MultiRegions::ExpListSharedPtr& pField):
        CellModel(pSession, pField)
    {
        m_nvar = 33;
        m_gates.push_back(1);
        m_gates.push_back(2);
        m_gates.push_back(3);
        m_gates.push_back(4);
        m_gates.push_back(5);
        m_gates.push_back(6);
        m_gates.push_back(7);
        m_concentrations.push_back(8);
        m_concentrations.push_back(9);
        m_concentrations.push_back(10);
        m_concentrations.push_back(11);
        m_concentrations.push_back(12);
        m_concentrations.push_back(13);
        m_concentrations.push_back(14);
        m_concentrations.push_back(15);
        m_concentrations.push_back(16);
        m_concentrations.push_back(17);
        m_concentrations.push_back(18);
        m_concentrations.push_back(19);
        m_gates.push_back(20);
        m_concentrations.push_back(21);
        m_concentrations.push_back(22);
        m_concentrations.push_back(23);
        m_concentrations.push_back(24);
        m_concentrations.push_back(25);
        m_concentrations.push_back(26);
        m_concentrations.push_back(27);
        m_concentrations.push_back(28);
        m_concentrations.push_back(29);
        m_concentrations.push_back(30);
        m_concentrations.push_back(31);
        m_concentrations.push_back(32);
    }
    
    
    
    
    
    void Winslow99::v_Update(
                     const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                           Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                                           const NekDouble time)
    {
        int nq = m_nq;
        for (unsigned int i = 0; i < nq; ++i)
        {
            
            // Inputs:
            // Time units: millisecond
            NekDouble var_chaste_interface__membrane__V = inarray[0][i];
            // Units: millivolt; Initial value: -96.1638
            NekDouble var_chaste_interface__fast_sodium_current_m_gate__m = inarray[1][i];
            // Units: dimensionless; Initial value: 0.0328302
            NekDouble var_chaste_interface__fast_sodium_current_h_gate__h = inarray[2][i];
            // Units: dimensionless; Initial value: 0.988354
            NekDouble var_chaste_interface__fast_sodium_current_j_gate__j = inarray[3][i];
            // Units: dimensionless; Initial value: 0.99254
            NekDouble var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = inarray[4][i];
            // Units: dimensionless; Initial value: 0.51
            NekDouble var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = inarray[5][i];
            // Units: dimensionless; Initial value: 0.264
            NekDouble var_chaste_interface__transient_outward_potassium_current_X_to1_gate__X_to1 = inarray[6][i];
            // Units: dimensionless; Initial value: 2.63
            NekDouble var_chaste_interface__transient_outward_potassium_current_Y_to1_gate__Y_to1 = inarray[7][i];
            // Units: dimensionless; Initial value: 0.99
            NekDouble var_chaste_interface__L_type_Ca_current__O = inarray[8][i];
            // Units: dimensionless; Initial value: 9.84546e-21
            NekDouble var_chaste_interface__L_type_Ca_current__O_Ca = inarray[9][i];
            // Units: dimensionless; Initial value: 0
            NekDouble var_chaste_interface__L_type_Ca_current__C0 = inarray[10][i];
            // Units: dimensionless; Initial value: 0.997208
            NekDouble var_chaste_interface__L_type_Ca_current__C1 = inarray[11][i];
            // Units: dimensionless; Initial value: 6.38897e-5
            NekDouble var_chaste_interface__L_type_Ca_current__C2 = inarray[12][i];
            // Units: dimensionless; Initial value: 1.535e-9
            NekDouble var_chaste_interface__L_type_Ca_current__C3 = inarray[13][i];
            // Units: dimensionless; Initial value: 1.63909e-14
            NekDouble var_chaste_interface__L_type_Ca_current__C4 = inarray[14][i];
            // Units: dimensionless; Initial value: 6.56337e-20
            NekDouble var_chaste_interface__L_type_Ca_current__C_Ca0 = inarray[15][i];
            // Units: dimensionless; Initial value: 0.00272826
            NekDouble var_chaste_interface__L_type_Ca_current__C_Ca1 = inarray[16][i];
            // Units: dimensionless; Initial value: 6.99215e-7
            NekDouble var_chaste_interface__L_type_Ca_current__C_Ca2 = inarray[17][i];
            // Units: dimensionless; Initial value: 6.71989e-11
            NekDouble var_chaste_interface__L_type_Ca_current__C_Ca3 = inarray[18][i];
            // Units: dimensionless; Initial value: 2.87031e-15
            NekDouble var_chaste_interface__L_type_Ca_current__C_Ca4 = inarray[19][i];
            // Units: dimensionless; Initial value: 4.59752e-20
            NekDouble var_chaste_interface__L_type_Ca_current_y_gate__y = inarray[20][i];
            // Units: dimensionless; Initial value: 0.798
            NekDouble var_chaste_interface__RyR_channel__P_O1 = inarray[21][i];
            // Units: dimensionless; Initial value: 0
            NekDouble var_chaste_interface__RyR_channel__P_O2 = inarray[22][i];
            // Units: dimensionless; Initial value: 0
            NekDouble var_chaste_interface__RyR_channel__P_C1 = inarray[23][i];
            // Units: dimensionless; Initial value: 0.47
            NekDouble var_chaste_interface__RyR_channel__P_C2 = inarray[24][i];
            // Units: dimensionless; Initial value: 0.53
            NekDouble var_chaste_interface__intracellular_Ca_fluxes__HTRPNCa = inarray[25][i];
            // Units: millimolar; Initial value: 0.98
            NekDouble var_chaste_interface__intracellular_Ca_fluxes__LTRPNCa = inarray[26][i];
            // Units: millimolar; Initial value: 0.078
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Nai = inarray[27][i];
            // Units: millimolar; Initial value: 10
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Cai = inarray[28][i];
            // Units: millimolar; Initial value: 0.00008
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ki = inarray[29][i];
            // Units: millimolar; Initial value: 157.8
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_ss = inarray[30][i];
            // Units: millimolar; Initial value: 0.00011
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_JSR = inarray[31][i];
            // Units: millimolar; Initial value: 0.257
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_NSR = inarray[32][i];
            // Units: millimolar; Initial value: 0.257


            // Mathematics
            NekDouble d_dt_chaste_interface__membrane__V;
            const NekDouble var_membrane__R = 8.314472; // joule_per_mole_kelvin
            const NekDouble var_membrane__T = 310.0; // kelvin
            const NekDouble var_membrane__F = 96.4853415; // coulomb_per_millimole
            const NekDouble var_fast_sodium_current__j = var_chaste_interface__fast_sodium_current_j_gate__j; // dimensionless
            const NekDouble var_fast_sodium_current__h = var_chaste_interface__fast_sodium_current_h_gate__h; // dimensionless
            const NekDouble var_fast_sodium_current__g_Na = 12.8; // milliS_per_microF
            const NekDouble var_fast_sodium_current__m = var_chaste_interface__fast_sodium_current_m_gate__m; // dimensionless
            const NekDouble var_fast_sodium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_fast_sodium_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_fast_sodium_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_standard_ionic_concentrations__Nao = 138.0; // millimolar
            const NekDouble var_fast_sodium_current__Nao = var_standard_ionic_concentrations__Nao; // millimolar
            const NekDouble var_fast_sodium_current__Nai = var_chaste_interface__intracellular_ion_concentrations__Nai; // millimolar
            const NekDouble var_fast_sodium_current__T = var_membrane__T; // kelvin
            const NekDouble var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Nao / var_fast_sodium_current__Nai); // millivolt
            const NekDouble var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na); // microA_per_microF
            const NekDouble var_L_type_Ca_current__O = var_chaste_interface__L_type_Ca_current__O; // dimensionless
            const NekDouble var_L_type_Ca_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_L_type_Ca_current__P_Ca = 0.0003125; // cm_per_second
            const NekDouble var_standard_ionic_concentrations__Cao = 2.0; // millimolar
            const NekDouble var_L_type_Ca_current__Cao = var_standard_ionic_concentrations__Cao; // millimolar
            const NekDouble var_L_type_Ca_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_L_type_Ca_current__T = var_membrane__T; // kelvin
            const NekDouble var_L_type_Ca_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_L_type_Ca_current__i_Ca_max = ((((var_L_type_Ca_current__P_Ca / (1.0 * 1.0)) * 4.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0) * 1000.0) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((0.001 * exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - (0.341 * var_L_type_Ca_current__Cao))) / (exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0); // microA_per_microF
            const NekDouble var_L_type_Ca_current__y = var_chaste_interface__L_type_Ca_current_y_gate__y; // dimensionless
            const NekDouble var_L_type_Ca_current__O_Ca = var_chaste_interface__L_type_Ca_current__O_Ca; // dimensionless
            const NekDouble var_L_type_Ca_current__i_Ca = var_L_type_Ca_current__i_Ca_max * var_L_type_Ca_current__y * (var_L_type_Ca_current__O + var_L_type_Ca_current__O_Ca); // microA_per_microF
            const NekDouble var_L_type_Ca_current__P_K = 5.79e-07; // cm_per_second
            const NekDouble var_L_type_Ca_current__i_Ca_half =  -0.265; // microA_per_microF
            const NekDouble var_L_type_Ca_current__p_prime_k = var_L_type_Ca_current__P_K / (1.0 + (var_L_type_Ca_current__i_Ca_max / var_L_type_Ca_current__i_Ca_half)); // cm_per_second
            const NekDouble var_standard_ionic_concentrations__Ko = 4.0; // millimolar
            const NekDouble var_L_type_Ca_current__Ko = var_standard_ionic_concentrations__Ko; // millimolar
            const NekDouble var_L_type_Ca_current__Ki = var_chaste_interface__intracellular_ion_concentrations__Ki; // millimolar
            const NekDouble var_L_type_Ca_current__i_Ca_K = ((((var_L_type_Ca_current__p_prime_k / (1.0 * 1.0)) * var_L_type_Ca_current__y * (var_L_type_Ca_current__O + var_L_type_Ca_current__O_Ca) * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__Ki * exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__Ko)) / (exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0); // microA_per_microF
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__Ko = var_standard_ionic_concentrations__Ko; // millimolar
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__f_Ko = sqrt(var_rapid_activating_delayed_rectifiyer_K_current__Ko / 4.0); // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__Ki = var_chaste_interface__intracellular_ion_concentrations__Ki; // millimolar
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__T = var_membrane__T; // kelvin
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__E_K = ((var_rapid_activating_delayed_rectifiyer_K_current__R * var_rapid_activating_delayed_rectifiyer_K_current__T) / var_rapid_activating_delayed_rectifiyer_K_current__F) * log(var_rapid_activating_delayed_rectifiyer_K_current__Ko / var_rapid_activating_delayed_rectifiyer_K_current__Ki); // millivolt
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__R_V = 1.0 / (1.0 + (1.4945 * exp(0.0446 * var_rapid_activating_delayed_rectifiyer_K_current__V))); // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__X_kr = var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr; // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__g_Kr = 0.0034; // milliS_per_microF
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__g_Kr * var_rapid_activating_delayed_rectifiyer_K_current__f_Ko * var_rapid_activating_delayed_rectifiyer_K_current__R_V * var_rapid_activating_delayed_rectifiyer_K_current__X_kr * (var_rapid_activating_delayed_rectifiyer_K_current__V - var_rapid_activating_delayed_rectifiyer_K_current__E_K); // microA_per_microF
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__g_Ks = 0.0027134; // milliS_per_microF
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__Ko = var_standard_ionic_concentrations__Ko; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__Nao = var_standard_ionic_concentrations__Nao; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__Ki = var_chaste_interface__intracellular_ion_concentrations__Ki; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__Nai = var_chaste_interface__intracellular_ion_concentrations__Nai; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__T = var_membrane__T; // kelvin
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__E_Ks = ((var_slow_activating_delayed_rectifiyer_K_current__R * var_slow_activating_delayed_rectifiyer_K_current__T) / var_slow_activating_delayed_rectifiyer_K_current__F) * log((var_slow_activating_delayed_rectifiyer_K_current__Ko + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Nao)) / (var_slow_activating_delayed_rectifiyer_K_current__Ki + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Nai))); // millivolt
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__X_ks = var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks; // dimensionless
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__g_Ks * pow(var_slow_activating_delayed_rectifiyer_K_current__X_ks, 2.0) * (var_slow_activating_delayed_rectifiyer_K_current__V - var_slow_activating_delayed_rectifiyer_K_current__E_Ks); // microA_per_microF
            const NekDouble var_transient_outward_potassium_current__X_to1 = var_chaste_interface__transient_outward_potassium_current_X_to1_gate__X_to1; // dimensionless
            const NekDouble var_transient_outward_potassium_current__Y_to1 = var_chaste_interface__transient_outward_potassium_current_Y_to1_gate__Y_to1; // dimensionless
            const NekDouble var_transient_outward_potassium_current__g_to1 = 0.23815; // milliS_per_microF
            const NekDouble var_transient_outward_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K; // millivolt
            const NekDouble var_transient_outward_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_transient_outward_potassium_current__i_to1 = var_transient_outward_potassium_current__g_to1 * var_transient_outward_potassium_current__X_to1 * var_transient_outward_potassium_current__Y_to1 * (var_transient_outward_potassium_current__V - var_transient_outward_potassium_current__E_K); // microA_per_microF
            const NekDouble var_time_independent_potassium_current__Ko = var_standard_ionic_concentrations__Ko; // millimolar
            const NekDouble var_time_independent_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K; // millivolt
            const NekDouble var_time_independent_potassium_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_time_independent_potassium_current_K1_gate__F = var_time_independent_potassium_current__F; // coulomb_per_millimole
            const NekDouble var_time_independent_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V; // millivolt
            const NekDouble var_time_independent_potassium_current__T = var_membrane__T; // kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__T = var_time_independent_potassium_current__T; // kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__E_K = var_time_independent_potassium_current__E_K; // millivolt
            const NekDouble var_time_independent_potassium_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__R = var_time_independent_potassium_current__R; // joule_per_mole_kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__K1_infinity_V = 1.0 / (2.0 + exp(((1.5 * var_time_independent_potassium_current_K1_gate__F) / (var_time_independent_potassium_current_K1_gate__R * var_time_independent_potassium_current_K1_gate__T)) * (var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K))); // dimensionless
            const NekDouble var_time_independent_potassium_current__K1_infinity_V = var_time_independent_potassium_current_K1_gate__K1_infinity_V; // dimensionless
            const NekDouble var_time_independent_potassium_current__g_K1 = 2.8; // milliS_per_microF
            const NekDouble var_time_independent_potassium_current__K_mK1 = 13.0; // millimolar
            const NekDouble var_time_independent_potassium_current__i_K1 = ((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity_V * var_time_independent_potassium_current__Ko) / (var_time_independent_potassium_current__Ko + var_time_independent_potassium_current__K_mK1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K); // microA_per_microF
            const NekDouble var_plateau_potassium_current__g_Kp = 0.002216; // milliS_per_microF
            const NekDouble var_plateau_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_plateau_potassium_current_Kp_gate__V = var_plateau_potassium_current__V; // millivolt
            const NekDouble var_plateau_potassium_current_Kp_gate__Kp_V = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current_Kp_gate__V) / 5.98)); // dimensionless
            const NekDouble var_plateau_potassium_current__Kp_V = var_plateau_potassium_current_Kp_gate__Kp_V; // dimensionless
            const NekDouble var_plateau_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K; // millivolt
            const NekDouble var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp_V * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_K); // microA_per_microF
            const NekDouble var_Na_Ca_exchanger__Nao = var_standard_ionic_concentrations__Nao; // millimolar
            const NekDouble var_Na_Ca_exchanger__K_sat = 0.2; // dimensionless
            const NekDouble var_Na_Ca_exchanger__K_mNa = 87.5; // millimolar
            const NekDouble var_Na_Ca_exchanger__Nai = var_chaste_interface__intracellular_ion_concentrations__Nai; // millimolar
            const NekDouble var_Na_Ca_exchanger__K_NaCa = 0.3; // microA_per_microF
            const NekDouble var_Na_Ca_exchanger__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_Na_Ca_exchanger__T = var_membrane__T; // kelvin
            const NekDouble var_Na_Ca_exchanger__Cao = var_standard_ionic_concentrations__Cao; // millimolar
            const NekDouble var_Na_Ca_exchanger__eta = 0.35; // dimensionless
            const NekDouble var_Na_Ca_exchanger__K_mCa = 1.38; // millimolar
            const NekDouble var_Na_Ca_exchanger__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_Na_Ca_exchanger__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_Na_Ca_exchanger__Cai = var_chaste_interface__intracellular_ion_concentrations__Cai; // millimolar
            const NekDouble var_Na_Ca_exchanger__i_NaCa = ((var_Na_Ca_exchanger__K_NaCa * 5000.0) / ((pow(var_Na_Ca_exchanger__K_mNa, 3.0) + pow(var_Na_Ca_exchanger__Nao, 3.0)) * (var_Na_Ca_exchanger__K_mCa + var_Na_Ca_exchanger__Cao) * (1.0 + (var_Na_Ca_exchanger__K_sat * exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)))))) * ((exp((var_Na_Ca_exchanger__eta * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Nai, 3.0) * var_Na_Ca_exchanger__Cao) - (exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Nao, 3.0) * var_Na_Ca_exchanger__Cai)); // microA_per_microF
            const NekDouble var_sodium_potassium_pump__Ko = var_standard_ionic_concentrations__Ko; // millimolar
            const NekDouble var_sodium_potassium_pump__K_mNai = 10.0; // millimolar
            const NekDouble var_sodium_potassium_pump__Nai = var_chaste_interface__intracellular_ion_concentrations__Nai; // millimolar
            const NekDouble var_sodium_potassium_pump__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_sodium_potassium_pump__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_sodium_potassium_pump__T = var_membrane__T; // kelvin
            const NekDouble var_sodium_potassium_pump__Nao = var_standard_ionic_concentrations__Nao; // millimolar
            const NekDouble var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Nao / 67.3) - 1.0); // dimensionless
            const NekDouble var_sodium_potassium_pump__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_sodium_potassium_pump__f_NaK = 1.0 / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__V * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__V) * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T)))); // dimensionless
            const NekDouble var_sodium_potassium_pump__I_NaK = 0.693; // microA_per_microF
            const NekDouble var_sodium_potassium_pump__K_mKo = 1.5; // millimolar
            const NekDouble var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__I_NaK * var_sodium_potassium_pump__f_NaK) / (1.0 + pow(var_sodium_potassium_pump__K_mNai / var_sodium_potassium_pump__Nai, 1.5))) * var_sodium_potassium_pump__Ko) / (var_sodium_potassium_pump__Ko + var_sodium_potassium_pump__K_mKo); // microA_per_microF
            const NekDouble var_sarcolemmal_calcium_pump__I_pCa = 0.05; // microA_per_microF
            const NekDouble var_sarcolemmal_calcium_pump__Cai = var_chaste_interface__intracellular_ion_concentrations__Cai; // millimolar
            const NekDouble var_sarcolemmal_calcium_pump__K_mpCa = 5e-05; // millimolar
            const NekDouble var_sarcolemmal_calcium_pump__i_p_Ca = (var_sarcolemmal_calcium_pump__I_pCa * var_sarcolemmal_calcium_pump__Cai) / (var_sarcolemmal_calcium_pump__K_mpCa + var_sarcolemmal_calcium_pump__Cai); // microA_per_microF
            const NekDouble var_calcium_background_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_calcium_background_current__Cai = var_chaste_interface__intracellular_ion_concentrations__Cai; // millimolar
            const NekDouble var_calcium_background_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_calcium_background_current__T = var_membrane__T; // kelvin
            const NekDouble var_calcium_background_current__Cao = var_standard_ionic_concentrations__Cao; // millimolar
            const NekDouble var_calcium_background_current__E_Ca = ((var_calcium_background_current__R * var_calcium_background_current__T) / (2.0 * var_calcium_background_current__F)) * log(var_calcium_background_current__Cao / var_calcium_background_current__Cai); // millivolt
            const NekDouble var_calcium_background_current__g_Cab = 0.0003842; // milliS_per_microF
            const NekDouble var_calcium_background_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_calcium_background_current__i_Ca_b = var_calcium_background_current__g_Cab * (var_calcium_background_current__V - var_calcium_background_current__E_Ca); // microA_per_microF
            const NekDouble var_sodium_background_current__g_Nab = 0.0031; // milliS_per_microF
            const NekDouble var_sodium_background_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_sodium_background_current__E_Na = var_fast_sodium_current__E_Na; // millivolt
            const NekDouble var_sodium_background_current__i_Na_b = var_sodium_background_current__g_Nab * (var_sodium_background_current__V - var_sodium_background_current__E_Na); // microA_per_microF
            const NekDouble var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_m_gate__beta_m = 80.0 * exp((-var_fast_sodium_current_m_gate__V) / 11.0); // per_second
            const NekDouble var_fast_sodium_current_m_gate__E0_m = var_fast_sodium_current_m_gate__V + 47.13; // millivolt
            const NekDouble var_fast_sodium_current_m_gate__alpha_m = (fabs(var_fast_sodium_current_m_gate__E0_m) < 1e-05) ? (1000.0 / (0.1 - (0.005 * var_fast_sodium_current_m_gate__E0_m))) : ((320.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp((-0.1) * var_fast_sodium_current_m_gate__E0_m))); // per_second
            const NekDouble var_fast_sodium_current_m_gate__m = var_fast_sodium_current__m; // dimensionless
            const NekDouble var_fast_sodium_current_m_gate__d_m_d_environment__time = (var_fast_sodium_current_m_gate__V >= (-90.0)) ? ((var_fast_sodium_current_m_gate__alpha_m * (1.0 - var_fast_sodium_current_m_gate__m)) - (var_fast_sodium_current_m_gate__beta_m * var_fast_sodium_current_m_gate__m)) : 0.0; // per_second
            const NekDouble var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time = var_fast_sodium_current_m_gate__d_m_d_environment__time; // per_second
            const NekDouble var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((3560.0 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.35 * var_fast_sodium_current_h_gate__V))) : (1000.0 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1))))); // per_second
            const NekDouble var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (135.0 * exp((80.0 + var_fast_sodium_current_h_gate__V) / (-6.8))) : 0.0; // per_second
            const NekDouble var_fast_sodium_current_h_gate__h = var_fast_sodium_current__h; // dimensionless
            const NekDouble var_fast_sodium_current_h_gate__d_h_d_environment__time = (var_fast_sodium_current_h_gate__alpha_h * (1.0 - var_fast_sodium_current_h_gate__h)) - (var_fast_sodium_current_h_gate__beta_h * var_fast_sodium_current_h_gate__h); // per_second
            const NekDouble var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time = var_fast_sodium_current_h_gate__d_h_d_environment__time; // per_second
            const NekDouble var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((1000.0 * (-((127140.0 * exp(0.2444 * var_fast_sodium_current_j_gate__V)) + (3.474e-05 * exp((-0.04391) * var_fast_sodium_current_j_gate__V)))) * (var_fast_sodium_current_j_gate__V + 37.78)) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0; // per_second
            const NekDouble var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((121.2 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((300.0 * exp((-2.535e-07) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0)))); // per_second
            const NekDouble var_fast_sodium_current_j_gate__j = var_fast_sodium_current__j; // dimensionless
            const NekDouble var_fast_sodium_current_j_gate__d_j_d_environment__time = (var_fast_sodium_current_j_gate__alpha_j * (1.0 - var_fast_sodium_current_j_gate__j)) - (var_fast_sodium_current_j_gate__beta_j * var_fast_sodium_current_j_gate__j); // per_second
            const NekDouble var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time = var_fast_sodium_current_j_gate__d_j_d_environment__time; // per_second
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_factor = 1.0; // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V = var_rapid_activating_delayed_rectifiyer_K_current__V; // millivolt
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__K21 = exp((-7.677) - (0.0128 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V)); // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__K12 = exp((-5.495) + (0.1691 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V)); // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_X_kr = (0.001 / (var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__K12 + var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__K21)) + (var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_factor * 0.027); // second
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = var_rapid_activating_delayed_rectifiyer_K_current__X_kr; // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr_inf = var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__K12 / (var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__K12 + var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__K21); // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time = (var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr_inf - var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr) / var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_X_kr; // per_second
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time = var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time; // per_second
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V = var_slow_activating_delayed_rectifiyer_K_current__V; // millivolt
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__tau_X_ks = 0.001 / (((7.19e-05 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) / (1.0 - exp((-0.148) * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)))) + ((0.000131 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) / (exp(0.0687 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) - 1.0))); // second
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks_infinity = 1.0 / (1.0 + exp((-(var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 24.7)) / 13.6)); // dimensionless
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = var_slow_activating_delayed_rectifiyer_K_current__X_ks; // dimensionless
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time = (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks_infinity - var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks) / var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__tau_X_ks; // per_second
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time = var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time; // per_second
            const NekDouble var_transient_outward_potassium_current_X_to1_gate__V = var_transient_outward_potassium_current__V; // millivolt
            const NekDouble var_transient_outward_potassium_current_X_to1_gate__alpha_X_to1 = 45.16 * exp(0.03577 * var_transient_outward_potassium_current_X_to1_gate__V); // per_second
            const NekDouble var_transient_outward_potassium_current_X_to1_gate__X_to1 = var_transient_outward_potassium_current__X_to1; // dimensionless
            const NekDouble var_transient_outward_potassium_current_X_to1_gate__beta_X_to1 = 98.9 * exp((-0.06237) * var_transient_outward_potassium_current_X_to1_gate__V); // per_second
            const NekDouble var_transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time = (var_transient_outward_potassium_current_X_to1_gate__alpha_X_to1 * (1.0 - var_transient_outward_potassium_current_X_to1_gate__X_to1)) - (var_transient_outward_potassium_current_X_to1_gate__beta_X_to1 * var_transient_outward_potassium_current_X_to1_gate__X_to1); // per_second
            const NekDouble var_transient_outward_potassium_current__transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time = var_transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time; // per_second
            const NekDouble var_transient_outward_potassium_current_Y_to1_gate__V = var_transient_outward_potassium_current__V; // millivolt
            const NekDouble var_transient_outward_potassium_current_Y_to1_gate__alpha_Y_to1 = (5.415 * exp((-(var_transient_outward_potassium_current_Y_to1_gate__V + 33.5)) / 5.0)) / (1.0 + (0.051335 * exp((-(var_transient_outward_potassium_current_Y_to1_gate__V + 33.5)) / 5.0))); // per_second
            const NekDouble var_transient_outward_potassium_current_Y_to1_gate__Y_to1 = var_transient_outward_potassium_current__Y_to1; // dimensionless
            const NekDouble var_transient_outward_potassium_current_Y_to1_gate__beta_Y_to1 = (5.415 * exp((var_transient_outward_potassium_current_Y_to1_gate__V + 33.5) / 5.0)) / (1.0 + (0.051335 * exp((var_transient_outward_potassium_current_Y_to1_gate__V + 33.5) / 5.0))); // per_second
            const NekDouble var_transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time = (var_transient_outward_potassium_current_Y_to1_gate__alpha_Y_to1 * (1.0 - var_transient_outward_potassium_current_Y_to1_gate__Y_to1)) - (var_transient_outward_potassium_current_Y_to1_gate__beta_Y_to1 * var_transient_outward_potassium_current_Y_to1_gate__Y_to1); // per_second
            const NekDouble var_transient_outward_potassium_current__transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time = var_transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time; // per_second
            const NekDouble var_L_type_Ca_current__alpha = 400.0 * exp((var_L_type_Ca_current__V + 2.0) / 10.0); // per_second
            const NekDouble var_L_type_Ca_current__beta = 50.0 * exp((-(var_L_type_Ca_current__V + 2.0)) / 13.0); // per_second
            const NekDouble var_L_type_Ca_current__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_L_type_Ca_current__gamma = (103.75 * var_L_type_Ca_current__Ca_ss) / 1.0; // per_second
            const NekDouble var_L_type_Ca_current__a = 2.0; // dimensionless
            const NekDouble var_L_type_Ca_current__alpha_a = var_L_type_Ca_current__alpha * var_L_type_Ca_current__a; // per_second
            const NekDouble var_L_type_Ca_current__b = 2.0; // dimensionless
            const NekDouble var_L_type_Ca_current__beta_b = var_L_type_Ca_current__beta / var_L_type_Ca_current__b; // per_second
            const NekDouble var_L_type_Ca_current__g = 2000.0; // per_second
            const NekDouble var_L_type_Ca_current__f = 300.0; // per_second
            const NekDouble var_L_type_Ca_current__gprime = 7000.0; // per_second
            const NekDouble var_L_type_Ca_current__fprime = 7.0; // per_second
            const NekDouble var_L_type_Ca_current__omega = 10.0; // per_second
            const NekDouble var_L_type_Ca_current__C0 = var_chaste_interface__L_type_Ca_current__C0; // dimensionless
            const NekDouble var_L_type_Ca_current__C1 = var_chaste_interface__L_type_Ca_current__C1; // dimensionless
            const NekDouble var_L_type_Ca_current__C2 = var_chaste_interface__L_type_Ca_current__C2; // dimensionless
            const NekDouble var_L_type_Ca_current__C3 = var_chaste_interface__L_type_Ca_current__C3; // dimensionless
            const NekDouble var_L_type_Ca_current__C4 = var_chaste_interface__L_type_Ca_current__C4; // dimensionless
            const NekDouble var_L_type_Ca_current__C_Ca0 = var_chaste_interface__L_type_Ca_current__C_Ca0; // dimensionless
            const NekDouble var_L_type_Ca_current__C_Ca1 = var_chaste_interface__L_type_Ca_current__C_Ca1; // dimensionless
            const NekDouble var_L_type_Ca_current__C_Ca2 = var_chaste_interface__L_type_Ca_current__C_Ca2; // dimensionless
            const NekDouble var_L_type_Ca_current__C_Ca3 = var_chaste_interface__L_type_Ca_current__C_Ca3; // dimensionless
            const NekDouble var_L_type_Ca_current__C_Ca4 = var_chaste_interface__L_type_Ca_current__C_Ca4; // dimensionless
            const NekDouble var_L_type_Ca_current__d_O_d_environment__time = (var_L_type_Ca_current__f * var_L_type_Ca_current__C4) - (var_L_type_Ca_current__g * var_L_type_Ca_current__O); // per_second
            const NekDouble var_L_type_Ca_current__d_O_Ca_d_environment__time = (var_L_type_Ca_current__fprime * var_L_type_Ca_current__C_Ca4) - (var_L_type_Ca_current__gprime * var_L_type_Ca_current__O_Ca); // per_second
            const NekDouble var_L_type_Ca_current__d_C0_d_environment__time = ((var_L_type_Ca_current__beta * var_L_type_Ca_current__C1) + (var_L_type_Ca_current__omega * var_L_type_Ca_current__C_Ca0)) - (((4.0 * var_L_type_Ca_current__alpha) + var_L_type_Ca_current__gamma) * var_L_type_Ca_current__C0); // per_second
            const NekDouble var_L_type_Ca_current__d_C1_d_environment__time = ((4.0 * var_L_type_Ca_current__alpha * var_L_type_Ca_current__C0) + (2.0 * var_L_type_Ca_current__beta * var_L_type_Ca_current__C2) + ((var_L_type_Ca_current__omega / var_L_type_Ca_current__b) * var_L_type_Ca_current__C_Ca1)) - ((var_L_type_Ca_current__beta + (3.0 * var_L_type_Ca_current__alpha) + (var_L_type_Ca_current__gamma * var_L_type_Ca_current__a)) * var_L_type_Ca_current__C1); // per_second
            const NekDouble var_L_type_Ca_current__d_C2_d_environment__time = ((3.0 * var_L_type_Ca_current__alpha * var_L_type_Ca_current__C1) + (3.0 * var_L_type_Ca_current__beta * var_L_type_Ca_current__C3) + ((var_L_type_Ca_current__omega / pow(var_L_type_Ca_current__b, 2.0)) * var_L_type_Ca_current__C_Ca2)) - (((var_L_type_Ca_current__beta * 2.0) + (2.0 * var_L_type_Ca_current__alpha) + (var_L_type_Ca_current__gamma * pow(var_L_type_Ca_current__a, 2.0))) * var_L_type_Ca_current__C2); // per_second
            const NekDouble var_L_type_Ca_current__d_C3_d_environment__time = ((2.0 * var_L_type_Ca_current__alpha * var_L_type_Ca_current__C2) + (4.0 * var_L_type_Ca_current__beta * var_L_type_Ca_current__C4) + ((var_L_type_Ca_current__omega / pow(var_L_type_Ca_current__b, 3.0)) * var_L_type_Ca_current__C_Ca3)) - (((var_L_type_Ca_current__beta * 3.0) + var_L_type_Ca_current__alpha + (var_L_type_Ca_current__gamma * pow(var_L_type_Ca_current__a, 3.0))) * var_L_type_Ca_current__C3); // per_second
            const NekDouble var_L_type_Ca_current__d_C4_d_environment__time = ((var_L_type_Ca_current__alpha * var_L_type_Ca_current__C3) + (var_L_type_Ca_current__g * var_L_type_Ca_current__O) + ((var_L_type_Ca_current__omega / pow(var_L_type_Ca_current__b, 4.0)) * var_L_type_Ca_current__C_Ca4)) - (((var_L_type_Ca_current__beta * 4.0) + var_L_type_Ca_current__f + (var_L_type_Ca_current__gamma * pow(var_L_type_Ca_current__a, 4.0))) * var_L_type_Ca_current__C4); // per_second
            const NekDouble var_L_type_Ca_current__d_C_Ca0_d_environment__time = ((var_L_type_Ca_current__beta_b * var_L_type_Ca_current__C_Ca1) + (var_L_type_Ca_current__gamma * var_L_type_Ca_current__C0)) - (((4.0 * var_L_type_Ca_current__alpha_a) + var_L_type_Ca_current__omega) * var_L_type_Ca_current__C_Ca0); // per_second
            const NekDouble var_L_type_Ca_current__d_C_Ca1_d_environment__time = ((4.0 * var_L_type_Ca_current__alpha_a * var_L_type_Ca_current__C_Ca0) + (2.0 * var_L_type_Ca_current__beta_b * var_L_type_Ca_current__C_Ca2) + (var_L_type_Ca_current__gamma * var_L_type_Ca_current__a * var_L_type_Ca_current__C1)) - ((var_L_type_Ca_current__beta_b + (3.0 * var_L_type_Ca_current__alpha_a) + (var_L_type_Ca_current__omega / var_L_type_Ca_current__b)) * var_L_type_Ca_current__C_Ca1); // per_second
            const NekDouble var_L_type_Ca_current__d_C_Ca2_d_environment__time = ((3.0 * var_L_type_Ca_current__alpha_a * var_L_type_Ca_current__C_Ca1) + (3.0 * var_L_type_Ca_current__beta_b * var_L_type_Ca_current__C_Ca3) + (var_L_type_Ca_current__gamma * pow(var_L_type_Ca_current__a, 2.0) * var_L_type_Ca_current__C2)) - (((var_L_type_Ca_current__beta_b * 2.0) + (2.0 * var_L_type_Ca_current__alpha_a) + (var_L_type_Ca_current__omega / pow(var_L_type_Ca_current__b, 2.0))) * var_L_type_Ca_current__C_Ca2); // per_second
            const NekDouble var_L_type_Ca_current__d_C_Ca3_d_environment__time = ((2.0 * var_L_type_Ca_current__alpha_a * var_L_type_Ca_current__C_Ca2) + (4.0 * var_L_type_Ca_current__beta_b * var_L_type_Ca_current__C_Ca4) + (var_L_type_Ca_current__gamma * pow(var_L_type_Ca_current__a, 3.0) * var_L_type_Ca_current__C3)) - (((var_L_type_Ca_current__beta_b * 3.0) + var_L_type_Ca_current__alpha_a + (var_L_type_Ca_current__omega / pow(var_L_type_Ca_current__b, 3.0))) * var_L_type_Ca_current__C_Ca3); // per_second
            const NekDouble var_L_type_Ca_current__d_C_Ca4_d_environment__time = ((var_L_type_Ca_current__alpha_a * var_L_type_Ca_current__C_Ca3) + (var_L_type_Ca_current__gprime * var_L_type_Ca_current__O_Ca) + (var_L_type_Ca_current__gamma * pow(var_L_type_Ca_current__a, 4.0) * var_L_type_Ca_current__C4)) - (((var_L_type_Ca_current__beta_b * 4.0) + var_L_type_Ca_current__fprime + (var_L_type_Ca_current__omega / pow(var_L_type_Ca_current__b, 4.0))) * var_L_type_Ca_current__C_Ca4); // per_second
            const NekDouble var_L_type_Ca_current_y_gate__y = var_L_type_Ca_current__y; // dimensionless
            const NekDouble var_L_type_Ca_current_y_gate__V = var_L_type_Ca_current__V; // millivolt
            const NekDouble var_L_type_Ca_current_y_gate__y_infinity = (0.8 / (1.0 + exp((var_L_type_Ca_current_y_gate__V + 12.5) / 5.0))) + 0.2; // dimensionless
            const NekDouble var_L_type_Ca_current_y_gate__tau_y = (20.0 + (600.0 / (1.0 + exp((var_L_type_Ca_current_y_gate__V + 20.0) / 9.5)))) / 1000.0; // second
            const NekDouble var_L_type_Ca_current_y_gate__d_y_d_environment__time = (var_L_type_Ca_current_y_gate__y_infinity - var_L_type_Ca_current_y_gate__y) / var_L_type_Ca_current_y_gate__tau_y; // per_second
            const NekDouble var_L_type_Ca_current__L_type_Ca_current_y_gate__d_y_d_environment__time = var_L_type_Ca_current_y_gate__d_y_d_environment__time; // per_second
            const NekDouble var_RyR_channel__P_O2 = var_chaste_interface__RyR_channel__P_O2; // dimensionless
            const NekDouble var_RyR_channel__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_RyR_channel__P_O1 = var_chaste_interface__RyR_channel__P_O1; // dimensionless
            const NekDouble var_RyR_channel__v1 = 1800.0; // per_second
            const NekDouble var_RyR_channel__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__Ca_JSR; // millimolar
            const NekDouble var_RyR_channel__J_rel = var_RyR_channel__v1 * (var_RyR_channel__P_O1 + var_RyR_channel__P_O2) * (var_RyR_channel__Ca_JSR - var_RyR_channel__Ca_ss); // millimolar_per_second
            const NekDouble var_RyR_channel__k_a_plus = 1.215e+13; // millimolar4_per_second
            const NekDouble var_RyR_channel__k_a_minus = 576.0; // per_second
            const NekDouble var_RyR_channel__k_b_plus = 4050000000.0; // millimolar3_per_second
            const NekDouble var_RyR_channel__k_b_minus = 1930.0; // per_second
            const NekDouble var_RyR_channel__k_c_plus = 100.0; // per_second
            const NekDouble var_RyR_channel__k_c_minus = 0.8; // per_second
            const NekDouble var_RyR_channel__P_C1 = var_chaste_interface__RyR_channel__P_C1; // dimensionless
            const NekDouble var_RyR_channel__P_C2 = var_chaste_interface__RyR_channel__P_C2; // dimensionless
            const NekDouble var_RyR_channel__n = 4.0; // dimensionless
            const NekDouble var_RyR_channel__m = 3.0; // dimensionless
            const NekDouble var_RyR_channel__d_P_O1_d_environment__time = ((var_RyR_channel__k_a_plus * pow(var_RyR_channel__Ca_ss, var_RyR_channel__n) * var_RyR_channel__P_C1) - ((var_RyR_channel__k_a_minus * var_RyR_channel__P_O1) + (var_RyR_channel__k_b_plus * pow(var_RyR_channel__Ca_ss, var_RyR_channel__m) * var_RyR_channel__P_O1) + (var_RyR_channel__k_c_plus * var_RyR_channel__P_O1))) + (var_RyR_channel__k_b_minus * var_RyR_channel__P_O2) + (var_RyR_channel__k_c_minus * var_RyR_channel__P_C2); // per_second
            const NekDouble var_RyR_channel__d_P_O2_d_environment__time = (var_RyR_channel__k_b_plus * pow(var_RyR_channel__Ca_ss, var_RyR_channel__m) * var_RyR_channel__P_O1) - (var_RyR_channel__k_b_minus * var_RyR_channel__P_O2); // per_second
            const NekDouble var_RyR_channel__d_P_C1_d_environment__time = ((-var_RyR_channel__k_a_plus) * pow(var_RyR_channel__Ca_ss, var_RyR_channel__n) * var_RyR_channel__P_C1) + (var_RyR_channel__k_a_minus * var_RyR_channel__P_O1); // per_second
            const NekDouble var_RyR_channel__d_P_C2_d_environment__time = (var_RyR_channel__k_c_plus * var_RyR_channel__P_O1) - (var_RyR_channel__k_c_minus * var_RyR_channel__P_C2); // per_second
            const NekDouble var_SERCA2a_pump__Cai = var_chaste_interface__intracellular_ion_concentrations__Cai; // millimolar
            const NekDouble var_SERCA2a_pump__N_fb = 1.2; // dimensionless
            const NekDouble var_SERCA2a_pump__K_fb = 0.000168; // millimolar
            const NekDouble var_SERCA2a_pump__fb = pow(var_SERCA2a_pump__Cai / var_SERCA2a_pump__K_fb, var_SERCA2a_pump__N_fb); // dimensionless
            const NekDouble var_SERCA2a_pump__Vmaxf = 0.0813; // millimolar_per_second
            const NekDouble var_SERCA2a_pump__K_SR = 1.0; // dimensionless
            const NekDouble var_SERCA2a_pump__Ca_NSR = var_chaste_interface__intracellular_ion_concentrations__Ca_NSR; // millimolar
            const NekDouble var_SERCA2a_pump__K_rb = 3.29; // millimolar
            const NekDouble var_SERCA2a_pump__N_rb = 1.0; // dimensionless
            const NekDouble var_SERCA2a_pump__rb = pow(var_SERCA2a_pump__Ca_NSR / var_SERCA2a_pump__K_rb, var_SERCA2a_pump__N_rb); // dimensionless
            const NekDouble var_SERCA2a_pump__Vmaxr = 0.318; // millimolar_per_second
            const NekDouble var_SERCA2a_pump__J_up = (var_SERCA2a_pump__K_SR * ((var_SERCA2a_pump__Vmaxf * var_SERCA2a_pump__fb) - (var_SERCA2a_pump__Vmaxr * var_SERCA2a_pump__rb))) / (1.0 + var_SERCA2a_pump__fb + var_SERCA2a_pump__rb); // millimolar_per_second
            const NekDouble var_intracellular_Ca_fluxes__Ca_NSR = var_chaste_interface__intracellular_ion_concentrations__Ca_NSR; // millimolar
            const NekDouble var_intracellular_Ca_fluxes__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__Ca_JSR; // millimolar
            const NekDouble var_intracellular_Ca_fluxes__tau_tr = 0.0005747; // second
            const NekDouble var_intracellular_Ca_fluxes__J_tr = (var_intracellular_Ca_fluxes__Ca_NSR - var_intracellular_Ca_fluxes__Ca_JSR) / var_intracellular_Ca_fluxes__tau_tr; // millimolar_per_second
            const NekDouble var_intracellular_Ca_fluxes__Cai = var_chaste_interface__intracellular_ion_concentrations__Cai; // millimolar
            const NekDouble var_intracellular_Ca_fluxes__tau_xfer = 0.0267; // second
            const NekDouble var_intracellular_Ca_fluxes__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_intracellular_Ca_fluxes__J_xfer = (var_intracellular_Ca_fluxes__Ca_ss - var_intracellular_Ca_fluxes__Cai) / var_intracellular_Ca_fluxes__tau_xfer; // millimolar_per_second
            const NekDouble var_intracellular_Ca_fluxes__k_htrpn_minus = 0.066; // per_second
            const NekDouble var_intracellular_Ca_fluxes__k_htrpn_plus = 20000.0; // per_millimolar_second
            const NekDouble var_intracellular_Ca_fluxes__HTRPNCa = var_chaste_interface__intracellular_Ca_fluxes__HTRPNCa; // millimolar
            const NekDouble var_intracellular_Ca_fluxes__d_HTRPNCa_d_environment__time = (var_intracellular_Ca_fluxes__k_htrpn_plus * var_intracellular_Ca_fluxes__Cai * (1.0 - var_intracellular_Ca_fluxes__HTRPNCa)) - (var_intracellular_Ca_fluxes__k_htrpn_minus * var_intracellular_Ca_fluxes__HTRPNCa); // 'millimole per litre per second'
            const NekDouble var_intracellular_Ca_fluxes__J_HTRPNCa = var_intracellular_Ca_fluxes__d_HTRPNCa_d_environment__time; // millimolar_per_second
            const NekDouble var_intracellular_Ca_fluxes__LTRPN_tot = 0.07; // dimensionless
            const NekDouble var_intracellular_Ca_fluxes__LTRPNCa = var_chaste_interface__intracellular_Ca_fluxes__LTRPNCa; // millimolar
            const NekDouble var_intracellular_Ca_fluxes__k_ltrpn_minus = 40.0; // per_second
            const NekDouble var_intracellular_Ca_fluxes__k_ltrpn_plus = 40000.0; // per_millimolar_second
            const NekDouble var_intracellular_Ca_fluxes__d_LTRPNCa_d_environment__time = (var_intracellular_Ca_fluxes__k_ltrpn_plus * var_intracellular_Ca_fluxes__Cai * (1.0 - var_intracellular_Ca_fluxes__LTRPNCa)) - (var_intracellular_Ca_fluxes__k_ltrpn_minus * var_intracellular_Ca_fluxes__LTRPNCa); // 'millimole per litre per second'
            const NekDouble var_intracellular_Ca_fluxes__J_LTRPNCa = var_intracellular_Ca_fluxes__d_LTRPNCa_d_environment__time; // millimolar_per_second
            const NekDouble var_intracellular_Ca_fluxes__HTRPN_tot = 0.14; // dimensionless
            const NekDouble var_intracellular_Ca_fluxes__J_trpn = (var_intracellular_Ca_fluxes__HTRPN_tot * var_intracellular_Ca_fluxes__J_HTRPNCa) + (var_intracellular_Ca_fluxes__LTRPN_tot * var_intracellular_Ca_fluxes__J_LTRPNCa); // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__Cai = var_chaste_interface__intracellular_ion_concentrations__Cai; // millimolar
            const NekDouble var_intracellular_ion_concentrations__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_intracellular_ion_concentrations__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__Ca_JSR; // millimolar
            const NekDouble var_intracellular_ion_concentrations__A_cap = 0.0001534; // cm2
            const NekDouble var_intracellular_ion_concentrations__V_myo = 2.584e-05; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__V_JSR = 1.6e-07; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__V_NSR = 2.1e-06; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__V_SS = 1.2e-09; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__K_mCMDN = 0.00238; // millimolar
            const NekDouble var_intracellular_ion_concentrations__K_mEGTA = 0.00015; // millimolar
            const NekDouble var_intracellular_ion_concentrations__K_mCSQN = 0.8; // millimolar
            const NekDouble var_intracellular_ion_concentrations__CMDN_tot = 0.05; // millimolar
            const NekDouble var_intracellular_ion_concentrations__EGTA_tot = 0.0; // millimolar
            const NekDouble var_intracellular_ion_concentrations__CSQN_tot = 15.0; // millimolar
            const NekDouble var_intracellular_ion_concentrations__beta_i = 1.0 / (1.0 + ((var_intracellular_ion_concentrations__CMDN_tot * var_intracellular_ion_concentrations__K_mCMDN) / pow(var_intracellular_ion_concentrations__K_mCMDN + var_intracellular_ion_concentrations__Cai, 2.0)) + ((var_intracellular_ion_concentrations__EGTA_tot * var_intracellular_ion_concentrations__K_mEGTA) / pow(var_intracellular_ion_concentrations__K_mEGTA + var_intracellular_ion_concentrations__Cai, 2.0))); // dimensionless
            const NekDouble var_intracellular_ion_concentrations__beta_SS = 1.0 / (1.0 + ((var_intracellular_ion_concentrations__CMDN_tot * var_intracellular_ion_concentrations__K_mCMDN) / pow(var_intracellular_ion_concentrations__K_mCMDN + var_intracellular_ion_concentrations__Ca_ss, 2.0)) + ((var_intracellular_ion_concentrations__EGTA_tot * var_intracellular_ion_concentrations__K_mEGTA) / pow(var_intracellular_ion_concentrations__K_mEGTA + var_intracellular_ion_concentrations__Ca_ss, 2.0))); // dimensionless
            const NekDouble var_intracellular_ion_concentrations__beta_JSR = 1.0 / (1.0 + ((var_intracellular_ion_concentrations__CSQN_tot * var_intracellular_ion_concentrations__K_mCSQN) / pow(var_intracellular_ion_concentrations__K_mCSQN + var_intracellular_ion_concentrations__Ca_JSR, 2.0))); // dimensionless
            const NekDouble var_intracellular_ion_concentrations__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_intracellular_ion_concentrations__i_Na = var_fast_sodium_current__i_Na; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_Ca = var_L_type_Ca_current__i_Ca; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_Na_b = var_sodium_background_current__i_Na_b; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_NaCa = var_Na_Ca_exchanger__i_NaCa; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_NaK = var_sodium_potassium_pump__i_NaK; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_Ca_K = var_L_type_Ca_current__i_Ca_K; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__i_Kr; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__i_Ks; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_K1 = var_time_independent_potassium_current__i_K1; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_Kp = var_plateau_potassium_current__i_Kp; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_to1 = var_transient_outward_potassium_current__i_to1; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__i_Ca_b = var_calcium_background_current__i_Ca_b; // microA_per_microF
            const NekDouble var_intracellular_ion_concentrations__J_up = var_SERCA2a_pump__J_up; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_rel = var_RyR_channel__J_rel; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_xfer = var_intracellular_Ca_fluxes__J_xfer; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_trpn = var_intracellular_Ca_fluxes__J_trpn; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_tr = var_intracellular_Ca_fluxes__J_tr; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__d_Nai_d_environment__time = ((-0.0) * (var_intracellular_ion_concentrations__i_Na + var_intracellular_ion_concentrations__i_Na_b + (var_intracellular_ion_concentrations__i_NaCa * 3.0) + (var_intracellular_ion_concentrations__i_NaK * 3.0)) * var_intracellular_ion_concentrations__A_cap * 1.0) / (var_intracellular_ion_concentrations__V_myo * var_intracellular_ion_concentrations__F); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Cai_d_environment__time = var_intracellular_ion_concentrations__beta_i * ((var_intracellular_ion_concentrations__J_xfer - (var_intracellular_ion_concentrations__J_up + var_intracellular_ion_concentrations__J_trpn)) + ((((2.0 * var_intracellular_ion_concentrations__i_NaCa) - (var_intracellular_ion_concentrations__i_p_Ca + var_intracellular_ion_concentrations__i_Ca_b)) * var_intracellular_ion_concentrations__A_cap * 1.0) / (2.0 * var_intracellular_ion_concentrations__V_myo * var_intracellular_ion_concentrations__F))); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ki_d_environment__time = ((-0.0) * (var_intracellular_ion_concentrations__i_Ca_K + var_intracellular_ion_concentrations__i_Kr + var_intracellular_ion_concentrations__i_Ks + var_intracellular_ion_concentrations__i_K1 + var_intracellular_ion_concentrations__i_Kp + var_intracellular_ion_concentrations__i_to1 + (var_intracellular_ion_concentrations__i_NaK * (-2.0))) * var_intracellular_ion_concentrations__A_cap * 1.0) / (var_intracellular_ion_concentrations__V_myo * var_intracellular_ion_concentrations__F); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_ss_d_environment__time = var_intracellular_ion_concentrations__beta_SS * ((((var_intracellular_ion_concentrations__J_rel * var_intracellular_ion_concentrations__V_JSR) / var_intracellular_ion_concentrations__V_SS) - ((var_intracellular_ion_concentrations__J_xfer * var_intracellular_ion_concentrations__V_myo) / var_intracellular_ion_concentrations__V_SS)) - ((var_intracellular_ion_concentrations__i_Ca * var_intracellular_ion_concentrations__A_cap * 1.0) / (2.0 * var_intracellular_ion_concentrations__V_SS * var_intracellular_ion_concentrations__F))); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_JSR_d_environment__time = var_intracellular_ion_concentrations__beta_JSR * (var_intracellular_ion_concentrations__J_tr - var_intracellular_ion_concentrations__J_rel); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_NSR_d_environment__time = ((var_intracellular_ion_concentrations__J_up * var_intracellular_ion_concentrations__V_myo) / var_intracellular_ion_concentrations__V_NSR) - ((var_intracellular_ion_concentrations__J_tr * var_intracellular_ion_concentrations__V_JSR) / var_intracellular_ion_concentrations__V_NSR); // 'millimole per litre per second'
            const NekDouble var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time_converter = var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time; // per_second
            const NekDouble var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time = 0.001 * var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time_converter = var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time; // per_second
            const NekDouble var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time = 0.001 * var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time_converter = var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time; // per_second
            const NekDouble var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time = 0.001 * var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time_converter = var_rapid_activating_delayed_rectifiyer_K_current__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time; // per_second
            const NekDouble var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time = 0.001 * var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time_converter = var_slow_activating_delayed_rectifiyer_K_current__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time; // per_second
            const NekDouble var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time = 0.001 * var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time_converter = var_transient_outward_potassium_current__transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time = 0.001 * var_chaste_interface__transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time_converter = var_transient_outward_potassium_current__transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time = 0.001 * var_chaste_interface__transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_O_d_environment__time_converter = var_L_type_Ca_current__d_O_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_O_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_O_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_O_Ca_d_environment__time_converter = var_L_type_Ca_current__d_O_Ca_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_O_Ca_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_O_Ca_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C0_d_environment__time_converter = var_L_type_Ca_current__d_C0_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C0_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C0_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C1_d_environment__time_converter = var_L_type_Ca_current__d_C1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C1_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C2_d_environment__time_converter = var_L_type_Ca_current__d_C2_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C2_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C2_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C3_d_environment__time_converter = var_L_type_Ca_current__d_C3_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C3_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C3_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C4_d_environment__time_converter = var_L_type_Ca_current__d_C4_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C4_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C4_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca0_d_environment__time_converter = var_L_type_Ca_current__d_C_Ca0_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca0_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C_Ca0_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca1_d_environment__time_converter = var_L_type_Ca_current__d_C_Ca1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca1_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C_Ca1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca2_d_environment__time_converter = var_L_type_Ca_current__d_C_Ca2_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca2_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C_Ca2_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca3_d_environment__time_converter = var_L_type_Ca_current__d_C_Ca3_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca3_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C_Ca3_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca4_d_environment__time_converter = var_L_type_Ca_current__d_C_Ca4_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current__d_C_Ca4_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current__d_C_Ca4_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_current_y_gate__d_y_d_environment__time_converter = var_L_type_Ca_current__L_type_Ca_current_y_gate__d_y_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_current_y_gate__d_y_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_current_y_gate__d_y_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__RyR_channel__d_P_O1_d_environment__time_converter = var_RyR_channel__d_P_O1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__RyR_channel__d_P_O1_d_environment__time = 0.001 * var_chaste_interface__RyR_channel__d_P_O1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__RyR_channel__d_P_O2_d_environment__time_converter = var_RyR_channel__d_P_O2_d_environment__time; // per_second
            const NekDouble var_chaste_interface__RyR_channel__d_P_O2_d_environment__time = 0.001 * var_chaste_interface__RyR_channel__d_P_O2_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__RyR_channel__d_P_C1_d_environment__time_converter = var_RyR_channel__d_P_C1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__RyR_channel__d_P_C1_d_environment__time = 0.001 * var_chaste_interface__RyR_channel__d_P_C1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__RyR_channel__d_P_C2_d_environment__time_converter = var_RyR_channel__d_P_C2_d_environment__time; // per_second
            const NekDouble var_chaste_interface__RyR_channel__d_P_C2_d_environment__time = 0.001 * var_chaste_interface__RyR_channel__d_P_C2_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__intracellular_Ca_fluxes__d_HTRPNCa_d_environment__time_converter = var_intracellular_Ca_fluxes__d_HTRPNCa_d_environment__time; // ___units_89
            const NekDouble var_chaste_interface__intracellular_Ca_fluxes__d_HTRPNCa_d_environment__time = 0.001 * var_chaste_interface__intracellular_Ca_fluxes__d_HTRPNCa_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_Ca_fluxes__d_LTRPNCa_d_environment__time_converter = var_intracellular_Ca_fluxes__d_LTRPNCa_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_Ca_fluxes__d_LTRPNCa_d_environment__time = 0.001 * var_chaste_interface__intracellular_Ca_fluxes__d_LTRPNCa_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Nai_d_environment__time_converter = var_intracellular_ion_concentrations__d_Nai_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Nai_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Nai_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Cai_d_environment__time_converter = var_intracellular_ion_concentrations__d_Cai_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Cai_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Cai_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ki_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ki_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ki_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ki_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ca_ss_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ca_JSR_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ca_NSR_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble d_dt_chaste_interface__fast_sodium_current_m_gate__m = var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__fast_sodium_current_h_gate__h = var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__fast_sodium_current_j_gate__j = var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__transient_outward_potassium_current_X_to1_gate__X_to1 = var_chaste_interface__transient_outward_potassium_current_X_to1_gate__d_X_to1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__transient_outward_potassium_current_Y_to1_gate__Y_to1 = var_chaste_interface__transient_outward_potassium_current_Y_to1_gate__d_Y_to1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__O = var_chaste_interface__L_type_Ca_current__d_O_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__O_Ca = var_chaste_interface__L_type_Ca_current__d_O_Ca_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C0 = var_chaste_interface__L_type_Ca_current__d_C0_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C1 = var_chaste_interface__L_type_Ca_current__d_C1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C2 = var_chaste_interface__L_type_Ca_current__d_C2_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C3 = var_chaste_interface__L_type_Ca_current__d_C3_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C4 = var_chaste_interface__L_type_Ca_current__d_C4_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C_Ca0 = var_chaste_interface__L_type_Ca_current__d_C_Ca0_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C_Ca1 = var_chaste_interface__L_type_Ca_current__d_C_Ca1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C_Ca2 = var_chaste_interface__L_type_Ca_current__d_C_Ca2_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C_Ca3 = var_chaste_interface__L_type_Ca_current__d_C_Ca3_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current__C_Ca4 = var_chaste_interface__L_type_Ca_current__d_C_Ca4_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_current_y_gate__y = var_chaste_interface__L_type_Ca_current_y_gate__d_y_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__RyR_channel__P_O1 = var_chaste_interface__RyR_channel__d_P_O1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__RyR_channel__P_O2 = var_chaste_interface__RyR_channel__d_P_O2_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__RyR_channel__P_C1 = var_chaste_interface__RyR_channel__d_P_C1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__RyR_channel__P_C2 = var_chaste_interface__RyR_channel__d_P_C2_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_Ca_fluxes__HTRPNCa = var_chaste_interface__intracellular_Ca_fluxes__d_HTRPNCa_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_Ca_fluxes__LTRPNCa = var_chaste_interface__intracellular_Ca_fluxes__d_LTRPNCa_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Nai = var_chaste_interface__intracellular_ion_concentrations__d_Nai_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Cai = var_chaste_interface__intracellular_ion_concentrations__d_Cai_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ki = var_chaste_interface__intracellular_ion_concentrations__d_Ki_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_NSR = var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time; // 'millimole per litre per millisecond'
            
            const NekDouble var_membrane__C_sc = 0.001; // microF_per_cm2
            const NekDouble var_membrane__i_Na = var_fast_sodium_current__i_Na; // microA_per_microF
            const NekDouble var_membrane__i_Ca = var_L_type_Ca_current__i_Ca; // microA_per_microF
            const NekDouble var_membrane__i_Ca_K = var_L_type_Ca_current__i_Ca_K; // microA_per_microF
            const NekDouble var_membrane__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__i_Kr; // microA_per_microF
            const NekDouble var_membrane__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__i_Ks; // microA_per_microF
            const NekDouble var_membrane__i_to1 = var_transient_outward_potassium_current__i_to1; // microA_per_microF
            const NekDouble var_membrane__i_K1 = var_time_independent_potassium_current__i_K1; // microA_per_microF
            const NekDouble var_membrane__i_Kp = var_plateau_potassium_current__i_Kp; // microA_per_microF
            const NekDouble var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa; // microA_per_microF
            const NekDouble var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK; // microA_per_microF
            const NekDouble var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca; // microA_per_microF
            const NekDouble var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b; // microA_per_microF
            const NekDouble var_membrane__i_Na_b = var_sodium_background_current__i_Na_b; // microA_per_microF
            const NekDouble var_chaste_interface__membrane__i_Stim = 0.0;
            const NekDouble var_membrane__i_Stim_converter = var_chaste_interface__membrane__i_Stim; // uA_per_cm2
            const NekDouble var_membrane__chaste_interface__chaste_membrane_capacitance = 1.0; // uF_per_cm2
            const NekDouble var_membrane__i_Stim = var_membrane__i_Stim_converter / var_membrane__chaste_interface__chaste_membrane_capacitance; // microA_per_microF
            const NekDouble var_membrane__d_V_d_environment__time = ((-1.0) * 1.0 * (var_membrane__i_Na + var_membrane__i_Ca + var_membrane__i_Ca_K + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_to1 + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_NaCa + var_membrane__i_NaK + var_membrane__i_p_Ca + var_membrane__i_Na_b + var_membrane__i_Ca_b + var_membrane__i_Stim)) / var_membrane__C_sc; // 'millivolt per second'
            const NekDouble var_chaste_interface__membrane__d_V_d_environment__time_converter = var_membrane__d_V_d_environment__time; // ___units_1
            const NekDouble var_chaste_interface__membrane__d_V_d_environment__time = 0.001 * var_chaste_interface__membrane__d_V_d_environment__time_converter; // 'millivolt per millisecond'
            d_dt_chaste_interface__membrane__V = var_chaste_interface__membrane__d_V_d_environment__time; // 'millivolt per millisecond'
            outarray[0][i] = d_dt_chaste_interface__membrane__V;
            outarray[1][i] = d_dt_chaste_interface__fast_sodium_current_m_gate__m;
            outarray[2][i] = d_dt_chaste_interface__fast_sodium_current_h_gate__h;
            outarray[3][i] = d_dt_chaste_interface__fast_sodium_current_j_gate__j;
            outarray[4][i] = d_dt_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr;
            outarray[5][i] = d_dt_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks;
            outarray[6][i] = d_dt_chaste_interface__transient_outward_potassium_current_X_to1_gate__X_to1;
            outarray[7][i] = d_dt_chaste_interface__transient_outward_potassium_current_Y_to1_gate__Y_to1;
            outarray[8][i] = d_dt_chaste_interface__L_type_Ca_current__O;
            outarray[9][i] = d_dt_chaste_interface__L_type_Ca_current__O_Ca;
            outarray[10][i] = d_dt_chaste_interface__L_type_Ca_current__C0;
            outarray[11][i] = d_dt_chaste_interface__L_type_Ca_current__C1;
            outarray[12][i] = d_dt_chaste_interface__L_type_Ca_current__C2;
            outarray[13][i] = d_dt_chaste_interface__L_type_Ca_current__C3;
            outarray[14][i] = d_dt_chaste_interface__L_type_Ca_current__C4;
            outarray[15][i] = d_dt_chaste_interface__L_type_Ca_current__C_Ca0;
            outarray[16][i] = d_dt_chaste_interface__L_type_Ca_current__C_Ca1;
            outarray[17][i] = d_dt_chaste_interface__L_type_Ca_current__C_Ca2;
            outarray[18][i] = d_dt_chaste_interface__L_type_Ca_current__C_Ca3;
            outarray[19][i] = d_dt_chaste_interface__L_type_Ca_current__C_Ca4;
            outarray[20][i] = d_dt_chaste_interface__L_type_Ca_current_y_gate__y;
            outarray[21][i] = d_dt_chaste_interface__RyR_channel__P_O1;
            outarray[22][i] = d_dt_chaste_interface__RyR_channel__P_O2;
            outarray[23][i] = d_dt_chaste_interface__RyR_channel__P_C1;
            outarray[24][i] = d_dt_chaste_interface__RyR_channel__P_C2;
            outarray[25][i] = d_dt_chaste_interface__intracellular_Ca_fluxes__HTRPNCa;
            outarray[26][i] = d_dt_chaste_interface__intracellular_Ca_fluxes__LTRPNCa;
            outarray[27][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Nai;
            outarray[28][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Cai;
            outarray[29][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ki;
            outarray[30][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_ss;
            outarray[31][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_JSR;
            outarray[32][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_NSR;
        }
        
    }

    /**
    *
    */
    void Winslow99::v_GenerateSummary(SummaryList& s)
    {
        SolverUtils::AddSummaryItem(s, "Cell model", "Winslow99");
    }


    /**
     *
     */
    void Winslow99::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, -96.1638,     m_cellSol[0],  1);
        Vmath::Fill(m_nq, 0.0328302,    m_cellSol[1],  1);
        Vmath::Fill(m_nq, 0.988354,     m_cellSol[2],  1);
        Vmath::Fill(m_nq, 0.99254,      m_cellSol[3],  1);
        Vmath::Fill(m_nq, 0.51,         m_cellSol[4],  1);
        Vmath::Fill(m_nq, 0.264,        m_cellSol[5],  1);
        Vmath::Fill(m_nq, 2.63,         m_cellSol[6],  1);
        Vmath::Fill(m_nq, 0.99,         m_cellSol[7],  1);
        Vmath::Fill(m_nq, 9.84546e-21,  m_cellSol[8],  1);
        Vmath::Fill(m_nq, 0.0,          m_cellSol[9],  1);
        Vmath::Fill(m_nq, 0.997208,     m_cellSol[10], 1);
        Vmath::Fill(m_nq, 6.38897e-5,   m_cellSol[11], 1);
        Vmath::Fill(m_nq, 1.535e-9,     m_cellSol[12], 1);
        Vmath::Fill(m_nq, 1.63909e-14,  m_cellSol[13], 1);
        Vmath::Fill(m_nq, 6.56337e-20,  m_cellSol[14], 1);
        Vmath::Fill(m_nq, 0.00272826,   m_cellSol[15], 1);
        Vmath::Fill(m_nq, 6.99215e-7,   m_cellSol[16], 1);
        Vmath::Fill(m_nq, 6.71989e-11,  m_cellSol[17], 1);
        Vmath::Fill(m_nq, 2.87031e-15,  m_cellSol[18], 1);
        Vmath::Fill(m_nq, 4.59752e-20,  m_cellSol[19], 1);
        Vmath::Fill(m_nq, 0.798,        m_cellSol[20], 1);
        Vmath::Fill(m_nq, 0.0,          m_cellSol[21],  1);
        Vmath::Fill(m_nq, 0.0,          m_cellSol[22], 1);
        Vmath::Fill(m_nq, 0.47,         m_cellSol[23], 1);
        Vmath::Fill(m_nq, 0.53,         m_cellSol[24], 1);
        Vmath::Fill(m_nq, 0.98,         m_cellSol[25], 1);
        Vmath::Fill(m_nq, 0.078,        m_cellSol[26], 1);
        Vmath::Fill(m_nq, 10.0,         m_cellSol[27], 1);
        Vmath::Fill(m_nq, 0.00008,      m_cellSol[28], 1);
        Vmath::Fill(m_nq, 157.8,        m_cellSol[29], 1);
        Vmath::Fill(m_nq, 0.00011,      m_cellSol[30], 1);
        Vmath::Fill(m_nq, 0.257,        m_cellSol[31], 1);
        Vmath::Fill(m_nq, 0.257,        m_cellSol[32], 1);
    }
}

