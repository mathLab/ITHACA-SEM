///////////////////////////////////////////////////////////////////////////////
//
// File Fox02.cpp
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
// Description: Fox 2002 cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
//#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <CardiacEPSolver/CellModels/Fox02.h>

namespace Nektar
{
    std::string Fox02::className
              = GetCellModelFactory().RegisterCreatorFunction(
                        "Fox02",
                        Fox02::create,
                         "Fox 2002 cell model.");
    
    
    /**
    *
    */
    Fox02::Fox02(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const MultiRegions::ExpListSharedPtr& pField):
        CellModel(pSession, pField)
    {
        pSession->LoadParameter("chi",       m_chi);
        pSession->LoadParameter("sigmai",    m_sigmai);

        m_nvar = 13;

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
        m_concentrations.push_back(11);
        m_concentrations.push_back(12);
    }
    
    
    
    
    
    void Fox02::v_Update(
                     const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                           Array<OneD,        Array<OneD, NekDouble> >&outarray,
                     const NekDouble time)
    {
        for (unsigned int i = 0; i < m_nq; ++i)
        {
            
            // Inputs:
            // Time units: millisecond
            NekDouble var_chaste_interface__membrane__V = inarray[0][i];
            // Units: millivolt; Initial value: -94.7
            NekDouble var_chaste_interface__fast_sodium_current_m_gate__m = inarray[1][i];
            // Units: dimensionless; Initial value: 0.00024676
            NekDouble var_chaste_interface__fast_sodium_current_h_gate__h = inarray[2][i];
            // Units: dimensionless; Initial value: 0.99869
            NekDouble var_chaste_interface__fast_sodium_current_j_gate__j = inarray[3][i];
            // Units: dimensionless; Initial value: 0.99887
            NekDouble var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = inarray[4][i];
            // Units: dimensionless; Initial value: 0.229
            NekDouble var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = inarray[5][i];
            // Units: dimensionless; Initial value: 0.0001
            NekDouble var_chaste_interface__transient_outward_potassium_current_X_to_gate__X_to = inarray[6][i];
            // Units: dimensionless; Initial value: 0.00003742
            NekDouble var_chaste_interface__transient_outward_potassium_current_Y_to_gate__Y_to = inarray[7][i];
            // Units: dimensionless; Initial value: 1
            NekDouble var_chaste_interface__L_type_Ca_current_f_gate__f = inarray[8][i];
            // Units: dimensionless; Initial value: 0.983
            NekDouble var_chaste_interface__L_type_Ca_current_d_gate__d = inarray[9][i];
            // Units: dimensionless; Initial value: 0.0001
            NekDouble var_chaste_interface__L_type_Ca_current_f_Ca_gate__f_Ca = inarray[10][i];
            // Units: dimensionless; Initial value: 0.942
            NekDouble var_chaste_interface__calcium_dynamics__Ca_i = inarray[11][i];
            // Units: micromolar; Initial value: 0.0472
            NekDouble var_chaste_interface__calcium_dynamics__Ca_SR = inarray[12][i];
            // Units: micromolar; Initial value: 320


            // Mathematics
            NekDouble d_dt_chaste_interface__membrane__V;
            const NekDouble var_membrane__R = 8.314; // joule_per_mole_kelvin
            const NekDouble var_membrane__T = 310.0; // kelvin
            const NekDouble var_membrane__F = 96.5; // coulomb_per_millimole
            const NekDouble var_fast_sodium_current__j = var_chaste_interface__fast_sodium_current_j_gate__j; // dimensionless
            const NekDouble var_fast_sodium_current__h = var_chaste_interface__fast_sodium_current_h_gate__h; // dimensionless
            const NekDouble var_fast_sodium_current__m = var_chaste_interface__fast_sodium_current_m_gate__m; // dimensionless
            const NekDouble var_fast_sodium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_standard_ionic_concentrations__Na_o = 138.0; // millimolar
            const NekDouble var_standard_ionic_concentrations__Na_i = 10.0; // millimolar
            const NekDouble var_L_type_Ca_current__f = var_chaste_interface__L_type_Ca_current_f_gate__f; // dimensionless
            const NekDouble var_L_type_Ca_current__f_Ca = var_chaste_interface__L_type_Ca_current_f_Ca_gate__f_Ca; // dimensionless
            const NekDouble var_L_type_Ca_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_L_type_Ca_current__P_Ca = 1.26e-05; // cm_per_millisecond
            const NekDouble var_L_type_Ca_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_L_type_Ca_current__C_sc = 1.0; // microF_per_cm2
            const NekDouble var_L_type_Ca_current__T = var_membrane__T; // kelvin
            const NekDouble var_standard_ionic_concentrations__Ca_o = 2000.0; // micromolar
            const NekDouble var_L_type_Ca_current__Ca_o = var_standard_ionic_concentrations__Ca_o; // micromolar
            const NekDouble var_L_type_Ca_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_L_type_Ca_current__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // micromolar
            const NekDouble var_L_type_Ca_current__i_Ca_max = ((((var_L_type_Ca_current__P_Ca / var_L_type_Ca_current__C_sc) * 4.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__Ca_i * exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - (0.341 * var_L_type_Ca_current__Ca_o))) / (exp((2.0 * var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0); // microA_per_microF
            const NekDouble var_L_type_Ca_current__d = var_chaste_interface__L_type_Ca_current_d_gate__d; // dimensionless
            const NekDouble var_L_type_Ca_current__i_Ca = var_L_type_Ca_current__i_Ca_max * var_L_type_Ca_current__f * var_L_type_Ca_current__d * var_L_type_Ca_current__f_Ca; // microA_per_microF
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__X_kr = var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr; // dimensionless
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__X_ks = var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks; // dimensionless
            const NekDouble var_transient_outward_potassium_current__Y_to = var_chaste_interface__transient_outward_potassium_current_Y_to_gate__Y_to; // dimensionless
            const NekDouble var_transient_outward_potassium_current__X_to = var_chaste_interface__transient_outward_potassium_current_X_to_gate__X_to; // dimensionless
            const NekDouble var_transient_outward_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_Na_Ca_exchanger__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_Na_Ca_exchanger__K_sat = 0.2; // dimensionless
            const NekDouble var_Na_Ca_exchanger__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // micromolar
            const NekDouble var_Na_Ca_exchanger__K_mNa = 87.5; // millimolar
            const NekDouble var_Na_Ca_exchanger__Ca_o = var_standard_ionic_concentrations__Ca_o; // micromolar
            const NekDouble var_Na_Ca_exchanger__K_NaCa = 1500.0; // microA_per_microF
            const NekDouble var_Na_Ca_exchanger__T = var_membrane__T; // kelvin
            const NekDouble var_Na_Ca_exchanger__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_Na_Ca_exchanger__Na_i = var_standard_ionic_concentrations__Na_i; // millimolar
            const NekDouble var_Na_Ca_exchanger__eta = 0.35; // dimensionless
            const NekDouble var_Na_Ca_exchanger__K_mCa = 1380.0; // micromolar
            const NekDouble var_Na_Ca_exchanger__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_Na_Ca_exchanger__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_Na_Ca_exchanger__i_NaCa = (var_Na_Ca_exchanger__K_NaCa / ((pow(var_Na_Ca_exchanger__K_mNa, 3.0) + pow(var_Na_Ca_exchanger__Na_o, 3.0)) * (var_Na_Ca_exchanger__K_mCa + var_Na_Ca_exchanger__Ca_o) * (1.0 + (var_Na_Ca_exchanger__K_sat * exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)))))) * ((exp((var_Na_Ca_exchanger__eta * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_i, 3.0) * var_Na_Ca_exchanger__Ca_o) - (exp(((var_Na_Ca_exchanger__eta - 1.0) * var_Na_Ca_exchanger__V * var_Na_Ca_exchanger__F) / (var_Na_Ca_exchanger__R * var_Na_Ca_exchanger__T)) * pow(var_Na_Ca_exchanger__Na_o, 3.0) * var_Na_Ca_exchanger__Ca_i)); // microA_per_microF
            const NekDouble var_sarcolemmal_calcium_pump__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // micromolar
            const NekDouble var_sarcolemmal_calcium_pump__i_pCa_max = 0.05; // microA_per_microF
            const NekDouble var_sarcolemmal_calcium_pump__K_mpCa = 0.05; // micromolar
            const NekDouble var_sarcolemmal_calcium_pump__i_p_Ca = (var_sarcolemmal_calcium_pump__i_pCa_max * var_sarcolemmal_calcium_pump__Ca_i) / (var_sarcolemmal_calcium_pump__K_mpCa + var_sarcolemmal_calcium_pump__Ca_i); // microA_per_microF
            const NekDouble var_calcium_background_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_calcium_background_current__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // micromolar
            const NekDouble var_calcium_background_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_calcium_background_current__T = var_membrane__T; // kelvin
            const NekDouble var_calcium_background_current__Ca_o = var_standard_ionic_concentrations__Ca_o; // micromolar
            const NekDouble var_calcium_background_current__E_Ca = ((var_calcium_background_current__R * var_calcium_background_current__T) / (2.0 * var_calcium_background_current__F)) * log(var_calcium_background_current__Ca_o / var_calcium_background_current__Ca_i); // millivolt
            const NekDouble var_calcium_background_current__g_Cab = 0.0003842; // milliS_per_microF
            const NekDouble var_calcium_background_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_calcium_background_current__i_Ca_b = var_calcium_background_current__g_Cab * (var_calcium_background_current__V - var_calcium_background_current__E_Ca); // microA_per_microF
            const NekDouble var_fast_sodium_current_m_gate__m = var_fast_sodium_current__m; // dimensionless
            const NekDouble var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_m_gate__E0_m = var_fast_sodium_current_m_gate__V + 47.13; // millivolt
            const NekDouble var_fast_sodium_current_m_gate__alpha_m = (0.32 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp((-0.1) * var_fast_sodium_current_m_gate__E0_m)); // per_millisecond
            const NekDouble var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_fast_sodium_current_m_gate__V) / 11.0); // per_millisecond
            const NekDouble var_fast_sodium_current_m_gate__d_m_d_environment__time = (var_fast_sodium_current_m_gate__alpha_m * (1.0 - var_fast_sodium_current_m_gate__m)) - (var_fast_sodium_current_m_gate__beta_m * var_fast_sodium_current_m_gate__m); // per_millisecond
            const NekDouble var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time = var_fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
            const NekDouble var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_h_gate__beta_h = 7.5 / (1.0 + exp((-0.1) * (var_fast_sodium_current_h_gate__V + 11.0))); // per_millisecond
            const NekDouble var_fast_sodium_current_h_gate__alpha_h = 0.135 * exp((var_fast_sodium_current_h_gate__V + 80.0) / (-6.8)); // per_millisecond
            const NekDouble var_fast_sodium_current_h_gate__h = var_fast_sodium_current__h; // dimensionless
            const NekDouble var_fast_sodium_current_h_gate__d_h_d_environment__time = (var_fast_sodium_current_h_gate__alpha_h * (1.0 - var_fast_sodium_current_h_gate__h)) - (var_fast_sodium_current_h_gate__beta_h * var_fast_sodium_current_h_gate__h); // per_millisecond
            const NekDouble var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time = var_fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
            const NekDouble var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_j_gate__alpha_j = (0.175 * exp((var_fast_sodium_current_j_gate__V + 100.0) / (-23.0))) / (1.0 + exp(0.15 * (var_fast_sodium_current_j_gate__V + 79.0))); // per_millisecond
            const NekDouble var_fast_sodium_current_j_gate__beta_j = 0.3 / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0))); // per_millisecond
            const NekDouble var_fast_sodium_current_j_gate__j = var_fast_sodium_current__j; // dimensionless
            const NekDouble var_fast_sodium_current_j_gate__d_j_d_environment__time = (var_fast_sodium_current_j_gate__alpha_j * (1.0 - var_fast_sodium_current_j_gate__j)) - (var_fast_sodium_current_j_gate__beta_j * var_fast_sodium_current_j_gate__j); // per_millisecond
            const NekDouble var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time = var_fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V = var_rapid_activating_delayed_rectifiyer_K_current__V; // millivolt
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_X_kr = 43.0 + (1.0 / (exp((-5.495) + (0.1691 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V)) + exp((-7.677) - (0.0128 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V)))); // millisecond
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = var_rapid_activating_delayed_rectifiyer_K_current__X_kr; // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr_inf = 1.0 / (1.0 + exp((-2.182) - (0.1819 * var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__V))); // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time = (var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr_inf - var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr) / var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__tau_X_kr; // per_millisecond
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time = var_rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time; // per_millisecond
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V = var_slow_activating_delayed_rectifiyer_K_current__V; // millivolt
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__tau_X_ks = 1.0 / (((7.19e-05 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) / (1.0 - exp((-0.148) * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)))) + ((0.000131 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) / (exp(0.0687 * (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 10.0)) - 1.0))); // millisecond
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks_infinity = 1.0 / (1.0 + exp((var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__V - 16.0) / (-13.6))); // dimensionless
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = var_slow_activating_delayed_rectifiyer_K_current__X_ks; // dimensionless
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time = (var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks_infinity - var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks) / var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__tau_X_ks; // per_millisecond
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time = var_slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time; // per_millisecond
            const NekDouble var_transient_outward_potassium_current_X_to_gate__V = var_transient_outward_potassium_current__V; // millivolt
            const NekDouble var_transient_outward_potassium_current_X_to_gate__alpha_X_to = 0.04516 * exp(0.03577 * var_transient_outward_potassium_current_X_to_gate__V); // per_millisecond
            const NekDouble var_transient_outward_potassium_current_X_to_gate__X_to = var_transient_outward_potassium_current__X_to; // dimensionless
            const NekDouble var_transient_outward_potassium_current_X_to_gate__beta_X_to = 0.0989 * exp((-0.06237) * var_transient_outward_potassium_current_X_to_gate__V); // per_millisecond
            const NekDouble var_transient_outward_potassium_current_X_to_gate__d_X_to_d_environment__time = (var_transient_outward_potassium_current_X_to_gate__alpha_X_to * (1.0 - var_transient_outward_potassium_current_X_to_gate__X_to)) - (var_transient_outward_potassium_current_X_to_gate__beta_X_to * var_transient_outward_potassium_current_X_to_gate__X_to); // per_millisecond
            const NekDouble var_transient_outward_potassium_current__transient_outward_potassium_current_X_to_gate__d_X_to_d_environment__time = var_transient_outward_potassium_current_X_to_gate__d_X_to_d_environment__time; // per_millisecond
            const NekDouble var_transient_outward_potassium_current_Y_to_gate__V = var_transient_outward_potassium_current__V; // millivolt
            const NekDouble var_transient_outward_potassium_current_Y_to_gate__beta_Y_to = (0.005415 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / 5.0)) / (1.0 + (0.051335 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / 5.0))); // per_millisecond
            const NekDouble var_transient_outward_potassium_current_Y_to_gate__alpha_Y_to = (0.005415 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / (-5.0))) / (1.0 + (0.051335 * exp((var_transient_outward_potassium_current_Y_to_gate__V + 33.5) / (-5.0)))); // per_millisecond
            const NekDouble var_transient_outward_potassium_current_Y_to_gate__Y_to = var_transient_outward_potassium_current__Y_to; // dimensionless
            const NekDouble var_transient_outward_potassium_current_Y_to_gate__d_Y_to_d_environment__time = (var_transient_outward_potassium_current_Y_to_gate__alpha_Y_to * (1.0 - var_transient_outward_potassium_current_Y_to_gate__Y_to)) - (var_transient_outward_potassium_current_Y_to_gate__beta_Y_to * var_transient_outward_potassium_current_Y_to_gate__Y_to); // per_millisecond
            const NekDouble var_transient_outward_potassium_current__transient_outward_potassium_current_Y_to_gate__d_Y_to_d_environment__time = var_transient_outward_potassium_current_Y_to_gate__d_Y_to_d_environment__time; // per_millisecond
            const NekDouble var_L_type_Ca_current_f_gate__V = var_L_type_Ca_current__V; // millivolt
            const NekDouble var_L_type_Ca_current_f_gate__tau_f = 30.0 + (200.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 20.0) / 9.5))); // millisecond
            const NekDouble var_L_type_Ca_current_f_gate__f = var_L_type_Ca_current__f; // dimensionless
            const NekDouble var_L_type_Ca_current_f_gate__f_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 12.5) / 5.0)); // dimensionless
            const NekDouble var_L_type_Ca_current_f_gate__d_f_d_environment__time = (var_L_type_Ca_current_f_gate__f_infinity - var_L_type_Ca_current_f_gate__f) / var_L_type_Ca_current_f_gate__tau_f; // per_millisecond
            const NekDouble var_L_type_Ca_current__L_type_Ca_current_f_gate__d_f_d_environment__time = var_L_type_Ca_current_f_gate__d_f_d_environment__time; // per_millisecond
            const NekDouble var_L_type_Ca_current_d_gate__V = var_L_type_Ca_current__V; // millivolt
            const NekDouble var_L_type_Ca_current_d_gate__E0_m = var_L_type_Ca_current_d_gate__V + 40.0; // millivolt
            const NekDouble var_L_type_Ca_current_d_gate__tau_d = 1.0 / (((0.25 * exp((-0.01) * var_L_type_Ca_current_d_gate__V)) / (1.0 + exp((-0.07) * var_L_type_Ca_current_d_gate__V))) + ((0.07 * exp((-0.05) * var_L_type_Ca_current_d_gate__E0_m)) / (1.0 + exp(0.05 * var_L_type_Ca_current_d_gate__E0_m)))); // millisecond
            const NekDouble var_L_type_Ca_current_d_gate__d = var_L_type_Ca_current__d; // dimensionless
            const NekDouble var_L_type_Ca_current_d_gate__d_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_current_d_gate__V + 10.0) / (-6.24))); // dimensionless
            const NekDouble var_L_type_Ca_current_d_gate__d_d_d_environment__time = (var_L_type_Ca_current_d_gate__d_infinity - var_L_type_Ca_current_d_gate__d) / var_L_type_Ca_current_d_gate__tau_d; // per_millisecond
            const NekDouble var_L_type_Ca_current__L_type_Ca_current_d_gate__d_d_d_environment__time = var_L_type_Ca_current_d_gate__d_d_d_environment__time; // per_millisecond
            const NekDouble var_L_type_Ca_current_f_Ca_gate__f_Ca = var_L_type_Ca_current__f_Ca; // dimensionless
            const NekDouble var_L_type_Ca_current_f_Ca_gate__K_mfCa = 0.18; // micromolar
            const NekDouble var_L_type_Ca_current_f_Ca_gate__Ca_i = var_L_type_Ca_current__Ca_i; // micromolar
            const NekDouble var_L_type_Ca_current_f_Ca_gate__f_Ca_infinity = 1.0 / (1.0 + pow(var_L_type_Ca_current_f_Ca_gate__Ca_i / var_L_type_Ca_current_f_Ca_gate__K_mfCa, 3.0)); // dimensionless
            const NekDouble var_L_type_Ca_current_f_Ca_gate__tau_f_Ca = 30.0; // millisecond
            const NekDouble var_L_type_Ca_current_f_Ca_gate__d_f_Ca_d_environment__time = (var_L_type_Ca_current_f_Ca_gate__f_Ca_infinity - var_L_type_Ca_current_f_Ca_gate__f_Ca) / var_L_type_Ca_current_f_Ca_gate__tau_f_Ca; // per_millisecond
            const NekDouble var_L_type_Ca_current__L_type_Ca_current_f_Ca_gate__d_f_Ca_d_environment__time = var_L_type_Ca_current_f_Ca_gate__d_f_Ca_d_environment__time; // per_millisecond
            const NekDouble var_calcium_dynamics__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // micromolar
            const NekDouble var_calcium_dynamics__CMDN_tot = 10.0; // micromolar
            const NekDouble var_calcium_dynamics__K_mCMDN = 2.0; // micromolar
            const NekDouble var_calcium_dynamics__beta_i = 1.0 / (1.0 + ((var_calcium_dynamics__CMDN_tot * var_calcium_dynamics__K_mCMDN) / pow(var_calcium_dynamics__K_mCMDN + var_calcium_dynamics__Ca_i, 2.0))); // dimensionless
            const NekDouble var_calcium_dynamics__V_myo = 2.584e-05; // microlitre
            const NekDouble var_calcium_dynamics__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_calcium_dynamics__C_sc = var_L_type_Ca_current__C_sc; // microF_per_cm2
            const NekDouble var_calcium_dynamics__A_Cap = 0.0001534; // cm2
            const NekDouble var_calcium_dynamics__Ca_SR = var_chaste_interface__calcium_dynamics__Ca_SR; // micromolar
            const NekDouble var_calcium_dynamics__P_rel = 6.0; // per_millisecond
            const NekDouble var_calcium_dynamics__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_calcium_dynamics__f_Ca = var_chaste_interface__L_type_Ca_current_f_Ca_gate__f_Ca; // dimensionless
            const NekDouble var_calcium_dynamics__gamma = 1.0 / (1.0 + pow(2000.0 / var_calcium_dynamics__Ca_SR, 3.0)); // dimensionless
            const NekDouble var_calcium_dynamics__d = var_chaste_interface__L_type_Ca_current_d_gate__d; // dimensionless
            const NekDouble var_calcium_dynamics__f = var_chaste_interface__L_type_Ca_current_f_gate__f; // dimensionless
            const NekDouble var_calcium_dynamics__J_rel = (var_calcium_dynamics__P_rel * var_calcium_dynamics__f * var_calcium_dynamics__d * var_calcium_dynamics__f_Ca * ((var_calcium_dynamics__gamma * var_calcium_dynamics__Ca_SR) - var_calcium_dynamics__Ca_i)) / (1.0 + (1.65 * exp(var_calcium_dynamics__V / 20.0))); // micromolar_per_millisecond
            const NekDouble var_calcium_dynamics__P_leak = 1e-06; // per_millisecond
            const NekDouble var_calcium_dynamics__J_leak = var_calcium_dynamics__P_leak * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_i); // micromolar_per_millisecond
            const NekDouble var_calcium_dynamics__K_mup = 0.32; // micromolar
            const NekDouble var_calcium_dynamics__V_up = 0.1; // micromolar_per_millisecond
            const NekDouble var_calcium_dynamics__J_up = var_calcium_dynamics__V_up / (1.0 + pow(var_calcium_dynamics__K_mup / var_calcium_dynamics__Ca_i, 2.0)); // micromolar_per_millisecond
            const NekDouble var_calcium_dynamics__i_Ca = var_L_type_Ca_current__i_Ca; // microA_per_microF
            const NekDouble var_calcium_dynamics__i_Ca_b = var_calcium_background_current__i_Ca_b; // microA_per_microF
            const NekDouble var_calcium_dynamics__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca; // microA_per_microF
            const NekDouble var_calcium_dynamics__i_NaCa = var_Na_Ca_exchanger__i_NaCa; // microA_per_microF
            const NekDouble var_calcium_dynamics__K_mCSQN = 600.0; // micromolar
            const NekDouble var_calcium_dynamics__CSQN_tot = 10000.0; // micromolar
            const NekDouble var_calcium_dynamics__V_SR = 2e-06; // microlitre
            const NekDouble var_calcium_dynamics__beta_SR = 1.0 / (1.0 + ((var_calcium_dynamics__CSQN_tot * var_calcium_dynamics__K_mCSQN) / pow(var_calcium_dynamics__K_mCSQN + var_calcium_dynamics__Ca_SR, 2.0))); // dimensionless
            const NekDouble var_calcium_dynamics__d_Ca_i_d_environment__time = var_calcium_dynamics__beta_i * (((var_calcium_dynamics__J_rel + var_calcium_dynamics__J_leak) - var_calcium_dynamics__J_up) - (((var_calcium_dynamics__A_Cap * var_calcium_dynamics__C_sc) / (2.0 * var_calcium_dynamics__F * var_calcium_dynamics__V_myo)) * ((var_calcium_dynamics__i_Ca + var_calcium_dynamics__i_Ca_b + var_calcium_dynamics__i_p_Ca) - (2.0 * var_calcium_dynamics__i_NaCa)))); // 'micromole per litre per millisecond'
            const NekDouble var_calcium_dynamics__d_Ca_SR_d_environment__time = (var_calcium_dynamics__beta_SR * ((var_calcium_dynamics__J_up - var_calcium_dynamics__J_leak) - var_calcium_dynamics__J_rel) * var_calcium_dynamics__V_myo) / var_calcium_dynamics__V_SR; // 'micromole per litre per millisecond'
            const NekDouble var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time = var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time = var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time = var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time = var_rapid_activating_delayed_rectifiyer_K_current__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time = var_slow_activating_delayed_rectifiyer_K_current__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__transient_outward_potassium_current_X_to_gate__d_X_to_d_environment__time = var_transient_outward_potassium_current__transient_outward_potassium_current_X_to_gate__d_X_to_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__transient_outward_potassium_current_Y_to_gate__d_Y_to_d_environment__time = var_transient_outward_potassium_current__transient_outward_potassium_current_Y_to_gate__d_Y_to_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__L_type_Ca_current_f_gate__d_f_d_environment__time = var_L_type_Ca_current__L_type_Ca_current_f_gate__d_f_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__L_type_Ca_current_d_gate__d_d_d_environment__time = var_L_type_Ca_current__L_type_Ca_current_d_gate__d_d_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__L_type_Ca_current_f_Ca_gate__d_f_Ca_d_environment__time = var_L_type_Ca_current__L_type_Ca_current_f_Ca_gate__d_f_Ca_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__calcium_dynamics__d_Ca_i_d_environment__time = var_calcium_dynamics__d_Ca_i_d_environment__time; // micromolar_per_millisecond
            const NekDouble var_chaste_interface__calcium_dynamics__d_Ca_SR_d_environment__time = var_calcium_dynamics__d_Ca_SR_d_environment__time; // micromolar_per_millisecond
            const NekDouble d_dt_chaste_interface__fast_sodium_current_m_gate__m = var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__fast_sodium_current_h_gate__h = var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__fast_sodium_current_j_gate__j = var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr = var_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__d_X_kr_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks = var_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__d_X_ks_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__transient_outward_potassium_current_X_to_gate__X_to = var_chaste_interface__transient_outward_potassium_current_X_to_gate__d_X_to_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__transient_outward_potassium_current_Y_to_gate__Y_to = var_chaste_interface__transient_outward_potassium_current_Y_to_gate__d_Y_to_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__L_type_Ca_current_f_gate__f = var_chaste_interface__L_type_Ca_current_f_gate__d_f_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__L_type_Ca_current_d_gate__d = var_chaste_interface__L_type_Ca_current_d_gate__d_d_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__L_type_Ca_current_f_Ca_gate__f_Ca = var_chaste_interface__L_type_Ca_current_f_Ca_gate__d_f_Ca_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__calcium_dynamics__Ca_i = var_chaste_interface__calcium_dynamics__d_Ca_i_d_environment__time; // 'micromole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__calcium_dynamics__Ca_SR = var_chaste_interface__calcium_dynamics__d_Ca_SR_d_environment__time; // 'micromole per litre per millisecond'
            
            const NekDouble var_fast_sodium_current__g_Na = 12.8; // milliS_per_microF
            const NekDouble var_fast_sodium_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_fast_sodium_current__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_fast_sodium_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_fast_sodium_current__T = var_membrane__T; // kelvin
            const NekDouble var_fast_sodium_current__Na_i = var_standard_ionic_concentrations__Na_i; // millimolar
            const NekDouble var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Na_o / var_fast_sodium_current__Na_i); // millivolt
            const NekDouble var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na); // microA_per_microF
            const NekDouble var_membrane__i_Na = var_fast_sodium_current__i_Na; // microA_per_microF
            const NekDouble var_membrane__i_Ca = var_L_type_Ca_current__i_Ca; // microA_per_microF
            const NekDouble var_L_type_Ca_current__P_CaK = 5.79e-07; // cm_per_millisecond
            const NekDouble var_standard_ionic_concentrations__K_o = 4.0; // millimolar
            const NekDouble var_L_type_Ca_current__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_L_type_Ca_current__i_Ca_half =  -0.265; // microA_per_microF
            const NekDouble var_standard_ionic_concentrations__K_i = 149.4; // millimolar
            const NekDouble var_L_type_Ca_current__K_i = var_standard_ionic_concentrations__K_i; // millimolar
            const NekDouble var_L_type_Ca_current__i_CaK = ((((((var_L_type_Ca_current__P_CaK / var_L_type_Ca_current__C_sc) * var_L_type_Ca_current__f * var_L_type_Ca_current__d * var_L_type_Ca_current__f_Ca) / (1.0 + (var_L_type_Ca_current__i_Ca_max / var_L_type_Ca_current__i_Ca_half))) * 1000.0 * var_L_type_Ca_current__V * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((var_L_type_Ca_current__K_i * exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__K_o)) / (exp((var_L_type_Ca_current__V * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0); // microA_per_microF
            const NekDouble var_membrane__i_CaK = var_L_type_Ca_current__i_CaK; // microA_per_microF
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__K_i = var_standard_ionic_concentrations__K_i; // millimolar
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__T = var_membrane__T; // kelvin
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__E_K = ((var_rapid_activating_delayed_rectifiyer_K_current__R * var_rapid_activating_delayed_rectifiyer_K_current__T) / var_rapid_activating_delayed_rectifiyer_K_current__F) * log(var_rapid_activating_delayed_rectifiyer_K_current__K_o / var_rapid_activating_delayed_rectifiyer_K_current__K_i); // millivolt
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__R_V = 1.0 / (1.0 + (2.5 * exp(0.1 * (var_rapid_activating_delayed_rectifiyer_K_current__V + 28.0)))); // dimensionless
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__g_Kr = 0.0136; // milliS_per_microF
            const NekDouble var_rapid_activating_delayed_rectifiyer_K_current__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__g_Kr * var_rapid_activating_delayed_rectifiyer_K_current__R_V * var_rapid_activating_delayed_rectifiyer_K_current__X_kr * sqrt(var_rapid_activating_delayed_rectifiyer_K_current__K_o / 4.0) * (var_rapid_activating_delayed_rectifiyer_K_current__V - var_rapid_activating_delayed_rectifiyer_K_current__E_K); // microA_per_microF
            const NekDouble var_membrane__i_Kr = var_rapid_activating_delayed_rectifiyer_K_current__i_Kr; // microA_per_microF
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__g_Ks = 0.0245; // milliS_per_microF
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__T = var_membrane__T; // kelvin
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__Na_i = var_standard_ionic_concentrations__Na_i; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__K_i = var_standard_ionic_concentrations__K_i; // millimolar
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__E_Ks = ((var_slow_activating_delayed_rectifiyer_K_current__R * var_slow_activating_delayed_rectifiyer_K_current__T) / var_slow_activating_delayed_rectifiyer_K_current__F) * log((var_slow_activating_delayed_rectifiyer_K_current__K_o + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Na_o)) / (var_slow_activating_delayed_rectifiyer_K_current__K_i + (0.01833 * var_slow_activating_delayed_rectifiyer_K_current__Na_i))); // millivolt
            const NekDouble var_slow_activating_delayed_rectifiyer_K_current__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__g_Ks * pow(var_slow_activating_delayed_rectifiyer_K_current__X_ks, 2.0) * (var_slow_activating_delayed_rectifiyer_K_current__V - var_slow_activating_delayed_rectifiyer_K_current__E_Ks); // microA_per_microF
            const NekDouble var_membrane__i_Ks = var_slow_activating_delayed_rectifiyer_K_current__i_Ks; // microA_per_microF
            const NekDouble var_transient_outward_potassium_current__g_to = 0.23815; // milliS_per_microF
            const NekDouble var_transient_outward_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K; // millivolt
            const NekDouble var_transient_outward_potassium_current__i_to = var_transient_outward_potassium_current__g_to * var_transient_outward_potassium_current__X_to * var_transient_outward_potassium_current__Y_to * (var_transient_outward_potassium_current__V - var_transient_outward_potassium_current__E_K); // microA_per_microF
            const NekDouble var_membrane__i_to = var_transient_outward_potassium_current__i_to; // microA_per_microF
            const NekDouble var_time_independent_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K; // millivolt
            const NekDouble var_time_independent_potassium_current__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_time_independent_potassium_current__g_K1 = 2.8; // milliS_per_microF
            const NekDouble var_time_independent_potassium_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_time_independent_potassium_current_K1_gate__F = var_time_independent_potassium_current__F; // coulomb_per_millimole
            const NekDouble var_time_independent_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_time_independent_potassium_current_K1_gate__V = var_time_independent_potassium_current__V; // millivolt
            const NekDouble var_time_independent_potassium_current__T = var_membrane__T; // kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__T = var_time_independent_potassium_current__T; // kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__E_K = var_time_independent_potassium_current__E_K; // millivolt
            const NekDouble var_time_independent_potassium_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__R = var_time_independent_potassium_current__R; // joule_per_mole_kelvin
            const NekDouble var_time_independent_potassium_current_K1_gate__K1_infinity = 1.0 / (2.0 + exp(((1.62 * var_time_independent_potassium_current_K1_gate__F) / (var_time_independent_potassium_current_K1_gate__R * var_time_independent_potassium_current_K1_gate__T)) * (var_time_independent_potassium_current_K1_gate__V - var_time_independent_potassium_current_K1_gate__E_K))); // dimensionless
            const NekDouble var_time_independent_potassium_current__K1_infinity = var_time_independent_potassium_current_K1_gate__K1_infinity; // dimensionless
            const NekDouble var_time_independent_potassium_current__K_mK1 = 13.0; // millimolar
            const NekDouble var_time_independent_potassium_current__i_K1 = ((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K1_infinity * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mK1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K); // microA_per_microF
            const NekDouble var_membrane__i_K1 = var_time_independent_potassium_current__i_K1; // microA_per_microF
            const NekDouble var_plateau_potassium_current__g_Kp = 0.002216; // milliS_per_microF
            const NekDouble var_plateau_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_plateau_potassium_current_Kp_gate__V = var_plateau_potassium_current__V; // millivolt
            const NekDouble var_plateau_potassium_current_Kp_gate__Kp_V = 1.0 / (1.0 + exp((7.488 - var_plateau_potassium_current_Kp_gate__V) / 5.98)); // dimensionless
            const NekDouble var_plateau_potassium_current__Kp_V = var_plateau_potassium_current_Kp_gate__Kp_V; // dimensionless
            const NekDouble var_plateau_potassium_current__E_K = var_rapid_activating_delayed_rectifiyer_K_current__E_K; // millivolt
            const NekDouble var_plateau_potassium_current__i_Kp = var_plateau_potassium_current__g_Kp * var_plateau_potassium_current__Kp_V * (var_plateau_potassium_current__V - var_plateau_potassium_current__E_K); // microA_per_microF
            const NekDouble var_membrane__i_Kp = var_plateau_potassium_current__i_Kp; // microA_per_microF
            const NekDouble var_membrane__i_NaCa = var_Na_Ca_exchanger__i_NaCa; // microA_per_microF
            const NekDouble var_sodium_potassium_pump__K_mNai = 10.0; // millimolar
            const NekDouble var_sodium_potassium_pump__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_sodium_potassium_pump__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_sodium_potassium_pump__T = var_membrane__T; // kelvin
            const NekDouble var_sodium_potassium_pump__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Na_o / 67.3) - 1.0); // dimensionless
            const NekDouble var_sodium_potassium_pump__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_sodium_potassium_pump__f_NaK = 1.0 / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__V * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__V) * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T)))); // dimensionless
            const NekDouble var_sodium_potassium_pump__i_NaK_max = 0.693; // microA_per_microF
            const NekDouble var_sodium_potassium_pump__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_sodium_potassium_pump__Na_i = var_standard_ionic_concentrations__Na_i; // millimolar
            const NekDouble var_sodium_potassium_pump__K_mKo = 1.5; // millimolar
            const NekDouble var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__f_NaK) / (1.0 + pow(var_sodium_potassium_pump__K_mNai / var_sodium_potassium_pump__Na_i, 1.5))) * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_o + var_sodium_potassium_pump__K_mKo); // microA_per_microF
            const NekDouble var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK; // microA_per_microF
            const NekDouble var_membrane__i_p_Ca = var_sarcolemmal_calcium_pump__i_p_Ca; // microA_per_microF
            const NekDouble var_membrane__i_Ca_b = var_calcium_background_current__i_Ca_b; // microA_per_microF
            const NekDouble var_sodium_background_current__g_Nab = 0.0031; // milliS_per_microF
            const NekDouble var_sodium_background_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_sodium_background_current__E_Na = var_fast_sodium_current__E_Na; // millivolt
            const NekDouble var_sodium_background_current__i_Na_b = var_sodium_background_current__g_Nab * (var_sodium_background_current__V - var_sodium_background_current__E_Na); // microA_per_microF
            const NekDouble var_membrane__i_Na_b = var_sodium_background_current__i_Na_b; // microA_per_microF
            const NekDouble var_chaste_interface__membrane__i_Stim = 0.0;
            const NekDouble var_membrane__i_Stim_converter = var_chaste_interface__membrane__i_Stim; // uA_per_cm2
            const NekDouble var_membrane__chaste_interface__chaste_membrane_capacitance = 1.0; // uF_per_cm2
            const NekDouble var_membrane__i_Stim = var_membrane__i_Stim_converter / var_membrane__chaste_interface__chaste_membrane_capacitance; // microA_per_microF
            const NekDouble var_membrane__d_V_d_environment__time = -(var_membrane__i_Na + var_membrane__i_Ca + var_membrane__i_CaK + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_to + var_membrane__i_K1 + var_membrane__i_Kp + var_membrane__i_NaCa + var_membrane__i_NaK + var_membrane__i_p_Ca + var_membrane__i_Na_b + var_membrane__i_Ca_b + var_membrane__i_Stim); // 'millivolt per millisecond'
            const NekDouble var_chaste_interface__membrane__d_V_d_environment__time = var_membrane__d_V_d_environment__time; // ___units_1
            d_dt_chaste_interface__membrane__V = var_chaste_interface__membrane__d_V_d_environment__time; // 'millivolt per millisecond'
            outarray[0][i] = d_dt_chaste_interface__membrane__V;
            outarray[1][i] = d_dt_chaste_interface__fast_sodium_current_m_gate__m;
            outarray[2][i] = d_dt_chaste_interface__fast_sodium_current_h_gate__h;
            outarray[3][i] = d_dt_chaste_interface__fast_sodium_current_j_gate__j;
            outarray[4][i] = d_dt_chaste_interface__rapid_activating_delayed_rectifiyer_K_current_X_kr_gate__X_kr;
            outarray[5][i] = d_dt_chaste_interface__slow_activating_delayed_rectifiyer_K_current_X_ks_gate__X_ks;
            outarray[6][i] = d_dt_chaste_interface__transient_outward_potassium_current_X_to_gate__X_to;
            outarray[7][i] = d_dt_chaste_interface__transient_outward_potassium_current_Y_to_gate__Y_to;
            outarray[8][i] = d_dt_chaste_interface__L_type_Ca_current_f_gate__f;
            outarray[9][i] = d_dt_chaste_interface__L_type_Ca_current_d_gate__d;
            outarray[10][i] = d_dt_chaste_interface__L_type_Ca_current_f_Ca_gate__f_Ca;
            outarray[11][i] = d_dt_chaste_interface__calcium_dynamics__Ca_i;
            outarray[12][i] = d_dt_chaste_interface__calcium_dynamics__Ca_SR;
        }
        
    }

    /**
    *
    */
    void Fox02::v_GenerateSummary(SummaryList& s)
    {
        SolverUtils::AddSummaryItem(s, "Cell model", "Fox02");
    }


    /**
     *
     */
    void Fox02::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, -94.7,       m_cellSol[0],  1);
        Vmath::Fill(m_nq, 0.00024676,  m_cellSol[1],  1);
        Vmath::Fill(m_nq, 0.99869,     m_cellSol[2],  1);
        Vmath::Fill(m_nq, 0.99887,     m_cellSol[3],  1);
        Vmath::Fill(m_nq, 0.229,       m_cellSol[4],  1);
        Vmath::Fill(m_nq, 0.0001,      m_cellSol[5],  1);
        Vmath::Fill(m_nq, 0.00003742,  m_cellSol[6],  1);
        Vmath::Fill(m_nq, 1.0,         m_cellSol[7],  1);
        Vmath::Fill(m_nq, 0.0472,      m_cellSol[8],  1);
        Vmath::Fill(m_nq, 0.983,       m_cellSol[9],  1);
        Vmath::Fill(m_nq, 0.0001,      m_cellSol[10], 1);
        Vmath::Fill(m_nq, 0.942,       m_cellSol[11], 1);
        Vmath::Fill(m_nq, 320.0,       m_cellSol[12], 1);
    }
}

