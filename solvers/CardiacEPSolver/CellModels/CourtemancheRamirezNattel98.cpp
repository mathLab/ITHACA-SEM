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
// Description: Courtemanche-Ramirez-Nattel ionic atrial cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <CardiacEPSolver/CellModels/CourtemancheRamirezNattel98.h>
namespace Nektar
{
    std::string CourtemancheRamirezNattel98::className
              = GetCellModelFactory().RegisterCreatorFunction(
                        "CourtemancheRamirezNattel98",
                        CourtemancheRamirezNattel98::create,
                         "Phenomological model of nerve cell electrophysiology.");
    
    
    /**
    *
    */
    CourtemancheRamirezNattel98::CourtemancheRamirezNattel98(
                const LibUtilities::SessionReaderSharedPtr& pSession, const int nq)
            : CellModel(pSession, nq)
    {
        ASSERTL0(pSession->GetVariables().size() == 22,
                 "Aliev-Panfilov cell model requires 22 variables.");

        m_nq = nq;
    }
    
    
    
    /**
    *
    */
    CourtemancheRamirezNattel98::~CourtemancheRamirezNattel98()
    {
        
    }
    
    
    
    /**
    * @param   inarray         Input array.
    * @param   outarray        Output array after addition of reaction terms.
    * @param   time            Current simulation time.
    */
    void CourtemancheRamirezNattel98::v_Update(
                     const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                           Array<OneD,        Array<OneD, NekDouble> >&outarray,
                     const NekDouble time)
    {
        int nvariables  = inarray.num_elements();
        int nq          = m_nq;

        for (unsigned int i = 0; i < nq; ++i)
        {
            // Inputs:
            // Time units: millisecond
            NekDouble var_chaste_interface__membrane__V = inarray[0][i];
            // Units: millivolt; Initial value: -81.18
            NekDouble var_chaste_interface__fast_sodium_current_m_gate__m = inarray[2][i];
            // Units: dimensionless; Initial value: 2.908e-3
            NekDouble var_chaste_interface__fast_sodium_current_h_gate__h = inarray[3][i];
            // Units: dimensionless; Initial value: 9.649e-1
            NekDouble var_chaste_interface__fast_sodium_current_j_gate__j = inarray[4][i];
            // Units: dimensionless; Initial value: 9.775e-1
            NekDouble var_chaste_interface__transient_outward_K_current_oa_gate__oa = inarray[5][i];
            // Units: dimensionless; Initial value: 3.043e-2
            NekDouble var_chaste_interface__transient_outward_K_current_oi_gate__oi = inarray[6][i];
            // Units: dimensionless; Initial value: 9.992e-1
            NekDouble var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ua_gate__ua = inarray[7][i];
            // Units: dimensionless; Initial value: 4.966e-3
            NekDouble var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ui_gate__ui = inarray[8][i];
            // Units: dimensionless; Initial value: 9.986e-1
            NekDouble var_chaste_interface__rapid_delayed_rectifier_K_current_xr_gate__xr = inarray[9][i];
            // Units: dimensionless; Initial value: 3.296e-5
            NekDouble var_chaste_interface__slow_delayed_rectifier_K_current_xs_gate__xs = inarray[10][i];
            // Units: dimensionless; Initial value: 1.869e-2
            NekDouble var_chaste_interface__L_type_Ca_channel_d_gate__d = inarray[11][i];
            // Units: dimensionless; Initial value: 1.367e-4
            NekDouble var_chaste_interface__L_type_Ca_channel_f_gate__f = inarray[12][i];
            // Units: dimensionless; Initial value: 9.996e-1
            NekDouble var_chaste_interface__L_type_Ca_channel_f_Ca_gate__f_Ca = inarray[13][i];
            // Units: dimensionless; Initial value: 7.755e-1
            NekDouble var_chaste_interface__Ca_release_current_from_JSR_u_gate__u = inarray[14][i];
            // Units: dimensionless; Initial value: 2.35e-112
            NekDouble var_chaste_interface__Ca_release_current_from_JSR_v_gate__v = inarray[15][i];
            // Units: dimensionless; Initial value: 1
            NekDouble var_chaste_interface__Ca_release_current_from_JSR_w_gate__w = inarray[16][i];
            // Units: dimensionless; Initial value: 0.9992
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Na_i = inarray[17][i];
            // Units: millimolar; Initial value: 1.117e+01
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_i = inarray[18][i];
            // Units: millimolar; Initial value: 1.013e-4
            NekDouble var_chaste_interface__intracellular_ion_concentrations__K_i = inarray[19][i];
            // Units: millimolar; Initial value: 1.39e+02
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_rel = inarray[20][i];
            // Units: millimolar; Initial value: 1.488
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_up = inarray[21][i];
            // Units: millimolar; Initial value: 1.488

            
            // Mathematics
            NekDouble d_dt_chaste_interface__membrane__V;
            const NekDouble var_membrane__R = 8.3143; // joule_per_mole_kelvin
            const NekDouble var_membrane__T = 310.0; // kelvin
            const NekDouble var_membrane__F = 96.4867; // coulomb_per_millimole
            const NekDouble var_membrane__Cm = 100.0; // picoF
            const NekDouble var_fast_sodium_current__j = var_chaste_interface__fast_sodium_current_j_gate__j; // dimensionless
            const NekDouble var_fast_sodium_current__h = var_chaste_interface__fast_sodium_current_h_gate__h; // dimensionless
            const NekDouble var_fast_sodium_current__g_Na = 7.8; // nanoS_per_picoF
            const NekDouble var_fast_sodium_current__m = var_chaste_interface__fast_sodium_current_m_gate__m; // dimensionless
            const NekDouble var_fast_sodium_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_fast_sodium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_fast_sodium_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_standard_ionic_concentrations__Na_o = 140.0; // millimolar
            const NekDouble var_fast_sodium_current__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_fast_sodium_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_fast_sodium_current__T = var_membrane__T; // kelvin
            const NekDouble var_fast_sodium_current__Na_i = var_chaste_interface__intracellular_ion_concentrations__Na_i; // millimolar
            const NekDouble var_fast_sodium_current__E_Na = ((var_fast_sodium_current__R * var_fast_sodium_current__T) / var_fast_sodium_current__F) * log(var_fast_sodium_current__Na_o / var_fast_sodium_current__Na_i); // millivolt
            const NekDouble var_fast_sodium_current__i_Na = var_fast_sodium_current__Cm * var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na); // picoA
            const NekDouble var_time_independent_potassium_current__g_K1 = 0.09; // nanoS_per_picoF
            const NekDouble var_time_independent_potassium_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_time_independent_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_standard_ionic_concentrations__K_o = 5.4; // millimolar
            const NekDouble var_time_independent_potassium_current__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_time_independent_potassium_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_time_independent_potassium_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_time_independent_potassium_current__K_i = var_chaste_interface__intracellular_ion_concentrations__K_i; // millimolar
            const NekDouble var_time_independent_potassium_current__T = var_membrane__T; // kelvin
            const NekDouble var_time_independent_potassium_current__E_K = ((var_time_independent_potassium_current__R * var_time_independent_potassium_current__T) / var_time_independent_potassium_current__F) * log(var_time_independent_potassium_current__K_o / var_time_independent_potassium_current__K_i); // millivolt
            const NekDouble var_time_independent_potassium_current__i_K1 = (var_time_independent_potassium_current__Cm * var_time_independent_potassium_current__g_K1 * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp(0.07 * (var_time_independent_potassium_current__V + 80.0))); // picoA
            const NekDouble var_transient_outward_K_current__oi = var_chaste_interface__transient_outward_K_current_oi_gate__oi; // dimensionless
            const NekDouble var_transient_outward_K_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_transient_outward_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_transient_outward_K_current__oa = var_chaste_interface__transient_outward_K_current_oa_gate__oa; // dimensionless
            const NekDouble var_transient_outward_K_current__g_to = 0.1652; // nanoS_per_picoF
            const NekDouble var_transient_outward_K_current__E_K = var_time_independent_potassium_current__E_K; // millivolt
            const NekDouble var_transient_outward_K_current__i_to = var_transient_outward_K_current__Cm * var_transient_outward_K_current__g_to * pow(var_transient_outward_K_current__oa, 3.0) * var_transient_outward_K_current__oi * (var_transient_outward_K_current__V - var_transient_outward_K_current__E_K); // picoA
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__g_Kur = 0.005 + (0.05 / (1.0 + exp((var_ultrarapid_delayed_rectifier_K_current__V - 15.0) / (-13.0)))); // nanoS_per_picoF
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__ua = var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ua_gate__ua; // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__E_K = var_time_independent_potassium_current__E_K; // millivolt
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__ui = var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ui_gate__ui; // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__i_Kur = var_ultrarapid_delayed_rectifier_K_current__Cm * var_ultrarapid_delayed_rectifier_K_current__g_Kur * pow(var_ultrarapid_delayed_rectifier_K_current__ua, 3.0) * var_ultrarapid_delayed_rectifier_K_current__ui * (var_ultrarapid_delayed_rectifier_K_current__V - var_ultrarapid_delayed_rectifier_K_current__E_K); // picoA
            const NekDouble var_rapid_delayed_rectifier_K_current__g_Kr = 0.029411765; // nanoS_per_picoF
            const NekDouble var_rapid_delayed_rectifier_K_current__E_K = var_time_independent_potassium_current__E_K; // millivolt
            const NekDouble var_rapid_delayed_rectifier_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_rapid_delayed_rectifier_K_current__xr = var_chaste_interface__rapid_delayed_rectifier_K_current_xr_gate__xr; // dimensionless
            const NekDouble var_rapid_delayed_rectifier_K_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_rapid_delayed_rectifier_K_current__i_Kr = (var_rapid_delayed_rectifier_K_current__Cm * var_rapid_delayed_rectifier_K_current__g_Kr * var_rapid_delayed_rectifier_K_current__xr * (var_rapid_delayed_rectifier_K_current__V - var_rapid_delayed_rectifier_K_current__E_K)) / (1.0 + exp((var_rapid_delayed_rectifier_K_current__V + 15.0) / 22.4)); // picoA
            const NekDouble var_slow_delayed_rectifier_K_current__E_K = var_time_independent_potassium_current__E_K; // millivolt
            const NekDouble var_slow_delayed_rectifier_K_current__xs = var_chaste_interface__slow_delayed_rectifier_K_current_xs_gate__xs; // dimensionless
            const NekDouble var_slow_delayed_rectifier_K_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_slow_delayed_rectifier_K_current__g_Ks = 0.12941176; // nanoS_per_picoF
            const NekDouble var_slow_delayed_rectifier_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_slow_delayed_rectifier_K_current__i_Ks = var_slow_delayed_rectifier_K_current__Cm * var_slow_delayed_rectifier_K_current__g_Ks * pow(var_slow_delayed_rectifier_K_current__xs, 2.0) * (var_slow_delayed_rectifier_K_current__V - var_slow_delayed_rectifier_K_current__E_K); // picoA
            const NekDouble var_L_type_Ca_channel__f = var_chaste_interface__L_type_Ca_channel_f_gate__f; // dimensionless
            const NekDouble var_L_type_Ca_channel__f_Ca = var_chaste_interface__L_type_Ca_channel_f_Ca_gate__f_Ca; // dimensionless
            const NekDouble var_L_type_Ca_channel__d = var_chaste_interface__L_type_Ca_channel_d_gate__d; // dimensionless
            const NekDouble var_L_type_Ca_channel__Cm = var_membrane__Cm; // picoF
            const NekDouble var_L_type_Ca_channel__g_Ca_L = 0.12375; // nanoS_per_picoF
            const NekDouble var_L_type_Ca_channel__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__Cm * var_L_type_Ca_channel__g_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f_Ca * (var_L_type_Ca_channel__V - 65.0); // picoA
            const NekDouble var_sarcolemmal_calcium_pump_current__i_CaP_max = 0.275; // picoA_per_picoF
            const NekDouble var_sarcolemmal_calcium_pump_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_sarcolemmal_calcium_pump_current__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_sarcolemmal_calcium_pump_current__i_CaP = (var_sarcolemmal_calcium_pump_current__Cm * var_sarcolemmal_calcium_pump_current__i_CaP_max * var_sarcolemmal_calcium_pump_current__Ca_i) / (0.0005 + var_sarcolemmal_calcium_pump_current__Ca_i); // picoA
            const NekDouble var_sodium_potassium_pump__Km_Na_i = 10.0; // millimolar
            const NekDouble var_sodium_potassium_pump__i_NaK_max = 0.59933874; // picoA_per_picoF
            const NekDouble var_sodium_potassium_pump__Km_K_o = 1.5; // millimolar
            const NekDouble var_sodium_potassium_pump__Cm = var_membrane__Cm; // picoF
            const NekDouble var_sodium_potassium_pump__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_sodium_potassium_pump__Na_i = var_chaste_interface__intracellular_ion_concentrations__Na_i; // millimolar
            const NekDouble var_sodium_potassium_pump__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_sodium_potassium_pump__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_sodium_potassium_pump__T = var_membrane__T; // kelvin
            const NekDouble var_sodium_potassium_pump__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_sodium_potassium_pump__sigma = (1.0 / 7.0) * (exp(var_sodium_potassium_pump__Na_o / 67.3) - 1.0); // dimensionless
            const NekDouble var_sodium_potassium_pump__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_sodium_potassium_pump__f_NaK = pow(1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__F * var_sodium_potassium_pump__V) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__F) * var_sodium_potassium_pump__V) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))), -1.0); // dimensionless
            const NekDouble var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__Cm * var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__f_NaK * 1.0) / (1.0 + pow(var_sodium_potassium_pump__Km_Na_i / var_sodium_potassium_pump__Na_i, 1.5))) * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_o + var_sodium_potassium_pump__Km_K_o); // picoA
            const NekDouble var_Na_Ca_exchanger_current__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_Na_Ca_exchanger_current__K_mNa = 87.5; // millimolar
            const NekDouble var_Na_Ca_exchanger_current__I_NaCa_max = 1600.0; // picoA_per_picoF
            const NekDouble var_Na_Ca_exchanger_current__T = var_membrane__T; // kelvin
            const NekDouble var_Na_Ca_exchanger_current__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_Na_Ca_exchanger_current__K_sat = 0.1; // dimensionless
            const NekDouble var_Na_Ca_exchanger_current__gamma = 0.35; // dimensionless
            const NekDouble var_standard_ionic_concentrations__Ca_o = 1.8; // millimolar
            const NekDouble var_Na_Ca_exchanger_current__Ca_o = var_standard_ionic_concentrations__Ca_o; // millimolar
            const NekDouble var_Na_Ca_exchanger_current__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_Na_Ca_exchanger_current__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_Na_Ca_exchanger_current__Cm = var_membrane__Cm; // picoF
            const NekDouble var_Na_Ca_exchanger_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_Na_Ca_exchanger_current__K_mCa = 1.38; // millimolar
            const NekDouble var_Na_Ca_exchanger_current__Na_i = var_chaste_interface__intracellular_ion_concentrations__Na_i; // millimolar
            const NekDouble var_Na_Ca_exchanger_current__i_NaCa = (var_Na_Ca_exchanger_current__Cm * var_Na_Ca_exchanger_current__I_NaCa_max * ((exp((var_Na_Ca_exchanger_current__gamma * var_Na_Ca_exchanger_current__F * var_Na_Ca_exchanger_current__V) / (var_Na_Ca_exchanger_current__R * var_Na_Ca_exchanger_current__T)) * pow(var_Na_Ca_exchanger_current__Na_i, 3.0) * var_Na_Ca_exchanger_current__Ca_o) - (exp(((var_Na_Ca_exchanger_current__gamma - 1.0) * var_Na_Ca_exchanger_current__F * var_Na_Ca_exchanger_current__V) / (var_Na_Ca_exchanger_current__R * var_Na_Ca_exchanger_current__T)) * pow(var_Na_Ca_exchanger_current__Na_o, 3.0) * var_Na_Ca_exchanger_current__Ca_i))) / ((pow(var_Na_Ca_exchanger_current__K_mNa, 3.0) + pow(var_Na_Ca_exchanger_current__Na_o, 3.0)) * (var_Na_Ca_exchanger_current__K_mCa + var_Na_Ca_exchanger_current__Ca_o) * (1.0 + (var_Na_Ca_exchanger_current__K_sat * exp(((var_Na_Ca_exchanger_current__gamma - 1.0) * var_Na_Ca_exchanger_current__V * var_Na_Ca_exchanger_current__F) / (var_Na_Ca_exchanger_current__R * var_Na_Ca_exchanger_current__T))))); // picoA
            const NekDouble var_background_currents__Cm = var_membrane__Cm; // picoF
            const NekDouble var_background_currents__E_Na = var_fast_sodium_current__E_Na; // millivolt
            const NekDouble var_background_currents__g_B_Na = 0.0006744375; // nanoS_per_picoF
            const NekDouble var_background_currents__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_background_currents__i_B_Na = var_background_currents__Cm * var_background_currents__g_B_Na * (var_background_currents__V - var_background_currents__E_Na); // picoA
            const NekDouble var_background_currents__R = var_membrane__R; // joule_per_mole_kelvin
            const NekDouble var_background_currents__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_background_currents__Ca_o = var_standard_ionic_concentrations__Ca_o; // millimolar
            const NekDouble var_background_currents__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_background_currents__T = var_membrane__T; // kelvin
            const NekDouble var_background_currents__E_Ca = ((var_background_currents__R * var_background_currents__T) / (2.0 * var_background_currents__F)) * log(var_background_currents__Ca_o / var_background_currents__Ca_i); // millivolt
            const NekDouble var_background_currents__g_B_Ca = 0.001131; // nanoS_per_picoF
            const NekDouble var_background_currents__i_B_Ca = var_background_currents__Cm * var_background_currents__g_B_Ca * (var_background_currents__V - var_background_currents__E_Ca); // picoA
            const NekDouble var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_m_gate__alpha_m = (var_fast_sodium_current_m_gate__V == (-47.13)) ? 3.2 : ((0.32 * (var_fast_sodium_current_m_gate__V + 47.13)) / (1.0 - exp((-0.1) * (var_fast_sodium_current_m_gate__V + 47.13)))); // per_millisecond
            const NekDouble var_fast_sodium_current_m_gate__beta_m = 0.08 * exp((-var_fast_sodium_current_m_gate__V) / 11.0); // per_millisecond
            const NekDouble var_fast_sodium_current_m_gate__tau_m = 1.0 / (var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m); // millisecond
            const NekDouble var_fast_sodium_current_m_gate__m_inf = var_fast_sodium_current_m_gate__alpha_m / (var_fast_sodium_current_m_gate__alpha_m + var_fast_sodium_current_m_gate__beta_m); // dimensionless
            const NekDouble var_fast_sodium_current_m_gate__m = var_fast_sodium_current__m; // dimensionless
            const NekDouble var_fast_sodium_current_m_gate__d_m_d_environment__time = (var_fast_sodium_current_m_gate__m_inf - var_fast_sodium_current_m_gate__m) / var_fast_sodium_current_m_gate__tau_m; // per_millisecond
            const NekDouble var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time = var_fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
            const NekDouble var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((3.56 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.35 * var_fast_sodium_current_h_gate__V))) : (1.0 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1))))); // per_millisecond
            const NekDouble var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.135 * exp((var_fast_sodium_current_h_gate__V + 80.0) / (-6.8))) : 0.0; // per_millisecond
            const NekDouble var_fast_sodium_current_h_gate__h_inf = var_fast_sodium_current_h_gate__alpha_h / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h); // dimensionless
            const NekDouble var_fast_sodium_current_h_gate__tau_h = 1.0 / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h); // millisecond
            const NekDouble var_fast_sodium_current_h_gate__h = var_fast_sodium_current__h; // dimensionless
            const NekDouble var_fast_sodium_current_h_gate__d_h_d_environment__time = (var_fast_sodium_current_h_gate__h_inf - var_fast_sodium_current_h_gate__h) / var_fast_sodium_current_h_gate__tau_h; // per_millisecond
            const NekDouble var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time = var_fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
            const NekDouble var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V; // millivolt
            const NekDouble var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? (((((-127140.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (3.474e-05 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0; // per_millisecond
            const NekDouble var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.1212 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.3 * exp((-2.535e-07) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0)))); // per_millisecond
            const NekDouble var_fast_sodium_current_j_gate__j_inf = var_fast_sodium_current_j_gate__alpha_j / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j); // dimensionless
            const NekDouble var_fast_sodium_current_j_gate__tau_j = 1.0 / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j); // millisecond
            const NekDouble var_fast_sodium_current_j_gate__j = var_fast_sodium_current__j; // dimensionless
            const NekDouble var_fast_sodium_current_j_gate__d_j_d_environment__time = (var_fast_sodium_current_j_gate__j_inf - var_fast_sodium_current_j_gate__j) / var_fast_sodium_current_j_gate__tau_j; // per_millisecond
            const NekDouble var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time = var_fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
            const NekDouble var_transient_outward_K_current__K_Q10 = 3.0; // dimensionless
            const NekDouble var_transient_outward_K_current_oa_gate__V = var_transient_outward_K_current__V; // millivolt
            const NekDouble var_transient_outward_K_current_oa_gate__beta_oa = 0.65 * pow(2.5 + exp(((var_transient_outward_K_current_oa_gate__V - (-10.0)) + 72.0) / 17.0), -1.0); // per_millisecond
            const NekDouble var_transient_outward_K_current_oa_gate__alpha_oa = 0.65 * pow(exp((var_transient_outward_K_current_oa_gate__V - (-10.0)) / (-8.5)) + exp(((var_transient_outward_K_current_oa_gate__V - (-10.0)) - 40.0) / (-59.0)), -1.0); // per_millisecond
            const NekDouble var_transient_outward_K_current_oa_gate__K_Q10 = var_transient_outward_K_current__K_Q10; // dimensionless
            const NekDouble var_transient_outward_K_current_oa_gate__tau_oa = pow(var_transient_outward_K_current_oa_gate__alpha_oa + var_transient_outward_K_current_oa_gate__beta_oa, -1.0) / var_transient_outward_K_current_oa_gate__K_Q10; // millisecond
            const NekDouble var_transient_outward_K_current_oa_gate__oa_infinity = pow(1.0 + exp(((var_transient_outward_K_current_oa_gate__V - (-10.0)) + 10.47) / (-17.54)), -1.0); // dimensionless
            const NekDouble var_transient_outward_K_current_oa_gate__oa = var_transient_outward_K_current__oa; // dimensionless
            const NekDouble var_transient_outward_K_current_oa_gate__d_oa_d_environment__time = (var_transient_outward_K_current_oa_gate__oa_infinity - var_transient_outward_K_current_oa_gate__oa) / var_transient_outward_K_current_oa_gate__tau_oa; // per_millisecond
            const NekDouble var_transient_outward_K_current__transient_outward_K_current_oa_gate__d_oa_d_environment__time = var_transient_outward_K_current_oa_gate__d_oa_d_environment__time; // per_millisecond
            const NekDouble var_transient_outward_K_current_oi_gate__V = var_transient_outward_K_current__V; // millivolt
            const NekDouble var_transient_outward_K_current_oi_gate__beta_oi = pow(35.56 + (1.0 * exp(((var_transient_outward_K_current_oi_gate__V - (-10.0)) - 8.74) / (-7.44))), -1.0); // per_millisecond
            const NekDouble var_transient_outward_K_current_oi_gate__alpha_oi = pow(18.53 + (1.0 * exp(((var_transient_outward_K_current_oi_gate__V - (-10.0)) + 103.7) / 10.95)), -1.0); // per_millisecond
            const NekDouble var_transient_outward_K_current_oi_gate__K_Q10 = var_transient_outward_K_current__K_Q10; // dimensionless
            const NekDouble var_transient_outward_K_current_oi_gate__tau_oi = pow(var_transient_outward_K_current_oi_gate__alpha_oi + var_transient_outward_K_current_oi_gate__beta_oi, -1.0) / var_transient_outward_K_current_oi_gate__K_Q10; // millisecond
            const NekDouble var_transient_outward_K_current_oi_gate__oi_infinity = pow(1.0 + exp(((var_transient_outward_K_current_oi_gate__V - (-10.0)) + 33.1) / 5.3), -1.0); // dimensionless
            const NekDouble var_transient_outward_K_current_oi_gate__oi = var_transient_outward_K_current__oi; // dimensionless
            const NekDouble var_transient_outward_K_current_oi_gate__d_oi_d_environment__time = (var_transient_outward_K_current_oi_gate__oi_infinity - var_transient_outward_K_current_oi_gate__oi) / var_transient_outward_K_current_oi_gate__tau_oi; // per_millisecond
            const NekDouble var_transient_outward_K_current__transient_outward_K_current_oi_gate__d_oi_d_environment__time = var_transient_outward_K_current_oi_gate__d_oi_d_environment__time; // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__K_Q10 = var_transient_outward_K_current__K_Q10; // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua = var_ultrarapid_delayed_rectifier_K_current__ua; // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__V = var_ultrarapid_delayed_rectifier_K_current__V; // millivolt
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__alpha_ua = 0.65 * pow(exp((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) / (-8.5)) + exp(((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) - 40.0) / (-59.0)), -1.0); // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__beta_ua = 0.65 * pow(2.5 + exp(((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) + 72.0) / 17.0), -1.0); // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__K_Q10 = var_ultrarapid_delayed_rectifier_K_current__K_Q10; // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__tau_ua = pow(var_ultrarapid_delayed_rectifier_K_current_ua_gate__alpha_ua + var_ultrarapid_delayed_rectifier_K_current_ua_gate__beta_ua, -1.0) / var_ultrarapid_delayed_rectifier_K_current_ua_gate__K_Q10; // millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua_infinity = pow(1.0 + exp(((var_ultrarapid_delayed_rectifier_K_current_ua_gate__V - (-10.0)) + 20.3) / (-9.6)), -1.0); // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ua_gate__d_ua_d_environment__time = (var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua_infinity - var_ultrarapid_delayed_rectifier_K_current_ua_gate__ua) / var_ultrarapid_delayed_rectifier_K_current_ua_gate__tau_ua; // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__ultrarapid_delayed_rectifier_K_current_ua_gate__d_ua_d_environment__time = var_ultrarapid_delayed_rectifier_K_current_ua_gate__d_ua_d_environment__time; // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui = var_ultrarapid_delayed_rectifier_K_current__ui; // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__V = var_ultrarapid_delayed_rectifier_K_current__V; // millivolt
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__alpha_ui = pow(21.0 + (1.0 * exp(((var_ultrarapid_delayed_rectifier_K_current_ui_gate__V - (-10.0)) - 195.0) / (-28.0))), -1.0); // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__beta_ui = 1.0 / exp(((var_ultrarapid_delayed_rectifier_K_current_ui_gate__V - (-10.0)) - 168.0) / (-16.0)); // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__K_Q10 = var_ultrarapid_delayed_rectifier_K_current__K_Q10; // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__tau_ui = pow(var_ultrarapid_delayed_rectifier_K_current_ui_gate__alpha_ui + var_ultrarapid_delayed_rectifier_K_current_ui_gate__beta_ui, -1.0) / var_ultrarapid_delayed_rectifier_K_current_ui_gate__K_Q10; // millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui_infinity = pow(1.0 + exp(((var_ultrarapid_delayed_rectifier_K_current_ui_gate__V - (-10.0)) - 109.45) / 27.48), -1.0); // dimensionless
            const NekDouble var_ultrarapid_delayed_rectifier_K_current_ui_gate__d_ui_d_environment__time = (var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui_infinity - var_ultrarapid_delayed_rectifier_K_current_ui_gate__ui) / var_ultrarapid_delayed_rectifier_K_current_ui_gate__tau_ui; // per_millisecond
            const NekDouble var_ultrarapid_delayed_rectifier_K_current__ultrarapid_delayed_rectifier_K_current_ui_gate__d_ui_d_environment__time = var_ultrarapid_delayed_rectifier_K_current_ui_gate__d_ui_d_environment__time; // per_millisecond
            const NekDouble var_rapid_delayed_rectifier_K_current_xr_gate__V = var_rapid_delayed_rectifier_K_current__V; // millivolt
            const NekDouble var_rapid_delayed_rectifier_K_current_xr_gate__alpha_xr = (fabs(var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1) < 1e-10) ? 0.0015 : ((0.0003 * (var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1)) / (1.0 - exp((var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1) / (-5.0)))); // per_millisecond
            const NekDouble var_rapid_delayed_rectifier_K_current_xr_gate__beta_xr = (fabs(var_rapid_delayed_rectifier_K_current_xr_gate__V - 3.3328) < 1e-10) ? 0.00037836118 : ((7.3898e-05 * (var_rapid_delayed_rectifier_K_current_xr_gate__V - 3.3328)) / (exp((var_rapid_delayed_rectifier_K_current_xr_gate__V - 3.3328) / 5.1237) - 1.0)); // per_millisecond
            const NekDouble var_rapid_delayed_rectifier_K_current_xr_gate__tau_xr = pow(var_rapid_delayed_rectifier_K_current_xr_gate__alpha_xr + var_rapid_delayed_rectifier_K_current_xr_gate__beta_xr, -1.0); // millisecond
            const NekDouble var_rapid_delayed_rectifier_K_current_xr_gate__xr_infinity = pow(1.0 + exp((var_rapid_delayed_rectifier_K_current_xr_gate__V + 14.1) / (-6.5)), -1.0); // dimensionless
            const NekDouble var_rapid_delayed_rectifier_K_current_xr_gate__xr = var_rapid_delayed_rectifier_K_current__xr; // dimensionless
            const NekDouble var_rapid_delayed_rectifier_K_current_xr_gate__d_xr_d_environment__time = (var_rapid_delayed_rectifier_K_current_xr_gate__xr_infinity - var_rapid_delayed_rectifier_K_current_xr_gate__xr) / var_rapid_delayed_rectifier_K_current_xr_gate__tau_xr; // per_millisecond
            const NekDouble var_rapid_delayed_rectifier_K_current__rapid_delayed_rectifier_K_current_xr_gate__d_xr_d_environment__time = var_rapid_delayed_rectifier_K_current_xr_gate__d_xr_d_environment__time; // per_millisecond
            const NekDouble var_slow_delayed_rectifier_K_current_xs_gate__V = var_slow_delayed_rectifier_K_current__V; // millivolt
            const NekDouble var_slow_delayed_rectifier_K_current_xs_gate__alpha_xs = (fabs(var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) < 1e-10) ? 0.00068 : ((4e-05 * (var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9)) / (1.0 - exp((var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) / (-17.0)))); // per_millisecond
            const NekDouble var_slow_delayed_rectifier_K_current_xs_gate__beta_xs = (fabs(var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) < 1e-10) ? 0.000315 : ((3.5e-05 * (var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9)) / (exp((var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) / 9.0) - 1.0)); // per_millisecond
            const NekDouble var_slow_delayed_rectifier_K_current_xs_gate__tau_xs = 0.5 * pow(var_slow_delayed_rectifier_K_current_xs_gate__alpha_xs + var_slow_delayed_rectifier_K_current_xs_gate__beta_xs, -1.0); // millisecond
            const NekDouble var_slow_delayed_rectifier_K_current_xs_gate__xs_infinity = pow(1.0 + exp((var_slow_delayed_rectifier_K_current_xs_gate__V - 19.9) / (-12.7)), -0.5); // dimensionless
            const NekDouble var_slow_delayed_rectifier_K_current_xs_gate__xs = var_slow_delayed_rectifier_K_current__xs; // dimensionless
            const NekDouble var_slow_delayed_rectifier_K_current_xs_gate__d_xs_d_environment__time = (var_slow_delayed_rectifier_K_current_xs_gate__xs_infinity - var_slow_delayed_rectifier_K_current_xs_gate__xs) / var_slow_delayed_rectifier_K_current_xs_gate__tau_xs; // per_millisecond
            const NekDouble var_slow_delayed_rectifier_K_current__slow_delayed_rectifier_K_current_xs_gate__d_xs_d_environment__time = var_slow_delayed_rectifier_K_current_xs_gate__d_xs_d_environment__time; // per_millisecond
            const NekDouble var_L_type_Ca_channel__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_L_type_Ca_channel_d_gate__d = var_L_type_Ca_channel__d; // dimensionless
            const NekDouble var_L_type_Ca_channel_d_gate__V = var_L_type_Ca_channel__V; // millivolt
            const NekDouble var_L_type_Ca_channel_d_gate__d_infinity = pow(1.0 + exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-8.0)), -1.0); // dimensionless
            const NekDouble var_L_type_Ca_channel_d_gate__tau_d = (fabs(var_L_type_Ca_channel_d_gate__V + 10.0) < 1e-10) ? (4.579 / (1.0 + exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-6.24)))) : ((1.0 - exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-6.24))) / (0.035 * (var_L_type_Ca_channel_d_gate__V + 10.0) * (1.0 + exp((var_L_type_Ca_channel_d_gate__V + 10.0) / (-6.24))))); // millisecond
            const NekDouble var_L_type_Ca_channel_d_gate__d_d_d_environment__time = (var_L_type_Ca_channel_d_gate__d_infinity - var_L_type_Ca_channel_d_gate__d) / var_L_type_Ca_channel_d_gate__tau_d; // per_millisecond
            const NekDouble var_L_type_Ca_channel__L_type_Ca_channel_d_gate__d_d_d_environment__time = var_L_type_Ca_channel_d_gate__d_d_d_environment__time; // per_millisecond
            const NekDouble var_L_type_Ca_channel_f_gate__f = var_L_type_Ca_channel__f; // dimensionless
            const NekDouble var_L_type_Ca_channel_f_gate__V = var_L_type_Ca_channel__V; // millivolt
            const NekDouble var_L_type_Ca_channel_f_gate__f_infinity = exp((-(var_L_type_Ca_channel_f_gate__V + 28.0)) / 6.9) / (1.0 + exp((-(var_L_type_Ca_channel_f_gate__V + 28.0)) / 6.9)); // dimensionless
            const NekDouble var_L_type_Ca_channel_f_gate__tau_f = 9.0 * pow((0.0197 * exp((-pow(0.0337, 2.0)) * pow(var_L_type_Ca_channel_f_gate__V + 10.0, 2.0))) + 0.02, -1.0); // millisecond
            const NekDouble var_L_type_Ca_channel_f_gate__d_f_d_environment__time = (var_L_type_Ca_channel_f_gate__f_infinity - var_L_type_Ca_channel_f_gate__f) / var_L_type_Ca_channel_f_gate__tau_f; // per_millisecond
            const NekDouble var_L_type_Ca_channel__L_type_Ca_channel_f_gate__d_f_d_environment__time = var_L_type_Ca_channel_f_gate__d_f_d_environment__time; // per_millisecond
            const NekDouble var_L_type_Ca_channel_f_Ca_gate__f_Ca = var_L_type_Ca_channel__f_Ca; // dimensionless
            const NekDouble var_L_type_Ca_channel_f_Ca_gate__tau_f_Ca = 2.0; // millisecond
            const NekDouble var_L_type_Ca_channel_f_Ca_gate__Ca_i = var_L_type_Ca_channel__Ca_i; // millimolar
            const NekDouble var_L_type_Ca_channel_f_Ca_gate__f_Ca_infinity = pow(1.0 + (var_L_type_Ca_channel_f_Ca_gate__Ca_i / 0.00035), -1.0); // dimensionless
            const NekDouble var_L_type_Ca_channel_f_Ca_gate__d_f_Ca_d_environment__time = (var_L_type_Ca_channel_f_Ca_gate__f_Ca_infinity - var_L_type_Ca_channel_f_Ca_gate__f_Ca) / var_L_type_Ca_channel_f_Ca_gate__tau_f_Ca; // per_millisecond
            const NekDouble var_L_type_Ca_channel__L_type_Ca_channel_f_Ca_gate__d_f_Ca_d_environment__time = var_L_type_Ca_channel_f_Ca_gate__d_f_Ca_d_environment__time; // per_millisecond
            const NekDouble var_background_currents__g_B_K = 0.0; // nanoS_per_picoF
            const NekDouble var_background_currents__E_K = var_time_independent_potassium_current__E_K; // millivolt
            const NekDouble var_background_currents__i_B_K = var_background_currents__Cm * var_background_currents__g_B_K * (var_background_currents__V - var_background_currents__E_K); // picoA
            const NekDouble var_Ca_release_current_from_JSR__K_rel = 30.0; // per_millisecond
            const NekDouble var_Ca_release_current_from_JSR__u = var_chaste_interface__Ca_release_current_from_JSR_u_gate__u; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR__w = var_chaste_interface__Ca_release_current_from_JSR_w_gate__w; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR__Ca_rel = var_chaste_interface__intracellular_ion_concentrations__Ca_rel; // millimolar
            const NekDouble var_Ca_release_current_from_JSR__v = var_chaste_interface__Ca_release_current_from_JSR_v_gate__v; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_Ca_release_current_from_JSR__i_rel = var_Ca_release_current_from_JSR__K_rel * pow(var_Ca_release_current_from_JSR__u, 2.0) * var_Ca_release_current_from_JSR__v * var_Ca_release_current_from_JSR__w * (var_Ca_release_current_from_JSR__Ca_rel - var_Ca_release_current_from_JSR__Ca_i); // millimolar_per_millisecond
            const NekDouble var_Ca_release_current_from_JSR__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_intracellular_ion_concentrations__V_cell = 20100.0; // micrometre_3
            const NekDouble var_intracellular_ion_concentrations__V_rel = 0.0048 * var_intracellular_ion_concentrations__V_cell; // micrometre_3
            const NekDouble var_Ca_release_current_from_JSR__V_rel = var_intracellular_ion_concentrations__V_rel; // micrometre_3
            const NekDouble var_Ca_release_current_from_JSR__i_NaCa = var_Na_Ca_exchanger_current__i_NaCa; // picoA
            const NekDouble var_Ca_release_current_from_JSR__i_Ca_L = var_L_type_Ca_channel__i_Ca_L; // picoA
            const NekDouble var_Ca_release_current_from_JSR__Fn = 1000.0 * ((1e-15 * var_Ca_release_current_from_JSR__V_rel * var_Ca_release_current_from_JSR__i_rel) - ((1e-15 / (2.0 * var_Ca_release_current_from_JSR__F)) * ((0.5 * var_Ca_release_current_from_JSR__i_Ca_L) - (0.2 * var_Ca_release_current_from_JSR__i_NaCa)))); // dimensionless
            const NekDouble var_Ca_release_current_from_JSR__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_Ca_release_current_from_JSR_u_gate__u = var_Ca_release_current_from_JSR__u; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_u_gate__tau_u = 8.0; // millisecond
            const NekDouble var_Ca_release_current_from_JSR_u_gate__Fn = var_Ca_release_current_from_JSR__Fn; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_u_gate__u_infinity = pow(1.0 + exp((-(var_Ca_release_current_from_JSR_u_gate__Fn - 3.4175e-13)) / 1.367e-15), -1.0); // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_u_gate__d_u_d_environment__time = (var_Ca_release_current_from_JSR_u_gate__u_infinity - var_Ca_release_current_from_JSR_u_gate__u) / var_Ca_release_current_from_JSR_u_gate__tau_u; // per_millisecond
            const NekDouble var_Ca_release_current_from_JSR__Ca_release_current_from_JSR_u_gate__d_u_d_environment__time = var_Ca_release_current_from_JSR_u_gate__d_u_d_environment__time; // per_millisecond
            const NekDouble var_Ca_release_current_from_JSR_v_gate__Fn = var_Ca_release_current_from_JSR__Fn; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_v_gate__tau_v = 1.91 + (2.09 * pow(1.0 + exp((-(var_Ca_release_current_from_JSR_v_gate__Fn - 3.4175e-13)) / 1.367e-15), -1.0)); // millisecond
            const NekDouble var_Ca_release_current_from_JSR_v_gate__v = var_Ca_release_current_from_JSR__v; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_v_gate__v_infinity = 1.0 - pow(1.0 + exp((-(var_Ca_release_current_from_JSR_v_gate__Fn - 6.835e-14)) / 1.367e-15), -1.0); // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_v_gate__d_v_d_environment__time = (var_Ca_release_current_from_JSR_v_gate__v_infinity - var_Ca_release_current_from_JSR_v_gate__v) / var_Ca_release_current_from_JSR_v_gate__tau_v; // per_millisecond
            const NekDouble var_Ca_release_current_from_JSR__Ca_release_current_from_JSR_v_gate__d_v_d_environment__time = var_Ca_release_current_from_JSR_v_gate__d_v_d_environment__time; // per_millisecond
            const NekDouble var_Ca_release_current_from_JSR_w_gate__w = var_Ca_release_current_from_JSR__w; // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_w_gate__V = var_Ca_release_current_from_JSR__V; // millivolt
            const NekDouble var_Ca_release_current_from_JSR_w_gate__tau_w = (fabs(var_Ca_release_current_from_JSR_w_gate__V - 7.9) < 1e-10) ? ((6.0 * 0.2) / 1.3) : ((6.0 * (1.0 - exp((-(var_Ca_release_current_from_JSR_w_gate__V - 7.9)) / 5.0))) / ((1.0 + (0.3 * exp((-(var_Ca_release_current_from_JSR_w_gate__V - 7.9)) / 5.0))) * 1.0 * (var_Ca_release_current_from_JSR_w_gate__V - 7.9))); // millisecond
            const NekDouble var_Ca_release_current_from_JSR_w_gate__w_infinity = 1.0 - pow(1.0 + exp((-(var_Ca_release_current_from_JSR_w_gate__V - 40.0)) / 17.0), -1.0); // dimensionless
            const NekDouble var_Ca_release_current_from_JSR_w_gate__d_w_d_environment__time = (var_Ca_release_current_from_JSR_w_gate__w_infinity - var_Ca_release_current_from_JSR_w_gate__w) / var_Ca_release_current_from_JSR_w_gate__tau_w; // per_millisecond
            const NekDouble var_Ca_release_current_from_JSR__Ca_release_current_from_JSR_w_gate__d_w_d_environment__time = var_Ca_release_current_from_JSR_w_gate__d_w_d_environment__time; // per_millisecond
            const NekDouble var_transfer_current_from_NSR_to_JSR__Ca_up = var_chaste_interface__intracellular_ion_concentrations__Ca_up; // millimolar
            const NekDouble var_transfer_current_from_NSR_to_JSR__Ca_rel = var_chaste_interface__intracellular_ion_concentrations__Ca_rel; // millimolar
            const NekDouble var_transfer_current_from_NSR_to_JSR__tau_tr = 180.0; // millisecond
            const NekDouble var_transfer_current_from_NSR_to_JSR__i_tr = (var_transfer_current_from_NSR_to_JSR__Ca_up - var_transfer_current_from_NSR_to_JSR__Ca_rel) / var_transfer_current_from_NSR_to_JSR__tau_tr; // millimolar_per_millisecond
            const NekDouble var_Ca_uptake_current_by_the_NSR__I_up_max = 0.005; // millimolar_per_millisecond
            const NekDouble var_Ca_uptake_current_by_the_NSR__K_up = 0.00092; // millimolar
            const NekDouble var_Ca_uptake_current_by_the_NSR__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_Ca_uptake_current_by_the_NSR__i_up = var_Ca_uptake_current_by_the_NSR__I_up_max / (1.0 + (var_Ca_uptake_current_by_the_NSR__K_up / var_Ca_uptake_current_by_the_NSR__Ca_i)); // millimolar_per_millisecond
            const NekDouble var_Ca_leak_current_by_the_NSR__Ca_up_max = 15.0; // millimolar
            const NekDouble var_Ca_leak_current_by_the_NSR__I_up_max = var_Ca_uptake_current_by_the_NSR__I_up_max; // millimolar_per_millisecond
            const NekDouble var_Ca_leak_current_by_the_NSR__Ca_up = var_chaste_interface__intracellular_ion_concentrations__Ca_up; // millimolar
            const NekDouble var_Ca_leak_current_by_the_NSR__i_up_leak = (var_Ca_leak_current_by_the_NSR__I_up_max * var_Ca_leak_current_by_the_NSR__Ca_up) / var_Ca_leak_current_by_the_NSR__Ca_up_max; // millimolar_per_millisecond
            const NekDouble var_Ca_buffers__CMDN_max = 0.05; // millimolar
            const NekDouble var_Ca_buffers__TRPN_max = 0.07; // millimolar
            const NekDouble var_Ca_buffers__CSQN_max = 10.0; // millimolar
            const NekDouble var_Ca_buffers__Km_CMDN = 0.00238; // millimolar
            const NekDouble var_Ca_buffers__Km_TRPN = 0.0005; // millimolar
            const NekDouble var_Ca_buffers__Km_CSQN = 0.8; // millimolar
            const NekDouble var_intracellular_ion_concentrations__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_intracellular_ion_concentrations__Ca_rel = var_chaste_interface__intracellular_ion_concentrations__Ca_rel; // millimolar
            const NekDouble var_intracellular_ion_concentrations__V_i = var_intracellular_ion_concentrations__V_cell * 0.68; // micrometre_3
            const NekDouble var_intracellular_ion_concentrations__V_up = 0.0552 * var_intracellular_ion_concentrations__V_cell; // micrometre_3
            const NekDouble var_intracellular_ion_concentrations__i_Ca_L = var_L_type_Ca_channel__i_Ca_L; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_rel = var_Ca_release_current_from_JSR__i_rel; // millimolar_per_millisecond
            const NekDouble var_intracellular_ion_concentrations__i_B_Ca = var_background_currents__i_B_Ca; // picoA
            const NekDouble var_intracellular_ion_concentrations__F = var_membrane__F; // coulomb_per_millimole
            const NekDouble var_intracellular_ion_concentrations__i_NaCa = var_Na_Ca_exchanger_current__i_NaCa; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_up_leak = var_Ca_leak_current_by_the_NSR__i_up_leak; // millimolar_per_millisecond
            const NekDouble var_intracellular_ion_concentrations__i_up = var_Ca_uptake_current_by_the_NSR__i_up; // millimolar_per_millisecond
            const NekDouble var_intracellular_ion_concentrations__i_CaP = var_sarcolemmal_calcium_pump_current__i_CaP; // picoA
            const NekDouble var_intracellular_ion_concentrations__B1 = (((2.0 * var_intracellular_ion_concentrations__i_NaCa) - (var_intracellular_ion_concentrations__i_CaP + var_intracellular_ion_concentrations__i_Ca_L + var_intracellular_ion_concentrations__i_B_Ca)) / (2.0 * var_intracellular_ion_concentrations__V_i * var_intracellular_ion_concentrations__F)) + (((var_intracellular_ion_concentrations__V_up * (var_intracellular_ion_concentrations__i_up_leak - var_intracellular_ion_concentrations__i_up)) + (var_intracellular_ion_concentrations__i_rel * var_intracellular_ion_concentrations__V_rel)) / var_intracellular_ion_concentrations__V_i); // millimolar_per_millisecond
            const NekDouble var_intracellular_ion_concentrations__Km_TRPN = var_Ca_buffers__Km_TRPN; // millimolar
            const NekDouble var_intracellular_ion_concentrations__CMDN_max = var_Ca_buffers__CMDN_max; // millimolar
            const NekDouble var_intracellular_ion_concentrations__Km_CMDN = var_Ca_buffers__Km_CMDN; // millimolar
            const NekDouble var_intracellular_ion_concentrations__TRPN_max = var_Ca_buffers__TRPN_max; // millimolar
            const NekDouble var_intracellular_ion_concentrations__B2 = 1.0 + ((var_intracellular_ion_concentrations__TRPN_max * var_intracellular_ion_concentrations__Km_TRPN) / pow(var_intracellular_ion_concentrations__Ca_i + var_intracellular_ion_concentrations__Km_TRPN, 2.0)) + ((var_intracellular_ion_concentrations__CMDN_max * var_intracellular_ion_concentrations__Km_CMDN) / pow(var_intracellular_ion_concentrations__Ca_i + var_intracellular_ion_concentrations__Km_CMDN, 2.0)); // dimensionless
            const NekDouble var_intracellular_ion_concentrations__i_NaK = var_sodium_potassium_pump__i_NaK; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_B_Na = var_background_currents__i_B_Na; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_Na = var_fast_sodium_current__i_Na; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_K1 = var_time_independent_potassium_current__i_K1; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_to = var_transient_outward_K_current__i_to; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_Kur = var_ultrarapid_delayed_rectifier_K_current__i_Kur; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_Kr = var_rapid_delayed_rectifier_K_current__i_Kr; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_Ks = var_slow_delayed_rectifier_K_current__i_Ks; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_B_K = var_background_currents__i_B_K; // picoA
            const NekDouble var_intracellular_ion_concentrations__i_tr = var_transfer_current_from_NSR_to_JSR__i_tr; // millimolar_per_millisecond
            const NekDouble var_intracellular_ion_concentrations__Km_CSQN = var_Ca_buffers__Km_CSQN; // millimolar
            const NekDouble var_intracellular_ion_concentrations__CSQN_max = var_Ca_buffers__CSQN_max; // millimolar
            const NekDouble var_intracellular_ion_concentrations__d_Na_i_d_environment__time = (((-3.0) * var_intracellular_ion_concentrations__i_NaK) - ((3.0 * var_intracellular_ion_concentrations__i_NaCa) + var_intracellular_ion_concentrations__i_B_Na + var_intracellular_ion_concentrations__i_Na)) / (var_intracellular_ion_concentrations__V_i * var_intracellular_ion_concentrations__F); // 'millimole per litre per millisecond'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_i_d_environment__time = var_intracellular_ion_concentrations__B1 / var_intracellular_ion_concentrations__B2; // 'millimole per litre per millisecond'
            const NekDouble var_intracellular_ion_concentrations__d_K_i_d_environment__time = ((2.0 * var_intracellular_ion_concentrations__i_NaK) - (var_intracellular_ion_concentrations__i_K1 + var_intracellular_ion_concentrations__i_to + var_intracellular_ion_concentrations__i_Kur + var_intracellular_ion_concentrations__i_Kr + var_intracellular_ion_concentrations__i_Ks + var_intracellular_ion_concentrations__i_B_K)) / (var_intracellular_ion_concentrations__V_i * var_intracellular_ion_concentrations__F); // 'millimole per litre per millisecond'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_rel_d_environment__time = (var_intracellular_ion_concentrations__i_tr - var_intracellular_ion_concentrations__i_rel) * pow(1.0 + ((var_intracellular_ion_concentrations__CSQN_max * var_intracellular_ion_concentrations__Km_CSQN) / pow(var_intracellular_ion_concentrations__Ca_rel + var_intracellular_ion_concentrations__Km_CSQN, 2.0)), -1.0); // 'millimole per litre per millisecond'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_up_d_environment__time = var_intracellular_ion_concentrations__i_up - (var_intracellular_ion_concentrations__i_up_leak + ((var_intracellular_ion_concentrations__i_tr * var_intracellular_ion_concentrations__V_rel) / var_intracellular_ion_concentrations__V_up)); // 'millimole per litre per millisecond'
            const NekDouble var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time = var_fast_sodium_current__fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time = var_fast_sodium_current__fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time = var_fast_sodium_current__fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__transient_outward_K_current_oa_gate__d_oa_d_environment__time = var_transient_outward_K_current__transient_outward_K_current_oa_gate__d_oa_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__transient_outward_K_current_oi_gate__d_oi_d_environment__time = var_transient_outward_K_current__transient_outward_K_current_oi_gate__d_oi_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ua_gate__d_ua_d_environment__time = var_ultrarapid_delayed_rectifier_K_current__ultrarapid_delayed_rectifier_K_current_ua_gate__d_ua_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ui_gate__d_ui_d_environment__time = var_ultrarapid_delayed_rectifier_K_current__ultrarapid_delayed_rectifier_K_current_ui_gate__d_ui_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__rapid_delayed_rectifier_K_current_xr_gate__d_xr_d_environment__time = var_rapid_delayed_rectifier_K_current__rapid_delayed_rectifier_K_current_xr_gate__d_xr_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__slow_delayed_rectifier_K_current_xs_gate__d_xs_d_environment__time = var_slow_delayed_rectifier_K_current__slow_delayed_rectifier_K_current_xs_gate__d_xs_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__L_type_Ca_channel_d_gate__d_d_d_environment__time = var_L_type_Ca_channel__L_type_Ca_channel_d_gate__d_d_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__L_type_Ca_channel_f_gate__d_f_d_environment__time = var_L_type_Ca_channel__L_type_Ca_channel_f_gate__d_f_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__L_type_Ca_channel_f_Ca_gate__d_f_Ca_d_environment__time = var_L_type_Ca_channel__L_type_Ca_channel_f_Ca_gate__d_f_Ca_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__Ca_release_current_from_JSR_u_gate__d_u_d_environment__time = var_Ca_release_current_from_JSR__Ca_release_current_from_JSR_u_gate__d_u_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__Ca_release_current_from_JSR_v_gate__d_v_d_environment__time = var_Ca_release_current_from_JSR__Ca_release_current_from_JSR_v_gate__d_v_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__Ca_release_current_from_JSR_w_gate__d_w_d_environment__time = var_Ca_release_current_from_JSR__Ca_release_current_from_JSR_w_gate__d_w_d_environment__time; // per_millisecond
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Na_i_d_environment__time = var_intracellular_ion_concentrations__d_Na_i_d_environment__time; // millimolar_per_millisecond
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_i_d_environment__time = var_intracellular_ion_concentrations__d_Ca_i_d_environment__time; // millimolar_per_millisecond
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_K_i_d_environment__time = var_intracellular_ion_concentrations__d_K_i_d_environment__time; // millimolar_per_millisecond
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_rel_d_environment__time = var_intracellular_ion_concentrations__d_Ca_rel_d_environment__time; // millimolar_per_millisecond
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_up_d_environment__time = var_intracellular_ion_concentrations__d_Ca_up_d_environment__time; // millimolar_per_millisecond
            const NekDouble d_dt_chaste_interface__fast_sodium_current_m_gate__m = var_chaste_interface__fast_sodium_current_m_gate__d_m_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__fast_sodium_current_h_gate__h = var_chaste_interface__fast_sodium_current_h_gate__d_h_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__fast_sodium_current_j_gate__j = var_chaste_interface__fast_sodium_current_j_gate__d_j_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__transient_outward_K_current_oa_gate__oa = var_chaste_interface__transient_outward_K_current_oa_gate__d_oa_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__transient_outward_K_current_oi_gate__oi = var_chaste_interface__transient_outward_K_current_oi_gate__d_oi_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__ultrarapid_delayed_rectifier_K_current_ua_gate__ua = var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ua_gate__d_ua_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__ultrarapid_delayed_rectifier_K_current_ui_gate__ui = var_chaste_interface__ultrarapid_delayed_rectifier_K_current_ui_gate__d_ui_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__rapid_delayed_rectifier_K_current_xr_gate__xr = var_chaste_interface__rapid_delayed_rectifier_K_current_xr_gate__d_xr_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__slow_delayed_rectifier_K_current_xs_gate__xs = var_chaste_interface__slow_delayed_rectifier_K_current_xs_gate__d_xs_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__L_type_Ca_channel_d_gate__d = var_chaste_interface__L_type_Ca_channel_d_gate__d_d_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__L_type_Ca_channel_f_gate__f = var_chaste_interface__L_type_Ca_channel_f_gate__d_f_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__L_type_Ca_channel_f_Ca_gate__f_Ca = var_chaste_interface__L_type_Ca_channel_f_Ca_gate__d_f_Ca_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__Ca_release_current_from_JSR_u_gate__u = var_chaste_interface__Ca_release_current_from_JSR_u_gate__d_u_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__Ca_release_current_from_JSR_v_gate__v = var_chaste_interface__Ca_release_current_from_JSR_v_gate__d_v_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__Ca_release_current_from_JSR_w_gate__w = var_chaste_interface__Ca_release_current_from_JSR_w_gate__d_w_d_environment__time; // per_millisecond
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Na_i = var_chaste_interface__intracellular_ion_concentrations__d_Na_i_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_i = var_chaste_interface__intracellular_ion_concentrations__d_Ca_i_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__K_i = var_chaste_interface__intracellular_ion_concentrations__d_K_i_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_rel = var_chaste_interface__intracellular_ion_concentrations__d_Ca_rel_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_up = var_chaste_interface__intracellular_ion_concentrations__d_Ca_up_d_environment__time; // 'millimole per litre per millisecond'
            
            const NekDouble var_chaste_interface__membrane__i_st = 0.0;
            const NekDouble var_membrane__i_st_converter = var_chaste_interface__membrane__i_st; // uA_per_cm2
            const NekDouble var_membrane__chaste_interface__chaste_membrane_capacitance = 1.0; // uF_per_cm2
            const NekDouble var_membrane__i_st = (var_membrane__i_st_converter * var_membrane__Cm) / var_membrane__chaste_interface__chaste_membrane_capacitance; // picoA
            const NekDouble var_membrane__i_Na = var_fast_sodium_current__i_Na; // picoA
            const NekDouble var_membrane__i_K1 = var_time_independent_potassium_current__i_K1; // picoA
            const NekDouble var_membrane__i_to = var_transient_outward_K_current__i_to; // picoA
            const NekDouble var_membrane__i_Kur = var_ultrarapid_delayed_rectifier_K_current__i_Kur; // picoA
            const NekDouble var_membrane__i_Kr = var_rapid_delayed_rectifier_K_current__i_Kr; // picoA
            const NekDouble var_membrane__i_Ks = var_slow_delayed_rectifier_K_current__i_Ks; // picoA
            const NekDouble var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L; // picoA
            const NekDouble var_membrane__i_CaP = var_sarcolemmal_calcium_pump_current__i_CaP; // picoA
            const NekDouble var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK; // picoA
            const NekDouble var_membrane__i_NaCa = var_Na_Ca_exchanger_current__i_NaCa; // picoA
            const NekDouble var_membrane__i_B_Na = var_background_currents__i_B_Na; // picoA
            const NekDouble var_membrane__i_B_Ca = var_background_currents__i_B_Ca; // picoA
            const NekDouble var_membrane__d_V_d_environment__time = (-(var_membrane__i_Na + var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kur + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_B_Na + var_membrane__i_B_Ca + var_membrane__i_NaK + var_membrane__i_CaP + var_membrane__i_NaCa + var_membrane__i_Ca_L + var_membrane__i_st)) / var_membrane__Cm; // 'millivolt per millisecond'
            const NekDouble var_chaste_interface__membrane__d_V_d_environment__time = var_membrane__d_V_d_environment__time; // ___units_1
            d_dt_chaste_interface__membrane__V = var_chaste_interface__membrane__d_V_d_environment__time; // 'millivolt per millisecond'
            outarray[0][i] = outarray[0][i]/var_membrane__Cm;
            outarray[0][i] += d_dt_chaste_interface__membrane__V;
            outarray[1][i] = 0.0;
            outarray[2][i] = d_dt_chaste_interface__fast_sodium_current_m_gate__m;
            outarray[3][i] = d_dt_chaste_interface__fast_sodium_current_h_gate__h;
            outarray[4][i] = d_dt_chaste_interface__fast_sodium_current_j_gate__j;
            outarray[5][i] = d_dt_chaste_interface__transient_outward_K_current_oa_gate__oa;
            outarray[6][i] = d_dt_chaste_interface__transient_outward_K_current_oi_gate__oi;
            outarray[7][i] = d_dt_chaste_interface__ultrarapid_delayed_rectifier_K_current_ua_gate__ua;
            outarray[8][i] = d_dt_chaste_interface__ultrarapid_delayed_rectifier_K_current_ui_gate__ui;
            outarray[9][i] = d_dt_chaste_interface__rapid_delayed_rectifier_K_current_xr_gate__xr;
            outarray[10][i] = d_dt_chaste_interface__slow_delayed_rectifier_K_current_xs_gate__xs;
            outarray[11][i] = d_dt_chaste_interface__L_type_Ca_channel_d_gate__d;
            outarray[12][i] = d_dt_chaste_interface__L_type_Ca_channel_f_gate__f;
            outarray[13][i] = d_dt_chaste_interface__L_type_Ca_channel_f_Ca_gate__f_Ca;
            outarray[14][i] = d_dt_chaste_interface__Ca_release_current_from_JSR_u_gate__u;
            outarray[15][i] = d_dt_chaste_interface__Ca_release_current_from_JSR_v_gate__v;
            outarray[16][i] = d_dt_chaste_interface__Ca_release_current_from_JSR_w_gate__w;
            outarray[17][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Na_i;
            outarray[18][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_i;
            outarray[19][i] = d_dt_chaste_interface__intracellular_ion_concentrations__K_i;
            outarray[20][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_rel;
            outarray[21][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_up;
        }
        
    }

    /**
    *
    */
    void CourtemancheRamirezNattel98::v_PrintSummary(std::ostream &out)
    {
        out << "	Cell model      : CourtemancheRamirezNattel98" << std::endl;
    }

}
