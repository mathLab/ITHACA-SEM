///////////////////////////////////////////////////////////////////////////////
//
// File PanditGilesDemir03.cpp
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
// Description: Pandit-Giles-Demir 2003 cell model
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <CardiacEPSolver/CellModels/PanditGilesDemir03.h>
namespace Nektar
{
    std::string PanditGilesDemir03::className
              = GetCellModelFactory().RegisterCreatorFunction(
                        "PanditGilesDemir03",
                        PanditGilesDemir03::create,
                         "Pandit-Giles-Demir 2003 cell model.");


    /**
    *
    */
    PanditGilesDemir03::PanditGilesDemir03(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const int nq): CellModel(pSession, nq)
    {
        m_nq   = nq;
    }


    void PanditGilesDemir03::v_Update(
                     const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                           Array<OneD,        Array<OneD, NekDouble> >&outarray,
                     const NekDouble time)
    {
        int nvariables  = inarray.size();
        int nq = m_nq;
        for (unsigned int i = 0; i < nq; ++i)
        {
            // Inputs:
            // Time units: millisecond
            NekDouble var_chaste_interface__membrane__V = inarray[0][i];
            // Units: millivolt; Initial value: -80.50146
            NekDouble var_chaste_interface__sodium_current_m_gate__m = inarray[2][i];
            // Units: dimensionless; Initial value: 0.004164108
            NekDouble var_chaste_interface__sodium_current_h_gate__h = inarray[3][i];
            // Units: dimensionless; Initial value: 0.6735613
            NekDouble var_chaste_interface__sodium_current_j_gate__j = inarray[4][i];
            // Units: dimensionless; Initial value: 0.6729362
            NekDouble var_chaste_interface__L_type_Ca_channel_d_gate__d = inarray[5][i];
            // Units: dimensionless; Initial value: 0.000002171081
            NekDouble var_chaste_interface__L_type_Ca_channel_f_11_gate__f_11 = inarray[6][i];
            // Units: dimensionless; Initial value: 0.9999529
            NekDouble var_chaste_interface__L_type_Ca_channel_f_12_gate__f_12 = inarray[7][i];
            // Units: dimensionless; Initial value: 0.9999529
            NekDouble var_chaste_interface__L_type_Ca_channel_Ca_inact_gate__Ca_inact = inarray[8][i];
            // Units: dimensionless; Initial value: 0.9913102
            NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__r = inarray[9][i];
            // Units: dimensionless; Initial value: 0.002191519
            NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__s = inarray[10][i];
            // Units: dimensionless; Initial value: 0.9842542
            NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__s_slow = inarray[11][i];
            // Units: dimensionless; Initial value: 0.6421196
            NekDouble var_chaste_interface__steady_state_outward_K_current_r_ss_gate__r_ss = inarray[12][i];
            // Units: dimensionless; Initial value: 0.002907171
            NekDouble var_chaste_interface__steady_state_outward_K_current_s_ss_gate__s_ss = inarray[13][i];
            // Units: dimensionless; Initial value: 0.3142767
            NekDouble var_chaste_interface__hyperpolarisation_activated_current_y_gate__y = inarray[14][i];
            // Units: dimensionless; Initial value: 0.003578708
            NekDouble var_chaste_interface__SR_Ca_release_channel__P_O1 = inarray[15][i];
            // Units: dimensionless; Initial value: 0.0004327548
            NekDouble var_chaste_interface__SR_Ca_release_channel__P_O2 = inarray[16][i];
            // Units: dimensionless; Initial value: 0.000000000606254
            NekDouble var_chaste_interface__SR_Ca_release_channel__P_C1 = inarray[17][i];
            // Units: dimensionless; Initial value: 0.6348229
            NekDouble var_chaste_interface__SR_Ca_release_channel__P_C2 = inarray[18][i];
            // Units: dimensionless; Initial value: 0.3647471
            NekDouble var_chaste_interface__intracellular_and_SR_Ca_fluxes__HTRPNCa = inarray[19][i];
            // Units: millimolar; Initial value: 1.394301e-1
            NekDouble var_chaste_interface__intracellular_and_SR_Ca_fluxes__LTRPNCa = inarray[20][i];
            // Units: millimolar; Initial value: 5.1619e-3
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Na_i = inarray[21][i];
            // Units: millimolar; Initial value: 10.73519
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_i = inarray[22][i];
            // Units: millimolar; Initial value: 0.00007901351
            NekDouble var_chaste_interface__intracellular_ion_concentrations__K_i = inarray[23][i];
            // Units: millimolar; Initial value: 139.2751
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_ss = inarray[24][i];
            // Units: millimolar; Initial value: 0.00008737212
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_JSR = inarray[25][i];
            // Units: millimolar; Initial value: 0.06607948
            NekDouble var_chaste_interface__intracellular_ion_concentrations__Ca_NSR = inarray[26][i];
            // Units: millimolar; Initial value: 0.06600742


            // Mathematics
            NekDouble d_dt_chaste_interface__membrane__V;
            const NekDouble var_membrane__R = 8314.5; // millijoule_per_mole_kelvin
            const NekDouble var_membrane__T = 295.0; // kelvin
            const NekDouble var_membrane__F = 96487.0; // coulomb_per_mole
            const NekDouble var_sodium_current__F = var_membrane__F; // coulomb_per_mole
            const NekDouble var_standard_ionic_concentrations__Na_o = 140.0; // millimolar
            const NekDouble var_sodium_current__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_sodium_current__T = var_membrane__T; // kelvin
            const NekDouble var_sodium_current__R = var_membrane__R; // millijoule_per_mole_kelvin
            const NekDouble var_sodium_current__Na_i = var_chaste_interface__intracellular_ion_concentrations__Na_i; // millimolar
            const NekDouble var_sodium_current__E_Na = ((var_sodium_current__R * var_sodium_current__T) / var_sodium_current__F) * log(var_sodium_current__Na_o / var_sodium_current__Na_i); // millivolt
            const NekDouble var_sodium_current__m = var_chaste_interface__sodium_current_m_gate__m; // dimensionless
            const NekDouble var_sodium_current__j = var_chaste_interface__sodium_current_j_gate__j; // dimensionless
            const NekDouble var_sodium_current__h = var_chaste_interface__sodium_current_h_gate__h; // dimensionless
            const NekDouble var_sodium_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_sodium_current__g_Na = 1.064; // microS
            const NekDouble var_sodium_current__i_Na = var_sodium_current__g_Na * pow(var_sodium_current__m, 3.0) * var_sodium_current__h * var_sodium_current__j * (var_sodium_current__V - var_sodium_current__E_Na); // nanoA
            const NekDouble var_L_type_Ca_channel__d = var_chaste_interface__L_type_Ca_channel_d_gate__d; // dimensionless
            const NekDouble var_L_type_Ca_channel__E_Ca_L = 65.0; // millivolt
            const NekDouble var_L_type_Ca_channel__g_Ca_L_normal = 0.0341; // microS
            const NekDouble var_membrane__Diabetes = 0.0; // dimensionless
            const NekDouble var_L_type_Ca_channel__Diabetes = var_membrane__Diabetes; // dimensionless
            const NekDouble var_L_type_Ca_channel__g_Ca_L = (var_L_type_Ca_channel__Diabetes == 0.0) ? var_L_type_Ca_channel__g_Ca_L_normal : (0.76 * var_L_type_Ca_channel__g_Ca_L_normal); // microS
            const NekDouble var_L_type_Ca_channel__f_12 = var_chaste_interface__L_type_Ca_channel_f_12_gate__f_12; // dimensionless
            const NekDouble var_L_type_Ca_channel__f_11 = var_chaste_interface__L_type_Ca_channel_f_11_gate__f_11; // dimensionless
            const NekDouble var_L_type_Ca_channel__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_L_type_Ca_channel__Ca_inact = var_chaste_interface__L_type_Ca_channel_Ca_inact_gate__Ca_inact; // dimensionless
            const NekDouble var_L_type_Ca_channel__i_Ca_L = var_L_type_Ca_channel__g_Ca_L * var_L_type_Ca_channel__d * (((0.9 + (var_L_type_Ca_channel__Ca_inact / 10.0)) * var_L_type_Ca_channel__f_11) + ((0.1 - (var_L_type_Ca_channel__Ca_inact / 10.0)) * var_L_type_Ca_channel__f_12)) * (var_L_type_Ca_channel__V - var_L_type_Ca_channel__E_Ca_L); // nanoA
            const NekDouble var_Ca_independent_transient_outward_K_current__Diabetes = var_membrane__Diabetes; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current__b = (var_Ca_independent_transient_outward_K_current__Diabetes == 0.0) ? 0.114 : 0.31; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current__a = 1.0 - var_Ca_independent_transient_outward_K_current__b; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current__g_t_normal = 0.04375; // microS
            const NekDouble var_Ca_independent_transient_outward_K_current__g_t = (var_Ca_independent_transient_outward_K_current__Diabetes == 0.0) ? var_Ca_independent_transient_outward_K_current__g_t_normal : (0.68 * var_Ca_independent_transient_outward_K_current__g_t_normal); // microS
            const NekDouble var_Ca_independent_transient_outward_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_Ca_independent_transient_outward_K_current__r = var_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__r; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current__s = var_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__s; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current__F = var_membrane__F; // coulomb_per_mole
            const NekDouble var_Ca_independent_transient_outward_K_current__T = var_membrane__T; // kelvin
            const NekDouble var_Ca_independent_transient_outward_K_current__K_i = var_chaste_interface__intracellular_ion_concentrations__K_i; // millimolar
            const NekDouble var_Ca_independent_transient_outward_K_current__R = var_membrane__R; // millijoule_per_mole_kelvin
            const NekDouble var_standard_ionic_concentrations__K_o = 5.4; // millimolar
            const NekDouble var_Ca_independent_transient_outward_K_current__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_Ca_independent_transient_outward_K_current__E_K = ((var_Ca_independent_transient_outward_K_current__R * var_Ca_independent_transient_outward_K_current__T) / var_Ca_independent_transient_outward_K_current__F) * log(var_Ca_independent_transient_outward_K_current__K_o / var_Ca_independent_transient_outward_K_current__K_i); // millivolt
            const NekDouble var_Ca_independent_transient_outward_K_current__s_slow = var_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__s_slow; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current__i_t = var_Ca_independent_transient_outward_K_current__g_t * var_Ca_independent_transient_outward_K_current__r * ((var_Ca_independent_transient_outward_K_current__a * var_Ca_independent_transient_outward_K_current__s) + (var_Ca_independent_transient_outward_K_current__b * var_Ca_independent_transient_outward_K_current__s_slow)) * (var_Ca_independent_transient_outward_K_current__V - var_Ca_independent_transient_outward_K_current__E_K); // nanoA
            const NekDouble var_steady_state_outward_K_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_steady_state_outward_K_current__s_ss = var_chaste_interface__steady_state_outward_K_current_s_ss_gate__s_ss; // dimensionless
            const NekDouble var_steady_state_outward_K_current__E_K = var_Ca_independent_transient_outward_K_current__E_K; // millivolt
            const NekDouble var_steady_state_outward_K_current__g_ss_normal = 0.0077; // microS
            const NekDouble var_steady_state_outward_K_current__Diabetes = var_membrane__Diabetes; // dimensionless
            const NekDouble var_steady_state_outward_K_current__g_ss = (var_steady_state_outward_K_current__Diabetes == 0.0) ? var_steady_state_outward_K_current__g_ss_normal : (0.77 * var_steady_state_outward_K_current__g_ss_normal); // microS
            const NekDouble var_steady_state_outward_K_current__r_ss = var_chaste_interface__steady_state_outward_K_current_r_ss_gate__r_ss; // dimensionless
            const NekDouble var_steady_state_outward_K_current__i_ss = var_steady_state_outward_K_current__g_ss * var_steady_state_outward_K_current__r_ss * var_steady_state_outward_K_current__s_ss * (var_steady_state_outward_K_current__V - var_steady_state_outward_K_current__E_K); // nanoA
            const NekDouble var_hyperpolarisation_activated_current__f_Na = 0.2; // dimensionless
            const NekDouble var_hyperpolarisation_activated_current__E_Na = var_sodium_current__E_Na; // millivolt
            const NekDouble var_hyperpolarisation_activated_current__y = var_chaste_interface__hyperpolarisation_activated_current_y_gate__y; // dimensionless
            const NekDouble var_hyperpolarisation_activated_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_hyperpolarisation_activated_current__g_f = 0.00145; // microS
            const NekDouble var_hyperpolarisation_activated_current__i_f_Na = var_hyperpolarisation_activated_current__g_f * var_hyperpolarisation_activated_current__y * var_hyperpolarisation_activated_current__f_Na * (var_hyperpolarisation_activated_current__V - var_hyperpolarisation_activated_current__E_Na); // nanoA
            const NekDouble var_hyperpolarisation_activated_current__f_K = 1.0 - var_hyperpolarisation_activated_current__f_Na; // dimensionless
            const NekDouble var_hyperpolarisation_activated_current__E_K = var_Ca_independent_transient_outward_K_current__E_K; // millivolt
            const NekDouble var_hyperpolarisation_activated_current__i_f_K = var_hyperpolarisation_activated_current__g_f * var_hyperpolarisation_activated_current__y * var_hyperpolarisation_activated_current__f_K * (var_hyperpolarisation_activated_current__V - var_hyperpolarisation_activated_current__E_K); // nanoA
            const NekDouble var_inward_rectifier__T = var_membrane__T; // kelvin
            const NekDouble var_inward_rectifier__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_inward_rectifier__g_K1 = 0.024; // microS
            const NekDouble var_inward_rectifier__R = var_membrane__R; // millijoule_per_mole_kelvin
            const NekDouble var_inward_rectifier__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_inward_rectifier__F = var_membrane__F; // coulomb_per_mole
            const NekDouble var_inward_rectifier__E_K = var_Ca_independent_transient_outward_K_current__E_K; // millivolt
            const NekDouble var_inward_rectifier__i_K1 = ((((48.0 / (exp((var_inward_rectifier__V + 37.0) / 25.0) + exp((var_inward_rectifier__V + 37.0) / (-25.0)))) + 10.0) * 0.0001) / (1.0 + exp((var_inward_rectifier__V - (var_inward_rectifier__E_K + 76.77)) / (-17.0)))) + ((var_inward_rectifier__g_K1 * (var_inward_rectifier__V - (var_inward_rectifier__E_K + 1.73))) / ((1.0 + exp((1.613 * var_inward_rectifier__F * (var_inward_rectifier__V - (var_inward_rectifier__E_K + 1.73))) / (var_inward_rectifier__R * var_inward_rectifier__T))) * (1.0 + exp((var_inward_rectifier__K_o - 0.9988) / (-0.124))))); // nanoA
            const NekDouble var_background_currents__R = var_membrane__R; // millijoule_per_mole_kelvin
            const NekDouble var_background_currents__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_standard_ionic_concentrations__Ca_o = 1.2; // millimolar
            const NekDouble var_background_currents__Ca_o = var_standard_ionic_concentrations__Ca_o; // millimolar
            const NekDouble var_background_currents__F = var_membrane__F; // coulomb_per_mole
            const NekDouble var_background_currents__T = var_membrane__T; // kelvin
            const NekDouble var_background_currents__E_Ca = ((0.5 * var_background_currents__R * var_background_currents__T) / var_background_currents__F) * log(var_background_currents__Ca_o / var_background_currents__Ca_i); // millivolt
            const NekDouble var_background_currents__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_background_currents__g_B_Ca_normal = 3.24e-05; // microS
            const NekDouble var_background_currents__Diabetes = var_membrane__Diabetes; // dimensionless
            const NekDouble var_background_currents__g_B_Ca = (var_background_currents__Diabetes == 0.0) ? var_background_currents__g_B_Ca_normal : (0.5 * var_background_currents__g_B_Ca_normal); // microS
            const NekDouble var_background_currents__i_B_Ca = var_background_currents__g_B_Ca * (var_background_currents__V - var_background_currents__E_Ca); // nanoA
            const NekDouble var_background_currents__g_B_K = 0.000138; // microS
            const NekDouble var_background_currents__E_K = var_Ca_independent_transient_outward_K_current__E_K; // millivolt
            const NekDouble var_background_currents__i_B_K = var_background_currents__g_B_K * (var_background_currents__V - var_background_currents__E_K); // nanoA
            const NekDouble var_background_currents__E_Na = var_sodium_current__E_Na; // millivolt
            const NekDouble var_background_currents__g_B_Na_normal = 8.015e-05; // microS
            const NekDouble var_background_currents__g_B_Na = (var_background_currents__Diabetes == 0.0) ? var_background_currents__g_B_Na_normal : (1.25 * var_background_currents__g_B_Na_normal); // microS
            const NekDouble var_background_currents__i_B_Na = var_background_currents__g_B_Na * (var_background_currents__V - var_background_currents__E_Na); // nanoA
            const NekDouble var_sodium_potassium_pump__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_sodium_potassium_pump__sigma = (exp(var_sodium_potassium_pump__Na_o / 67.3) - 1.0) / 7.0; // dimensionless
            const NekDouble var_sodium_potassium_pump__K_m_Na = 10.0; // millimolar
            const NekDouble var_sodium_potassium_pump__K_m_K = 1.5; // millimolar
            const NekDouble var_sodium_potassium_pump__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_sodium_potassium_pump__T = var_membrane__T; // kelvin
            const NekDouble var_sodium_potassium_pump__R = var_membrane__R; // millijoule_per_mole_kelvin
            const NekDouble var_sodium_potassium_pump__Na_i = var_chaste_interface__intracellular_ion_concentrations__Na_i; // millimolar
            const NekDouble var_sodium_potassium_pump__F = var_membrane__F; // coulomb_per_mole
            const NekDouble var_sodium_potassium_pump__i_NaK_max_normal = 0.08; // nanoA
            const NekDouble var_sodium_potassium_pump__Diabetes = var_membrane__Diabetes; // dimensionless
            const NekDouble var_sodium_potassium_pump__i_NaK_max = (var_sodium_potassium_pump__Diabetes == 0.0) ? var_sodium_potassium_pump__i_NaK_max_normal : (0.63 * var_sodium_potassium_pump__i_NaK_max_normal); // nanoA
            const NekDouble var_sodium_potassium_pump__K_o = var_standard_ionic_concentrations__K_o; // millimolar
            const NekDouble var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump__V * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))) + (0.0365 * var_sodium_potassium_pump__sigma * exp(((-var_sodium_potassium_pump__V) * var_sodium_potassium_pump__F) / (var_sodium_potassium_pump__R * var_sodium_potassium_pump__T))))) * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_o + var_sodium_potassium_pump__K_m_K)) / (1.0 + pow(var_sodium_potassium_pump__K_m_Na / var_sodium_potassium_pump__Na_i, 1.5)); // nanoA
            const NekDouble var_Na_Ca_ion_exchanger_current__Na_o = var_standard_ionic_concentrations__Na_o; // millimolar
            const NekDouble var_Na_Ca_ion_exchanger_current__Na_i = var_chaste_interface__intracellular_ion_concentrations__Na_i; // millimolar
            const NekDouble var_Na_Ca_ion_exchanger_current__gamma_NaCa = 0.5; // dimensionless
            const NekDouble var_Na_Ca_ion_exchanger_current__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_Na_Ca_ion_exchanger_current__Ca_o = var_standard_ionic_concentrations__Ca_o; // millimolar
            const NekDouble var_Na_Ca_ion_exchanger_current__K_NaCa = 9.984e-06; // nanoA_millimolar_4
            const NekDouble var_Na_Ca_ion_exchanger_current__d_NaCa = 0.0001; // millimolar_4
            const NekDouble var_Na_Ca_ion_exchanger_current__V = var_chaste_interface__membrane__V; // millivolt
            const NekDouble var_Na_Ca_ion_exchanger_current__i_NaCa = (var_Na_Ca_ion_exchanger_current__K_NaCa * ((pow(var_Na_Ca_ion_exchanger_current__Na_i, 3.0) * var_Na_Ca_ion_exchanger_current__Ca_o * exp(0.03743 * var_Na_Ca_ion_exchanger_current__V * var_Na_Ca_ion_exchanger_current__gamma_NaCa)) - (pow(var_Na_Ca_ion_exchanger_current__Na_o, 3.0) * var_Na_Ca_ion_exchanger_current__Ca_i * exp(0.03743 * var_Na_Ca_ion_exchanger_current__V * (var_Na_Ca_ion_exchanger_current__gamma_NaCa - 1.0))))) / (1.0 + (var_Na_Ca_ion_exchanger_current__d_NaCa * ((var_Na_Ca_ion_exchanger_current__Ca_i * pow(var_Na_Ca_ion_exchanger_current__Na_o, 3.0)) + (var_Na_Ca_ion_exchanger_current__Ca_o * pow(var_Na_Ca_ion_exchanger_current__Na_i, 3.0))))); // nanoA
            const NekDouble var_sarcolemmal_calcium_pump_current__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_sarcolemmal_calcium_pump_current__i_Ca_P_max = 0.004; // nanoA
            const NekDouble var_sarcolemmal_calcium_pump_current__i_Ca_P = (var_sarcolemmal_calcium_pump_current__i_Ca_P_max * var_sarcolemmal_calcium_pump_current__Ca_i) / (var_sarcolemmal_calcium_pump_current__Ca_i + 0.0004); // nanoA
            const NekDouble var_sodium_current_m_gate__V = var_sodium_current__V; // millivolt
            const NekDouble var_sodium_current_m_gate__m_infinity = 1.0 / (1.0 + exp((var_sodium_current_m_gate__V + 45.0) / (-6.5))); // dimensionless
            const NekDouble var_sodium_current_m_gate__tau_m = 0.00136 / (((0.32 * (var_sodium_current_m_gate__V + 47.13)) / (1.0 - exp((-0.1) * (var_sodium_current_m_gate__V + 47.13)))) + (0.08 * exp((-var_sodium_current_m_gate__V) / 11.0))); // second
            const NekDouble var_sodium_current_m_gate__m = var_sodium_current__m; // dimensionless
            const NekDouble var_sodium_current_m_gate__d_m_d_environment__time = (var_sodium_current_m_gate__m_infinity - var_sodium_current_m_gate__m) / var_sodium_current_m_gate__tau_m; // per_second
            const NekDouble var_sodium_current__sodium_current_m_gate__d_m_d_environment__time = var_sodium_current_m_gate__d_m_d_environment__time; // per_second
            const NekDouble var_sodium_current_h_gate__V = var_sodium_current__V; // millivolt
            const NekDouble var_sodium_current_h_gate__h_infinity = 1.0 / (1.0 + exp((var_sodium_current_h_gate__V + 76.1) / 6.07)); // dimensionless
            const NekDouble var_sodium_current_h_gate__h = var_sodium_current__h; // dimensionless
            const NekDouble var_sodium_current_h_gate__tau_h = (var_sodium_current_h_gate__V >= (-40.0)) ? (0.0004537 * (1.0 + exp((-(var_sodium_current_h_gate__V + 10.66)) / 11.1))) : (0.00349 / ((0.135 * exp((-(var_sodium_current_h_gate__V + 80.0)) / 6.8)) + (3.56 * exp(0.079 * var_sodium_current_h_gate__V)) + (310000.0 * exp(0.35 * var_sodium_current_h_gate__V)))); // second
            const NekDouble var_sodium_current_h_gate__d_h_d_environment__time = (var_sodium_current_h_gate__h_infinity - var_sodium_current_h_gate__h) / var_sodium_current_h_gate__tau_h; // per_second
            const NekDouble var_sodium_current__sodium_current_h_gate__d_h_d_environment__time = var_sodium_current_h_gate__d_h_d_environment__time; // per_second
            const NekDouble var_sodium_current_j_gate__V = var_sodium_current__V; // millivolt
            const NekDouble var_sodium_current_j_gate__j_infinity = 1.0 / (1.0 + exp((var_sodium_current_j_gate__V + 76.1) / 6.07)); // dimensionless
            const NekDouble var_sodium_current_j_gate__j = var_sodium_current__j; // dimensionless
            const NekDouble var_sodium_current_j_gate__tau_j = (var_sodium_current_j_gate__V >= (-40.0)) ? ((0.01163 * (1.0 + exp((-0.1) * (var_sodium_current_j_gate__V + 32.0)))) / exp((-2.535e-07) * var_sodium_current_j_gate__V)) : (0.00349 / ((((var_sodium_current_j_gate__V + 37.78) / (1.0 + exp(0.311 * (var_sodium_current_j_gate__V + 79.23)))) * (((-127140.0) * exp(0.2444 * var_sodium_current_j_gate__V)) - (3.474e-05 * exp((-0.04391) * var_sodium_current_j_gate__V)))) + ((0.1212 * exp((-0.01052) * var_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_sodium_current_j_gate__V + 40.14)))))); // second
            const NekDouble var_sodium_current_j_gate__d_j_d_environment__time = (var_sodium_current_j_gate__j_infinity - var_sodium_current_j_gate__j) / var_sodium_current_j_gate__tau_j; // per_second
            const NekDouble var_sodium_current__sodium_current_j_gate__d_j_d_environment__time = var_sodium_current_j_gate__d_j_d_environment__time; // per_second
            const NekDouble var_L_type_Ca_channel__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_L_type_Ca_channel_d_gate__d = var_L_type_Ca_channel__d; // dimensionless
            const NekDouble var_L_type_Ca_channel_d_gate__V = var_L_type_Ca_channel__V; // millivolt
            const NekDouble var_L_type_Ca_channel_d_gate__d_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_channel_d_gate__V + 15.3) / (-5.0))); // dimensionless
            const NekDouble var_L_type_Ca_channel_d_gate__tau_d = (0.00305 * exp((-0.0045) * pow(var_L_type_Ca_channel_d_gate__V + 7.0, 2.0))) + (0.00105 * exp((-0.002) * pow(var_L_type_Ca_channel_d_gate__V - 18.0, 2.0))) + 0.00025; // second
            const NekDouble var_L_type_Ca_channel_d_gate__d_d_d_environment__time = (var_L_type_Ca_channel_d_gate__d_infinity - var_L_type_Ca_channel_d_gate__d) / var_L_type_Ca_channel_d_gate__tau_d; // per_second
            const NekDouble var_L_type_Ca_channel__L_type_Ca_channel_d_gate__d_d_d_environment__time = var_L_type_Ca_channel_d_gate__d_d_d_environment__time; // per_second
            const NekDouble var_L_type_Ca_channel_f_11_gate__f_11 = var_L_type_Ca_channel__f_11; // dimensionless
            const NekDouble var_L_type_Ca_channel_f_11_gate__Diabetes = var_L_type_Ca_channel__Diabetes; // dimensionless
            const NekDouble var_L_type_Ca_channel_f_11_gate__V = var_L_type_Ca_channel__V; // millivolt
            const NekDouble var_L_type_Ca_channel_f_11_gate__tau_f_11 = (var_L_type_Ca_channel_f_11_gate__Diabetes == 0.0) ? ((0.105 * exp(-pow((var_L_type_Ca_channel_f_11_gate__V + 45.0) / 12.0, 2.0))) + (0.04 / (1.0 + exp(((-var_L_type_Ca_channel_f_11_gate__V) + 25.0) / 25.0))) + (0.015 / (1.0 + exp((var_L_type_Ca_channel_f_11_gate__V + 75.0) / 25.0))) + 0.0017) : ((0.9 * 0.105 * exp(-pow((var_L_type_Ca_channel_f_11_gate__V + 45.0) / 12.0, 2.0))) + (0.04 / (1.0 + exp(((-var_L_type_Ca_channel_f_11_gate__V) + 25.0) / 25.0))) + (0.015 / (1.0 + exp((var_L_type_Ca_channel_f_11_gate__V + 75.0) / 25.0))) + 0.0017); // second
            const NekDouble var_L_type_Ca_channel_f_11_gate__f_11_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_channel_f_11_gate__V + 26.7) / 5.4)); // dimensionless
            const NekDouble var_L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time = (var_L_type_Ca_channel_f_11_gate__f_11_infinity - var_L_type_Ca_channel_f_11_gate__f_11) / var_L_type_Ca_channel_f_11_gate__tau_f_11; // per_second
            const NekDouble var_L_type_Ca_channel__L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time = var_L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time; // per_second
            const NekDouble var_L_type_Ca_channel_f_12_gate__f_12 = var_L_type_Ca_channel__f_12; // dimensionless
            const NekDouble var_L_type_Ca_channel_f_12_gate__V = var_L_type_Ca_channel__V; // millivolt
            const NekDouble var_L_type_Ca_channel_f_12_gate__f_12_infinity = 1.0 / (1.0 + exp((var_L_type_Ca_channel_f_12_gate__V + 26.7) / 5.4)); // dimensionless
            const NekDouble var_L_type_Ca_channel_f_12_gate__tau_f_12 = (0.041 * exp(-pow((var_L_type_Ca_channel_f_12_gate__V + 47.0) / 12.0, 2.0))) + (0.08 / (1.0 + exp((var_L_type_Ca_channel_f_12_gate__V + 55.0) / (-5.0)))) + (0.015 / (1.0 + exp((var_L_type_Ca_channel_f_12_gate__V + 75.0) / 25.0))) + 0.0017; // second
            const NekDouble var_L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time = (var_L_type_Ca_channel_f_12_gate__f_12_infinity - var_L_type_Ca_channel_f_12_gate__f_12) / var_L_type_Ca_channel_f_12_gate__tau_f_12; // per_second
            const NekDouble var_L_type_Ca_channel__L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time = var_L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time; // per_second
            const NekDouble var_L_type_Ca_channel_Ca_inact_gate__Ca_ss = var_L_type_Ca_channel__Ca_ss; // millimolar
            const NekDouble var_L_type_Ca_channel_Ca_inact_gate__Ca_inact_infinity = 1.0 / (1.0 + (var_L_type_Ca_channel_Ca_inact_gate__Ca_ss / 0.01)); // dimensionless
            const NekDouble var_L_type_Ca_channel_Ca_inact_gate__tau_Ca_inact = 0.009; // second
            const NekDouble var_L_type_Ca_channel_Ca_inact_gate__Ca_inact = var_L_type_Ca_channel__Ca_inact; // dimensionless
            const NekDouble var_L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time = (var_L_type_Ca_channel_Ca_inact_gate__Ca_inact_infinity - var_L_type_Ca_channel_Ca_inact_gate__Ca_inact) / var_L_type_Ca_channel_Ca_inact_gate__tau_Ca_inact; // per_second
            const NekDouble var_L_type_Ca_channel__L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time = var_L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time; // per_second
            const NekDouble var_Ca_independent_transient_outward_K_current_r_gate__V = var_Ca_independent_transient_outward_K_current__V; // millivolt
            const NekDouble var_Ca_independent_transient_outward_K_current_r_gate__r_infinity = 1.0 / (1.0 + exp((var_Ca_independent_transient_outward_K_current_r_gate__V + 10.6) / (-11.42))); // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current_r_gate__r = var_Ca_independent_transient_outward_K_current__r; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current_r_gate__tau_r = 1.0 / ((45.16 * exp(0.03577 * (var_Ca_independent_transient_outward_K_current_r_gate__V + 50.0))) + (98.9 * exp((-0.1) * (var_Ca_independent_transient_outward_K_current_r_gate__V + 38.0)))); // second
            const NekDouble var_Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time = (var_Ca_independent_transient_outward_K_current_r_gate__r_infinity - var_Ca_independent_transient_outward_K_current_r_gate__r) / var_Ca_independent_transient_outward_K_current_r_gate__tau_r; // per_second
            const NekDouble var_Ca_independent_transient_outward_K_current__Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time = var_Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time; // per_second
            const NekDouble var_Ca_independent_transient_outward_K_current_s_gate__V = var_Ca_independent_transient_outward_K_current__V; // millivolt
            const NekDouble var_Ca_independent_transient_outward_K_current_s_gate__s_infinity = 1.0 / (1.0 + exp((var_Ca_independent_transient_outward_K_current_s_gate__V + 45.3) / 6.8841)); // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current_s_gate__tau_s = (0.35 * exp(-pow((var_Ca_independent_transient_outward_K_current_s_gate__V + 70.0) / 15.0, 2.0))) + 0.035; // second
            const NekDouble var_Ca_independent_transient_outward_K_current_s_gate__s = var_Ca_independent_transient_outward_K_current__s; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time = (var_Ca_independent_transient_outward_K_current_s_gate__s_infinity - var_Ca_independent_transient_outward_K_current_s_gate__s) / var_Ca_independent_transient_outward_K_current_s_gate__tau_s; // per_second
            const NekDouble var_Ca_independent_transient_outward_K_current__Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time = var_Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time; // per_second
            const NekDouble var_Ca_independent_transient_outward_K_current_s_slow_gate__V = var_Ca_independent_transient_outward_K_current__V; // millivolt
            const NekDouble var_Ca_independent_transient_outward_K_current_s_slow_gate__s_slow_infinity = 1.0 / (1.0 + exp((var_Ca_independent_transient_outward_K_current_s_slow_gate__V + 45.3) / 6.8841)); // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current_s_slow_gate__s_slow = var_Ca_independent_transient_outward_K_current__s_slow; // dimensionless
            const NekDouble var_Ca_independent_transient_outward_K_current_s_slow_gate__tau_s_slow = (3.7 * exp(-pow((var_Ca_independent_transient_outward_K_current_s_slow_gate__V + 70.0) / 30.0, 2.0))) + 0.035; // second
            const NekDouble var_Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time = (var_Ca_independent_transient_outward_K_current_s_slow_gate__s_slow_infinity - var_Ca_independent_transient_outward_K_current_s_slow_gate__s_slow) / var_Ca_independent_transient_outward_K_current_s_slow_gate__tau_s_slow; // per_second
            const NekDouble var_Ca_independent_transient_outward_K_current__Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time = var_Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time; // per_second
            const NekDouble var_steady_state_outward_K_current_r_ss_gate__V = var_steady_state_outward_K_current__V; // millivolt
            const NekDouble var_steady_state_outward_K_current_r_ss_gate__r_ss_infinity = 1.0 / (1.0 + exp((var_steady_state_outward_K_current_r_ss_gate__V + 11.5) / (-11.82))); // dimensionless
            const NekDouble var_steady_state_outward_K_current_r_ss_gate__r_ss = var_steady_state_outward_K_current__r_ss; // dimensionless
            const NekDouble var_steady_state_outward_K_current_r_ss_gate__tau_r_ss = 10.0 / ((45.16 * exp(0.03577 * (var_steady_state_outward_K_current_r_ss_gate__V + 50.0))) + (98.9 * exp((-0.1) * (var_steady_state_outward_K_current_r_ss_gate__V + 38.0)))); // second
            const NekDouble var_steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time = (var_steady_state_outward_K_current_r_ss_gate__r_ss_infinity - var_steady_state_outward_K_current_r_ss_gate__r_ss) / var_steady_state_outward_K_current_r_ss_gate__tau_r_ss; // per_second
            const NekDouble var_steady_state_outward_K_current__steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time = var_steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time; // per_second
            const NekDouble var_steady_state_outward_K_current_s_ss_gate__s_ss = var_steady_state_outward_K_current__s_ss; // dimensionless
            const NekDouble var_steady_state_outward_K_current_s_ss_gate__tau_s_ss = 2.1; // second
            const NekDouble var_steady_state_outward_K_current_s_ss_gate__V = var_steady_state_outward_K_current__V; // millivolt
            const NekDouble var_steady_state_outward_K_current_s_ss_gate__s_ss_infinity = 1.0 / (1.0 + exp((var_steady_state_outward_K_current_s_ss_gate__V + 87.5) / 10.3)); // dimensionless
            const NekDouble var_steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time = (var_steady_state_outward_K_current_s_ss_gate__s_ss_infinity - var_steady_state_outward_K_current_s_ss_gate__s_ss) / var_steady_state_outward_K_current_s_ss_gate__tau_s_ss; // per_second
            const NekDouble var_steady_state_outward_K_current__steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time = var_steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time; // per_second
            const NekDouble var_hyperpolarisation_activated_current_y_gate__V = var_hyperpolarisation_activated_current__V; // millivolt
            const NekDouble var_hyperpolarisation_activated_current_y_gate__tau_y = 1.0 / ((0.11885 * exp((var_hyperpolarisation_activated_current_y_gate__V + 80.0) / 28.37)) + (0.5623 * exp((var_hyperpolarisation_activated_current_y_gate__V + 80.0) / (-14.19)))); // second
            const NekDouble var_hyperpolarisation_activated_current_y_gate__y = var_hyperpolarisation_activated_current__y; // dimensionless
            const NekDouble var_hyperpolarisation_activated_current_y_gate__y_infinity = 1.0 / (1.0 + exp((var_hyperpolarisation_activated_current_y_gate__V + 138.6) / 10.48)); // dimensionless
            const NekDouble var_hyperpolarisation_activated_current_y_gate__d_y_d_environment__time = (var_hyperpolarisation_activated_current_y_gate__y_infinity - var_hyperpolarisation_activated_current_y_gate__y) / var_hyperpolarisation_activated_current_y_gate__tau_y; // per_second
            const NekDouble var_hyperpolarisation_activated_current__hyperpolarisation_activated_current_y_gate__d_y_d_environment__time = var_hyperpolarisation_activated_current_y_gate__d_y_d_environment__time; // per_second
            const NekDouble var_SR_Ca_release_channel__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_SR_Ca_release_channel__v1 = 1800.0; // per_second
            const NekDouble var_SR_Ca_release_channel__P_O2 = var_chaste_interface__SR_Ca_release_channel__P_O2; // dimensionless
            const NekDouble var_SR_Ca_release_channel__P_O1 = var_chaste_interface__SR_Ca_release_channel__P_O1; // dimensionless
            const NekDouble var_SR_Ca_release_channel__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__Ca_JSR; // millimolar
            const NekDouble var_SR_Ca_release_channel__J_rel = var_SR_Ca_release_channel__v1 * (var_SR_Ca_release_channel__P_O1 + var_SR_Ca_release_channel__P_O2) * (var_SR_Ca_release_channel__Ca_JSR - var_SR_Ca_release_channel__Ca_ss); // millimolar_per_second
            const NekDouble var_SR_Ca_release_channel__k_a_plus = 1.215e+13; // per_second
            const NekDouble var_SR_Ca_release_channel__k_a_minus = 576.0; // per_second
            const NekDouble var_SR_Ca_release_channel__k_b_plus = 4050000000.0; // per_second
            const NekDouble var_SR_Ca_release_channel__k_b_minus = 1930.0; // per_second
            const NekDouble var_SR_Ca_release_channel__k_c_plus = 100.0; // per_second
            const NekDouble var_SR_Ca_release_channel__k_c_minus = 0.8; // per_second
            const NekDouble var_SR_Ca_release_channel__P_C1 = var_chaste_interface__SR_Ca_release_channel__P_C1; // dimensionless
            const NekDouble var_SR_Ca_release_channel__P_C2 = var_chaste_interface__SR_Ca_release_channel__P_C2; // dimensionless
            const NekDouble var_SR_Ca_release_channel__n = 4.0; // dimensionless
            const NekDouble var_SR_Ca_release_channel__m = 3.0; // dimensionless
            const NekDouble var_SR_Ca_release_channel__d_P_O1_d_environment__time = ((var_SR_Ca_release_channel__k_a_plus * pow(var_SR_Ca_release_channel__Ca_ss / 1.0, var_SR_Ca_release_channel__n) * var_SR_Ca_release_channel__P_C1) - ((var_SR_Ca_release_channel__k_a_minus * var_SR_Ca_release_channel__P_O1) + (var_SR_Ca_release_channel__k_b_plus * pow(var_SR_Ca_release_channel__Ca_ss / 1.0, var_SR_Ca_release_channel__m) * var_SR_Ca_release_channel__P_O1) + (var_SR_Ca_release_channel__k_c_plus * var_SR_Ca_release_channel__P_O1))) + (var_SR_Ca_release_channel__k_b_minus * var_SR_Ca_release_channel__P_O2) + (var_SR_Ca_release_channel__k_c_minus * var_SR_Ca_release_channel__P_C2); // per_second
            const NekDouble var_SR_Ca_release_channel__d_P_O2_d_environment__time = (var_SR_Ca_release_channel__k_b_plus * pow(var_SR_Ca_release_channel__Ca_ss / 1.0, var_SR_Ca_release_channel__m) * var_SR_Ca_release_channel__P_O1) - (var_SR_Ca_release_channel__k_b_minus * var_SR_Ca_release_channel__P_O2); // per_second
            const NekDouble var_SR_Ca_release_channel__d_P_C1_d_environment__time = ((-var_SR_Ca_release_channel__k_a_plus) * pow(var_SR_Ca_release_channel__Ca_ss / 1.0, var_SR_Ca_release_channel__n) * var_SR_Ca_release_channel__P_C1) + (var_SR_Ca_release_channel__k_a_minus * var_SR_Ca_release_channel__P_O1); // per_second
            const NekDouble var_SR_Ca_release_channel__d_P_C2_d_environment__time = (var_SR_Ca_release_channel__k_c_plus * var_SR_Ca_release_channel__P_O1) - (var_SR_Ca_release_channel__k_c_minus * var_SR_Ca_release_channel__P_C2); // per_second
            const NekDouble var_SERCA2a_pump__N_fb = 1.2; // dimensionless
            const NekDouble var_SERCA2a_pump__K_fb = 0.000168; // millimolar
            const NekDouble var_SERCA2a_pump__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_SERCA2a_pump__fb = pow(var_SERCA2a_pump__Ca_i / var_SERCA2a_pump__K_fb, var_SERCA2a_pump__N_fb); // dimensionless
            const NekDouble var_SERCA2a_pump__Diabetes = var_membrane__Diabetes; // dimensionless
            const NekDouble var_SERCA2a_pump__Vmaxf_normal = 0.04; // millimolar_per_second
            const NekDouble var_SERCA2a_pump__Vmaxf = (var_SERCA2a_pump__Diabetes == 0.0) ? var_SERCA2a_pump__Vmaxf_normal : (0.8 * var_SERCA2a_pump__Vmaxf_normal); // millimolar_per_second
            const NekDouble var_SERCA2a_pump__K_SR_normal = 1.0; // dimensionless
            const NekDouble var_SERCA2a_pump__K_SR = (var_SERCA2a_pump__Diabetes == 0.0) ? var_SERCA2a_pump__K_SR_normal : (0.55 * var_SERCA2a_pump__K_SR_normal); // dimensionless
            const NekDouble var_SERCA2a_pump__Ca_NSR = var_chaste_interface__intracellular_ion_concentrations__Ca_NSR; // millimolar
            const NekDouble var_SERCA2a_pump__K_rb = 3.29; // millimolar
            const NekDouble var_SERCA2a_pump__N_rb = 1.0; // dimensionless
            const NekDouble var_SERCA2a_pump__rb = pow(var_SERCA2a_pump__Ca_NSR / var_SERCA2a_pump__K_rb, var_SERCA2a_pump__N_rb); // dimensionless
            const NekDouble var_SERCA2a_pump__Vmaxr = 0.9; // millimolar_per_second
            const NekDouble var_SERCA2a_pump__J_up = (var_SERCA2a_pump__K_SR * ((var_SERCA2a_pump__Vmaxf * var_SERCA2a_pump__fb) - (var_SERCA2a_pump__Vmaxr * var_SERCA2a_pump__rb))) / (1.0 + var_SERCA2a_pump__fb + var_SERCA2a_pump__rb); // millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__tau_tr = 0.0005747; // second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__Ca_JSR; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__Ca_NSR = var_chaste_interface__intracellular_ion_concentrations__Ca_NSR; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__J_tr = (var_intracellular_and_SR_Ca_fluxes__Ca_NSR - var_intracellular_and_SR_Ca_fluxes__Ca_JSR) / var_intracellular_and_SR_Ca_fluxes__tau_tr; // millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__tau_xfer = 0.0267; // second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__J_xfer = (var_intracellular_and_SR_Ca_fluxes__Ca_ss - var_intracellular_and_SR_Ca_fluxes__Ca_i) / var_intracellular_and_SR_Ca_fluxes__tau_xfer; // millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__k_htrpn_plus = 200000.0; // per_millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__k_htrpn_minus = 0.066; // per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__HTRPNCa = var_chaste_interface__intracellular_and_SR_Ca_fluxes__HTRPNCa; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__HTRPN_tot = 0.14; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__J_HTRPNCa = (var_intracellular_and_SR_Ca_fluxes__k_htrpn_plus * var_intracellular_and_SR_Ca_fluxes__Ca_i * (var_intracellular_and_SR_Ca_fluxes__HTRPN_tot - var_intracellular_and_SR_Ca_fluxes__HTRPNCa)) - (var_intracellular_and_SR_Ca_fluxes__k_htrpn_minus * var_intracellular_and_SR_Ca_fluxes__HTRPNCa); // millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__LTRPNCa = var_chaste_interface__intracellular_and_SR_Ca_fluxes__LTRPNCa; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__LTRPN_tot = 0.07; // millimolar
            const NekDouble var_intracellular_and_SR_Ca_fluxes__k_ltrpn_plus = 40000.0; // per_millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__k_ltrpn_minus = 40.0; // per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__J_LTRPNCa = (var_intracellular_and_SR_Ca_fluxes__k_ltrpn_plus * var_intracellular_and_SR_Ca_fluxes__Ca_i * (var_intracellular_and_SR_Ca_fluxes__LTRPN_tot - var_intracellular_and_SR_Ca_fluxes__LTRPNCa)) - (var_intracellular_and_SR_Ca_fluxes__k_ltrpn_minus * var_intracellular_and_SR_Ca_fluxes__LTRPNCa); // millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__J_trpn = var_intracellular_and_SR_Ca_fluxes__J_HTRPNCa + var_intracellular_and_SR_Ca_fluxes__J_LTRPNCa; // millimolar_per_second
            const NekDouble var_intracellular_and_SR_Ca_fluxes__d_HTRPNCa_d_environment__time = (var_intracellular_and_SR_Ca_fluxes__k_htrpn_plus * var_intracellular_and_SR_Ca_fluxes__Ca_i * (var_intracellular_and_SR_Ca_fluxes__HTRPN_tot - var_intracellular_and_SR_Ca_fluxes__HTRPNCa)) - (var_intracellular_and_SR_Ca_fluxes__k_htrpn_minus * var_intracellular_and_SR_Ca_fluxes__HTRPNCa); // 'millimole per litre per second'
            const NekDouble var_intracellular_and_SR_Ca_fluxes__d_LTRPNCa_d_environment__time = (var_intracellular_and_SR_Ca_fluxes__k_ltrpn_plus * var_intracellular_and_SR_Ca_fluxes__Ca_i * (var_intracellular_and_SR_Ca_fluxes__LTRPN_tot - var_intracellular_and_SR_Ca_fluxes__LTRPNCa)) - (var_intracellular_and_SR_Ca_fluxes__k_ltrpn_minus * var_intracellular_and_SR_Ca_fluxes__LTRPNCa); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__Ca_i = var_chaste_interface__intracellular_ion_concentrations__Ca_i; // millimolar
            const NekDouble var_intracellular_ion_concentrations__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__Ca_ss; // millimolar
            const NekDouble var_intracellular_ion_concentrations__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__Ca_JSR; // millimolar
            const NekDouble var_intracellular_ion_concentrations__V_myo = 9.36e-06; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__V_JSR = 5.6e-07; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__V_NSR = 5.04e-06; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__V_SS = 1.2e-08; // micro_litre
            const NekDouble var_intracellular_ion_concentrations__K_mCMDN = 0.00238; // millimolar
            const NekDouble var_intracellular_ion_concentrations__K_mCSQN = 0.8; // millimolar
            const NekDouble var_intracellular_ion_concentrations__K_mEGTA = 0.00015; // millimolar
            const NekDouble var_intracellular_ion_concentrations__CMDN_tot = 0.05; // millimolar
            const NekDouble var_intracellular_ion_concentrations__CSQN_tot = 15.0; // millimolar
            const NekDouble var_intracellular_ion_concentrations__EGTA_tot = 10.0; // millimolar
            const NekDouble var_intracellular_ion_concentrations__beta_i = 1.0 / (1.0 + ((var_intracellular_ion_concentrations__CMDN_tot * var_intracellular_ion_concentrations__K_mCMDN) / pow(var_intracellular_ion_concentrations__K_mCMDN + var_intracellular_ion_concentrations__Ca_i, 2.0)) + ((var_intracellular_ion_concentrations__EGTA_tot * var_intracellular_ion_concentrations__K_mEGTA) / pow(var_intracellular_ion_concentrations__K_mEGTA + var_intracellular_ion_concentrations__Ca_i, 2.0))); // dimensionless
            const NekDouble var_intracellular_ion_concentrations__beta_SS = 1.0 / (1.0 + ((var_intracellular_ion_concentrations__CMDN_tot * var_intracellular_ion_concentrations__K_mCMDN) / pow(var_intracellular_ion_concentrations__K_mCMDN + var_intracellular_ion_concentrations__Ca_ss, 2.0))); // dimensionless
            const NekDouble var_intracellular_ion_concentrations__beta_JSR = 1.0 / (1.0 + ((var_intracellular_ion_concentrations__CSQN_tot * var_intracellular_ion_concentrations__K_mCSQN) / pow(var_intracellular_ion_concentrations__K_mCSQN + var_intracellular_ion_concentrations__Ca_JSR, 2.0))); // dimensionless
            const NekDouble var_intracellular_ion_concentrations__F = var_membrane__F; // coulomb_per_mole
            const NekDouble var_intracellular_ion_concentrations__i_Na = var_sodium_current__i_Na; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_Ca_L = var_L_type_Ca_channel__i_Ca_L; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_B_Na = var_background_currents__i_B_Na; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_NaCa = var_Na_Ca_ion_exchanger_current__i_NaCa; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_NaK = var_sodium_potassium_pump__i_NaK; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_f_Na = var_hyperpolarisation_activated_current__i_f_Na; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_f_K = var_hyperpolarisation_activated_current__i_f_K; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_B_K = var_background_currents__i_B_K; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_K1 = var_inward_rectifier__i_K1; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_t = var_Ca_independent_transient_outward_K_current__i_t; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_ss = var_steady_state_outward_K_current__i_ss; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_Ca_P = var_sarcolemmal_calcium_pump_current__i_Ca_P; // nanoA
            const NekDouble var_intracellular_ion_concentrations__i_B_Ca = var_background_currents__i_B_Ca; // nanoA
            const NekDouble var_intracellular_ion_concentrations__J_up = var_SERCA2a_pump__J_up; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_rel = var_SR_Ca_release_channel__J_rel; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_xfer = var_intracellular_and_SR_Ca_fluxes__J_xfer; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_trpn = var_intracellular_and_SR_Ca_fluxes__J_trpn; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__J_tr = var_intracellular_and_SR_Ca_fluxes__J_tr; // millimolar_per_second
            const NekDouble var_intracellular_ion_concentrations__d_Na_i_d_environment__time = (-(var_intracellular_ion_concentrations__i_Na + var_intracellular_ion_concentrations__i_B_Na + (var_intracellular_ion_concentrations__i_NaCa * 3.0) + (var_intracellular_ion_concentrations__i_NaK * 3.0) + var_intracellular_ion_concentrations__i_f_Na)) / (var_intracellular_ion_concentrations__V_myo * var_intracellular_ion_concentrations__F); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_i_d_environment__time = var_intracellular_ion_concentrations__beta_i * (var_intracellular_ion_concentrations__J_xfer - (var_intracellular_ion_concentrations__J_up + var_intracellular_ion_concentrations__J_trpn + (((var_intracellular_ion_concentrations__i_B_Ca - (2.0 * var_intracellular_ion_concentrations__i_NaCa)) + var_intracellular_ion_concentrations__i_Ca_P) / (2.0 * var_intracellular_ion_concentrations__V_myo * var_intracellular_ion_concentrations__F)))); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_K_i_d_environment__time = (-(var_intracellular_ion_concentrations__i_ss + var_intracellular_ion_concentrations__i_B_K + var_intracellular_ion_concentrations__i_t + var_intracellular_ion_concentrations__i_K1 + var_intracellular_ion_concentrations__i_f_K + (var_intracellular_ion_concentrations__i_NaK * (-2.0)))) / (var_intracellular_ion_concentrations__V_myo * var_intracellular_ion_concentrations__F); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_ss_d_environment__time = var_intracellular_ion_concentrations__beta_SS * ((((var_intracellular_ion_concentrations__J_rel * var_intracellular_ion_concentrations__V_JSR) / var_intracellular_ion_concentrations__V_SS) - ((var_intracellular_ion_concentrations__J_xfer * var_intracellular_ion_concentrations__V_myo) / var_intracellular_ion_concentrations__V_SS)) - (var_intracellular_ion_concentrations__i_Ca_L / (2.0 * var_intracellular_ion_concentrations__V_SS * var_intracellular_ion_concentrations__F))); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_JSR_d_environment__time = var_intracellular_ion_concentrations__beta_JSR * (var_intracellular_ion_concentrations__J_tr - var_intracellular_ion_concentrations__J_rel); // 'millimole per litre per second'
            const NekDouble var_intracellular_ion_concentrations__d_Ca_NSR_d_environment__time = ((var_intracellular_ion_concentrations__J_up * var_intracellular_ion_concentrations__V_myo) / var_intracellular_ion_concentrations__V_NSR) - ((var_intracellular_ion_concentrations__J_tr * var_intracellular_ion_concentrations__V_JSR) / var_intracellular_ion_concentrations__V_NSR); // 'millimole per litre per second'
            const NekDouble var_chaste_interface__sodium_current_m_gate__d_m_d_environment__time_converter = var_sodium_current__sodium_current_m_gate__d_m_d_environment__time; // per_second
            const NekDouble var_chaste_interface__sodium_current_m_gate__d_m_d_environment__time = 0.001 * var_chaste_interface__sodium_current_m_gate__d_m_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__sodium_current_h_gate__d_h_d_environment__time_converter = var_sodium_current__sodium_current_h_gate__d_h_d_environment__time; // per_second
            const NekDouble var_chaste_interface__sodium_current_h_gate__d_h_d_environment__time = 0.001 * var_chaste_interface__sodium_current_h_gate__d_h_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__sodium_current_j_gate__d_j_d_environment__time_converter = var_sodium_current__sodium_current_j_gate__d_j_d_environment__time; // per_second
            const NekDouble var_chaste_interface__sodium_current_j_gate__d_j_d_environment__time = 0.001 * var_chaste_interface__sodium_current_j_gate__d_j_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_channel_d_gate__d_d_d_environment__time_converter = var_L_type_Ca_channel__L_type_Ca_channel_d_gate__d_d_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_channel_d_gate__d_d_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_channel_d_gate__d_d_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time_converter = var_L_type_Ca_channel__L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time_converter = var_L_type_Ca_channel__L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time_converter = var_L_type_Ca_channel__L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time; // per_second
            const NekDouble var_chaste_interface__L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time = 0.001 * var_chaste_interface__L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time_converter = var_Ca_independent_transient_outward_K_current__Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time; // per_second
            const NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time = 0.001 * var_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time_converter = var_Ca_independent_transient_outward_K_current__Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time; // per_second
            const NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time = 0.001 * var_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time_converter = var_Ca_independent_transient_outward_K_current__Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time; // per_second
            const NekDouble var_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time = 0.001 * var_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time_converter = var_steady_state_outward_K_current__steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time; // per_second
            const NekDouble var_chaste_interface__steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time = 0.001 * var_chaste_interface__steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time_converter = var_steady_state_outward_K_current__steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time; // per_second
            const NekDouble var_chaste_interface__steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time = 0.001 * var_chaste_interface__steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__hyperpolarisation_activated_current_y_gate__d_y_d_environment__time_converter = var_hyperpolarisation_activated_current__hyperpolarisation_activated_current_y_gate__d_y_d_environment__time; // per_second
            const NekDouble var_chaste_interface__hyperpolarisation_activated_current_y_gate__d_y_d_environment__time = 0.001 * var_chaste_interface__hyperpolarisation_activated_current_y_gate__d_y_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_O1_d_environment__time_converter = var_SR_Ca_release_channel__d_P_O1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_O1_d_environment__time = 0.001 * var_chaste_interface__SR_Ca_release_channel__d_P_O1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_O2_d_environment__time_converter = var_SR_Ca_release_channel__d_P_O2_d_environment__time; // per_second
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_O2_d_environment__time = 0.001 * var_chaste_interface__SR_Ca_release_channel__d_P_O2_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_C1_d_environment__time_converter = var_SR_Ca_release_channel__d_P_C1_d_environment__time; // per_second
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_C1_d_environment__time = 0.001 * var_chaste_interface__SR_Ca_release_channel__d_P_C1_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_C2_d_environment__time_converter = var_SR_Ca_release_channel__d_P_C2_d_environment__time; // per_second
            const NekDouble var_chaste_interface__SR_Ca_release_channel__d_P_C2_d_environment__time = 0.001 * var_chaste_interface__SR_Ca_release_channel__d_P_C2_d_environment__time_converter; // 'per millisecond'
            const NekDouble var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_HTRPNCa_d_environment__time_converter = var_intracellular_and_SR_Ca_fluxes__d_HTRPNCa_d_environment__time; // ___units_26
            const NekDouble var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_HTRPNCa_d_environment__time = 0.001 * var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_HTRPNCa_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_LTRPNCa_d_environment__time_converter = var_intracellular_and_SR_Ca_fluxes__d_LTRPNCa_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_LTRPNCa_d_environment__time = 0.001 * var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_LTRPNCa_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Na_i_d_environment__time_converter = var_intracellular_ion_concentrations__d_Na_i_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Na_i_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Na_i_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_i_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ca_i_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_i_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ca_i_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_K_i_d_environment__time_converter = var_intracellular_ion_concentrations__d_K_i_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_K_i_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_K_i_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ca_ss_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ca_JSR_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time_converter = var_intracellular_ion_concentrations__d_Ca_NSR_d_environment__time; // millimole_per_litre_per_second
            const NekDouble var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time = 0.001 * var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time_converter; // 'millimolar per millisecond'
            const NekDouble d_dt_chaste_interface__sodium_current_m_gate__m = var_chaste_interface__sodium_current_m_gate__d_m_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__sodium_current_h_gate__h = var_chaste_interface__sodium_current_h_gate__d_h_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__sodium_current_j_gate__j = var_chaste_interface__sodium_current_j_gate__d_j_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_channel_d_gate__d = var_chaste_interface__L_type_Ca_channel_d_gate__d_d_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_channel_f_11_gate__f_11 = var_chaste_interface__L_type_Ca_channel_f_11_gate__d_f_11_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_channel_f_12_gate__f_12 = var_chaste_interface__L_type_Ca_channel_f_12_gate__d_f_12_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__L_type_Ca_channel_Ca_inact_gate__Ca_inact = var_chaste_interface__L_type_Ca_channel_Ca_inact_gate__d_Ca_inact_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__r = var_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__d_r_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__s = var_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__d_s_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__s_slow = var_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__d_s_slow_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__steady_state_outward_K_current_r_ss_gate__r_ss = var_chaste_interface__steady_state_outward_K_current_r_ss_gate__d_r_ss_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__steady_state_outward_K_current_s_ss_gate__s_ss = var_chaste_interface__steady_state_outward_K_current_s_ss_gate__d_s_ss_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__hyperpolarisation_activated_current_y_gate__y = var_chaste_interface__hyperpolarisation_activated_current_y_gate__d_y_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__SR_Ca_release_channel__P_O1 = var_chaste_interface__SR_Ca_release_channel__d_P_O1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__SR_Ca_release_channel__P_O2 = var_chaste_interface__SR_Ca_release_channel__d_P_O2_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__SR_Ca_release_channel__P_C1 = var_chaste_interface__SR_Ca_release_channel__d_P_C1_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__SR_Ca_release_channel__P_C2 = var_chaste_interface__SR_Ca_release_channel__d_P_C2_d_environment__time; // 'per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_and_SR_Ca_fluxes__HTRPNCa = var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_HTRPNCa_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_and_SR_Ca_fluxes__LTRPNCa = var_chaste_interface__intracellular_and_SR_Ca_fluxes__d_LTRPNCa_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Na_i = var_chaste_interface__intracellular_ion_concentrations__d_Na_i_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_i = var_chaste_interface__intracellular_ion_concentrations__d_Ca_i_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__K_i = var_chaste_interface__intracellular_ion_concentrations__d_K_i_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_ss = var_chaste_interface__intracellular_ion_concentrations__d_Ca_ss_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_JSR = var_chaste_interface__intracellular_ion_concentrations__d_Ca_JSR_d_environment__time; // 'millimole per litre per millisecond'
            const NekDouble d_dt_chaste_interface__intracellular_ion_concentrations__Ca_NSR = var_chaste_interface__intracellular_ion_concentrations__d_Ca_NSR_d_environment__time; // 'millimole per litre per millisecond'

            const NekDouble var_membrane__Cm = 0.0001; // microF
            const NekDouble var_membrane__i_Na = var_sodium_current__i_Na; // nanoA
            const NekDouble var_membrane__i_Ca_L = var_L_type_Ca_channel__i_Ca_L; // nanoA
            const NekDouble var_membrane__i_t = var_Ca_independent_transient_outward_K_current__i_t; // nanoA
            const NekDouble var_membrane__i_ss = var_steady_state_outward_K_current__i_ss; // nanoA
            const NekDouble var_hyperpolarisation_activated_current__i_f = var_hyperpolarisation_activated_current__i_f_Na + var_hyperpolarisation_activated_current__i_f_K; // nanoA
            const NekDouble var_membrane__i_f = var_hyperpolarisation_activated_current__i_f; // nanoA
            const NekDouble var_membrane__i_K1 = var_inward_rectifier__i_K1; // nanoA
            const NekDouble var_background_currents__i_B = var_background_currents__i_B_Na + var_background_currents__i_B_Ca + var_background_currents__i_B_K; // nanoA
            const NekDouble var_membrane__i_B = var_background_currents__i_B; // nanoA
            const NekDouble var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK; // nanoA
            const NekDouble var_membrane__i_NaCa = var_Na_Ca_ion_exchanger_current__i_NaCa; // nanoA
            const NekDouble var_membrane__i_Ca_P = var_sarcolemmal_calcium_pump_current__i_Ca_P; // nanoA
            const NekDouble var_chaste_interface__membrane__i_Stim = 0.0;
            const NekDouble var_membrane__i_Stim_converter = var_chaste_interface__membrane__i_Stim; // uA_per_cm2
            const NekDouble var_membrane__chaste_interface__chaste_membrane_capacitance = 1.0; // uF_per_cm2
            const NekDouble var_membrane__i_Stim = 1000.0 * ((var_membrane__i_Stim_converter * var_membrane__Cm) / var_membrane__chaste_interface__chaste_membrane_capacitance); // nanoA
            const NekDouble var_membrane__d_V_d_environment__time = (-(var_membrane__i_Na + var_membrane__i_Ca_L + var_membrane__i_t + var_membrane__i_ss + var_membrane__i_f + var_membrane__i_K1 + var_membrane__i_B + var_membrane__i_NaK + var_membrane__i_NaCa + var_membrane__i_Ca_P + var_membrane__i_Stim)) / var_membrane__Cm; // 'millivolt per second'
            const NekDouble var_chaste_interface__membrane__d_V_d_environment__time_converter = var_membrane__d_V_d_environment__time; // ___units_1
            const NekDouble var_chaste_interface__membrane__d_V_d_environment__time = 0.001 * var_chaste_interface__membrane__d_V_d_environment__time_converter; // 'millivolt per millisecond'
            d_dt_chaste_interface__membrane__V = var_chaste_interface__membrane__d_V_d_environment__time; // 'millivolt per millisecond'
            outarray[0][i] = d_dt_chaste_interface__membrane__V;
            outarray[1][i] = 0.0;
            outarray[2][i] = d_dt_chaste_interface__sodium_current_m_gate__m;
            outarray[3][i] = d_dt_chaste_interface__sodium_current_h_gate__h;
            outarray[4][i] = d_dt_chaste_interface__sodium_current_j_gate__j;
            outarray[5][i] = d_dt_chaste_interface__L_type_Ca_channel_d_gate__d;
            outarray[6][i] = d_dt_chaste_interface__L_type_Ca_channel_f_11_gate__f_11;
            outarray[7][i] = d_dt_chaste_interface__L_type_Ca_channel_f_12_gate__f_12;
            outarray[8][i] = d_dt_chaste_interface__L_type_Ca_channel_Ca_inact_gate__Ca_inact;
            outarray[9][i] = d_dt_chaste_interface__Ca_independent_transient_outward_K_current_r_gate__r;
            outarray[10][i] = d_dt_chaste_interface__Ca_independent_transient_outward_K_current_s_gate__s;
            outarray[11][i] = d_dt_chaste_interface__Ca_independent_transient_outward_K_current_s_slow_gate__s_slow;
            outarray[12][i] = d_dt_chaste_interface__steady_state_outward_K_current_r_ss_gate__r_ss;
            outarray[13][i] = d_dt_chaste_interface__steady_state_outward_K_current_s_ss_gate__s_ss;
            outarray[14][i] = d_dt_chaste_interface__hyperpolarisation_activated_current_y_gate__y;
            outarray[15][i] = d_dt_chaste_interface__SR_Ca_release_channel__P_O1;
            outarray[16][i] = d_dt_chaste_interface__SR_Ca_release_channel__P_O2;
            outarray[17][i] = d_dt_chaste_interface__SR_Ca_release_channel__P_C1;
            outarray[18][i] = d_dt_chaste_interface__SR_Ca_release_channel__P_C2;
            outarray[19][i] = d_dt_chaste_interface__intracellular_and_SR_Ca_fluxes__HTRPNCa;
            outarray[20][i] = d_dt_chaste_interface__intracellular_and_SR_Ca_fluxes__LTRPNCa;
            outarray[21][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Na_i;
            outarray[22][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_i;
            outarray[23][i] = d_dt_chaste_interface__intracellular_ion_concentrations__K_i;
            outarray[24][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_ss;
            outarray[25][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_JSR;
            outarray[26][i] = d_dt_chaste_interface__intracellular_ion_concentrations__Ca_NSR;
        }

    }

    /**
    *
    */
    void PanditGilesDemir03::v_GenerateSummary(SummaryList& s)
    {
        SolverUtils::AddSummaryItem(s, "Cell model", "PanditGilesDemir03");
    }

}

