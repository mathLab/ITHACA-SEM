///////////////////////////////////////////////////////////////////////////////
//
// File TenTusscher06.cpp
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
// Description: ten Tusscher 2006 Epicardium cell model
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
//#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <CardiacEPSolver/CellModels/TenTusscher06.h>

namespace Nektar
{

std::string TenTusscher06::className
= GetCellModelFactory().RegisterCreatorFunction(
        "TenTusscher06",
        TenTusscher06::create,
        "ten Tusscher 2006.");

std::string TenTusscher06::lookupIds[4] = {
        LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                "Epicardium", TenTusscher06::eEpicardium),
        LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                "Endocardium", TenTusscher06::eEndocardium),
        LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                "Mid", TenTusscher06::eMid),
        LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                "Ischemia", TenTusscher06::eIschemia)
};

std::string TenTusscher06::def =
        LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "CellModelVariant", "Epicardium");

/**
 *
 */
TenTusscher06::TenTusscher06(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const MultiRegions::ExpListSharedPtr& pField):
                    CellModel(pSession, pField)
{
    model_variant = pSession->GetSolverInfoAsEnum<
            TenTusscher06::Variants>("CellModelVariant");

    switch (model_variant) {
    case eEpicardium:
        g_to = 0.294;
        g_Ks = 0.392;
        s_inf_factor = 20.0;
        s_tau_f1 = 85.0;
        s_tau_f2 = 45.0;
        s_tau_f3 = 320.0;
        s_tau_f4 = 0.0;
        s_tau_f5 = 1.0;
        k_0 = 5.4;
        break;
    case eEndocardium:
        g_to = 0.073;
        g_Ks = 0.392;
        s_inf_factor = 28.0;
        s_tau_f1 = 1000.0;
        s_tau_f2 = 67.0;
        s_tau_f3 = 1000.0;
        s_tau_f4 = 8.0;
        s_tau_f5 = 0.0;
        k_0 = 5.4;
        break;
    case eMid:
        g_to = 0.294;
        g_Ks = 0.098;
        s_inf_factor = 20.0;
        s_tau_f1 = 85.0;
        s_tau_f2 = 45.0;
        s_tau_f3 = 320.0;
        s_tau_f4 = 0.0;
        s_tau_f5 = 1.0;
        k_0 = 5.4;
        break;
    case eIschemia:
        g_to = 0.294;
        g_Ks = 0.392;
        s_inf_factor = 20.0;
        s_tau_f1 = 85.0;
        s_tau_f2 = 45.0;
        s_tau_f3 = 320.0;
        s_tau_f4 = 0.0;
        s_tau_f5 = 1.0;
        k_0 = 9.0;
        break;
    }
    m_nvar = 19;

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
    m_concentrations.push_back(13);
    m_concentrations.push_back(14);
    m_concentrations.push_back(15);
    m_concentrations.push_back(16);
    m_concentrations.push_back(17);
    m_concentrations.push_back(18);
}


/**
 *
 */
void TenTusscher06::v_Update(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
        Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
{
    for (unsigned int i = 0; i < m_nq; ++i)
    {
        // Inputs:
        // Time units: millisecond
        NekDouble var_chaste_interface__membrane__V = inarray[0][i];
        // Units: millivolt; Initial value: -85.23
        NekDouble var_chaste_interface__rapid_time_dependent_potassium_current_Xr1_gate__Xr1 = inarray[1][i];
        // Units: dimensionless; Initial value: 0.00621
        NekDouble var_chaste_interface__rapid_time_dependent_potassium_current_Xr2_gate__Xr2 = inarray[2][i];
        // Units: dimensionless; Initial value: 0.4712
        NekDouble var_chaste_interface__slow_time_dependent_potassium_current_Xs_gate__Xs = inarray[3][i];
        // Units: dimensionless; Initial value: 0.0095
        NekDouble var_chaste_interface__fast_sodium_current_m_gate__m = inarray[4][i];
        // Units: dimensionless; Initial value: 0.00172
        NekDouble var_chaste_interface__fast_sodium_current_h_gate__h = inarray[5][i];
        // Units: dimensionless; Initial value: 0.7444
        NekDouble var_chaste_interface__fast_sodium_current_j_gate__j = inarray[6][i];
        // Units: dimensionless; Initial value: 0.7045
        NekDouble var_chaste_interface__L_type_Ca_current_d_gate__d = inarray[7][i];
        // Units: dimensionless; Initial value: 3.373e-5
        NekDouble var_chaste_interface__L_type_Ca_current_f_gate__f = inarray[8][i];
        // Units: dimensionless; Initial value: 0.7888
        NekDouble var_chaste_interface__L_type_Ca_current_f2_gate__f2 = inarray[9][i];
        // Units: dimensionless; Initial value: 0.9755
        NekDouble var_chaste_interface__L_type_Ca_current_fCass_gate__fCass = inarray[10][i];
        // Units: dimensionless; Initial value: 0.9953
        NekDouble var_chaste_interface__transient_outward_current_s_gate__s = inarray[11][i];
        // Units: dimensionless; Initial value: 0.999998
        NekDouble var_chaste_interface__transient_outward_current_r_gate__r = inarray[12][i];
        // Units: dimensionless; Initial value: 2.42e-8
        NekDouble var_chaste_interface__calcium_dynamics__Ca_i = inarray[13][i];
        // Units: millimolar; Initial value: 0.000126
        NekDouble var_chaste_interface__calcium_dynamics__Ca_SR = inarray[14][i];
        // Units: millimolar; Initial value: 3.64
        NekDouble var_chaste_interface__calcium_dynamics__Ca_ss = inarray[15][i];
        // Units: millimolar; Initial value: 0.00036
        NekDouble var_chaste_interface__calcium_dynamics__R_prime = inarray[16][i];
        // Units: dimensionless; Initial value: 0.9073
        NekDouble var_chaste_interface__sodium_dynamics__Na_i = inarray[17][i];
        // Units: millimolar; Initial value: 8.604
        NekDouble var_chaste_interface__potassium_dynamics__K_i = inarray[18][i];
        // Units: millimolar; Initial value: 136.89


        // Mathematics
        NekDouble d_dt_chaste_interface__membrane__V;
        const NekDouble var_membrane__R = 8314.472; // joule_per_mole_kelvin
        const NekDouble var_membrane__T = 310.0; // kelvin
        const NekDouble var_membrane__F = 96485.3415; // coulomb_per_millimole
        const NekDouble var_membrane__Cm = 0.185; // microF
        const NekDouble var_membrane__V_c = 0.016404; // micrometre3
        const NekDouble var_inward_rectifier_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_reversal_potentials__K_i = var_chaste_interface__potassium_dynamics__K_i; // millimolar
        const NekDouble var_reversal_potentials__R = var_membrane__R; // joule_per_mole_kelvin
        const NekDouble var_reversal_potentials__T = var_membrane__T; // kelvin
        const NekDouble var_reversal_potentials__F = var_membrane__F; // coulomb_per_millimole
        const NekDouble var_potassium_dynamics__K_o = k_0; // millimolar
        const NekDouble var_reversal_potentials__K_o = var_potassium_dynamics__K_o; // millimolar
        const NekDouble var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i); // millivolt
        const NekDouble var_inward_rectifier_potassium_current__E_K = var_reversal_potentials__E_K; // millivolt
        const NekDouble var_inward_rectifier_potassium_current__beta_K1 = ((3.0 * exp(0.0002 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) + 100.0))) + exp(0.1 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) - 10.0))) / (1.0 + exp((-0.5) * (var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K))); // dimensionless
        const NekDouble var_inward_rectifier_potassium_current__alpha_K1 = 0.1 / (1.0 + exp(0.06 * ((var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K) - 200.0))); // dimensionless
        const NekDouble var_inward_rectifier_potassium_current__xK1_inf = var_inward_rectifier_potassium_current__alpha_K1 / (var_inward_rectifier_potassium_current__alpha_K1 + var_inward_rectifier_potassium_current__beta_K1); // dimensionless
        const NekDouble var_inward_rectifier_potassium_current__K_o = var_potassium_dynamics__K_o; // millimolar
        const NekDouble var_inward_rectifier_potassium_current__g_K1 = 5.405; // nanoS_per_picoF
        const NekDouble var_inward_rectifier_potassium_current__i_K1 = var_inward_rectifier_potassium_current__g_K1 * var_inward_rectifier_potassium_current__xK1_inf * sqrt(var_inward_rectifier_potassium_current__K_o / 5.4) * (var_inward_rectifier_potassium_current__V - var_inward_rectifier_potassium_current__E_K); // picoA_per_picoF
        const NekDouble var_transient_outward_current__s = var_chaste_interface__transient_outward_current_s_gate__s; // dimensionless
        const NekDouble var_transient_outward_current__r = var_chaste_interface__transient_outward_current_r_gate__r; // dimensionless
        const NekDouble var_transient_outward_current__g_to = g_to; // nanoS_per_picoF
        const NekDouble var_transient_outward_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_transient_outward_current__E_K = var_reversal_potentials__E_K; // millivolt
        const NekDouble var_transient_outward_current__i_to = var_transient_outward_current__g_to * var_transient_outward_current__r * var_transient_outward_current__s * (var_transient_outward_current__V - var_transient_outward_current__E_K); // picoA_per_picoF
        const NekDouble var_rapid_time_dependent_potassium_current__Xr1 = var_chaste_interface__rapid_time_dependent_potassium_current_Xr1_gate__Xr1; // dimensionless
        const NekDouble var_rapid_time_dependent_potassium_current__Xr2 = var_chaste_interface__rapid_time_dependent_potassium_current_Xr2_gate__Xr2; // dimensionless
        const NekDouble var_rapid_time_dependent_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_rapid_time_dependent_potassium_current__K_o = var_potassium_dynamics__K_o; // millimolar
        const NekDouble var_rapid_time_dependent_potassium_current__E_K = var_reversal_potentials__E_K; // millivolt
        const NekDouble var_rapid_time_dependent_potassium_current__g_Kr = 0.153; // nanoS_per_picoF
        const NekDouble var_rapid_time_dependent_potassium_current__i_Kr = var_rapid_time_dependent_potassium_current__g_Kr * sqrt(var_rapid_time_dependent_potassium_current__K_o / 5.4) * var_rapid_time_dependent_potassium_current__Xr1 * var_rapid_time_dependent_potassium_current__Xr2 * (var_rapid_time_dependent_potassium_current__V - var_rapid_time_dependent_potassium_current__E_K); // picoA_per_picoF
        const NekDouble var_slow_time_dependent_potassium_current__g_Ks = g_Ks; // nanoS_per_picoF
        const NekDouble var_sodium_dynamics__Na_o = 140.0; // millimolar
        const NekDouble var_reversal_potentials__Na_o = var_sodium_dynamics__Na_o; // millimolar
        const NekDouble var_reversal_potentials__Na_i = var_chaste_interface__sodium_dynamics__Na_i; // millimolar
        const NekDouble var_reversal_potentials__P_kna = 0.03; // dimensionless
        const NekDouble var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i))); // millivolt
        const NekDouble var_slow_time_dependent_potassium_current__E_Ks = var_reversal_potentials__E_Ks; // millivolt
        const NekDouble var_slow_time_dependent_potassium_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_slow_time_dependent_potassium_current__Xs = var_chaste_interface__slow_time_dependent_potassium_current_Xs_gate__Xs; // dimensionless
        const NekDouble var_slow_time_dependent_potassium_current__i_Ks = var_slow_time_dependent_potassium_current__g_Ks * pow(var_slow_time_dependent_potassium_current__Xs, 2.0) * (var_slow_time_dependent_potassium_current__V - var_slow_time_dependent_potassium_current__E_Ks); // picoA_per_picoF
        const NekDouble var_L_type_Ca_current__Ca_ss = var_chaste_interface__calcium_dynamics__Ca_ss; // millimolar
        const NekDouble var_L_type_Ca_current__g_CaL = 3.98e-05; // litre_per_farad_second
        const NekDouble var_L_type_Ca_current__f = var_chaste_interface__L_type_Ca_current_f_gate__f; // dimensionless
        const NekDouble var_L_type_Ca_current__d = var_chaste_interface__L_type_Ca_current_d_gate__d; // dimensionless
        const NekDouble var_L_type_Ca_current__F = var_membrane__F; // coulomb_per_millimole
        const NekDouble var_L_type_Ca_current__f2 = var_chaste_interface__L_type_Ca_current_f2_gate__f2; // dimensionless
        const NekDouble var_L_type_Ca_current__fCass = var_chaste_interface__L_type_Ca_current_fCass_gate__fCass; // dimensionless
        const NekDouble var_L_type_Ca_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_L_type_Ca_current__T = var_membrane__T; // kelvin
        const NekDouble var_calcium_dynamics__Ca_o = 2.0; // millimolar
        const NekDouble var_L_type_Ca_current__Ca_o = var_calcium_dynamics__Ca_o; // millimolar
        const NekDouble var_L_type_Ca_current__R = var_membrane__R; // joule_per_mole_kelvin
        const NekDouble var_L_type_Ca_current__i_CaL = (((var_L_type_Ca_current__g_CaL * var_L_type_Ca_current__d * var_L_type_Ca_current__f * var_L_type_Ca_current__f2 * var_L_type_Ca_current__fCass * 4.0 * (var_L_type_Ca_current__V - 15.0) * pow(var_L_type_Ca_current__F, 2.0)) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) * ((0.25 * var_L_type_Ca_current__Ca_ss * exp((2.0 * (var_L_type_Ca_current__V - 15.0) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T))) - var_L_type_Ca_current__Ca_o)) / (exp((2.0 * (var_L_type_Ca_current__V - 15.0) * var_L_type_Ca_current__F) / (var_L_type_Ca_current__R * var_L_type_Ca_current__T)) - 1.0); // picoA_per_picoF
        const NekDouble var_sodium_potassium_pump_current__Na_i = var_chaste_interface__sodium_dynamics__Na_i; // millimolar
        const NekDouble var_sodium_potassium_pump_current__R = var_membrane__R; // joule_per_mole_kelvin
        const NekDouble var_sodium_potassium_pump_current__T = var_membrane__T; // kelvin
        const NekDouble var_sodium_potassium_pump_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_sodium_potassium_pump_current__K_mk = 1.0; // millimolar
        const NekDouble var_sodium_potassium_pump_current__P_NaK = 2.724; // picoA_per_picoF
        const NekDouble var_sodium_potassium_pump_current__K_mNa = 40.0; // millimolar
        const NekDouble var_sodium_potassium_pump_current__F = var_membrane__F; // coulomb_per_millimole
        const NekDouble var_sodium_potassium_pump_current__K_o = var_potassium_dynamics__K_o; // millimolar
        const NekDouble var_sodium_potassium_pump_current__i_NaK = ((((var_sodium_potassium_pump_current__P_NaK * var_sodium_potassium_pump_current__K_o) / (var_sodium_potassium_pump_current__K_o + var_sodium_potassium_pump_current__K_mk)) * var_sodium_potassium_pump_current__Na_i) / (var_sodium_potassium_pump_current__Na_i + var_sodium_potassium_pump_current__K_mNa)) / (1.0 + (0.1245 * exp(((-0.1) * var_sodium_potassium_pump_current__V * var_sodium_potassium_pump_current__F) / (var_sodium_potassium_pump_current__R * var_sodium_potassium_pump_current__T))) + (0.0353 * exp(((-var_sodium_potassium_pump_current__V) * var_sodium_potassium_pump_current__F) / (var_sodium_potassium_pump_current__R * var_sodium_potassium_pump_current__T)))); // picoA_per_picoF
        const NekDouble var_fast_sodium_current__j = var_chaste_interface__fast_sodium_current_j_gate__j; // dimensionless
        const NekDouble var_fast_sodium_current__h = var_chaste_interface__fast_sodium_current_h_gate__h; // dimensionless
        const NekDouble var_fast_sodium_current__g_Na = 14.838; // nanoS_per_picoF
        const NekDouble var_fast_sodium_current__m = var_chaste_interface__fast_sodium_current_m_gate__m; // dimensionless
        const NekDouble var_fast_sodium_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i); // millivolt
        const NekDouble var_fast_sodium_current__E_Na = var_reversal_potentials__E_Na; // millivolt
        const NekDouble var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * var_fast_sodium_current__j * (var_fast_sodium_current__V - var_fast_sodium_current__E_Na); // picoA_per_picoF
        const NekDouble var_sodium_background_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_sodium_background_current__E_Na = var_reversal_potentials__E_Na; // millivolt
        const NekDouble var_sodium_background_current__g_bna = 0.00029; // nanoS_per_picoF
        const NekDouble var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na); // picoA_per_picoF
        const NekDouble var_sodium_calcium_exchanger_current__alpha = 2.5; // dimensionless
        const NekDouble var_sodium_calcium_exchanger_current__gamma = 0.35; // dimensionless
        const NekDouble var_sodium_calcium_exchanger_current__K_sat = 0.1; // dimensionless
        const NekDouble var_sodium_calcium_exchanger_current__Km_Ca = 1.38; // millimolar
        const NekDouble var_sodium_calcium_exchanger_current__K_NaCa = 1000.0; // picoA_per_picoF
        const NekDouble var_sodium_calcium_exchanger_current__F = var_membrane__F; // coulomb_per_millimole
        const NekDouble var_sodium_calcium_exchanger_current__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // millimolar
        const NekDouble var_sodium_calcium_exchanger_current__Ca_o = var_calcium_dynamics__Ca_o; // millimolar
        const NekDouble var_sodium_calcium_exchanger_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_sodium_calcium_exchanger_current__R = var_membrane__R; // joule_per_mole_kelvin
        const NekDouble var_sodium_calcium_exchanger_current__Km_Nai = 87.5; // millimolar
        const NekDouble var_sodium_calcium_exchanger_current__Na_o = var_sodium_dynamics__Na_o; // millimolar
        const NekDouble var_sodium_calcium_exchanger_current__Na_i = var_chaste_interface__sodium_dynamics__Na_i; // millimolar
        const NekDouble var_sodium_calcium_exchanger_current__T = var_membrane__T; // kelvin
        const NekDouble var_sodium_calcium_exchanger_current__i_NaCa = (var_sodium_calcium_exchanger_current__K_NaCa * ((exp((var_sodium_calcium_exchanger_current__gamma * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)) * pow(var_sodium_calcium_exchanger_current__Na_i, 3.0) * var_sodium_calcium_exchanger_current__Ca_o) - (exp(((var_sodium_calcium_exchanger_current__gamma - 1.0) * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T)) * pow(var_sodium_calcium_exchanger_current__Na_o, 3.0) * var_sodium_calcium_exchanger_current__Ca_i * var_sodium_calcium_exchanger_current__alpha))) / ((pow(var_sodium_calcium_exchanger_current__Km_Nai, 3.0) + pow(var_sodium_calcium_exchanger_current__Na_o, 3.0)) * (var_sodium_calcium_exchanger_current__Km_Ca + var_sodium_calcium_exchanger_current__Ca_o) * (1.0 + (var_sodium_calcium_exchanger_current__K_sat * exp(((var_sodium_calcium_exchanger_current__gamma - 1.0) * var_sodium_calcium_exchanger_current__V * var_sodium_calcium_exchanger_current__F) / (var_sodium_calcium_exchanger_current__R * var_sodium_calcium_exchanger_current__T))))); // picoA_per_picoF
        const NekDouble var_reversal_potentials__Ca_o = var_calcium_dynamics__Ca_o; // millimolar
        const NekDouble var_reversal_potentials__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // millimolar
        const NekDouble var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i); // millivolt
        const NekDouble var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca; // millivolt
        const NekDouble var_calcium_background_current__g_bca = 0.000592; // nanoS_per_picoF
        const NekDouble var_calcium_background_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca); // picoA_per_picoF
        const NekDouble var_potassium_pump_current__g_pK = 0.0146; // nanoS_per_picoF
        const NekDouble var_potassium_pump_current__V = var_chaste_interface__membrane__V; // millivolt
        const NekDouble var_potassium_pump_current__E_K = var_reversal_potentials__E_K; // millivolt
        const NekDouble var_potassium_pump_current__i_p_K = (var_potassium_pump_current__g_pK * (var_potassium_pump_current__V - var_potassium_pump_current__E_K)) / (1.0 + exp((25.0 - var_potassium_pump_current__V) / 5.98)); // picoA_per_picoF
        const NekDouble var_calcium_pump_current__K_pCa = 0.0005; // millimolar
        const NekDouble var_calcium_pump_current__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // millimolar
        const NekDouble var_calcium_pump_current__g_pCa = 0.1238; // picoA_per_picoF
        const NekDouble var_calcium_pump_current__i_p_Ca = (var_calcium_pump_current__g_pCa * var_calcium_pump_current__Ca_i) / (var_calcium_pump_current__Ca_i + var_calcium_pump_current__K_pCa); // picoA_per_picoF
        const NekDouble var_chaste_interface__membrane__i_Stim = 0.0;
        const NekDouble var_rapid_time_dependent_potassium_current_Xr1_gate__V = var_rapid_time_dependent_potassium_current__V; // millivolt
        const NekDouble var_rapid_time_dependent_potassium_current_Xr1_gate__alpha_xr1 = 450.0 / (1.0 + exp(((-45.0) - var_rapid_time_dependent_potassium_current_Xr1_gate__V) / 10.0)); // dimensionless
        const NekDouble var_rapid_time_dependent_potassium_current_Xr1_gate__beta_xr1 = 6.0 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr1_gate__V + 30.0) / 11.5)); // dimensionless
        const NekDouble var_rapid_time_dependent_potassium_current_Xr1_gate__tau_xr1 = 1.0 * var_rapid_time_dependent_potassium_current_Xr1_gate__alpha_xr1 * var_rapid_time_dependent_potassium_current_Xr1_gate__beta_xr1; // millisecond
        const NekDouble var_rapid_time_dependent_potassium_current_Xr1_gate__xr1_inf = 1.0 / (1.0 + exp(((-26.0) - var_rapid_time_dependent_potassium_current_Xr1_gate__V) / 7.0)); // dimensionless
        const NekDouble var_rapid_time_dependent_potassium_current_Xr2_gate__V = var_rapid_time_dependent_potassium_current__V; // millivolt
        const NekDouble var_rapid_time_dependent_potassium_current_Xr2_gate__alpha_xr2 = 3.0 / (1.0 + exp(((-60.0) - var_rapid_time_dependent_potassium_current_Xr2_gate__V) / 20.0)); // dimensionless
        const NekDouble var_rapid_time_dependent_potassium_current_Xr2_gate__beta_xr2 = 1.12 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr2_gate__V - 60.0) / 20.0)); // dimensionless
        const NekDouble var_rapid_time_dependent_potassium_current_Xr2_gate__tau_xr2 = 1.0 * var_rapid_time_dependent_potassium_current_Xr2_gate__alpha_xr2 * var_rapid_time_dependent_potassium_current_Xr2_gate__beta_xr2; // millisecond
        const NekDouble var_rapid_time_dependent_potassium_current_Xr2_gate__xr2_inf = 1.0 / (1.0 + exp((var_rapid_time_dependent_potassium_current_Xr2_gate__V + 88.0) / 24.0)); // dimensionless
        const NekDouble var_slow_time_dependent_potassium_current_Xs_gate__V = var_slow_time_dependent_potassium_current__V; // millivolt
        const NekDouble var_slow_time_dependent_potassium_current_Xs_gate__beta_xs = 1.0 / (1.0 + exp((var_slow_time_dependent_potassium_current_Xs_gate__V - 35.0) / 15.0)); // dimensionless
        const NekDouble var_slow_time_dependent_potassium_current_Xs_gate__alpha_xs = 1400.0 / sqrt(1.0 + exp((5.0 - var_slow_time_dependent_potassium_current_Xs_gate__V) / 6.0)); // dimensionless
        const NekDouble var_slow_time_dependent_potassium_current_Xs_gate__tau_xs = (1.0 * var_slow_time_dependent_potassium_current_Xs_gate__alpha_xs * var_slow_time_dependent_potassium_current_Xs_gate__beta_xs) + 80.0; // millisecond
        const NekDouble var_slow_time_dependent_potassium_current_Xs_gate__xs_inf = 1.0 / (1.0 + exp(((-5.0) - var_slow_time_dependent_potassium_current_Xs_gate__V) / 14.0)); // dimensionless
        const NekDouble var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V; // millivolt
        const NekDouble var_fast_sodium_current_m_gate__alpha_m = 1.0 / (1.0 + exp(((-60.0) - var_fast_sodium_current_m_gate__V) / 5.0)); // dimensionless
        const NekDouble var_fast_sodium_current_m_gate__beta_m = (0.1 / (1.0 + exp((var_fast_sodium_current_m_gate__V + 35.0) / 5.0))) + (0.1 / (1.0 + exp((var_fast_sodium_current_m_gate__V - 50.0) / 200.0))); // dimensionless
        const NekDouble var_fast_sodium_current_m_gate__tau_m = 1.0 * var_fast_sodium_current_m_gate__alpha_m * var_fast_sodium_current_m_gate__beta_m; // millisecond
        const NekDouble var_fast_sodium_current_m_gate__m_inf = 1.0 / pow(1.0 + exp(((-56.86) - var_fast_sodium_current_m_gate__V) / 9.03), 2.0); // dimensionless
        const NekDouble var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V; // millivolt
        const NekDouble var_fast_sodium_current_h_gate__h_inf = 1.0 / pow(1.0 + exp((var_fast_sodium_current_h_gate__V + 71.55) / 7.43), 2.0); // dimensionless
        const NekDouble var_fast_sodium_current_h_gate__beta_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? ((2.7 * exp(0.079 * var_fast_sodium_current_h_gate__V)) + (310000.0 * exp(0.3485 * var_fast_sodium_current_h_gate__V))) : (0.77 / (0.13 * (1.0 + exp((var_fast_sodium_current_h_gate__V + 10.66) / (-11.1))))); // per_millisecond
        const NekDouble var_fast_sodium_current_h_gate__alpha_h = (var_fast_sodium_current_h_gate__V < (-40.0)) ? (0.057 * exp((-(var_fast_sodium_current_h_gate__V + 80.0)) / 6.8)) : 0.0; // per_millisecond
        const NekDouble var_fast_sodium_current_h_gate__tau_h = 1.0 / (var_fast_sodium_current_h_gate__alpha_h + var_fast_sodium_current_h_gate__beta_h); // millisecond
        const NekDouble var_fast_sodium_current_j_gate__V = var_fast_sodium_current__V; // millivolt
        const NekDouble var_fast_sodium_current_j_gate__j_inf = 1.0 / pow(1.0 + exp((var_fast_sodium_current_j_gate__V + 71.55) / 7.43), 2.0); // dimensionless
        const NekDouble var_fast_sodium_current_j_gate__alpha_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((((((-25428.0) * exp(0.2444 * var_fast_sodium_current_j_gate__V)) - (6.948e-06 * exp((-0.04391) * var_fast_sodium_current_j_gate__V))) * (var_fast_sodium_current_j_gate__V + 37.78)) / 1.0) / (1.0 + exp(0.311 * (var_fast_sodium_current_j_gate__V + 79.23)))) : 0.0; // per_millisecond
        const NekDouble var_fast_sodium_current_j_gate__beta_j = (var_fast_sodium_current_j_gate__V < (-40.0)) ? ((0.02424 * exp((-0.01052) * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1378) * (var_fast_sodium_current_j_gate__V + 40.14)))) : ((0.6 * exp(0.057 * var_fast_sodium_current_j_gate__V)) / (1.0 + exp((-0.1) * (var_fast_sodium_current_j_gate__V + 32.0)))); // per_millisecond
        const NekDouble var_fast_sodium_current_j_gate__tau_j = 1.0 / (var_fast_sodium_current_j_gate__alpha_j + var_fast_sodium_current_j_gate__beta_j); // millisecond
        const NekDouble var_L_type_Ca_current_d_gate__V = var_L_type_Ca_current__V; // millivolt
        const NekDouble var_L_type_Ca_current_d_gate__alpha_d = (1.4 / (1.0 + exp(((-35.0) - var_L_type_Ca_current_d_gate__V) / 13.0))) + 0.25; // dimensionless
        const NekDouble var_L_type_Ca_current_d_gate__gamma_d = 1.0 / (1.0 + exp((50.0 - var_L_type_Ca_current_d_gate__V) / 20.0)); // millisecond
        const NekDouble var_L_type_Ca_current_d_gate__beta_d = 1.4 / (1.0 + exp((var_L_type_Ca_current_d_gate__V + 5.0) / 5.0)); // dimensionless
        const NekDouble var_L_type_Ca_current_d_gate__tau_d = (1.0 * var_L_type_Ca_current_d_gate__alpha_d * var_L_type_Ca_current_d_gate__beta_d) + var_L_type_Ca_current_d_gate__gamma_d; // millisecond
        const NekDouble var_L_type_Ca_current_d_gate__d_inf = 1.0 / (1.0 + exp(((-8.0) - var_L_type_Ca_current_d_gate__V) / 7.5)); // dimensionless
        const NekDouble var_L_type_Ca_current_f_gate__V = var_L_type_Ca_current__V; // millivolt
        const NekDouble var_L_type_Ca_current_f_gate__tau_f = (1102.5 * exp((-pow(var_L_type_Ca_current_f_gate__V + 27.0, 2.0)) / 225.0)) + (200.0 / (1.0 + exp((13.0 - var_L_type_Ca_current_f_gate__V) / 10.0))) + (180.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 30.0) / 10.0))) + 20.0; // millisecond
        const NekDouble var_L_type_Ca_current_f_gate__f_inf = 1.0 / (1.0 + exp((var_L_type_Ca_current_f_gate__V + 20.0) / 7.0)); // dimensionless
        const NekDouble var_L_type_Ca_current_f2_gate__V = var_L_type_Ca_current__V; // millivolt
        const NekDouble var_L_type_Ca_current_f2_gate__f2_inf = (0.67 / (1.0 + exp((var_L_type_Ca_current_f2_gate__V + 35.0) / 7.0))) + 0.33; // dimensionless
        const NekDouble var_L_type_Ca_current_f2_gate__tau_f2 = (562.0 * exp((-pow(var_L_type_Ca_current_f2_gate__V + 27.0, 2.0)) / 240.0)) + (31.0 / (1.0 + exp((25.0 - var_L_type_Ca_current_f2_gate__V) / 10.0))) + (80.0 / (1.0 + exp((var_L_type_Ca_current_f2_gate__V + 30.0) / 10.0))); // millisecond
        const NekDouble var_L_type_Ca_current_fCass_gate__Ca_ss = var_L_type_Ca_current__Ca_ss; // millimolar
        const NekDouble var_L_type_Ca_current_fCass_gate__tau_fCass = (80.0 / (1.0 + pow(var_L_type_Ca_current_fCass_gate__Ca_ss / 0.05, 2.0))) + 2.0; // millisecond
        const NekDouble var_L_type_Ca_current_fCass_gate__fCass_inf = (0.6 / (1.0 + pow(var_L_type_Ca_current_fCass_gate__Ca_ss / 0.05, 2.0))) + 0.4; // dimensionless
        const NekDouble var_transient_outward_current_s_gate__V = var_transient_outward_current__V; // millivolt
        const NekDouble var_transient_outward_current_s_gate__s_inf = 1.0 / (1.0 + exp((var_transient_outward_current_s_gate__V + s_inf_factor) / 5.0)); // dimensionless
        const NekDouble var_transient_outward_current_s_gate__tau_s = (s_tau_f1 * exp((-pow(var_transient_outward_current_s_gate__V + s_tau_f2, 2.0)) / s_tau_f3)) + s_tau_f4 + s_tau_f5*((5.0 / (1.0 + exp((var_transient_outward_current_s_gate__V - 20.0) / 5.0))) + 3.0); // millisecond
        const NekDouble var_transient_outward_current_r_gate__V = var_transient_outward_current__V; // millivolt
        const NekDouble var_transient_outward_current_r_gate__r_inf = 1.0 / (1.0 + exp((20.0 - var_transient_outward_current_r_gate__V) / 6.0)); // dimensionless
        const NekDouble var_transient_outward_current_r_gate__tau_r = (9.5 * exp((-pow(var_transient_outward_current_r_gate__V + 40.0, 2.0)) / 1800.0)) + 0.8; // millisecond
        const NekDouble var_calcium_dynamics__Ca_i = var_chaste_interface__calcium_dynamics__Ca_i; // millimolar
        const NekDouble var_calcium_dynamics__Ca_SR = var_chaste_interface__calcium_dynamics__Ca_SR; // millimolar
        const NekDouble var_calcium_dynamics__Ca_ss = var_chaste_interface__calcium_dynamics__Ca_ss; // millimolar
        const NekDouble var_calcium_dynamics__V_rel = 0.102; // per_millisecond
        const NekDouble var_calcium_dynamics__R_prime = var_chaste_interface__calcium_dynamics__R_prime; // dimensionless
        const NekDouble var_calcium_dynamics__k1_prime = 0.15; // per_millimolar2_per_millisecond
        const NekDouble var_calcium_dynamics__max_sr = 2.5; // dimensionless
        const NekDouble var_calcium_dynamics__EC = 1.5; // millimolar
        const NekDouble var_calcium_dynamics__min_sr = 1.0; // dimensionless
        const NekDouble var_calcium_dynamics__kcasr = var_calcium_dynamics__max_sr - ((var_calcium_dynamics__max_sr - var_calcium_dynamics__min_sr) / (1.0 + pow(var_calcium_dynamics__EC / var_calcium_dynamics__Ca_SR, 2.0))); // dimensionless
        const NekDouble var_calcium_dynamics__k1 = var_calcium_dynamics__k1_prime / var_calcium_dynamics__kcasr; // per_millimolar2_per_millisecond
        const NekDouble var_calcium_dynamics__k3 = 0.06; // per_millisecond
        const NekDouble var_calcium_dynamics__O = (var_calcium_dynamics__k1 * pow(var_calcium_dynamics__Ca_ss, 2.0) * var_calcium_dynamics__R_prime) / (var_calcium_dynamics__k3 + (var_calcium_dynamics__k1 * pow(var_calcium_dynamics__Ca_ss, 2.0))); // dimensionless
        const NekDouble var_calcium_dynamics__i_rel = var_calcium_dynamics__V_rel * var_calcium_dynamics__O * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_ss); // millimolar_per_millisecond
        const NekDouble var_calcium_dynamics__Vmax_up = 0.006375; // millimolar_per_millisecond
        const NekDouble var_calcium_dynamics__K_up = 0.00025; // millimolar
        const NekDouble var_calcium_dynamics__i_up = var_calcium_dynamics__Vmax_up / (1.0 + (pow(var_calcium_dynamics__K_up, 2.0) / pow(var_calcium_dynamics__Ca_i, 2.0))); // millimolar_per_millisecond
        const NekDouble var_calcium_dynamics__V_leak = 0.00036; // per_millisecond
        const NekDouble var_calcium_dynamics__i_leak = var_calcium_dynamics__V_leak * (var_calcium_dynamics__Ca_SR - var_calcium_dynamics__Ca_i); // millimolar_per_millisecond
        const NekDouble var_calcium_dynamics__V_xfer = 0.0038; // per_millisecond
        const NekDouble var_calcium_dynamics__i_xfer = var_calcium_dynamics__V_xfer * (var_calcium_dynamics__Ca_ss - var_calcium_dynamics__Ca_i); // millimolar_per_millisecond
        const NekDouble var_calcium_dynamics__k2_prime = 0.045; // per_millimolar_per_millisecond
        const NekDouble var_calcium_dynamics__k2 = var_calcium_dynamics__k2_prime * var_calcium_dynamics__kcasr; // per_millimolar_per_millisecond
        const NekDouble var_calcium_dynamics__k4 = 0.005; // per_millisecond
        const NekDouble var_calcium_dynamics__Buf_c = 0.2; // millimolar
        const NekDouble var_calcium_dynamics__K_buf_c = 0.001; // millimolar
        const NekDouble var_calcium_dynamics__Ca_i_bufc = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_c * var_calcium_dynamics__K_buf_c) / pow(var_calcium_dynamics__Ca_i + var_calcium_dynamics__K_buf_c, 2.0))); // dimensionless
        const NekDouble var_calcium_dynamics__K_buf_sr = 0.3; // millimolar
        const NekDouble var_calcium_dynamics__Buf_sr = 10.0; // millimolar
        const NekDouble var_calcium_dynamics__Ca_sr_bufsr = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_sr * var_calcium_dynamics__K_buf_sr) / pow(var_calcium_dynamics__Ca_SR + var_calcium_dynamics__K_buf_sr, 2.0))); // dimensionless
        const NekDouble var_calcium_dynamics__Buf_ss = 0.4; // millimolar
        const NekDouble var_calcium_dynamics__K_buf_ss = 0.00025; // millimolar
        const NekDouble var_calcium_dynamics__Ca_ss_bufss = 1.0 / (1.0 + ((var_calcium_dynamics__Buf_ss * var_calcium_dynamics__K_buf_ss) / pow(var_calcium_dynamics__Ca_ss + var_calcium_dynamics__K_buf_ss, 2.0))); // dimensionless
        const NekDouble var_calcium_dynamics__V_sr = 0.001094; // micrometre3
        const NekDouble var_calcium_dynamics__V_ss = 5.468e-05; // micrometre3
        const NekDouble var_calcium_dynamics__V_c = var_membrane__V_c; // micrometre3
        const NekDouble var_calcium_dynamics__F = var_membrane__F; // coulomb_per_millimole
        const NekDouble var_calcium_dynamics__Cm = var_membrane__Cm; // microF
        const NekDouble var_calcium_dynamics__i_CaL = var_L_type_Ca_current__i_CaL; // picoA_per_picoF
        const NekDouble var_calcium_dynamics__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa; // picoA_per_picoF
        const NekDouble var_calcium_dynamics__i_p_Ca = var_calcium_pump_current__i_p_Ca; // picoA_per_picoF
        const NekDouble var_calcium_dynamics__i_b_Ca = var_calcium_background_current__i_b_Ca; // picoA_per_picoF
        const NekDouble var_calcium_dynamics__d_Ca_i_d_environment__time = var_calcium_dynamics__Ca_i_bufc * (((((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_c) + var_calcium_dynamics__i_xfer) - ((1.0 * ((var_calcium_dynamics__i_b_Ca + var_calcium_dynamics__i_p_Ca) - (2.0 * var_calcium_dynamics__i_NaCa)) * var_calcium_dynamics__Cm) / (2.0 * 1.0 * var_calcium_dynamics__V_c * var_calcium_dynamics__F))); // 'millimole per litre per millisecond'
        const NekDouble var_calcium_dynamics__d_Ca_SR_d_environment__time = var_calcium_dynamics__Ca_sr_bufsr * (var_calcium_dynamics__i_up - (var_calcium_dynamics__i_rel + var_calcium_dynamics__i_leak)); // 'millimole per litre per millisecond'
        const NekDouble var_calcium_dynamics__d_Ca_ss_d_environment__time = var_calcium_dynamics__Ca_ss_bufss * (((((-1.0) * var_calcium_dynamics__i_CaL * var_calcium_dynamics__Cm) / (2.0 * 1.0 * var_calcium_dynamics__V_ss * var_calcium_dynamics__F)) + ((var_calcium_dynamics__i_rel * var_calcium_dynamics__V_sr) / var_calcium_dynamics__V_ss)) - ((var_calcium_dynamics__i_xfer * var_calcium_dynamics__V_c) / var_calcium_dynamics__V_ss)); // 'millimole per litre per millisecond'
        const NekDouble var_calcium_dynamics__d_R_prime_d_environment__time = ((-var_calcium_dynamics__k2) * var_calcium_dynamics__Ca_ss * var_calcium_dynamics__R_prime) + (var_calcium_dynamics__k4 * (1.0 - var_calcium_dynamics__R_prime)); // per_millisecond
        const NekDouble var_sodium_dynamics__F = var_membrane__F; // coulomb_per_millimole
        const NekDouble var_sodium_dynamics__Cm = var_membrane__Cm; // microF
        const NekDouble var_sodium_dynamics__V_c = var_membrane__V_c; // micrometre3
        const NekDouble var_sodium_dynamics__i_Na = var_fast_sodium_current__i_Na; // picoA_per_picoF
        const NekDouble var_sodium_dynamics__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa; // picoA_per_picoF
        const NekDouble var_sodium_dynamics__i_NaK = var_sodium_potassium_pump_current__i_NaK; // picoA_per_picoF
        const NekDouble var_sodium_dynamics__i_b_Na = var_sodium_background_current__i_b_Na; // picoA_per_picoF
        const NekDouble var_sodium_dynamics__d_Na_i_d_environment__time = (((-1.0) * (var_sodium_dynamics__i_Na + var_sodium_dynamics__i_b_Na + (3.0 * var_sodium_dynamics__i_NaK) + (3.0 * var_sodium_dynamics__i_NaCa))) / (1.0 * var_sodium_dynamics__V_c * var_sodium_dynamics__F)) * var_sodium_dynamics__Cm; // 'millimole per litre per millisecond'
        const NekDouble var_potassium_dynamics__F = var_membrane__F; // coulomb_per_millimole
        const NekDouble var_potassium_dynamics__Cm = var_membrane__Cm; // microF
        const NekDouble var_potassium_dynamics__V_c = var_membrane__V_c; // micrometre3
        const NekDouble var_potassium_dynamics__i_K1 = var_inward_rectifier_potassium_current__i_K1; // picoA_per_picoF
        const NekDouble var_potassium_dynamics__i_to = var_transient_outward_current__i_to; // picoA_per_picoF
        const NekDouble var_potassium_dynamics__i_NaK = var_sodium_potassium_pump_current__i_NaK; // picoA_per_picoF
        const NekDouble var_potassium_dynamics__i_Kr = var_rapid_time_dependent_potassium_current__i_Kr; // picoA_per_picoF
        const NekDouble var_potassium_dynamics__i_Ks = var_slow_time_dependent_potassium_current__i_Ks; // picoA_per_picoF
        const NekDouble var_potassium_dynamics__i_p_K = var_potassium_pump_current__i_p_K; // picoA_per_picoF
        const NekDouble var_potassium_dynamics__chaste_interface__chaste_membrane_capacitance = 1.0; // uF_per_cm2
        const NekDouble var_potassium_dynamics__i_Stim_converter = var_chaste_interface__membrane__i_Stim; // uA_per_cm2
        const NekDouble var_potassium_dynamics__i_Stim = var_potassium_dynamics__i_Stim_converter / var_potassium_dynamics__chaste_interface__chaste_membrane_capacitance; // picoA_per_picoF
        const NekDouble var_potassium_dynamics__d_K_i_d_environment__time = (((-1.0) * ((var_potassium_dynamics__i_K1 + var_potassium_dynamics__i_to + var_potassium_dynamics__i_Kr + var_potassium_dynamics__i_Ks + var_potassium_dynamics__i_p_K + var_potassium_dynamics__i_Stim) - (2.0 * var_potassium_dynamics__i_NaK))) / (1.0 * var_potassium_dynamics__V_c * var_potassium_dynamics__F)) * var_potassium_dynamics__Cm; // 'millimole per litre per millisecond'
        const NekDouble var_chaste_interface__calcium_dynamics__d_Ca_i_d_environment__time = var_calcium_dynamics__d_Ca_i_d_environment__time; // millimolar_per_millisecond
        const NekDouble var_chaste_interface__calcium_dynamics__d_Ca_SR_d_environment__time = var_calcium_dynamics__d_Ca_SR_d_environment__time; // millimolar_per_millisecond
        const NekDouble var_chaste_interface__calcium_dynamics__d_Ca_ss_d_environment__time = var_calcium_dynamics__d_Ca_ss_d_environment__time; // millimolar_per_millisecond
        const NekDouble var_chaste_interface__calcium_dynamics__d_R_prime_d_environment__time = var_calcium_dynamics__d_R_prime_d_environment__time; // per_millisecond
        const NekDouble var_chaste_interface__sodium_dynamics__d_Na_i_d_environment__time = var_sodium_dynamics__d_Na_i_d_environment__time; // millimolar_per_millisecond
        const NekDouble var_chaste_interface__potassium_dynamics__d_K_i_d_environment__time = var_potassium_dynamics__d_K_i_d_environment__time; // millimolar_per_millisecond
        const NekDouble d_dt_chaste_interface__calcium_dynamics__Ca_i = var_chaste_interface__calcium_dynamics__d_Ca_i_d_environment__time; // 'millimole per litre per millisecond'
        const NekDouble d_dt_chaste_interface__calcium_dynamics__Ca_SR = var_chaste_interface__calcium_dynamics__d_Ca_SR_d_environment__time; // 'millimole per litre per millisecond'
        const NekDouble d_dt_chaste_interface__calcium_dynamics__Ca_ss = var_chaste_interface__calcium_dynamics__d_Ca_ss_d_environment__time; // 'millimole per litre per millisecond'
        const NekDouble d_dt_chaste_interface__calcium_dynamics__R_prime = var_chaste_interface__calcium_dynamics__d_R_prime_d_environment__time; // per_millisecond
        const NekDouble d_dt_chaste_interface__sodium_dynamics__Na_i = var_chaste_interface__sodium_dynamics__d_Na_i_d_environment__time; // 'millimole per litre per millisecond'
        const NekDouble d_dt_chaste_interface__potassium_dynamics__K_i = var_chaste_interface__potassium_dynamics__d_K_i_d_environment__time; // 'millimole per litre per millisecond'

        const NekDouble var_membrane__i_K1 = var_inward_rectifier_potassium_current__i_K1; // picoA_per_picoF
        const NekDouble var_membrane__i_to = var_transient_outward_current__i_to; // picoA_per_picoF
        const NekDouble var_membrane__i_Kr = var_rapid_time_dependent_potassium_current__i_Kr; // picoA_per_picoF
        const NekDouble var_membrane__i_Ks = var_slow_time_dependent_potassium_current__i_Ks; // picoA_per_picoF
        const NekDouble var_membrane__i_CaL = var_L_type_Ca_current__i_CaL; // picoA_per_picoF
        const NekDouble var_membrane__i_NaK = var_sodium_potassium_pump_current__i_NaK; // picoA_per_picoF
        const NekDouble var_membrane__i_Na = var_fast_sodium_current__i_Na; // picoA_per_picoF
        const NekDouble var_membrane__i_b_Na = var_sodium_background_current__i_b_Na; // picoA_per_picoF
        const NekDouble var_membrane__i_NaCa = var_sodium_calcium_exchanger_current__i_NaCa; // picoA_per_picoF
        const NekDouble var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca; // picoA_per_picoF
        const NekDouble var_membrane__i_p_K = var_potassium_pump_current__i_p_K; // picoA_per_picoF
        const NekDouble var_membrane__i_p_Ca = var_calcium_pump_current__i_p_Ca; // picoA_per_picoF
        const NekDouble var_membrane__i_Stim_converter = var_chaste_interface__membrane__i_Stim; // uA_per_cm2
        const NekDouble var_membrane__chaste_interface__chaste_membrane_capacitance = 1.0; // uF_per_cm2
        const NekDouble var_membrane__i_Stim = var_membrane__i_Stim_converter / var_membrane__chaste_interface__chaste_membrane_capacitance; // picoA_per_picoF
        const NekDouble var_membrane__d_V_d_environment__time = ((-1.0) / 1.0) * (var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_CaL + var_membrane__i_NaK + var_membrane__i_Na + var_membrane__i_b_Na + var_membrane__i_NaCa + var_membrane__i_b_Ca + var_membrane__i_p_K + var_membrane__i_p_Ca + var_membrane__i_Stim); // 'millivolt per millisecond'
        const NekDouble var_chaste_interface__membrane__d_V_d_environment__time = var_membrane__d_V_d_environment__time; // ___units_1
        d_dt_chaste_interface__membrane__V = var_chaste_interface__membrane__d_V_d_environment__time; // 'millivolt per millisecond'
        outarray[0][i] = d_dt_chaste_interface__membrane__V;
        outarray[1][i] = var_rapid_time_dependent_potassium_current_Xr1_gate__xr1_inf;
        m_gates_tau[0][i] = var_rapid_time_dependent_potassium_current_Xr1_gate__tau_xr1;
        outarray[2][i] = var_rapid_time_dependent_potassium_current_Xr2_gate__xr2_inf;
        m_gates_tau[1][i] = var_rapid_time_dependent_potassium_current_Xr2_gate__tau_xr2;
        outarray[3][i] = var_slow_time_dependent_potassium_current_Xs_gate__xs_inf;
        m_gates_tau[2][i] = var_slow_time_dependent_potassium_current_Xs_gate__tau_xs;
        outarray[4][i] = var_fast_sodium_current_m_gate__m_inf;
        m_gates_tau[3][i] = var_fast_sodium_current_m_gate__tau_m;
        outarray[5][i] = var_fast_sodium_current_h_gate__h_inf;
        m_gates_tau[4][i] = var_fast_sodium_current_h_gate__tau_h;
        outarray[6][i] = var_fast_sodium_current_j_gate__j_inf;
        m_gates_tau[5][i] = var_fast_sodium_current_j_gate__tau_j;
        outarray[7][i] = var_L_type_Ca_current_d_gate__d_inf;
        m_gates_tau[6][i] = var_L_type_Ca_current_d_gate__tau_d;
        outarray[8][i] = var_L_type_Ca_current_f_gate__f_inf;
        m_gates_tau[7][i] = var_L_type_Ca_current_f_gate__tau_f;
        outarray[9][i] = var_L_type_Ca_current_f2_gate__f2_inf;
        m_gates_tau[8][i] = var_L_type_Ca_current_f2_gate__tau_f2;
        outarray[10][i] = var_L_type_Ca_current_fCass_gate__fCass_inf;
        m_gates_tau[9][i] = var_L_type_Ca_current_fCass_gate__tau_fCass;
        outarray[11][i] = var_transient_outward_current_s_gate__s_inf;
        m_gates_tau[10][i] = var_transient_outward_current_s_gate__tau_s;
        outarray[12][i] = var_transient_outward_current_r_gate__r_inf;
        m_gates_tau[11][i] = var_transient_outward_current_r_gate__tau_r;
        outarray[13][i] = d_dt_chaste_interface__calcium_dynamics__Ca_i;
        outarray[14][i] = d_dt_chaste_interface__calcium_dynamics__Ca_SR;
        outarray[15][i] = d_dt_chaste_interface__calcium_dynamics__Ca_ss;
        outarray[16][i] = d_dt_chaste_interface__calcium_dynamics__R_prime;
        outarray[17][i] = d_dt_chaste_interface__sodium_dynamics__Na_i;
        outarray[18][i] = d_dt_chaste_interface__potassium_dynamics__K_i;
    }

}


/**
 *
 */
void TenTusscher06::v_GenerateSummary(SummaryList& s)
{
    SolverUtils::AddSummaryItem(s, "Cell model", "TenTusscher06");
}


/**
 *
 */
void TenTusscher06::v_SetInitialConditions()
{
    // Values taken from CellML website
    switch (model_variant)
    {
    case eEpicardium:
        Vmath::Fill(m_nq, -85.23,       m_cellSol[0],   1);
        Vmath::Fill(m_nq, 0.00621,      m_cellSol[1],   1);
        Vmath::Fill(m_nq, 0.4712,       m_cellSol[2],   1);
        Vmath::Fill(m_nq, 0.0095,       m_cellSol[3],   1);
        Vmath::Fill(m_nq, 0.00172,      m_cellSol[4],   1);
        Vmath::Fill(m_nq, 0.7444,       m_cellSol[5],   1);
        Vmath::Fill(m_nq, 0.7045,       m_cellSol[6],   1);
        Vmath::Fill(m_nq, 3.373e-5,     m_cellSol[7],   1);
        Vmath::Fill(m_nq, 0.7888,       m_cellSol[8],   1);
        Vmath::Fill(m_nq, 0.9755,       m_cellSol[9],   1);
        Vmath::Fill(m_nq, 0.9953,       m_cellSol[10],  1);
        Vmath::Fill(m_nq, 0.999998,     m_cellSol[11],  1);
        Vmath::Fill(m_nq, 2.42e-8,      m_cellSol[12],  1);
        Vmath::Fill(m_nq, 0.000126,     m_cellSol[13],  1);
        Vmath::Fill(m_nq, 3.64,         m_cellSol[14],  1);
        Vmath::Fill(m_nq, 0.00036,      m_cellSol[15],  1);
        Vmath::Fill(m_nq, 0.9073,       m_cellSol[16],  1);
        Vmath::Fill(m_nq, 8.604,        m_cellSol[17],  1);
        Vmath::Fill(m_nq, 136.89,       m_cellSol[18],  1);
        break;
    case eEndocardium:
        Vmath::Fill(m_nq, -86.709,      m_cellSol[0],   1);
        Vmath::Fill(m_nq, 0.00448,      m_cellSol[1],   1);
        Vmath::Fill(m_nq, 0.476,        m_cellSol[2],   1);
        Vmath::Fill(m_nq, 0.0087,       m_cellSol[3],   1);
        Vmath::Fill(m_nq, 0.00155,      m_cellSol[4],   1);
        Vmath::Fill(m_nq, 0.7573,       m_cellSol[5],   1);
        Vmath::Fill(m_nq, 0.7225,       m_cellSol[6],   1);
        Vmath::Fill(m_nq, 3.164e-5,     m_cellSol[7],   1);
        Vmath::Fill(m_nq, 0.8009,       m_cellSol[8],   1);
        Vmath::Fill(m_nq, 0.9778,       m_cellSol[9],   1);
        Vmath::Fill(m_nq, 0.9953,       m_cellSol[10],  1);
        Vmath::Fill(m_nq, 0.3212,       m_cellSol[11],  1);
        Vmath::Fill(m_nq, 2.235e-8,     m_cellSol[12],  1);
        Vmath::Fill(m_nq, 0.00013,      m_cellSol[13],  1);
        Vmath::Fill(m_nq, 3.715,        m_cellSol[14],  1);
        Vmath::Fill(m_nq, 0.00036,      m_cellSol[15],  1);
        Vmath::Fill(m_nq, 0.9068,       m_cellSol[16],  1);
        Vmath::Fill(m_nq, 10.355,       m_cellSol[17],  1);
        Vmath::Fill(m_nq, 138.4,        m_cellSol[18],  1);
        break;
    case eMid:
        Vmath::Fill(m_nq, -85.423,      m_cellSol[0],   1);
        Vmath::Fill(m_nq, 0.0165,       m_cellSol[1],   1);
        Vmath::Fill(m_nq, 0.473,        m_cellSol[2],   1);
        Vmath::Fill(m_nq, 0.0174,       m_cellSol[3],   1);
        Vmath::Fill(m_nq, 0.00165,      m_cellSol[4],   1);
        Vmath::Fill(m_nq, 0.749,        m_cellSol[5],   1);
        Vmath::Fill(m_nq, 0.6788,       m_cellSol[6],   1);
        Vmath::Fill(m_nq, 3.288e-5,     m_cellSol[7],   1);
        Vmath::Fill(m_nq, 0.7026,       m_cellSol[8],   1);
        Vmath::Fill(m_nq, 0.9526,       m_cellSol[9],   1);
        Vmath::Fill(m_nq, 0.9942,       m_cellSol[10],  1);
        Vmath::Fill(m_nq, 0.999998,     m_cellSol[11],  1);
        Vmath::Fill(m_nq, 2.347e-8,     m_cellSol[12],  1);
        Vmath::Fill(m_nq, 0.000153,     m_cellSol[13],  1);
        Vmath::Fill(m_nq, 4.272,        m_cellSol[14],  1);
        Vmath::Fill(m_nq, 0.00042,      m_cellSol[15],  1);
        Vmath::Fill(m_nq, 0.8978,       m_cellSol[16],  1);
        Vmath::Fill(m_nq, 10.132,       m_cellSol[17],  1);
        Vmath::Fill(m_nq, 138.52,       m_cellSol[18],  1);
        break;
    case eIschemia:
        Vmath::Fill(m_nq, -85.23,       m_cellSol[0],   1);
        Vmath::Fill(m_nq, 0.00621,      m_cellSol[1],   1);
        Vmath::Fill(m_nq, 0.4712,       m_cellSol[2],   1);
        Vmath::Fill(m_nq, 0.0095,       m_cellSol[3],   1);
        Vmath::Fill(m_nq, 0.00172,      m_cellSol[4],   1);
        Vmath::Fill(m_nq, 0.7444,       m_cellSol[5],   1);
        Vmath::Fill(m_nq, 0.7045,       m_cellSol[6],   1);
        Vmath::Fill(m_nq, 3.373e-5,     m_cellSol[7],   1);
        Vmath::Fill(m_nq, 0.7888,       m_cellSol[8],   1);
        Vmath::Fill(m_nq, 0.9755,       m_cellSol[9],   1);
        Vmath::Fill(m_nq, 0.9953,       m_cellSol[10],  1);
        Vmath::Fill(m_nq, 0.999998,     m_cellSol[11],  1);
        Vmath::Fill(m_nq, 2.42e-8,      m_cellSol[12],  1);
        Vmath::Fill(m_nq, 0.000126,     m_cellSol[13],  1);
        Vmath::Fill(m_nq, 3.64,         m_cellSol[14],  1);
        Vmath::Fill(m_nq, 0.00036,      m_cellSol[15],  1);
        Vmath::Fill(m_nq, 0.9073,       m_cellSol[16],  1);
        Vmath::Fill(m_nq, 8.604,        m_cellSol[17],  1);
        Vmath::Fill(m_nq, 136.89,       m_cellSol[18],  1);
        break;
    }

}

}

