///////////////////////////////////////////////////////////////////////////////
//
// File FentonKarma.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) n6 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Fenton-Karma cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <CardiacEPSolver/CellModels/FentonKarma.h>

using namespace std;

namespace Nektar
{
    // Register cell model
    std::string FentonKarma::className =
            GetCellModelFactory().RegisterCreatorFunction(
                                                    "FentonKarma",
                                                    FentonKarma::create,
                                                    "Phenomenological Model.");

    // Register cell model variants
    std::string FentonKarma::lookupIds[] = {
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "BR",   FentonKarma::eBR),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "MBR",  FentonKarma::eMBR),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "MLR1", FentonKarma::eMLR1),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "GP",   FentonKarma::eGP),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "CF1",  FentonKarma::eCF1),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "CF2a",  FentonKarma::eCF2a),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "CF2b",  FentonKarma::eCF2b),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "CF2c",  FentonKarma::eCF2c),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "CF3a",  FentonKarma::eCF3a),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "CF3b",  FentonKarma::eCF3b),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set1",  FentonKarma::eFC2002Set1a),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set1a",  FentonKarma::eFC2002Set1a),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set1b",  FentonKarma::eFC2002Set1b),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set1c",  FentonKarma::eFC2002Set1c),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set1d",  FentonKarma::eFC2002Set1d),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set1e",  FentonKarma::eFC2002Set1e),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set2",  FentonKarma::eFC2002Set2),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set4",  FentonKarma::eFC2002Set4a),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set4a",  FentonKarma::eFC2002Set4a),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set4b",  FentonKarma::eFC2002Set4b),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set4c",  FentonKarma::eFC2002Set4c),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set4d",  FentonKarma::eFC2002Set4d),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set5",  FentonKarma::eFC2002Set5),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set6",  FentonKarma::eFC2002Set6),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set7",  FentonKarma::eFC2002Set7),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set8",  FentonKarma::eFC2002Set8),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "FC2002Set9",  FentonKarma::eFC2002Set9),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "Lawson",  FentonKarma::eLawson),
            LibUtilities::SessionReader::RegisterEnumValue("CellModelVariant",
                    "CAF",  FentonKarma::eCAF)
    };

    // Register default variant
    std::string FentonKarma::def =
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                    "CellModelVariant","BR");
    
    /**
     * Initialise Fenton-Karma model parameters.
     * k1 is k in the original model.
     * k2 is an additional parameter introduced in Cherry&Fenton 2004.
     * u_r and u_fi are introduced by Cherry&Fenton 2004 and is the same as
     * u_c in the original model.
     */
    FentonKarma::FentonKarma(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const MultiRegions::ExpListSharedPtr& pField)
        : CellModel(pSession, pField)
    {
        model_variant = pSession->GetSolverInfoAsEnum<FentonKarma::Variants>(
                                                            "CellModelVariant");

        C_m  =  1; // picoF
        V_0  = -85;
        V_fi =  15;

        switch (model_variant) {
            case eBR:   // (Fenton, Karma, Chaos 8(1), 1998)
                g_fi_max     = 4;
                tau_r        = 33.33;
                tau_si       = 29;
                tau_0        = 12.5;
                tau_v_plus   = 3.33;
                tau_v1_minus = 1250;
                tau_v2_minus = 19.6;
                tau_w_plus   = 870;
                tau_w_minus  = 41;
                u_c          = 0.13;
                u_v          = 0.04;
                u_r          = 0.13;
                u_fi         = 0.13;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0.0;
                break;
            case eMBR:  // (Fenton, Karma, Chaos 8(1), 1998)
                g_fi_max     = 4;
                tau_r        = 50;
                tau_si       = 44.84;
                tau_0        = 8.3;
                tau_v_plus   = 3.33;
                tau_v1_minus = 1000;
                tau_v2_minus = 19.2;
                tau_w_plus   = 667;
                tau_w_minus  = 11;
                u_c          = 0.13;
                u_v          = 0.04;
                u_r          = 0.13;
                u_fi         = 0.13;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0.0;
                break;
            case eMLR1: // (Fenton, Karma, Chaos 8(1), 1998)
                g_fi_max     = 5.8;
                tau_r        = 130;
                tau_si       = 127;
                tau_0        = 12.5;
                tau_v_plus   = 10;
                tau_v1_minus = 18.2;
                tau_v2_minus = 18.2;
                tau_w_plus   = 1020;
                tau_w_minus  = 80;
                u_c          = 0.13;
                u_v          = 0.0;
                u_r          = 0.13;
                u_fi         = 0.13;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0.0;
                break;
            case eGP:   // (Fenton, Karma, Chaos 8(1), 1998)
                g_fi_max     = 8.7;
                tau_r        = 25;
                tau_si       = 22.22;
                tau_0        = 12.5;
                tau_v_plus   = 10;
                tau_v1_minus = 333;
                tau_v2_minus = 40;
                tau_w_plus   = 1000;
                tau_w_minus  = 65;
                u_c          = 0.13;
                u_v          = 0.025;
                u_r          = 0.13;
                u_fi         = 0.13;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0.0;
                break;
            case eCF1:  // (Cherry, Fenton, Am J Physiol Heart 286, 2004)
                g_fi_max     = 6.6666;
                tau_r        = 12.5;
                tau_si       = 10;
                tau_0        = 1.5;
                tau_v_plus   = 10;
                tau_v1_minus = 350;
                tau_v2_minus = 80;
                tau_w_plus   = 562;
                tau_w_minus  = 48.5;
                u_c          = 0.25;
                u_v          = 0.001;
                u_r          = 0.25;
                u_fi         = 0.15;
                u_csi        = 0.2;
                k1           = 15;
                k2           = 0;
                break;
            case eCF2a: // (Cherry, Fenton, Am J Physiol Heart 286, 2004)
                g_fi_max     = 6.6666;
                tau_r        = 31;
                tau_si       = 26.5;
                tau_0        = 1.5;
                tau_v_plus   = 10;
                tau_v1_minus = 20;
                tau_v2_minus = 20;
                tau_w_plus   = 800;
                tau_w_minus  = 45;
                u_c          = 0.25;
                u_v          = 0.05;
                u_r          = 0.6;
                u_fi         = 0.11;
                u_csi        = 0.7;
                k1           = 10;
                k2           = 1;
                break;
            case eCF2b: // (Cherry, Fenton, Am J Physiol Heart 286, 2004)
                g_fi_max     = 6.6666;
                tau_r        = 31;
                tau_si       = 26.5;
                tau_0        = 1.5;
                tau_v_plus   = 10;
                tau_v1_minus = 100;
                tau_v2_minus = 20;
                tau_w_plus   = 800;
                tau_w_minus  = 45;
                u_c          = 0.25;
                u_v          = 0.05;
                u_r          = 0.6;
                u_fi         = 0.11;
                u_csi        = 0.7;
                k1           = 10;
                k2           = 1;
                break;
            case eCF2c: // (Cherry, Fenton, Am J Physiol Heart 286, 2004)
                g_fi_max     = 6.6666;
                tau_r        = 31;
                tau_si       = 26.5;
                tau_0        = 1.5;
                tau_v_plus   = 10;
                tau_v1_minus = 150;
                tau_v2_minus = 20;
                tau_w_plus   = 800;
                tau_w_minus  = 45;
                u_c          = 0.25;
                u_v          = 0.05;
                u_r          = 0.6;
                u_fi         = 0.11;
                u_csi        = 0.7;
                k1           = 10;
                k2           = 1;
                break;
            case eCF3a: // (Cherry, Fenton, Am J Physiol Heart 286, 2004)
                g_fi_max     = 13.3333;
                tau_r        = 38;
                tau_si       = 127;
                tau_0        = 8.3;
                tau_v_plus   = 3.33;
                tau_v1_minus = 45;
                tau_v2_minus = 300;
                tau_w_plus   = 600;
                tau_w_minus  = 40;
                tau_y_plus   = 1000;
                tau_y_minus  = 230;
                u_c          = 0.25;
                u_v          = 0.5;
                u_r          = 0.25;
                u_fi         = 0.25;
                u_csi        = 0.7;
                k1           = 60;
                k2           = 0;
                break;
            case eCF3b: // (Cherry, Fenton, Am J Physiol Heart 286, 2004)
                g_fi_max     = 13.3333;
                tau_r        = 38;
                tau_si       = 127;
                tau_0        = 8.3;
                tau_v_plus   = 3.33;
                tau_v1_minus = 20;
                tau_v2_minus = 300;
                tau_w_plus   = 600;
                tau_w_minus  = 40;
                tau_y_plus   = 1000;
                tau_y_minus  = 230;
                u_c          = 0.25;
                u_v          = 0.5;
                u_r          = 0.25;
                u_fi         = 0.25;
                u_csi        = 0.7;
                k1           = 60;
                k2           = 0;
                break;
            case eFC2002Set1a:  // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1/0.41; //2.44 - Fig 1A
                tau_r        = 50;
                tau_si       = 45;
                tau_0        = 1/0.12; //8.33
                tau_v_plus   = 3.33;
                tau_v1_minus = 1000;
                tau_v2_minus = 19.6;
                tau_w_plus   = 667;
                tau_w_minus  = 11;
                u_c          = 0.13;
                u_v          = 0.055;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set1b:  // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1/0.392; //2.55 - Fig 1B
                tau_r        = 50;
                tau_si       = 45;
                tau_0        = 1/0.12; //8.33
                tau_v_plus   = 3.33;
                tau_v1_minus = 1000;
                tau_v2_minus = 19.6;
                tau_w_plus   = 667;
                tau_w_minus  = 11;
                u_c          = 0.13;
                u_v          = 0.055;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set1c:  // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1/0.381; //2.62 - Fig 1C
                tau_r        = 50;
                tau_si       = 45;
                tau_0        = 1/0.12; //8.33
                tau_v_plus   = 3.33;
                tau_v1_minus = 1000;
                tau_v2_minus = 19.6;
                tau_w_plus   = 667;
                tau_w_minus  = 11;
                u_c          = 0.13;
                u_v          = 0.055;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set1d:  // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1/0.36; //2.77 - Fig 1D
                tau_r        = 50;
                tau_si       = 45;
                tau_0        = 1/0.12; //8.33
                tau_v_plus   = 3.33;
                tau_v1_minus = 1000;
                tau_v2_minus = 19.6;
                tau_w_plus   = 667;
                tau_w_minus  = 11;
                u_c          = 0.13;
                u_v          = 0.055;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set1e:  // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1/0.25;  // 4 - Fig 1E
                tau_r        = 50;
                tau_si       = 45;
                tau_0        = 1/0.12; //8.33
                tau_v_plus   = 3.33;
                tau_v1_minus = 1000;
                tau_v2_minus = 19.6;
                tau_w_plus   = 667;
                tau_w_minus  = 11;
                u_c          = 0.13;
                u_v          = 0.055;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set2:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 4;
                tau_r        = 190;
                tau_si       = 100000;
                tau_0        = 10;
                tau_v_plus   = 10;
                tau_v1_minus = 10;
                tau_v2_minus = 10;
                tau_w_plus   = 100000;
                tau_w_minus  = 100000;
                u_c          = 0.13;
                u_v          = 100;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 100;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set4a:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1.0/0.407;
                tau_r        = 34;
                tau_si       = 26.5;
                tau_0        = 9;
                tau_v_plus   = 3.33;
                tau_v1_minus = 5;       // Opposite convention
                tau_v2_minus = 15.6;    // to the paper
                tau_w_plus   = 350;
                tau_w_minus  = 80;
                u_c          = 0.15;
                u_v          = 0.04;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.45;
                k1           = 15;
                k2           = 0;
                break;
            case eFC2002Set4b:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1.0/0.41;
                tau_r        = 34;
                tau_si       = 26.5;
                tau_0        = 9;
                tau_v_plus   = 3.33;
                tau_v1_minus = 5;       // Opposite convention
                tau_v2_minus = 15.6;    // to the paper
                tau_w_plus   = 350;
                tau_w_minus  = 80;
                u_c          = 0.15;
                u_v          = 0.04;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.45;
                k1           = 15;
                k2           = 0;
                break;
            case eFC2002Set4c:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1.0/0.405;
                tau_r        = 34;
                tau_si       = 26.5;
                tau_0        = 9;
                tau_v_plus   = 3.33;
                tau_v1_minus = 5;       // Opposite convention
                tau_v2_minus = 15.6;    // to the paper
                tau_w_plus   = 350;
                tau_w_minus  = 80;
                u_c          = 0.15;
                u_v          = 0.04;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.45;
                k1           = 15;
                k2           = 0;
                break;
            case eFC2002Set4d:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 1.0/0.4;
                tau_r        = 34;
                tau_si       = 26.5;
                tau_0        = 9;
                tau_v_plus   = 3.33;
                tau_v1_minus = 5;       // Opposite convention
                tau_v2_minus = 15.6;    // to the paper
                tau_w_plus   = 350;
                tau_w_minus  = 80;
                u_c          = 0.15;
                u_v          = 0.04;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.45;
                k1           = 15;
                k2           = 0;
                break;
            case eFC2002Set5:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 2.76;
                tau_r        = 33.33;
                tau_si       = 29;
                tau_0        = 5;
                tau_v_plus   = 3.33;
                tau_v1_minus = 2;       // Opposite convention
                tau_v2_minus = 12;      // to the paper
                tau_w_plus   = 1000;
                tau_w_minus  = 100;
                u_c          = 0.13;
                u_v          = 0.14;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.70;
                k1           = 15;
                k2           = 0;
                break;
            case eFC2002Set6:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 2.53;
                tau_r        = 33.33;
                tau_si       = 29;
                tau_0        = 9;
                tau_v_plus   = 3.33;
                tau_v1_minus = 8;       // Opposite convention
                tau_v2_minus = 9;       // to the paper
                tau_w_plus   = 250;
                tau_w_minus  = 60;
                u_c          = 0.13;
                u_v          = 0.04;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.50;
                k1           = 15;
                k2           = 0;
                break;
            case eFC2002Set7:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 4;
                tau_r        = 100;
                tau_si       = 100000;
                tau_0        = 12;
                tau_v_plus   = 10;
                tau_v1_minus = 7;       // Opposite convention
                tau_v2_minus = 7;       // to the paper
                tau_w_plus   = 100000;
                tau_w_minus  = 100000;
                u_c          = 0.13;
                u_v          = 100;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 100;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set8:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 2.2;
                tau_r        = 33.25;
                tau_si       = 29;
                tau_0        = 12.5;
                tau_v_plus   = 13.03;
                tau_v1_minus = 1250;    // Opposite convention
                tau_v2_minus = 19.6;    // to the paper
                tau_w_plus   = 800;
                tau_w_minus  = 40;
                u_c          = 0.13;
                u_v          = 0.04;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
            case eFC2002Set9:   // (Fenton, Cherry, Chaos 12(852), 2002)
                g_fi_max     = 4;
                tau_r        = 28;
                tau_si       = 29;
                tau_0        = 12.5;
                tau_v_plus   = 3.33;
                tau_v1_minus = 2;       // Opposite convention
                tau_v2_minus = 15;      // to the paper
                tau_w_plus   = 670;
                tau_w_minus  = 61;
                u_c          = 0.13;
                u_v          = 0.05;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.45;
                k1           = 10;
                k2           = 0;
                break;
            case eLawson:       // (Lawson, Front. Physiol, 2018)
                g_fi_max     = 3;
                tau_r        = 1/0.02;
                tau_si       = 1/0.0223;
                tau_0        = 1/0.12;
                tau_v_plus   = 3.33;
                tau_v1_minus = 1000;
                tau_v2_minus = 19.6;
                tau_w_plus   = 667;
                tau_w_minus  = 11;
                u_c          = 0.13;
                u_v          = 0.055;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
            case eCAF:
                g_fi_max     = 8;
                tau_r        = 70;
                tau_si       = 114;
                tau_0        = 32.5; //8.33
                tau_v_plus   = 5.75;
                tau_v1_minus = 60;
                tau_v2_minus = 82.5;
                tau_w_plus   = 300;
                tau_w_minus  = 400;
                u_c          = 0.16;
                u_v          = 0.04;
                u_r          = u_c;
                u_fi         = u_c;
                u_csi        = 0.85;
                k1           = 10;
                k2           = 0;
                break;
        }

        tau_d  = C_m/g_fi_max;
        
        isCF3 = (model_variant == eCF3a || model_variant == eCF3b);

        // List gates and concentrations
        m_gates.push_back(1);
        m_gates.push_back(2);
        m_nvar = 3;

        // Cherry-Fenton 2004 Model 3 has extra gating variable
        if (isCF3)
        {
            m_gates.push_back(3);
            m_nvar = 4;
        }

        ASSERTL0(!isCF3, "Cherry-Fenton model 3 not implemented yet.");
    }
    
    
    
    /**
     *
     */
    FentonKarma::~FentonKarma()
    {
        
    }
    
    
    
    void FentonKarma::v_Update(
                                 const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                 Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                 const NekDouble time )
    {
        ASSERTL0(inarray.get() != outarray.get(),
                 "Must have different arrays for input and output.");
        
        // Variables
        //  0   u    membrane potential
        //  1   v    v gate
        //  2   w    w gate
        int n = m_nq;
        int i = 0;

        // Declare pointers
        const NekDouble *u = &inarray[0][0];
        const NekDouble *v = &inarray[1][0];
        const NekDouble *w = &inarray[2][0];
        const NekDouble *y = isCF3 ? &inarray[3][0] : 0;
        NekDouble *u_new   = &outarray[0][0];
        NekDouble *v_new   = &outarray[1][0];
        NekDouble *w_new   = &outarray[2][0];
        //NekDouble *y_new   = isCF3 ? &outarray[3][0] : 0;
        NekDouble *v_tau   = &m_gates_tau[0][0];
        NekDouble *w_tau   = &m_gates_tau[1][0];
        //NekDouble *y_tau   = isCF3 ? &m_gates_tau[2][0] : 0;

        // Temporary variables
        NekDouble J_fi, J_so, J_si, h1, h2, h3, alpha, beta;

        double V;
        // Compute rates for each point in domain
        for (i = 0; i < n; ++i)
        {
            // *u is dimensional
            // V is non-dimensional for what follows
            V = (*u - V_0)/(V_fi - V_0);

            // Heavyside functions
            h1 = (V < u_c) ? 0.0 : 1.0;
            h2 = (V < u_v) ? 0.0 : 1.0;
            h3 = (V < u_r) ? 0.0 : 1.0;

            // w-gate
            alpha = (1-h1)/tau_w_minus;
            beta = h1/tau_w_plus;
            *w_tau = 1.0 / (alpha + beta);
            *w_new = alpha * (*w_tau);

            // v-gate
            alpha = (1-h1)/(h2*tau_v1_minus + (1-h2)*tau_v2_minus);
            beta = h1/tau_v_plus;
            *v_tau = 1.0 / (alpha + beta);
            *v_new = alpha * (*v_tau);

            // y-gate
            if (isCF3)
            {
                // TODO: implementation for y_tau and y_new
            }

            // J_fi
            J_fi = -(*v)*h1*(1 - V)*(V - u_c)/tau_d;

            // J_so
            // added extra (1-k2*v) term from Cherry&Fenton 2004
            J_so = V*(1-h3)*(1-k2*(*v))/tau_0 +
                    (isCF3 ? h3*V*(*y)/tau_r : h3/tau_r);

            // J_si
            J_si = -(*w)*(1 + tanh(k1*(V - u_csi)))/(2.0*tau_si);

            // u
            *u_new = -J_fi - J_so - J_si;
            *u_new *= C_m*(V_fi - V_0);

            ++u, ++v, ++w, ++u_new, ++v_new, ++w_new, ++v_tau, ++w_tau;
        }
    }
    
    void FentonKarma::v_GenerateSummary(SummaryList& s)
    {
        SolverUtils::AddSummaryItem(s, "Cell model", "Fenton Karma");
        SolverUtils::AddSummaryItem(s, "Cell model var.", lookupIds[model_variant]);
    }
    
    
    void FentonKarma::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, 0.0,  m_cellSol[0],  1);
        Vmath::Fill(m_nq, 1.0,  m_cellSol[1],  1);
        Vmath::Fill(m_nq, 1.0,  m_cellSol[2],  1);
        if (isCF3)
        {
            Vmath::Fill(m_nq, 0.1,  m_cellSol[2],  1);
        }
        
    }

    std::string FentonKarma::v_GetCellVarName(unsigned int idx)
    {
        switch (idx)
        {
            case 0:  return "u";
            case 1:  return "v";
            case 2:  return "w";
            case 3:  return "y";
            default: return "unknown";
        }
    }
}
