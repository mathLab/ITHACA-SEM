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
    std::string FentonKarma::lookupIds[9] = {
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
                    "CF3",  FentonKarma::eCF3)
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

        switch (model_variant) {
            case eBR:
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
            case eMBR:
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
            case eMLR1:
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
            case eGP:
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
            case eCF1:
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
            case eCF2a:
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
            case eCF2b:
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
            case eCF2c:
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
            case eCF3:
                g_fi_max     = 13.3333;
                tau_r        = 38;
                tau_si       = 127;
                tau_0        = 8.3;
                tau_v_plus   = 3.33;
                tau_v1_minus = 45;
                tau_v2_minus = 300;
                tau_w_plus   = 600;
                tau_w_minus  = 40;
                u_c          = 0.25;
                u_v          = 0.5;
                u_r          = 0.25;
                u_fi         = 0.25;
                u_csi        = 0.7;
                k1           = 60;
                k2           = 0;
                break;
        }

        tau_d  = C_m/g_fi_max;
        
        m_nvar = 3;
        
        // List gates and concentrations
        m_gates.push_back(1);
        m_gates.push_back(2);
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
        const NekDouble *u;
        const NekDouble *v;
        const NekDouble *w;
        NekDouble *u_new;
        NekDouble *v_new;
        NekDouble *w_new;
        NekDouble *v_tau;
        NekDouble *w_tau;

        // Temporary variables
        NekDouble J_fi, J_so, J_si, h1, h2, h3, alpha, beta;

        // Compute rates for each point in domain
        for (i = 0, u = &inarray[0][0], v = &inarray[1][0], w = &inarray[2][0],
                    u_new = &outarray[0][0], v_new = &outarray[1][0], w_new = &outarray[2][0],
                    v_tau = &m_gates_tau[0][0], w_tau = &m_gates_tau[1][0];
                i < n; ++i, ++u, ++v, ++w)
        {
            // Heavyside functions
            h1 = (*u < u_c) ? 0.0 : 1.0;
            h2 = (*u < u_v) ? 0.0 : 1.0;
            h3 = (*u < u_r) ? 0.0 : 1.0;

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

            // J_fi
            J_fi = -(*v)*h1*(1 - *u)*(*u - u_fi)/tau_d;

            // J_so
            // added extra (1-k2*v) term from Cherry&Fenton 2004
            J_so = (*u)*(1-h3)*(1-k2*(*v))/tau_0 + h3/tau_r;

            // J_si
            J_si = -(*w)*(1 + tanh(k1*(*u - u_csi)))/(2.0*tau_si);

            // u
            *u_new = -J_fi - J_so - J_si;
        }
    }
    
    void FentonKarma::v_PrintSummary(std::ostream &out)
    {
        out << "	Cell model      : FentonKarma" << std::endl;
        out << "    Cell model var. : " << lookupIds[model_variant] << std::endl;
    }
    
    
    void FentonKarma::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, 0.0,  m_cellSol[0],  1);
        Vmath::Fill(m_nq, 1.0,  m_cellSol[1],  1);
        Vmath::Fill(m_nq, 1.0,  m_cellSol[2],  1);
      
        
    }
}
