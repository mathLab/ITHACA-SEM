///////////////////////////////////////////////////////////////////////////////
//
// File CourtemancheRamirezNattel.h
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
// Description: Courtemanche-Ramirez-Nattel cell model
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_COURTEMANCHE_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_COURTEMANCHE_H

#include <CardiacEPSolver/CellModels/CellModel.h>

namespace Nektar
{
    class CourtemancheRamirezNattel98: public CellModel
    {

    public:
        /// Creates an instance of this class
        static CellModelSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField)
        {
            return MemoryManager<CourtemancheRamirezNattel98>
                                        ::AllocateSharedPtr(pSession, pField);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        CourtemancheRamirezNattel98(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField);

        /// Destructor
        virtual ~CourtemancheRamirezNattel98();

    protected:
        /// Computes the reaction terms $f(u,v)$ and $g(u,v)$.
        virtual void v_Update(
                const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        /// Prints a summary of the model parameters.
        virtual void v_PrintSummary(std::ostream &out);

        virtual void v_SetInitialConditions();

        virtual std::string v_GetCellVarName(unsigned int idx);

    private:
        NekDouble C_m;
        NekDouble g_Na;
        NekDouble g_K1;
        NekDouble g_to;
        NekDouble g_Kr;
        NekDouble g_Ks;
        NekDouble g_b_Na;
        NekDouble g_b_Ca;
        NekDouble g_Ca_L;
        NekDouble R;
        NekDouble T;
        NekDouble F;
        NekDouble Na_o;
        NekDouble K_o;
        NekDouble sigma;
        NekDouble K_i;
        NekDouble K_m_Na_i;
        NekDouble I_Na_K_max;
        NekDouble I_NaCa_max;
        NekDouble gamma;
        NekDouble Ca_o;
        NekDouble K_m_Na;
        NekDouble K_m_Ca;
        NekDouble K_sat;
        NekDouble I_p_Ca_max;
        NekDouble Trpn_max;
        NekDouble Km_Trpn;
        NekDouble Cmdn_max;
        NekDouble Csqn_max;
        NekDouble Km_Cmdn;
        NekDouble Km_Csqn;
        NekDouble NSR_I_up_max;
        NekDouble NSR_I_Ca_max;
        NekDouble NSR_K_up;
        NekDouble JSR_K_rel;
        NekDouble JSR_V_cell;
        NekDouble JSR_V_rel;
        NekDouble JSR_V_up;
        NekDouble tau_tr;
        NekDouble K_Q10;
        NekDouble V_i;
    };

}

#endif
