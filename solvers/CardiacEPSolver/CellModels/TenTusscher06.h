///////////////////////////////////////////////////////////////////////////////
//
// File TenTusscher06Epi.h
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

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_TENTUSSCHER_PANFILOV_2006_CELL_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_TENTUSSCHER_PANFILOV_2006_CELL_H

#include <CardiacEPSolver/CellModels/CellModel.h>

namespace Nektar
{
    class TenTusscher06 : public CellModel
    {

    public:
        /// Creates an instance of this class
        static CellModelSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField)
        {
            return MemoryManager<TenTusscher06>::AllocateSharedPtr(pSession, pField);
        }

        /// Name of class
        static std::string className;

        /// Constructor
        TenTusscher06(const LibUtilities::SessionReaderSharedPtr& pSession, const MultiRegions::ExpListSharedPtr& pField);

        /// Desctructor
        virtual ~TenTusscher06() {}

    protected:
        virtual void v_Update(
               const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                     Array<OneD,        Array<OneD, NekDouble> >&outarray,
               const NekDouble time);

        /// Prints a summary of the model parameters.
        virtual void v_GenerateSummary(SummaryList& s);

        virtual void v_SetInitialConditions();

        NekDouble g_to;
        NekDouble g_Ks;
        NekDouble s_inf_factor;
        NekDouble s_tau_f1;
        NekDouble s_tau_f2;
        NekDouble s_tau_f3;
        NekDouble s_tau_f4;
        NekDouble s_tau_f5;
        NekDouble k_0;

        enum Variants {
        	eEpicardium,
        	eEndocardium,
        	eMid,
        	eIschemia
        };
        enum Variants model_variant;

        static std::string lookupIds[];
        static std::string def;
    };
}

#endif
