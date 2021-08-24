///////////////////////////////////////////////////////////////////////////////
//
// File FilterCellHistoryPoints.h
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
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FilterCellHistoryPoints_H
#define NEKTAR_SOLVERUTILS_FILTERS_FilterCellHistoryPoints_H

#include <SolverUtils/Filters/FilterHistoryPoints.h>
#include <CardiacEPSolver/CellModels/CellModel.h>

namespace Nektar
{

class FilterCellHistoryPoints : public SolverUtils::FilterHistoryPoints
{
    public:
        friend class MemoryManager<FilterCellHistoryPoints>;

        /// Creates an instance of this class
        static SolverUtils::FilterSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr         &pSession,
            const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
            const ParamMap &pParams)
        {
            SolverUtils::FilterSharedPtr p
                        = MemoryManager<FilterCellHistoryPoints>
                            ::AllocateSharedPtr(pSession, pEquation, pParams);
            return p;
        }

        ///Name of the class
        static std::string className;

        FilterCellHistoryPoints(
            const LibUtilities::SessionReaderSharedPtr         &pSession,
            const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
            const ParamMap &pParams);
        ~FilterCellHistoryPoints();

        void SetCellModel(CellModelSharedPtr &pCellModel)
        {
            m_cell = pCellModel;
        }

    protected:
        virtual void v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);

        CellModelSharedPtr m_cell;

};

}

#endif /* NEKTAR_SOLVERUTILS_FILTERS_FILTERCHECKPOINT_H */
