///////////////////////////////////////////////////////////////////////////////
//
// File FilterCellHistoryPoints.cpp
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

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <CardiacEPSolver/Filters/FilterCellHistoryPoints.h>

using namespace std;

namespace Nektar
{

std::string FilterCellHistoryPoints::className
    = SolverUtils::GetFilterFactory().RegisterCreatorFunction(
            "CellHistoryPoints", FilterCellHistoryPoints::create);

/**
 *
 */
FilterCellHistoryPoints::FilterCellHistoryPoints(
    const LibUtilities::SessionReaderSharedPtr         &pSession,
    const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
    const ParamMap &pParams) :
    FilterHistoryPoints(pSession, pEquation, pParams)
{
}


/**
 *
 */
FilterCellHistoryPoints::~FilterCellHistoryPoints()
{

}


/**
 *
 */
void FilterCellHistoryPoints::v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    // Only output every m_outputFrequency.
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }

    int j         = 0;
    int k         = 0;
    int numPoints = m_historyPoints.size();
    int numFields = m_cell->GetNumCellVariables();
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    Array<OneD, NekDouble> data(numPoints*numFields, 0.0);
    Array<OneD, NekDouble> gloCoord(3, 0.0);
    Array<OneD, NekDouble> physvals;
    Array<OneD, NekDouble> locCoord;
    int expId;
    int nppp = 0; // Number of points per plane

    // Pull out data values field by field
    for (j = 0; j < numFields; ++j)
    {
        k = 0;
        if(m_isHomogeneous1D)
        {
            for (auto &x : m_historyList)
            {
                locCoord = x.second;
                expId    = x.first->GetVid();
                nppp     = pFields[0]->GetPlane(0)->GetTotPoints();

                physvals = m_cell->GetCellSolution(j) + m_outputPlane*nppp
                            + pFields[j]->GetPhys_Offset(expId);

                // interpolate point can do with zero plane methods
                data[m_historyLocalPointMap[k++]*numFields+j]
                    = pFields[0]->GetExp(expId)->StdPhysEvaluate(
                                                    locCoord,physvals);
            }
        }
        else
        {
            for (auto &x : m_historyList)
            {
                locCoord = x.second;
                expId    = x.first->GetVid();

                physvals = m_cell->GetCellSolution(j)
                            + pFields[0]->GetPhys_Offset(expId);

                // interpolate point
                data[m_historyLocalPointMap[k++]*numFields+j]
                    = pFields[0]->GetExp(expId)->StdPhysEvaluate(
                                                    locCoord,physvals);
            }
        }
    }

    // Exchange history data
    // This could be improved to reduce communication but works for now
    vComm->AllReduce(data, LibUtilities::ReduceSum);

    // Only the root process writes out history data
    if (vComm->GetRank() == 0)
    {

        // Write data values point by point
        for (k = 0; k < m_historyPoints.size(); ++k)
        {
            m_outputStream.width(8);
            m_outputStream << setprecision(6) << time;
            for (int j = 0; j < numFields; ++j)
            {
                m_outputStream.width(25);
                m_outputStream << setprecision(16) << data[k*numFields+j];
            }
            m_outputStream << endl;
        }
    }
}

}
