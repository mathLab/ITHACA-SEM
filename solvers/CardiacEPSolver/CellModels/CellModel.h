///////////////////////////////////////////////////////////////////////////////
//
// File CellModel.h
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
// Description: Cell model base class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_CELLMODEL
#define NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_CELLMODEL

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SpatialDomains/SpatialData.h>

namespace Nektar
{
    // Forward declaration
    class CellModel;

    /// A shared pointer to an EquationSystem object
    typedef boost::shared_ptr<CellModel> CellModelSharedPtr;
    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the EquationSystem class.
    typedef LibUtilities::NekFactory< std::string, CellModel,
                const LibUtilities::SessionReaderSharedPtr&, int> CellModelFactory;
    CellModelFactory& GetCellModelFactory();

    /// Cell model base class.
    class CellModel
    {
    public:
        CellModel(const LibUtilities::SessionReaderSharedPtr& pSession, const int nq);
        virtual ~CellModel() {}

        virtual void Update(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time) = 0;

        virtual void v_PrintSummary(std::ostream &out) = 0;

    protected:
        /// Spatially varying parameters.
        SpatialDomains::SpatialParametersSharedPtr  m_spatialParameters;
        /// Number of physical points.
        int m_nq;
    };

}

#endif /* CELLMODEL_H_ */
