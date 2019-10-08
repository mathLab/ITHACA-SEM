///////////////////////////////////////////////////////////////////////////////
//
// File FilterModalEnergy.h
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
// Description: Outputs values at the modal energy.

///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERMODALENERGY_H
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERMODALENERGY_H

#include <SolverUtils/Filters/Filter.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>

#include <MultiRegions/ExpList2D.h>     // for ExpList2D, etc
#include <MultiRegions/ExpList3D.h>     // for ExpList3D
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>


namespace Nektar
{
namespace SolverUtils
{
class FilterModalEnergy : public Filter
{
public:
    friend class MemoryManager<FilterModalEnergy>;

    // Creates an instance of this class
    static FilterSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr  &pSession,
        const std::weak_ptr<EquationSystem>       &pEquation,
        const ParamMap                              &pParams)
    {
        FilterSharedPtr p = MemoryManager<FilterModalEnergy>::
                                AllocateSharedPtr(pSession, pEquation, pParams);
        return p;
    }

    // Name of the class
    static std::string className;

    SOLVER_UTILS_EXPORT FilterModalEnergy(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem>      &pEquation,
        const ParamMap                             &pParams);
    SOLVER_UTILS_EXPORT virtual ~FilterModalEnergy();

protected:
    virtual void v_Initialise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble                           &time);
    virtual void v_Update(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble                           &time);
    virtual void v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble                           &time);
    virtual bool v_IsTimeDependent();
    NekDouble L2Error (
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        unsigned int                              field,
        const NekDouble                           &time);
    void SetUpBaseFields(
        SpatialDomains::MeshGraphSharedPtr        &mesh);
    void ImportFldBase(
        std::string                               pInfile);

private:
    enum MultiRegions::ProjectionType           m_projectionType;
    Array<OneD, MultiRegions::ExpListSharedPtr> m_base;
    LibUtilities::FieldIOSharedPtr              m_fld;

    // ID's of boundary regions where we want the forces
    std::vector<unsigned int>    m_boundaryRegionsIdList;
    // Determines if a given Boundary Region is in
    // m_boundaryRegionsIdList
    std::vector<bool>            m_boundaryRegionIsInList;
    unsigned int                 m_index;
    unsigned int                 m_outputFrequency;
    // plane to take history point from if using a homogeneous1D
    // expansion
    unsigned int                 m_outputPlane;
    bool                         m_isHomogeneous1D;
    bool                         m_isHomogeneous2D;
    bool                         m_PertEnergy;
    int                          m_npointsZ;
    int                          m_nproc;
    std::string                  m_outputFile;
    std::string                  m_EqTypeStr;
    std::ofstream                m_outputStream;
    LibUtilities::BasisSharedPtr m_homogeneousBasis;
    std::string                  m_BoundaryString;
    int                          m_nplanes;
    int                          m_NumQuadPointsError;
    bool                         m_SingleMode;
    bool                         m_HalfMode;
    bool                         m_MultipleModes;
    bool                         m_useFFT;
    NekDouble                    m_LhomZ;
    bool                         m_homogen_dealiasing;
};

}
}

#endif
