///////////////////////////////////////////////////////////////////////////////
//
// File LinearElasticSystem.h
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
// Description: Linear elasticity equation system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_LINEARELASTICSOLVER_EQUATIONSYSTEMS_LINELASSYSTEM_H
#define NEKTAR_SOLVERS_LINEARELASTICSOLVER_EQUATIONSYSTEMS_LINELASSYSTEM_H

#include <SolverUtils/EquationSystem.h>
#include <LinearElasticSolver/EquationSystems/CoupledAssemblyMap.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

/**
 * @brief Base class for linear elastic system.
 */
class LinearElasticSystem : public EquationSystem
{
public:
    /// Class may only be instantiated through the MemoryManager.
    friend class MemoryManager<LinearElasticSystem>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        EquationSystemSharedPtr p = MemoryManager<
            LinearElasticSystem>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    void BuildMatrixSystem();

    void SetStaticCondBlock(
        const int                              n,
        const LocalRegions::ExpansionSharedPtr exp,
        Array<TwoD, DNekMatSharedPtr>         &mat);

    DNekMatSharedPtr BuildLaplacianIJMatrix(
        const int                        k1,
        const int                        k2,
        const NekDouble                  scale,
        LocalRegions::ExpansionSharedPtr exp);

protected:
    LinearElasticSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph);

    /// Poisson ratio.
    NekDouble m_nu;
    /// Young's modulus.
    NekDouble m_E;
    /// Parameter dictating amount of thermal stress to add.
    NekDouble m_beta;
    /// Assembly map for the coupled (u,v,w) system.
    CoupledAssemblyMapSharedPtr m_assemblyMap;
    /// Schur complement boundary-boundary matrix.
    DNekScalBlkMatSharedPtr m_schurCompl;
    /// Matrix of elemental \f$ B^{-1}D \f$ components
    DNekScalBlkMatSharedPtr m_BinvD;
    /// Matrix of elemental \f$ C \f$ components
    DNekScalBlkMatSharedPtr m_C;
    /// Matrix of elemental \f$ D^{-1} \f$ components
    DNekScalBlkMatSharedPtr m_Dinv;
    /// Boundary maps for each of the fields.
    Array<OneD, Array<OneD, unsigned int> > m_bmap;
    /// Interior maps for each of the fields.
    Array<OneD, Array<OneD, unsigned int> > m_imap;
    /// Storage for the temperature terms.
    Array<OneD, Array<OneD, NekDouble> > m_temperature;
    /// Storage for the thermal stress terms.
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_stress;

    virtual void v_InitObject();
    virtual void v_GenerateSummary(SolverUtils::SummaryList& s);
    virtual void v_DoSolve();
    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables);
};

}

#endif
