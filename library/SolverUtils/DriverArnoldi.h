///////////////////////////////////////////////////////////////////////////////
//
// File DriverArnoldi.h
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
// Description: Base Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERARNOLDI_H
#define NEKTAR_SOLVERUTILS_DRIVERARNOLDI_H

#include <SolverUtils/Driver.h>

namespace Nektar
{
namespace SolverUtils
{

/// Base class for the development of solvers.
class DriverArnoldi: public Driver
{
public:
    friend class MemoryManager<DriverArnoldi>;

    SOLVER_UTILS_EXPORT void ArnoldiSummary(std::ostream &out);

    SOLVER_UTILS_EXPORT inline const Array<OneD, const NekDouble> &GetMaskCoeff() const;

    SOLVER_UTILS_EXPORT inline const Array<OneD, const NekDouble> &GetMaskPhys() const;

protected:
    int       m_kdim;  /// Dimension of Krylov subspace
    int       m_nvec;  /// Number of vectors to test
    int       m_nits;  /// Maxmum number of iterations
    NekDouble m_evtol; /// Tolerance of iteratiosn
    NekDouble m_period;/// Period of time stepping algorithm
    bool      m_timeSteppingAlgorithm; /// underlying operator is time stepping

    int       m_infosteps; /// interval to dump information if required.

    int       m_nfields;
    NekDouble m_realShift;
    NekDouble m_imagShift;
    int       m_negatedOp;   /// Operator in solve call is negated

    Array<OneD, NekDouble> m_real_evl;
    Array<OneD, NekDouble> m_imag_evl;

    bool m_useMask;
    Array<OneD, NekDouble> m_maskCoeffs;
    Array<OneD, NekDouble> m_maskPhys;

    /// Constructor
    DriverArnoldi(const LibUtilities::SessionReaderSharedPtr pSession,
                  const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    virtual ~DriverArnoldi();

    /// Copy Arnoldi storage to fields.
    void CopyArnoldiArrayToField(Array<OneD, NekDouble> &array);

    /// Copy fields to Arnoldi storage.
    void CopyFieldToArnoldiArray(Array<OneD, NekDouble> &array);

    /// Copy the  forward field to the adjoint system in transient growth 
    /// calculations
    void CopyFwdToAdj();

    /// Write coefficients to file
    void WriteFld(std::string file, 
                  std::vector<Array<OneD, NekDouble> > coeffs);

    void WriteFld(std::string file, Array<OneD, NekDouble> coeffs);

    void WriteEvs(std::ostream &evlout, const int k,
                  const NekDouble real, const NekDouble imag,
                  NekDouble resid = NekConstants::kNekUnsetDouble,
                  bool DumpInverse = true);

    /// init mask
    void MaskInit();

    void GetUnmaskFunction(std::vector<std::vector<LibUtilities::EquationSharedPtr> > & unmaskfun);

    virtual void v_InitObject(std::ostream &out = std::cout);

    virtual  Array<OneD, NekDouble> v_GetRealEvl(void)
    {
        return m_real_evl;
    }

    virtual Array<OneD, NekDouble> v_GetImagEvl(void)
    {
        return m_imag_evl;
    }

};

inline const Array<OneD, const NekDouble> & DriverArnoldi::GetMaskCoeff() const
{
    return m_maskCoeffs;
}

inline const Array<OneD, const NekDouble> & DriverArnoldi::GetMaskPhys() const
{
    return m_maskPhys;
}
}
} //end of namespace

#endif //NEKTAR_SOLVERUTILS_DRIVERARNOLDI_H
