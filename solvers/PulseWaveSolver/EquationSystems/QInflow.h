///////////////////////////////////////////////////////////////////////////////
//
// File QInflow.h
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
// Description: QInflow header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_QINFLOW_H
#define NEKTAR_QINFLOW_H

#include <string>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <PulseWaveSolver/EquationSystems/PulseWaveBoundary.h>

namespace Nektar
{
    // Forward declarations
    class QInflow;

    /// Pointer to a PulseWaveOutflow object.
    typedef std::shared_ptr<QInflow> QInflowSharedPtr;
    
    /// A global linear system.
    class QInflow : public PulseWaveBoundary
    {
    public:
        /// Creates an instance of this class
      static PulseWaveBoundarySharedPtr create(Array<OneD, MultiRegions::ExpListSharedPtr>& pVessel, 
                                               const LibUtilities::SessionReaderSharedPtr& pSession,
                                               PulseWavePressureAreaSharedPtr& pressureArea)
        {
            return MemoryManager<QInflow>::AllocateSharedPtr(pVessel,pSession,pressureArea);
        }

        /// Name of class
        static std::string className;
        
        QInflow(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel, 
                const LibUtilities::SessionReaderSharedPtr pSession,
                PulseWavePressureAreaSharedPtr pressureArea);

        virtual ~QInflow();
    protected:
        virtual void v_DoBoundary(
            const Array<OneD,const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &A_0,
            Array<OneD, Array<OneD, NekDouble> > &beta,
            const NekDouble time,
            int omega,int offset,int n);

    private:
        void Q_RiemannSolver(NekDouble Q,NekDouble A_r,NekDouble u_r,
                             NekDouble A_0, NekDouble beta,NekDouble &Au,
                             NekDouble &uu);
    };
}

#endif
