///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveBoundary.h
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
// Description: PulseWaveBoundary header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_PULSEWAVEBOUNDARY_H
#define NEKTAR_PULSEWAVEBOUNDARY_H

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <PulseWaveSolver/EquationSystems/PulseWavePressureArea.h>

namespace Nektar
{
    class PulseWaveBoundary;
    typedef std::shared_ptr<PulseWaveBoundary>  PulseWaveBoundarySharedPtr;
    
    static PulseWaveBoundarySharedPtr NullPulseWaveBoundarySharedPtr;

    typedef LibUtilities::NekFactory< std::string, 
        PulseWaveBoundary, 
        Array<OneD, MultiRegions::ExpListSharedPtr>&, 
        const LibUtilities::SessionReaderSharedPtr&, 
        PulseWavePressureAreaSharedPtr& > BoundaryFactory;
    BoundaryFactory& GetBoundaryFactory();
    
    class PulseWaveBoundary
    {
    public:
        PulseWaveBoundary(Array<OneD, MultiRegions::ExpListSharedPtr> &pVessel,
                          const LibUtilities::SessionReaderSharedPtr &pSession,
                          PulseWavePressureAreaSharedPtr & pressureArea);

        virtual ~PulseWaveBoundary();

        inline void DoBoundary(
            const Array<OneD,const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &A_0,
            Array<OneD, Array<OneD, NekDouble> > &beta,
            const NekDouble time, 
            int omega, int offset,int n);

    protected:
        virtual void v_DoBoundary(
            const Array<OneD,const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &A_0,
            Array<OneD, Array<OneD, NekDouble> > &beta,
            const NekDouble time,
            int omega,int offset,int n) = 0;

        Array<OneD, MultiRegions::ExpListSharedPtr> m_vessels;
	LibUtilities::SessionReaderSharedPtr m_session;
        PulseWavePressureAreaSharedPtr m_pressureArea;

        NekDouble m_pext;
        NekDouble m_pout;
        NekDouble m_rho;


    private:
    };

    /**
     *
     */
    inline void PulseWaveBoundary::DoBoundary(
        const Array<OneD,const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &A_0,
        Array<OneD, Array<OneD, NekDouble> > &beta,
        const NekDouble time,
        int omega,int offset,int n)
    {
        v_DoBoundary(inarray,A_0,beta,time,omega,offset,n);
    }



}
#endif
