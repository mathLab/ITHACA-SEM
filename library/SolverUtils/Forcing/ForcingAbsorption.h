///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingAbsorption.h
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
// Description: Absorption layer forcing.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGABSORPTION
#define NEKTAR_SOLVERUTILS_FORCINGABSORPTION

#include <string>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace SolverUtils
{
    class ForcingAbsorption : public Forcing
    {
        public:
            friend class MemoryManager<ForcingAbsorption> ;

            /// Creates an instance of this class
            SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::weak_ptr<EquationSystem>      &pEquation,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields,
                    const TiXmlElement* pForce)
            {
                ForcingSharedPtr p = MemoryManager<ForcingAbsorption>::
                                        AllocateSharedPtr(pSession, pEquation);
                p->InitObject(pFields, pNumForcingFields, pForce);
                return p;
            }

            ///Name of the class
            static std::string className;

        protected:

            typedef bg::model::point<NekDouble, 3, bg::cs::cartesian> BPoint;
            typedef std::pair<BPoint, unsigned int>                   BPointPair;
            typedef bgi::rtree<BPointPair, bgi::rstar<16> >           BRTree;

            bool                                    m_hasRefFlow;
            bool                                    m_hasRefFlowTime;
            Array<OneD, Array<OneD, NekDouble> >    m_Absorption;
            Array<OneD, Array<OneD, NekDouble> >    m_Refflow;
            std::string                             m_funcNameTime;
            std::vector<unsigned int>               m_bRegions;
            std::shared_ptr<BRTree>               m_rtree;
  
            SOLVER_UTILS_EXPORT virtual void v_InitObject(
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                    const unsigned int& pNumForcingFields,
                    const TiXmlElement* pForce);

            SOLVER_UTILS_EXPORT virtual void v_Apply(
                    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                    const Array<OneD, Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble &time);

        private:
            ForcingAbsorption(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const std::weak_ptr<EquationSystem>      &pEquation);
            virtual ~ForcingAbsorption(void){};

            void CalcAbsorption(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const TiXmlElement *pForce);
    };
}
}
// Hui XU  21 Jul 2013 Created 
// Yumnah Mohamied May 2014 Modified and generalised. 
#endif
