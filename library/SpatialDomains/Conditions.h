////////////////////////////////////////////////////////////////////////////////
//
//  File:  Conditions.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H
#define NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

#include <string>
#include <map>
#include <iostream>
#include <sstream>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

class TiXmlElement;
class TiXmlDocument;

class MeshGraph;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum BoundaryConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin,
            ePeriodic
        };

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(BoundaryConditionType type):
                m_boundaryConditionType(type)
            {
            }

            BoundaryConditionBase(BoundaryConditionType type, LibUtilities::Equation userDefined):
                m_boundaryConditionType(type),  m_userDefined(userDefined)
            {
            }

            virtual ~BoundaryConditionBase()
            {};

            BoundaryConditionType GetBoundaryConditionType() const
            {
                return m_boundaryConditionType;
            }

            void SetUserDefined(LibUtilities::Equation equation)
            {
                m_userDefined = equation;
            }

            LibUtilities::Equation GetUserDefined() const
            {
                return m_userDefined;
            }


        protected:
            BoundaryConditionType m_boundaryConditionType;
            LibUtilities::Equation m_userDefined;
        };


        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {

             DirichletBoundaryCondition(const std::string &eqn, const std::string &userDefined = std::string("NoUserDefined")):
            BoundaryConditionBase(eDirichlet, userDefined),
                m_dirichletCondition(eqn)
                {
                }

             LibUtilities::Equation m_dirichletCondition;
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
            NeumannBoundaryCondition(const std::string &eqn, const std::string &userDefined = std::string("NoUserDefined")):
                BoundaryConditionBase(eNeumann, userDefined),
                m_neumannCondition(eqn)
            {
            }

            LibUtilities::Equation m_neumannCondition;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            RobinBoundaryCondition( const std::string &a, const std::string &b, const std::string &userDefined = std::string("NoUserDefined")):
                BoundaryConditionBase(eRobin, userDefined),
                m_robinFunction(a), m_robinPrimitiveCoeff(b)
            {
            }
                // \frac{\partial {u}}{\partial{n}} +
                // m_robinPrimativeCoeff(x,y,z)*u = m_robinFunction(x,y,z)
            LibUtilities::Equation m_robinFunction;
            LibUtilities::Equation m_robinPrimitiveCoeff;
        };


        struct PeriodicBoundaryCondition : public BoundaryConditionBase
        {
            PeriodicBoundaryCondition(const unsigned int n):
                BoundaryConditionBase(ePeriodic),
                m_connectedBoundaryRegion(n)
            {
            }

            unsigned int m_connectedBoundaryRegion;
        };

        typedef std::map<int, Composite> BoundaryRegion;
        typedef boost::shared_ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef boost::shared_ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::vector<BoundaryRegionShPtr> BoundaryRegionCollection;

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef boost::shared_ptr<DirichletBoundaryCondition> DirichletBCShPtr;
        typedef boost::shared_ptr<NeumannBoundaryCondition>   NeumannBCShPtr;
        typedef boost::shared_ptr<RobinBoundaryCondition>     RobinBCShPtr;
        typedef std::map<std::string,BoundaryConditionShPtr>  BoundaryConditionMap;
        typedef boost::shared_ptr<BoundaryConditionMap>  BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;

        const static Array<OneD, BoundaryConditionShPtr> NullBoundaryConditionShPtrArray;

        class BoundaryConditions
        {
        public:
            SPATIAL_DOMAINS_EXPORT BoundaryConditions(LibUtilities::SessionReaderSharedPtr pSession, const MeshGraph *meshGraph);
            SPATIAL_DOMAINS_EXPORT ~BoundaryConditions();

            BoundaryRegionCollection &GetBoundaryRegions(void)
            {
                return m_boundaryRegions;
            }

            BoundaryConditionCollection &GetBoundaryConditions(void)
            {
                return m_boundaryConditions;
            }

            const std::string GetVariable(unsigned int indx)
            {
                return m_session->GetVariable(indx);
            }

        protected:
            /// The mesh graph to use for referencing geometry info.
            const MeshGraph                        *m_meshGraph;
            LibUtilities::SessionReaderSharedPtr    m_session;

            BoundaryRegionCollection                m_boundaryRegions;
            BoundaryConditionCollection             m_boundaryConditions;

        private:
            BoundaryConditions();

            /// Read segments (and general MeshGraph) given TiXmlDocument.
            void Read(TiXmlElement *conditions);

            void ReadBoundaryRegions(TiXmlElement *regions);
            void ReadBoundaryConditions(TiXmlElement *conditions);
        };

        typedef boost::shared_ptr<BoundaryConditions> BoundaryConditionsSharedPtr;
    }
}

#endif //NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

