////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/BoundaryConditions.h,v $
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

#include <SpatialDomains/Equation.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SpatialDomains/MeshGraph.h>

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
            eRobin
        };

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(BoundaryConditionType type):
                m_BoundaryConditionType(type)
            {
            };

            virtual ~BoundaryConditionBase(){};

            BoundaryConditionType GetBoundaryConditionType()
            {
                return m_BoundaryConditionType;
            };

        protected:
            BoundaryConditionType m_BoundaryConditionType;
        };

        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {
            DirichletBoundaryCondition(const std::string &eqn):
                m_DirichletCondition(eqn),
                BoundaryConditionBase(eDirichlet)
            {
            };

            Equation<NekDouble> m_DirichletCondition;
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
            NeumannBoundaryCondition(const std::string &eqn):
                m_NeumannCondition(eqn),
                BoundaryConditionBase(eNeumann)
            {
            };

            Equation<NekDouble> m_NeumannCondition;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            RobinBoundaryCondition(const std::string &a, const std::string &b):
                m_a(a), m_b(b),
                BoundaryConditionBase(eRobin)
            {
            }

            // u = a(x,y,z) + b(x,y,z)*\frac{\partial{u}}{\partial{n}}
            Equation<NekDouble> m_a;
            Equation<NekDouble> m_b;
        };

        typedef std::map<std::string, NekDouble> ParamMap;
        typedef std::vector<std::string> Variable;
        typedef std::vector<Composite> BoundaryRegion;
        typedef ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::vector<BoundaryRegionShPtr> BoundaryRegionCollection;
        typedef ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef std::map<std::string,BoundaryConditionShPtr> BoundaryConditionMap;
        typedef ptr<BoundaryConditionMap> BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;
        typedef Equation<NekDouble> ForcingFunction;
        typedef ptr<ForcingFunction> ForcingFunctionShPtr;
        typedef ptr<const ForcingFunction> ConstForcingFunctionShPtr;
        typedef std::map<std::string, ForcingFunctionShPtr> ForcingFunctionsMap;
        typedef Equation<NekDouble> InitialCondition;
        typedef ptr<InitialCondition> InitialConditionShPtr;
        typedef ptr<const InitialCondition> ConstInitialConditionShPtr;
        typedef std::map<std::string, InitialConditionShPtr> InitialConditionsMap;

        class BoundaryConditions
        {
        public:
            BoundaryConditions(const MeshGraph *meshGraph);
            ~BoundaryConditions();

            void Read(std::string &infilename);
            void Read(TiXmlDocument &doc);

            NekDouble GetParameter(const std::string &parmName);

            BoundaryRegionCollection &GetBoundaryRegions(void) 
            {
                return m_BoundaryRegions;
            }

            BoundaryConditionCollection &GetBoundaryConditions(void)
            {
                return m_BoundaryConditions;
            }

            /// Get forcing function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            ConstForcingFunctionShPtr GetForcingFunction(int indx);

            /// Get forcing function based on name of variable.
            ConstForcingFunctionShPtr GetForcingFunction(const string &var);

            /// Get initial condition function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            ConstInitialConditionShPtr GetInitialCondition(int indx);

            /// Get initial condition function based on name of variable.
            ConstInitialConditionShPtr GetInitialCondition(const string &var);

        protected:
            void ReadParameters(TiXmlElement *conditions);
            void ReadVariables(TiXmlElement *conditions);
            void ReadBoundaryRegions(TiXmlElement *conditions);
            void ReadExpansionTypes(TiXmlElement *conditions);
            void ReadBoundaryConditions(TiXmlElement *conditions);
            void ReadForcingFunctions(TiXmlElement *conditions);
            void ReadInitialConditions(TiXmlElement *conditions);

            // Containers to hold conditions and associated data
            ParamMap m_Parameters;
            Variable m_Variables;
            BoundaryRegionCollection m_BoundaryRegions;
            BoundaryConditionCollection m_BoundaryConditions;
            ForcingFunctionsMap m_ForcingFunctions;
            InitialConditionsMap m_InitialConditions;

            /// The mesh graph to use for referencing geometry info.
            const MeshGraph *m_MeshGraph;

        private:
            BoundaryConditions();
        };
    };
}

#endif NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H
