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

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <SpatialDomains/MeshGraph.h>


namespace Nektar
{
    struct OneD;

    namespace SpatialDomains
    {
        enum BoundaryConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin,
            ePeriodic,
            eJunction,
            eBifurcation,
            eMerging
        };

        enum BndUserDefinedType
        {
            eI,
            eMG,
            eHigh,
            eWall,
            eWALL,
            eWallViscous,
            eWALLVISCOUS,
            eSymmetry,
            eRinglebFlow,
            eTimeDependent,
            eRadiation,
            eIsentropicVortex,
            eCalcBC,
            eQinflow,
            eTerminal,
            eRterminal,
            eCRterminal,
            eRCRterminal,
            eInflowCFE,
            eOutflowCFE,
            eNoUserDefined
        };

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(
                    BoundaryConditionType type,
                    const std::string &userDefined = std::string("NoUserDefined")):
                m_boundaryConditionType(type)
            {
                std::map<const std::string, BndUserDefinedType>  known_type;
                known_type["H"] = eHigh;
                known_type["I"] = eI;
                known_type["MG"] = eMG;
                known_type["Wall"] = eWall;
                known_type["WallViscous"] = eWallViscous;
                known_type["Q-inflow"] = eQinflow;
                known_type["Terminal"] = eTerminal;
                known_type["R-terminal"] = eRterminal;
                known_type["CR-terminal"] = eCRterminal;
                known_type["RCR-terminal"] = eRCRterminal;
                known_type["WALL"] = eWALL;
                known_type["WALLVISCOUS"] = eWALLVISCOUS;
                known_type["CalcBC"] = eCalcBC;
                known_type["RinglebFlow"] = eRinglebFlow;
                known_type["Symmetry"] = eSymmetry;
                known_type["TimeDependent"] = eTimeDependent;
                known_type["Radiation"] = eRadiation;
                known_type["IsentropicVortex"] = eIsentropicVortex;
                known_type["NoUserDefined"] = eNoUserDefined;

                std::map<const std::string, BndUserDefinedType>::const_iterator it = known_type.find(userDefined);
                if (it != known_type.end())
                {
                    m_userDefined = it->second;
                }
                else
                {
                    //ASSERTL0(false, std::string("Unknown boundary condition user defined type [") + userDefined + std::string("]"));
                    m_userDefined = eNoUserDefined;
                }
            }

            virtual ~BoundaryConditionBase()
            {};

            BoundaryConditionType GetBoundaryConditionType() const
            {
                return m_boundaryConditionType;
            }

            void SetUserDefined(BndUserDefinedType type)
            {
                m_userDefined = type;
            }

            BndUserDefinedType GetUserDefined() const
            {
                return m_userDefined;
            }

            int m_parent;
            int m_daughter1;
            int m_daughter2;

            void SetJunction(int P, int D1)
            {
                m_parent = P;
				m_daughter1 = D1;
            }
			
			void SetBifurcation(int P, int D1, int D2)
            {
                m_parent = P;
				m_daughter1 = D1;
				m_daughter2 = D2;
            }
			
			void SetMerging(int P, int D1, int D2)
            {
                m_parent = P;
				m_daughter1 = D1;
				m_daughter2 = D2;
            }
			
			int GetParent() const
            {
                return m_parent;
            }
			
			int GetDaughter1() const
            {
                return m_daughter1;
            }
			
			int GetDaughter2() const
            {
                return m_daughter2;
            }


        protected:
            BoundaryConditionType m_boundaryConditionType;
            BndUserDefinedType    m_userDefined;
        };


        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {

             DirichletBoundaryCondition(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::string& eqn,
                    const std::string& userDefined = std::string("NoUserDefined"),
                    const std::string& filename=std::string("")):
                BoundaryConditionBase(eDirichlet, userDefined),
                m_dirichletCondition(pSession, eqn),
                m_filename(filename)
                {
                }

            LibUtilities::Equation m_dirichletCondition;
            std::string m_filename;
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
            NeumannBoundaryCondition(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::string& eqn,
                    const std::string& userDefined = std::string("NoUserDefined"),
                    const std::string& filename=std::string("")):
                BoundaryConditionBase(eNeumann, userDefined),
                m_neumannCondition(pSession, eqn),
                m_filename(filename)
            {
            }

            LibUtilities::Equation m_neumannCondition;
            std::string m_filename;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            RobinBoundaryCondition(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const std::string &a,
                    const std::string &b,
                    const std::string &userDefined = std::string("NoUserDefined"),
                    const std::string& filename=std::string("")):
                BoundaryConditionBase(eRobin, userDefined),
                m_robinFunction(pSession, a),
                m_robinPrimitiveCoeff(pSession, b),
                m_filename(filename)
            {
            }
                // \frac{\partial {u}}{\partial{n}} +
                // m_robinPrimativeCoeff(x,y,z)*u = m_robinFunction(x,y,z)
            LibUtilities::Equation m_robinFunction;
            LibUtilities::Equation m_robinPrimitiveCoeff;
            std::string m_filename;
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

        struct JunctionBoundaryCondition : public BoundaryConditionBase
        {
            JunctionBoundaryCondition( const int &P, const int &D1, const std::string &userDefined = std::string("NoUserDefined")):
            BoundaryConditionBase(eJunction, userDefined),
            m_parent(P), m_daughter1(D1)
            {
                SetJunction(P, D1);
            }
            int m_parent;
            int m_daughter1;
        };

        struct BifurcationBoundaryCondition : public BoundaryConditionBase
        {
            BifurcationBoundaryCondition( const int &P, const int &D1, const int &D2, const std::string &userDefined = std::string("NoUserDefined")):
            BoundaryConditionBase(eBifurcation, userDefined),
            m_parent(P), m_daughter1(D1), m_daughter2(D2)
            {
                SetBifurcation(P, D1, D2);
            }
            int m_parent;
            int m_daughter1;
            int m_daughter2;
        };

        struct MergingBoundaryCondition : public BoundaryConditionBase
        {
            MergingBoundaryCondition( const int &P, const int &D1, const int &D2, const std::string &userDefined = std::string("NoUserDefined")):
            BoundaryConditionBase(eMerging, userDefined),
            m_parent(P), m_daughter1(D1), m_daughter2(D2)
            {
                SetMerging(P, D1, D2);
            }
            int m_parent;
            int m_daughter1;
            int m_daughter2;
        };


        typedef std::map<int, Composite> BoundaryRegion;
        typedef boost::shared_ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef boost::shared_ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::vector<BoundaryRegionShPtr> BoundaryRegionCollection;

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef boost::shared_ptr<DirichletBoundaryCondition> DirichletBCShPtr;
        typedef boost::shared_ptr<NeumannBoundaryCondition>   NeumannBCShPtr;
        typedef boost::shared_ptr<RobinBoundaryCondition>     RobinBCShPtr;
        typedef boost::shared_ptr<JunctionBoundaryCondition>  JunctionBCShPtr;
        typedef boost::shared_ptr<BifurcationBoundaryCondition>  BifurcationBCShPtr;
        typedef boost::shared_ptr<MergingBoundaryCondition>   MergingBCShPtr;
        typedef std::map<std::string,BoundaryConditionShPtr>  BoundaryConditionMap;
        typedef boost::shared_ptr<BoundaryConditionMap>  BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;

        const static Array<OneD, BoundaryConditionShPtr> NullBoundaryConditionShPtrArray;

        class BoundaryConditions
        {
        public:
            SPATIAL_DOMAINS_EXPORT BoundaryConditions(const LibUtilities::SessionReaderSharedPtr &pSession, const MeshGraphSharedPtr &meshGraph);
            SPATIAL_DOMAINS_EXPORT ~BoundaryConditions();

            const BoundaryRegionCollection &GetBoundaryRegions(void) const
            {
                return m_boundaryRegions;
            }

            const BoundaryConditionCollection &GetBoundaryConditions(void) const
            {
                return m_boundaryConditions;
            }

            const std::string GetVariable(unsigned int indx)
            {
                return m_session->GetVariable(indx);
            }

        protected:
            /// The mesh graph to use for referencing geometry info.
            MeshGraphSharedPtr                      m_meshGraph;
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

