///////////////////////////////////////////////////////////////////////////////
//
// File: CollectionOptimisation.cpp
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
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <Collections/CollectionOptimisation.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{
namespace Collections
{

// static manager for Operator ImplementationMap
map<OpImpTimingKey,OperatorImpMap> CollectionOptimisation::m_opImpMap;

CollectionOptimisation::CollectionOptimisation(
        LibUtilities::SessionReaderSharedPtr pSession,
        ImplementationType defaultType)
{
    int i;
    map<ElmtOrder, ImplementationType> defaults, defaultsPhysDeriv;
    bool verbose  = (pSession.get()) &&
                    (pSession->DefinesCmdLineArgument("verbose")) &&
                    (pSession->GetComm()->GetRank() == 0);

    m_setByXml    = false;
    m_autotune    = false;
    m_maxCollSize = 0;
    m_defaultType = defaultType == eNoImpType ? eIterPerExp : defaultType;

    map<string, LibUtilities::ShapeType> elTypes;
    elTypes["S"] = LibUtilities::eSegment;
    elTypes["T"] = LibUtilities::eTriangle;
    elTypes["Q"] = LibUtilities::eQuadrilateral;
    elTypes["A"] = LibUtilities::eTetrahedron;
    elTypes["P"] = LibUtilities::ePyramid;
    elTypes["R"] = LibUtilities::ePrism;
    elTypes["H"] = LibUtilities::eHexahedron;

    // Set defaults for all element types.
    for (auto &it2 : elTypes)
    {
        defaults          [ElmtOrder(it2.second, -1)] = m_defaultType;
        defaultsPhysDeriv [ElmtOrder(it2.second, -1)] = m_defaultType;
    }

    if (defaultType == eNoImpType)
    {
        for (auto &it2 : elTypes)
        {
            // For 1<=N<=5 use StdMat otherwise IterPerExp or given default type
            for (int i = 1; i < 5; ++i)
            {
                defaults[ElmtOrder(it2.second, i)] = eStdMat;
            }

            // For 1<=N<=3 use SumFac otherwise NoCollection. Note that
            // default is not currently overwritten by given default
            // type
            defaultsPhysDeriv [ElmtOrder(it2.second, -1)] = eNoCollection;
            for (int i = 1; i < 3; ++i)
            {
                defaultsPhysDeriv[ElmtOrder(it2.second, i)] = eSumFac;
            }
        }
    }

    map<string, OperatorType> opTypes;
    for (i = 0; i < SIZE_OperatorType; ++i)
    {
        opTypes[OperatorTypeMap[i]] = (OperatorType)i;
        switch ((OperatorType)i)
        {
            case ePhysDeriv:
                m_global[(OperatorType)i] = defaultsPhysDeriv;
                break;
            default:
                m_global[(OperatorType)i] = defaults;
        }
    }

    map<string, ImplementationType> impTypes;
    for (i = 0; i < SIZE_ImplementationType; ++i)
    {
        impTypes[ImplementationTypeMap[i]] = (ImplementationType)i;
    }

    if(pSession.get()) // turn off file reader if dummy pointer is given
    {
        TiXmlDocument &doc = pSession->GetDocument();
        TiXmlHandle docHandle(&doc);
        TiXmlElement *master = docHandle.FirstChildElement("NEKTAR").Element();
        ASSERTL0(master, "Unable to find NEKTAR tag in file.");

        TiXmlElement *xmlCol = master->FirstChildElement("COLLECTIONS");

        // Check if user has specified some options
        if (xmlCol)
        {
            // Set the maxsize and default implementation type if provided
            const char *maxSize = xmlCol->Attribute("MAXSIZE");
            m_maxCollSize = (maxSize ? atoi(maxSize) : 0);

            const char *defaultImpl = xmlCol->Attribute("DEFAULT");
            m_defaultType = defaultType;

            // If user has specified a default impl type, autotuning
            // and set this default across all operators.
            if (defaultType == eNoImpType && defaultImpl)
            {
                const std::string collinfo = string(defaultImpl);
                m_autotune = boost::iequals(collinfo, "auto");

                if (!m_autotune)
                {
                    for(i = 1; i < Collections::SIZE_ImplementationType; ++i)
                    {
                        if(boost::iequals(collinfo,
                                Collections::ImplementationTypeMap[i]))
                        {
                            m_defaultType = (Collections::ImplementationType) i;
                            break;
                        }
                    }

                    ASSERTL0(i != Collections::SIZE_ImplementationType,
                         "Unknown default collection scheme: "+collinfo);

                    defaults.clear();
                    // Override default types
                    for (auto &it2 : elTypes)
                    {
                        defaults[ElmtOrder(it2.second, -1)] = m_defaultType;
                    }

                    for (i = 0; i < SIZE_OperatorType; ++i)
                    {
                        m_global[(OperatorType)i] = defaults;
                    }
                }
            }

            // Now process operator-specific implementation selections
            TiXmlElement *elmt = xmlCol->FirstChildElement();
            while (elmt)
            {
                m_setByXml = true;

                string tagname = elmt->ValueStr();

                ASSERTL0(boost::iequals(tagname, "OPERATOR"),
                        "Only OPERATOR tags are supported inside the "
                        "COLLECTIONS tag.");

                const char *attr = elmt->Attribute("TYPE");
                ASSERTL0(attr, "Missing TYPE in OPERATOR tag.");
                string opType(attr);

                ASSERTL0(opTypes.count(opType) > 0,
                        "Unknown OPERATOR type " + opType + ".");

                OperatorType ot = opTypes[opType];

                TiXmlElement *elmt2 = elmt->FirstChildElement();

                while (elmt2)
                {
                    string tagname = elmt2->ValueStr();
                    ASSERTL0(boost::iequals(tagname, "ELEMENT"),
                            "Only ELEMENT tags are supported inside the "
                            "OPERATOR tag.");

                    const char *attr = elmt2->Attribute("TYPE");
                    ASSERTL0(attr, "Missing TYPE in ELEMENT tag.");

                    string elType(attr);
                    auto it2 = elTypes.find(elType);
                    ASSERTL0(it2 != elTypes.end(),
                            "Unknown element type "+elType+" in ELEMENT "
                            "tag");

                    const char *attr2 = elmt2->Attribute("IMPTYPE");
                    ASSERTL0(attr2, "Missing IMPTYPE in ELEMENT tag.");
                    string impType(attr2);
                    ASSERTL0(impTypes.count(impType) > 0,
                            "Unknown IMPTYPE type " + impType + ".");

                    const char *attr3 = elmt2->Attribute("ORDER");
                    ASSERTL0(attr3, "Missing ORDER in ELEMENT tag.");
                    string order(attr3);

                    if (order == "*")
                    {
                        m_global[ot][ElmtOrder(it2->second, -1)]
                                     = impTypes[impType];
                    }
                    else
                    {
                        vector<unsigned int> orders;
                        ParseUtils::GenerateSeqVector(order, orders);

                        for (int i = 0; i < orders.size(); ++i)
                        {
                            m_global[ot][ElmtOrder(it2->second, orders[i])]
                                         = impTypes[impType];
                        }
                    }

                    elmt2 = elmt2->NextSiblingElement();
                }

                elmt = elmt->NextSiblingElement();
            }

            // Print out operator map
            if (verbose)
            {
                if (!m_setByXml && !m_autotune)
                {
                    cout << "Setting Collection optimisation using: "
                         << Collections::ImplementationTypeMap[m_defaultType]
                         << endl;
                }

                if (m_setByXml)
                {
                    for (auto &mIt : m_global)
                    {
                        cout << "Operator " << OperatorTypeMap[mIt.first]
                             << ":" << endl;

                        for (auto &eIt : mIt.second)
                        {
                            cout << "- "
                                 << LibUtilities::ShapeTypeMap[eIt.first.first]
                                 << " order " << eIt.first.second << " -> "
                                 << ImplementationTypeMap[eIt.second] << endl;
                        }
                    }
                }
            }
        }
    }
}

OperatorImpMap  CollectionOptimisation::GetOperatorImpMap(
        StdRegions::StdExpansionSharedPtr pExp)
{
    OperatorImpMap ret;
    ElmtOrder searchKey(pExp->DetShapeType(),
                        pExp->GetBasisNumModes(0));
    ElmtOrder defSearch(pExp->DetShapeType(), -1);

    for (auto &it : m_global)
    {
        ImplementationType impType;

        auto it2 = it.second.find(searchKey);

        if (it2 == it.second.end())
        {
            it2 = it.second.find(defSearch);
            if (it2 == it.second.end())
            {
                // Shouldn't be able to reach here.
                impType = eNoCollection;
            }
            else
            {
                impType = it2->second;
            }
        }
        else
        {
            impType = it2->second;
        }

        ret[it.first] = impType;
    }

    return ret;
}

OperatorImpMap CollectionOptimisation::SetWithTimings(
        vector<StdRegions::StdExpansionSharedPtr> pCollExp,
        OperatorImpMap &impTypes,
        bool verbose )
{
    boost::ignore_unused(impTypes);

    OperatorImpMap ret;

    StdRegions::StdExpansionSharedPtr pExp = pCollExp[0];

    // check to see if already defined for this expansion
    OpImpTimingKey OpKey(pExp,pCollExp.size(),pExp->GetNumBases());
    if(m_opImpMap.count(OpKey) != 0)
    {
        ret = m_opImpMap[OpKey];
        return ret;
    }

    int maxsize = pCollExp.size()*max(pExp->GetNcoeffs(),pExp->GetTotPoints());
    Array<OneD, NekDouble> inarray(maxsize,1.0);
    Array<OneD, NekDouble> outarray1(maxsize);
    Array<OneD, NekDouble> outarray2(maxsize);
    Array<OneD, NekDouble> outarray3(maxsize);

    LibUtilities::Timer t;

    if(verbose)
    {
        cout << "Collection Implemenation for "
             << LibUtilities::ShapeTypeMap[pExp->DetShapeType()] << " ( ";
        for(int i = 0; i < pExp->GetNumBases(); ++i)
        {
            cout << pExp->GetBasis(i)->GetNumModes() <<" ";
        }
        cout << ")" <<  " for ngeoms = " << pCollExp.size() << endl;
    }
    // set  up an array of collections
    CollectionVector coll;
    for(int imp = 1; imp < SIZE_ImplementationType; ++imp)
    {
        ImplementationType impType = (ImplementationType)imp;
        OperatorImpMap impTypes;
        for (int i = 0; i < SIZE_OperatorType; ++i)
        {
            OperatorType opType = (OperatorType)i;
            OperatorKey opKey(pCollExp[0]->DetShapeType(), opType, impType,
                              pCollExp[0]->IsNodalNonTensorialExp());

            if (GetOperatorFactory().ModuleExists(opKey))
            {
                impTypes[opType] = impType;
            }
            else
            {
                cout << "Note: Implementation does not exist: " << opKey << endl;
            }
        }

        Collection collloc(pCollExp,impTypes);
        coll.push_back(collloc);
    }

    // Determine the number of tests to do in one second
    Array<OneD, int> Ntest(SIZE_OperatorType);
    for(int i = 0; i < SIZE_OperatorType; ++i)
    {
        OperatorType OpType = (OperatorType)i;

        t.Start();
        coll[0].ApplyOperator(OpType,
                           inarray,
                           outarray1,
                           outarray2,
                           outarray3);
        t.Stop();

        NekDouble oneTest = t.TimePerTest(1);

        Ntest[i] = max((int)(0.25/oneTest),1);
    }

    Array<OneD, NekDouble> timing(SIZE_ImplementationType);
    // loop over all operators and determine fastest implementation
    for(int i = 0; i < SIZE_OperatorType; ++i)
    {
        OperatorType OpType = (OperatorType)i;

        // call collection implementation in thorugh ExpList.
        for (int imp = 0; imp < coll.size(); ++imp)
        {
            if (coll[imp].HasOperator(OpType))
            {
                t.Start();
                for(int n = 0; n < Ntest[i]; ++n)
                {
                    coll[imp].ApplyOperator(OpType,
                                      inarray,
                                      outarray1,
                                      outarray2,
                                      outarray3);
                }
                t.Stop();
                timing[imp] = t.TimePerTest(Ntest[i]);
            }
            else
            {
                timing[imp] = 1000.0;
            }
        }
        // determine optimal implementation. Note +1 to
        // remove NoImplementationType flag
        int minImp = Vmath::Imin(coll.size(),timing,1)+1;

        if(verbose)
        {
            cout << "\t " << OperatorTypeMap[i] << ": "
                 << ImplementationTypeMap[minImp] << "\t (";
            for(int j = 0; j < coll.size(); ++j)
            {
                if (timing[j] > 999.0)
                {
                    cout << "-";
                }
                else
                {
                    cout << timing[j] ;
                }
                if(j != coll.size()-1)
                {
                    cout <<", ";
                }
            }
            cout << ")" <<endl;
        }
        // could reset global map if reusing  method?
        //m_global[OpType][pExp->DetShapeType()] = (ImplementationType)minImp;
        // set up new map
        ret[OpType] = (ImplementationType)minImp;
    }

    // store map for use by another expansion.
    m_opImpMap[OpKey] = ret;
    return ret;
}

}
}
