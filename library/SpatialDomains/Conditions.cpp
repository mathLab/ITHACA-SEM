////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditions.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/Conditions.h>
#include <tinyxml.h>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{
/**
 * Constructor - collective on the session's communicator.
 */
BoundaryConditions::BoundaryConditions(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const MeshGraphSharedPtr &meshGraph)
    : m_meshGraph(meshGraph), m_session(pSession)

{
    Read(m_session->GetElement("Nektar/Conditions"));
}

BoundaryConditions::BoundaryConditions(void)
{
}

BoundaryConditions::~BoundaryConditions(void)
{
}

/**
 * Helper that turns a set into an array.
 */
Array<OneD, int> ToArray(const std::set<int> &set)
{
    Array<OneD, int> ans(set.size());
    auto it = set.begin(), end = set.end();
    for (int i = 0; it != end; ++it, ++i)
    {
        ans[i] = *it;
    }
    return ans;
}

/*
 * Helper function that effectively does an MPI_Allreduce for the
 * sets of boundary region IDs.
 *
 * Can't actually use an MPI_Allreduce because the sizes of the input
 * sets and output set are (in general) different.
 *
 * Instead, use a simple binary tree reduction and two MPI_Bcast calls.
 */
std::set<int> ShareAllBoundaryIDs(
    const BoundaryRegionCollection &boundaryRegions,
    LibUtilities::CommSharedPtr comm)
{
    // Turn the keys of boundaryRegions into set.
    std::set<int> ids;
    auto it = boundaryRegions.begin(), end = boundaryRegions.end();
    int i = 0;
    for (; it != end; ++it, ++i)
        ids.insert(it->first);

    int np = comm->GetSize();
    int ip = comm->GetRank();

    int half_size = 1;
    bool involved = true;
    while (involved && half_size < np)
    {
        if (ip & half_size)
        {
            // I'm sender

            // The receiver rank
            int receiver = ip - half_size;

            Array<OneD, int> idsArray = ToArray(ids);
            // Send my size (to alloc the reciever array)
            Array<OneD, int> sender_size(1);
            sender_size[0] = idsArray.size();
            comm->Send(receiver, sender_size);

            // Send my data
            comm->Send(receiver, idsArray);

            // Once we've sent, we're no longer involved.
            involved = false;
        }
        else
        {
            // I'm receiver

            // The sender rank
            int sender = ip + half_size;

            if (sender < np)
            {
                // Receive the size
                Array<OneD, int> sender_size(1);
                comm->Recv(sender, sender_size);

                // Receive the data
                Array<OneD, int> other_ids(sender_size[0]);
                comm->Recv(sender, other_ids);

                // Merge
                ids.insert(other_ids.begin(), other_ids.end());
            }
        }
        half_size *= 2;
    }

    // Bcast the size
    int nIds;
    if (ip == 0)
        nIds = ids.size();

    comm->Bcast(nIds, 0);

    // Bcast the data
    Array<OneD, int> idsArray;
    if (ip == 0)
        idsArray = ToArray(ids);
    else
        idsArray = Array<OneD, int>(nIds);

    comm->Bcast(idsArray, 0);

    return std::set<int>(idsArray.begin(), idsArray.end());
}

/**
 * Create a new communicator for each boundary region.
 * Collective on the session's communicator.
 */
void BoundaryConditions::CreateBoundaryComms()
{
    LibUtilities::CommSharedPtr comm = m_session->GetComm();

    if (comm->IsSerial())
    {
        // Do not try and generate a communicator if we have a serial
        // communicator. Arises with a FieldConvert communicator when
        // using --nparts in FieldConvert. Just set communicator to comm
        // in this case.
        for (auto &it : m_boundaryRegions)
        {
            m_boundaryCommunicators[it.first] = comm;
        }
        return;
    }

    std::set<int> allids = ShareAllBoundaryIDs(m_boundaryRegions, comm);

    for (auto &it : allids)
    {
        auto reg_it                = m_boundaryRegions.find(it);
        int this_rank_participates = (reg_it != m_boundaryRegions.end());
        LibUtilities::CommSharedPtr comm_region =
            comm->CommCreateIf(this_rank_participates);

        ASSERTL0(bool(comm_region) == bool(this_rank_participates),
                 "Rank should be in communicator but wasn't or is in "
                 "communicator but shouldn't be.");

        if (this_rank_participates)
        {
            m_boundaryCommunicators[reg_it->first] = comm_region;
        }
    }
}

/**
 * Collective on the session's communicator.
 */
void BoundaryConditions::Read(TiXmlElement *conditions)
{
    ASSERTL0(conditions, "Unable to find CONDITIONS tag in file.");

    TiXmlElement *boundaryRegions =
        conditions->FirstChildElement("BOUNDARYREGIONS");

    if (boundaryRegions)
    {
        ReadBoundaryRegions(conditions);
        CreateBoundaryComms();
        ReadBoundaryConditions(conditions);
    }
}

/**
 *
 */
void BoundaryConditions::ReadBoundaryRegions(TiXmlElement *conditions)
{
    // ensure boundary regions only read once per class definition
    if (m_boundaryRegions.size() != 0)
    {
        return;
    }

    TiXmlElement *boundaryRegions =
        conditions->FirstChildElement("BOUNDARYREGIONS");
    ASSERTL0(boundaryRegions, "Unable to find BOUNDARYREGIONS block.");

    // See if we have boundary regions defined.
    TiXmlElement *boundaryRegionsElement =
        boundaryRegions->FirstChildElement("B");

    while (boundaryRegionsElement)
    {
        /// All elements are of the form: "<B ID="#"> ... </B>", with
        /// ? being the element type.
        int indx;
        int err = boundaryRegionsElement->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");

        TiXmlNode *boundaryRegionChild = boundaryRegionsElement->FirstChild();
        // This is primarily to skip comments that may be present.
        // Comments appear as nodes just like elements.
        // We are specifically looking for text in the body
        // of the definition.
        while (boundaryRegionChild &&
               boundaryRegionChild->Type() != TiXmlNode::TINYXML_TEXT)
        {
            boundaryRegionChild = boundaryRegionChild->NextSibling();
        }

        ASSERTL0(boundaryRegionChild,
                 "Unable to read variable definition body.");
        std::string boundaryRegionStr =
            boundaryRegionChild->ToText()->ValueStr();

        std::string::size_type indxBeg =
            boundaryRegionStr.find_first_of('[') + 1;
        std::string::size_type indxEnd =
            boundaryRegionStr.find_last_of(']') - 1;

        ASSERTL0(indxBeg <= indxEnd,
                 (std::string("Error reading boundary region definition:") +
                  boundaryRegionStr)
                     .c_str());

        std::string indxStr =
            boundaryRegionStr.substr(indxBeg, indxEnd - indxBeg + 1);

        if (!indxStr.empty())
        {
            // Extract the composites from the string and return them in a list.
            BoundaryRegionShPtr boundaryRegion(
                MemoryManager<BoundaryRegion>::AllocateSharedPtr());

            ASSERTL0(m_boundaryRegions.count(indx) == 0,
                     "Boundary region " + indxStr +
                     " defined more than "
                     "once!");

            m_meshGraph->GetCompositeList(indxStr, *boundaryRegion);
            if (boundaryRegion->size() > 0)
            {
                m_boundaryRegions[indx] = boundaryRegion;
            }
        }

        boundaryRegionsElement =
            boundaryRegionsElement->NextSiblingElement("B");
    }
}

/**
 *
 */
void BoundaryConditions::ReadBoundaryConditions(TiXmlElement *conditions)
{
    // Protect against multiple reads.
    if (m_boundaryConditions.size() != 0)
    {
        return;
    }

    // Read REGION tags
    TiXmlElement *boundaryConditionsElement =
        conditions->FirstChildElement("BOUNDARYCONDITIONS");
    ASSERTL0(boundaryConditionsElement,
             "Boundary conditions must be specified.");

    TiXmlElement *regionElement =
        boundaryConditionsElement->FirstChildElement("REGION");

    // Read R (Robin), D (Dirichlet), N (Neumann), P (Periodic) C(Cauchy) tags
    while (regionElement)
    {
        BoundaryConditionMapShPtr boundaryConditions =
            MemoryManager<BoundaryConditionMap>::AllocateSharedPtr();

        int boundaryRegionID;
        int err = regionElement->QueryIntAttribute("REF", &boundaryRegionID);
        ASSERTL0(err == TIXML_SUCCESS,
                 "Error reading boundary region reference.");

        ASSERTL0(m_boundaryConditions.count(boundaryRegionID) == 0,
                 "Boundary region '" +
                 boost::lexical_cast<std::string>(boundaryRegionID) +
                 "' appears multiple times.");

        // Find the boundary region corresponding to this ID.
        std::string boundaryRegionIDStr;
        std::ostringstream boundaryRegionIDStrm(boundaryRegionIDStr);
        boundaryRegionIDStrm << boundaryRegionID;

        if (m_boundaryRegions.count(boundaryRegionID) == 0)
        {
            regionElement = regionElement->NextSiblingElement("REGION");
            continue;
        }

        ASSERTL0(m_boundaryRegions.count(boundaryRegionID) == 1,
                 "Boundary region " +
                 boost::lexical_cast<string>(boundaryRegionID) +
                 " not found");

        // Find the communicator that belongs to this ID
        LibUtilities::CommSharedPtr boundaryRegionComm =
            m_boundaryCommunicators[boundaryRegionID];

        TiXmlElement *conditionElement = regionElement->FirstChildElement();
        std::vector<std::string> vars  = m_session->GetVariables();

        while (conditionElement)
        {
            // Check type.
            std::string conditionType = conditionElement->Value();
            std::string attrData;
            bool isTimeDependent = false;

            // All have var specified, or else all variables are zero.
            TiXmlAttribute *attr = conditionElement->FirstAttribute();

            std::vector<std::string>::iterator iter;
            std::string attrName;

            attrData = conditionElement->Attribute("VAR");

            if (!attrData.empty())
            {
                iter = std::find(vars.begin(), vars.end(), attrData);
                ASSERTL0(
                    iter != vars.end(),
                    (std::string("Cannot find variable: ") + attrData).c_str());
            }

            if (conditionType == "N")
            {
                if (attrData.empty())
                {
                    // All variables are Neumann and are set to zero.
                    for (auto &varIter : vars)
                    {
                        BoundaryConditionShPtr neumannCondition(
                            MemoryManager<NeumannBoundaryCondition>::
                            AllocateSharedPtr(m_session, "00.0"));
                        (*boundaryConditions)[varIter] = neumannCondition;
                    }
                }
                else
                {
                    if (attr)
                    {
                        std::string equation, userDefined, filename;

                        while (attr)
                        {

                            attrName = attr->Name();

                            if (attrName == "VAR")
                            {
                                // if VAR do nothing
                            }
                            else if (attrName == "USERDEFINEDTYPE")
                            {
                                // Do stuff for the user defined attribute
                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "USERDEFINEDTYPE attribute must have "
                                         "associated value.");

                                // Suppose to go here?
                                m_session->SubstituteExpressions(attrData);

                                userDefined = attrData;
                                isTimeDependent =
                                    boost::iequals(attrData, "TimeDependent");
                            }
                            else if (attrName == "VALUE")
                            {
                                ASSERTL0(attrName == "VALUE",
                                         (std::string("Unknown attribute: ") +
                                          attrName)
                                             .c_str());

                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "VALUE attribute must be specified.");

                                m_session->SubstituteExpressions(attrData);

                                equation = attrData;
                            }
                            else if (attrName == "FILE")
                            {
                                ASSERTL0(attrName == "FILE",
                                         (std::string("Unknown attribute: ") +
                                          attrName)
                                             .c_str());

                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "FILE attribute must be specified.");

                                m_session->SubstituteExpressions(attrData);

                                filename = attrData;
                            }
                            else
                            {
                                ASSERTL0(false,
                                         (std::string("Unknown boundary "
                                                      "condition attribute: ") +
                                          attrName)
                                             .c_str());
                            }
                            attr = attr->Next();
                        }

                        BoundaryConditionShPtr neumannCondition(
                            MemoryManager<NeumannBoundaryCondition>::
                            AllocateSharedPtr(m_session, equation,
                                              userDefined, filename,
                                              boundaryRegionComm));
                        neumannCondition->SetIsTimeDependent(isTimeDependent);
                        (*boundaryConditions)[*iter] = neumannCondition;
                    }
                    else
                    {
                        // This variable's condition is zero.
                        BoundaryConditionShPtr neumannCondition(
                            MemoryManager<NeumannBoundaryCondition>::
                            AllocateSharedPtr(m_session, "0"));
                        (*boundaryConditions)[*iter] = neumannCondition;
                    }
                }
            }
            else if (conditionType == "D")
            {
                if (attrData.empty())
                {
                    // All variables are Dirichlet and are set to zero.
                    for (auto &varIter : vars)
                    {
                        BoundaryConditionShPtr dirichletCondition(
                            MemoryManager<DirichletBoundaryCondition>::
                            AllocateSharedPtr(m_session, "0"));
                        (*boundaryConditions)[varIter] = dirichletCondition;
                    }
                }
                else
                {
                    if (attr)
                    {
                        std::string equation, userDefined, filename;

                        while (attr)
                        {

                            attrName = attr->Name();

                            if (attrName == "VAR")
                            {
                                // if VAR do nothing
                            }
                            else if (attrName == "USERDEFINEDTYPE")
                            {

                                // Do stuff for the user defined attribute
                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "USERDEFINEDTYPE attribute must have "
                                         "associated value.");

                                m_session->SubstituteExpressions(attrData);

                                userDefined = attrData;
                                isTimeDependent =
                                    boost::iequals(attrData, "TimeDependent");
                            }
                            else if (attrName == "VALUE")
                            {
                                ASSERTL0(attrName == "VALUE",
                                         (std::string("Unknown attribute: ") +
                                          attrName)
                                             .c_str());

                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "VALUE attribute must have associated "
                                         "value.");

                                m_session->SubstituteExpressions(attrData);

                                equation = attrData;
                            }
                            else if (attrName == "FILE")
                            {
                                ASSERTL0(attrName == "FILE",
                                         (std::string("Unknown attribute: ") +
                                          attrName)
                                             .c_str());

                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "FILE attribute must be specified.");

                                m_session->SubstituteExpressions(attrData);

                                filename = attrData;
                            }
                            else
                            {
                                ASSERTL0(false,
                                         (std::string("Unknown boundary "
                                                      "condition attribute: ") +
                                          attrName)
                                             .c_str());
                            }
                            attr = attr->Next();
                        }

                        BoundaryConditionShPtr dirichletCondition(
                            MemoryManager<DirichletBoundaryCondition>::
                            AllocateSharedPtr(m_session, equation,
                                              userDefined, filename,
                                              boundaryRegionComm));
                        dirichletCondition->SetIsTimeDependent(isTimeDependent);
                        (*boundaryConditions)[*iter] = dirichletCondition;
                    }
                    else
                    {
                        // This variable's condition is zero.
                        BoundaryConditionShPtr dirichletCondition(
                            MemoryManager<DirichletBoundaryCondition>::
                            AllocateSharedPtr(m_session, "0"));
                        (*boundaryConditions)[*iter] = dirichletCondition;
                    }
                }
            }
            else if (conditionType == "R") // Read du/dn +  PRIMCOEFF u = VALUE
            {
                if (attrData.empty())
                {
                    // All variables are Robin and are set to zero.
                    for (auto &varIter : vars)
                    {
                        BoundaryConditionShPtr robinCondition(
                            MemoryManager<RobinBoundaryCondition>::
                            AllocateSharedPtr(m_session, "0", "0"));
                        (*boundaryConditions)[varIter] = robinCondition;
                    }
                }
                else
                {

                    if (attr)
                    {
                        std::string equation1, equation2, userDefined;
                        std::string filename;

                        bool primcoeffset = false;

                        while (attr)
                        {

                            attrName = attr->Name();

                            if (attrName == "VAR")
                            {
                                // if VAR do nothing
                            }
                            else if (attrName == "USERDEFINEDTYPE")
                            {

                                // Do stuff for the user defined attribute
                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "USERDEFINEDTYPE attribute must have "
                                         "associated value.");

                                m_session->SubstituteExpressions(attrData);
                                userDefined = attrData;
                                isTimeDependent =
                                    boost::iequals(attrData, "TimeDependent");
                            }
                            else if (attrName == "VALUE")
                            {

                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "VALUE attributes must have "
                                         "associated values.");

                                m_session->SubstituteExpressions(attrData);

                                equation1 = attrData;
                            }
                            else if (attrName == "PRIMCOEFF")
                            {

                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "PRIMCOEFF attributes must have "
                                         "associated values.");

                                m_session->SubstituteExpressions(attrData);

                                equation2 = attrData;

                                primcoeffset = true;
                            }
                            else if (attrName == "FILE")
                            {
                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "FILE attribute must be specified.");

                                m_session->SubstituteExpressions(attrData);

                                filename = attrData;
                            }
                            else
                            {
                                ASSERTL0(false,
                                         (std::string("Unknown boundary "
                                                      "condition attribute: ") +
                                          attrName)
                                             .c_str());
                            }
                            attr = attr->Next();
                        }

                        if (primcoeffset == false)
                        {
                            ASSERTL0(false, "PRIMCOEFF must be specified in a "
                                            "Robin boundary condition");
                        }
                        BoundaryConditionShPtr robinCondition(
                            MemoryManager<RobinBoundaryCondition>::
                            AllocateSharedPtr(
                                m_session, equation1, equation2,
                                userDefined, filename, boundaryRegionComm));
                        (*boundaryConditions)[*iter] = robinCondition;
                    }
                    else
                    {
                        // This variable's condition is zero.
                        BoundaryConditionShPtr robinCondition(
                            MemoryManager<RobinBoundaryCondition>::
                            AllocateSharedPtr(m_session, "0", "0"));
                        robinCondition->SetIsTimeDependent(isTimeDependent);
                        (*boundaryConditions)[*iter] = robinCondition;
                    }
                }
            }
            else if (conditionType == "P")
            {
                if (attrData.empty())
                {
                    ASSERTL0(false, "Periodic boundary conditions should "
                                    "be explicitely defined");
                }
                else
                {

                    if (attr)
                    {
                        std::string userDefined;
                        vector<unsigned int> periodicBndRegionIndex;
                        while (attr)
                        {
                            attrName = attr->Name();

                            if (attrName == "VAR")
                            {
                                // if VAR do nothing
                            }
                            else if (attrName == "USERDEFINEDTYPE")
                            {
                                // Do stuff for the user defined attribute
                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "USERDEFINEDTYPE attribute must have "
                                         "associated value.");

                                m_session->SubstituteExpressions(attrData);
                                userDefined = attrData;
                                isTimeDependent =
                                    boost::iequals(attrData, "TimeDependent");
                            }
                            else if (attrName == "VALUE")
                            {
                                attrData = attr->Value();
                                ASSERTL0(!attrData.empty(),
                                         "VALUE attribute must have associated "
                                         "value.");

                                int beg = attrData.find_first_of("[");
                                int end = attrData.find_first_of("]");
                                std::string periodicBndRegionIndexStr =
                                    attrData.substr(beg + 1, end - beg - 1);
                                ASSERTL0(
                                    beg < end,
                                    (std::string("Error reading periodic "
                                                 "boundary region definition "
                                                 "for boundary region: ") +
                                     boundaryRegionIDStrm.str())
                                        .c_str());

                                bool parseGood = ParseUtils::GenerateSeqVector(
                                    periodicBndRegionIndexStr.c_str(),
                                    periodicBndRegionIndex);

                                ASSERTL0(
                                    parseGood &&
                                    (periodicBndRegionIndex.size() == 1),
                                    (std::string(
                                        "Unable to read periodic boundary "
                                        "condition for boundary region: ") +
                                     boundaryRegionIDStrm.str())
                                        .c_str());
                            }
                            attr = attr->Next();
                        }
                        BoundaryConditionShPtr periodicCondition(
                            MemoryManager<PeriodicBoundaryCondition>::
                            AllocateSharedPtr(periodicBndRegionIndex[0],
                                              userDefined,
                                              boundaryRegionComm));
                        (*boundaryConditions)[*iter] = periodicCondition;
                    }
                    else
                    {
                        ASSERTL0(false, "Periodic boundary conditions should "
                                        "be explicitely defined");
                    }
                }
            }
            else if (conditionType == "C")
            {
                ASSERTL0(false, "Cauchy type boundary conditions not "
                                "implemented.");
            }

            conditionElement = conditionElement->NextSiblingElement();
        }

        m_boundaryConditions[boundaryRegionID] = boundaryConditions;
        regionElement = regionElement->NextSiblingElement("REGION");
    }
}
} // namespace SpatialDomains
} // namespace Nektar
