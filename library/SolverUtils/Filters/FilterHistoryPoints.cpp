///////////////////////////////////////////////////////////////////////////////
//
// File FilterHistoryPoints.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Outputs values at specific points during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <SolverUtils/Filters/FilterHistoryPoints.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterHistoryPoints::className = GetFilterFactory().RegisterCreatorFunction("HistoryPoints", FilterHistoryPoints::create);

        /**
         *
         */
        FilterHistoryPoints::FilterHistoryPoints(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            if (pParams.find("OutputFile") == pParams.end())
            {
                m_outputFile = m_session->GetSessionName();
            }
            else
            {
                ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                         "Missing parameter 'OutputFile'.");
                m_outputFile = pParams.find("OutputFile")->second;
            }
            if (!(m_outputFile.length() >= 4 
                  && m_outputFile.substr(m_outputFile.length() - 4) == ".his"))
            {
                m_outputFile += ".his";
            }

            if (pParams.find("OutputFrequency") == pParams.end())
            {
                m_outputFrequency = 1;
            }
            else
            {
                m_outputFrequency = atoi(pParams.find("OutputFrequency")->second.c_str());
            }


            m_session->MatchSolverInfo("Homogeneous","1D",m_isHomogeneous1D,false);
            
            if(m_isHomogeneous1D)
            {
                if (pParams.find("OutputPlane") == pParams.end())
                {
                    m_outputPlane = 0;
                }
                else
                {
                    m_outputPlane = atoi(pParams.find("OutputPlane")->second.c_str());
                }
            }

            ASSERTL0(pParams.find("Points") != pParams.end(),
                     "Missing parameter 'Points'.");
            m_historyPointStream.str(pParams.find("Points")->second);
            m_index = 0;
        }


        /**
         *
         */
        FilterHistoryPoints::~FilterHistoryPoints()
        {

        }


        /**
         *
         */
        void FilterHistoryPoints::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            ASSERTL0(!m_historyPointStream.fail(),
                     "No history points in stream.");

            m_index = 0;
            m_historyList.clear();

            // Read history points
            Array<OneD, NekDouble>  gloCoord(3,0.0);
            int dim = pFields[0]->GetGraph()->GetSpaceDimension();
            int i = 0;
            while (!m_historyPointStream.fail())
            {
                m_historyPointStream >> gloCoord[0] >> gloCoord[1] >> gloCoord[2];
                if(m_isHomogeneous1D) // overwrite with plane z
                {
                    NekDouble Z = (pFields[0]->GetHomogeneousBasis()->GetZ())[m_outputPlane];
                    if(fabs(gloCoord[2] - Z) > NekConstants::kVertexTheSameDouble)
                    {
                        cout <<"Reseting History point from " << gloCoord[2] << 
                            " to " << Z << endl;
                    }
                    gloCoord[2] = Z;
                }
                
                if (!m_historyPointStream.fail())
                {
                    SpatialDomains::VertexComponentSharedPtr vert
                        = MemoryManager<SpatialDomains::VertexComponent>
                        ::AllocateSharedPtr(dim, i, gloCoord[0],
                                            gloCoord[1], gloCoord[2]);

                    m_historyPoints.push_back(vert);
                    ++i;
                }
            }


            // Determine the unique process responsible for each history point
            // For points on a partition boundary, must select a single process
            LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
            int vRank = vComm->GetRank();
            Array<OneD, int>        idList(m_historyPoints.size());
            Array<OneD, int>        procList(m_historyPoints.size(), -1);
            for (i = 0; i < m_historyPoints.size(); ++i)
            {
                m_historyPoints[i]->GetCoords(  gloCoord[0],
                                                gloCoord[1],
                                                gloCoord[2]);

                idList[i] = pFields[0]->GetExpIndex(gloCoord);
                if (idList[i] != -1) 
                {
                    if(m_isHomogeneous1D)
                    {
                        int j;
                        Array<OneD, const unsigned int> IDs = pFields[0]->GetZIDs();
                        for(j = 0; j < IDs.num_elements(); ++j)
                        {
                            if(IDs[j] == m_outputPlane)
                            {
                                break;
                            }
                        }
                        
                        if(j != IDs.num_elements())
                        {
                            m_outputPlane = j;
                            procList[i] = vRank;
                        }
                    }
                    else
                    {
                        procList[i] = vRank;
                    }
                }
            }
            vComm->AllReduce(procList, LibUtilities::ReduceMax);

            // Determine the element in which each history point resides.
            // If point is not in mesh (on this process), id is -1.
            for (i = 0; i < m_historyPoints.size(); ++i)
            {
                // If point lies on partition boundary, only the proc with max
                // rank retains posession.
                if (procList[i] != vRank)
                {
                    idList[i] = -1;
                }
                if (idList[i] != -1)
                {
                    m_historyLocalPointMap[m_historyList.size()] = i;
                    m_historyList.push_back(
                        std::pair<SpatialDomains::VertexComponentSharedPtr, int>(
                            m_historyPoints[i], idList[i]));
                }
            }

            // Collate the element ID list across processes and check each
            // history point is allocated to a process
            vComm->AllReduce(idList, LibUtilities::ReduceMax);
            if (vComm->GetRank() == 0)
            {
                for (i = 0; i < m_historyPoints.size(); ++i)
                {
                    m_historyPoints[i]->GetCoords(  gloCoord[0],
                                                    gloCoord[1],
                                                    gloCoord[2]);
                    ASSERTL0(idList[i] != -1, "History point " +
                             boost::lexical_cast<std::string>(gloCoord[0]) + ", " +
                             boost::lexical_cast<std::string>(gloCoord[1]) + ", " +
                             boost::lexical_cast<std::string>(gloCoord[2]) +
                             " cannot be found in the mesh.");
                }

                // Open output stream
                m_outputStream.open(m_outputFile.c_str());
                m_outputStream << "# History data for variables (:";

                for (i = 0; i < pFields.num_elements(); ++i)
                {
                    m_outputStream << m_session->GetVariable(i) <<",";
                }

                if(m_isHomogeneous1D)
                {
                    m_outputStream << ") at points:";
                }
                else
                {
                    m_outputStream << ") at points:" << endl;
                }

                for (i = 0; i < m_historyPoints.size(); ++i)
                {
                    m_historyPoints[i]->GetCoords(  gloCoord[0],
                                                    gloCoord[1],
                                                    gloCoord[2]);

                    m_outputStream << "# \t" << i;
                    m_outputStream.width(8);
                    m_outputStream << gloCoord[0];
                    m_outputStream.width(8);
                    m_outputStream << gloCoord[1];
                    m_outputStream.width(8);
                    m_outputStream << gloCoord[2];
                    m_outputStream << endl;
                }

                if(m_isHomogeneous1D)
                {
                    m_outputStream << "(in Wavespace)" << endl;
                }
            }
            v_Update(pFields, time);
        }


        /**
         *
         */
        void FilterHistoryPoints::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            // Only output every m_outputFrequency.
            if ((m_index++) % m_outputFrequency)
            {
                return;
            }

            int j         = 0;
            int k         = 0;
            int numPoints = m_historyPoints.size();
            int numFields = pFields.num_elements();
            LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
            Array<OneD, NekDouble> data(numPoints*numFields, 0.0);
            Array<OneD, NekDouble> gloCoord(3, 0.0);
            std::list<std::pair<SpatialDomains::VertexComponentSharedPtr, int> >::iterator x;

            // Pull out data values field by field
            for (j = 0; j < numFields; ++j)
            {
                if(m_isHomogeneous1D)
                {
                    if(pFields[j]->GetPhysState() == false)
                    {
                        pFields[j]->GetPlane(m_outputPlane)->BwdTrans(pFields[j]->GetPlane(m_outputPlane)->GetCoeffs(),pFields[j]->GetPlane(m_outputPlane)->UpdatePhys());
                    }
                    pFields[j]->GetPlane(m_outputPlane)->PutPhysInToElmtExp();
                    for (k = 0, x = m_historyList.begin(); x != m_historyList.end(); ++x, ++k)
                    {
                        (*x).first->GetCoords(gloCoord[0], gloCoord[1], gloCoord[2]);
                        data[m_historyLocalPointMap[k]*numFields+j] = pFields[j]->GetPlane(m_outputPlane)->GetExp((*x).second)->PhysEvaluate(gloCoord);
                    }
                }
                else
                {
                    if(pFields[j]->GetPhysState() == false)
                    {
                        pFields[j]->BwdTrans(pFields[j]->GetCoeffs(),pFields[j]->UpdatePhys());
                    }
                    pFields[j]->PutPhysInToElmtExp();
                    for (k = 0, x = m_historyList.begin(); x != m_historyList.end(); ++x, ++k)
                    {
                        (*x).first->GetCoords(gloCoord[0], gloCoord[1], gloCoord[2]);
                        data[m_historyLocalPointMap[k]*numFields+j] = pFields[j]->GetExp((*x).second)->PhysEvaluate(gloCoord);
                    }
                }
            }

            // Exchange history data
            // This could be improved to reduce communication but works for now
            vComm->AllReduce(data, LibUtilities::ReduceSum);

            // Only the root process writes out history data
            if (vComm->GetRank() == 0)
            {

                // Write data values point by point
                for (k = 0; k < m_historyPoints.size(); ++k)
                {
                    m_outputStream.width(8);
                    m_outputStream << setprecision(6) << time;
                    for (int j = 0; j < numFields; ++j)
                    {
                        m_outputStream.width(25);
                        m_outputStream << setprecision(16) << data[k*numFields+j];
                    }
                    m_outputStream << endl;
                }
            }
        }


        /**
         *
         */
        void FilterHistoryPoints::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            if (pFields[0]->GetComm()->GetRank() == 0)
            {
                m_outputStream.close();
            }
        }


        /**
         *
         */
        bool FilterHistoryPoints::v_IsTimeDependent()
        {
            return true;
        }
    }
}
