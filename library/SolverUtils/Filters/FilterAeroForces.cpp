///////////////////////////////////////////////////////////////////////////////
//
// File FilterAeroForces.cpp
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
// Description: Output values of aerodynamic forces during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <SolverUtils/Filters/FilterAeroForces.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterAeroForces::className =
                GetFilterFactory().RegisterCreatorFunction("AeroForces",
                                                    FilterAeroForces::create);

        /**
         *
         */
        FilterAeroForces::FilterAeroForces(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            // Load name of output file
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
                  && m_outputFile.substr(m_outputFile.length() - 4) == ".fce"))
            {
                m_outputFile += ".fce";
            }

            // Load frequency (in time-steps) of output
            if (pParams.find("OutputFrequency") == pParams.end())
            {
                m_outputFrequency = 1;
            }
            else
            {
                m_outputFrequency =
                    atoi(pParams.find("OutputFrequency")->second.c_str());
            }

            // Load time after which we need to calculate the forces
            if (pParams.find("StartTime") == pParams.end())
            {
                m_startTime = 0;
            }
            else
            {
                m_startTime = atof(pParams.find("StartTime")->second.c_str());
            }

            m_session->MatchSolverInfo("Homogeneous", "1D",
                                       m_isHomogeneous1D, false);

            // For 3DH1D, determine if we should calculate the forces in each
            //    plane, or just the average force
            if(m_isHomogeneous1D)
            {
                if (pParams.find("OutputAllPlanes") == pParams.end())
                {
                    m_outputAllPlanes = false;
                }
                else
                {
                    std::string sOption =
                            pParams.find("OutputAllPlanes")->second.c_str();
                    m_outputAllPlanes = ( boost::iequals(sOption,"true")) ||
                                        ( boost::iequals(sOption,"yes"));
                }
            }
            else
            {
                m_outputAllPlanes = false;
            }

            //specify the boundary to calculate the forces
            if (pParams.find("Boundary") == pParams.end())
            {
                ASSERTL0(false, "Missing parameter 'Boundary'.");
            }
            else
            {
                ASSERTL0(!(pParams.find("Boundary")->second.empty()),
                         "Missing parameter 'Boundary'.");
                m_BoundaryString = pParams.find("Boundary")->second;
            }

            //Get directions on which the forces will be calculated
            //    default is the coordinate axes orientation

            // Allocate m_directions
            m_directions = Array<OneD, Array<OneD, NekDouble> > (3);
            //Initialise directions to default values (ex, ey, ez)
            for (int i = 0; i < 3; ++i)
            {
                m_directions[i] = Array<OneD, NekDouble>(3, 0.0);
                m_directions[i][i] = 1.0;
            }
            std::stringstream       directionStream;
            std::string             directionString;
            //Override with input from xml file (if defined)
            for (int i = 0; i < 3; ++i)
            {
                std::stringstream tmp;
                tmp << i+1;
                std::string dir = "Direction" + tmp.str();
                if ( pParams.find(dir) != pParams.end() )
                {
                    ASSERTL0(!(pParams.find(dir)->second.empty()),
                             "Missing parameter '"+dir+"'.");
                    directionStream.str(pParams.find(dir)->second);
                    // Guarantee the stream is in its start position
                    //      before extracting
                    directionStream.clear();
                    // normalisation factor
                    NekDouble norm = 0.0;
                    for (int j = 0; j < 3; j++)
                    {
                        directionStream >> directionString;
                        if (!directionString.empty())
                        {
                            try
                            {
                                LibUtilities::Equation expression(
                                        pSession, directionString);
                                m_directions[i][j] = expression.Evaluate();
                                norm += m_directions[i][j]*m_directions[i][j];
                            }
                            catch (const std::runtime_error &)
                            {
                                ASSERTL0(false,
                                        "Error evaluating parameter expression"
                                        " '" + directionString + "'." );
                            }
                        }
                    }
                    //Normalise direction
                    for( int j = 0; j < 3; j++)
                    {
                        m_directions[i][j] /= sqrt(norm);
                    }
                }
            }
        }


        /**
         *
         */
        FilterAeroForces::~FilterAeroForces()
        {

        }


        /**
         *
         */
        void FilterAeroForces::v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            // Parse the boundary regions into a list.
            std::string::size_type FirstInd =
                                    m_BoundaryString.find_first_of('[') + 1;
            std::string::size_type LastInd =
                                    m_BoundaryString.find_last_of(']') - 1;

            ASSERTL0(FirstInd <= LastInd,
                    (std::string("Error reading boundary region definition:") +
                     m_BoundaryString).c_str());

            std::string IndString =
                    m_BoundaryString.substr(FirstInd, LastInd - FirstInd + 1);
            bool parseGood = ParseUtils::GenerateSeqVector(IndString.c_str(),
                                                       m_boundaryRegionsIdList);
            ASSERTL0(parseGood && !m_boundaryRegionsIdList.empty(),
                     (std::string("Unable to read boundary regions index "
                      "range for FilterAeroForces: ") + IndString).c_str());

            // determine what boundary regions need to be considered
            int cnt;
            unsigned int numBoundaryRegions =
                                pFields[0]->GetBndConditions().num_elements();
            m_boundaryRegionIsInList.insert(m_boundaryRegionIsInList.end(),
                                            numBoundaryRegions, 0);

            SpatialDomains::BoundaryConditions bcs(m_session,
                                                    pFields[0]->GetGraph());
            const SpatialDomains::BoundaryRegionCollection &bregions =
                                                    bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryRegionCollection::const_iterator it;

            for (cnt = 0, it = bregions.begin(); it != bregions.end();
                    ++it, cnt++)
            {
                if ( std::find(m_boundaryRegionsIdList.begin(),
                               m_boundaryRegionsIdList.end(), it->first) !=
                        m_boundaryRegionsIdList.end() )
                {
                    m_boundaryRegionIsInList[cnt] = 1;
                }
            }

            // Create map for element and edge/face of each boundary expansion
            if(m_isHomogeneous1D)
            {
                pFields[0]->GetPlane(0)->GetBoundaryToElmtMap
                                                (m_BCtoElmtID,m_BCtoTraceID);
            }
            else
            {
                pFields[0]->GetBoundaryToElmtMap(m_BCtoElmtID,m_BCtoTraceID);
            }

            // Define number of planes  to calculate the forces 
            //     in the Homogeneous direction ( if m_outputAllPlanes is false,
            //      consider only first plane in wave space)
            // If flow has no Homogeneous direction, use 1 to make code general
            if(m_isHomogeneous1D && m_outputAllPlanes)
            {
                m_nPlanes = pFields[0]->GetHomogeneousBasis()->
                                                    GetZ().num_elements();
            }
            else
            {
                m_nPlanes = 1;
            }
            
            // Create map for Planes ID for Homogeneous direction
            //    If flow has no Homogeneous direction, create trivial map  
            int j;
            m_planesID = Array<OneD, int> (m_nPlanes,-1);
            if(m_isHomogeneous1D)
            {
                Array<OneD, const unsigned int> IDs = pFields[0]->GetZIDs();
                //Loop through all planes
                for(cnt = 0; cnt < m_nPlanes; cnt++)
                {
                    for(j = 0; j < IDs.num_elements(); ++j)
                    {
                        //Look for current plane ID in this process
                        if(IDs[j] == cnt)
                        {
                            break;
                        }
                    }
                    // Assign ID to planesID
                    // If plane is not found in this process, leave it with -1
                    if(j != IDs.num_elements())
                    {
                        m_planesID[cnt] = j;
                    }
                }
            }
            else
            {
                m_planesID[0] = 0;
            }            

            LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

            // Write header
            int expdim = pFields[0]->GetGraph()->GetMeshDimension();
            if (vComm->GetRank() == 0)
            {
                // Open output stream
                m_outputStream.open(m_outputFile.c_str());
                m_outputStream << "# Forces acting on bodies" << endl;
                for( int i = 0; i < expdim; i++ )
                {
                    m_outputStream << "#" << " Direction" << i+1 << " = (";
                    m_outputStream.width(8);
                    m_outputStream << setprecision(4) << m_directions[i][0];
                    m_outputStream.width(8);
                    m_outputStream << setprecision(4) << m_directions[i][1];
                    m_outputStream.width(8);
                    m_outputStream << setprecision(4) << m_directions[i][2];
                    m_outputStream << ")" << endl;
                }
                m_outputStream << "# Boundary regions: " << IndString.c_str() << endl;
                m_outputStream << "#";
                m_outputStream.width(7);
                m_outputStream << "Time";
                for( int i = 1; i <= expdim; i++ )
                {
                    m_outputStream.width(8);
                    m_outputStream <<  "F" << i << "-press";
                    m_outputStream.width(9);
                    m_outputStream <<  "F" << i << "-visc";
                    m_outputStream.width(8);
                    m_outputStream <<  "F" << i << "-total";
                }
                if( m_outputAllPlanes )
                {
                    m_outputStream.width(10);
                    m_outputStream << "Plane";                   
                }
                if (m_session->DefinesSolverInfo("HomoStrip"))
                {
                    ASSERTL0(m_outputAllPlanes==false,
                            "Output forces on all planes not compatible with HomoStrips");
                    m_outputStream.width(10);
                    m_outputStream << "Strip";                    
                }
                
                m_outputStream << endl;
            }

            m_index = 0;
            v_Update(pFields, time);
        }


        /**
         *
         */
        void FilterAeroForces::v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            // Only output every m_outputFrequency.
            if ((m_index++) % m_outputFrequency  || (time < m_startTime))
            {
                return;
            }

            int i, j, k, n, cnt, elmtid, nq, offset, boundary, plane;
            
            // Get number of quadrature points and dimensions
            int physTot = pFields[0]->GetNpoints();
            int expdim = pFields[0]->GetGraph()->GetMeshDimension();
            int nVel = expdim;
            if( m_isHomogeneous1D )
            {
                nVel = nVel + 1;
            }
            
            StdRegions::StdExpansionSharedPtr elmt;
            
            // Fields used to calculate forces (a single plane for 3DH1D)
            Array<OneD, MultiRegions::ExpListSharedPtr>  
                                            fields( pFields.num_elements() );
            
            // Arrays of variables in the element
            Array<OneD, Array<OneD, NekDouble> >       velocity(expdim);
            Array<OneD, NekDouble>                     P(physTot);
            
            // Velocity gradient
            Array<OneD, Array<OneD, NekDouble> >       grad( expdim*expdim);
            
            // Values at the boundary
            Array<OneD, NekDouble>                     Pb; 
            Array<OneD, Array<OneD, NekDouble> >       gradb( expdim*expdim);

            // Communicators to exchange results
            LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
            LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
            LibUtilities::CommSharedPtr colComm = 
                                    m_session->DefinesSolverInfo("HomoStrip") ?
                                        vComm->GetColumnComm()->GetColumnComm():
                                        vComm->GetColumnComm();
            
            // Arrays with forces in each plane
            Array<OneD, Array<OneD, NekDouble> > Fpplane (expdim);
            Array<OneD, Array<OneD, NekDouble> > Fvplane (expdim);
            Array<OneD, Array<OneD, NekDouble> > Ftplane (expdim);
            for( i = 0; i < expdim; i++)
            {
                Fpplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
                Fvplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
                Ftplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
            }
            
            // Arrays with force
            Array<OneD, NekDouble>  Fp(expdim,0.0);
            Array<OneD, NekDouble>  Fv(expdim,0.0);
            Array<OneD, NekDouble>  Ft(expdim,0.0);
            
            // Forces per element length in a boundary
            Array<OneD, Array<OneD, NekDouble> >       fp( expdim );
            Array<OneD, Array<OneD, NekDouble> >       fv( expdim );

            // Get viscosity
            NekDouble rho = (m_session->DefinesParameter("rho"))
                    ? (m_session->GetParameter("rho"))
                    : 1;
            NekDouble mu = rho*m_session->GetParameter("Kinvis");
            
            // Perform BwdTrans: when we only want the mean force in a 3DH1D
            //     we work in wavespace, otherwise we use physical space            
            for(i = 0; i < pFields.num_elements(); ++i)
            {
                if (m_isHomogeneous1D && m_outputAllPlanes)
                {
                    pFields[i]->SetWaveSpace(false);
                }
                pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                                     pFields[i]->UpdatePhys());
                pFields[i]->SetPhysState(true);
            }
            
            // Define boundary expansions
            Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
            Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
            if(m_isHomogeneous1D)
            {
                BndConds = pFields[0]->GetPlane(0)->GetBndConditions();
                BndExp   = pFields[0]->GetPlane(0)->GetBndCondExpansions();                
            }
            else
            {
                BndConds = pFields[0]->GetBndConditions();
                BndExp   = pFields[0]->GetBndCondExpansions();                
            }
            
            // For Homogeneous, calculate force on each 2D plane
            // Otherwise, m_nPlanes = 1, and loop only runs once
            for(plane = 0; plane < m_nPlanes; plane++ )
            {
                // Check if plane is in this proc
                if( m_planesID[plane] != -1 )
                {
                    // For Homogeneous, consider the 2D expansion
                    //      on the current plane
                    if(m_isHomogeneous1D)
                    {
                        for(n = 0; n < pFields.num_elements(); n++)
                        {
                            fields[n] = pFields[n]->GetPlane(m_planesID[plane]);
                        }
                    }
                    else
                    {
                        for(n = 0; n < pFields.num_elements(); n++)
                        {
                            fields[n] = pFields[n];
                        }
                    } 
                    
                    //Loop all the Boundary Regions
                    for( cnt = n = 0; n < BndConds.num_elements(); n++)
                    {
                        if(m_boundaryRegionIsInList[n] == 1)
                        {
                            for (i=0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                            {
                                elmtid = m_BCtoElmtID[cnt];
                                elmt   = fields[0]->GetExp(elmtid);
                                nq     = elmt->GetTotPoints();
                                offset = fields[0]->GetPhys_Offset(elmtid);
                                
                                // Extract  fields on this element
                                for( j=0; j<expdim; j++)
                                {
                                    velocity[j] = fields[j]->GetPhys() + offset;
                                }
                                P = fields[nVel]->GetPhys() + offset;
                                
                                // Compute the velocity gradients
                                for (j=0; j<expdim; j++)
                                {
                                    for (k=0; k<expdim; k++)
                                    {
                                        grad[j*expdim+k] = 
                                                Array<OneD, NekDouble>(nq,0.0);
                                        elmt->PhysDeriv(k,velocity[j],
                                                grad[j*expdim+k]);
                                    }
                                }
                                
                                // identify boundary of element
                                boundary = m_BCtoTraceID[cnt];
                                
                                // Dimension specific part for obtaining values
                                //   at boundary and normal vector
                                Array<OneD, Array<OneD, NekDouble> > normals;
                                int nbc;
                                switch(expdim)
                                {
                                    case 2:
                                    {
                                        // Get expansion on boundary
                                        LocalRegions::Expansion1DSharedPtr bc;
                                        bc =  BndExp[n]->GetExp(i)->
                                               as<LocalRegions::Expansion1D> ();
                                        
                                        // Get number of points on the boundary
                                        nbc = bc->GetTotPoints();
                                        
                                        // Get normals
                                        normals = elmt->GetEdgeNormal(boundary);
                                        
                                        // Extract values at boundary
                                        Pb = Array<OneD, NekDouble>(nbc,0.0);
                                        elmt->GetEdgePhysVals(boundary,bc,P,Pb);
                                        for(int j = 0; j < expdim*expdim; ++j)
                                        {
                                            gradb[j] = Array<OneD, NekDouble>
                                                            (nbc,0.0);
                                            elmt->GetEdgePhysVals(boundary,
                                                           bc,grad[j],gradb[j]);
                                        }
                                    }
                                    break;
                                    
                                    case 3:
                                    {
                                        // Get expansion on boundary
                                        LocalRegions::Expansion2DSharedPtr bc;
                                        bc =  BndExp[n]->GetExp(i)->
                                               as<LocalRegions::Expansion2D> ();
 
                                        // Get number of points on the boundary
                                        nbc = bc->GetTotPoints();
                                        
                                        // Get normals
                                        normals = elmt->GetFaceNormal(boundary);
                                        
                                        // Extract values at boundary
                                        Pb = Array<OneD, NekDouble>(nbc,0.0);
                                        elmt->GetFacePhysVals(boundary,bc,P,Pb);
                                        for(int j = 0; j < expdim*expdim; ++j)
                                        {
                                            gradb[j] = Array<OneD, NekDouble>
                                                            (nbc,0.0);
                                            elmt->GetFacePhysVals(boundary,
                                                           bc,grad[j],gradb[j]);
                                        }                                        
                                    }
                                    break;
                                    
                                    default:
                                        ASSERTL0(false,
                                            "Expansion not supported by FilterForces");
                                    break;
                                }
                                
                                // Calculate forces per unit length
                                
                                // Pressure component: fp[j] = p*n[j]
                                for ( j = 0; j < expdim; j++)
                                {
                                    fp[j] = Array<OneD, NekDouble> (nbc,0.0);
                                    Vmath::Vmul (nbc, Pb, 1, 
                                                      normals[j], 1, 
                                                      fp[j], 1);
                                }
                                
                                // Viscous component:
                                //     fv[j] = -mu*{(grad[k,j]+grad[j,k]) *n[k]}
                                for ( j = 0; j < expdim; j++ )
                                {
                                    fv[j] = Array<OneD, NekDouble> (nbc,0.0);
                                    for ( k = 0; k < expdim; k++ )
                                    {
                                        Vmath::Vvtvp (nbc, gradb[k*expdim+j], 1,
                                                           normals[k], 1,
                                                           fv[j], 1,
                                                           fv[j], 1);
                                        Vmath::Vvtvp (nbc, gradb[j*expdim+k], 1,
                                                           normals[k], 1,
                                                           fv[j], 1, 
                                                           fv[j], 1);                                                
                                    }
                                    Vmath::Smul(nbc, -mu, fv[j], 1, fv[j], 1);
                                }

                                // Integrate to obtain force
                                for ( j = 0; j < expdim; j++)
                                {
                                    Fpplane[j][plane] += BndExp[n]->GetExp(i)->
                                                            Integral(fp[j]);
                                    Fvplane[j][plane] += BndExp[n]->GetExp(i)->
                                                            Integral(fv[j]);
                                }    
                            }
                        }
                        else
                        {
                            cnt += BndExp[n]->GetExpSize();
                        }
                    }
                }
            }
            
            // Combine contributions from different processes
            //    this is split between row and col comm because of
            //      homostrips case, which only keeps its own strip
            for( i = 0; i < expdim; i++)
            {
                rowComm->AllReduce(Fpplane[i], LibUtilities::ReduceSum);
                colComm->AllReduce(Fpplane[i], LibUtilities::ReduceSum);
                
                rowComm->AllReduce(Fvplane[i], LibUtilities::ReduceSum);
                colComm->AllReduce(Fvplane[i], LibUtilities::ReduceSum);
            }
            
            // Project results to new directions
            for(plane = 0; plane < m_nPlanes; plane++)
            {
                Array< OneD, NekDouble> tmpP(expdim, 0.0);
                Array< OneD, NekDouble> tmpV(expdim, 0.0);
                for( i = 0; i < expdim; i++)
                {   
                    for( j = 0; j < expdim; j++ )
                    {
                        tmpP[i] += Fpplane[j][plane]*m_directions[i][j];
                        tmpV[i] += Fvplane[j][plane]*m_directions[i][j];
                    }
                }
                // Copy result
                for( i = 0; i < expdim; i++)
                {
                    Fpplane[i][plane] = tmpP[i];
                    Fvplane[i][plane] = tmpV[i];
                }
            }
            
            // Sum viscous and pressure components
            for(plane = 0; plane < m_nPlanes; plane++)
            {
                for( i = 0; i < expdim; i++)
                {
                    Ftplane[i][plane] = Fpplane[i][plane] + Fvplane[i][plane];
                }
            }
            // Calculate forces including all planes
            for( i = 0; i < expdim; i++)
            {
                Fp[i] = Vmath::Vsum(m_nPlanes, Fpplane[i], 1) / m_nPlanes;
                Fv[i] = Vmath::Vsum(m_nPlanes, Fvplane[i], 1) / m_nPlanes;
                Ft[i] = Fp[i] + Fv[i];
            }            
            
            //Write Results
            if (vComm->GetRank() == 0)
            {
                // Write result in each plane
                if( m_outputAllPlanes)
                {
                    for( plane = 0; plane < m_nPlanes; plane++)
                    {
                        // Write time
                        m_outputStream.width(8);
                        m_outputStream << setprecision(6) << time;
                        // Write forces
                        for( i = 0; i < expdim; i++ )
                        {
                            m_outputStream.width(15);
                            m_outputStream << setprecision(8) 
                                           << Fpplane[i][plane];
                            m_outputStream.width(15);
                            m_outputStream << setprecision(8)
                                           << Fvplane[i][plane];
                            m_outputStream.width(15);
                            m_outputStream << setprecision(8)  
                                           << Ftplane[i][plane];
                        }
                        m_outputStream.width(10);
                        m_outputStream << plane;
                        m_outputStream << endl;
                    }                       
                }
                // Output average (or total) force
                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                for( i = 0; i < expdim; i++)
                {
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Fp[i];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Fv[i];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) << Ft[i];
                }
                if( m_outputAllPlanes)
                {
                    m_outputStream.width(10);
                    m_outputStream << "average";
                }
                
                if( m_session->DefinesSolverInfo("HomoStrip"))
                {
                    // The result we already wrote is for strip 0
                    m_outputStream.width(10);
                    m_outputStream << 0;
                    m_outputStream << endl;
                    
                    // Now get result from other strips
                    int nstrips;
                    m_session->LoadParameter("Strip_Z", nstrips);
                    for(i = 1; i<nstrips; i++)
                    {
                        vComm->GetColumnComm()->Recv(i, Fp);
                        vComm->GetColumnComm()->Recv(i, Fv);
                        vComm->GetColumnComm()->Recv(i, Ft);
                        
                        m_outputStream.width(8);
                        m_outputStream << setprecision(6) << time;
                        for( j = 0; j < expdim; j++)
                        {
                            m_outputStream.width(15);
                            m_outputStream << setprecision(8) << Fp[j];
                            m_outputStream.width(15);
                            m_outputStream << setprecision(8) << Fv[j];
                            m_outputStream.width(15);
                            m_outputStream << setprecision(8) << Ft[j];
                        }
                    m_outputStream.width(10);
                    m_outputStream << i;
                    m_outputStream << endl;                            
                    }                    
                }
                else
                {
                    m_outputStream << endl;
                }
            }
            else
            {
                // In the homostrips case, we need to send result to root
                if (m_session->DefinesSolverInfo("HomoStrip") &&
                        (rowComm->GetRank() == 0) )
                {
                        vComm->GetColumnComm()->Send(0, Fp);
                        vComm->GetColumnComm()->Send(0, Fv);
                        vComm->GetColumnComm()->Send(0, Ft);                    
                }
            }
            
            // Put results back to wavespace, if necessary
            if( m_isHomogeneous1D && m_outputAllPlanes )
            {
                for (i = 0; i < pFields.num_elements(); ++i)
                {
                    pFields[i]->SetWaveSpace(true);
                    pFields[i]->HomogeneousFwdTrans(pFields[i]->GetPhys(),
                                                    pFields[i]->UpdatePhys());
                }
            }
        }


        /**
         *
         */
        void FilterAeroForces::v_Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, 
            const NekDouble &time)
        {
            if (pFields[0]->GetComm()->GetRank() == 0)
            {
                m_outputStream.close();
            }
        }


        /**
         *
         */
        bool FilterAeroForces::v_IsTimeDependent()
        {
            return true;
        }
    }
}
