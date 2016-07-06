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
#include <MultiRegions/ExpList2D.h>     
#include <MultiRegions/ExpList3D.h>    
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <SolverUtils/Filters/FilterAeroForces.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterAeroForces::className =
        GetFilterFactory().RegisterCreatorFunction(
                "AeroForces", FilterAeroForces::create);

/**
 *
 */
FilterAeroForces::FilterAeroForces(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const ParamMap &pParams) :
    Filter(pSession)
{
    ParamMap::const_iterator it;

    // OutputFile
    it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    if (!(m_outputFile.length() >= 4
          && m_outputFile.substr(m_outputFile.length() - 4) == ".fce"))
    {
        m_outputFile += ".fce";
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_outputFrequency = floor(equ.Evaluate());
    }
    
    // Time after which we need to calculate the forces
    it = pParams.find("StartTime");
    if (it == pParams.end())
    {
        m_startTime = 0;
    }
    else
    {
        LibUtilities::Equation equ(m_session, it->second);
        m_startTime = equ.Evaluate();
    }

    // For 3DH1D, OutputAllPlanes or just average forces?
    m_session->MatchSolverInfo("Homogeneous", "1D",
                               m_isHomogeneous1D, false);
    if(m_isHomogeneous1D)
    {
        it = pParams.find("OutputAllPlanes");
        if (it == pParams.end())
        {
            m_outputAllPlanes = false;
        }
        else
        {
            std::string sOption =
                            it->second.c_str();
            m_outputAllPlanes = ( boost::iequals(sOption,"true")) ||
                                ( boost::iequals(sOption,"yes"));
        }
    }
    else
    {
        m_outputAllPlanes = false;
    }

    // Boundary (to calculate forces on)
    it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(),   "Missing parameter 'Boundary");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Boundary'.");
    m_BoundaryString = it->second;
    
    //
    // Directions (to project forces)
    //
    
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
        it = pParams.find(dir);
        if ( it != pParams.end() )
        {
            ASSERTL0(!(it->second.empty()),
                     "Missing parameter '"+dir+"'.");
            directionStream.str(it->second);
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
                    LibUtilities::Equation equ(m_session, directionString);
                    m_directions[i][j] = equ.Evaluate();
                    norm += m_directions[i][j]*m_directions[i][j];
                }
            }
            // Normalise direction
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
    // Load mapping
    m_mapping = GlobalMapping::Mapping::Load(m_session, pFields);

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
    if(m_isHomogeneous1D &&(m_outputAllPlanes || m_mapping->IsDefined()))
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

    m_lastTime = -1;
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
    // Calculate the forces
    CalculateForces(pFields, time);

    // Calculate forces including all planes
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    Array<OneD, NekDouble>  Fp(expdim,0.0);
    Array<OneD, NekDouble>  Fv(expdim,0.0);
    Array<OneD, NekDouble>  Ft(expdim,0.0);
    for( int i = 0; i < expdim; i++)
    {
        Fp[i] = Vmath::Vsum(m_nPlanes, m_Fpplane[i], 1) / m_nPlanes;
        Fv[i] = Vmath::Vsum(m_nPlanes, m_Fvplane[i], 1) / m_nPlanes;
        Ft[i] = Fp[i] + Fv[i];
    }   

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();                   

    //Write Results
    if (vComm->GetRank() == 0)
    {
        // Write result in each plane
        if( m_outputAllPlanes)
        {
            for( int plane = 0; plane < m_nPlanes; plane++)
            {
                // Write time
                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                // Write forces
                for( int i = 0; i < expdim; i++ )
                {
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8) 
                                   << m_Fpplane[i][plane];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8)
                                   << m_Fvplane[i][plane];
                    m_outputStream.width(15);
                    m_outputStream << setprecision(8)  
                                   << m_Ftplane[i][plane];
                }
                m_outputStream.width(10);
                m_outputStream << plane;
                m_outputStream << endl;
            }                       
        }
        // Output average (or total) force
        m_outputStream.width(8);
        m_outputStream << setprecision(6) << time;
        for( int i = 0; i < expdim; i++)
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
            for(int i = 1; i<nstrips; i++)
            {
                vComm->GetColumnComm()->Recv(i, Fp);
                vComm->GetColumnComm()->Recv(i, Fv);
                vComm->GetColumnComm()->Recv(i, Ft);

                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                for( int j = 0; j < expdim; j++)
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
                (vComm->GetRowComm()->GetRank() == 0) )
        {
                vComm->GetColumnComm()->Send(0, Fp);
                vComm->GetColumnComm()->Send(0, Fv);
                vComm->GetColumnComm()->Send(0, Ft);                    
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

/**
 *     This function outputs the force on all planes of the current
 *          process, in the format required by ForcingMovingBody
 */        
void FilterAeroForces::GetForces(
                    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
                    Array<OneD, NekDouble> &Aeroforces,
                    const NekDouble &time)
{
    // Calculate forces if the result we have is not up-to-date
    if(time > m_lastTime)
    {
        CalculateForces(pFields, time);
    }
    // Get information to write result
    Array<OneD, unsigned int> ZIDs = pFields[0]->GetZIDs();
    int local_planes = ZIDs.num_elements();
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();

    // Copy results to Aeroforces
    if (m_outputAllPlanes)
    {
        for(int plane = 0 ; plane < local_planes; plane++)
        {
            for(int dir =0; dir < expdim; dir++)
            {
                Aeroforces[plane + dir*local_planes] = 
                        m_Ftplane[dir][ZIDs[plane]];
            }
        }     
    }
    else
    {
        for(int plane = 0 ; plane < local_planes; plane++)
        {
            for(int dir =0; dir < expdim; dir++)
            {
                Aeroforces[plane + dir*local_planes] = 
                        m_Ftplane[dir][0];
            }
        }                  
    }
}

/**
 *     This function calculates the forces
 */        
void FilterAeroForces::CalculateForces(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    // Store time so we can check if result is up-to-date
    m_lastTime = time;

    // If a mapping is defined, call specific function
    //   Note: CalculateForcesMapping should work without a mapping,
    //         but since it is not very efficient the way it is now,
    //         it is only used when actually required
    if (m_mapping->IsDefined())
    {
        CalculateForcesMapping( pFields, time);
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
    m_Fpplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Fvplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Ftplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    for( i = 0; i < expdim; i++)
    {
        m_Fpplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Fvplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Ftplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
    }

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
                            m_Fpplane[j][plane] += BndExp[n]->GetExp(i)->
                                                    Integral(fp[j]);
                            m_Fvplane[j][plane] += BndExp[n]->GetExp(i)->
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
        rowComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);

        rowComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
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
                tmpP[i] += m_Fpplane[j][plane]*m_directions[i][j];
                tmpV[i] += m_Fvplane[j][plane]*m_directions[i][j];
            }
        }
        // Copy result
        for( i = 0; i < expdim; i++)
        {
            m_Fpplane[i][plane] = tmpP[i];
            m_Fvplane[i][plane] = tmpV[i];
        }
    }

    // Sum viscous and pressure components
    for(plane = 0; plane < m_nPlanes; plane++)
    {
        for( i = 0; i < expdim; i++)
        {
            m_Ftplane[i][plane] = m_Fpplane[i][plane] + m_Fvplane[i][plane];
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
 *     This function calculates the forces when we have a mapping
 *         defining a coordinate system transformation
 */        
void FilterAeroForces::CalculateForcesMapping(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble &time)
{
    int i, j, k, n, cnt, elmtid, offset, boundary, plane;
    // Get number of quadrature points and dimensions
    int physTot = pFields[0]->GetNpoints();
    int expdim = pFields[0]->GetGraph()->GetMeshDimension();
    int nVel = expdim;
    if( m_isHomogeneous1D )
    {
        nVel = nVel + 1;
    }

    StdRegions::StdExpansionSharedPtr elmt;

    // Pressure stress tensor 
    //    (global, in a plane, in element and boundary)
    Array<OneD, MultiRegions::ExpListSharedPtr>  P      ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  PPlane ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         PElmt  ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         PBnd  ( nVel*nVel);
    // Velocity gradient
    Array<OneD, MultiRegions::ExpListSharedPtr>  grad     ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  gradPlane( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         gradElmt ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         gradBnd  ( nVel*nVel);

    // Transformation to cartesian system
    Array<OneD, MultiRegions::ExpListSharedPtr>  C     ( nVel*nVel);
    Array<OneD, MultiRegions::ExpListSharedPtr>  CPlane( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         CElmt ( nVel*nVel);
    Array<OneD, Array<OneD, NekDouble> >         CBnd  ( nVel*nVel);            

    // Jacobian
    MultiRegions::ExpListSharedPtr  Jac;
    MultiRegions::ExpListSharedPtr  JacPlane;
    Array<OneD, NekDouble>          JacElmt;
    Array<OneD, NekDouble>          JacBnd;               

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
    LibUtilities::CommSharedPtr colComm = 
                            m_session->DefinesSolverInfo("HomoStrip") ?
                                vComm->GetColumnComm()->GetColumnComm():
                                vComm->GetColumnComm();

    // Arrays with forces in each plane
    m_Fpplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Fvplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    m_Ftplane = Array<OneD, Array<OneD, NekDouble> >  (expdim);
    for( i = 0; i < expdim; i++)
    {
        m_Fpplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Fvplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
        m_Ftplane[i] = Array<OneD, NekDouble>(m_nPlanes,0.0);
    }

    // Forces per element length in a boundary
    Array<OneD, Array<OneD, NekDouble> >       fp( nVel );
    Array<OneD, Array<OneD, NekDouble> >       fv( nVel );

    // Get viscosity
    NekDouble rho = (m_session->DefinesParameter("rho"))
            ? (m_session->GetParameter("rho"))
            : 1;
    NekDouble mu = rho*m_session->GetParameter("Kinvis");

    // Perform BwdTrans: for case with mapping, we can never work
    //                   in wavespace
    for(i = 0; i < pFields.num_elements(); ++i)
    {
        if (m_isHomogeneous1D)
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

    //
    // Calculate pressure stress tensor, velocity gradient
    //      and get informations about the mapping

    // Initialise variables
    switch (expdim)
    {
        case 2:
        {
            if (m_isHomogeneous1D)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
                Exp3DH1 = boost::dynamic_pointer_cast
                                <MultiRegions::ExpList3DHomogeneous1D>
                                                    (pFields[0]);
                for(i = 0; i < nVel*nVel; i++)
                {
                    grad[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);

                    P[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);

                    C[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);
                }
                Jac = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);                        
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = boost::dynamic_pointer_cast
                                <MultiRegions::ExpList2D>
                                                    (pFields[0]);
                for(i = 0; i < nVel*nVel; i++)
                {
                    grad[i] = MemoryManager<MultiRegions::ExpList2D>::
                                AllocateSharedPtr(*Exp2D);

                    P[i] = MemoryManager<MultiRegions::ExpList2D>::
                                AllocateSharedPtr(*Exp2D);

                    C[i] = MemoryManager<MultiRegions::ExpList2D>::
                                AllocateSharedPtr(*Exp2D);
                }
                Jac = MemoryManager<MultiRegions::ExpList2D>::
                                AllocateSharedPtr(*Exp2D);                         
            }
            break;
        }
        case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = boost::dynamic_pointer_cast
                            <MultiRegions::ExpList3D>
                                                (pFields[0]);
            for(i = 0; i < nVel*nVel; i++)
            {
                grad[i] = MemoryManager<MultiRegions::ExpList3D>::
                            AllocateSharedPtr(*Exp3D);

                P[i] = MemoryManager<MultiRegions::ExpList3D>::
                            AllocateSharedPtr(*Exp3D);

                C[i] = MemoryManager<MultiRegions::ExpList3D>::
                            AllocateSharedPtr(*Exp3D);
            }
            Jac = MemoryManager<MultiRegions::ExpList3D>::
                            AllocateSharedPtr(*Exp3D);

            break;
        }
        default:
            ASSERTL0(false,"Expansion dimension not supported by FilterAeroForces");
            break;
    }


    // Get g^ij
    Array<OneD, Array<OneD, NekDouble> >        tmp( nVel*nVel );
    m_mapping->GetInvMetricTensor(tmp);

    // Calculate P^ij = g^ij*p
    for (i = 0; i < nVel*nVel; i++)
    {
        Vmath::Vmul(physTot, tmp[i], 1,
                            pFields[nVel]->GetPhys(), 1,
                            P[i]->UpdatePhys(), 1);
    }

    // Calculate covariant derivatives of U = u^i_,k
    Array<OneD, Array<OneD, NekDouble> >        wk( nVel );
    for (i=0; i<nVel; i++)
    {
        wk[i] = Array<OneD, NekDouble>(physTot, 0.0);
        Vmath::Vcopy(physTot, pFields[i]->GetPhys(), 1,
                            wk[i], 1);                    
    }
    m_mapping->ApplyChristoffelContravar(wk, tmp);        
    for (i=0; i< nVel; i++)
    {
        for (k=0; k< nVel; k++)
        {
            pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                   wk[i], grad[i*nVel+k]->UpdatePhys());

            Vmath::Vadd(physTot,tmp[i*nVel+k],1,
                                grad[i*nVel+k]->UpdatePhys(), 1,
                                grad[i*nVel+k]->UpdatePhys(), 1);               
        }
    }   
    // Raise index to obtain Grad^ij = g^jk u^i_,k
    for (i=0; i< nVel; i++)
    {
        for (k=0; k< nVel; k++)
        {
            Vmath::Vcopy(physTot, grad[i*nVel+k]->GetPhys(), 1,
                                  wk[k], 1);
        }
        m_mapping->RaiseIndex(wk, wk);
        for (j=0; j<nVel; j++)
        {
            Vmath::Vcopy(physTot, wk[j], 1,
                                  grad[i*nVel+j]->UpdatePhys(), 1);
        }
    } 

    // Get Jacobian
    m_mapping->GetJacobian( Jac->UpdatePhys());

    // Get transformation to Cartesian system
    for (i=0; i< nVel; i++)
    {
        // Zero wk storage
        for (k=0; k< nVel; k++)
        {
            wk[k] = Array<OneD, NekDouble>(physTot, 0.0);
        }
        // Make wk[i] = 1
        wk[i] = Array<OneD, NekDouble>(physTot, 1.0);
        // Transform wk to Cartesian
        m_mapping->ContravarToCartesian(wk,wk);
        // Copy result to a column in C
        for (k=0; k< nVel; k++)
        {
            Vmath::Vcopy(physTot, wk[k], 1,
                                  C[k*nVel+i]->UpdatePhys(), 1);
        }                
    }           

    //
    // Calculate forces
    //

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
                for(n = 0; n < nVel*nVel; n++)
                {
                    PPlane[n]    = P[n]->GetPlane(m_planesID[plane]);
                    gradPlane[n] = grad[n]->GetPlane(m_planesID[plane]);
                    CPlane[n]    = C[n]->GetPlane(m_planesID[plane]);
                }
                JacPlane = Jac->GetPlane(m_planesID[plane]);
            }
            else
            {
                for(n = 0; n < nVel*nVel; n++)
                {
                    PPlane[n]    = P[n];
                    gradPlane[n] = grad[n];
                    CPlane[n] = C[n];
                }
                JacPlane = Jac;
            } 

            //Loop all the Boundary Regions
            for( cnt = n = 0; n < BndConds.num_elements(); n++)
            {
                if(m_boundaryRegionIsInList[n] == 1)
                {
                    for (i=0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                    {
                        elmtid = m_BCtoElmtID[cnt];
                        elmt   = PPlane[0]->GetExp(elmtid);
                        offset = PPlane[0]->GetPhys_Offset(elmtid);

                        // Extract  fields on this element
                        for( j=0; j<nVel*nVel; j++)
                        {
                            PElmt[j]    = PPlane[j]->GetPhys() 
                                        + offset;
                            gradElmt[j] = gradPlane[j]->GetPhys() 
                                        + offset;
                            CElmt[j]    = CPlane[j]->GetPhys() 
                                        + offset;
                        }
                        JacElmt = JacPlane->GetPhys() + offset;

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
                                for(int j = 0; j < nVel*nVel; ++j)
                                {
                                    gradBnd[j] = Array<OneD, NekDouble>
                                                    (nbc,0.0);
                                    elmt->GetEdgePhysVals(boundary,
                                                 bc,gradElmt[j],
                                                    gradBnd[j]);

                                    PBnd[j] = Array<OneD, NekDouble>
                                                    (nbc,0.0);
                                    elmt->GetEdgePhysVals(boundary,
                                                 bc,PElmt[j],
                                                    PBnd[j]);
                                    CBnd[j] = Array<OneD, NekDouble>
                                                    (nbc,0.0);
                                    elmt->GetEdgePhysVals(boundary,
                                                 bc,CElmt[j],
                                                    CBnd[j]);
                                }
                                JacBnd = Array<OneD, NekDouble>
                                                    (nbc,0.0);
                                elmt->GetEdgePhysVals(boundary,
                                                 bc,JacElmt,
                                                    JacBnd);
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
                                for(int j = 0; j < nVel*nVel; ++j)
                                {
                                    gradBnd[j] = Array<OneD, NekDouble>
                                                    (nbc,0.0);
                                    elmt->GetFacePhysVals(boundary,
                                                 bc,gradElmt[j],
                                                    gradBnd[j]);

                                    PBnd[j] = Array<OneD, NekDouble>
                                                    (nbc,0.0);
                                    elmt->GetFacePhysVals(boundary,
                                                 bc,PElmt[j],
                                                    PBnd[j]);          

                                    CBnd[j] = Array<OneD, NekDouble>
                                                    (nbc,0.0);
                                    elmt->GetFacePhysVals(boundary,
                                                 bc,CElmt[j],
                                                    CBnd[j]);
                                }
                                JacBnd = Array<OneD, NekDouble>
                                                (nbc,0.0);
                                elmt->GetFacePhysVals(boundary,
                                             bc,JacElmt,
                                                JacBnd);                                        
                            }
                            break;

                            default:
                                ASSERTL0(false,
                                    "Expansion not supported by FilterForces");
                            break;
                        }

                        // Calculate forces per unit length

                        // Pressure component: fp[j] = P[j,k]*n[k]
                        for ( j = 0; j < nVel; j++)
                        {
                            fp[j] = Array<OneD, NekDouble> (nbc,0.0);
                            // Normals only has expdim elements
                            for ( k = 0; k < expdim; k++)
                            {
                                Vmath::Vvtvp (nbc, PBnd[ j*nVel + k], 1,
                                                   normals[k], 1,
                                                   fp[j], 1, 
                                                   fp[j], 1);
                            }
                        }

                        // Viscous component:
                        //     fv[j] = -mu*{(grad[k,j]+grad[j,k]) *n[k]}
                        for ( j = 0; j < nVel; j++ )
                        {
                            fv[j] = Array<OneD, NekDouble> (nbc,0.0);
                            for ( k = 0; k < expdim; k++ )
                            {
                                Vmath::Vvtvp (nbc,gradBnd[k*nVel+j],1,
                                                   normals[k], 1,
                                                   fv[j], 1,
                                                   fv[j], 1);
                                Vmath::Vvtvp (nbc,gradBnd[j*nVel+k],1,
                                                   normals[k], 1,
                                                   fv[j], 1, 
                                                   fv[j], 1);                                                
                            }
                            Vmath::Smul(nbc, -mu, fv[j], 1, fv[j], 1);
                        }

                        // Multiply by Jacobian
                        for ( k = 0; k < nVel; k++ )
                        {
                            Vmath::Vmul(nbc, JacBnd, 1, fp[k], 1,
                                                        fp[k], 1);
                            Vmath::Vmul(nbc, JacBnd, 1, fv[k], 1,
                                                        fv[k], 1);
                        }                                

                        // Convert to cartesian system
                        for ( k = 0; k < nVel; k++ )
                        {
                            wk[k] = Array<OneD, NekDouble>(physTot,0.0);
                            for ( j = 0; j < nVel; j++ )
                            {
                                Vmath::Vvtvp(nbc, CBnd[k*nVel+j], 1,
                                                    fp[j], 1,
                                                    wk[k], 1,
                                                    wk[k], 1);
                            }
                        }
                        for ( k = 0; k < nVel; k++ )
                        {
                            Vmath::Vcopy(nbc, wk[k], 1, fp[k], 1);
                        }

                        for ( k = 0; k < nVel; k++ )
                        {
                            wk[k] = Array<OneD, NekDouble>(physTot,0.0);
                            for ( j = 0; j < nVel; j++ )
                            {
                                Vmath::Vvtvp(nbc, CBnd[k*nVel+j], 1,
                                                    fv[j], 1,
                                                    wk[k], 1,
                                                    wk[k], 1);
                            }
                        }
                        for ( k = 0; k < nVel; k++ )
                        {
                            Vmath::Vcopy(nbc, wk[k], 1, fv[k], 1);
                        }

                        // Integrate to obtain force
                        for ( j = 0; j < expdim; j++)
                        {
                            m_Fpplane[j][plane] += BndExp[n]->GetExp(i)->
                                                    Integral(fp[j]);
                            m_Fvplane[j][plane] += BndExp[n]->GetExp(i)->
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
        rowComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fpplane[i], LibUtilities::ReduceSum);

        rowComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
        colComm->AllReduce(m_Fvplane[i], LibUtilities::ReduceSum);
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
                tmpP[i] += m_Fpplane[j][plane]*m_directions[i][j];
                tmpV[i] += m_Fvplane[j][plane]*m_directions[i][j];
            }
        }
        // Copy result
        for( i = 0; i < expdim; i++)
        {
            m_Fpplane[i][plane] = tmpP[i];
            m_Fvplane[i][plane] = tmpV[i];
        }
    }

    // Sum viscous and pressure components
    for(plane = 0; plane < m_nPlanes; plane++)
    {
        for( i = 0; i < expdim; i++)
        {
            m_Ftplane[i][plane] = m_Fpplane[i][plane] + m_Fvplane[i][plane];
        }
    }

    // Put results back to wavespace, if necessary
    if( m_isHomogeneous1D)
    {
        for (i = 0; i < pFields.num_elements(); ++i)
        {
            pFields[i]->SetWaveSpace(true);
            pFields[i]->HomogeneousFwdTrans(pFields[i]->GetPhys(),
                                            pFields[i]->UpdatePhys());
        }
    }            
}

}
}
