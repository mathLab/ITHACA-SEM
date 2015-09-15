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
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <SolverUtils/Filters/FilterAeroForces.h>

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

    // OutputPlane
    m_session->MatchSolverInfo("Homogeneous", "1D",
                               m_isHomogeneous1D, false);
    if(m_isHomogeneous1D)
    {
        it = pParams.find("OutputPlane");
        if (it == pParams.end())
        {
            m_outputPlane = 0;
        }
        else
        {
            LibUtilities::Equation equ(m_session, it->second);
            m_outputPlane = floor(equ.Evaluate());
        }
    }

    // Boundary (to calculate forces on)
    it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(),   "Missing parameter 'Boundary");
    ASSERTL0(it->second.length() > 0, "Empty parameter 'Boundary'.");
    m_BoundaryString = it->second;
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

    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    if (vComm->GetRank() == 0)
    {
        // Open output stream
        m_outputStream.open(m_outputFile.c_str());
        m_outputStream << "#";
        m_outputStream.width(7);
        m_outputStream << "Time";
        m_outputStream.width(25);
        m_outputStream << "Fx (press)";
        m_outputStream.width(25);
        m_outputStream << "Fx (visc)";
        m_outputStream.width(25);
        m_outputStream << "Fx (tot)";
        m_outputStream.width(25);
        m_outputStream << "Fy (press)";
        m_outputStream.width(25);
        m_outputStream << "Fy (visc)";
        m_outputStream.width(25);
        m_outputStream << "Fy (tot)";
        m_outputStream.width(25);
        m_outputStream << "Fz (press)";
        m_outputStream.width(25);
        m_outputStream << "Fz (visc)";
        m_outputStream.width(25);
        m_outputStream << "Fz (tot)";
        m_outputStream << endl;
    }

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
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }

    int n, cnt, elmtid, nq, offset, nt, boundary;
    nt = pFields[0]->GetNpoints();
    int dim = pFields.num_elements()-1;

    StdRegions::StdExpansionSharedPtr elmt;
    Array<OneD, int> BoundarytoElmtID;
    Array<OneD, int> BoundarytoTraceID;
    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

    Array<OneD, const NekDouble> P(nt);
    Array<OneD, const NekDouble> U(nt);
    Array<OneD, const NekDouble> V(nt);
    Array<OneD, const NekDouble> W(nt);

    Array<OneD, Array<OneD, NekDouble> > gradU(dim);
    Array<OneD, Array<OneD, NekDouble> > gradV(dim);
    Array<OneD, Array<OneD, NekDouble> > gradW(dim);

    Array<OneD, Array<OneD, NekDouble> > fgradU(dim);
    Array<OneD, Array<OneD, NekDouble> > fgradV(dim);
    Array<OneD, Array<OneD, NekDouble> > fgradW(dim);

    Array<OneD, NekDouble> values;
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    NekDouble Fx,Fy,Fz,Fxp,Fxv,Fyp,Fyv,Fzp,Fzv;

    Fxp = 0.0; // x-component of the force due to pressure difference
    Fxv = 0.0; // x-component of the force due to viscous stress
    Fx = 0.0;  // x-component of the force (total) Fx = Fxp + Fxv (Drag)

    Fyp = 0.0; // y-component of the force due to pressure difference
    Fyv = 0.0; // y-component of the force due to viscous stress
    Fy = 0.0;  // y-component of the force (total) Fy = Fyp + Fyv (Lift)

    Fzp = 0.0; // z-component of the force due to pressure difference
    Fzv = 0.0; // z-component of the force due to viscous stress
    Fz = 0.0;  // z-component of the force (total) Fz = Fzp + Fzv (Side)

    NekDouble rho = (m_session->DefinesParameter("rho"))
            ? (m_session->GetParameter("rho"))
            : 1;
    NekDouble mu = rho*m_session->GetParameter("Kinvis");

    for(int i = 0; i < pFields.num_elements(); ++i)
    {
        pFields[i]->SetWaveSpace(false);
        pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                             pFields[i]->UpdatePhys());
        pFields[i]->SetPhysState(true);
    }

    // Homogeneous 1D case  Compute forces on all WALL boundaries
    // This only has to be done on the zero (mean) Fourier mode.
    if(m_isHomogeneous1D)
    {
        if(vComm->GetColumnComm()->GetRank() == 0)
        {
            pFields[0]->GetPlane(0)->GetBoundaryToElmtMap(
                                 BoundarytoElmtID,BoundarytoTraceID);
            BndExp = pFields[0]->GetPlane(0)->GetBndCondExpansions();
            StdRegions::StdExpansion1DSharedPtr bc;

            // loop over the types of boundary conditions
            for(cnt = n = 0; n < BndExp.num_elements(); ++n)
            {
                if(m_boundaryRegionIsInList[n] == 1)
                {
                    for(int i = 0; i <  BndExp[n]->GetExpSize();
                            ++i, cnt++)
                    {
                        // find element of this expansion.
                        elmtid = BoundarytoElmtID[cnt];
                        elmt   = pFields[0]->GetPlane(0)->GetExp(elmtid);
                        nq     = elmt->GetTotPoints();
                        offset = pFields[0]->GetPlane(0)->GetPhys_Offset(elmtid);

                        // Initialise local arrays for the velocity
                        // gradients size of total number of quadrature
                        // points for each element (hence local).
                        for(int j = 0; j < dim; ++j)
                        {
                            gradU[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradV[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradW[j] = Array<OneD, NekDouble>(nq,0.0);
                        }

                        // identify boundary of element
                        boundary = BoundarytoTraceID[cnt];

                        // Extract  fields
                        U = pFields[0]->GetPlane(0)->GetPhys() + offset;
                        V = pFields[1]->GetPlane(0)->GetPhys() + offset;
                        P = pFields[3]->GetPlane(0)->GetPhys() + offset;

                        // compute the gradients
                        elmt->PhysDeriv(U,gradU[0],gradU[1]);
                        elmt->PhysDeriv(V,gradV[0],gradV[1]);

                        // Get face 1D expansion from element expansion
                        bc =  BndExp[n]->GetExp(i)
                                ->as<LocalRegions::Expansion1D> ();

                        // number of points on the boundary
                        int nbc = bc->GetTotPoints();

                        // several vectors for computing the forces
                        Array<OneD, NekDouble> Pb(nbc,0.0);

                        for(int j = 0; j < dim; ++j)
                        {
                            fgradU[j] = Array<OneD, NekDouble>(nbc,0.0);
                            fgradV[j] = Array<OneD, NekDouble>(nbc,0.0);
                        }

                        Array<OneD, NekDouble>  drag_t(nbc,0.0);
                        Array<OneD, NekDouble>  lift_t(nbc,0.0);
                        Array<OneD, NekDouble>  drag_p(nbc,0.0);
                        Array<OneD, NekDouble>  lift_p(nbc,0.0);
                        Array<OneD, NekDouble>  temp(nbc,0.0);
                        Array<OneD, NekDouble>  temp2(nbc,0.0);

                        // identify boundary of element .
                        boundary = BoundarytoTraceID[cnt];

                        // extraction of the pressure and wss on the
                        // boundary of the element
                        elmt->GetEdgePhysVals(boundary,bc,P,Pb);

                        for(int j = 0; j < dim; ++j)
                        {
                            elmt->GetEdgePhysVals(boundary,bc,gradU[j],
                                                  fgradU[j]);
                            elmt->GetEdgePhysVals(boundary,bc,gradV[j],
                                                  fgradV[j]);
                        }

                        //normals of the element
                        const Array<OneD, Array<OneD, NekDouble> > &normals
                                    = elmt->GetEdgeNormal(boundary);

                        //
                        // Compute viscous tractive forces on wall from
                        //
                        //  t_i  = - T_ij * n_j  (minus sign for force
                        //                        exerted BY fluid ON wall),
                        //
                        // where
                        //
                        //  T_ij = viscous stress tensor (here in Cartesian
                        //         coords)
                        //                          dU_i    dU_j
                        //       = RHO * KINVIS * ( ----  + ---- ) .
                        //                          dx_j    dx_i

                        //a) DRAG TERMS
                        //-rho*kinvis*(2*du/dx*nx+(du/dy+dv/dx)*ny

                        Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,drag_t,1);
                        Vmath::Vmul(nbc,drag_t,1,normals[1],1,drag_t,1);

                        Vmath::Smul(nbc,2.0,fgradU[0],1,fgradU[0],1);
                        Vmath::Vmul(nbc,fgradU[0],1,normals[0],1,temp2,1);
                        Vmath::Smul(nbc,0.5,fgradU[0],1,fgradU[0],1);

                        Vmath::Vadd(nbc,temp2,1,drag_t,1,drag_t,1);
                        Vmath::Smul(nbc,-mu,drag_t,1,drag_t,1);

                        //zero temporary storage vector
                        Vmath::Zero(nbc,temp,0);
                        Vmath::Zero(nbc,temp2,0);


                        //b) LIFT TERMS
                        //-rho*kinvis*(2*dv/dy*nx+(du/dy+dv/dx)*nx

                        Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,lift_t,1);
                        Vmath::Vmul(nbc,lift_t,1,normals[0],1,lift_t,1);

                        Vmath::Smul(nbc,2.0,fgradV[1],1,fgradV[1],1);
                        Vmath::Vmul(nbc,fgradV[1],1,normals[1],1,temp2,1);
                        Vmath::Smul(nbc,-0.5,fgradV[1],1,fgradV[1],1);


                        Vmath::Vadd(nbc,temp2,1,lift_t,1,lift_t,1);
                        Vmath::Smul(nbc,-mu,lift_t,1,lift_t,1);

                        // Compute normal tractive forces on all WALL
                        // boundaries

                        Vmath::Vvtvp(nbc,Pb,1,normals[0],1,
                                        drag_p,1,drag_p, 1);
                        Vmath::Vvtvp(nbc,Pb,1,normals[1],1,
                                        lift_p,1,lift_p,1);

                        //integration over the boundary
                        Fxv += bc->Integral(drag_t);
                        Fyv += bc->Integral(lift_t);

                        Fxp += bc->Integral(drag_p);
                        Fyp += bc->Integral(lift_p);
                    }
                }
                else
                {
                        cnt += BndExp[n]->GetExpSize();
                }
            }
        }

        for(int i = 0; i < pFields.num_elements(); ++i)
        {
            pFields[i]->SetWaveSpace(true);
            pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                                 pFields[i]->UpdatePhys());
            pFields[i]->SetPhysState(false);
        }
    }
    //3D WALL case
    else if(dim==3 && !m_isHomogeneous1D)
    {
        pFields[0]->GetBoundaryToElmtMap(BoundarytoElmtID,
                                         BoundarytoTraceID);
        BndExp = pFields[0]->GetBndCondExpansions();
        LocalRegions::Expansion2DSharedPtr bc;

        // loop over the types of boundary conditions
        for(cnt = n = 0; n < BndExp.num_elements(); ++n)
        {
            if(m_boundaryRegionIsInList[n] == 1)
            {
                for(int i = 0; i <  BndExp[n]->GetExpSize(); ++i, cnt++)
                {
                    // find element of this expansion.
                    elmtid = BoundarytoElmtID[cnt];
                    elmt   = pFields[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = pFields[0]->GetPhys_Offset(elmtid);

                    // Initialise local arrays for the velocity
                    // gradients size of total number of quadrature
                    // points for each element (hence local).
                    for(int j = 0; j < dim; ++j)
                    {
                        gradU[j] = Array<OneD, NekDouble>(nq,0.0);
                        gradV[j] = Array<OneD, NekDouble>(nq,0.0);
                        gradW[j] = Array<OneD, NekDouble>(nq,0.0);
                    }

                    //identify boundary of element
                    boundary = BoundarytoTraceID[cnt];

                    //Extract  fields
                    U = pFields[0]->GetPhys() + offset;
                    V = pFields[1]->GetPhys() + offset;
                    W = pFields[2]->GetPhys() + offset;
                    P = pFields[3]->GetPhys() + offset;

                    //compute the gradients
                    elmt->PhysDeriv(U,gradU[0],gradU[1],gradU[2]);
                    elmt->PhysDeriv(V,gradV[0],gradV[1],gradV[2]);
                    elmt->PhysDeriv(W,gradW[0],gradW[1],gradW[2]);

                    // Get face 2D expansion from element expansion
                    bc =  BndExp[n]->GetExp(i)
                            ->as<LocalRegions::Expansion2D> ();

                    //number of points on the boundary
                    int nbc = bc->GetTotPoints();

                    //several vectors for computing the forces
                    Array<OneD, NekDouble> Pb(nbc,0.0);

                    for(int j = 0; j < dim; ++j)
                    {
                        fgradU[j] = Array<OneD, NekDouble>(nbc,0.0);
                        fgradV[j] = Array<OneD, NekDouble>(nbc,0.0);
                        fgradW[j] = Array<OneD, NekDouble>(nbc,0.0);

                    }

                    Array<OneD, NekDouble>  drag_t(nbc,0.0);
                    Array<OneD, NekDouble>  lift_t(nbc,0.0);
                    Array<OneD, NekDouble>  side_t(nbc,0.0);
                    Array<OneD, NekDouble>  drag_p(nbc,0.0);
                    Array<OneD, NekDouble>  lift_p(nbc,0.0);
                    Array<OneD, NekDouble>  side_p(nbc,0.0);
                    Array<OneD, NekDouble>  temp(nbc,0.0);
                    Array<OneD, NekDouble>  temp2(nbc,0.0);

                    // identify boundary of element .
                    boundary = BoundarytoTraceID[cnt];

                    // extraction of the pressure and wss on the
                    // boundary of the element
                    elmt->GetFacePhysVals(boundary,bc,P,Pb);

                    for(int j = 0; j < dim; ++j)
                    {
                        elmt->GetFacePhysVals(boundary,bc,gradU[j],
                                              fgradU[j]);
                        elmt->GetFacePhysVals(boundary,bc,gradV[j],
                                              fgradV[j]);
                        elmt->GetFacePhysVals(boundary,bc,gradW[j],
                                              fgradW[j]);
                    }

                    // normals of the element
                    const Array<OneD, Array<OneD, NekDouble> > &normals
                                        = elmt->GetFaceNormal(boundary);

                    //
                    // Compute viscous tractive forces on wall from
                    //
                    //  t_i  = - T_ij * n_j  (minus sign for force
                    //                        exerted BY fluid ON wall),
                    //
                    // where
                    //
                    //  T_ij = viscous stress tensor (here in Cartesian
                    //         coords)
                    //                          dU_i    dU_j
                    //       = RHO * KINVIS * ( ----  + ---- ) .
                    //                          dx_j    dx_i

                    //a) DRAG TERMS
                    //-rho*kinvis*
                    //    (2*du/dx*nx+(du/dy+dv/dx)*ny+(du/dz+dw/dx)*nz)
                    Vmath::Vadd(nbc,fgradU[2],1,fgradW[0],1,temp,1);
                    Vmath::Neg(nbc,temp,1);
                    Vmath::Vmul(nbc,temp,1,normals[2],1,temp,1);

                    Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,drag_t,1);
                    Vmath::Neg(nbc,drag_t,1);
                    Vmath::Vmul(nbc,drag_t,1,normals[1],1,drag_t,1);

                    Vmath::Smul(nbc,-2.0,fgradU[0],1,fgradU[0],1);
                    Vmath::Vmul(nbc,fgradU[0],1,normals[0],1,temp2,1);
                    Vmath::Smul(nbc,-0.5,fgradU[0],1,fgradU[0],1);

                    Vmath::Vadd(nbc,temp,1,temp2,1,temp,1);
                    Vmath::Vadd(nbc,temp,1,drag_t,1,drag_t,1);
                    Vmath::Smul(nbc,mu,drag_t,1,drag_t,1);

                    //zero temporary storage vector
                    Vmath::Zero(nbc,temp,0);
                    Vmath::Zero(nbc,temp2,0);


                    //b) LIFT TERMS
                    //-rho*kinvis*
                    //    (2*dv/dy*nx+(du/dy+dv/dx)*nx+(dv/dz+dw/dy)*nz)
                    Vmath::Vadd(nbc,fgradV[2],1,fgradW[1],1,temp,1);
                    Vmath::Neg(nbc,temp,1);
                    Vmath::Vmul(nbc,temp,1,normals[2],1,temp,1);

                    Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,lift_t,1);
                    Vmath::Neg(nbc,lift_t,1);
                    Vmath::Vmul(nbc,lift_t,1,normals[0],1,lift_t,1);

                    Vmath::Smul(nbc,-2.0,fgradV[1],1,fgradV[1],1);
                    Vmath::Vmul(nbc,fgradV[1],1,normals[1],1,temp2,1);
                    Vmath::Smul(nbc,-0.5,fgradV[1],1,fgradV[1],1);

                    Vmath::Vadd(nbc,temp,1,temp2,1,temp,1);
                    Vmath::Vadd(nbc,temp,1,lift_t,1,lift_t,1);
                    Vmath::Smul(nbc,mu,lift_t,1,lift_t,1);

                    //zero temporary storage vector
                    Vmath::Zero(nbc,temp,0);
                    Vmath::Zero(nbc,temp2,0);

                    //b) SIDE TERMS
                    //-rho*kinvis*
                    //    (2*dv/dy*nx+(du/dy+dv/dx)*nx+(dv/dz+dw/dy)*nz)
                    Vmath::Vadd(nbc,fgradV[2],1,fgradW[1],1,temp,1);
                    Vmath::Neg(nbc,temp,1);
                    Vmath::Vmul(nbc,temp,1,normals[1],1,temp,1);

                    Vmath::Vadd(nbc,fgradU[2],1,fgradW[0],1,side_t,1);
                    Vmath::Neg(nbc,side_t,1);
                    Vmath::Vmul(nbc,side_t,1,normals[0],1,side_t,1);

                    Vmath::Smul(nbc,-2.0,fgradW[2],1,fgradW[2],1);
                    Vmath::Vmul(nbc,fgradW[2],1,normals[2],1,temp2,1);
                    Vmath::Smul(nbc,-0.5,fgradW[2],1,fgradW[2],1);

                    Vmath::Vadd(nbc,temp,1,temp2,1,temp,1);
                    Vmath::Vadd(nbc,temp,1,side_t,1,side_t,1);
                    Vmath::Smul(nbc,mu,side_t,1,side_t,1);


                    // Compute normal tractive forces on all WALL
                    // boundaries
                    Vmath::Vvtvp(nbc,Pb,1,normals[0],1,
                                     drag_p,1,drag_p,1);
                    Vmath::Vvtvp(nbc,Pb,1,normals[1],1,
                                     lift_p,1,lift_p,1);
                    Vmath::Vvtvp(nbc,Pb,1,normals[2],1,
                                 side_p,1,side_p,1);

                    //integration over the boundary
                    Fxv += bc->Expansion::Integral(drag_t);
                    Fyv += bc->Expansion::Integral(lift_t);
                    Fzv += bc->Expansion::Integral(side_t);

                    Fxp += bc->Expansion::Integral(drag_p);
                    Fyp += bc->Expansion::Integral(lift_p);
                    Fzp += bc->Expansion::Integral(side_p);
                }
            }
            else
            {
                cnt += BndExp[n]->GetExpSize();
            }
        }
    }
    //2D WALL Condition
    else
    {
        pFields[0]->GetBoundaryToElmtMap(BoundarytoElmtID,
                                         BoundarytoTraceID);
        BndExp = pFields[0]->GetBndCondExpansions();
        StdRegions::StdExpansion1DSharedPtr bc;

        // loop over the types of boundary conditions
        for(cnt = n = 0; n < BndExp.num_elements(); ++n)
        {
            if(m_boundaryRegionIsInList[n] == 1)
            {
                for(int i = 0; i <  BndExp[n]->GetExpSize(); ++i, cnt++)
                {

                    elmtid = BoundarytoElmtID[cnt];
                    elmt   = pFields[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = pFields[0]->GetPhys_Offset(elmtid);

                    for(int j = 0; j < dim; ++j)
                    {
                        gradU[j] = Array<OneD, NekDouble>(nq,0.0);
                        gradV[j] = Array<OneD, NekDouble>(nq,0.0);
                    }

                    boundary = BoundarytoTraceID[cnt];

                    U = pFields[0]->GetPhys() + offset;
                    V = pFields[1]->GetPhys() + offset;
                    P = pFields[2]->GetPhys() + offset;

                    elmt->PhysDeriv(U,gradU[0],gradU[1]);
                    elmt->PhysDeriv(V,gradV[0],gradV[1]);

                    bc =  BndExp[n]->GetExp(i)
                            ->as<LocalRegions::Expansion1D> ();

                    int nbc = bc->GetTotPoints();
                    Array<OneD, NekDouble> Pb(nbc,0.0);

                    Array<OneD, NekDouble>  drag_t(nbc,0.0);
                    Array<OneD, NekDouble>  lift_t(nbc,0.0);
                    Array<OneD, NekDouble>  drag_p(nbc,0.0);
                    Array<OneD, NekDouble>  lift_p(nbc,0.0);
                    Array<OneD, NekDouble>  temp(nbc,0.0);

                    boundary = BoundarytoTraceID[cnt];

                    elmt->GetEdgePhysVals(boundary,bc,P,Pb);

                    for(int j = 0; j < dim; ++j)
                    {
                        fgradU[j] = Array<OneD, NekDouble>(nbc,0.0);
                        fgradV[j] = Array<OneD, NekDouble>(nbc,0.0);

                    }

                    for(int j = 0; j < dim; ++j)
                    {
                        elmt->GetEdgePhysVals(boundary,bc,gradU[j],
                                              fgradU[j]);
                        elmt->GetEdgePhysVals(boundary,bc,gradV[j],
                                              fgradV[j]);
                    }

                    const Array<OneD, Array<OneD, NekDouble> > &normals
                                        = elmt->GetEdgeNormal(boundary);

                    Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,drag_t,1);
                    Vmath::Neg(nbc,drag_t,1);
                    Vmath::Vmul(nbc,drag_t,1,normals[1],1,drag_t,1);

                    Vmath::Smul(nbc,-2.0,fgradU[0],1,fgradU[0],1);
                    Vmath::Vmul(nbc,fgradU[0],1,normals[0],1,temp,1);
                    Vmath::Vadd(nbc,temp,1,drag_t,1,drag_t,1);
                    Vmath::Smul(nbc,mu,drag_t,1,drag_t,1);

                    Vmath::Vadd(nbc,fgradU[1],1,fgradV[0],1,lift_t,1);
                    Vmath::Neg(nbc,lift_t,1);
                    Vmath::Vmul(nbc,lift_t,1,normals[0],1,lift_t,1);
                    Vmath::Smul(nbc,-2.0,fgradV[1],1,fgradV[1],1);
                    Vmath::Vmul(nbc,fgradV[1],1,normals[1],1,temp,1);
                    Vmath::Vadd(nbc,temp,1,lift_t,1,lift_t,1);
                    Vmath::Smul(nbc,mu,lift_t,1,lift_t,1);

                    Vmath::Vvtvp(nbc,Pb,1,normals[0],1,
                                     drag_p,1,drag_p,1);
                    Vmath::Vvtvp(nbc,Pb,1,normals[1],1,
                                     lift_p,1,lift_p,1);

                    Fxp += bc->Integral(drag_p);
                    Fyp += bc->Integral(lift_p);

                    Fxv += bc->Integral(drag_t);
                    Fyv += bc->Integral(lift_t);
                }
            }
            else
            {
                cnt += BndExp[n]->GetExpSize();
            }

        }

    }

    vComm->AllReduce(Fxp, LibUtilities::ReduceSum);
    vComm->AllReduce(Fxv, LibUtilities::ReduceSum);
    Fx = Fxp + Fxv;

    vComm->AllReduce(Fyp, LibUtilities::ReduceSum);
    vComm->AllReduce(Fyv, LibUtilities::ReduceSum);
    Fy = Fyp + Fyv;

    vComm->AllReduce(Fzp, LibUtilities::ReduceSum);
    vComm->AllReduce(Fzv, LibUtilities::ReduceSum);
    Fz = Fzp + Fzv;


    if (vComm->GetRank() == 0)
    {
        m_outputStream.width(8);
        m_outputStream << setprecision(6) << time;

        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fxp;
        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fxv;
        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fx;

        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fyp;
        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fyv;
        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fy;

        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fzp;
        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fzv;
        m_outputStream.width(25);
        m_outputStream << setprecision(8) << Fz;

        m_outputStream << endl;
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
