///////////////////////////////////////////////////////////////////////////////
//
// File Expansion3D.cpp
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
// Description: File for Expansion3D routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion3D.h>
#include <SpatialDomains/Geometry3D.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        Expansion3D::Expansion3D(){}
		
        DNekMatSharedPtr Expansion3D::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            DNekMatSharedPtr returnval;
            ASSERTL0(false,"This matrix type is not set up");
            
            return returnval;
        }

        void Expansion3D::SetFaceExp(const int face, Expansion2DSharedPtr &f)
        {
            int nFaces = GetNfaces();
            ASSERTL1(face < nFaces, "Face is out of range.");
            if (m_faceExp.size() < nFaces)
            {
                m_faceExp.resize(nFaces);
            }
            m_faceExp[face] = f;
        }
        
        Expansion2DSharedPtr Expansion3D::GetFaceExp(const int face)
        {
            return m_faceExp[face].lock();
        }
        
        void Expansion3D::v_AddFaceNormBoundaryInt(
            const int                            face,
            StdRegions::StdExpansion2DSharedPtr &FaceExp,
            const Array<OneD, const NekDouble>  &Fx,
            const Array<OneD, const NekDouble>  &Fy,
            const Array<OneD, const NekDouble>  &Fz,
                  Array<OneD,       NekDouble>  &outarray)
        {
            const Array<OneD, const Array<OneD, NekDouble> > normals
                                    = GetFaceNormal(face);
            int nquad_f = normals[0].num_elements();

            Vmath::Zero (nquad_f,FaceExp->UpdatePhys(),1);
            Vmath::Vmul (nquad_f,normals[0],1,Fx,1,FaceExp->UpdatePhys(),1);
            Vmath::Vvtvp(nquad_f,normals[1],1,Fy,1,FaceExp->GetPhys(),1,FaceExp->UpdatePhys(),1);
            Vmath::Vvtvp(nquad_f,normals[2],1,Fz,1,FaceExp->GetPhys(),1,FaceExp->UpdatePhys(),1);

            LocalRegions::Expansion2DSharedPtr locExp = 
                boost::dynamic_pointer_cast<
                    LocalRegions::Expansion2D>(FaceExp);
            
            if (locExp->GetRightAdjacentElementFace() != -1)
            {
                if (GetGeom3D()->GetFid(locExp->GetRightAdjacentElementFace()) ==
                    locExp->GetRightAdjacentElementExp()->GetGeom3D()->
                    GetFid(locExp->GetRightAdjacentElementFace()))
                {
                    Vmath::Neg(nquad_f,FaceExp->UpdatePhys(),1);
                }
            }
            
            AddFaceNormBoundaryInt(face, FaceExp, FaceExp->GetPhys(), outarray);
        }

        void Expansion3D::v_AddFaceNormBoundaryInt(
            const int                            face,
            StdRegions::StdExpansion2DSharedPtr &FaceExp,
            const Array<OneD, const NekDouble>  &Fn,
                  Array<OneD,       NekDouble>  &outarray)
        {
            int                         i;
            Array<OneD, unsigned int>   map;
            Array<OneD, int>            sign;
            StdRegions::FaceOrientation facedir = GetFaceOrient(face);
            
            GetFaceToElementMap(face,facedir,map,sign);
            int order_e = map.num_elements(); // Order of the element
            int n_coeffs = (FaceExp->GetCoeffs()).num_elements(); // Order of the trace

            if(n_coeffs!=order_e) // Going to orthogonal space
            {
                ASSERTL0(false, "Variable order not supported in 3D.");
            }
            else
            {
                FaceExp->IProductWRTBase(Fn,FaceExp->UpdateCoeffs());
                
                LocalRegions::Expansion2DSharedPtr locExp = 
                    boost::dynamic_pointer_cast<
                        LocalRegions::Expansion2D>(FaceExp);
                
                /*
                 * Coming into this routine, the velocity V will have been
                 * multiplied by the trace normals to give the input vector
                 * Vn. By convention, these normals are inwards facing for
                 * elements which have FaceExp as their right-adjacent face.
                 * This conditional statement therefore determines whether the
                 * normals must be negated, since the integral being performed
                 * here requires an outwards facing normal.
                 */ 
                if (locExp->GetRightAdjacentElementFace() != -1)
                {
                    if (GetGeom3D()->GetFid(locExp->GetRightAdjacentElementFace()) ==
                        locExp->GetRightAdjacentElementExp()->GetGeom3D()->
                        GetFid(locExp->GetRightAdjacentElementFace()))
                    {
                        Vmath::Neg(order_e,FaceExp->UpdateCoeffs(),1);
                    }
                }
            }
            
            for(i = 0; i < order_e; ++i)
            {
                outarray[map[i]] += sign[i]*FaceExp->GetCoeff(i);
            }
        }


#if 0 //needs m_faceMap to be defined and setupin Expansion3D similar to 2D case
        void Expansion3D::AddRobinMassMatrix(const int face, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
        {
            int i,j;
            int id1,id2;
            int order_e = m_faceExp[face]->GetNcoeffs();                    
            // Checks to see if this is a boundary interior
            // decomposed expansion - Routine not appropriate otherwise
            int nbndry = v_NumBndryCoeffs(); 

            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;
            
            LocalRegions::MatrixKey mkey(StdRegions::eMass,StdRegions::eSegment, *m_edgeExp[edge], primCoeffs);
            DNekScalMat &edgemat = *m_edgeExp[edge]->GetLocMatrix(mkey);

            v_GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);
            
            // Need to reset mapping to boundary rather than 
            // elemental mapping
            if(IsBoundaryMatrix == true) 
            {

                Array<OneD,unsigned int> bmap(nbndry);
                v_GetBoundaryMap(bmap);
                
                for(i = 0; i < order_e; ++i)
                {
                    for(j = 0; j < nbndry; ++j)
                    {
                        if(map[i] == bmap[j])
                        {
                            map[i] = j;
                            break;
                        }
                    }
                    ASSERTL1(j != nbndry,"Did not find number in map");
                }
            }

            for(i = 0; i < order_e; ++i)
            {
                id1 = map[i];
                for(j = 0; j < order_e; ++j)
                {
                    id2 = map[j];
                    (*inoutmat)(id1,id2) +=  edgemat(i,j)*sign[i]*sign[j];
                }
            }
        }

        void Expansion3D::v_AddRobinMassMatrix(const int faceid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
        {
            AddRobinMassMatrix(edgeid,primCoeffs,inoutmat);
        }
#endif

    } //end of namespace
} //end of namespace

/** 
 *    $Log: Expansion3D.cpp,v $
 *    Revision 1.2  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
