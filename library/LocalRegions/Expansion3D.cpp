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

namespace Nektar
{
    namespace LocalRegions 
    {
        
        DNekMatSharedPtr Expansion3D::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {

            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
                default:
                    ASSERTL0(false,"This matrix type cannot is not set up");
                    break;
            }
            
            return returnval;
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
