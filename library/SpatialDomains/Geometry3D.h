////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/Geometry3D.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H


#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/Geometry.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry3D: public Geometry
        {
        public:
            Geometry3D();
            Geometry3D(const int coordim);
            ~Geometry3D();

            void AddElmtConnected(int gvo_id, int locid)
            {
                v_AddElmtConnected(gvo_id,locid);
            }

            int  NumElmtConnected() const
            {
                return v_NumElmtConnected();
            }

            bool IsElmtConnected(int gvo_id, int locid) const
            {
                return v_IsElmtConnected(gvo_id,locid);
            }
            
            inline int GetEid() const
            {
                return v_GetEid();
            }

            inline int GetCoordDim() const 
            {
                return v_GetCoordDim();
            }

            inline const LibUtilities::BasisSharedPtr GetBasis(const int i, const int j)
            {
                return v_GetBasis(i,j);
            }

            inline Array<OneD,NekDouble> &UpdatePhys(const int i)
            {
                return v_UpdatePhys(i);
            }

            NekDouble GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
            {
                return v_GetCoord(i,Lcoord);
            }

        protected:

        private:
            virtual void v_AddElmtConnected(int gvo_id, int locid)
            {  
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");            
            }
            
            virtual int v_NumElmtConnected() const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }
            
            virtual bool v_IsElmtConnected(int gvo_id, int locid) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return false;
            }
            
            virtual int v_GetEid() const 
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual int v_GetCoordDim() const 
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0;
            }

            virtual const LibUtilities::BasisSharedPtr v_GetBasis(const int i, const int j)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                LibUtilities::BasisSharedPtr returnval;        
                return returnval;    
            }

            virtual Array<OneD,NekDouble> &v_UpdatePhys(const int i)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return NullNekDouble1DArray;
            }

            virtual NekDouble v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for shape type geometries");
                return 0.0;
            }
        };
        // shorthand for boost pointer
        typedef boost::shared_ptr<Geometry3D> Geometry3DSharedPtr;
        typedef std::vector< Geometry3DSharedPtr > Geometry3DVector;
        typedef std::vector< Geometry3DSharedPtr >::iterator Geometry3DVectorIter;

    }; //end of namespace
}; //end of namespace


#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY3D_H

//
// $Log: Geometry3D.h,v $
// Revision 1.3  2008/04/02 22:19:03  pvos
// Update for 2D local to global mapping
//
// Revision 1.2  2008/01/31 11:02:55  ehan
// Added constructor and destructor.
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.14  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.13  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.12  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.11  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

