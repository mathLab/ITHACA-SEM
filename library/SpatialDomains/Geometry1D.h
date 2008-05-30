////////////////////////////////////////////////////////////////////////////////
//
//  File: Geometry1D.h
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
#ifndef NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H
#define NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H

#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/Geometry.h>
#include <SpatialDomains/MeshComponents.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class Geometry1D: public Geometry
        {
        public:
            Geometry1D();
            Geometry1D(const int coordim);
            ~Geometry1D();

            // Wrappers around virtual Functions (for the EdgeComponent class)
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

            inline const LibUtilities::BasisSharedPtr GetBasis(const int i, const int j)
            {
                return v_GetBasis(i,j);
            }

            inline const StdRegions::StdExpansion1DSharedPtr &GetXmap(const int i)
            {
                return v_GetXmap(i);
            }

            inline Array<OneD,NekDouble> &UpdatePhys(const int i)
            {
                return v_UpdatePhys(i);
            }

            inline VertexComponentSharedPtr GetVertex(const int i) const
            {
                return v_GetVertex(i);
            }

            NekDouble GetCoord(const int i, const Array<OneD, const NekDouble> &Lcoord)
            {
                return v_GetCoord(i,Lcoord);
            }

            inline int GetVid(int i) const
            {
                return v_GetVid(i);
            }

            void FillGeom()
            {
                v_FillGeom();
            }
            
            StdRegions::ExpansionType DetExpansionType() const
            {
                return v_DetExpansionType();
            }
            
            void GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
            {
                v_GetLocCoords(coords,Lcoords);
            }            
            
            inline void SetOwnData()
            {
                v_SetOwnData();
            }
            
            void WriteToFile(std::ofstream &outfile, const int dumpVar)
            {
                v_WriteToFile(outfile,dumpVar);
            }
            
        protected:

        private:
            virtual void v_AddElmtConnected(int gvo_id, int locid)
            {  
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");            
            }
            
            virtual int v_NumElmtConnected() const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                return 0;
            }
            
            virtual bool v_IsElmtConnected(int gvo_id, int locid) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                return false;
            }
            
            virtual int v_GetEid() const 
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                return 0;
            }

            virtual const LibUtilities::BasisSharedPtr v_GetBasis(const int i, const int j)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                LibUtilities::BasisSharedPtr returnval;        
                return returnval;    
            }

            virtual const StdRegions::StdExpansion1DSharedPtr &v_GetXmap(const int i)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                static StdRegions::StdExpansion1DSharedPtr returnval;        
                return returnval; 
            }

            virtual Array<OneD,NekDouble> &v_UpdatePhys(const int i)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                return NullNekDouble1DArray;
            }

            virtual VertexComponentSharedPtr v_GetVertex(const int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                VertexComponentSharedPtr returnval;        
                return returnval;
            }

            virtual NekDouble v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                return 0.0;
            }

            virtual int v_GetVid(int i) const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                return 0;
            }
            
            virtual void v_FillGeom()
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
            }
            
            virtual StdRegions::ExpansionType v_DetExpansionType() const
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
                return StdRegions::eNoExpansionType;
            }
            
            virtual void v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
            }            
            
            virtual  void v_SetOwnData()
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
            }
            
            virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
            {
                NEKERROR(ErrorUtil::efatal,
                         "This function is only valid for expansion type geometries");
            }
        };
        // shorthand for boost pointer
        typedef boost::shared_ptr<Geometry1D> Geometry1DSharedPtr;
        typedef std::vector< Geometry1DSharedPtr > Geometry1DVector;
        typedef std::vector< Geometry1DSharedPtr >::iterator Geometry1DVectorIter;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_GEOMETRY1D_H

//
// $Log: Geometry1D.h,v $
// Revision 1.6  2008/04/06 06:00:37  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.5  2008/04/02 22:19:03  pvos
// Update for 2D local to global mapping
//
// Revision 1.4  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.3  2007/03/20 09:17:40  kirby
//
// GeomFactors now added; metricinfo used instead of minfo; styles updated
//
// Revision 1.2  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.1  2006/05/04 18:58:59  kirby
// *** empty log message ***
//
// Revision 1.19  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.18  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.17  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.16  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
