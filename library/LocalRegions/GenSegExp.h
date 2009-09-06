///////////////////////////////////////////////////////////////////////////////
//
// File GenSegExp.h
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
// Description: Header file for GenSegExp  (Generalised SegExp) routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef GENSEGEXP_H
#define GENSEGEXP_H

#include <SpatialDomains/Geometry2D.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {
    
        class GenSegExp: public SegExp
        {

        public:
            
            GenSegExp(const LibUtilities::BasisKey &Ba, 
                      const SpatialDomains::Geometry1DSharedPtr &geom);

            ///Copy Constructor
            GenSegExp(const GenSegExp &S);


            ///Destructor
            ~GenSegExp();

            
            /**  \brief Given a 2D Expansion, \a exp2d, set normals to
                 the same as the values along edge \a edge. By
                 definition the normal will be defined as outwards
                 with respect to the two-dimensional element if the
                 edge is orientated \e eForwards otherwise it will be
                 inwards facing.

                 \param exp2D is the 2D Element with respect to the
                 normal which will be defined

                 \param edge is the edge id where the normal should be set. 
            **/
            void SetUpPhysNormals(const StdRegions::StdExpansionSharedPtr &exp2d, const int edge);

            void SetPhysNormals(Array<OneD, const NekDouble> &normal)
            {
                m_physNormal = normal;
            }

            void SetPhysBiNormals(Array<OneD,const NekDouble> &binormal)
            {
                m_physBiNormal = binormal;
            }
            
            const Array<OneD, const NekDouble>& GetPhysNormals(void)
            {
                ASSERTL1(m_physNormal.num_elements(),"normals not set");
                return  m_physNormal;
            }

            inline const Array<OneD, const NekDouble>& GetPhysBiNormals()
            {
                ASSERTL1(m_physBiNormal.num_elements(),"binormals not set");
                return  m_physBiNormal;
            }

            void NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx, const Array<OneD, const NekDouble> &Fy, Array< OneD, NekDouble> &outarray, bool NegateNormal = false);

        protected:

        private:
            Array<OneD, NekDouble> m_physNormal;
            Array<OneD, NekDouble> m_physBiNormal;
            
            GenSegExp();

            virtual const Array<OneD, const NekDouble>& v_GetPhysNormals(void)
            {
                return GetPhysNormals(); 
            }


            virtual void v_SetPhysNormals(Array<OneD, const NekDouble> &normals)
            {
                SetPhysNormals(normals); 
            }

            virtual void v_SetUpPhysNormals(const StdRegions::StdExpansionSharedPtr &exp2d, const int edge)
            {
                SetUpPhysNormals(exp2d,edge);
            }

            void v_NormVectorIProductWRTBase(const Array<OneD, const NekDouble> &Fx, const Array<OneD, const NekDouble> &Fy, Array< OneD, NekDouble> &outarray,
                                             bool NegateNorm = false)
            {
                NormVectorIProductWRTBase(Fx,Fy,outarray, NegateNorm);
            }
        };
        
        // type defines for use of SegExp in a boost vector
        typedef boost::shared_ptr<GenSegExp>      GenSegExpSharedPtr;
        typedef std::vector< GenSegExpSharedPtr > GenSegExpVector;
    } //end of namespace
} //end of namespace

#endif // GENSEGEXP_H

//
// $Log: GenSegExp.h,v $
// Revision 1.2  2008/09/09 15:05:09  sherwin
// Updates related to cuved geometries. Normals have been removed from m_metricinfo and replaced with a direct evaluation call. Interp methods have been moved to LibUtilities
//
// Revision 1.1  2008/07/29 22:24:49  sherwin
// Generalised Segment expansion which include a normal and binormal at the physical quadrature points
//
