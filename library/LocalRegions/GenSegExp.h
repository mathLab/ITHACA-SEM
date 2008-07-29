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

            
            void SetUpPhysNormals(const SpatialDomains::Geometry2DSharedPtr& Geom, const int edge, bool NegateNormals = false);

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
        };
        
        // type defines for use of SegExp in a boost vector
        typedef boost::shared_ptr<GenSegExp>      GenSegExpSharedPtr;
        typedef std::vector< GenSegExpSharedPtr > GenSegExpVector;
    } //end of namespace
} //end of namespace

#endif // GENSEGEXP_H

//
// $Log: GenSegExp.h,v $
