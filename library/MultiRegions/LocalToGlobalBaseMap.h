///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalBaseMap.h
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
// Description: Local to Global base mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_LOCALTOGLOBALBASEMAP_H
#define MULTIREGIONS_LOCALTOGLOBALBASEMAP_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        class LocalToGlobalBaseMap
        {
        public:
            LocalToGlobalBaseMap(); 
            
            ~LocalToGlobalBaseMap(); 

            int GetLocalToGlobalBndMapSize() const
            {
                return m_localToGlobalBndMap.num_elements();
            }
                
            int GetLocalToGlobalBndMap(const int i) const
            {
                ASSERTL1((i >= 0) && (i < m_localToGlobalBndMap.num_elements()),
                         "index is out of range ");
                return m_localToGlobalBndMap[i];
            }
            
            inline Array<OneD,const int>  GetLocalToGlobalBndMap(void)
            {
                return m_localToGlobalBndMap;
            }

            int GetLocalToGlobalBndSign(const int i) const
            {
                ASSERTL1((i >= 0) && (i < m_localToGlobalBndSign.num_elements()),
                         "index is out of range ");
                return m_localToGlobalBndSign[i];
            }
             
            inline Array<OneD,const int>  GetLocalToGlobalBndSign(void)
            {
                return m_localToGlobalBndSign;
            }

            inline int GetNumDirichletBndCoeffs() const
            {
                return m_numDirichletBndCoeffs;
            }


            inline int GetNumGlobalBndCoeffs() const
            {
                return m_numGlobalBndCoeffs;
            }
        protected:
            int m_numGlobalBndCoeffs;     // Total number of global dofs
            int m_numDirichletBndCoeffs;  // Number of Dirichlet Boundary Coefficient
            Array<OneD,int> m_localToGlobalBndMap;  //< integer map of boundary/trace space 
            Array<OneD,int> m_localToGlobalBndSign; //< integer sign of boundary/trace space 
        private:
        };
        typedef boost::shared_ptr<LocalToGlobalBaseMap>  LocalToGlobalBaseMapSharedPtr;
        
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP_H


/** $Log: LocalToGlobalBaseMap.h,v $
 */

