////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNek.h
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
//  Description: Nektar file format converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_PREPROCESSING_MESHCONVERT_INPUTNEK
#define UTILITIES_PREPROCESSING_MESHCONVERT_INPUTNEK

#include "Module.h"

namespace Nektar
{
    namespace Utilities
    {
        enum NekCurve
        {
            eFile,
            eRecon
        };

        typedef HOTriangle<NodeSharedPtr> HOSurf;
        typedef boost::shared_ptr<HOSurf> HOSurfSharedPtr;

        /**
         * Hash class for high-order surfaces.
         */
        struct HOSurfHash : std::unary_function<HOSurfSharedPtr, std::size_t>
        {
            /** 
             * Calculate hash of a given high-order surface p by taking
             * successive hashes of the vertex IDs.
             */
            std::size_t operator()(HOSurfSharedPtr const& p) const
            {
                std::size_t seed = 0;
                std::vector<int> ids = p->vertId;
                std::sort(ids.begin(), ids.end());
                for (int i = 0; i < ids.size(); ++i)
                {
                    boost::hash_combine(seed, ids[i]);
                }
                return seed;
            }
        };

        bool operator==(HOSurfSharedPtr const &p1, HOSurfSharedPtr const &p2);
        
        typedef boost::unordered_set<HOSurfSharedPtr, HOSurfHash> HOSurfSet;

        /**
         * Converter class for Nektar session files.
         */
        class InputNek : public InputModule
        {
        public:
            InputNek(MeshSharedPtr p_m);
            virtual ~InputNek();
            virtual void Process();
            
            /// Creates an instance of this class.
            static ModuleSharedPtr create(MeshSharedPtr m) {
                return MemoryManager<InputNek>::AllocateSharedPtr(m);
            }
            /// %ModuleKey for class.
            static ModuleKey className;
            
        private:
            void LoadHOSurfaces();
            int  GetNnodes     (LibUtilities::ShapeType elType);

            /**
             * Maps a curve tag to a filename containing surface information.
             */
            std::map<string, pair<NekCurve, string> > curveTags;

            /**
             * Maps a curve tag to the high-order surface data for that tag.
             */
            std::map<string, HOSurfSet> hoData;

            /**
             * Maps ordering of hsf standard element to Nektar++ ordering.
             */
            std::map<int, int> hoMap;
        };
    }
}

#endif
