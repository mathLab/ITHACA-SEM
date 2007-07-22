///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalMap1D.h
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
// Description: Local to Global mapping routines in 1D, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP1D_H
#define NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP1D_H

#include <MultiRegions/LocalToGlobalMap.h>
#include <SpatialDomains/MeshGraph1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
    
    class LocalToGlobalMap1D: 
            public LocalToGlobalMap
    {
        public:
            LocalToGlobalMap1D(){};
        LocalToGlobalMap1D(const int loclen, 
                               const StdRegions::StdExpansionVector &locexp, 
                               const SpatialDomains::Composite &cmps);

            virtual ~LocalToGlobalMap1D();

            void ResetMapping(const int NumDirichlet, 
                              SpatialDomains::BoundaryConditions &bcs);
        
        protected:
        
        private:
    };
    
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP1D_H


/** $Log: LocalToGlobalMap1D.h,v $
/** Revision 1.8  2007/07/20 02:04:13  bnelson
/** Replaced boost::shared_ptr with Nektar::ptr
/**
/** Revision 1.7  2007/07/10 08:54:30  pvos
/** Updated ContField1D constructor
/**
/** Revision 1.6  2007/07/06 18:39:34  pvos
/** ContField1D constructor updates
/**
/** Revision 1.5  2007/05/28 16:15:00  sherwin
/** Updated files in MultiRegions to make 1D demos work
/**
/** Revision 1.4  2007/05/27 16:09:43  bnelson
/** Update to new Array type.
/**
/** Revision 1.3  2007/04/26 15:00:16  sherwin
/** SJS compiling working version using SHaredArrays
/**
/** Revision 1.2  2007/03/20 16:58:42  sherwin
/** Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
/**
/** Revision 1.1  2006/07/02 17:16:17  sherwin
/**
/** Modifications to make MultiRegions work for a connected domain in 2D (Tris)
/**
/** Revision 1.3  2006/06/05 00:14:33  bnelson
/** Fixed a compiler error (couldn't find boost::shared_ptr<) and a couple of formatting updates for the standard.
/** */

