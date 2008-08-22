///////////////////////////////////////////////////////////////////////////////
//
// File GenExpList1D.h
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
// Description: Generalised Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef GENEXPLIST1D_H
#define GENEXPLIST1D_H

#include <LocalRegions/GenSegExp.h>
#include <MultiRegions/ExpList1D.h>

namespace Nektar
{
    namespace MultiRegions
    {     

        class GenExpList1D: 
        public ExpList1D
        {
        public:

            /**
             * \brief The default constructor.  
             */  
            GenExpList1D();
            
            // constructor for trace space in connection with DisContField2D.cpp
            GenExpList1D(const Array<OneD,const MultiRegions::ExpList1DSharedPtr> &bndConstraint,  
                         const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond, 
                         const StdRegions::StdExpansionVector &locexp, 
                         SpatialDomains::MeshGraph2D &graph2D,
                         const map<int,int> &periodicEdges);
            
            /**
             * \brief The copy constructor.
             */  
            GenExpList1D(const GenExpList1D &In);   
            
            
            /**
             * \brief The default destructor.
             */  
            ~GenExpList1D();
            
            void Upwind(Array<OneD,Array<OneD, const NekDouble> > &Vec, 
                        Array<OneD, const NekDouble> &Fwd, 
                        Array<OneD, const NekDouble> &Bwd, 
                        Array<OneD, NekDouble> &Upwind);

	    void GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals); 

        protected:

        private:

        };

        typedef boost::shared_ptr<GenExpList1D>     GenExpList1DSharedPtr;
        typedef std::vector<GenExpList1DSharedPtr>   GenExpList1DVector;
        typedef std::vector<GenExpList1DSharedPtr>::iterator GenExpList1DVectorIter;

    } //end of namespace
} //end of namespace

#endif//GENEXPLIST1D_H

/**
 * $Log: GenExpList1D.h,v $
 * Revision 1.2  2008/08/14 22:15:51  sherwin
 * Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
 *
 * Revision 1.1  2008/07/29 22:26:35  sherwin
 * Generalised 1D Segment list which includes a normal direction at physical points
 *
 **/
