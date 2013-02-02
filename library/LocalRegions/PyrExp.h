///////////////////////////////////////////////////////////////////////////////
//
// File PyrExp.h
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
// Description: Header file for PyrExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef PYREXP_H
#define PYREXP_H

#include <StdRegions/StdPyrExp.h>
#include <SpatialDomains/PyrGeom.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar
{
    namespace LocalRegions
    {

        class PyrExp: virtual public StdRegions::StdPyrExp, 
                      virtual public Expansion3D
        {
        public:

            /** \brief Constructor using BasisKey class for quadrature points and order
                definition */
            LOCAL_REGIONS_EXPORT PyrExp(const LibUtilities::BasisKey &Ba,
                                        const LibUtilities::BasisKey &Bb,
                                        const LibUtilities::BasisKey &Bc,
                                        const SpatialDomains::PyrGeomSharedPtr &geom);

            LOCAL_REGIONS_EXPORT PyrExp(const PyrExp &T);

            LOCAL_REGIONS_EXPORT ~PyrExp();

        protected:
            //-------------------------------
            // Integration Methods
            //-------------------------------
            LOCAL_REGIONS_EXPORT virtual NekDouble v_Integral(
                const Array<OneD, const NekDouble>& inarray);            
            
            
            //----------------------------
            // Differentiation Methods
            //----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble>& inarray, 
                      Array<OneD,       NekDouble>& out_d0,
                      Array<OneD,       NekDouble>& out_d1,
                      Array<OneD,       NekDouble>& out_d2);


            //---------------------------------------
            // Transforms
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual void v_FwdTrans(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            
            
            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble>& inarray, 
                      Array<OneD,       NekDouble>& outarray);


            //---------------------------------------
            // Evaluation functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual void v_GetCoords(
                Array<OneD, NekDouble> &coords_0,
                Array<OneD, NekDouble> &coords_1 = NullNekDouble1DArray,
                Array<OneD, NekDouble> &coords_2 = NullNekDouble1DArray);
            
            LOCAL_REGIONS_EXPORT virtual void v_GetCoord(
                const Array<OneD, const NekDouble> &Lcoords, 
                      Array<OneD,       NekDouble> &coords);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble> &coord);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble>& coord,
                const Array<OneD, const NekDouble>& physvals);


            //---------------------------------------
            // Helper functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT void v_WriteToFile(
                std::ofstream &outfile, 
                OutputFormat   format, 
                const bool     dumpVar = true, 
                std::string    var = "v");
            LOCAL_REGIONS_EXPORT virtual const SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const;
            LOCAL_REGIONS_EXPORT virtual const SpatialDomains::GeometrySharedPtr v_GetGeom() const;
            LOCAL_REGIONS_EXPORT virtual const SpatialDomains::Geometry3DSharedPtr& v_GetGeom3D() const;
            LOCAL_REGIONS_EXPORT virtual int v_GetCoordim();

            
            //---------------------------------------
            // Matrix creation functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                const StdRegions::StdMatrixKey &mkey);
            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                const StdRegions::StdMatrixKey &mkey);
            LOCAL_REGIONS_EXPORT virtual DNekScalMatSharedPtr v_GetLocMatrix(
                const MatrixKey &mkey);
            LOCAL_REGIONS_EXPORT virtual DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(
                const MatrixKey &mkey);
            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr CreateMatrix(
                const MatrixKey &mkey);
            LOCAL_REGIONS_EXPORT DNekScalBlkMatSharedPtr CreateStaticCondMatrix(
                const MatrixKey &mkey);
            
        private:
            SpatialDomains::Geometry3DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;
        };

        // type defines for use of PyrExp in a boost vector
        typedef boost::shared_ptr<PyrExp> PyrExpSharedPtr;
        typedef std::vector<PyrExpSharedPtr> PyrExpVector;
        typedef std::vector<PyrExpSharedPtr>::iterator PyrExpVectorIter;
    } //end of namespace
} //end of namespace

#endif  //PYREXP_H
