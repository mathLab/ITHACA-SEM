///////////////////////////////////////////////////////////////////////////////
//
// File Expansion.cpp
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
// Description: File for Expansion routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion.h>
#include <LocalRegions/MatrixKey.h>

#include <SpatialDomains/MeshComponents.h>
namespace Nektar
{
    namespace LocalRegions 
    {
        Expansion::Expansion(SpatialDomains::GeometrySharedPtr pGeom) :
                    m_geom(pGeom),
                    m_metricinfo(m_geom->GetGeomFactors(m_base))
        {
        }
        
        Expansion::Expansion(const Expansion &pSrc) :
                m_geom(pSrc.m_geom),
                m_metricinfo(pSrc.m_metricinfo)
        {

        }

        Expansion::~Expansion()
        {
        }

        DNekScalMatSharedPtr Expansion::GetLocMatrix(const LocalRegions::MatrixKey &mkey)
        {
            return v_GetLocMatrix(mkey);
        }

        DNekScalMatSharedPtr Expansion::GetLocMatrix(const StdRegions::MatrixType mtype,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeffs)
        {
            MatrixKey mkey(mtype, DetExpansionType(), *this, factors, varcoeffs);
            return GetLocMatrix(mkey);
        }

        SpatialDomains::GeometrySharedPtr Expansion::GetGeom() const
        {
            return m_geom;
        }

        const SpatialDomains::GeomFactorsSharedPtr& Expansion::v_GetMetricInfo() const
        {
            return m_metricinfo;
        }


        DNekScalMatSharedPtr Expansion::v_GetLocMatrix(const LocalRegions::MatrixKey &mkey)
        {
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            return NullDNekScalMatSharedPtr;
        }

        void Expansion::v_MultiplyByQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
                                                 Array<OneD, NekDouble> &outarray)
        {
            const int nqtot = GetTotPoints();

            if (m_metrics.count(MetricQuadrature) == 0)
            {
                ComputeQuadratureMetric();
            }

            Vmath::Vmul(nqtot, m_metrics[MetricQuadrature], 1, inarray, 1, outarray, 1);
        }

        void Expansion::ComputeLaplacianMetric()
        {
            if (m_metrics.count(MetricQuadrature) == 0)
            {
                ComputeQuadratureMetric();
            }

            const SpatialDomains::GeomType type = m_metricinfo->GetGtype();
            const unsigned int nqtot = GetTotPoints();
            const unsigned int dim = GetShapeDimension();
            const unsigned int pts =
                        (type == SpatialDomains::eRegular ||
                         type == SpatialDomains::eMovingRegular) ? 1 : nqtot;
            const MetricType m[3][3] = { {MetricLaplacian00, MetricLaplacian01, MetricLaplacian02},
                                       {MetricLaplacian01, MetricLaplacian11, MetricLaplacian12},
                                       {MetricLaplacian02, MetricLaplacian12, MetricLaplacian22}
            };

            for (unsigned int i = 0; i < dim; ++i)
            {
                for (unsigned int j = i; j < dim; ++j)
                {
                    m_metrics[m[i][j]] = Array<OneD, NekDouble>(nqtot);
                    Vmath::Vcopy(nqtot, &m_metricinfo->GetGmat()[i*dim+j][0], 1,
                                        &m_metrics[m[i][j]][0], 1);
                    MultiplyByQuadratureMetric(m_metrics[m[i][j]],
                                               m_metrics[m[i][j]]);

                }
            }
        }

        void Expansion::ComputeQuadratureMetric()
        {
            unsigned int nqtot = GetTotPoints();
            SpatialDomains::GeomType type = m_metricinfo->GetGtype();
            if (type == SpatialDomains::eRegular ||
                   type == SpatialDomains::eMovingRegular)
            {
                m_metrics[MetricQuadrature] = Array<OneD, NekDouble>(nqtot, m_metricinfo->GetJac()[0]);
            }
            else
            {
                m_metrics[MetricQuadrature] = m_metricinfo->GetJac();
            }

            MultiplyByStdQuadratureMetric(m_metrics[MetricQuadrature],
                                                   m_metrics[MetricQuadrature]);
        }

    } //end of namespace
} //end of namespace

