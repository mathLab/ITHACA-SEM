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

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/Interp.h>
#include <LocalRegions/Expansion.h>
#include <LocalRegions/MatrixKey.h>

using namespace std;

namespace Nektar
{
    namespace LocalRegions
    {
        Expansion::Expansion(SpatialDomains::GeometrySharedPtr pGeom) :
            m_IndexMapManager
               (std::bind(&Expansion::CreateIndexMap,this, std::placeholders::_1),
                std::string("ExpansionIndexMap")),
            m_geom(pGeom),
            m_metricinfo(m_geom->GetGeomFactors()),
            m_elementTraceLeft(-1),
            m_elementTraceRight(-1)
        {
            if (!m_metricinfo)
            {
                return;
            }

            if (!m_metricinfo->IsValid())
            {
                int nDim = m_base.size();
                string type = "regular";
                if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                {
                    type = "deformed";
                }

                stringstream err;
                err << nDim << "D " << type << " Jacobian not positive "
                    << "(element ID = " << m_geom->GetGlobalID() << ") "
                    << "(first vertex ID = " << m_geom->GetVid(0) << ")";
                NEKERROR(ErrorUtil::ewarning, err.str());
            }
        }

        Expansion::Expansion(const Expansion &pSrc) :
                StdExpansion(pSrc),
                m_IndexMapManager(pSrc.m_IndexMapManager),
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

        DNekMatSharedPtr Expansion::BuildTransformationMatrix(
            const DNekScalMatSharedPtr &r_bnd,
            const StdRegions::MatrixType matrixType)
        {
            return v_BuildTransformationMatrix(r_bnd,matrixType);
        }


        DNekMatSharedPtr Expansion::BuildVertexMatrix(
            const DNekScalMatSharedPtr &r_bnd)
        {
            return v_BuildVertexMatrix(r_bnd);
        }


        void Expansion::ExtractDataToCoeffs(
            const NekDouble *data,
            const std::vector<unsigned int > &nummodes,
            const int nmodes_offset,
            NekDouble *coeffs,
            std::vector<LibUtilities::BasisType> &fromType)
        {
            v_ExtractDataToCoeffs(data,nummodes,nmodes_offset,coeffs,fromType);
        }

        void Expansion::AddEdgeNormBoundaryInt(
            const int                           edge,
            const std::shared_ptr<Expansion>   &EdgeExp,
            const Array<OneD, const NekDouble> &Fx,
            const Array<OneD, const NekDouble> &Fy,
                  Array<OneD,       NekDouble> &outarray)
        {
            v_AddEdgeNormBoundaryInt(edge, EdgeExp, Fx, Fy, outarray);
        }

        void Expansion::AddEdgeNormBoundaryInt(
            const int                           edge,
            const std::shared_ptr<Expansion>   &EdgeExp,
            const Array<OneD, const NekDouble> &Fn,
                  Array<OneD,       NekDouble> &outarray)
        {
            v_AddEdgeNormBoundaryInt(edge, EdgeExp, Fn, outarray);
        }

        void Expansion::AddFaceNormBoundaryInt(
            const int                           face,
            const std::shared_ptr<Expansion>   &FaceExp,
            const Array<OneD, const NekDouble> &Fn,
                  Array<OneD,       NekDouble> &outarray)
        {
            v_AddFaceNormBoundaryInt(face, FaceExp, Fn, outarray);
        }

        void Expansion::DGDeriv(
            const int                                   dir,
            const Array<OneD, const NekDouble>&         inarray,
                  Array<OneD, ExpansionSharedPtr>      &EdgeExp,
                  Array<OneD, Array<OneD, NekDouble> > &coeffs,
                  Array<OneD,             NekDouble>   &outarray)
        {
            v_DGDeriv(dir, inarray, EdgeExp, coeffs, outarray);
        }

        NekDouble Expansion::VectorFlux(
            const Array<OneD, Array<OneD, NekDouble > > &vec)
        {
            return v_VectorFlux(vec);
        }

        DNekScalMatSharedPtr Expansion::GetLocMatrix(const StdRegions::MatrixType mtype,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varcoeffs)
        {
            MatrixKey mkey(mtype, DetShapeType(), *this, factors, varcoeffs);
            return GetLocMatrix(mkey);
        }

        SpatialDomains::GeometrySharedPtr Expansion::GetGeom() const
        {
            return m_geom;
        }

        void Expansion::Reset()
        {
            // Clear metrics
            m_metrics.clear();

            // Regenerate geometry factors
            m_metricinfo = m_geom->GetGeomFactors();
        }

        IndexMapValuesSharedPtr Expansion::CreateIndexMap(const IndexMapKey &ikey)
        {
            IndexMapValuesSharedPtr returnval;

            IndexMapType itype = ikey.GetIndexMapType();

            int entity = ikey.GetIndexEntity();

            StdRegions::Orientation orient = ikey.GetIndexOrientation();

            Array<OneD,unsigned int>     map;
            Array<OneD,int>             sign;

            switch(itype)
            {
            case eEdgeToElement:
                {
                    GetTraceToElementMap(entity,map,sign,orient);
                }
                break;
            case eFaceToElement:
                {
                    GetTraceToElementMap(entity,map,sign,orient);
                }
                break;
            case eEdgeInterior:
                {
                    ASSERTL0(false,"Boundary Index Map not implemented yet.");
                    //v_GetEdgeInteriorMap(entity,orient,map,sign);
                }
                break;
            case eFaceInterior:
                {
                    ASSERTL0(false,"Boundary Index Map not implemented yet.");
                    //v_GetFaceInteriorMap(entity,orient,map,sign);
                }
                break;
            case eBoundary:
                {
                    ASSERTL0(false,"Boundary Index Map not implemented yet.");
                }
                break;
            case eVertex:
                {
                    ASSERTL0(false,"Vertex Index Map not implemented yet.");
                }
                break;
            default:
                {
                    ASSERTL0(false,"The Index Map you are requiring "
                             "is not between the possible options.");
                }
            }
            
            returnval = MemoryManager<IndexMapValues>::AllocateSharedPtr(map.size());
            
            for(int i = 0; i < map.size(); i++)
            {
                (*returnval)[i].index =  map[i];
                (*returnval)[i].sign  =  sign[i];
            }

            return returnval;
        }

        const SpatialDomains::GeomFactorsSharedPtr& Expansion::GetMetricInfo() const
        {
            return m_metricinfo;
        }


        DNekScalMatSharedPtr Expansion::v_GetLocMatrix(const LocalRegions::MatrixKey &mkey)
        {
            boost::ignore_unused(mkey);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            return NullDNekScalMatSharedPtr;
        }

        void Expansion::v_MultiplyByQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
                                                 Array<OneD, NekDouble> &outarray)
        {
            const int nqtot = GetTotPoints();

            if (m_metrics.count(eMetricQuadrature) == 0)
            {
                ComputeQuadratureMetric();
            }

            Vmath::Vmul(nqtot, m_metrics[eMetricQuadrature], 1, inarray, 1, outarray, 1);
        }

        void Expansion::v_DivideByQuadratureMetric(
            const Array<OneD, 
            const NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            const int nqtot = GetTotPoints();

            if (m_metrics.count(eMetricQuadrature) == 0)
            {
                ComputeQuadratureMetric();
            }

            Vmath::Vdiv(nqtot, inarray, 
                        1, m_metrics[eMetricQuadrature], 
                        1, outarray, 1);
        }

        void Expansion::ComputeLaplacianMetric()
        {
            v_ComputeLaplacianMetric();
        }

        void Expansion::ComputeQuadratureMetric()
        {
            unsigned int nqtot = GetTotPoints();
            SpatialDomains::GeomType type = m_metricinfo->GetGtype();
            LibUtilities::PointsKeyVector p = GetPointsKeys();
            if (type == SpatialDomains::eRegular ||
                   type == SpatialDomains::eMovingRegular)
            {
                m_metrics[eMetricQuadrature] = Array<OneD, NekDouble>(nqtot, m_metricinfo->GetJac(p)[0]);
            }
            else
            {
                m_metrics[eMetricQuadrature] = m_metricinfo->GetJac(p);
            }

            v_MultiplyByStdQuadratureMetric(m_metrics[eMetricQuadrature],
                                                   m_metrics[eMetricQuadrature]);
        }

        void Expansion::v_GetCoords(
            Array<OneD, NekDouble> &coords_0,
            Array<OneD, NekDouble> &coords_1,
            Array<OneD, NekDouble> &coords_2)
        {
            ASSERTL1(m_geom, "m_geom not defined");

            // get physical points defined in Geom
            m_geom->FillGeom();

            const int expDim = m_base.size();
            int       nqGeom = 1;
            bool      doCopy = true;

            Array<OneD, LibUtilities::BasisSharedPtr> CBasis(expDim);
            Array<OneD, Array<OneD, NekDouble> > tmp(3);

            for (int i = 0; i < expDim; ++i)
            {
                CBasis[i] = m_geom->GetXmap()->GetBasis(i);
                nqGeom   *= CBasis[i]->GetNumPoints();
                doCopy    = doCopy && m_base[i]->GetBasisKey().SamePoints(
                                                      CBasis[i]->GetBasisKey());
            }

            tmp[0] = coords_0;
            tmp[1] = coords_1;
            tmp[2] = coords_2;

            if (doCopy)
            {
                for (int i = 0; i < m_geom->GetCoordim(); ++i)
                {
                    m_geom->GetXmap()->BwdTrans(m_geom->GetCoeffs(i), tmp[i]);
                }
            }
            else
            {
                for (int i = 0; i < m_geom->GetCoordim(); ++i)
                {
                    Array<OneD, NekDouble> tmpGeom(nqGeom);
                    m_geom->GetXmap()->BwdTrans(m_geom->GetCoeffs(i), tmpGeom);

                    switch (expDim)
                    {
                        case 1:
                        {
                            LibUtilities::Interp1D(
                                CBasis[0]->GetPointsKey(), &tmpGeom[0],
                                m_base[0]->GetPointsKey(), &tmp[i][0]);
                            break;
                        }
                        case 2:
                        {
                            LibUtilities::Interp2D(
                                CBasis[0]->GetPointsKey(),
                                CBasis[1]->GetPointsKey(),
                                &tmpGeom[0],
                                m_base[0]->GetPointsKey(),
                                m_base[1]->GetPointsKey(),
                                &tmp[i][0]);
                            break;
                        }
                        case 3:
                        {
                            LibUtilities::Interp3D(
                                CBasis[0]->GetPointsKey(),
                                CBasis[1]->GetPointsKey(),
                                CBasis[2]->GetPointsKey(),
                                &tmpGeom[0],
                                m_base[0]->GetPointsKey(),
                                m_base[1]->GetPointsKey(),
                                m_base[2]->GetPointsKey(),
                                &tmp[i][0]);
                            break;
                        }
                    }
                }
            }
        }

        void Expansion::ComputeGmatcdotMF(
                const Array<TwoD, const NekDouble> &df,
                const Array<OneD, const NekDouble> &direction,
                Array<OneD, Array<OneD, NekDouble> > &dfdir)
        {
            int shapedim = dfdir.size();
            int coordim = GetCoordim();
            int nqtot = direction.size()/coordim;

            for(int j = 0; j < shapedim; j++)
            {
                dfdir[j] = Array<OneD, NekDouble>(nqtot,0.0);
                for (int k = 0; k < coordim; k++)
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vvtvp(nqtot,
                                     &df[shapedim*k+j][0], 1,
                                     &direction[k*nqtot],  1,
                                     &dfdir[j][0],         1,
                                     &dfdir[j][0],         1);
                    }
                    else
                    {
                        Vmath::Svtvp(nqtot,
                                     df[shapedim*k+j][0],
                                     &direction[k*nqtot],  1,
                                     &dfdir[j][0],         1,
                                     &dfdir[j][0],         1);
                    }
                }
            }
        }

        // Get Moving frames
        Array<OneD, NekDouble> Expansion::v_GetMF(
                const int dir,
                const int shapedim,
                const StdRegions::VarCoeffMap   &varcoeffs)
        {
            int coordim = GetCoordim();

            int nquad0, nquad1, nquad2;
            int nqtot=1;
            switch(shapedim)
            {
                case 2:
                {
                    nquad0  = m_base[0]->GetNumPoints();
                    nquad1  = m_base[1]->GetNumPoints();
                    nqtot   = nquad0*nquad1;
                    break;
                }
                case 3:
                {
                    nquad0 = m_base[0]->GetNumPoints();
                    nquad1 = m_base[1]->GetNumPoints();
                    nquad2 = m_base[2]->GetNumPoints();
                    nqtot   = nquad0 * nquad1 * nquad2;
                    break;
                }
                default:
                    break;
            }

            StdRegions::VarCoeffType MMFCoeffs[15] =
            {
                StdRegions::eVarCoeffMF1x,
                StdRegions::eVarCoeffMF1y,
                StdRegions::eVarCoeffMF1z,
                StdRegions::eVarCoeffMF1Div,
                StdRegions::eVarCoeffMF1Mag,
                StdRegions::eVarCoeffMF2x,
                StdRegions::eVarCoeffMF2y,
                StdRegions::eVarCoeffMF2z,
                StdRegions::eVarCoeffMF2Div,
                StdRegions::eVarCoeffMF2Mag,
                StdRegions::eVarCoeffMF3x,
                StdRegions::eVarCoeffMF3y,
                StdRegions::eVarCoeffMF3z,
                StdRegions::eVarCoeffMF3Div,
                StdRegions::eVarCoeffMF3Mag
            };

            Array<OneD, NekDouble> MF(coordim*nqtot);
            Array<OneD, NekDouble> tmp(nqtot);
            for (int k = 0; k < coordim; k++)
            {
                StdRegions::VarCoeffMap::const_iterator MFdir =
                        varcoeffs.find(MMFCoeffs[5*dir+k]);
                tmp = MFdir->second;

                Vmath::Vcopy(nqtot, &tmp[0], 1, &MF[k*nqtot], 1);
            }

            return MF;
        }

        // Get magnitude of MF
        Array<OneD, NekDouble> Expansion::v_GetMFDiv(
                const int dir,
                const StdRegions::VarCoeffMap   &varcoeffs)
        {
            int indxDiv = 3;

            StdRegions::VarCoeffType MMFCoeffs[15] =
            {
                StdRegions::eVarCoeffMF1x,
                StdRegions::eVarCoeffMF1y,
                StdRegions::eVarCoeffMF1z,
                StdRegions::eVarCoeffMF1Div,
                StdRegions::eVarCoeffMF1Mag,
                StdRegions::eVarCoeffMF2x,
                StdRegions::eVarCoeffMF2y,
                StdRegions::eVarCoeffMF2z,
                StdRegions::eVarCoeffMF2Div,
                StdRegions::eVarCoeffMF2Mag,
                StdRegions::eVarCoeffMF3x,
                StdRegions::eVarCoeffMF3y,
                StdRegions::eVarCoeffMF3z,
                StdRegions::eVarCoeffMF3Div,
                StdRegions::eVarCoeffMF3Mag
            };

            StdRegions::VarCoeffMap::const_iterator MFdir =
                    varcoeffs.find(MMFCoeffs[5*dir+indxDiv]);
            Array<OneD, NekDouble> MFDiv = MFdir->second;

            return MFDiv;
        }

        // Get magnitude of MF
        Array<OneD, NekDouble> Expansion::v_GetMFMag(
                const int dir,
                const StdRegions::VarCoeffMap   &varcoeffs)
        {
            int indxMag = 4;

            StdRegions::VarCoeffType MMFCoeffs[15] =
            {
                StdRegions::eVarCoeffMF1x,
                StdRegions::eVarCoeffMF1y,
                StdRegions::eVarCoeffMF1z,
                StdRegions::eVarCoeffMF1Div,
                StdRegions::eVarCoeffMF1Mag,
                StdRegions::eVarCoeffMF2x,
                StdRegions::eVarCoeffMF2y,
                StdRegions::eVarCoeffMF2z,
                StdRegions::eVarCoeffMF2Div,
                StdRegions::eVarCoeffMF2Mag,
                StdRegions::eVarCoeffMF3x,
                StdRegions::eVarCoeffMF3y,
                StdRegions::eVarCoeffMF3z,
                StdRegions::eVarCoeffMF3Div,
                StdRegions::eVarCoeffMF3Mag
            };

            StdRegions::VarCoeffMap::const_iterator MFdir = varcoeffs.find(MMFCoeffs[5*dir+indxMag]);
            Array<OneD, NekDouble> MFmag = MFdir->second;

            return MFmag;
        }


        DNekMatSharedPtr Expansion::v_BuildTransformationMatrix(
            const DNekScalMatSharedPtr &r_bnd,
            const StdRegions::MatrixType matrixType)
        {
            boost::ignore_unused(r_bnd, matrixType);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            return NullDNekMatSharedPtr;
        }

        DNekMatSharedPtr Expansion::v_BuildVertexMatrix(
            const DNekScalMatSharedPtr &r_bnd)
        {
            boost::ignore_unused(r_bnd);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            return NullDNekMatSharedPtr;
        }

        void Expansion::v_ExtractDataToCoeffs(
            const NekDouble *data,
            const std::vector<unsigned int > &nummodes,
            const int nmodes_offset,
            NekDouble *coeffs,
            std::vector<LibUtilities::BasisType> &fromType)
        {
            boost::ignore_unused(data, nummodes, nmodes_offset,
                                 coeffs, fromType);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
        }

        void Expansion::v_AddEdgeNormBoundaryInt(
            const int                           edge,
            const std::shared_ptr<Expansion>   &EdgeExp,
            const Array<OneD, const NekDouble> &Fx,
            const Array<OneD, const NekDouble> &Fy,
                  Array<OneD,       NekDouble> &outarray)
        {
            boost::ignore_unused(edge, EdgeExp, Fx, Fy, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
        }

        void Expansion::v_AddEdgeNormBoundaryInt(
            const int                           edge,
            const std::shared_ptr<Expansion>   &EdgeExp,
            const Array<OneD, const NekDouble> &Fn,
                  Array<OneD,       NekDouble> &outarray)
        {
            boost::ignore_unused(edge, EdgeExp, Fn, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
        }

        void Expansion::v_AddFaceNormBoundaryInt(
            const int                           face,
            const std::shared_ptr<Expansion>   &FaceExp,
            const Array<OneD, const NekDouble> &Fn,
                  Array<OneD,       NekDouble> &outarray)
        {
            boost::ignore_unused(face, FaceExp, Fn, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
        }

        void Expansion::v_DGDeriv(
            const int                                   dir,
            const Array<OneD, const NekDouble>&         inarray,
                  Array<OneD, ExpansionSharedPtr>      &EdgeExp,
                  Array<OneD, Array<OneD, NekDouble> > &coeffs,
                  Array<OneD,             NekDouble>   &outarray)
        {
            boost::ignore_unused(dir, inarray, EdgeExp, coeffs, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
        }

        NekDouble Expansion::v_VectorFlux(
            const Array<OneD, Array<OneD, NekDouble > > &vec)
        {
            boost::ignore_unused(vec);
            NEKERROR(ErrorUtil::efatal, "This function is only valid for LocalRegions");
            return 0.0;
        }

        StdRegions::Orientation Expansion::v_GetTraceOrient(int trace)
        {
            boost::ignore_unused(trace);
            return StdRegions::eForwards;
        }

        void Expansion::v_SetCoeffsToOrientation(StdRegions::Orientation dir,
                                                 Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(dir, inarray, outarray);
            NEKERROR(ErrorUtil::efatal, "This function is not defined for this shape");
        }

        void Expansion::v_GetTraceQFactors(const int trace,
                                          Array<OneD, NekDouble> &outarray)
        {
            boost::ignore_unused(trace, outarray);
            NEKERROR(ErrorUtil::efatal,
                     "Method does not exist for this shape or library");
        }
        
        void Expansion::v_GetTracePhysVals
                    (const int trace,
                     const StdRegions::StdExpansionSharedPtr &TraceExp,
                     const Array<OneD, const NekDouble> &inarray,
                     Array<OneD,       NekDouble> &outarray,
                     StdRegions::Orientation  orient)
        {
            boost::ignore_unused(trace,TraceExp,inarray,outarray,orient);
            NEKERROR(ErrorUtil::efatal,
                     "Method does not exist for this shape or library" );
        }

        void Expansion::v_GetTracePhysMap(const int  edge,
                                         Array<OneD, int>   &outarray)
        {
            boost::ignore_unused(edge, outarray);
            NEKERROR(ErrorUtil::efatal,
                     "Method does not exist for this shape or library" );
        }
        
        void Expansion::v_ReOrientTracePhysMap
                                (const StdRegions::Orientation orient,
                                 Array<OneD, int> &idmap,
                                 const int nq0,  const int nq1)
        {
            boost::ignore_unused(orient,idmap,nq0,nq1);
            NEKERROR(ErrorUtil::efatal,
                     "Method does not exist for this shape or library" );
        }
 
        const NormalVector & Expansion::v_GetTraceNormal(const int id) const
        {
            boost::ignore_unused(id);
            ASSERTL0(false, "Cannot get trace normals for this expansion.");
            static NormalVector result;
            return result;
        }

        void Expansion::v_ComputeTraceNormal(const int id)
        {
            boost::ignore_unused(id);
            ASSERTL0(false, "Cannot compute trace normal for this expansion.");
        }

        const Array<OneD, const NekDouble>& Expansion::v_GetPhysNormals(void)
        {
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
            return NullNekDouble1DArray;
        }


        void Expansion::v_SetPhysNormals(Array<OneD, const NekDouble> &normal)
        {
            boost::ignore_unused(normal);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        void Expansion::v_SetUpPhysNormals(const int edge)
        {
            boost::ignore_unused(edge);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        void Expansion::v_AddRobinMassMatrix(
            const int                            trace,
            const Array<OneD, const NekDouble > &primCoeffs,
                  DNekMatSharedPtr              &inoutmat)
        {
            boost::ignore_unused(trace,primCoeffs,inoutmat);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

        void Expansion::v_AddRobinTraceContribution(
                const int traceid,
                const Array<OneD, const NekDouble> &primCoeffs,
                const Array<OneD, NekDouble> &incoeffs,
                Array<OneD, NekDouble> &coeffs)
        {
            boost::ignore_unused(traceid,primCoeffs,incoeffs, coeffs);
            NEKERROR(ErrorUtil::efatal, "This function is not valid for this class");
        }

    } //end of namespace
} //end of namespace

