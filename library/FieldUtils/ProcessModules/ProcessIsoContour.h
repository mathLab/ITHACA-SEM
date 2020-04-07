///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessIsoContour.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: Generate isocontours from field data.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSISOCONTOUR
#define FIELDUTILS_PROCESSISOCONTOUR

#include "../Module.h"
#include "ProcessEquiSpacedOutput.h"

namespace Nektar
{
namespace FieldUtils
{

class Iso
{
    public:
        void  Condense(void);
        void  GlobalCondense(std::vector<std::shared_ptr<Iso> > &iso, bool verbose);
        void  SeparateRegions(std::vector<std::shared_ptr<Iso> > &iso, int minsize, bool verbose);

        void  Smooth(int n_iter, NekDouble lambda, NekDouble mu);

        int  GetNVert(void)
        {
            return m_nvert;
        }

        void SetNVert(int n)
        {
            m_nvert = n;
        }

        int  GetNTris(void)
        {
            return m_ntris;
        }

        void SetNTris(int n)
        {
            m_ntris = n;
        }

        void SetFields(const int loc,
                        const Array<OneD,Array<OneD, NekDouble> > &intfields,
                        const int j)
        {
            m_x[loc] = intfields[0][j];
            m_y[loc] = intfields[1][j];
            m_z[loc] = intfields[2][j];

            for(int i = 0; i < intfields.size()-3; ++i)
            {
                m_fields[i][loc] = intfields[i+3][j];
            }
        }

        NekDouble GetFields(const int i, const int j)
        {
            return m_fields[i][j];
        }

        void SetX(int loc, NekDouble val)
        {
            m_x[loc] = val;
        }

        void SetY(int loc, NekDouble val)
        {
            m_y[loc] = val;
        }

        void SetZ(int loc, NekDouble val)
        {
            m_z[loc] = val;
        }

        NekDouble GetX(int loc)
        {
            return m_x[loc];
        }

        NekDouble GetY(int loc)
        {
            return m_y[loc];
        }

        NekDouble GetZ(int loc)
        {
            return m_z[loc];
        }

        int  GetVId(int i)
        {
            return m_vid[i];
        }

        void ResizeVId(int nconn)
        {
            m_vid = Array<OneD, int>(nconn);
        }

        void SetVId(int i, int j)
        {
            m_vid[i] = j;
        }

        void ResizeFields(int size)
        {
            if(size > m_x.size()) // add 1000 element to vectors
            {
                m_x.resize(size+100);
                m_y.resize(size+100);
                m_z.resize(size+100);;
                for(int i = 0; i < m_fields.size(); ++i)
                {
                    m_fields[i].resize(size+1000);
                }

            }
            m_nvert = size;
        }

        Iso(int nfields)
        {
            m_condensed = false;
            m_nvert     = 0;
            m_fields.resize(nfields);
            // set up initial vectors to be 1000 long
            m_x.resize(1000);
            m_y.resize(1000);
            m_z.resize(1000);
            for(int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i].resize(1000);
            }
        };

        ~Iso(void)
        {
        }

    private:
        bool                                 m_condensed;
        int                                  m_nvert;           // number of vertices
        int                                  m_ntris;           // number of triangles introduced.
        std::vector<NekDouble>               m_x;
        std::vector<NekDouble>               m_y;
        std::vector<NekDouble>               m_z;
        std::vector<std::vector<NekDouble> > m_fields;
        Array<OneD, int>                     m_vid; // used when condensing field

};

typedef std::shared_ptr<Iso> IsoSharedPtr;

class IsoVertex
{
    public:
        friend class Iso;

        IsoVertex (void)
        {
            m_id = -1;
            m_x = m_y = m_z = -99999;
        }

        ~IsoVertex(){};

        int get_iso_id()
        {
            return m_iso_id;
        }

        int get_iso_vert_id()
        {
            return m_iso_vert_id;
        }

        friend bool operator == (const IsoVertex& x, const IsoVertex& y);
        friend bool operator != (const IsoVertex& x, const IsoVertex& y);

    private:
        int                      m_id;
        int                      m_iso_id;
        int                      m_iso_vert_id;
        NekDouble                m_x, m_y, m_z;
        std::vector<NekDouble >  m_fields;

};

/**
 * @brief This processing module extracts an isocontour
 */
class ProcessIsoContour : public ProcessModule
{
    public:
        /// Creates an instance of this class
        static std::shared_ptr<Module> create(FieldSharedPtr f)
        {
            return MemoryManager<ProcessIsoContour>::AllocateSharedPtr(f);
        }
        static ModuleKey className;

        ProcessIsoContour(FieldSharedPtr f);
        virtual ~ProcessIsoContour();

        /// Write mesh to output file.
        virtual void Process(po::variables_map &vm);

        virtual std::string GetModuleName()
        {
            return "ProcessIsoContour";
        }

        virtual std::string GetModuleDescription()
        {
            return "Extracting contour";
        }

        virtual ModulePriority GetModulePriority()
        {
            return eModifyPts;
        }

    protected:
        ProcessIsoContour(){};
        void ResetFieldPts(std::vector<IsoSharedPtr> &iso);
        void SetupIsoFromFieldPts(std::vector<IsoSharedPtr> &isovec);

    private:

        std::vector<IsoSharedPtr> ExtractContour(
            const int fieldid,
            const NekDouble val);
};

}
}

#endif
