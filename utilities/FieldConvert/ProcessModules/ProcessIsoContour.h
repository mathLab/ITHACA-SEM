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
//  Description: Generate isocontours from field data.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_PREPROCESSING_FIELDCONVERT_PROCESSISOCONTOUR
#define UTILITIES_PREPROCESSING_FIELDCONVERT_PROCESSISOCONTOUR

#include "../Module.h"
#include "ProcessEquiSpacedOutput.h"

namespace Nektar
{
namespace Utilities
{

class Iso
{
    public:
        void  condense(void);
        void  globalcondense(vector<boost::shared_ptr<Iso> > &iso);
        void  separate_regions(vector<boost::shared_ptr<Iso> > &iso, int minsize);

        void  smooth(int n_iter, NekDouble lambda, NekDouble mu);

        int  get_nvert(void)
        {
            return m_nvert;
        }

        void set_nvert(int n)
        {
            m_nvert = n;
        }

        int  get_ntris(void)
        {
            return m_ntris;
        }

        void set_ntris(int n)
        {
            m_ntris = n;
        }

        void set_fields(const int loc,
                        const Array<OneD,Array<OneD, NekDouble> > &intfields,
                        const int j)
        {
            m_x[loc] = intfields[0][j];
            m_y[loc] = intfields[1][j];
            m_z[loc] = intfields[2][j];

            for(int i = 0; i < intfields.num_elements()-3; ++i)
            {
                m_fields[i][loc] = intfields[i+3][j];
            }
        }

        NekDouble get_fields(const int i, const int j)
        {
            return m_fields[i][j];
        }
        void set_x(int loc, NekDouble val)
        {
            m_x[loc] = val;
        }

        void set_y(int loc, NekDouble val)
        {
            m_y[loc] = val;
        }

        void set_z(int loc, NekDouble val)
        {
            m_z[loc] = val;
        }

        NekDouble get_x(int loc)
        {
            return m_x[loc];
        }

        NekDouble get_y(int loc)
        {
            return m_y[loc];
        }

        NekDouble get_z(int loc)
        {
            return m_z[loc];
        }

        int  get_vid(int i)
        {
            return m_vid[i];
        }

        void resize_vid(int nconn)
        {
            m_vid = Array<OneD, int>(nconn);
        }

        void set_vid(int i, int j)
        {
            m_vid[i] = j;
        }

        void resize_fields(int size)
        {
            if(size > m_x.size()) // add 1000 element to vectors
            {
                m_x.resize(size+1000);
                m_y.resize(size+1000);
                m_z.resize(size+1000);;
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
            // set up initial vectors to be 10000 long
            m_x.resize(10000);
            m_y.resize(10000);
            m_z.resize(10000);
            for(int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i].resize(10000);
            }
        };

        ~Iso(void)
        {
        }

    private:
        bool                       m_condensed;
        int                        m_nvert;           // number of vertices
        int                        m_ntris;           // number of triangles introduced.
        vector<NekDouble>          m_x;
        vector<NekDouble>          m_y;
        vector<NekDouble>          m_z;
        vector<vector<NekDouble> > m_fields;
        Array<OneD, int>           m_vid; // used when condensing field

};

typedef boost::shared_ptr<Iso> IsoSharedPtr;

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
        int                 m_id;
        int                 m_iso_id;
        int                 m_iso_vert_id;
        NekDouble           m_x, m_y, m_z;
        vector<NekDouble >  m_fields;

};

/**
 * @brief This processing module extracts an isocontour
 */
class ProcessIsoContour : public ProcessEquiSpacedOutput
{
    public:
        /// Creates an instance of this class
        static boost::shared_ptr<Module> create(FieldSharedPtr f)
        {
            return MemoryManager<ProcessIsoContour>::AllocateSharedPtr(f);
        }
        static ModuleKey className;

        ProcessIsoContour(FieldSharedPtr f);
        virtual ~ProcessIsoContour();

        /// Write mesh to output file.
        virtual void Process(po::variables_map &vm);

    protected:
        ProcessIsoContour(){};
        void ResetFieldPts(vector<IsoSharedPtr> &iso);

    private:

        vector<IsoSharedPtr> ExtractContour(
            const int fieldid,
            const NekDouble val);
};

}
}

#endif
